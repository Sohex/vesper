#!/usr/bin/env python3
"""The land water ledger: what balances, over what control volume, and the
check that fails when it does not.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a declared graph of control volumes, the water transfers
between them, and the ownership rules that decide which component may advance
each store.

    python hydrography/scripts/land_water_ledger.py            # status, exit 0
    python hydrography/scripts/land_water_ledger.py --strict   # exit 1 while the ledger does not close

This module is the enforcement for `hydrography/config/land_water_ledger.yaml`.
It does four things, in this order, because each depends on the one before it:

1. **The graph.** Every term moves water from one declared node to one declared
   node. A dangling destination, a term drawing on a terminal node, a reservoir
   with no outflow, a phase change with no process named, a phase change
   carrying the wrong latent heat, or an evaporative component debited from two
   stores is an error.

2. **Closure, on fixtures that can fail.** A ledger that only ever runs on real
   inputs has no failing case, so the balance is exercised on reduced synthetic
   domains where the answer is known in advance. Ten of the eleven fixtures are
   built to be WRONG in a named way and are required to be rejected; a fixture
   that does not get the verdict it was built for is a defect in this checker
   and exits non-zero even when nothing else is wrong.

3. **Ownership.** Each store has one owner and one symbol, and every symbol in
   another component that presently holds an independent copy of the same store
   is declared as a shadow copy. A shadow copy is not a second store; it is the
   same water counted twice, and it is what lets two models evaporate one
   millimetre.

4. **Timing.** Each term declares the shortest interval its producing artifact
   can supply it on. A term reaching a store that changes inside the reporting
   interval, produced only as an annual mean, cannot partition that store.

The report goes to `hydrography/analysis/land_water_ledger_report.json`.
The contract is `hydrography/notes/land-water-ledger.md`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import ANALYSIS, PROJECT_ROOT  # noqa: F401  (adds lib/ to sys.path)

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = COMPONENT_ROOT / "config" / "land_water_ledger.yaml"
REPORT = ANALYSIS / "land_water_ledger_report.json"

UNDECLARED = "undeclared"

# The three phases one substance takes. The ledger books H2O and nothing else,
# so a phase is a property of a term and never a different quantity.
PHASES = ("liquid", "ice", "vapour")

# What a term may declare it does to phase, and the latent heat that change
# actually consumes or releases. `none` moves water unchanged; every other
# transform MUST change the phase and MUST name the matching latent heat.
# A melt term carrying the latent heat of vaporisation conserves water and
# loses energy, which is the class of error this mapping exists to catch.
TRANSFORM_LATENT = {
    "none": "none",
    "melting": "fusion",
    "freezing": "fusion",
    "vaporisation": "vaporisation",
    "condensation": "vaporisation",
    "sublimation": "sublimation",
    "deposition": "sublimation",
}

# The phase pairs each transform is allowed to move between.
TRANSFORM_PHASES = {
    "melting": ("ice", "liquid"),
    "freezing": ("liquid", "ice"),
    "vaporisation": ("liquid", "vapour"),
    "condensation": ("vapour", "liquid"),
    "sublimation": ("ice", "vapour"),
    "deposition": ("vapour", "ice"),
}

# The shortest interval a producing artifact can supply a term on.
INTERVAL_FLOORS = ("instantaneous", "daily", "monthly", "annual", UNDECLARED)


def load(path: Path = DECLARATION) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# 1. The graph
# ---------------------------------------------------------------------------

def check_graph(decl: dict) -> list[str]:
    """Structural failures in the declared ledger. Empty means the graph closes.

    Every one of these has a failing case that is a real modelling error rather
    than a typo, which is the test for whether a check is worth having.
    """
    nodes = decl["nodes"]
    bad: list[str] = []
    inflow: dict[str, int] = {name: 0 for name in nodes}
    outflow: dict[str, int] = {name: 0 for name in nodes}
    # et_component -> [(term, source node)]. One component, one debit.
    debits: dict[str, list[tuple[str, str]]] = {}

    for name, term in decl["terms"].items():
        src, dst = term.get("from"), term.get("to")
        for end, where in ((src, "from"), (dst, "to")):
            if end not in nodes:
                bad.append(f"term {name}: {where} names {end!r}, which is not a declared node")
        if src not in nodes or dst not in nodes:
            continue
        outflow[src] += 1
        inflow[dst] += 1

        if src == dst:
            bad.append(f"term {name}: source and destination are the same node")
        if nodes[src]["kind"] == "terminal":
            bad.append(f"term {name}: draws on {src}, which is terminal. "
                       "A terminal node is where this ledger stops")
        if nodes[dst]["kind"] == "boundary" and not nodes[dst].get("reentrant"):
            bad.append(f"term {name}: delivers to {dst}, which is a boundary and "
                       "is not declared reentrant. Water leaving the domain has a "
                       "terminal node to go to")

        pf, pt = term.get("phase_from"), term.get("phase_to")
        for phase, where in ((pf, "phase_from"), (pt, "phase_to")):
            if phase not in PHASES:
                bad.append(f"term {name}: {where} names {phase!r}, which is not a declared phase")
        transform = term.get("transform")
        if transform not in TRANSFORM_LATENT:
            bad.append(f"term {name}: unknown transform {transform!r}")
        else:
            if transform == "none" and pf != pt:
                bad.append(f"term {name}: transform is none but the phase changes "
                           f"{pf} to {pt}. A phase change is a process, it costs "
                           "energy, and it needs naming")
            if transform != "none":
                if pf == pt:
                    bad.append(f"term {name}: transform is {transform} and the phase "
                               "does not change")
                elif TRANSFORM_PHASES[transform] != (pf, pt):
                    want = TRANSFORM_PHASES[transform]
                    bad.append(f"term {name}: transform {transform} moves "
                               f"{want[0]} to {want[1]}, not {pf} to {pt}")
            want_latent = TRANSFORM_LATENT[transform]
            if term.get("latent") != want_latent:
                bad.append(f"term {name}: transform {transform} releases or consumes "
                           f"the latent heat of {want_latent}, and the term declares "
                           f"{term.get('latent')!r}. A phase term carrying the wrong "
                           "latent heat conserves water and loses energy")

        floor = term.get("interval_floor")
        if floor not in INTERVAL_FLOORS:
            bad.append(f"term {name}: interval_floor {floor!r} is not one of "
                       f"{', '.join(INTERVAL_FLOORS)}")

        component = term.get("et_component")
        if component is not None:
            debits.setdefault(component, []).append((name, src))

    for component, entries in sorted(debits.items()):
        if len(entries) > 1:
            where = ", ".join(f"{t} from {s}" for t, s in entries)
            bad.append(f"et_component {component}: debited by {len(entries)} terms "
                       f"({where}). One evaporative component, one debit, one store")

    for name, node in nodes.items():
        kind = node["kind"]
        if kind == "reservoir":
            if outflow[name] == 0:
                bad.append(f"node {name}: a reservoir with no outflow is an accumulator")
            if inflow[name] == 0:
                bad.append(f"node {name}: a reservoir with no inflow cannot be filled")
        if kind == "boundary" and inflow[name] and not node.get("reentrant"):
            bad.append(f"node {name}: a boundary that is not declared reentrant "
                       "receives nothing from inside the domain")
        if node.get("area_basis") is None:
            bad.append(f"node {name}: no area_basis. A per-area depth has no meaning "
                       "without one, and three different denominators appear here")
        if kind in ("reservoir", "interface"):
            if not node.get("owner"):
                bad.append(f"node {name}: a store with no owner. Exactly one component "
                           "advances each store")
            if node.get("shadow_copies") is None:
                bad.append(f"node {name}: shadow_copies is not declared. An empty list "
                           "is a claim; a missing key is silence")

    heat = decl["latent_heat"]
    residual = (heat["sublimation_j_per_kg"]
                - heat["vaporisation_j_per_kg"] - heat["fusion_j_per_kg"])
    if abs(residual) > float(heat["identity_tolerance_j_per_kg"]):
        bad.append(f"latent_heat: sublimation - vaporisation - fusion = "
                   f"{residual:.6g} J/kg, outside the declared tolerance. "
                   "ExoPlaSim derives fusion as the difference, so a fusion "
                   "constant that disagrees is a fourth number")
    return bad


def undeclared_terms(decl: dict) -> list[tuple[str, str]]:
    """(term, owning issue) for every term still carrying the sentinel."""
    return [(name, term.get("owner", "unowned"))
            for name, term in decl["terms"].items()
            if term.get("flux") == UNDECLARED]


def shadow_copies(decl: dict) -> list[dict]:
    """Stores presently held independently by more than one component."""
    out = []
    for name, node in decl["nodes"].items():
        copies = node.get("shadow_copies") or []
        if copies:
            out.append({"node": name, "owner": node.get("owner"),
                        "symbol": node.get("symbol"), "shadow_copies": copies})
    return out


def interval_refusals(decl: dict) -> list[str]:
    """Terms whose producer cannot resolve the store they deliver to.

    An annual mean is a total, not a partition. A term produced only annually
    that feeds a store which empties and refills inside the year sets that
    store's timing from an artifact that has none.
    """
    fast = set(decl["sub_annual_stores"])
    out = []
    for name, term in decl["terms"].items():
        if term.get("to") in fast and term.get("interval_floor") == "annual":
            out.append(f"term {name}: delivers to {term['to']}, which changes inside "
                       f"the reporting interval, and its producer "
                       f"({term.get('producer')}) supplies it only as an annual "
                       "mean. A total is not a partition")
    return out


# ---------------------------------------------------------------------------
# 2. Closure, on fixtures that can fail
# ---------------------------------------------------------------------------

class LedgerError(RuntimeError):
    """A transfer the ledger refuses to book."""


class Ledger:
    """A reduced ledger over a declared node set, in kg of water.

    Everything is booked in ABSOLUTE mass. A per-area depth is converted by the
    area of the node it is stated on, and that area is an argument rather than
    an assumption, so a depth generated over the land fraction and summed over
    the whole cell produces a residual instead of a silent rescale.
    """

    def __init__(self, decl: dict):
        self.kinds = {name: node["kind"] for name, node in decl["nodes"].items()}
        self.stock = {name: 0.0 for name in self.kinds}
        self.external_in = 0.0
        self.terminal_out = 0.0
        # Every evaporative withdrawal booked, by component, so a second debit
        # of the same component is visible even when mass happens to close.
        self.et_debits: dict[str, list[str]] = {}

    def seed(self, node: str, kg: float) -> None:
        """Put water in a store before the interval starts."""
        if self.kinds.get(node) not in ("reservoir", "interface"):
            raise LedgerError(f"cannot seed {node}: only a store holds an initial stock")
        self.stock[node] += kg

    def book(self, frm: str, to: str, kg: float, *,
             phase_from: str = "liquid", phase_to: str = "liquid",
             transform: str = "none", et_component: str | None = None) -> None:
        if frm not in self.kinds or to not in self.kinds:
            raise LedgerError(f"transfer {frm} -> {to} names a node the ledger does not have")
        if self.kinds[frm] == "terminal":
            raise LedgerError(f"transfer draws on {frm}, which is terminal")
        if kg < 0:
            raise LedgerError("a transfer carries a non-negative mass; reverse the direction instead")
        if transform == "none" and phase_from != phase_to:
            raise LedgerError(
                f"transfer {frm} -> {to} changes phase {phase_from} to {phase_to} "
                "with no process named. A phase change costs energy")
        if transform != "none" and TRANSFORM_PHASES.get(transform) != (phase_from, phase_to):
            raise LedgerError(
                f"transfer {frm} -> {to} declares {transform} and moves "
                f"{phase_from} to {phase_to}")
        if et_component is not None:
            seen = self.et_debits.setdefault(et_component, [])
            if seen:
                raise LedgerError(
                    f"evaporative component {et_component} is debited a second time, "
                    f"from {frm}, having already been debited from {seen[0]}. "
                    "One component, one debit, one store")
            seen.append(frm)

        if self.kinds[frm] == "boundary":
            self.external_in += kg
        else:
            if self.stock[frm] - kg < -1e-12:
                raise LedgerError(
                    f"transfer of {kg:.6g} kg leaves {frm} at "
                    f"{self.stock[frm] - kg:.6g} kg, below zero")
            self.stock[frm] -= kg

        if self.kinds[to] == "terminal":
            self.terminal_out += kg
        else:
            self.stock[to] += kg

    def handoff(self, frm: str, to: str, kg_out: float, kg_in: float, **kw) -> None:
        """A transfer whose two sides are stated SEPARATELY.

        This is what a cross-component transfer actually looks like: one model
        writes a rate on its own denominator and its own interval, and another
        reads it and multiplies by an area or a duration of its own. When the
        two disagree the mass changes and nothing in either file says so, so the
        ledger books the debit and the credit as given and lets the difference
        appear as a residual.
        """
        self.book(frm, to, kg_out, **kw)
        if self.kinds[to] == "terminal":
            self.terminal_out += kg_in - kg_out
        else:
            self.stock[to] += kg_in - kg_out

    def residual(self, seeded_kg: float) -> float:
        held = sum(v for name, v in self.stock.items()
                   if self.kinds[name] != "terminal")
        return (held + self.terminal_out) - (seeded_kg + self.external_in)


# One land cell, one Earth year, chosen so the fixtures read in kg without
# arithmetic. A T42 land cell is of order 1e10 m2; 1 mm over it is 1e7 kg.
CELL_AREA_M2 = 1.0e10
MM_KG = CELL_AREA_M2                      # kg of water per mm depth over the cell
EARTH_YEAR_S = 365.2422 * 86400.0


def _fixture_balanced(decl):
    """Rain and snow in, through canopy and pack and soil, out at the coast.

    Every store is left with a stock, every evaporative component is debited
    once, and the domain closes exactly.
    """
    led = Ledger(decl)
    led.seed("soil_liquid", 100.0 * MM_KG)
    led.seed("snowpack", 20.0 * MM_KG)
    led.seed("groundwater", 5000.0 * MM_KG)
    led.seed("open_water", 10.0 * MM_KG)
    seeded = 5130.0 * MM_KG

    led.book("atmosphere", "canopy", 400.0 * MM_KG)
    led.book("canopy", "atmosphere", 30.0 * MM_KG, phase_to="vapour",
             transform="vaporisation", et_component="canopy_evaporation")
    led.book("canopy", "soil_liquid", 370.0 * MM_KG)
    led.book("atmosphere", "snowpack", 80.0 * MM_KG,
             phase_from="ice", phase_to="ice")
    led.book("snowpack", "atmosphere", 6.0 * MM_KG, phase_from="ice",
             phase_to="vapour", transform="sublimation",
             et_component="snow_sublimation")
    led.book("snowpack", "soil_liquid", 80.0 * MM_KG, phase_from="ice",
             transform="melting")
    led.book("soil_liquid", "atmosphere", 120.0 * MM_KG, phase_to="vapour",
             transform="vaporisation", et_component="soil_evaporation")
    led.book("soil_liquid", "vegetation", 200.0 * MM_KG)
    led.book("vegetation", "atmosphere", 200.0 * MM_KG, phase_to="vapour",
             transform="vaporisation", et_component="transpiration")
    led.book("soil_liquid", "open_water", 90.0 * MM_KG)
    led.book("soil_liquid", "groundwater", 60.0 * MM_KG)
    led.book("groundwater", "open_water", 40.0 * MM_KG)
    led.book("open_water", "atmosphere", 25.0 * MM_KG, phase_to="vapour",
             transform="vaporisation", et_component="open_water_evaporation")
    led.book("open_water", "ocean_export", 100.0 * MM_KG)
    return abs(led.residual(seeded))


def _fixture_shadow_store_evaporation(decl):
    """Two models evaporating one soil, and only one of them owns a store.

    Finding 4. ExoPlaSim withdraws its own bare-soil evaporation from `dwatc`
    and the atmosphere is credited it. LPJ-GUESS transpires against `wcont`,
    which is not a node in the coupled domain, and the atmosphere is credited
    that too. The domain is debited once and credited twice, so the ledger
    reports the mass the coupled system does not have.
    """
    led = Ledger(decl)
    led.seed("soil_liquid", 300.0 * MM_KG)
    led.book("soil_liquid", "atmosphere", 120.0 * MM_KG, phase_to="vapour",
             transform="vaporisation", et_component="soil_evaporation")
    # LPJ's withdrawal: the credit is stated, the debit falls on a shadow store.
    led.handoff("soil_liquid", "atmosphere", kg_out=0.0, kg_in=180.0 * MM_KG,
                phase_to="vapour", transform="vaporisation")
    return abs(led.residual(300.0 * MM_KG))


def _fixture_double_debited_component(decl):
    """One evaporative component, two owners.

    The mass would close: both withdrawals come out of the same store. What
    fails is ownership, and it fails at the second booking rather than in the
    residual, because a component debited twice is wrong even when the
    arithmetic is not.
    """
    led = Ledger(decl)
    led.seed("soil_liquid", 300.0 * MM_KG)
    led.seed("canopy", 50.0 * MM_KG)
    led.book("soil_liquid", "atmosphere", 100.0 * MM_KG, phase_to="vapour",
             transform="vaporisation", et_component="transpiration")
    led.book("canopy", "atmosphere", 20.0 * MM_KG, phase_to="vapour",
             transform="vaporisation", et_component="transpiration")
    return 0.0


def _fixture_pme_triple_ownership(decl):
    """Positive P - E spent three times.

    Finding 3, exactly as the pipeline stands. One interval's atmospheric
    surplus over one cell is routed as ExoPlaSim's bucket overflow, taken again
    by the surface-water builder as catchment runoff, and taken a third time by
    the groundwater builder as recharge. The store holds one surplus, so the
    second withdrawal empties it and the third is refused.
    """
    led = Ledger(decl)
    surplus = 200.0 * MM_KG
    led.seed("soil_liquid", surplus)
    led.book("soil_liquid", "open_water", surplus)       # drunoff, routed
    led.book("soil_liquid", "open_water", surplus)       # max(P - E, 0), routed again
    led.book("soil_liquid", "groundwater", surplus)      # the same quantity as recharge
    return 0.0


def _fixture_area_basis_mismatch(decl):
    """A depth generated on one denominator and summed on another.

    Lake evaporation is quoted per square metre of OPEN WATER and the cell is
    mostly land; a consumer that multiplies it by the land-cell area gets a
    different mass and nothing in the units says so. This is the same shape as
    the cell-mean `dwmax` of finding 7, where a capacity built on the land
    fraction and a capacity built on the lake fraction are averaged into one
    number.
    """
    led = Ledger(decl)
    open_water_area_m2 = 0.02 * CELL_AREA_M2
    depth_m = 1.2                                  # over the open water it evaporates from
    led.seed("open_water", 5.0e10)
    led.handoff("open_water", "atmosphere",
                kg_out=depth_m * open_water_area_m2 * 1000.0,
                kg_in=depth_m * CELL_AREA_M2 * 1000.0,
                phase_to="vapour", transform="vaporisation",
                et_component="open_water_evaporation")
    return abs(led.residual(5.0e10))


def _fixture_interval_mismatch(decl):
    """A rate integrated over two different durations.

    ExoPlaSim reports `drunoff` in m/s at the model timestep; the hydrography
    path works in depth per Earth year. A consumer that integrates the
    instantaneous rate over the wrong interval changes the mass, and the
    ledger books the debit and the credit as stated so the difference appears.
    """
    led = Ledger(decl)
    rate_m_s = 2.0e-8
    led.seed("soil_liquid", 1.0e11)
    led.handoff("soil_liquid", "open_water",
                kg_out=rate_m_s * EARTH_YEAR_S * CELL_AREA_M2 * 1000.0,
                kg_in=rate_m_s * 360.0 * 86400.0 * CELL_AREA_M2 * 1000.0)
    return abs(led.residual(1.0e11))


def _fixture_unnamed_phase_change(decl):
    """Melt booked as a move.

    Snow leaves the pack as ice and arrives in the soil as liquid. Booked with
    no process named, water closes and the latent heat of fusion is neither
    absorbed nor released, which is how two models can melt the same snow on
    two different energy budgets.
    """
    led = Ledger(decl)
    led.seed("snowpack", 50.0 * MM_KG)
    led.book("snowpack", "soil_liquid", 10.0 * MM_KG,
             phase_from="ice", phase_to="liquid")
    return 0.0


def _fixture_wrong_latent_heat(decl):
    """Sublimation booked as evaporation off ice.

    The phases are right and the transform is not: a term taking ice straight
    to vapour with `vaporisation` conserves water and underpays the atmosphere
    by the latent heat of fusion on every kilogram.
    """
    led = Ledger(decl)
    led.seed("snowpack", 50.0 * MM_KG)
    led.book("snowpack", "atmosphere", 5.0 * MM_KG, phase_from="ice",
             phase_to="vapour", transform="vaporisation")
    return 0.0


def _fixture_terminal_source(decl):
    """Water drawn back out of the ocean."""
    led = Ledger(decl)
    led.book("ocean_export", "open_water", 1.0 * MM_KG)
    return 0.0


def _fixture_store_below_zero(decl):
    """A withdrawal larger than the store holds.

    This is the shape `landmod.f90`'s floor at zero would have if it were
    reachable: the net flux takes the bucket past zero and the model raises it
    back, crediting water while the latent heat that removed it has already
    been paid. `bucket_floor_bound` shows the climate column cannot get there,
    which is why the crossing is a declared absence rather than a term. The
    fixture stays because the shape is general -- any coupling that draws on a
    store harder than it holds has it -- and refusing it is what the ledger
    does instead of clipping.
    """
    led = Ledger(decl)
    led.seed("soil_liquid", 4.0 * MM_KG)
    led.book("soil_liquid", "atmosphere", 7.0 * MM_KG, phase_to="vapour",
             transform="vaporisation", et_component="soil_evaporation")
    return 0.0


def _fixture_negative_transfer(decl):
    """A transfer with a negative mass, standing in for a reversed sign.

    ExoPlaSim's `devap` is negative upward and hydrography's evaporation is
    positive, so a sign convention crossing this boundary is a live hazard.
    """
    led = Ledger(decl)
    led.seed("soil_liquid", 10.0 * MM_KG)
    led.book("soil_liquid", "open_water", -3.0 * MM_KG)
    return 0.0


# name, what it is built to do, tolerance on the residual when it should close
FIXTURES = (
    ("balanced", "closes", 1e-3, _fixture_balanced),
    ("shadow_store_evaporation", "rejects", None, _fixture_shadow_store_evaporation),
    ("double_debited_component", "rejects", None, _fixture_double_debited_component),
    ("pme_triple_ownership", "rejects", None, _fixture_pme_triple_ownership),
    ("area_basis_mismatch", "rejects", None, _fixture_area_basis_mismatch),
    ("interval_mismatch", "rejects", None, _fixture_interval_mismatch),
    ("unnamed_phase_change", "rejects", None, _fixture_unnamed_phase_change),
    ("wrong_latent_heat", "rejects", None, _fixture_wrong_latent_heat),
    ("terminal_source", "rejects", None, _fixture_terminal_source),
    ("store_below_zero", "rejects", None, _fixture_store_below_zero),
    ("negative_transfer", "rejects", None, _fixture_negative_transfer),
)


def run_fixtures(decl: dict) -> tuple[list[dict], list[str]]:
    """Every fixture, and the defects found in this checker by running them.

    A `rejects` fixture that produces no residual and raises nothing has been
    accepted by a ledger that should have refused it, which means the ledger is
    broken and not the fixture.
    """
    results, defects = [], []
    for name, expect, tol, fn in FIXTURES:
        row = {"fixture": name, "expects": expect}
        try:
            residual = fn(decl)
        except LedgerError as exc:
            row["outcome"] = "refused"
            row["reason"] = str(exc)
            if expect == "closes":
                defects.append(f"fixture {name} should close and was refused: {exc}")
        else:
            row["residual_kg"] = residual
            if expect == "closes":
                row["outcome"] = "closed" if residual <= tol else "residual"
                if residual > tol:
                    defects.append(f"fixture {name} should close and left a "
                                   f"residual of {residual:.6g} kg")
            else:
                row["outcome"] = "residual" if residual > 0 else "accepted"
                if residual == 0:
                    defects.append(f"fixture {name} should be rejected and the "
                                   "ledger accepted it without a residual")
        results.append(row)
    return results, defects


# ---------------------------------------------------------------------------
# 3. The graph check, against declarations built wrong
# ---------------------------------------------------------------------------
#
# `check_graph` returning nothing on the real declaration proves nothing on its
# own: a check that cannot fail and a declaration that is correct look the same
# from here. Each mutation below breaks the declaration in one named way and
# the graph check is required to catch it. A mutation the check passes is a
# defect in the check, not in the declaration.

def _mutate_dangling_destination(decl):
    decl["terms"]["snowmelt"]["to"] = "permafrost"


def _mutate_terminal_source(decl):
    decl["terms"]["coastal_discharge"]["from"] = "ocean_export"


def _mutate_unnamed_phase_change(decl):
    decl["terms"]["snowmelt"]["transform"] = "none"
    decl["terms"]["snowmelt"]["latent"] = "none"


def _mutate_wrong_latent_heat(decl):
    decl["terms"]["snowmelt"]["latent"] = "vaporisation"


def _mutate_reversed_phase_transform(decl):
    decl["terms"]["snowmelt"]["transform"] = "freezing"
    decl["terms"]["snowmelt"]["latent"] = "fusion"


def _mutate_double_debited_component(decl):
    decl["terms"]["soil_evaporation"]["et_component"] = "transpiration"


def _mutate_accumulating_reservoir(decl):
    for name in ("open_water_evaporation", "coastal_discharge"):
        decl["terms"].pop(name)


def _mutate_unfillable_reservoir(decl):
    decl["terms"].pop("soil_freezing")


def _mutate_store_without_owner(decl):
    decl["nodes"]["soil_liquid"].pop("owner")


def _mutate_missing_area_basis(decl):
    decl["nodes"]["open_water"].pop("area_basis")


def _mutate_fusion_constant(decl):
    decl["latent_heat"]["fusion_j_per_kg"] = 3.34e5


def _mutate_boundary_receives(decl):
    decl["terms"]["surface_runoff"]["to"] = "unowned_source"


GRAPH_MUTATIONS = (
    ("dangling_destination", _mutate_dangling_destination),
    ("terminal_source", _mutate_terminal_source),
    ("unnamed_phase_change", _mutate_unnamed_phase_change),
    ("wrong_latent_heat", _mutate_wrong_latent_heat),
    ("reversed_phase_transform", _mutate_reversed_phase_transform),
    ("double_debited_component", _mutate_double_debited_component),
    ("accumulating_reservoir", _mutate_accumulating_reservoir),
    ("unfillable_reservoir", _mutate_unfillable_reservoir),
    ("store_without_owner", _mutate_store_without_owner),
    ("missing_area_basis", _mutate_missing_area_basis),
    ("fusion_constant_disagrees", _mutate_fusion_constant),
    ("boundary_receives_from_inside", _mutate_boundary_receives),
)


def run_graph_mutations(path: Path = DECLARATION) -> tuple[list[dict], list[str]]:
    """Every mutation, and the defects in `check_graph` that running them found."""
    rows, defects = [], []
    for name, mutate in GRAPH_MUTATIONS:
        broken = load(path)
        mutate(broken)
        caught = check_graph(broken)
        rows.append({"mutation": name, "failures": len(caught),
                     "caught": bool(caught),
                     "first": caught[0] if caught else None})
        if not caught:
            defects.append(f"graph mutation {name} was not caught. The check "
                           "passes a declaration built wrong in that way")
    return rows, defects


# ---------------------------------------------------------------------------
# 4. The unowned crossing, measured
# ---------------------------------------------------------------------------
#
# `bucket_floor_creation` is the one term booked from the `unowned_source`
# boundary: `landmod.f90`'s land water step adds the net surface flux to the
# store and then raises it back to zero, which credits the store water no node
# was debited for, while the latent heat that removed that water has already
# been paid to the air column.
#
# THIS MEASURES IT RATHER THAN DESIGNING FOR IT, and the answer is that the
# crossing does not exist in exact arithmetic. Three properties of the model's
# own code close it, and none of them is the snow path the declaration first
# named:
#
#   1. `fluxmod.f90`'s humidity solve takes AMAX1 against the level's own
#      humidity before anything else, so the surface evaporation is never
#      negative and the land column cannot gain water from the air here.
#   2. The line below it limits that same quantity to the store divided by the
#      timestep, on every land cell, whatever the snowpack is doing. So the
#      withdrawal charged to the store over one step is at most the store.
#   3. Every other contribution to the land surface flux adds. Snowfall is a
#      SUBSET of the total precipitation in `rainmod.f90`, so precipitation
#      minus snowfall is the rain and is non-negative; melt and the snow depth
#      cap both return water; and where the snowpack absorbs the evaporation
#      the remainder charged to the store is smaller still.
#
# The snow-exhaustion path the declaration named is the wrong suspect. The
# limit is applied unconditionally, so a snowpack can only take part of a
# withdrawal that was already capped against the soil store, and exhausting it
# leaves the store FULLER than the bare-ground case rather than emptier.
#
# What is left is the rounding of a difference that cancels exactly, and the
# sweep below bounds it. The two arms after it are what make the sweep a
# check: each removes one of the properties above and each is required to
# produce a crossing many orders larger. A sweep that cannot report a crossing
# is not evidence that there is none.

# The model's dry gas constant is derived from the declared composition at run
# time; this is the value config/planet.yaml records that derivation at. The
# bound is insensitive to it: it enters the limit and the flux at the same
# place and cancels.
DRY_GAS_CONSTANT_J_KG_K = 287.017

# `landmod.f90`'s wetness knee: the fraction of capacity above which the land
# column evaporates at the potential rate.
WETNESS_KNEE = 0.4

WATER_DENSITY_KG_M3 = 1000.0

EARTH_RADIUS_M = 6371000.0

# The swept axes. Deliberately wider than this world reaches, because the
# result is a BOUND. In order: the mass of the lowest model layer, the surface
# pressure, the surface temperature, the turbulent transfer coefficient, the
# surface saturation humidity, the level humidity as a share of it, the bucket
# depth from the driest cell of the soil map to the model's uniform default,
# and the store as a share of the withdrawal the unlimited solve would take,
# which straddles the point the limit starts binding rather than sampling
# around it.
_SWEEP_LAYER_MASS_KG_M2 = (200.0, 500.0, 780.0, 1200.0, 2000.0)
_SWEEP_SURFACE_PRESSURE_PA = (0.5e5, 1.0e5, 1.5e5)
_SWEEP_SURFACE_TEMPERATURE_K = (230.0, 280.0, 320.0, 340.0)
_SWEEP_TRANSFER_M_S = (1e-4, 1e-3, 1e-2, 5e-2, 2e-1)
_SWEEP_SURFACE_HUMIDITY = (1e-5, 1e-3, 1e-2, 5e-2, 1.5e-1)
_SWEEP_LEVEL_SHARE = (0.0, 0.1, 0.5, 0.9, 0.999)
_SWEEP_BUCKET_M = (0.010, 0.05, 0.2, 0.5)
_SWEEP_STORE_SHARE = (0.0, 1e-12, 1e-6, 0.25, 0.5, 0.9, 0.99, 0.999,
                      1.0, 1.0001, 1.01, 2.0)
# (large-scale rain, convective rain, large-scale snow, convective snow), m/s.
# Split this way because `rainmod.f90` builds the total and the snowfall rate
# out of the same four condensate sums, which is property 3 above.
_SWEEP_PRECIP_M_S = ((0.0, 0.0, 0.0, 0.0),
                     (1e-8, 0.0, 0.0, 0.0),
                     (0.0, 1e-7, 0.0, 0.0),
                     (0.0, 0.0, 1e-7, 0.0),
                     (0.0, 0.0, 0.0, 1e-6),
                     (1e-8, 0.0, 1e-7, 0.0),
                     (1e-6, 1e-6, 1e-6, 1e-6))
_SWEEP_SNOWPACK_M = (0.0, 1e-9, 1e-6, 1e-4, 1e-2, 1.0)


def _sweep(deltsec: float, ga: float, limiter: bool = True,
           snow_is_subset: bool = True) -> dict:
    """One timestep of the land water chain over every swept combination.

    `fluxmod.f90` vdiff's humidity solve and evaporation limit, then
    `landmod.f90` tands' snow partition, then `landcolumn.f90` bucket_step's
    floor, in the order the model runs them and in the model's own real kind.
    The branches that depend only on the precipitation and the snowpack are
    taken outside the array expression, so each is evaluated as the model
    evaluates it rather than as a blended `where`.

    `limiter` and `snow_is_subset` are the falsification arms and are true for
    the model as it stands.
    """
    import itertools

    import numpy as np

    grid = np.meshgrid(*(np.asarray(axis, dtype=float) for axis in (
        _SWEEP_LAYER_MASS_KG_M2, _SWEEP_SURFACE_PRESSURE_PA,
        _SWEEP_SURFACE_TEMPERATURE_K, _SWEEP_TRANSFER_M_S,
        _SWEEP_SURFACE_HUMIDITY, _SWEEP_LEVEL_SHARE, _SWEEP_BUCKET_M,
        _SWEEP_STORE_SHARE)), indexing="ij")
    layer_mass, dp, ts, dtransh, dq_lep, share, dwmax, store_share = (
        axis.ravel() for axis in grid)
    dsigma = layer_mass * ga / dp
    dq_lev = share * dq_lep
    keep = dsigma <= 0.5
    layer_mass, dp, ts, dtransh = (a[keep] for a in (layer_mass, dp, ts, dtransh))
    dq_lep, dwmax, store_share = (a[keep] for a in (dq_lep, dwmax, store_share))
    dsigma, dq_lev = dsigma[keep], dq_lev[keep]

    worst, worst_index, worst_key = 0.0, None, None
    cases = bound_cases = exhausted_cases = 0
    for deltsec2, precip, dsnowz in itertools.product(
            (deltsec, 2.0 * deltsec), _SWEEP_PRECIP_M_S, _SWEEP_SNOWPACK_M):
        rain_l, rain_c, snow_l, snow_c = precip
        dprl, dprc = rain_l + snow_l, rain_c + snow_c
        dprs = snow_l + snow_c
        if not snow_is_subset:
            dprs = dprs + dprl + dprc

        zkonst1 = ga * deltsec2 / (DRY_GAS_CONSTANT_J_KG_K * dsigma)
        zkonst2 = dsigma / deltsec2 / ga
        # The store, placed against the withdrawal the unlimited solve takes.
        zk_full = zkonst1 * dtransh / ts
        zqn_full = np.maximum(dq_lev,
                              (dq_lev + zk_full * dq_lep) / (1.0 + zk_full))
        full = dp * zkonst2 / 1000.0 * (zqn_full - dq_lev) * deltsec
        dwatc = store_share * full
        live = (dwatc >= 0.0) & (dwatc <= dwmax)

        drhs = np.where(dwmax <= 0.0, 1.0,
                        np.minimum(1.0, dwatc / (WETNESS_KNEE * dwmax)))
        zkdiff = drhs * zkonst1 * dtransh / ts
        zqn = np.maximum(dq_lev, (dq_lev + zkdiff * dq_lep) / (1.0 + zkdiff))
        unlimited = zqn
        if limiter:
            zqn = np.minimum(zqn, dwatc / deltsec * 1000.0 / dp / zkonst2 + dq_lev)
        devap = -dp * zkonst2 / 1000.0 * (zqn - dq_lev)

        zdsnowz = np.full_like(devap, dprs if dprs > 0.0 else 0.0)
        if dsnowz > 0.0:
            zdsnowz = zdsnowz + devap
        zsnowz = np.maximum(0.0, dsnowz + zdsnowz * deltsec)
        zdsnowz = (zsnowz - dsnowz) / deltsec
        dwater = devap + dprl + dprc - zdsnowz
        zw = dwatc + deltsec * dwater
        created = np.where(live, np.maximum(0.0, -zw), 0.0)

        cases += int(live.sum())
        bound_cases += int((live & (zqn < unlimited)).sum())
        exhausted_cases += int((live & (zsnowz == 0.0)
                                & (dsnowz > 0.0 or dprs > 0.0)).sum())
        here = int(created.argmax())
        if created[here] > worst:
            worst = float(created[here])
            worst_index = here
            worst_key = (deltsec2, precip, dsnowz)

    case = None
    if worst_index is not None:
        deltsec2, precip, dsnowz = worst_key
        i = worst_index
        case = {"lowest_layer_mass_kg_m2": float(layer_mass[i]),
                "surface_pressure_pa": float(dp[i]),
                "surface_temperature_k": float(ts[i]),
                "transfer_coefficient_m_s": float(dtransh[i]),
                "surface_humidity": float(dq_lep[i]),
                "level_humidity": float(dq_lev[i]),
                "bucket_depth_m": float(dwmax[i]),
                "leapfrog_seconds": float(deltsec2),
                "precipitation_m_s": list(precip),
                "snowpack_m": float(dsnowz)}
    return {"cases": cases, "limiter_bound": bound_cases,
            "snowpack_exhausted": exhausted_cases,
            "worst_created_m": worst, "worst_case": case}


def bucket_floor_bound() -> dict:
    """What `bucket_floor_creation` is worth, in kilograms.

    The depth bound comes from the sweep and is a property of the model's
    arithmetic alone. Turning it into a mass needs an area and an interval, and
    the ledger's own rule is that the area is an argument: the land cells and
    their areas come from the soil map of the configured build on the
    configured grid, so the denominator is the land cell the crossing is stated
    on rather than a reconstruction of one.
    """
    import numpy as np
    import yaml as _yaml
    from numpy.polynomial.legendre import leggauss

    import builds
    import orbit

    config = _yaml.safe_load((PROJECT_ROOT / "config" / "planet.yaml")
                             .read_text(encoding="utf-8"))
    ga = float(config["planet"]["gravity_m_s2"])
    radius_m = float(config["planet"]["radius_earth"]) * EARTH_RADIUS_M
    model = config["model"]
    nlat, nlon = int(model["latitudes"]), int(model["longitudes"])
    deltsec = float(model["timestep_minutes"]) * 60.0
    steps_per_orbit = int(round(orbit.orbital_year_days(config)
                                * float(config["planet"]["rotation_hours"])
                                * 3600.0 / deltsec))

    measured = _sweep(deltsec, ga)

    # The falsification arms. Each removes ONE of the properties the bound
    # rests on and is required to produce a crossing far above the rounding
    # residue. The bar is one micrometre of water in one step: twelve orders
    # above the residue, and still far below anything the ledger would report
    # as a term. An arm that stays at the residue means the sweep is measuring
    # nothing, which is a defect in the check and not a result about the model.
    macroscopic_m = 1e-6
    arms = []
    for name, kwargs, why in (
        ("evaporation_limit_removed", {"limiter": False},
         "fluxmod's cap of the withdrawal at the store divided by the "
         "timestep is what empties the store exactly rather than past empty"),
        ("snowfall_not_a_subset_of_precipitation", {"snow_is_subset": False},
         "rainmod builds the snowfall rate out of the same condensate sums "
         "the total is built from, so precipitation minus snowfall is the "
         "rain and cannot be negative"),
    ):
        arm = _sweep(deltsec, ga, **kwargs)
        arms.append({"arm": name, "why_it_matters": why,
                     "worst_created_m": arm["worst_created_m"],
                     "reached_macroscopic": arm["worst_created_m"] >= macroscopic_m})

    weights = leggauss(nlat)[1][::-1]
    latitudes = np.rad2deg(np.arcsin(leggauss(nlat)[0][::-1]))
    area_by_latitude = weights * 2.0 * np.pi * radius_m * radius_m / nlon

    # The capacity the MODEL installs, which is the land column property
    # contract's, not the soil map's `awc`. Those are two numbers for one soil
    # and WORLD-OF6N settled which one the world has; a stock computed from the
    # other would not be the stock the bucket floor is being weighed against.
    states = builds.land_column_states(config)
    lines = states.read_text(encoding="utf-8").splitlines()
    header = lines[0].split()
    lat_column, awc_column = header.index("Lat"), header.index("awc_mm")
    rows = [line.split() for line in lines[1:] if line.strip()]
    cell_lat = np.array([float(row[lat_column]) for row in rows])
    awc_mm = np.array([float(row[awc_column]) for row in rows])
    nearest = np.abs(cell_lat[:, None] - latitudes[None, :]).argmin(axis=1)
    cell_area = area_by_latitude[nearest]

    worst_m = measured["worst_created_m"]
    per_cell_step_kg = worst_m * WATER_DENSITY_KG_M3 * float(cell_area.max())
    ceiling_kg = (worst_m * WATER_DENSITY_KG_M3 * float(cell_area.sum())
                  * steps_per_orbit)
    stock_kg = float((awc_mm / 1000.0 * cell_area).sum()) * WATER_DENSITY_KG_M3

    defects = [f"bucket floor arm {arm['arm']} did not reach a macroscopic "
               "crossing, so the sweep cannot report one and its bound on the "
               "unmutated model is not evidence"
               for arm in arms if not arm["reached_macroscopic"]]

    return {
        "term": "bucket_floor_creation",
        "verdict": (
            "no crossing in exact arithmetic. The land water step's floor at "
            "zero is unreachable: the surface evaporation is non-negative, it "
            "is capped at the store divided by the timestep on every land "
            "cell whatever the snowpack holds, and every other contribution "
            "to the surface flux adds. What is left is the rounding of a "
            "difference that cancels exactly"),
        "snow_exhaustion_path": (
            "not the mechanism. The cap is applied unconditionally, so a "
            "snowpack absorbs part of a withdrawal already limited against "
            "the soil store and leaves the store fuller than bare ground "
            "would, not emptier"),
        "rests_on": [
            "fluxmod.f90 vdiff: AMAX1 against the level humidity, so the "
            "surface evaporation never adds water to the land surface flux",
            "fluxmod.f90 vdiff: the withdrawal is capped at dwatc/deltsec on "
            "every land cell",
            "rainmod.f90: the snowfall rate is a subset of the total "
            "precipitation rate",
            "landmod.f90 soilini and the floor itself: the store enters every "
            "step at or above zero, which restart_schema.py carries as the "
            "lower bound on dwatc and dwatcl",
        ],
        "swept_cases": measured["cases"],
        "cases_with_the_limit_binding": measured["limiter_bound"],
        "cases_with_the_snowpack_exhausted": measured["snowpack_exhausted"],
        "worst_created_m_per_cell_step": worst_m,
        "worst_created_kg_per_cell_step": per_cell_step_kg,
        "worst_case": measured["worst_case"],
        "orbit_ceiling_kg": ceiling_kg,
        "orbit_ceiling_basis": (
            "every land cell of the configured build, at every timestep of "
            "one orbit, at the worst crossing the sweep found. Nothing can "
            "reach it: the floor is touched at all only where the cap binds"),
        "land_soil_water_stock_kg": stock_kg,
        "orbit_ceiling_over_stock": ceiling_kg / stock_kg if stock_kg else None,
        "land_cells": len(rows),
        "steps_per_orbit": steps_per_orbit,
        "land_column_states": str(states.relative_to(PROJECT_ROOT)),
        "falsification_arms": arms,
        "checker_defects": defects,
    }


# ---------------------------------------------------------------------------

def build_report(decl: dict) -> dict:
    graph = check_graph(decl)
    fixtures, defects = run_fixtures(decl)
    mutations, mutation_defects = run_graph_mutations()
    floor = bucket_floor_bound()
    defects = defects + mutation_defects + floor["checker_defects"]
    return {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": "hydrography/config/land_water_ledger.yaml",
        "declaration_sha256": hashlib.sha256(DECLARATION.read_bytes()).hexdigest(),
        "version": decl["version"],
        "units": decl["units"],
        "control_volume": decl["control_volume"],
        "interval_earth_days": decl["interval_earth_days"],
        "latent_heat": decl["latent_heat"],
        "graph_failures": graph,
        "graph_mutations": mutations,
        "closure_fixtures": fixtures,
        "checker_defects": defects,
        "shadow_copies": shadow_copies(decl),
        "interval_refusals": interval_refusals(decl),
        "undeclared_terms": [{"term": t, "owner": o} for t, o in undeclared_terms(decl)],
        "absences": [{"absence": name, "owner": body.get("owner"),
                      "what": " ".join(body["what"].split())}
                     for name, body in decl["absences"].items()],
        "bucket_floor_bound": floor,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 while the ledger does not close: any term "
                             "undeclared, any store with a shadow copy, or any "
                             "term whose producer cannot resolve its destination")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    decl = load()
    report = build_report(decl)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2, sort_keys=False) + "\n",
                      encoding="utf-8")

    hard = report["graph_failures"] + report["checker_defects"]
    open_items = (report["undeclared_terms"] + report["shadow_copies"]
                  + report["interval_refusals"])
    if not args.quiet:
        print(f"land water ledger v{report['version']}: {len(decl['nodes'])} nodes, "
              f"{len(decl['terms'])} terms, {len(decl['absences'])} named absences")
        for row in report["closure_fixtures"]:
            print(f"  fixture {row['fixture']:<26} expects {row['expects']:<8} "
                  f"-> {row['outcome']}")
        caught = sum(1 for row in report["graph_mutations"] if row["caught"])
        print(f"  graph mutations            {caught} of "
              f"{len(report['graph_mutations'])} declarations built wrong were caught")
        for row in report["shadow_copies"]:
            print(f"  SHADOW  {row['node']:<14} owned by {row['owner']} as "
                  f"{row['symbol']}, also held as {', '.join(row['shadow_copies'])}")
        for line in report["interval_refusals"]:
            print(f"  TIMING  {line}")
        floor = report["bucket_floor_bound"]
        print(f"  BOUND   {floor['term']}: "
              f"{floor['worst_created_kg_per_cell_step']:.3g} kg per land cell "
              f"and timestep at worst over {floor['swept_cases']} swept cases, "
              f"{floor['orbit_ceiling_kg']:.3g} kg over one orbit if every "
              f"land cell clipped at every step, "
              f"{floor['orbit_ceiling_over_stock']:.3g} of the land soil "
              "water stock")
        print(f"  {len(report['undeclared_terms'])} of {len(decl['terms'])} terms "
              "undeclared; the ledger is defined and does not close")
        for problem in hard:
            print(f"  FAIL {problem}")
        print(f"wrote {REPORT.relative_to(PROJECT_ROOT)}")

    if hard:
        return 1
    if args.strict and open_items:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
