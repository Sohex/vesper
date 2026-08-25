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
   domains where the answer is known in advance. Eight of the nine fixtures are
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

    This is the shape of `wandr`'s clip: the net flux takes the bucket past
    zero and the model raises it back, creating water while the latent heat
    that removed it has already been paid. The ledger refuses instead, which is
    what makes the clip a term with an owner rather than a rounding.
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

def build_report(decl: dict) -> dict:
    graph = check_graph(decl)
    fixtures, defects = run_fixtures(decl)
    mutations, mutation_defects = run_graph_mutations()
    defects = defects + mutation_defects
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
