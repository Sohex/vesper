#!/usr/bin/env python3
"""The abiotic nutrient ledger: what balances, over what control volume, and the
check that fails when it does not.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a declared graph of control volumes, the mass transfers
between them, and the screens that decide which sources are allowed to book
into it.

    python biosphere/scripts/abiotic_nutrient_ledger.py            # status, exit 0
    python biosphere/scripts/abiotic_nutrient_ledger.py --strict   # exit 1 while any term is undeclared

This module is the enforcement for `biosphere/config/abiotic_nutrients.yaml`.
It does four things, in this order, because each depends on the one before it:

1. **The graph.** Every term moves one element from one declared node to one
   declared node. A dangling destination, a term drawing on a terminal node, a
   reservoir with no outflow, a transfer that changes element, or a refractory
   particle reaching the root zone without a dissolution term is an error.

2. **Closure, on fixtures that can fail.** A ledger that only ever runs on real
   inputs has no failing case, so the balance is exercised on reduced synthetic
   domains where the answer is known in advance. Six of the seven fixtures are
   built to be WRONG in a named way and are required to be rejected; a fixture
   that does not get the verdict it was built for is a defect in this checker
   and exits non-zero even when nothing else is wrong.

3. **The ANUT-7 screen.** Every registered candidate has a verdict, a rule that
   produced it, and a rate provenance. A candidate carried on an aerosol optical
   depth is refused by name whatever its verdict: an optical depth folds in
   refractive index, size distribution and water uptake, and no elemental mass
   comes back out of it.

4. **The ANUT-8 adequacy screen.** A critical runoff per element and per
   lithology, from bounds whose directions are declared in the config before any
   result is seen. It can return `adequate`; it can return `not_settled`, which
   is a bracket and not a verdict of inadequate; and `not_adequate` has a stated
   meaning it could reach.

The report goes to `biosphere/generated/abiotic_nutrient_ledger_report.json`.
The contract is `biosphere/notes/abiotic-nutrient-ledger.md`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import GENERATED, PROJECT_ROOT  # noqa: F401  (adds lib/ to sys.path)

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = COMPONENT_ROOT / "config" / "abiotic_nutrients.yaml"
REPORT = GENERATED / "abiotic_nutrient_ledger_report.json"

UNDECLARED = "undeclared"

ELEMENTS = ("N", "P", "K", "Ca", "Mg", "S", "Fe", "Na", "Cl", "Si")

# Phases a term may carry. The distinction that matters is the last two: a
# refractory particle is not a nutrient until something dissolves it, and the
# CNP fork's `snow_pinput`/`pmass_add` path would make every deposited apatite
# grain immediately available if a deposition flux were booked straight into
# the labile system.
PHASES = ("dissolved_inorganic", "dissolved_organic",
          "particulate_reactive", "particulate_refractory", "gas")

# The transform names a term may declare. `none` moves material unchanged;
# everything else must change the phase, and none of them may change the
# element.
TRANSFORMS = ("none", "dissolution", "sorption", "desorption", "occlusion")

# Equivalent weights, g per equivalent, for the Meybeck table's ionic columns.
# Sulfur is taken from the sulfate column at two equivalents per sulfur, which
# is exactly the conversion the ledger's element basis exists to make explicit.
EQUIVALENT_WEIGHT = {"Ca": 40.078 / 2, "Mg": 24.305 / 2, "K": 39.0983,
                     "Na": 22.98977, "Cl": 35.453, "S": 32.06 / 2}
MEYBECK_COLUMN = {"Ca": "ca_ueq_l", "Mg": "mg_ueq_l", "K": "k_ueq_l",
                  "Na": "na_ueq_l", "Cl": "cl_ueq_l", "S": "so4_ueq_l"}


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
                       "is not declared reentrant. Mass leaving the domain has a "
                       "terminal node to go to")

        elements = term.get("elements") or []
        if not elements:
            bad.append(f"term {name}: carries no element")
        for element in elements:
            if element not in ELEMENTS:
                bad.append(f"term {name}: unknown element {element!r}")

        pf, pt = term.get("phase_from"), term.get("phase_to")
        for phase, where in ((pf, "phase_from"), (pt, "phase_to")):
            if phase not in PHASES:
                bad.append(f"term {name}: {where} names {phase!r}, which is not a declared phase")
        transform = term.get("transform")
        if transform not in TRANSFORMS:
            bad.append(f"term {name}: unknown transform {transform!r}")
        elif transform == "none" and pf != pt:
            bad.append(f"term {name}: transform is none but the phase changes "
                       f"{pf} to {pt}. A phase change is a process and needs naming")
        elif transform != "none" and pf == pt:
            bad.append(f"term {name}: transform is {transform} and the phase does not change")

        if (pt == "dissolved_inorganic" and pf in ("particulate_refractory",
                                                   "particulate_reactive")
                and transform not in ("dissolution", "desorption")):
            bad.append(f"term {name}: puts {pf} material into a dissolved pool "
                       "without a dissolution or desorption term")
        if dst == "soil_solution" and pt != "dissolved_inorganic":
            bad.append(f"term {name}: delivers {pt} to soil_solution. "
                       "A root takes from solution and nothing else")

    for name, node in nodes.items():
        if node["kind"] == "reservoir":
            if outflow[name] == 0:
                bad.append(f"node {name}: a reservoir with no outflow is an accumulator, "
                           "and endorheic retention is exactly the error that shape hides")
            if inflow[name] == 0:
                bad.append(f"node {name}: a reservoir with no inflow cannot be filled")
        if node["kind"] == "boundary" and inflow[name] and not node.get("reentrant"):
            bad.append(f"node {name}: a boundary that is not declared reentrant "
                       "receives nothing from inside the domain")
        if node.get("area_basis") is None:
            bad.append(f"node {name}: no area_basis. A per-area rate has no meaning without one")
    return bad


def undeclared_terms(decl: dict) -> list[tuple[str, str]]:
    """(term, owning issue) for every term still carrying the sentinel."""
    return [(name, term.get("owner", "unowned"))
            for name, term in decl["terms"].items()
            if term.get("flux") == UNDECLARED]


# ---------------------------------------------------------------------------
# 2. Closure, on fixtures that can fail
# ---------------------------------------------------------------------------

class LedgerError(RuntimeError):
    """A transfer the ledger refuses to book."""


class Ledger:
    """A reduced ledger over a declared node set, in kg of one element.

    Everything is booked in ABSOLUTE mass. A per-area rate is converted by the
    area of the node it is stated on, and that area is an argument rather than
    an assumption, so quoting a rate on the source area and summing it on the
    receiving area produces a residual instead of a silent rescale.
    """

    def __init__(self, decl: dict, element: str):
        self.kinds = {name: node["kind"] for name, node in decl["nodes"].items()}
        self.element = element
        self.stock = {name: 0.0 for name in self.kinds}
        self.external_in = 0.0
        self.terminal_out = 0.0

    def seed(self, node: str, kg: float) -> None:
        """Put mass in a reservoir before the interval starts."""
        if self.kinds.get(node) != "reservoir":
            raise LedgerError(f"cannot seed {node}: only a reservoir holds an initial stock")
        self.stock[node] += kg

    def book(self, frm: str, to: str, kg: float, *,
             element: str | None = None) -> None:
        if element is not None and element != self.element:
            raise LedgerError(
                f"transfer carries {element} into a {self.element} ledger. "
                "A term never changes element")
        if frm not in self.kinds or to not in self.kinds:
            raise LedgerError(f"transfer {frm} -> {to} names a node the ledger does not have")
        if self.kinds[frm] == "terminal":
            raise LedgerError(f"transfer draws on {frm}, which is terminal")
        if kg < 0:
            raise LedgerError("a transfer carries a non-negative mass; reverse the direction instead")

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

    def handoff(self, frm: str, to: str, kg_out: float, kg_in: float) -> None:
        """A transfer whose two sides are stated SEPARATELY.

        This is what a cross-component transfer actually looks like: one script
        writes a per-area rate and another reads it and multiplies by an area of
        its own. When the two denominators differ the mass changes and nothing
        in either file says so, so the ledger books the debit and the credit as
        given and lets the difference appear as a residual.
        """
        self.book(frm, to, kg_out)
        if self.kinds[to] == "terminal":
            self.terminal_out += kg_in - kg_out
        else:
            self.stock[to] += kg_in - kg_out

    def residual(self, seeded_kg: float) -> float:
        held = sum(v for name, v in self.stock.items()
                   if self.kinds[name] != "terminal")
        return (held + self.terminal_out) - (seeded_kg + self.external_in)


def _fixture_balanced(decl):
    """Import, weather, take up, leach, export. Closes exactly."""
    led = Ledger(decl, "Ca")
    led.seed("parent_material", 1000.0)
    led.book("parent_material", "soil_solution", 40.0)
    led.book("atmosphere_import", "atmosphere_column", 5.0)
    led.book("atmosphere_column", "soil_solution", 5.0)
    led.book("soil_solution", "plant_and_litter", 30.0)
    led.book("plant_and_litter", "soil_solution", 28.0)
    led.book("soil_solution", "surface_water", 12.0)
    led.book("surface_water", "coastal_export", 12.0)
    return abs(led.residual(1000.0))


def _fixture_area_basis(decl):
    """The same delivery, quoted on the source area and summed on the receiving
    area. This is the error `phosphorus_budget.py` makes geometrically when it
    concentrates a relative score over a basin floor: nothing in the units says
    which denominator is meant, and the mass silently changes."""
    led = Ledger(decl, "P")
    source_area, receiving_area = 1.0e9, 2.0e7
    rate = 3.0e-6                      # kg/m2 over the SOURCE area
    led.seed("regolith_mineral", 1.0e6)
    led.handoff("regolith_mineral", "surface_water",
                kg_out=rate * source_area, kg_in=rate * receiving_area)
    return abs(led.residual(1.0e6))


def _fixture_dangling_destination(decl):
    led = Ledger(decl, "P")
    led.seed("soil_solution", 10.0)
    led.book("soil_solution", "ocean_productivity", 1.0)
    return 0.0


def _fixture_terminal_source(decl):
    led = Ledger(decl, "P")
    led.book("coastal_export", "surface_water", 1.0)
    return 0.0


def _fixture_element_change(decl):
    led = Ledger(decl, "P")
    led.seed("soil_solution", 10.0)
    led.book("soil_solution", "plant_and_litter", 1.0, element="N")
    return 0.0


def _fixture_stock_below_zero(decl):
    """A source drawn harder than it holds. This is what a weathering flux does
    when it is applied to a parent stock nothing renews, and it is why
    `geomorphic_renewal` has to be a term rather than an assumption."""
    led = Ledger(decl, "P")
    led.seed("parent_material", 5.0)
    led.book("parent_material", "soil_solution", 7.0)
    return 0.0


def _fixture_negative_transfer(decl):
    led = Ledger(decl, "P")
    led.seed("soil_solution", 10.0)
    led.book("soil_solution", "surface_water", -3.0)
    return 0.0


# name, what it is built to do, tolerance on the residual when it should close
FIXTURES = (
    ("balanced", "closes", 1e-9, _fixture_balanced),
    ("area_basis_mismatch", "rejects", None, _fixture_area_basis),
    ("dangling_destination", "rejects", None, _fixture_dangling_destination),
    ("terminal_source", "rejects", None, _fixture_terminal_source),
    ("element_change", "rejects", None, _fixture_element_change),
    ("stock_below_zero", "rejects", None, _fixture_stock_below_zero),
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
# 3. The ANUT-7 screen
# ---------------------------------------------------------------------------

def check_screen(decl: dict) -> tuple[list[dict], list[str]]:
    screen = decl["screen"]
    forbidden = tuple(str(c).lower() for c in screen["forbidden_carriers"])
    terms = decl["terms"]
    rows, refusals = [], []
    for name, cand in screen["candidates"].items():
        verdict = cand.get("verdict")
        row = {"candidate": name, "verdict": verdict,
               "rule": cand.get("verdict_rule"),
               "rate_source": cand.get("rate_source"),
               "pulse_timing": cand.get("pulse_timing"),
               "books_into": cand.get("books_into")}
        if verdict not in ("retain", "register_only", "out_of_scope"):
            refusals.append(f"screen {name}: verdict {verdict!r} is not one of "
                            "retain, register_only, out_of_scope")
        if verdict in ("retain", "register_only"):
            booked = cand.get("books_into") or []
            if isinstance(booked, str):
                booked = [booked]
            if not booked:
                refusals.append(f"screen {name}: {verdict} and books into nothing. "
                                "A retained source with nowhere to book is a list, not a register")
            for target in booked:
                if target not in terms:
                    refusals.append(f"screen {name}: {verdict} and books into "
                                    f"{target!r}, which is not a ledger term. "
                                    "A retained source with nowhere to book is a list, not a register")
            if not cand.get("rate_source"):
                refusals.append(f"screen {name}: {verdict} with no rate_source. "
                                "Orogen has no time axis, so where a rate comes from is "
                                "part of the verdict")
            if not cand.get("pulse_timing"):
                refusals.append(f"screen {name}: {verdict} with no pulse_timing verdict")
        carrier = str(cand.get("carrier", "")).lower()
        if carrier and carrier != UNDECLARED:
            for token in forbidden:
                if token in carrier:
                    refusals.append(
                        f"screen {name}: carried on {cand['carrier']!r}. An optical "
                        "depth is a mass-and-optics quantity and no elemental mass "
                        "comes back out of it")
                    row["carrier_refused"] = True
        rows.append(row)
    return rows, refusals


# ---------------------------------------------------------------------------
# 4. The ANUT-8 adequacy screen
# ---------------------------------------------------------------------------

def standing_pool_g_per_m2(decl: dict) -> dict[str, float]:
    """Upper-bound standing circulating pool, g of element per m2 of land.

    Direction: UP, at every step. The maximum over a compilation rather than a
    mean, the maximum sulfur content rather than a central one, and the upper
    end of the below-ground and exchangeable bracket. A larger pool takes longer
    to build, so every choice here makes the screen harder to pass.
    """
    cfg = decl["adequacy"]
    multiplier = float(cfg["belowground_and_exchangeable_multiplier"])
    pool = {element: kg_ha / 10.0 * multiplier          # kg/ha -> g/m2
            for element, kg_ha in cfg["standing_pool_aboveground_kg_ha"].items()}
    s = cfg["sulfur_from_biomass"]
    dry_g_m2 = float(s["max_biomass_t_ha"]) * 100.0     # t/ha -> g/m2
    carbon_g_m2 = dry_g_m2 * float(s["carbon_fraction_of_dry_mass"])
    pool["S"] = carbon_g_m2 / float(s["min_c_to_s_mass"]) * multiplier
    return pool


def release_mg_per_litre(table: dict, element: str) -> dict[str, float]:
    """Per-lithology dissolved release concentration, mg of element per litre.

    Meybeck (1987) Table 2C is the representative stream analysis for each rock
    type, so concentration times runoff is the areal release flux -- which is
    the paper's own model form and not a construction added here.
    """
    columns = table["table_2c"]["columns"]
    index = columns.index(MEYBECK_COLUMN[element])
    weight = EQUIVALENT_WEIGHT[element]
    out = {}
    for rock, row in table["table_2c"]["rows"].items():
        ueq = row[index]
        if ueq is None:
            continue
        out[rock] = ueq * weight * 1e-3
    return out


def adequacy_screen(decl: dict) -> dict:
    cfg = decl["adequacy"]
    table = json.loads((PROJECT_ROOT / cfg["release_table"]).read_text(encoding="utf-8"))
    residence = float(cfg["residence_time_earth_years"])
    pool = standing_pool_g_per_m2(decl)

    elements = {}
    for element, spec in cfg["elements"].items():
        if not spec.get("screened"):
            elements[element] = {"screened": False, "refusal": spec["refusal"].strip(),
                                 "elements": spec.get("elements")}
            continue
        concentrations = release_mg_per_litre(table, element)
        q = pool[element]
        # F(g/m2/yr) = c(mg/L) * R(m/yr), so the runoff at which weathering
        # release alone builds the whole standing pool in one soil residence
        # time is Q / (c * T). Below it the screen is not settled; above it the
        # element is adequate on rock weathering alone, with every other source
        # counted as zero.
        critical = {rock: q / (c * residence) * 1000.0     # m/yr -> mm/yr
                    for rock, c in concentrations.items()}
        ordered = sorted(critical.items(), key=lambda kv: kv[1])
        elements[element] = {
            "screened": True,
            "standing_pool_g_per_m2": q,
            "critical_runoff_mm_per_earth_year": critical,
            "easiest_lithology": ordered[0][0],
            "hardest_lithology": ordered[-1][0],
            "verdict": "adequate where local runoff exceeds the critical value "
                       "for the cell's lithology; not_settled below it",
        }
    return {
        "residence_time_earth_years": residence,
        "residence_time_bracket_earth_years": cfg["residence_time_bracket_earth_years"],
        "belowground_and_exchangeable_bracket": cfg["belowground_and_exchangeable_bracket"],
        "release_table": cfg["release_table"],
        "release_table_sha256": table["source"]["sha256"],
        "model_boundary_declaration": cfg["model_boundary_declaration"].strip(),
        "elements": elements,
    }


# ---------------------------------------------------------------------------

def build_report(decl: dict) -> dict:
    graph = check_graph(decl)
    fixtures, defects = run_fixtures(decl)
    screen, refusals = check_screen(decl)
    return {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": "biosphere/config/abiotic_nutrients.yaml",
        "declaration_sha256": hashlib.sha256(
            DECLARATION.read_bytes()).hexdigest(),
        "version": decl["version"],
        "units": decl["units"],
        "root_zone_depth_m": decl["root_zone_depth_m"],
        "graph_failures": graph,
        "closure_fixtures": fixtures,
        "checker_defects": defects,
        "screen": screen,
        "screen_refusals": refusals,
        "undeclared_terms": [{"term": t, "owner": o} for t, o in undeclared_terms(decl)],
        "adequacy": adequacy_screen(decl),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--strict", action="store_true",
                        help="exit 1 while any ledger term is still undeclared")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    decl = load()
    report = build_report(decl)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2, sort_keys=False) + "\n",
                      encoding="utf-8")

    hard = (report["graph_failures"] + report["checker_defects"]
            + report["screen_refusals"])
    if not args.quiet:
        print(f"ledger v{report['version']}: {len(decl['nodes'])} nodes, "
              f"{len(decl['terms'])} terms")
        for row in report["closure_fixtures"]:
            print(f"  fixture {row['fixture']:<22} expects {row['expects']:<8} "
                  f"-> {row['outcome']}")
        for row in report["screen"]:
            print(f"  screen  {row['candidate']:<26} {row['verdict']}")
        for element, res in report["adequacy"]["elements"].items():
            if not res["screened"]:
                print(f"  {element:<6} not screened")
                continue
            crit = res["critical_runoff_mm_per_earth_year"]
            print(f"  {element:<6} pool {res['standing_pool_g_per_m2']:8.1f} g/m2, "
                  f"critical runoff {min(crit.values()):7.1f} to "
                  f"{max(crit.values()):9.1f} mm/Earth yr")
        undeclared = report["undeclared_terms"]
        print(f"  {len(undeclared)} of {len(decl['terms'])} terms undeclared; "
              "the ledger is defined and does not close")
        for problem in hard:
            print(f"  FAIL {problem}")
        print(f"wrote {REPORT.relative_to(PROJECT_ROOT)}")

    if hard:
        return 1
    if args.strict and report["undeclared_terms"]:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
