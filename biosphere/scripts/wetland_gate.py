"""The wetland, peat and methane activation gate: what has to be declared and
what has to be REPAIRED before the simulated biosphere's peatland and methane
modules may be switched on, and what an activated run has to emit before its
output is accepted.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a vegetation model's dormant peatland and methane modules,
the declaration file that governs them, and the output tables they would write.

`run_peatland 1` plus `ifmethane 1` reads like two switches and is four models:
where the simulated wetlands are, how water reaches and leaves them, how peat
carbon and redox make and consume CH4, and what an atmosphere does with the
flux. The argument is in `biosphere/notes/wetland-activation-contract.md`; the
finding is `biosphere/notes/wetlands-peat-methane-audit.md`.

This module is the enforcement, and it refuses on two different grounds. It
reads `biosphere/config/wetlands.yaml` and refuses, by name, for every
precondition still carrying the `undeclared` sentinel. It ALSO reads the
vendored source and refuses for every hydrology defect still present there --
water created rather than routed, a constant runon, an annual water-table
average on an ordinal the calendar never reaches, and prognostic peat hydrology
absent from the serializer. A declaration that says a defect is closed while the
probe still finds it is reported as a contradiction rather than believed.

One of the grounds it refuses on is another component's declaration. A
wetland extent's saturated non-inundated mineral class needed a saturated
fraction, and the saturated-area closure that would have supplied one is
withdrawn in `hydrography/config/topographic_index.yaml`. This module reads
that file on every invocation rather than carrying a copy of the verdict, so a
revived closure refuses here until what it may key has been decided again.

It is fail-closed in one direction only. A biosphere run with peat and methane
OFF is a correct run, so with `requested: false` the gate reports and exits 0.
Only a request to activate can be refused.

A fixture set runs on every invocation. Most are mutations of the declaration,
built to be wrong in one named way each; one is the declaration as it stands,
which must be refused for some reason; six mutate the EVIDENCE instead, so
that what the gate does when a declaration disagrees with the vendored source,
or with what another component declares and publishes, stays checkable after
the repair or the withdrawal that made them agree; and the last is a met
declaration against a repaired source, which must be granted, because a gate
nothing can satisfy is a wall refusing for a reason that is never written down. A fixture that does not get
the verdict it was built for is a defect in this checker, and the gate exits
non-zero on it.

    python biosphere/scripts/wetland_gate.py           # status, exit 0
    python biosphere/scripts/wetland_gate.py --check-run biosphere/runs/<run_id>

`run_lpj_guess.py` calls `granted()` and writes `run_peatland 0`, `ifmethane 0`,
`ifsaturatewetlands 0` and `wetland_runon 0` whenever activation has not been
granted, which is every case that exists today.

LPJ-GUESS does not build on this tree until a baseline climatology exists, so
every probe below is a STATIC read of the vendored source. Nothing here compiles
or runs the model.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import CONFIG, GENERATED, GUESS_SOURCE, PROJECT_ROOT
from paths import rel  # noqa: E402

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = COMPONENT_ROOT / "config" / "wetlands.yaml"
REPORT = GENERATED / "wetland_gate_report.json"

UNDECLARED = "undeclared"

# The four instruction-file parameters activation moves. All four are top-level
# `declareitem` parameters that plib requires to be parsed, so an instruction
# file that omits one fails rather than defaults; `run_lpj_guess.py` writes all
# four rather than inheriting them.
SWITCHES = ("run_peatland", "ifmethane", "ifsaturatewetlands", "wetland_runon")

# The residual a closed ledger is allowed, as a fraction of the corresponding
# gross flux, and the LOOSEST the declaration may set it to. A water or C-N-P
# ledger that closes is an identity, so this is a floating-point tolerance and
# not a physical allowance; anything looser admits a real leak wearing the name
# of a rounding error. Fixed here, before any run exists.
RESIDUAL_TOLERANCE_CEILING = 1.0e-6

# The declaration fields whose sentinel refuses, and the bd issue that supplies
# each. The order is the order the refusals print in: the hydrology repair
# first, because saturation cannot be an outcome of a ledger that does not
# exist, then extent, then peat, then traits, then process, then the ledger,
# then the atmosphere it would close through.
PRECONDITIONS = (
    ("hydrology", "free_water_repair",
     "an infiltration path bounded by the water that actually arrived",
     "WET-3"),
    ("hydrology", "water_ledger",
     "the one daily exchange closing precipitation, routing, groundwater, "
     "uptake, storage and runoff", "WET-3"),
    ("hydrology", "runon_source",
     "what supplies runon, in place of the fixed wetland_runon scalar",
     "WET-3"),
    ("hydrology", "annual_water_table_update",
     "the annual water-table average, on an ordinal the calendar reaches",
     "WET-3"),
    ("hydrology", "restart_continuity_fixture",
     "a fixture requiring a resumed run to reproduce the uninterrupted one",
     "WET-3"),
    ("extent", "classification",
     "the mutually exclusive surface classification, as fractions", "WET-2"),
    ("extent", "downscaling_rule",
     "how a climate-grid fraction reaches a consumer that wants finer",
     "WORLD-D9U4"),
    ("extent", "regime_selector",
     "what selects the wetland process regime, which may not be latitude",
     "WET-2"),
    ("peat", "accumulation_rate_source",
     "where an accumulation rate comes from, given no time axis", "WET-4"),
    ("peat", "age_bracket", "both ends of the peat age bracket", "WET-4"),
    ("peat", "depth_bracket", "both ends of the peat depth bracket", "WET-4"),
    ("peat", "methane_active_column",
     "the methane-active column, kept apart from the total inventory",
     "WET-4"),
    ("traits", "registry",
     "the covarying plant trait registry the wetland traits live in",
     "PCAR-5"),
    ("traits", "coupled_traits",
     "roots, water-table tolerance, exudation, litter and gas transport, "
     "coupled through the registry", "WET-5"),
    ("traits", "transplanted_earth_bracket",
     "wetlandpfts.ins retained as one named transplanted-Earth bracket",
     "WET-5"),
    ("traits", "gas_transport_traits",
     "the aerenchyma, tiller-geometry and phenology traits gas transport "
     "reads", "WET-5"),
    ("process", "substrate_interface",
     "vertically resolved substrate, temperature, saturation and redox",
     "WET-6"),
    ("process", "oxidation_separate",
     "oxidation exposed as its own flux rather than folded into an emission "
     "factor", "WET-6"),
    ("process", "planet_gas_boundary",
     "local gravity, local pressure and the declared atmospheric gas state",
     "WET-7"),
    ("process", "transport_bracket",
     "winter and thaw pulses, dead-tiller venting, and bubbles crossing an "
     "unsaturated layer", "WET-7"),
    ("ledger", "surfaces",
     "the complete labelled surface ledger, on the same fractions as extent",
     "WET-8"),
    ("ledger", "dry_soil_sink",
     "the near-surface aerobic methane sink, which the module does not have",
     "WET-8"),
    ("ledger", "aquatic_sources",
     "lakes, rivers, reservoirs and playas, on the same fractions", "WET-8"),
    ("atmosphere", "oxidant_interface",
     "the reduced oxidant and lifetime bracket the ledger closes through",
     "WET-10"),
    ("acceptance", "declared_bracket_factor",
     "the declared uncertainty brackets on extent, flux, sink and production "
     "ratio", "WET-11"),
)

# The five classes WET-2's classification partitions the land into, and the
# closed set of fields a declaration may name as the source of each. A source
# named in prose cannot be checked against what hydrography publishes, so the
# set is closed and anything outside it is refused rather than believed.
#
# `topographic_index:f_sat` is in the set for one purpose: so that an extent
# keyed on a saturated fraction is refused BY NAME, with the withdrawal as the
# reason. Nothing may declare it.
EXTENT_CLASSES = ("peat_forming", "saturated_mineral", "seasonal_inundation",
                  "open_water", "dry_mineral")

EXTENT_SOURCES = {
    "surface_water:lake": (
        "surface_water.nc",
        "the solved equilibrium lake surface, one label per mesh region"),
    "surface_water:lake_cycle_fraction": (
        "surface_water.nc",
        "the per-bin inundated share of a closed basin's regions. It is the "
        "closed-basin part of seasonal inundation and no more: a floodplain's "
        "inundated area needs a height-above-nearest-drainage CDF and a "
        "routing model, and a seasonally saturated SOIL is a water content "
        "rather than an area"),
    "water_table:at_surface": (
        "water_table.nc",
        "where the solved water table meets the surface, a groundwater "
        "DISCHARGE extent per region, which needs no closure. It carries a "
        "regime licence: the depth field it comes from scores below no "
        "discrimination on one of the two Earth bore sets, so a cell mixing "
        "the licensed regime with the unlicensed one has no verdict attached"),
    "wetness:class_shares": (
        "wetness_*.nc",
        "the resolved wetness class shares over a climate-grid cell's land"),
    "lpj_guess:peatland_stand": (
        None,
        "the simulated peatland stand's own state, which is a model output "
        "rather than a hydrography field"),
    "residual": (
        None,
        "the land the resolved classes leave, which is a class only where "
        "every other class in the partition has a source"),
    "topographic_index:f_sat": (
        "topographic_index_*.nc",
        "WITHDRAWN. No consumer gets a saturated fraction from hydrography"),
}
WITHDRAWN_SOURCE = "topographic_index:f_sat"
# The one class that source would have formed, named in
# `hydrography/config/wetness.yaml` as the class the withdrawal leaves absent.
# It is here so that the refusal a reader actually meets carries the
# withdrawal, rather than a blank field's sentinel.
WITHDRAWN_CLASS = "saturated_mineral"
CLOSURE_DECLARATION = (PROJECT_ROOT / "hydrography" / "config" /
                       "topographic_index.yaml")

# One hydrology defect, one probe, one declaration field. A probe that still
# matches refuses while the field is undeclared, and CONTRADICTS the field once
# it is declared: a declaration cannot close a defect the source still carries.
HYDROLOGY_DEFECTS = (
    ("free_water_repair", "FREE-WATER",
     "modules/soilwater.cpp adds each layer's full saturation deficit to a "
     "low-latitude wetland stand whether or not rain_melt could supply it, so "
     "the simulated wetland creates water"),
    ("runon_source", "CONSTANT-RUNON",
     "modules/soil.cpp assigns soiltype.runon = wetland_runon, a namelist "
     "scalar in mm/day applied to every wetland stand and carrying no "
     "catchment"),
    ("annual_water_table_update", "UNREACHABLE-ANNUAL-UPDATE",
     "modules/soil.cpp guards the annual water-table average on "
     "date.day == Date::MAX_YEAR_LENGTH, which Date::next() never produces, "
     "so awtp holds its initial 0.0 and pins the simulated acrotelm CO2 to "
     "the pore-water value"),
)


class Refusal:
    """One named reason activation cannot be granted.

    `code` is what the error says and what a later report is grepped for; it is
    stable and does not move when the prose around it is rewritten.
    """

    def __init__(self, code: str, detail: str, issue: str | None = None):
        self.code = code
        self.detail = detail
        self.issue = issue

    def __str__(self) -> str:
        tail = f"  [{self.issue}]" if self.issue else ""
        return f"{self.code}: {self.detail}{tail}"

    def as_dict(self) -> dict:
        return {"code": self.code, "detail": self.detail, "issue": self.issue}


def read_declaration(path: Path = DECLARATION) -> dict:
    if not path.is_file():
        # The declaration's absence is itself a refusal, not an empty default.
        raise SystemExit(f"WET-DECLARATION-MISSING: {rel(path)} does not "
                         f"exist. There is no default declaration.")
    return yaml.safe_load(path.read_text())


def _get(declaration: dict, section: str, field: str):
    return (declaration.get(section) or {}).get(field, UNDECLARED)


def _declared(value) -> bool:
    return value is not None and value != UNDECLARED


def _read(path: Path) -> str:
    return path.read_text(errors="replace") if path.is_file() else ""


def _serialized_members(soil_source: str) -> set[str]:
    """Every member name `Soil::serialize` streams into the archive.

    Parsed rather than grepped over the whole file, because the point of the
    check is what the RESTART carries: a member named anywhere else in
    soil.cpp is not restart state.
    """
    start = soil_source.find("void Soil::serialize(")
    if start < 0:
        return set()
    depth = 0
    end = start
    for index in range(soil_source.find("{", start), len(soil_source)):
        if soil_source[index] == "{":
            depth += 1
        elif soil_source[index] == "}":
            depth -= 1
            if depth == 0:
                end = index
                break
    body = soil_source[start:end]
    return set(re.findall(r"&\s*([A-Za-z_]\w*)", body))


def probes() -> dict:
    """Read-only evidence, and the source facts the refusals are taken from.

    A probe never grants: the declaration decides, and a probe that stops
    matching because someone reworded a line must not be able to grant what a
    declaration still refuses. What the probes are for is putting the current
    state of the vendored source next to the refusal, so the reason is
    checkable against the artifact rather than against this file.
    """
    found: dict = {}

    soilwater = _read(GUESS_SOURCE / "modules" / "soilwater.cpp")
    if soilwater:
        # The low-latitude block clamps rain_melt to zero when it is smaller
        # than the total deficit, then adds the whole deficit anyway. Both
        # halves are required before the defect is claimed.
        clamps = bool(re.search(
            r"if\s*\(\s*soil\.rain_melt\s*<\s*total_potential\s*\)\s*"
            r"soil\.rain_melt\s*=\s*0\.0\s*;", soilwater))
        adds_full = bool(re.search(
            r"water_input_ly\s*=\s*potential_layer\[\s*ly\s*\]\s*;",
            soilwater))
        records_only = bool(re.search(
            r"if\s*\(\s*ifsaturatewetlands\s*\)\s*\n\s*"
            r"patch\.wetland_water_added_today\s*=\s*total_potential\s*;",
            soilwater))
        found["free_water_in_the_wetland_infiltration_path"] = {
            "rain_melt_clamped_to_zero": clamps,
            "full_deficit_added_regardless": adds_full,
            "switch_only_records_the_created_water": records_only,
            "defect_present": clamps and adds_full,
        }

    soil = _read(GUESS_SOURCE / "modules" / "soil.cpp")
    if soil:
        found["constant_wetland_runon"] = bool(re.search(
            r"soiltype\.runon\s*=\s*wetland_runon\s*;", soil))
        # The defect is the ORDINAL, not the constant. `day` runs 0 through
        # MAX_YEAR_LENGTH - 1, so a guard on MAX_YEAR_LENGTH itself never fires
        # and a guard on MAX_YEAR_LENGTH - 1 fires on the last simulated day of
        # the year. The lookahead is what tells those two apart.
        found["annual_water_table_guard"] = bool(re.search(
            r"date\.day\s*==\s*Date::MAX_YEAR_LENGTH(?!\s*-)", soil))
        # The whole set, not a count: `_source_refusals` reads it from here so
        # that every refusal is a function of the declaration and this evidence
        # and of nothing else, which is what makes the gate testable against a
        # repaired source without a repaired source to hand.
        found["serialized_soil_members"] = sorted(_serialized_members(soil))

    guessh = _read(GUESS_SOURCE / "framework" / "guess.h")
    if guessh:
        # Why the guard above is unreachable: `day` resets to 0 on the last day
        # of the last month, so it never reaches MAX_YEAR_LENGTH. Probed rather
        # than asserted, because the calendar is this project's own port.
        found["day_resets_to_zero_at_year_end"] = bool(re.search(
            r"if\s*\(islastmonth\)\s*\{[^}]*\bday\s*=\s*0\s*;", guessh,
            re.DOTALL))

    guesscpp = _read(GUESS_SOURCE / "framework" / "guess.cpp")
    if guesscpp:
        found["latitude_selects_the_wetland_regime"] = bool(re.search(
            r"landcover\s*==\s*PEATLAND\s*&&\s*lat\s*[<>]=?\s*"
            r"PEATLAND_WETLAND_LATITUDE_LIMIT", guesscpp))

    soilh = _read(GUESS_SOURCE / "modules" / "soil.h")
    if soilh:
        ratios = {}
        for name in ("CH4toCO2_peat", "CH4toCO2_inundated"):
            match = re.search(rf"^const\s+double\s+{name}\s*=\s*([0-9.eE+-]+)",
                              soilh, re.MULTILINE)
            if match:
                ratios[name] = float(match.group(1))
        found["production_ratio_constants"] = ratios
        limit = re.search(
            r"const\s+double\s+PEATLAND_WETLAND_LATITUDE_LIMIT\s*=\s*"
            r"([0-9.]+)", soilh)
        found["peatland_wetland_latitude_limit"] = (
            float(limit.group(1)) if limit else None)

    # Every output-table parameter commonoutput.cpp declares to plib.
    #
    # This is what decides whether a retained table can be ASKED FOR at all.
    # `run_lpj_guess.py` derives an instruction row from each retained
    # filename's stem -- `mwtp.out` becomes `file_mwtp "mwtp.out"` -- and plib
    # rejects an instruction file naming a parameter nothing declared. So a
    # retained output with no emitter does not produce an empty table; it
    # produces a run that aborts while parsing, naming one unknown parameter and
    # nothing about why it is unknown.
    common = _read(GUESS_SOURCE / "modules" / "commonoutput.cpp")
    if common:
        found["declared_output_parameters"] = sorted(set(re.findall(
            r'declare_parameter\s*\(\s*"(file_\w+)"', common)))

    wetpfts = GUESS_SOURCE / "data" / "ins" / "wetlandpfts.ins"
    if wetpfts.is_file():
        text = wetpfts.read_text(errors="replace")
        found["transplanted_earth_wetland_pfts"] = {
            "pfts": len(re.findall(r'^\s*pft\s*"', text, re.MULTILINE)),
            "with_aerenchyma": len(re.findall(r"has_aerenchyma\s+1", text)),
        }

    # ANOTHER COMPONENT'S DECLARATION, READ RATHER THAN RESTATED. Whether a
    # saturated fraction exists is hydrography's to say, and a refusal here
    # carrying its own copy of the verdict would go on refusing after a revival
    # and would grant after a second withdrawal, with nothing objecting either
    # time. Unreadable is its own answer and not a default.
    closure = {"declaration": str(rel(CLOSURE_DECLARATION)), "readable": False}
    if CLOSURE_DECLARATION.is_file():
        try:
            block = (yaml.safe_load(CLOSURE_DECLARATION.read_text())
                     or {}).get("closure") or {}
        except yaml.YAMLError:
            block = {}
        if block:
            closure.update(readable=True,
                           status=str(block.get("status", "active")),
                           withdrawn_on=str(block.get("withdrawn_on")),
                           withdrawn_because=block.get("withdrawn_because"))
    found["saturated_area_closure"] = closure

    # Which of the artifacts the closed source set names exist on any build. A
    # field nothing has produced is not a source, whatever the declaration says
    # the extent is keyed on.
    data = PROJECT_ROOT / "hydrography" / "data"
    found["hydrography_published"] = {
        token: sorted({p.parent.name for p in data.glob(f"*/{artifact}")})
        for token, (artifact, _what) in EXTENT_SOURCES.items()
        if artifact is not None}
    return found


def evaluate(declaration: dict, planet: dict, evidence: dict | None = None
             ) -> list[Refusal]:
    """Every reason activation cannot be granted, in declaration order."""
    refusals: list[Refusal] = []
    evidence = probes() if evidence is None else evidence

    for section, field, what, issue in PRECONDITIONS:
        if not _declared(_get(declaration, section, field)):
            refusals.append(Refusal(
                f"WET-UNDECLARED-{section.upper()}-"
                f"{field.upper().replace('_', '-')}",
                f"{what} is undeclared", issue))

    refusals.extend(_source_refusals(declaration, evidence))
    refusals.extend(_extent_refusals(declaration, evidence))
    refusals.extend(_peat_refusals(declaration))
    refusals.extend(_trait_refusals(declaration))
    refusals.extend(_acceptance_refusals(declaration, evidence))

    # The simulated atmosphere's methane is prescribed from a photochemical
    # calculation that holds modern Earth's biogenic surface fluxes fixed.
    # Driving an atmospheric abundance from this world's own simulated wetlands
    # while that stays fixed is not one atmospheric state, and the
    # inconsistency is the thing that must not be combined silently.
    if not (declaration.get("atmosphere") or {}).get("trace_gas_state_closed"):
        prescribed = [k for k in ("pCH4_bar", "ozone")
                      if k in (planet.get("atmosphere") or {})]
        if prescribed:
            refusals.append(Refusal(
                "WET-TRACE-GAS-STATE-NOT-CLOSED",
                "config/planet.yaml prescribes " + ", ".join(prescribed) +
                " from a fixed-Earth-flux photochemical calculation, so a "
                "wetland-derived methane source cannot be combined with it as "
                "one atmospheric state", "WET-10 / BVOC-6"))
    return refusals


def _source_refusals(declaration: dict, evidence: dict) -> list[Refusal]:
    """The hydrology defects still present in the vendored source.

    These are the half of the gate that can fail on the artifact rather than on
    a sentinel. Each has an answer that is wrong rather than merely different:
    either the source still adds water it did not receive, or it does not.
    """
    refusals: list[Refusal] = []
    free_water = evidence.get("free_water_in_the_wetland_infiltration_path")
    present = {
        "free_water_repair": bool(free_water and free_water["defect_present"]),
        "runon_source": bool(evidence.get("constant_wetland_runon")),
        "annual_water_table_update": bool(
            evidence.get("annual_water_table_guard")
            and evidence.get("day_resets_to_zero_at_year_end")),
    }
    for field, code, detail in HYDROLOGY_DEFECTS:
        if not present[field]:
            continue
        if _declared(_get(declaration, "hydrology", field)):
            refusals.append(Refusal(
                "WET-DECLARATION-CONTRADICTED",
                f"hydrology.{field} is declared closed while the source still "
                f"carries the defect: {detail}", "WET-3"))
        else:
            refusals.append(Refusal(f"WET-SOURCE-{code}", detail, "WET-3"))

    # Restart state. A member absent from Soil::serialize is not restart state,
    # whatever else the module does with it, so an arbitrary-day restart
    # resumes a different simulated state than the run it continues.
    required = (declaration.get("hydrology") or {}).get(
        "required_restart_members") or []
    carried = evidence.get("serialized_soil_members")
    if carried is not None and required:
        missing = [name for name in required if name not in carried]
        if missing:
            refusals.append(Refusal(
                "WET-SOURCE-RESTART-UNSERIALIZED",
                f"Soil::serialize does not carry {', '.join(missing)}, so a "
                "resumed run does not continue the simulated peat hydrology it "
                "stopped", "WET-3"))

    # Latitude selecting physics. It may locate a cell for geometry and
    # diagnostics; it may not decide which hydrology, decomposition and methane
    # model a stand runs, and on signed latitude a cell at 60 S takes the
    # low-latitude path.
    if evidence.get("latitude_selects_the_wetland_regime"):
        limit = evidence.get("peatland_wetland_latitude_limit")
        refusals.append(Refusal(
            "WET-SOURCE-LATITUDE-REGIME",
            "framework/guess.cpp selects the wetland hydrology, decomposition "
            f"and methane model on signed latitude against a fixed limit of "
            f"{limit}, so a simulated cell in the southern hemisphere takes "
            "the low-latitude path whatever its climate", "WET-2"))
    return refusals


def _withdrawal(closure: dict) -> str:
    """Why no saturated fraction is published, and what would reopen one.

    Read from hydrography's own declaration on every invocation, so a revival
    or a second withdrawal reaches the refusal text without anything here being
    edited. The reopening route is stated positively and its terms belong to
    `hydrography/notes/subgrid-water-table.md` section 7, which is cited rather
    than restated.
    """
    where = closure.get("declaration",
                        "hydrography/config/topographic_index.yaml")
    return (
        f"{where} carries closure.status {closure.get('status')!r} "
        f"({closure.get('withdrawn_because')}), "
        "hydrography/scripts/build_topographic_index.py implements no "
        "closure, and hydrography/scripts/build_wetness.py forms no saturated "
        "class and refuses to run if that status ever says anything else. "
        "hydrography/notes/subgrid-water-table.md section 7 measures why no "
        "narrower criterion can license one: the gain the terrain half "
        "carries over the cell-mean depth is at or below the scatter its own "
        "support puts on it on every arm of both bore sets, and a narrower "
        "criterion has less support rather than more. What reopens it is a "
        "score whose SUPPORT can resolve a gain of a few thousandths of an "
        "AUC -- more cells carrying observations than the 34 and 37 two "
        "continents of bores give, or an areal saturation or inundation "
        "observation in place of point depths -- with its criterion declared "
        "before those "
        "observations are in hand. The other route is a reduced form declared "
        "under WET-12 with what dropping the class costs in claims")


def _extent_refusals(declaration: dict, evidence: dict) -> list[Refusal]:
    """What may key an extent, and whether the fields it names exist."""
    refusals: list[Refusal] = []
    extent = declaration.get("extent") or {}

    support = extent.get("saturated_fraction_support", UNDECLARED)
    if _declared(support) and str(support) != "climate_grid":
        refusals.append(Refusal(
            "WET-EXTENT-SUPPORT-NOT-GRID",
            f"extent.saturated_fraction_support is {support!r}. A saturated "
            "fraction is a fraction of a population, the only population this "
            "project holds is the mesh regions inside a climate-grid cell, and "
            "a region has no sub-population: the sub-region hypsometry a "
            "native-mesh fraction would need does not exist and is not "
            "recoverable by a finer generation", "WORLD-D9U4 / GW-6"))

    # THE SATURATED-AREA CLOSURE IS WITHDRAWN, and this is where a reader of a
    # refusal meets that. The support statement above stays true and is no
    # longer the operative one: hydrography publishes no saturated fraction at
    # ANY support, so an extent keyed on one has no source rather than the
    # wrong support. The verdict is read from hydrography's own declaration on
    # every invocation, so a revival cannot pass unnoticed in either direction.
    closure = evidence.get("saturated_area_closure") or {}
    where = closure.get("declaration",
                        "hydrography/config/topographic_index.yaml")
    if not closure.get("readable"):
        refusals.append(Refusal(
            "WET-EXTENT-CLOSURE-UNREADABLE",
            f"{where} carries no readable closure block, so whether "
            "hydrography publishes a saturated fraction cannot be "
            "established here. It is refused rather than assumed either way",
            "WET-2 / GW-26"))
    elif closure.get("status") != "withdrawn":
        refusals.append(Refusal(
            "WET-EXTENT-CLOSURE-REVIVED",
            f"{where} closure.status is {closure.get('status')!r}, and every "
            "extent refusal here is written on the withdrawal. A revived "
            "closure is a new closure with a new criterion, and what that "
            "criterion may key in a wetland extent is a decision to be taken "
            "again before any fraction is read", "WET-2 / GW-26"))

    # One source per class, from the closed set, and each source checked
    # against what hydrography has actually published. The classification is a
    # PARTITION, so a class with no named source is land the partition cannot
    # account for, and a class declared absent is a reduced form rather than an
    # omission.
    sources = extent.get("class_sources")
    if not isinstance(sources, dict):
        refusals.append(Refusal(
            "WET-EXTENT-CLASS-SOURCES-UNDECLARED",
            "extent.class_sources is not a mapping from each of " +
            ", ".join(EXTENT_CLASSES) + " to the field it is taken from",
            "WET-2"))
        sources = {}
    published = evidence.get("hydrography_published") or {}
    for name in EXTENT_CLASSES:
        value = sources.get(name, UNDECLARED)
        if isinstance(value, dict):
            if not (value.get("absent_because") and value.get("licensed_by")):
                refusals.append(Refusal(
                    "WET-EXTENT-CLASS-ABSENT-UNLICENSED",
                    f"extent.class_sources.{name} declares the class absent "
                    "without both `absent_because` and `licensed_by`. An "
                    "extent that drops one of its classes is a reduced form, "
                    "and WET-12 is where a reduced form is declared together "
                    "with what taking it costs in claims", "WET-12"))
            continue
        if not _declared(value):
            withdrawn = closure.get("status") == "withdrawn"
            if name == WITHDRAWN_CLASS and withdrawn:
                refusals.append(Refusal(
                    "WET-EXTENT-CLASS-SOURCE-WITHDRAWN",
                    f"extent.class_sources.{name} is undeclared and has "
                    "nothing to declare: hydrography publishes no saturated "
                    f"fraction at any support. {_withdrawal(closure)}",
                    "WET-2 / GW-26 / WET-12"))
            else:
                refusals.append(Refusal(
                    "WET-EXTENT-CLASS-SOURCE-UNDECLARED",
                    f"extent.class_sources.{name} is undeclared", "WET-2"))
            continue
        token = str(value)
        if token == WITHDRAWN_SOURCE or re.search(
                r"f_sat\b|saturated[ _]fraction", token, re.IGNORECASE):
            refusals.append(Refusal(
                "WET-EXTENT-SATURATED-FRACTION-WITHDRAWN",
                f"extent.class_sources.{name} is {token!r} and hydrography "
                f"publishes no saturated fraction: {_withdrawal(closure)}",
                "WET-2 / GW-26"))
            continue
        if token not in EXTENT_SOURCES:
            refusals.append(Refusal(
                "WET-EXTENT-CLASS-SOURCE-UNKNOWN",
                f"extent.class_sources.{name} is {token!r}, which is not one "
                "of " + ", ".join(sorted(EXTENT_SOURCES)) + ". The set is "
                "closed because a source named in prose cannot be checked "
                "against what hydrography publishes", "WET-2"))
            continue
        artifact, _what = EXTENT_SOURCES[token]
        if artifact is not None and not published.get(token):
            refusals.append(Refusal(
                "WET-EXTENT-CLASS-SOURCE-NOT-PUBLISHED",
                f"extent.class_sources.{name} is {token!r} and no build under "
                f"hydrography/data carries {artifact}. A field nothing has "
                "produced is not a source", "WET-2"))

    # The compound topographic index transports as a rank statistic and not as
    # an absolute threshold: `a` carries a length, so the whole distribution
    # shifts with the mesh, by about ln 2 between this project's own two builds.
    if str(extent.get("absolute_index_threshold", "forbidden")) != "forbidden":
        refusals.append(Refusal(
            "WET-EXTENT-ABSOLUTE-INDEX-PERMITTED",
            "extent.absolute_index_threshold is not `forbidden`, and an "
            "absolute compound-topographic-index cut is a statement about two "
            "indices' scales rather than about this world", "WET-2"))
    keyed = " ".join(str(extent.get(k, "")) for k in
                     ("classification", "downscaling_rule"))
    if re.search(r"cti_\w*crit|absolute (?:cti|index) (?:threshold|cut)"
                 r"|index\s*[<>]=?\s*[0-9]", keyed, re.IGNORECASE):
        refusals.append(Refusal(
            "WET-EXTENT-ABSOLUTE-INDEX",
            "the declared extent model keys on an absolute compound "
            "topographic index value; the index shifts with the mesh by about "
            "ln 2 and only its rank statistics transport", "WET-2"))

    selector = extent.get("regime_selector", UNDECLARED)
    if _declared(selector) and re.search(r"latitude|\blat\b", str(selector),
                                         re.IGNORECASE):
        refusals.append(Refusal(
            "WET-REGIME-SELECTOR-LATITUDE",
            f"extent.regime_selector is {selector!r}. Latitude may locate a "
            "simulated cell for geometry and diagnostics and may not select an "
            "ecological process regime", "WET-2"))
    return refusals


def _peat_refusals(declaration: dict) -> list[Refusal]:
    """Peat age and depth are brackets, because this world has no time axis."""
    refusals: list[Refusal] = []
    peat = declaration.get("peat") or {}
    for field in ("age_bracket", "depth_bracket"):
        value = peat.get(field, UNDECLARED)
        if not _declared(value):
            continue
        if not isinstance(value, (list, tuple)) or len(value) != 2:
            refusals.append(Refusal(
                "WET-PEAT-NOT-BRACKETED",
                f"peat.{field} is {value!r} rather than two ends. Peat depth "
                "is a rate over a duration and a duration is undefined here, "
                "so a scalar is a claim to a history this project does not "
                "keep", "WET-4"))
        elif float(value[0]) >= float(value[1]):
            refusals.append(Refusal(
                "WET-PEAT-BRACKET-DEGENERATE",
                f"peat.{field} has ends {value[0]} and {value[1]}. A bracket "
                "whose ends agree is a central value with two names", "WET-4"))
    return refusals


def _trait_refusals(declaration: dict) -> list[Refusal]:
    """The strategies the registry has to be able to tell apart."""
    refusals: list[Refusal] = []
    traits = declaration.get("traits") or {}
    required = set(traits.get("required_strategies") or [])
    # `non_vegetated_inundation` is a state with an area and no plants. It is
    # required by name because leaving it out is exactly how inundated area
    # acquires a vegetation type by default.
    for name in ("bog_nutrient_strategy", "fen_nutrient_strategy",
                 "aerenchymatous_emergent", "non_vegetated_inundation"):
        if name not in required:
            refusals.append(Refusal(
                "WET-STRATEGY-MISSING",
                f"traits.required_strategies does not require {name}, so the "
                "registry could satisfy this declaration without being able to "
                "tell it apart from anything else", "WET-5"))
    return refusals


def _acceptance_refusals(declaration: dict, evidence: dict) -> list[Refusal]:
    """The bracket floors, the residual tolerance, and the registered arms."""
    refusals: list[Refusal] = []
    acceptance = declaration.get("acceptance") or {}

    floors = acceptance.get("minimum_bracket_factor") or {}
    declared = acceptance.get("declared_bracket_factor")
    if isinstance(declared, dict):
        for name, floor in floors.items():
            value = declared.get(name)
            if value is None:
                refusals.append(Refusal(
                    "WET-BRACKET-MISSING",
                    f"declared_bracket_factor has no entry for {name}",
                    "WET-11"))
            elif float(value) < float(floor):
                refusals.append(Refusal(
                    "WET-BRACKET-TOO-NARROW",
                    f"{name} declared at factor {value}, below the declared "
                    f"floor of {floor}", "WET-11"))

    # The production-ratio floor is not a literature range: it is the spread
    # between the two constants the vendored fork itself carries. If either
    # constant moves, the floor stops being what it says it is, and a floor
    # that has quietly stopped matching its own derivation is worse than none.
    ratios = evidence.get("production_ratio_constants") or {}
    floor = floors.get("production_ratio")
    if floor and len(ratios) == 2:
        peat = ratios.get("CH4toCO2_peat")
        inundated = ratios.get("CH4toCO2_inundated")
        if peat and inundated:
            actual = max(peat, inundated) / min(peat, inundated)
            if abs(actual - float(floor)) > 0.01 * float(floor):
                refusals.append(Refusal(
                    "WET-BRACKET-FLOOR-DRIFT",
                    f"the production_ratio floor is {floor}, and the fork's "
                    f"own CH4:CO2 constants now stand at {peat} and "
                    f"{inundated}, a spread of {actual:.3f}", "WET-11"))

    # Every retained table, against the emitter that would have to write it.
    #
    # `retained_output_status` says, per table, whether the vendored fork can
    # emit it and which issue supplies the one that cannot. Both halves are
    # checked against commonoutput.cpp, so neither a table that quietly lost its
    # emitter nor one that quietly gained one can sit here unnoticed. Activation
    # refuses while any is missing, which turns a plib parse error deep inside a
    # run into a refusal that names each table and what it waits on.
    status = acceptance.get("retained_output_status")
    declared_parameters = evidence.get("declared_output_parameters")
    retained = list(acceptance.get("retained_outputs") or [])
    if status is None:
        refusals.append(Refusal(
            "WET-OUTPUT-STATUS-UNDECLARED",
            "acceptance.retained_output_status is absent, so nothing says "
            "which retained tables the fork can emit and a run would find out "
            "by failing to parse its own instruction file", "WORLD-E5T5"))
    elif declared_parameters is not None:
        for name in retained:
            parameter = f"file_{name.split('.')[0]}"
            entry = (status or {}).get(name)
            has_emitter = parameter in declared_parameters
            if entry is None:
                refusals.append(Refusal(
                    "WET-OUTPUT-STATUS-MISSING",
                    f"{name} is retained and retained_output_status says "
                    "nothing about it", "WORLD-E5T5"))
                continue
            claims_emitter = bool(entry.get("emitted"))
            if claims_emitter and not has_emitter:
                refusals.append(Refusal(
                    "WET-OUTPUT-NO-EMITTER",
                    f"{name} is declared emitted and commonoutput.cpp declares "
                    f"no {parameter}. run_lpj_guess.py derives that row from "
                    "the filename, so an activated run would write an "
                    "instruction file plib rejects", "WORLD-E5T5"))
            elif not claims_emitter and has_emitter:
                refusals.append(Refusal(
                    "WET-OUTPUT-STATUS-STALE",
                    f"{name} is declared unemitted and commonoutput.cpp now "
                    f"declares {parameter}. A table that has become available "
                    "and is still recorded as waiting is one nobody will think "
                    "to ask for", "WORLD-E5T5"))
            elif not claims_emitter:
                if not entry.get("waits_on"):
                    refusals.append(Refusal(
                        "WET-OUTPUT-UNCLAIMED",
                        f"{name} has no emitter and names no issue that would "
                        "supply one", "WORLD-E5T5"))
                else:
                    refusals.append(Refusal(
                        "WET-OUTPUT-ABSENT",
                        f"{name} has no emitter in the vendored fork: "
                        f"{entry.get('why', 'no reason recorded')}",
                        str(entry.get("waits_on"))))

    tolerance = acceptance.get("residual_tolerance")
    if tolerance is None:
        refusals.append(Refusal(
            "WET-RESIDUAL-TOLERANCE-UNDECLARED",
            "acceptance.residual_tolerance is absent, so a ledger has no bar "
            "to close against", "WET-11"))
    elif float(tolerance) > RESIDUAL_TOLERANCE_CEILING:
        refusals.append(Refusal(
            "WET-RESIDUAL-TOLERANCE-TOO-LOOSE",
            f"acceptance.residual_tolerance is {tolerance}, above the "
            f"{RESIDUAL_TOLERANCE_CEILING} ceiling. A closed ledger is an "
            "identity, so the tolerance is a floating-point tolerance and not "
            "a physical allowance", "WET-11"))

    # The arms are registered before any of them can be run, and an arm chosen
    # after its result is seen is not an arm. Requiring both extent ends by name
    # is the specific case the audit's finding 11 is about: a total that agrees
    # while extent and production compensate.
    arms = set(acceptance.get("registered_arms") or [])
    for name in ("no_methane", "extent_low_end", "extent_high_end",
                 "transport_partition"):
        if name not in arms:
            refusals.append(Refusal(
                "WET-ARM-UNREGISTERED",
                f"acceptance.registered_arms does not register {name}, so an "
                "aggregate agreement could not be told apart from compensating "
                "errors in the components", "WET-11"))
    return refusals


def granted(declaration: dict | None = None,
            planet: dict | None = None) -> tuple[bool, list[Refusal]]:
    """Whether an activated run may proceed, and why not.

    `run_lpj_guess.py` calls this. A declaration that does not request
    activation is granted nothing and refuses nothing: the run writes the four
    switches at zero and is a correct run.
    """
    declaration = declaration if declaration is not None else read_declaration()
    planet = planet if planet is not None else yaml.safe_load(CONFIG.read_text())
    refusals = evaluate(declaration, planet)
    if not declaration.get("requested"):
        return False, refusals
    return not refusals, refusals


def require(declaration: dict | None = None,
            planet: dict | None = None) -> bool:
    """`granted`, but a refused REQUEST exits rather than returning False."""
    declaration = declaration if declaration is not None else read_declaration()
    ok, refusals = granted(declaration, planet)
    if declaration.get("requested") and not ok:
        lines = "\n  ".join(str(r) for r in refusals)
        raise SystemExit(
            f"{rel(DECLARATION)} requests wetland, peat and methane activation "
            f"and the gate refuses it for {len(refusals)} reasons:\n  {lines}"
            "\n\nEach is a declaration to make in that file, or a repair to "
            "make in the vendored source, with the bd issue that supplies it. "
            "Nothing here has a default value.")
    return ok


def switches(active: bool) -> dict[str, int | float]:
    """What `run_lpj_guess.py` writes into the instruction file.

    All four, always, whether or not activation was granted. plib takes the
    later declaration, so these override whatever the imported files carry, and
    an explicit zero is the difference between a decision and an inheritance.
    """
    if not active:
        return {"run_peatland": 0, "ifmethane": 0, "ifsaturatewetlands": 0,
                "wetland_runon": 0.0}
    # Granted activation still writes ifsaturatewetlands 0 and wetland_runon 0:
    # both belong to the path the hydrology repair replaces, and neither is a
    # source of water once the ledger owns it.
    return {"run_peatland": 1, "ifmethane": 1, "ifsaturatewetlands": 0,
            "wetland_runon": 0.0}


# --------------------------------------------------------------------------
# Acceptance, for output that does not exist yet.

def check_run(run_dir: Path, declaration: dict | None = None) -> list[Refusal]:
    """Reject a run's wetland output rather than accept it by default.

    Every check has an answer that is wrong rather than merely different: a
    simulated cell present in the productivity output and absent here is a
    missing cell, a negative area or a negative production is impossible, a
    monthly table without one column per simulated month is not this world's
    year, a residual above the declared tolerance is an unclosed ledger, and a
    run whose manifest carries no forcing hash cannot be attributed to a
    climate.
    """
    declaration = declaration if declaration is not None else read_declaration()
    acceptance = declaration.get("acceptance") or {}
    expected = list(acceptance.get("retained_outputs") or [])
    tolerance = float(acceptance.get("residual_tolerance")
                      or RESIDUAL_TOLERANCE_CEILING)

    refusals: list[Refusal] = []
    run_dir = Path(run_dir)

    missing = [name for name in expected if not (run_dir / name).is_file()]
    if missing:
        refusals.append(Refusal(
            "WET-OUTPUT-MISSING",
            f"{len(missing)} of {len(expected)} retained tables absent from "
            f"{rel(run_dir)}: " + ", ".join(missing[:4]) +
            ("..." if len(missing) > 4 else ""), "WET-11"))

    reference = run_dir / "anpp.out"
    reference_cells: set[tuple[str, str]] = set()
    if reference.is_file():
        for line in reference.read_text().splitlines()[1:]:
            parts = line.split()
            if len(parts) >= 2:
                reference_cells.add((parts[0], parts[1]))
    else:
        refusals.append(Refusal(
            "WET-REFERENCE-MISSING",
            f"{rel(run_dir)}/anpp.out is absent, so the simulated cell set the "
            "wetland output has to cover is unknown"))

    months = _simulated_months()
    for name in expected:
        path = run_dir / name
        if not path.is_file():
            continue
        lines = path.read_text().splitlines()
        if not lines:
            refusals.append(Refusal("WET-OUTPUT-EMPTY", f"{name} is empty"))
            continue
        header = lines[0].split()
        residual = name.startswith("aresidual_")
        cells: set[tuple[str, str]] = set()
        for line in lines[1:]:
            parts = line.split()
            if len(parts) < 3:
                continue
            cells.add((parts[0], parts[1]))
            for token in parts[2:]:
                try:
                    value = float(token)
                except ValueError:
                    refusals.append(Refusal(
                        "WET-OUTPUT-UNREADABLE",
                        f"{name} carries a non-numeric value {token!r}"))
                    break
                if math.isnan(value) or math.isinf(value):
                    refusals.append(Refusal(
                        "WET-OUTPUT-NOT-FINITE",
                        f"{name} carries a non-finite value"))
                    break
                if residual:
                    if abs(value) > tolerance:
                        refusals.append(Refusal(
                            "WET-LEDGER-NOT-CLOSED",
                            f"{name} carries a residual of {value:.3g}, above "
                            f"the declared tolerance of {tolerance:g}",
                            "WET-11"))
                        break
                elif value < 0.0:
                    # An area, a stock and a gross pathway flux are all
                    # non-negative by construction; the signed net flux lives in
                    # the residual tables, which is why they are exempt above.
                    refusals.append(Refusal(
                        "WET-OUTPUT-IMPOSSIBLE",
                        f"{name} carries a negative area, stock or gross flux"))
                    break
        if reference_cells and cells != reference_cells:
            refusals.append(Refusal(
                "WET-OUTPUT-INCOMPLETE",
                f"{name} covers {len(cells)} simulated cells against "
                f"{len(reference_cells)} in anpp.out"))
        if name.startswith("m"):
            if months is None:
                refusals.append(Refusal(
                    "WET-CALENDAR-UNKNOWN",
                    "biosphere/generated/vesper.h is absent, so how many "
                    "months the simulated year has is unknown and a monthly "
                    "table cannot be accepted"))
            elif len(header) >= 2 and len(header) - 2 < months:
                refusals.append(Refusal(
                    "WET-OUTPUT-CALENDAR",
                    f"{name} has {len(header) - 2} data columns against "
                    f"{months} months in the simulated year"))

    manifest = run_dir / "run_manifest.json"
    if not manifest.is_file():
        refusals.append(Refusal(
            "WET-PROVENANCE-MISSING",
            f"{rel(run_dir)}/run_manifest.json is absent, so the forcing this "
            "methane flux came from is unrecorded"))
    else:
        record = json.loads(manifest.read_text())
        if not record.get("source_build"):
            refusals.append(Refusal(
                "WET-PROVENANCE-NO-BUILD",
                "the run manifest carries no source_build"))
        if not (record.get("inputs") or {}).get("driver", {}).get("sha256"):
            refusals.append(Refusal(
                "WET-PROVENANCE-NO-FORCING-HASH",
                "the run manifest carries no hash of the climate forcing"))
        arm = (record.get("physical") or {}).get("wetland_arm")
        registered = set(acceptance.get("registered_arms") or [])
        if arm is not None and arm not in registered:
            refusals.append(Refusal(
                "WET-ARM-NOT-REGISTERED",
                f"the run manifest names arm {arm!r}, which the declaration "
                "did not register before the run", "WET-11"))
    return refusals


def _simulated_months() -> int | None:
    """Months in the simulated year, from the generated header.

    Delegated to the BVOC gate rather than copied. There is one answer to how
    long this world's year is and one place that reads it; the note in
    `_paths.py` about four components keeping their own copy of a path resolver
    and drifting apart is the reason this is not a second copy.
    """
    import bvoc_gate
    return bvoc_gate._simulated_months()


# --------------------------------------------------------------------------
# The fixtures. Each is a declaration that MUST be refused by a named code, and
# each asserts that code appears. The first is the declaration as it stands,
# which has to be refused for the reason the contract says it is refused for;
# the rest are mutations built to be wrong in one way each. A fixture that does
# not get its verdict is a defect in this checker, not in the declaration.

def _fixtures(declaration: dict, planet: dict, evidence: dict) -> list[dict]:

    def mutate(fn):
        candidate = copy.deepcopy(declaration)
        fn(candidate)
        return candidate

    def put(section, field, value):
        return lambda d: d.setdefault(section, {}).__setitem__(field, value)

    def drop(section, field, item):
        return lambda d: d[section][field].remove(item)

    # `None` means "must be refused for at least one reason", which is the only
    # durable expectation for the live declaration: the specific code it earns
    # moves as fields are declared and repairs land, and a fixture that fails
    # because someone made progress is a fixture nobody will keep.
    cases = [
        ("the declaration as it stands", declaration, None),
        ("a saturated fraction taken on the native mesh",
         mutate(put("extent", "saturated_fraction_support", "native_mesh")),
         "WET-EXTENT-SUPPORT-NOT-GRID"),
        ("the saturated mineral class left with nothing to declare",
         mutate(lambda d: d["extent"]["class_sources"].__setitem__(
             "saturated_mineral", "undeclared")),
         "WET-EXTENT-CLASS-SOURCE-WITHDRAWN"),
        ("an extent class keyed on the withdrawn saturated fraction",
         mutate(lambda d: d["extent"]["class_sources"].__setitem__(
             "saturated_mineral", "topographic_index:f_sat")),
         "WET-EXTENT-SATURATED-FRACTION-WITHDRAWN"),
        ("an extent class keyed on a source named in prose",
         mutate(lambda d: d["extent"]["class_sources"].__setitem__(
             "open_water", "the lakes, obviously")),
         "WET-EXTENT-CLASS-SOURCE-UNKNOWN"),
        ("a class dropped from the extent with no reduced form to license it",
         mutate(lambda d: d["extent"]["class_sources"].__setitem__(
             "saturated_mineral", {"absent_because": "there is no closure"})),
         "WET-EXTENT-CLASS-ABSENT-UNLICENSED"),
        ("an extent keyed on an absolute topographic index cut",
         mutate(put("extent", "classification",
                    "TOPMODEL fraction with cti_mean_crit 5.5")),
         "WET-EXTENT-ABSOLUTE-INDEX"),
        ("the absolute-threshold prohibition relaxed",
         mutate(put("extent", "absolute_index_threshold", "allowed")),
         "WET-EXTENT-ABSOLUTE-INDEX-PERMITTED"),
        ("latitude selecting the wetland process regime",
         mutate(put("extent", "regime_selector", "latitude >= 40")),
         "WET-REGIME-SELECTOR-LATITUDE"),
        ("a scalar peat age",
         mutate(put("peat", "age_bracket", 8000)), "WET-PEAT-NOT-BRACKETED"),
        ("a peat depth bracket whose ends agree",
         mutate(put("peat", "depth_bracket", [1.5, 1.5])),
         "WET-PEAT-BRACKET-DEGENERATE"),
        ("a wetland extent bracket narrower than the WETCHIMP spread",
         mutate(put("acceptance", "declared_bracket_factor",
                    {"wetland_extent": 1.5, "wetland_methane_flux": 2.0,
                     "dry_soil_sink": 6.0, "production_ratio": 4.0})),
         "WET-BRACKET-TOO-NARROW"),
        ("a declared bracket set with an entry missing",
         mutate(put("acceptance", "declared_bracket_factor",
                    {"wetland_methane_flux": 2.0, "dry_soil_sink": 6.0,
                     "production_ratio": 4.0})),
         "WET-BRACKET-MISSING"),
        ("a residual tolerance loose enough to admit a leak",
         mutate(put("acceptance", "residual_tolerance", 0.01)),
         "WET-RESIDUAL-TOLERANCE-TOO-LOOSE"),
        ("the low-extent arm dropped from the registered set",
         mutate(drop("acceptance", "registered_arms", "extent_low_end")),
         "WET-ARM-UNREGISTERED"),
        ("non-vegetated inundation dropped from the required strategies",
         mutate(drop("traits", "required_strategies",
                     "non_vegetated_inundation")),
         "WET-STRATEGY-MISSING"),
        ("a production-ratio floor that no longer matches the fork constants",
         mutate(lambda d: d["acceptance"]["minimum_bracket_factor"]
                .__setitem__("production_ratio", 9.0)),
         "WET-BRACKET-FLOOR-DRIFT"),
        ("no statement at all about which retained tables the fork can emit",
         mutate(lambda d: d["acceptance"].pop("retained_output_status")),
         "WET-OUTPUT-STATUS-UNDECLARED"),
        ("a retained table the emitter statement says nothing about",
         mutate(lambda d: d["acceptance"]["retained_output_status"]
                .pop("mwtp.out")),
         "WET-OUTPUT-STATUS-MISSING"),
        ("a table with no emitter and no issue that would supply one",
         mutate(lambda d: d["acceptance"]["retained_output_status"]
                .__setitem__("apeat_stock.out", {"emitted": False})),
         "WET-OUTPUT-UNCLAIMED"),
        # Both directions against the LIVE source, which is what makes them
        # tests: the fork declares file_mwtp and does not declare
        # file_mch4_production, so each mutation is a claim the source refutes.
        ("a table declared emitted that commonoutput.cpp cannot be asked for",
         mutate(lambda d: d["acceptance"]["retained_output_status"]
                .__setitem__("mch4_production.out",
                             {"emitted": True, "quantity": "wishful"})),
         "WET-OUTPUT-NO-EMITTER"),
        ("a table still recorded as waiting after its emitter arrived",
         mutate(lambda d: d["acceptance"]["retained_output_status"]
                .__setitem__("mwtp.out",
                             {"emitted": False, "waits_on": "WET-3",
                              "why": "stale"})),
         "WET-OUTPUT-STATUS-STALE"),
    ]

    # Mutations of the EVIDENCE rather than the declaration. These carry their
    # own source facts because the property they check is what the gate does
    # when a declaration and the source disagree, and that cannot be shown
    # against a source where the defect is absent. They stay meaningful after a
    # repair lands, which is exactly what a fixture written against live
    # evidence stops doing the moment the repair it was built for arrives.
    # What the gate does when the PRODUCER'S declaration moves. The withdrawal
    # is read on every invocation, so both directions of a change to it have to
    # be checkable from here rather than from a copy of the verdict.
    revived = copy.deepcopy(evidence)
    revived["saturated_area_closure"] = dict(
        evidence.get("saturated_area_closure") or {},
        readable=True, status="active")
    unreadable = copy.deepcopy(evidence)
    unreadable["saturated_area_closure"] = {
        "declaration": "hydrography/config/topographic_index.yaml",
        "readable": False}
    unbuilt = dict(evidence)
    unbuilt["hydrography_published"] = dict(
        evidence.get("hydrography_published") or {})
    unbuilt["hydrography_published"]["surface_water:lake"] = []

    defective = _defective(evidence)
    cases += [
        ("the saturated-area closure revived in hydrography's declaration",
         declaration, "WET-EXTENT-CLOSURE-REVIVED", revived),
        ("hydrography's closure declaration unreadable",
         declaration, "WET-EXTENT-CLOSURE-UNREADABLE", unreadable),
        ("an extent class keyed on an artifact no build carries",
         _satisfied(declaration), "WET-EXTENT-CLASS-SOURCE-NOT-PUBLISHED",
         unbuilt),
        ("the free-water repair declared while the source still creates water",
         mutate(put("hydrology", "free_water_repair", "done, honestly")),
         "WET-DECLARATION-CONTRADICTED", defective),
        ("the runon source declared while the source still assigns the scalar",
         mutate(put("hydrology", "runon_source", "done, honestly")),
         "WET-DECLARATION-CONTRADICTED", defective),
        ("a serializer missing one required peat-hydrology member",
         declaration, "WET-SOURCE-RESTART-UNSERIALIZED",
         _unserialized(evidence, "Wtot")),
    ]

    results = []
    for case in cases:
        label, candidate, expect = case[0], case[1], case[2]
        facts = case[3] if len(case) > 3 else evidence
        codes = {r.code for r in evaluate(candidate, planet, facts)}
        ok = bool(codes) if expect is None else expect in codes
        results.append({"fixture": label,
                        "expected": expect or "refused for some reason",
                        "refusals": len(codes), "pass": ok})

    # The last fixture is the other direction, and it is the one that says this
    # is a gate rather than a wall: a declaration with every precondition met,
    # evaluated against evidence of a repaired source, must come back with NO
    # refusals. A gate nothing can satisfy refuses for a reason that is never
    # written down.
    codes = {r.code for r in
             evaluate(_satisfied(declaration), planet,
                      _repaired(evidence,
                                declaration["acceptance"]["retained_outputs"]))}
    results.append({
        "fixture": "a met declaration against a repaired source is granted",
        "expected": "no refusal", "refusals": len(codes),
        "pass": not codes, "codes": sorted(codes)})
    return results


def _satisfied(declaration: dict) -> dict:
    """The declaration with every precondition met, for the fixture above.

    Placeholders, not proposals: nothing here is a value anyone should copy into
    biosphere/config/wetlands.yaml. It exists so the gate can be shown to be
    satisfiable in principle.
    """
    candidate = copy.deepcopy(declaration)
    for section, field, _what, _issue in PRECONDITIONS:
        candidate.setdefault(section, {})[field] = "declared, for the fixture"
    # One source per class. The saturated mineral class is where this fixture
    # has to say something rather than fill a blank: it has no source at all
    # and cannot acquire one from a declaration, so the only route that
    # satisfies the gate is a reduced form declared under WET-12. That the
    # other four resolve says nothing about whether any of them IS a wetland
    # extent, which is WET-2's decision and not this fixture's.
    candidate["extent"]["class_sources"] = {
        "peat_forming": "lpj_guess:peatland_stand",
        "saturated_mineral": {
            "absent_because": "the saturated-area closure is withdrawn",
            "licensed_by": "WET-12, for the fixture"},
        "seasonal_inundation": "surface_water:lake_cycle_fraction",
        "open_water": "surface_water:lake",
        "dry_mineral": "residual"}
    candidate["extent"]["regime_selector"] = "saturation and its persistence"
    candidate["peat"]["age_bracket"] = [1, 2]
    candidate["peat"]["depth_bracket"] = [0.1, 10.0]
    candidate["atmosphere"]["trace_gas_state_closed"] = True
    floors = candidate["acceptance"]["minimum_bracket_factor"]
    candidate["acceptance"]["declared_bracket_factor"] = dict(floors)
    # Every retained table declared emitted, which `_repaired` matches by giving
    # the source the parameters to go with it. Both halves have to move together
    # or the fixture would pass on a declaration that lies about the fork.
    candidate["acceptance"]["retained_output_status"] = {
        name: {"emitted": True, "quantity": "declared, for the fixture"}
        for name in candidate["acceptance"]["retained_outputs"]}
    return candidate


def _defective(evidence: dict) -> dict:
    """The evidence the vendored source produced before the WET-3 repairs.

    The inverse of `_repaired`, and it exists for the same reason: a fixture
    that reads live evidence can only exercise the branch the source is
    currently on, so the contradiction branch became untestable the moment the
    repair landed. Both directions are named here so both stay checkable.
    """
    defective = copy.deepcopy(evidence)
    defective["free_water_in_the_wetland_infiltration_path"] = {
        "rain_melt_clamped_to_zero": True,
        "full_deficit_added_regardless": True,
        "switch_only_records_the_created_water": True,
        "defect_present": True}
    defective["constant_wetland_runon"] = True
    defective["annual_water_table_guard"] = True
    defective["day_resets_to_zero_at_year_end"] = True
    defective["serialized_soil_members"] = sorted(
        set(evidence.get("serialized_soil_members") or [])
        - {"Wtot", "wtd", "stand_water", "mwtp", "Frac_ice", "rootfrac"})
    return defective


def _unserialized(evidence: dict, member: str) -> dict:
    """The evidence a serializer missing exactly one required member gives."""
    without = copy.deepcopy(evidence)
    without["serialized_soil_members"] = sorted(
        set(evidence.get("serialized_soil_members") or []) - {member})
    return without


def _repaired(evidence: dict, retained: list[str] | None = None) -> dict:
    """The evidence a repaired vendored source would produce."""
    repaired = copy.deepcopy(evidence)
    repaired["free_water_in_the_wetland_infiltration_path"] = {
        "defect_present": False}
    repaired["constant_wetland_runon"] = False
    repaired["annual_water_table_guard"] = False
    repaired["latitude_selects_the_wetland_regime"] = False
    repaired["serialized_soil_members"] = sorted(
        set(evidence.get("serialized_soil_members") or [])
        | {"Wtot", "wtd", "stand_water", "mwtp", "Frac_ice", "rootfrac"})
    repaired["declared_output_parameters"] = sorted(
        set(evidence.get("declared_output_parameters") or [])
        | {f"file_{name.split('.')[0]}" for name in (retained or [])})
    return repaired


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-run", type=Path, default=None,
                        help="accept or reject one run's wetland output")
    parser.add_argument("--declaration", type=Path, default=DECLARATION)
    args = parser.parse_args()

    declaration = read_declaration(args.declaration)
    planet = yaml.safe_load(CONFIG.read_text())
    evidence = probes()
    refusals = evaluate(declaration, planet, evidence)
    requested = bool(declaration.get("requested"))
    fixtures = _fixtures(declaration, planet, evidence)

    run_refusals: list[Refusal] = []
    if args.check_run is not None:
        run_refusals = check_run(args.check_run, declaration)

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "biosphere/scripts/wetland_gate.py",
        "declaration": str(rel(args.declaration)),
        "declaration_sha256": hashlib.sha256(
            args.declaration.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "requested": requested,
        "granted": requested and not refusals,
        "switches_written": switches(requested and not refusals),
        "residual_tolerance_ceiling": RESIDUAL_TOLERANCE_CEILING,
        "months_in_the_simulated_year": _simulated_months(),
        "refusals": [r.as_dict() for r in refusals],
        "run_checked": str(args.check_run) if args.check_run else None,
        "run_refusals": [r.as_dict() for r in run_refusals],
        "fixtures": fixtures,
        "evidence": evidence,
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    state = "REQUESTED" if requested else "not requested"
    print(f"wetland, peat and methane activation {state}; "
          f"{len(refusals)} preconditions unmet")
    for refusal in refusals:
        print(f"  {refusal}")
    broken = [f for f in fixtures if not f["pass"]]
    print(f"\nfixtures: {len(fixtures) - len(broken)} of {len(fixtures)} got "
          "the verdict they were built for")
    for case in broken:
        print(f"  BROKEN: {case['fixture']}: expected {case['expected']}")
    print(f"\nwrote instruction settings {switches(False)} while refused")
    if args.check_run is not None:
        print(f"\n{rel(args.check_run)}: {len(run_refusals)} acceptance "
              "failures")
        for refusal in run_refusals:
            print(f"  {refusal}")
    print(f"\nwrote {rel(REPORT)}")

    if broken:
        print("\nA fixture did not get the verdict it was built for. That is a "
              "defect in this checker, not in the declaration.")
        sys.exit(1)
    if requested and refusals:
        sys.exit(1)
    if run_refusals:
        sys.exit(1)


if __name__ == "__main__":
    main()
