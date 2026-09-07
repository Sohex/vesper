#!/usr/bin/env python3
"""The ocean tier gate: what the adopted ocean is declared to be, and whether the model still agrees.

WORLDBUILDING CONTEXT, stated first because this file borrows vocabulary from a
real discipline: Vesper is an invented planet and everything below is about the
simulation of it -- a climate model's slab ocean, its sea ice, the sea surface
its atmosphere sees, and the file that declares them. Nothing here is Earth
science.

WHICH TIER THIS GUARDS. The ocean that runs is ExoPlaSim's slab.
`oceanmod.f90` integrates one mixed layer, `icemod.f90` runs the sea-ice
thermodynamics and owns the sea water, `seamod.f90` is the surface the
atmosphere sees. `vendor/cgenie` is a CANDIDATE under OCN-3 that nothing reads
and that does not build where it stands, so nothing here looks at it: a lint on
a configuration no step in `config/pipeline.yaml` produces would be built before
anyone had asked whether it should exist.

WHAT IT IS FOR. `notes/audits/ocean-tier-implicit-earth.md` asked, of the
candidate tiers, where an Earth number reaches this world with nothing saying it
was chosen. Its headline was two independent statements of one quantity that
nothing compared, and its preconditions were a check that they agree, a
fail-open closed, and a column saying of every constant whether a run can reach
it. `exoplasim/config/ocean_tier.yaml` is that question answered for the tier
that actually runs, and this module is the enforcement.

It can fail:

  duplicate    a quantity two modules act on that is no longer stated once and
               handed on: an owner's declaration gone, a receiver's assignment
               gone, a handoff gone, or a second compile-time copy back beside
               the handoff. All four are checked, because a handoff whose
               assignment has gone leaves the receiver on its own default and a
               second copy beside a live handoff is the split re-created
  drift        a constant whose declaration line the source no longer runs, or
               whose value is not what that line's own initialiser evaluates to
  reach        a constant whose declared reachability the source contradicts: one
               declared `namelist` that is not in its module's own namelist
               block, or one declared `compile_time` or `handed` that is. This is
               the column the audit added, and it decides how a constant fails:
               a namelist key can be moved by a run, a compiled one only by a
               source edit
  range        the sea-ice albedo ramp leaving the hull of its own endpoints,
               swept over the whole driver domain and past both ends of it, at
               every corner of the endpoint bracket a run can produce. The
               fractional form cannot leave it; the absolute-slope form the model
               used to run did, by 0.004 at the melting point
  invariant    a line something written down elsewhere rests on that the model no
               longer runs
  refusal      a fail-closed guard that has gone. Each is a refusal the model
               makes rather than running on an Earth number or an unset field

Reduced fixtures run on every invocation, all but one built to be wrong in a
named way. A fixture that does not get the verdict it was built for is a defect
in this checker rather than in the declaration.

    python exoplasim/scripts/ocean_tier_gate.py            # status, exit 0
    python exoplasim/scripts/ocean_tier_gate.py --strict   # refuses on the
                                                           # named residual

`--strict` refuses on exactly one thing: a constant of the modelled sea or its
ice that is declared `sourced: false`, meaning an Earth measurement standing
where nothing says it was chosen for this world. That set is not empty and is
not meant to be: it is the residue the audit's inventory left, and each entry
has to be settled before the tier's numbers are Vesper's rather than inherited.

The default arm reports and exits 0. A run on the declared tier is a correct run
of a declared model boundary, and the boundary is what this file states.

The Fortran reader is `lib/fortran_source.py`, shared with
`config/cgenie_calibration.yaml`'s gate so that the two declarations cannot
disagree about what a source line says or about which symbols a run can
reach.

NOTHING HERE IS VERIFIED BY EXECUTION, and that is stated rather than left to be
inferred. Every statement is against the source and the declaration. What the
slab DELIVERS is verified separately and by running the model, in
`exoplasim/scripts/verify_ocean_flux_channel.py`.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import COMPONENT_ROOT, PROJECT_ROOT  # noqa: F401  (also puts lib/ on sys.path)

from fortran_source import (  # noqa: E402  from lib/, via _paths
    initialiser, namelist_blocks, squash, tight)
from paths import rel  # noqa: E402  from lib/, via _paths

DECLARATION = COMPONENT_ROOT / "config" / "ocean_tier.yaml"
REPORT = COMPONENT_ROOT / "analysis" / "ocean_tier_gate_report.json"

# How finely the albedo ramp's drivers are swept. Both ramps are piecewise
# linear in the temperature and exactly linear in the ice fraction, so a maximum
# found on a grid this dense is the maximum exactly; the density is there so the
# sweep lands on both breakpoints and inside every branch.
TEMPERATURE_SAMPLES = 601
FRACTION_SAMPLES = 51

# How far past each end of the ramp the temperature sweep runs, in units of the
# ramp width. The clipped branches are where a mis-declared form shows up, and
# the driver reaches them: the modelled sea surface goes well below the melting
# point under ice and above it at the margin.
OVERRUN_WIDTHS = 3.0

# What a constant's reachability may say. The audit's added column.
REACHES = ("namelist", "handed", "compile_time")
# Where a number came from.
ORIGINS = ("earth", "declared", "derived")
# Which ramp forms this model has had. `fractional` is what it runs.
RAMP_FORMS = ("fractional", "absolute_slope")

# Floating-point slack on the range statement. The hull check is exact
# arithmetic on a convex combination, so this only absorbs binary rounding.
RANGE_EPS = 1e-12


# ---------------------------------------------------------------------------
# The two ramps. Reproduced rather than approximated: the range statement bounds
# these functions, and a paraphrase of them would bound something else.
# ---------------------------------------------------------------------------

def ramp_fraction_fractional(driver_k: float, warm_anchor_k: float,
                             width_k: float) -> float:
    """seamod's `zsicf`: AMAX1(0.0, AMIN1(1.0, (tmelt - dts)/dicealbdt))."""
    if width_k == 0.0:
        return 1.0 if driver_k < warm_anchor_k else 0.0
    return max(0.0, min(1.0, (warm_anchor_k - driver_k) / width_k))


def sea_albedo(fraction_of_ice: float, driver_k: float, open_water: float,
               ice_min: float, ice_max: float, *, form: str,
               warm_anchor_k: float, width_k: float,
               anchor_k: float = 0.0, slope_per_k: float = 0.0) -> float:
    """The albedo seastep writes into one band, for one cell.

    `fractional` is what seamod runs: the ice albedo is the endpoint pair
    blended by a FRACTION of the ramp, and the cell is that blended against the
    open-water value by the ice cover. Two nested convex combinations, so the
    result lies in the hull of the three endpoints whatever the drivers do.

    `absolute_slope` is what upstream ran and is here so the range statement has
    something that can fail. It ramps the albedo by `slope_per_k` per K below
    `anchor_k` and clips at the maximum, which reaches the declared minimum only
    where the span is exactly slope times width; anywhere else it either
    saturates early or, at the warm end, returns a value BELOW the minimum it
    declares.
    """
    if form == "fractional":
        z = ramp_fraction_fractional(driver_k, warm_anchor_k, width_k)
        ice = ice_min + z * (ice_max - ice_min)
    elif form == "absolute_slope":
        ice = min(ice_max, ice_min + slope_per_k * (anchor_k - driver_k))
    else:
        raise ValueError(form)
    return open_water * (1.0 - fraction_of_ice) + fraction_of_ice * ice


def _linspace(lo: float, hi: float, n: int) -> list[float]:
    if n < 2:
        return [lo]
    step = (hi - lo) / (n - 1)
    return [lo + step * i for i in range(n)]


def _corners(bracket: dict) -> list[tuple[float, float, float]]:
    """Every corner of the endpoint bracket.

    The endpoints are not the compiled ones: radmod re-derives the ocean and
    sea-ice albedos from the stellar spectrum at init, so the statement being
    checked has to hold for every endpoint set a run can produce, and the
    corners of the declared bracket are what bound that.
    """
    out = []
    for w in bracket["open_water"]:
        for mn in bracket["ice_min"]:
            for mx in bracket["ice_max"]:
                out.append((float(w), float(mn), float(mx)))
    return out


def ramp_call(entry: dict) -> str:
    """The statement the source has to contain, built from the declaration.

    Built and not declared, so an expression changed in the source and one
    changed in the declaration fail the same check. Nothing here can move alone.
    """
    return (f"{entry['fraction']}(:)=AMAX1(0.0,AMIN1(1.0,"
            f"({entry['warm_anchor']}-{entry['driver']}(:))/{entry['width']}))")


# ---------------------------------------------------------------------------
# The checks.
# ---------------------------------------------------------------------------

def check(declaration: dict, sources: dict[str, str]) -> list[dict]:
    """Every check, as a list of findings. An empty list is a clean gate."""
    findings: list[dict] = []

    def bad(kind: str, what: str, detail: str) -> None:
        findings.append({"kind": kind, "what": what, "detail": detail})

    def source_of(entry: dict, kind: str, what: str) -> str | None:
        name = entry.get("source")
        if name not in sources:
            bad(kind, what, f"names source {name!r}, which this gate does not read")
            return None
        return sources[name]

    # -- one quantity, one declaration -------------------------------------
    for entry in declaration.get("single_declaration", []):
        what = entry.get("id", "?")
        if not entry.get("why"):
            bad("duplicate", what, "a handoff saying nothing about what splitting it would cost")
        owner = entry.get("owner") or {}
        text = source_of(owner, "duplicate", what)
        if text is not None:
            for field in ("declaration", "via"):
                form = owner.get(field)
                if form and tight(form) not in tight(text):
                    bad("duplicate", what,
                        f"the owner's {field} {form!r} is not in {owner['source']}, so "
                        "the quantity has no single statement to hand on")
            if not (owner.get("declaration") or owner.get("via")):
                bad("duplicate", what, "an owner that neither declares the quantity nor names where it comes from")
        receivers = entry.get("receivers") or []
        if not receivers:
            bad("duplicate", what, "a handoff with no receiver, which is not a handoff")
        for receiver in receivers:
            text = source_of(receiver, "duplicate", what)
            if text is None:
                continue
            flat = tight(text)
            everywhere = "".join(tight(t) for t in sources.values())
            for field in ("declaration", "assignment"):
                form = receiver.get(field)
                if not form:
                    bad("duplicate", what, f"a receiver naming no {field}")
                elif tight(form) not in flat:
                    bad("duplicate", what,
                        f"{receiver['source']} no longer runs the receiver's {field} "
                        f"{form!r}, so it holds its own value and nothing compares them")
            # The handoff is a CALL, and the call site is in the module that
            # hands the value over rather than the one that takes it, so it is
            # looked for across the tier and not in the receiver alone.
            form = receiver.get("handoff")
            if not form:
                bad("duplicate", what, "a receiver naming no handoff")
            elif tight(form) not in everywhere:
                bad("duplicate", what,
                    f"nothing in the tier still runs {form!r}, so the receiver is "
                    "never handed the value and keeps its own")
        forbidden = entry.get("forbidden") or []
        if not forbidden:
            bad("duplicate", what,
                "a handoff naming no form that must stay out of the source. Without "
                "one, a second compile-time copy beside the handoff would pass")
        for item in forbidden:
            text = source_of(item, "duplicate", what)
            if text is None:
                continue
            if not item.get("why"):
                bad("duplicate", what, f"a forbidden form {item.get('form')!r} saying nothing about what it would cost")
            if tight(item["form"]) in tight(text):
                bad("duplicate", what,
                    f"{item['source']} carries {item['form']!r} again, which is a "
                    "second statement of the quantity beside the handoff")

    # -- the constants ------------------------------------------------------
    block = declaration.get("constants", {})
    reach_by_name: dict[str, str] = {}
    for item in block.get("values", []):
        name = item["name"]
        what = f"{item.get('source', '?')}:{name}"
        text = source_of(item, "drift", what)
        if text is None:
            continue
        line = item["declaration"]
        if tight(line) not in tight(text):
            bad("drift", what, f"declared as {line!r}, and that source does not run that line")
        derived = initialiser(line)
        if derived is None:
            bad("drift", what, f"the declared line {line!r} has no readable initialiser")
        elif abs(derived - float(item["value"])) > 1e-9 * max(1.0, abs(derived)):
            bad("drift", what,
                f"declared {item['value']}, and its own line evaluates to {derived:.10g}")
        if item.get("origin") not in ORIGINS:
            bad("drift", what, f"unknown origin {item.get('origin')!r}")
        if not item.get("bites"):
            bad("drift", what, "a constant saying nothing about what it changes")
        if item.get("sourced") and not item.get("source_of_value"):
            bad("drift", what, "declared sourced and naming no source")

        reach = item.get("reach")
        reach_by_name[what] = reach
        if reach not in REACHES:
            bad("reach", what,
                f"unknown reach {reach!r}. A constant has to say whether a run can "
                "move it, because that is what decides how it fails")
            continue
        blocks = namelist_blocks(text)
        everywhere = set().union(*blocks.values()) if blocks else set()
        lowered = name.lower()
        if reach == "namelist":
            wanted = item.get("namelist")
            if not wanted:
                bad("reach", what, "declared namelist-reachable and naming no namelist")
            elif lowered not in blocks.get(wanted, set()):
                bad("reach", what,
                    f"declared reachable through {wanted}, and that block does not "
                    "carry it. A run that set it would be refused and a run that "
                    "did not would take the compiled value with nothing saying so")
        elif lowered in everywhere:
            bad("reach", what,
                f"declared {reach}, and its own source has it in a namelist block. "
                "Its reachability is what the audit's inventory turns on")

    # -- the sea-ice albedo ramp -------------------------------------------
    for name, entry in (declaration.get("ramps") or {}).items():
        text = source_of(entry, "drift", name)
        if text is None:
            continue
        form = entry.get("form")
        if form not in RAMP_FORMS:
            bad("drift", name, f"unknown ramp form {form!r}")
            continue
        if form == "fractional":
            call = ramp_call(entry)
            found = tight(text).count(tight(call))
            wanted = int(entry["occurrences"])
            if found != wanted:
                bad("drift", name,
                    f"declared as {call!r} appearing {wanted} time(s), and "
                    f"{entry['source']} runs it {found}. seaini and seastep each "
                    "write it, and one moving without the other is a cold start "
                    "that disagrees with every step after it")
        for item in entry.get("forbidden_forms", []) or []:
            if not item.get("why"):
                bad("drift", name, f"a forbidden form {item.get('form')!r} saying nothing about what it would mean")
            if tight(item["form"]) in tight(text):
                bad("drift", name,
                    f"{entry['source']} carries {item['form']!r}, which belongs to "
                    "the superseded ramp")

        superseded = entry.get("superseded") or {}
        anchor = float(superseded.get("anchor", 0.0))
        slope = float(superseded.get("slope", 0.0))
        # The endpoints of the ramp itself come from the declared constants, so
        # the sweep runs on the widths the model actually holds.
        width = None
        warm = None
        for item in block.get("values", []):
            if item["name"] == entry["width"]:
                width = float(item["value"])
            if item["name"] == entry["warm_anchor"]:
                warm = float(item["value"])
        if width is None or warm is None:
            bad("range", name,
                f"ramps on {entry['warm_anchor']} and {entry['width']}, and at least "
                "one of those is not a declared constant")
            continue
        lo = warm - OVERRUN_WIDTHS * abs(width)
        hi = warm + abs(width)
        worst: tuple[float, float, tuple] | None = None
        for open_water, ice_min, ice_max in _corners(entry["endpoint_bracket"]):
            allowed_lo = min(open_water, ice_min, ice_max)
            allowed_hi = max(open_water, ice_min, ice_max)
            for driver in _linspace(lo, hi, TEMPERATURE_SAMPLES):
                for cover in _linspace(0.0, 1.0, FRACTION_SAMPLES):
                    a = sea_albedo(cover, driver, open_water, ice_min, ice_max,
                                   form=form, warm_anchor_k=warm, width_k=width,
                                   anchor_k=anchor, slope_per_k=slope)
                    if a < allowed_lo - RANGE_EPS or a > allowed_hi + RANGE_EPS:
                        excess = max(allowed_lo - a, a - allowed_hi)
                        if worst is None or excess > worst[0]:
                            worst = (excess, a, (open_water, ice_min, ice_max,
                                                 driver, cover))
        if worst is not None:
            excess, value, (w, mn, mx, driver, cover) = worst
            bad("range", name,
                f"reaches {value:.6g} at open water {w}, ice {mn} to {mx}, driver "
                f"{driver:.4g} K, cover {cover:.3g} -- outside the [{min(w, mn, mx):.6g}, "
                f"{max(w, mn, mx):.6g}] its own three endpoints allow, by {excess:.6g}. "
                "The blend is two nested convex combinations and cannot leave that hull")

    # -- lines the argument rests on ---------------------------------------
    for entry in declaration.get("invariants", []) or []:
        what = entry.get("name", "?")
        text = source_of(entry, "invariant", what)
        if text is None:
            continue
        if not entry.get("why"):
            bad("invariant", what, "an invariant saying nothing about what rests on it")
        if tight(entry["line"]) not in tight(text):
            bad("invariant", what,
                f"rests on {entry['line']!r}, and {entry['source']} does not run it")

    # -- the fail-closed guards --------------------------------------------
    for entry in declaration.get("refusals", []) or []:
        what = entry.get("name", "?")
        text = source_of(entry, "refusal", what)
        if text is None:
            continue
        if not entry.get("guards"):
            bad("refusal", what, "a guard saying nothing about what it refuses to run on")
        if squash(entry["form"]) not in squash(text):
            bad("refusal", what,
                f"{entry['source']} no longer refuses on {entry['form']!r}. Without "
                "it the run proceeds on a compiled default or an unset field")
    return findings


# ---------------------------------------------------------------------------
# The fixtures. The first is the declaration as it stands and has to come back
# clean; every other is built to be wrong in a named way. A fixture that does
# not get its verdict is a defect in the checker.
# ---------------------------------------------------------------------------

def _fixtures(declaration: dict, sources: dict[str, str]) -> list[dict]:

    def mutate(fn):
        d = copy.deepcopy(declaration)
        fn(d)
        return d

    def _constant(d, name):
        for item in d["constants"]["values"]:
            if item["name"] == name:
                return item
        raise KeyError(name)

    def constant_claim(name, field, value):
        return lambda d: _constant(d, name).__setitem__(field, value)

    def _handoff(d, entry_id):
        for entry in d["single_declaration"]:
            if entry["id"] == entry_id:
                return entry
        raise KeyError(entry_id)

    def handoff_claim(entry_id, field, value):
        return lambda d: _handoff(d, entry_id).__setitem__(field, value)

    def ramp_claim(field, value):
        return lambda d: d["ramps"]["sea_ice_albedo"].__setitem__(field, value)

    def invariant_claim(name, field, value):
        def apply(d):
            for entry in d["invariants"]:
                if entry["name"] == name:
                    entry[field] = value
                    return
            raise KeyError(name)
        return apply

    def refusal_claim(name, field, value):
        def apply(d):
            for entry in d["refusals"]:
                if entry["name"] == name:
                    entry[field] = value
                    return
            raise KeyError(name)
        return apply

    cases = [
        ("the declaration as it stands", declaration, None),

        ("a handoff whose receiver no longer takes the value",
         mutate(handoff_claim("melting_point", "receivers", [
             {"source": "ice",
              "declaration": "      real :: tmelt = 273.16            ! melting temp. for snow (0 deg C)",
              "assignment": "tmelt = no_such_argument",
              "handoff": "call iceini"}])),
         "duplicate"),
        ("a handoff whose owner no longer states the quantity",
         mutate(handoff_claim("sea_ice_density", "owner",
                              {"source": "ice",
                               "declaration": "      parameter(CRHOI = 917.)"})),
         "duplicate"),
        # The form here follows icemod's OWNER LINE, which world-04ok changed
        # from `parameter(CRHOI = 920.)` to a namelist-reachable `real ::`. The
        # fixture tests that the owner's own line, wrongly listed as forbidden,
        # is caught -- so it has to name whatever that line currently is, or it
        # matches nothing and tests nothing.
        ("a second compile-time copy back beside a live handoff",
         mutate(handoff_claim("sea_ice_density", "forbidden",
                              [{"source": "ice", "form": "real :: CRHOI",
                                "why": "the owner's own line, which is not a second copy"}])),
         "duplicate"),
        ("a handoff naming no forbidden form",
         mutate(handoff_claim("sea_water", "forbidden", [])),
         "duplicate"),

        ("a constant whose declaration line the source does not run",
         mutate(constant_claim("CPS", "declaration",
                               "      real :: CPS       = 4180.0   ! specific heat of sea water (J/(kg*K))")),
         "drift"),
        ("a constant whose value is not what its own line evaluates to",
         mutate(constant_claim("CLFSN", "value", 3.337)),
         "drift"),
        ("a constant saying nothing about what it changes",
         mutate(constant_claim("CKAPI", "bites", "")),
         "drift"),
        # `crhosn` rather than CRHOI, which world-04ok made a real icemod_nl key:
        # a fixture claiming namelist reach has to name a constant the namelist
        # genuinely does not carry, and crhosn is handed from landmod.
        ("a constant declared reachable through a namelist that does not carry it",
         mutate(constant_claim("crhosn", "reach", "namelist")),
         "reach"),
        ("a namelist key declared compile-time",
         mutate(constant_claim("hlead", "reach", "compile_time")),
         "reach"),
        ("a constant declared reachable and naming no namelist",
         mutate(constant_claim("thicec", "namelist", None)),
         "reach"),

        ("the superseded absolute-slope ramp declared as the one that runs",
         mutate(ramp_claim("form", "absolute_slope")),
         "range"),
        ("a ramp whose expression the source does not contain",
         mutate(ramp_claim("fraction", "zsicg")),
         "drift"),
        ("a ramp called from fewer places than declared",
         mutate(ramp_claim("occurrences", 3)),
         "drift"),
        ("a ramp whose forbidden form is in fact what the source runs",
         mutate(ramp_claim("forbidden_forms",
                           [{"form": "AMAX1(0.0,AMIN1(1.0,",
                             "why": "the live ramp, which is not a superseded form"}])),
         "drift"),

        ("an invariant the source no longer runs",
         mutate(invariant_claim("ocean_timestep_is_the_planets_day", "line",
                                "dtmix = 86400.0 / real(ntspd)")),
         "invariant"),
        ("an invariant saying nothing about what rests on it",
         mutate(invariant_claim("ocean_radius_is_the_planets", "why", "")),
         "invariant"),

        ("a fail-closed guard the source no longer makes",
         mutate(refusal_claim("cold_start_must_be_declared", "form",
                              "icemod: a cold start needs no declaration")),
         "refusal"),
        ("a guard saying nothing about what it refuses to run on",
         mutate(refusal_claim("ocean_diffusion_needs_a_radius", "guards", "")),
         "refusal"),
    ]

    results = []
    for label, candidate, expect in cases:
        findings = check(candidate, sources)
        kinds = {f["kind"] for f in findings}
        ok = (not findings) if expect is None else (expect in kinds)
        results.append({"fixture": label, "expected": expect or "clean",
                        "found": sorted(kinds), "pass": ok})
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--strict", action="store_true",
                        help="refuse while an Earth constant of the tier has no source")
    parser.add_argument("--json", action="store_true", help="report only")
    args = parser.parse_args()

    declaration = yaml.safe_load(DECLARATION.read_text(encoding="utf-8"))
    sources = {name: (PROJECT_ROOT / path).read_text(encoding="utf-8")
               for name, path in declaration["sources"].items()}

    fixtures = _fixtures(declaration, sources)
    findings = check(declaration, sources)

    unsourced = [{"name": item["name"], "source": item.get("source", "?"),
                  "reach": item.get("reach", "?"), "bites": item.get("bites", "")}
                 for item in declaration["constants"]["values"]
                 if item.get("origin") == "earth" and not item.get("sourced")]
    by_reach: dict[str, int] = {}
    for item in declaration["constants"]["values"]:
        by_reach[item.get("reach", "?")] = by_reach.get(item.get("reach", "?"), 0) + 1

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": rel(DECLARATION),
        "tier": "ExoPlaSim slab; vendor/cgenie is a candidate nothing reads",
        "sources": {k: rel(PROJECT_ROOT / v) for k, v in declaration["sources"].items()},
        "findings": findings,
        "constants_by_reach": by_reach,
        "handoffs": [e["id"] for e in declaration.get("single_declaration", [])],
        "refusals": [e["name"] for e in declaration.get("refusals", [])],
        "unsourced_earth_constants": unsourced,
        "execution_verified": False,
        "fixtures": fixtures,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    broken = [f for f in fixtures if not f["pass"]]
    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print("The adopted ocean tier, checked against the source the model compiles.\n")
        print(f"  fixtures: {len(fixtures) - len(broken)} of {len(fixtures)} got their verdict")
        for case in broken:
            print(f"    BROKEN: {case['fixture']}: expected {case['expected']}, found {case['found']}")
        if findings:
            print(f"\n  {len(findings)} finding(s):")
            for finding in findings:
                print(f"    [{finding['kind']}] {finding['what']}: {finding['detail']}")
        else:
            print("  every handoff stands, every constant matches its own line, the")
            print("  sea-ice albedo ramp holds the hull of its endpoints, and every")
            print("  fail-closed guard is still in the source")
        print(f"\n  {len(report['handoffs'])} quantity(ies) stated once and handed on:")
        for item in report["handoffs"]:
            print(f"    {item}")
        print(f"\n  reachability of {sum(by_reach.values())} declared constant(s): "
              + ", ".join(f"{n} {k}" for k, n in sorted(by_reach.items())))
        print(f"  {len(report['refusals'])} fail-closed guard(s) in the source")
        if unsourced:
            print(f"\n  {len(unsourced)} Earth constant(s) of the modelled sea with no")
            print("  source, which is what --strict refuses on:")
            for item in unsourced:
                print(f"    {item['name']:10s} [{item['reach']}]  {item['bites'].strip().splitlines()[0]}")
        print("\n  NOTHING HERE IS VERIFIED BY EXECUTION. Every statement is against the")
        print("  source and the declaration; what the slab DELIVERS is verified by")
        print("  running the model, in verify_ocean_flux_channel.py.")
        print(f"\n  report: {rel(REPORT)}")

    if broken:
        print(f"\n{len(broken)} fixture(s) did not get the verdict they were built for.",
              file=sys.stderr)
        print("That is a defect in this checker, not in the declaration.", file=sys.stderr)
        return 2
    if findings:
        return 1
    if args.strict and unsourced:
        print("\n--strict: refused, on exactly what is still unsourced and nothing else.",
              file=sys.stderr)
        print("  Earth constant(s) of the modelled sea with no source: "
              + ", ".join(item["name"] for item in unsourced), file=sys.stderr)
        print("Each has to be settled before this tier's numbers are Vesper's.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
