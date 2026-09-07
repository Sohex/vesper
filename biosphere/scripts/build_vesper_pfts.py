"""Apply the biosphere time-base unit contract to the PFT instruction file.

LPJ-GUESS's plant functional types are Earth's, and keeping them is a declared
choice: this is an Earth-analogue biosphere, not a prediction of alien
physiology. But their parameters are Earth calibrations, and several of them are
per *year* -- an Earth year of 365.2569 days, not this world's orbit of roughly
half that. This script is the boundary between the two: `global.ins` is in EARTH
units, and the file it writes is in MODEL units, where one year is one orbit.
Every parameter it carries is therefore in one of four classes, and the class
decides the conversion. The classes, the evidence and the reader-by-reader
inventory are in `biosphere/notes/time-base-unit-contract.md`.

**ANNUAL_SUM** -- a threshold on something accumulated over one year.
`gdd5min_est` is the minimum growing degree-days above 5 C a PFT needs in a year
to establish; `greff_min` is a threshold on annual net production per unit leaf
area. Degree-days and carbon accumulate per day at the same rate on both worlds,
so a shorter year reaches a smaller total for an identical climate, and the
threshold scales DOWN with the year. Left alone, Earth thresholds exclude nearly
every tree PFT for reasons that have nothing to do with the climate. That is not
hypothetical: the patched model run on Earth's own demo data collapses boreal
needleleaf and temperate broadleaf to grass.

**YEAR_COUNT** -- a duration counted in years. A simulation year is about half an
Earth year, so representing the same absolute span needs about twice as many of
them, and these scale UP by the reciprocal. `nyear_spinup 500` reads like an
absolute statement and is not.

**ANNUAL_RATE** -- a FRACTION applied once per year. A fraction does not scale
linearly: an Earth-calibrated `r` per Earth year becomes `1 - (1-r)**f` per
orbit, which leaves 1.0 at 1.0 where a linear multiplier would silently change
what "all of it" means.

**UNSCALED** -- physiology, within-season accumulations, sentinels, and anything
counted in seasonal cycles rather than absolute time. Named explicitly, because a
parameter that is merely absent from the lists reads as an oversight.

Confusing the directions would be worse than doing neither, so the class of every
parameter this script touches is written into the generated header and into the
provenance JSON beside it.

The scale factor is derived from the configured orbit, not written down, because
the year length moves with the stellar flux. It uses the MODELLED year -- whole
24-hour steps, which is what the model integrates over -- and not the true
orbital period, which the model cannot represent.

    python biosphere/scripts/build_vesper_pfts.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import CONFIG, GENERATED, GUESS_SOURCE, PROJECT_ROOT
from paths import rel  # noqa: E402

import lpj_pfts
import orbit
import run_lengths

EARTH_YEAR_DAYS = orbit.EARTH_SIDEREAL_YEAR_DAYS

# A threshold on a quantity accumulated over one year. Scales DOWN by the year
# ratio, linearly, because the accumulation itself is linear in time.
#
#   gdd5min_est  growing degree-days above 5 C needed in a year to establish
#   gdd5min      the same limit under its other spelling
#   gdd0_min     degree-day window on a 0 C base; sentinels in global.ins
#   gdd0_max
#   greff_min    kgC/m2 leaf per year below which growth suppression kills
ANNUAL_SUM = ("gdd5min_est", "gdd5min", "gdd0_min", "gdd0_max", "greff_min")

# A duration counted in years. Scales UP by the reciprocal.
#
# This is the same trap gdd5min was, one level up, and it is easy to miss because
# "500 years of spin-up" reads like an absolute statement and is not. At
# nyear_spinup 500 this world would get about 247 Earth years of soil and
# vegetation development where Earth practice assumes 500.
#
#   nyear_spinup   time for vegetation and soil pools to reach steady state
#   distinterval   mean return time of generic patch-destroying disturbance
#   freenyears     time allowed to build an N pool before N limitation bites
#   longevity      age at which 0.1% of a cohort survives; compared against
#                  Individual::age, which counts simulation years
#   leaflong       leaf lifespan; Pft::initsla converts it back to the absolute
#                  months Reich et al. (1992) regressed SLA and leaf C:N against
YEAR_COUNT = ("nyear_spinup", "distinterval", "freenyears", "longevity", "leaflong")

# A fraction of a pool moved once per year. Converted as a rate, not multiplied.
#
#   turnover_leaf  leaf C moved to litter per year
#   turnover_root  fine root C moved to litter per year
#   turnover_sap   sapwood C converted to heartwood per year
ANNUAL_RATE = ("turnover_leaf", "turnover_root", "turnover_sap")

# Named so the decision not to scale them is explicit and reviewable, rather
# than an omission someone later reads as an oversight.
DELIBERATELY_UNSCALED = (
    "phengdd5ramp",   # within-season accumulation, already absolute time
    "tcmin_surv", "tcmin_est", "tcmax_est", "twmin_est", "twminusc",
    "k_chilla", "k_chillb", "k_chillk",  # chill-day sums, within-season
    # Counted in growing seasons rather than in absolute time. One simulation
    # year is one seasonal cycle on this world just as on Earth, so these are
    # already in the right unit.
    "estinterval",
    "est_max",        # saplings per m2 per growing season, not per Earth year
)

# Declared as int by plib, so a fractional value would not parse. distinterval is
# NOT here: parameters.cpp declares it double, and rounding a disturbance return
# time to a whole number of orbits throws away precision for nothing.
WHOLE_NUMBER = ("nyear_spinup", "freenyears")

# Above this a value is a "no restriction" sentinel rather than a limit, and
# scaling it would quietly turn it into one.
SENTINEL_ABOVE = 1e4

# The plib bounds declared in vendor/lpj-guess/framework/parameters.cpp for the
# parameters this script rewrites. A conversion that leaves a value outside its
# own declared range is a generation-time failure rather than a run-time one:
# plib would reject the file, but only after a build and a launch, and only for
# the first offending line. A large enough change of orbit can do it -- longevity
# 500 leaves the 3000 ceiling below a sixth of an Earth year -- so this is a
# guard against a future flux, not a hypothetical.
DECLARED_BOUNDS = {
    "nyear_spinup": (1, 10000),
    "estinterval": (1, 10),
    "distinterval": (1.0, 1.0e10),
    "freenyears": (0, 1000),
    "phengdd5ramp": (0.0, 1000.0),
    "turnover_leaf": (0.0, 1.0),
    "turnover_root": (0.0, 1.0),
    "turnover_sap": (0.0, 1.0),
    "gdd5min_est": (0.0, 5000.0),
    "est_max": (1.0e-4, 1.0),
    "longevity": (0.0, 3000.0),
    "greff_min": (0.0, 1.0),
    "leaflong": (0.1, 100.0),
    "gdd0_min": (0.0, 100000.0),
    "gdd0_max": (0.0, 100000.0),
}


UNIT = {
    "ANNUAL_SUM": "per simulation year",
    "YEAR_COUNT": "simulation years",
    "ANNUAL_RATE": "fraction per simulation year",
}

CLASS_OF = {}
for _name in ANNUAL_SUM:
    CLASS_OF[_name] = "ANNUAL_SUM"
for _name in YEAR_COUNT:
    CLASS_OF[_name] = "YEAR_COUNT"
for _name in ANNUAL_RATE:
    CLASS_OF[_name] = "ANNUAL_RATE"


def convert(kind: str, value: float, factor: float) -> float:
    """The one place each class's arithmetic lives."""
    if kind == "ANNUAL_SUM":
        return value * factor
    if kind == "YEAR_COUNT":
        return value / factor
    if kind == "ANNUAL_RATE":
        # A fraction, so a rate conversion and not a multiplier. Leaves 1.0 at
        # 1.0, which is what a summergreen shedding every leaf every year means.
        return 1.0 - (1.0 - value) ** factor
    raise ValueError(kind)


def render(name: str, value: float) -> str:
    if name in WHOLE_NUMBER:
        return f"{value:.0f}"
    return f"{value:.6g}"


def self_check(factor: float, changes: list[dict]) -> None:
    """Checks with a right answer, run every time the artifact is generated.

    Each one can fail. That is the point: an identity, a declared range and a
    round trip, rather than a comparison against whatever the last run produced.
    """
    both = set(ANNUAL_SUM) | set(YEAR_COUNT) | set(ANNUAL_RATE)
    overlap = both & set(DELIBERATELY_UNSCALED)
    if overlap:
        raise SystemExit(
            f"{', '.join(sorted(overlap))} is both scaled and deliberately "
            f"unscaled. A parameter has exactly one class.")

    # The Reich identity. Pft::initsla and Pft::init_cton_min multiply leaflong
    # by VESPER_EARTH_MONTHS_PER_ORBIT, which is 12*factor, so rescaling leaflong
    # by 1/factor has to leave the product at the Earth value. If this fails, SLA
    # and leaf C:N have silently moved off their Earth calibration.
    earth_months_per_orbit = 12.0 * factor
    for change in changes:
        if change["parameter"] != "leaflong":
            continue
        earth = change["from"] * 12.0
        model = change["to"] * earth_months_per_orbit
        if abs(model - earth) > 1e-4 * earth:
            raise SystemExit(
                f"leaflong {change['from']} -> {change['to']} does not preserve "
                f"the Reich regression argument: {model:.6f} against "
                f"{earth:.6f} absolute months.")

    for change in changes:
        name = change["parameter"]
        lo, hi = DECLARED_BOUNDS[name]
        if not (lo <= change["to"] <= hi):
            raise SystemExit(
                f"{name} {change['from']} -> {change['to']} falls outside the "
                f"{lo} to {hi} range parameters.cpp declares for it. The orbit "
                f"has moved far enough that this parameter needs a decision, "
                f"not a rescale.")
        # The conversion itself, checked against the class table rather than
        # against whatever the last run produced.
        exact = convert(change["class"], change["from"], factor)
        if abs(exact - change["to"]) > 1e-6 * max(abs(exact), 1e-9):
            raise SystemExit(
                f"{name} recorded {change['to']} where its class gives {exact}.")
        # And the rendering, which is where precision is actually lost. A whole
        # number parameter can lose half a unit and nothing else may lose more
        # than the six significant figures render() writes.
        written = float(render(name, change["to"]))
        tolerance = 0.5 if name in WHOLE_NUMBER else 1e-5 * max(abs(exact), 1e-9)
        if abs(written - exact) > tolerance:
            raise SystemExit(
                f"{name} renders {exact} as {written}, losing more than "
                f"{tolerance}. Widen render().")


def _derived_spinup_cycles() -> dict:
    """The spin-up floor this world's ecology needs, and which derivation gave it.

    Reads `lib/run_lengths.py`, which reads the acceptance artifacts and states no
    number of its own. Returns the convention it would replace alongside it, so
    the provenance records a substitution rather than a bare value.

    THE CONVENTION IS THE FALLBACK OF LAST RESORT and is reached only where there
    is no acceptance artifact to read at all. Where there is one, a spin-up is
    always derived: at a measured relaxation time where one is admissible, and
    otherwise at the minimax over every relaxation time, which needs no
    measurement. `spinup_basis` names which, and it travels into the provenance.
    """
    convention = None
    try:
        source = (GUESS_SOURCE / "data" / "ins" / "global.ins").read_text()
        match = re.search(r"^[ \t]*\bnyear_spinup\b[ \t]+([0-9.]+)", source,
                          re.MULTILINE)
        if match:
            convention = float(match.group(1))
    except OSError:
        pass
    try:
        derived = run_lengths.ecological_run_cycles(PROJECT_ROOT)
    except RuntimeError as exc:
        return {"cycles": None, "source": "Earth convention, rescaled",
                "reason": str(exc), "convention_cycles": convention}
    # The RECORD floor and the SPIN-UP floor are independent and are recorded
    # independently, so a reader can see what the run has to RETAIN beside what
    # has to precede it.
    return {"cycles": int(math.ceil(derived["spinup_cycles"])),
            "source": "lib/run_lengths.py:ecological_run_cycles",
            "basis": derived["spinup_basis"],
            "is_minimax_over_relaxation_time": derived["spinup_is_minimax"],
            "reason": derived["spinup_reason"],
            "is_a_floor": derived["is_a_floor"],
            "floor_because": derived["floor_because"],
            "record_cycles": int(math.ceil(derived["record_cycles"])),
            "convention_cycles": convention,
            "brackets": derived["brackets"]}


NATIVE_PFTS = PROJECT_ROOT / "biosphere" / "config" / "native_pfts.yaml"


def _native_parameters(owner: str, spec: dict, arm: str | None):
    """One block's parameter lines, its resolved values and its bracket ends."""
    lines: list[str] = []
    values: dict[str, float] = {}
    ends: dict[str, str] = {}
    for key, entry in spec.get("parameters", {}).items():
        if "basis" not in entry:
            raise SystemExit(f"native {owner}.{key} states no basis")
        if "bracket" in entry:
            bracket = entry["bracket"]
            chosen = arm or bracket["default"]
            if chosen == "default" or chosen not in bracket:
                raise SystemExit(
                    f"native {owner}.{key} has no bracket end named {chosen!r}; "
                    f"it declares {sorted(k for k in bracket if k != 'default')}")
            value = float(bracket[chosen])
            ends[key] = chosen
            note = (f"NATIVE, BRACKETED: end {chosen!r} of "
                    f"{[bracket[k] for k in bracket if k != 'default']}")
        else:
            value = entry["value"]
            note = "NATIVE"
        kind = CLASS_OF.get(key)
        if kind and entry.get("units") != "model":
            raise SystemExit(
                f"native {owner}.{key} is a {kind} parameter, so a reader cannot "
                "tell a native value from an Earth value the conversion missed. "
                "Declare units: model to state that it is already in model units.")
        if kind:
            note += f"; {UNIT[kind]}, already in model units, not rescaled"
        text = f'"{value}"' if isinstance(value, str) else f"{value:g}"
        lines.append(f"\t{key} {text}\t! {note}")
        if not isinstance(value, str):
            values[key] = float(value)
    return lines, values, ends


def native_blocks(declared: dict, arm: str | None, source_text: str) -> tuple[str, list[dict]]:
    """Render the Vesper-native groups and types, and refuse the ways they go wrong.

    THESE ARE NOT RESCALED, and that is the whole reason they are rendered here
    rather than appended to the source before the conversion pass. A block
    declared in this file has no Earth calibration behind it; its values are in
    model units already. Passing them through the conversion would apply the
    year-length factor to a number that is already in the target unit.

    The hazard that creates is that `leaflong 0.2392` written natively is
    indistinguishable by eye from an Earth value the conversion missed. So any
    native parameter whose name falls in a conversion class must state
    `units: model`, and the emitted line carries that statement too.

    GROUPS COME FIRST AND EXIST FOR ONE REASON: two native types sharing a
    derived trait must not each state it. A shared value written twice is two
    declarations that drift, so the group holds what both derive and each type
    holds only what distinguishes it.
    """
    existing = set(re.findall(r'^\s*(?:pft|group)\s+"([\w.]+)"', source_text,
                              re.MULTILINE))
    rendered: list[str] = []
    recorded: list[dict] = []
    native_names: set[str] = set()

    for kind_name, blocks in (("group", declared.get("groups", {}) or {}),
                              ("pft", declared.get("types", {}) or {})):
        for name, spec in blocks.items():
            if name in existing:
                raise SystemExit(
                    f"native {kind_name} {name!r} collides with a type or group "
                    "already in the source instruction file; a native block must "
                    "not shadow one whose parameters came through the Earth "
                    "conversion")
            inherits = spec["inherits"]
            if inherits not in existing and inherits not in native_names:
                raise SystemExit(
                    f"native {kind_name} {name!r} inherits {inherits!r}, which is "
                    "neither in the source instruction file nor declared above it "
                    "here")
            native_names.add(name)
            head = [f'{kind_name} "{name}" (',
                    f"\t! {spec['title']} -- NATIVE to Vesper, invented biology.",
                    f"\t! Derivation: {spec['derivation']}; issue {spec['issue']}.",
                    "\t! Declared in biosphere/config/native_pfts.yaml and NOT",
                    "\t! rescaled: its values are in model units already.",
                    f"\t{inherits}"]
            if kind_name == "pft":
                head.append(
                    f"\tinclude {int(spec['include'])}"
                    f"\t! NATIVE; "
                    + ("reaches runs" if spec["include"] else
                       "declared but not instantiated"))
            body, values, ends = _native_parameters(name, spec, arm)
            rendered.append("\n".join(head + body + [")"]))
            recorded.append({"kind": kind_name, "name": name,
                             "inherits": inherits,
                             "include": int(spec["include"]) if kind_name == "pft" else None,
                             "derivation": spec["derivation"], "issue": spec["issue"],
                             "bracket_ends": ends, "arm_requested": arm,
                             "values": values,
                             "basis": {k: v["basis"]
                                       for k, v in spec.get("parameters", {}).items()}})
    if not rendered:
        return "", []
    banner = ("\n\n"
              "!///////////////////////////////////////////////////////////////////////////////\n"
              "!// NATIVE TO VESPER. Everything below is invented biology for this world,\n"
              "!// declared in biosphere/config/native_pfts.yaml and appended AFTER the Earth\n"
              "!// conversion, so none of it is rescaled: these values are in model units\n"
              "!// already. A block states only what its derivation changes and inherits the\n"
              "!// rest; a group holds what two types share, so no derived value is written\n"
              "!// twice.\n"
              "!///////////////////////////////////////////////////////////////////////////////\n\n")
    return banner + "\n\n".join(rendered) + "\n", recorded


def _selftest() -> int:
    """Exercise the native-block guards on fixtures, writing nothing.

    Each guard exists because the failure it catches is invisible in the
    generated file: a native value in a conversion class reads exactly like an
    Earth value the conversion missed, and a native type shadowing a shipped one
    reads exactly like the shipped one.
    """
    source = 'group "C3G" (\n\tgrass\n)\n\npft "C3G" (\n\tC3G\n)\n'

    def spec(**over):
        base = {"title": "t", "derivation": "d.md",
                "issue": "world-x", "include": 0, "inherits": "C3G",
                "parameters": {"pstemp_low": {"value": 15.0, "units": "absolute",
                                              "basis": "b"}}}
        name = over.pop("name", "VPE")
        base.update(over)
        return {"types": {name: base}}

    def refuses(declared, arm=None) -> bool:
        try:
            native_blocks(declared, arm, source)
        except SystemExit:
            return True
        return False

    checks: list[tuple[str, bool]] = []
    text, records = native_blocks(spec(), None, source)
    checks.append(("a well-formed native type renders and records",
                   'pft "VPE"' in text and records[0]["values"]["pstemp_low"] == 15.0))
    checks.append(("the block names the type it inherits from",
                   "\n\tC3G\n" in text))
    checks.append(("a native type shadowing a shipped one is refused",
                   refuses(spec(name="C3G"))))
    checks.append(("inheriting a type the source does not declare is refused",
                   refuses(spec(inherits="NOSUCH"))))
    checks.append(("a parameter with no basis is refused",
                   refuses(spec(parameters={"pstemp_low": {"value": 1.0,
                                                           "units": "absolute"}}))))
    # THE ONE THAT MATTERS. leaflong is YEAR_COUNT, so a native value sitting
    # beside a converted one is indistinguishable by eye.
    checks.append(("a conversion-class parameter not declared in model units is refused",
                   refuses(spec(parameters={"leaflong": {"value": 0.24,
                                                         "units": "absolute",
                                                         "basis": "b"}}))))
    ok_text, _ = native_blocks(
        spec(parameters={"leaflong": {"value": 0.24, "units": "model",
                                      "basis": "b"}}), None, source)
    checks.append(("the same parameter declared in model units states so on its line",
                   "already in model units, not rescaled" in ok_text))
    bracketed = spec(parameters={"ltor_max": {"bracket": {"low": 1.0, "high": 2.0,
                                                          "default": "low"},
                                              "units": "absolute", "basis": "b"}})
    low, low_rec = native_blocks(bracketed, None, source)
    high, high_rec = native_blocks(bracketed, "high", source)
    checks.append(("a bracket takes its declared default when no arm is named",
                   low_rec[0]["values"]["ltor_max"] == 1.0))
    checks.append(("a named arm overrides the default",
                   high_rec[0]["values"]["ltor_max"] == 2.0))
    checks.append(("both renderings say which end they took",
                   "end 'low'" in low and "end 'high'" in high))
    checks.append(("an arm the bracket does not declare is refused",
                   refuses(bracketed, "middle")))
    checks.append(("no type declared renders nothing rather than an empty banner",
                   native_blocks({"types": {}}, None, source) == ("", [])))

    for label, ok in checks:
        print(f"[{'ok' if ok else 'FAIL'}] {label}")
    failed = sum(1 for _, ok in checks if not ok)
    print(f"\n{failed} failures")
    return 1 if failed else 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path,
                        default=GUESS_SOURCE / "data" / "ins" / "global.ins")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--native-arm", default=None,
                        help="which end of a bracketed native parameter to emit; "
                             "the config's declared default is used when absent, "
                             "and whichever was taken is written into the header "
                             "and the provenance either way")
    parser.add_argument("--self-test", action="store_true",
                        help="exercise the native-block guards on fixtures and "
                             "write nothing")
    args = parser.parse_args()

    if args.self_test:
        raise SystemExit(_selftest())

    config = yaml.safe_load(CONFIG.read_text())
    orbital_days = orbit.orbital_year_days(config)
    model_year_days = orbit.model_year_days(config)
    # The MODELLED year, because that is the span the model integrates a degree
    # day sum over and counts a simulation year as. It differs from the true
    # orbit by the rounding to whole days.
    factor = orbit.earth_years_per_model_year(config)

    text = args.source.read_text()
    changes = []

    # Longest first so no name can match as the prefix of another.
    names = sorted(CLASS_OF, key=len, reverse=True)
    pattern = re.compile(
        r"^([ \t]*)\b(" + "|".join(names) + r")\b([ \t]+)([0-9.]+)[^\n]*$",
        re.MULTILINE)

    def rescale(match: re.Match) -> str:
        indent, name, gap, raw = match.groups()
        value = float(raw)
        kind = CLASS_OF[name]
        # Sentinels meaning "no limit" stay sentinels: scaling 0 is a no-op but
        # scaling a large "no restriction" value would quietly become a limit.
        if value <= 0.0 or value >= SENTINEL_ABOVE:
            return match.group(0)
        new = convert(kind, value, factor)
        # The exact converted value, not the rendered one: the provenance has to
        # answer "is this artifact still current?" against the conversion, and
        # the rendering is checked separately in self_check.
        changes.append({"parameter": name, "class": kind,
                        "from": value, "to": new,
                        "written": float(render(name, new))})
        return (f"{indent}{name}{gap}{render(name, new)}"
                f"\t! {UNIT[kind]}; Earth calibration {value:g}")

    rescaled = pattern.sub(rescale, text)
    self_check(factor, changes)

    # nyear_spinup IS NOT AN EARTH DURATION TO RESCALE, once anything has
    # measured this world's own relaxation time. The rescale above converts the
    # shipped 500 to 998 simulation years, which is Earth's CONVENTION carried
    # across a correct conversion -- and this world's woody types and slow pools
    # relax on hundreds of cycles, so the convention is short. The derived floor
    # is already in simulation years and must NOT be scaled again; scaling it
    # would be the exact YEAR_COUNT trap this file exists to avoid, applied to a
    # number that is not in the unit the trap assumes.
    #
    # BEST AVAILABLE, NOT FIRST AVAILABLE. On a tree with no assessed run there
    # is no measured relaxation time and the rescaled convention is genuinely the
    # best there is; once a run has been assessed the floor is, and pinning this
    # to the convention afterwards would converge to the wrong place quietly.
    # Which one was used is recorded rather than inferred.
    spinup = _derived_spinup_cycles()
    if spinup.get("cycles"):
        rescaled, applied = re.subn(
            r"^([ \t]*)\bnyear_spinup\b([ \t]+)([0-9.]+)[^\n]*$",
            lambda m: (f"{m.group(1)}nyear_spinup{m.group(2)}{spinup['cycles']}"
                       "\t! simulation years; DERIVED floor from the acceptance "
                       "contract, not Earth's convention"),
            rescaled, flags=re.MULTILINE)
        spinup["applied"] = bool(applied)
        if applied:
            # The rescale pass already recorded a YEAR_COUNT conversion for this
            # parameter and the override supersedes it. Leaving both would file
            # the value under a class it is no longer in, and the header's class
            # list is read as the statement of what happened to each parameter.
            changes[:] = [c for c in changes if c["parameter"] != "nyear_spinup"]
            changes.append({"parameter": "nyear_spinup", "class": "DERIVED",
                            "from": float(spinup["convention_cycles"]),
                            "to": float(spinup["cycles"]),
                            "written": float(spinup["cycles"])})

    # THE RESPIRATION PATH IS NAMED, NOT INHERITED. `global.ins` ships
    # `acclimated_respiration 1` and `global_p.ins` ships 0, so which one this
    # file carries was decided by which source it happened to import. Nothing
    # ran on the wrong path -- `run_lpj_guess.py` writes 0 into the run
    # instruction and `plib` takes the later declaration -- but an inherited 1
    # that happens to be overridden is not a stated choice, and a reader of this
    # artifact could not tell the difference. WORLD-VNBC.
    #
    # The DECISION is `biosphere/config/respiration_acclimation.yaml`'s
    # `path.runs`, taken under world-uhfh, and `acclimation_gate.py` already
    # holds `run_lpj_guess.py` to the same line. Reading it here puts the
    # generated file on that one declaration too, so the three cannot disagree.
    declared = yaml.safe_load(
        (PROJECT_ROOT / "biosphere" / "config" / "respiration_acclimation.yaml")
        .read_text(encoding="utf-8"))
    runs_acclimated = bool(declared["path"]["runs"])
    rescaled, applied_resp = re.subn(
        r"^([ \t]*)\bacclimated_respiration\b([ \t]+)([0-9]+)[^\n]*$",
        lambda m: (f"{m.group(1)}acclimated_respiration{m.group(2)}"
                   f"{int(runs_acclimated)}"
                   "\t! NAMED from biosphere/config/respiration_acclimation.yaml"
                   " path.runs, not inherited from the source instruction file"),
        rescaled, flags=re.MULTILINE)
    if not applied_resp:
        raise SystemExit(
            "the source instruction file declares no acclimated_respiration, so "
            "the respiration path cannot be named here. "
            "biosphere/config/respiration_acclimation.yaml expects to govern it; "
            "world-vnbc.")

    # THE NATIVE TYPES GO ON AFTER EVERYTHING THE CONVERSION TOUCHES. They carry
    # no Earth calibration, so their values are in model units already and the
    # conversion must not reach them. Appending rather than merging is what makes
    # that structural instead of remembered.
    native_declared = yaml.safe_load(NATIVE_PFTS.read_text(encoding="utf-8"))
    native_text, native_records = native_blocks(native_declared, args.native_arm,
                                                text)
    rescaled = rescaled + native_text

    # The shipped file's own inline comments on rescaled lines are replaced,
    # because several of them state the Earth unit and would now be wrong. The
    # upstream annotation is in the vendored source the header names.
    by_class = {}
    for change in changes:
        by_class.setdefault(change["class"], set()).add(change["parameter"])
    class_lines = "\n".join(
        f"!// {kind:12s} {', '.join(sorted(by_class[kind]))}"
        for kind in ("ANNUAL_SUM", "YEAR_COUNT", "ANNUAL_RATE", "DERIVED")
        if kind in by_class) or "!// nothing rescaled"
    if "DERIVED" in by_class:
        class_lines += (
            "\n!//"
            "\n!// DERIVED is not a rescale. nyear_spinup is not an Earth duration"
            "\n!// converted into simulation years; it is already in simulation"
            "\n!// years, read from lib/run_lengths.py:ecological_run_cycles,"
            "\n!// which reads the acceptance artifacts. Rescaling it would apply"
            "\n!// the YEAR_COUNT factor to a number already in the target unit."
            f"\n!// Basis: {spinup.get('basis', '')}"
            f"\n!// It is a FLOOR: {spinup.get('floor_because', '')}")

    if native_records:
        native_summary = "\n".join(
            [f"!// NATIVE to Vesper, appended after the conversion and NOT rescaled:"]
            + [f"!//   {r['name']:6s} inherits {r['inherits']:6s} include "
               f"{r['include']}  "
               + (", ".join(f"{k} at {v}" for k, v in r["bracket_ends"].items())
                  or "no bracket")
               + f"  ({r['derivation']})"
               for r in native_records]
            + ["!// Declared in biosphere/config/native_pfts.yaml. Their values are in",
               "!// model units already; rescaling one would apply the year-length",
               "!// factor to a number that is already in the target unit."])
    else:
        native_summary = "!// no Vesper-native types are declared"

    header = f"""!///////////////////////////////////////////////////////////////////////////////
!// GENERATED by biosphere/scripts/build_vesper_pfts.py. Do not edit.
!//
!// {args.source.name}, which is in EARTH units, converted to MODEL units, where
!// one year is one orbit. The registry that decides each parameter's class is
!// biosphere/notes/time-base-unit-contract.md.
!//
!// orbit        {orbital_days:.4f} Earth days at {config['orbit']['baseline_flux_earth']} S-Earth
!// model year   {model_year_days} steps of 24 h
!// Earth year   {EARTH_YEAR_DAYS:.4f} days
!// factor       {factor:.6f} Earth years per simulation year
!//
{class_lines}
!// NOT rescaled {", ".join(DELIBERATELY_UNSCALED)}
!//
!// ANNUAL_SUM scales down with the year, YEAR_COUNT up by the reciprocal, and
!// ANNUAL_RATE converts as 1-(1-r)^factor rather than multiplying, so that a
!// fraction of 1.0 stays 1.0. phengdd5ramp is deliberately untouched: it is a
!// within-season accumulation, already in absolute time, and scaling it would be
!// a real error.
!//
!// Inline comments on rescaled lines are regenerated, because the shipped ones
!// state the Earth unit. The upstream annotation is in the source named above.
!//
{native_summary}
!//
!// generated {datetime.now(timezone.utc).isoformat(timespec="seconds")}
!///////////////////////////////////////////////////////////////////////////////

"""
    output = args.output or (GENERATED / "vesper_pfts.ins")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(header + rescaled)

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": str(args.source),
        "source_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        # The parsed config, so that a later "is this artifact still current?"
        # can be answered per KEY. `config_sha256` alone cannot: it moves for an
        # edited comment exactly as it does for an edited parameter, and
        # `check_consistency.py` reported these artifacts stale on a comment
        # change until it had this to read. Same field, same purpose as
        # `source_config` in an ExoPlaSim run manifest.
        "source_config": config,
        "orbital_year_earth_days": orbital_days,
        "model_year_days": model_year_days,
        "earth_year_days": EARTH_YEAR_DAYS,
        "scale_factor": factor,
        "contract": "biosphere/notes/time-base-unit-contract.md",
        "annual_sum_scaled_down": ANNUAL_SUM,
        "year_count_scaled_up": YEAR_COUNT,
        "annual_rate_converted": ANNUAL_RATE,
        "deliberately_unscaled": DELIBERATELY_UNSCALED,
        "changes": changes,
        "spinup": spinup,
        "native_pfts": {
            "declaration": rel(NATIVE_PFTS),
            "declaration_sha256": hashlib.sha256(NATIVE_PFTS.read_bytes()).hexdigest(),
            "arm_requested": args.native_arm,
            "types": native_records,
        },
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    report["generator"] = "biosphere/scripts/build_vesper_pfts.py"
    report_path = output.with_name(output.stem + "_provenance.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")

    print(f"orbit   {orbital_days:.4f} Earth days, model year {model_year_days} "
          f"-> factor {factor:.6f}")
    print(f"changed {len(changes)} parameter values:")
    for change in changes:
        print(f"   {change['class']:12s} {change['parameter']:14s} "
              f"{change['from']:>8g} -> {change['written']:g}")
    for record in native_records:
        print(f"native  {record['kind']:5s} {record['name']:12s} "
              f"inherits {record['inherits']:12s} "
              + (f"include {record['include']}, " if record["kind"] == "pft" else "shared,      ")
              + (", ".join(f"{k} at the {v} end"
                           for k, v in record["bracket_ends"].items())
                 or "no bracket") + ": "
              + ", ".join(f"{k} {v:g}" for k, v in record["values"].items()))

    # THE CONTROL, and it reads the artifact rather than the intent: resolve the
    # emitted type through the reader every consumer uses and require the
    # declared values to come back. A block that inherits from the wrong place,
    # or whose override is shadowed by a later group reference, passes every
    # check above and fails here.
    for record in native_records:
        resolved = lpj_pfts.parameters(record["name"], path=output)
        for key, want in record["values"].items():
            got = resolved.get(key)
            if got is None or abs(got - want) > 1e-9:
                raise SystemExit(
                    f"{record['name']}.{key} was declared {want:g} but resolves "
                    f"to {got!r} in the file just written")
        if record["kind"] == "pft" and resolved.get("include") != record["include"]:
            raise SystemExit(
                f"{record['name']} was declared include {record['include']} but "
                f"resolves to {resolved.get('include')!r}")
    if native_records:
        print(f"        every native value resolves back out of the written file")

    print(f"\nwrote {rel(output)}")
    print(f"      {report_path.name}")


if __name__ == "__main__":
    main()
