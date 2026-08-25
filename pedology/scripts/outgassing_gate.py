"""Whether this world's prescribed pCO2 is backed by a plausible outgassing rate.

`weathering_fluxes.py` already computes the outgassing the prescribed CO2
requires and compares it against a mass-scaling estimate of what the planet could
supply. That comparison had every input a gate needs -- a requirement, a supply,
a margin -- and none of the machinery. There was no threshold, no branch, no
refusal and no exit code, and the verdict was a constant string: had the margin
come back at 0.3, the same sentence would have been printed.

This turns it into a gate on the shape `biosphere/scripts/bvoc_gate.py`
established. `pedology/config/outgassing.yaml` is the declaration; with
`requested: false` this reports what is undeclared and exits 0, and with
`requested: true` -- which is what claiming anywhere that the prescribed CO2 is
SUPPORTED by outgassing means -- it refuses by name for every unmet precondition.

    python pedology/scripts/outgassing_gate.py
    python pedology/scripts/outgassing_gate.py --strict
    python pedology/scripts/outgassing_gate.py --json

The subject throughout is the SIMULATED planet: its modelled silicate weathering,
its prescribed atmosphere, and the outgassing rate that atmosphere would need.

VOLC-8.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import ANALYSIS, CONFIG, OUTGASSING, PROJECT_ROOT

from paths import rel  # noqa: E402  from lib/, via _paths.

REPORT = ANALYSIS / "outgassing_gate_report.json"
UNDECLARED = "undeclared"

# Section, field, what it is, and the issue that owns the decision. Supply comes
# first because a supply side that cannot state its own scaling makes every later
# refinement unusable: the margin is a ratio and one half of it is asserted.
PRECONDITIONS = (
    ("supply", "mass_scaling_exponent",
     "the power of planet mass the outgassing supply estimate claims", "volc-8"),
    ("supply", "mass_scaling_source",
     "where that exponent comes from", "volc-8"),
    ("supply", "per_km_segment_bracket",
     "the per-kilometre boundary CO2 flux bracket, spanning the observed range", "volc-8"),
    ("supply", "per_km_segment_source",
     "where that bracket comes from", "volc-8"),
    ("supply", "source_partition",
     "how the global requirement splits across ridge, arc, plume/rift and metamorphic", "volc-8"),
    ("supply", "degassing_efficiency",
     "the mantle/slab carbon and degassing efficiencies behind that partition", "volc-8"),
    ("supply", "atmosphere_versus_ocean_split",
     "which fraction reaches the atmosphere and which the ocean", "ocn-13"),
)

# The four source classes the partition has to name. A partition missing one is
# not a partition of the same thing.
PARTITION_CLASSES = ("ridge", "arc_recycled", "plume_rift", "metamorphic")

# Every declared bracket, and the value it brackets where there is one. Written
# out rather than discovered by key name: a check that guesses which keys are
# brackets stops checking the moment someone names one differently, and stops
# silently.
BRACKETS = (
    ("supply.mass_scaling_exponent",
     ("supply", "mass_scaling_exponent"),
     ("supply", "mass_scaling_exponent_bracket")),
    ("supply.per_km_segment_bracket",
     None,
     ("supply", "per_km_segment_bracket")),
    ("supply.atmosphere_versus_ocean_split.ocean_fraction",
     ("supply", "atmosphere_versus_ocean_split", "ocean_fraction"),
     ("supply", "atmosphere_versus_ocean_split", "ocean_fraction_bracket")),
    ("supply.degassing_efficiency.arc_recycled.subducted_carbon_returned_to_surface",
     None,
     ("supply", "degassing_efficiency", "arc_recycled",
      "subducted_carbon_returned_to_surface")),
) + tuple(
    (f"supply.source_partition.{cls}.flux_mol_per_year",
     ("supply", "source_partition", cls, "flux_mol_per_year"),
     ("supply", "source_partition", cls, "bracket_mol_per_year"))
    for cls in PARTITION_CLASSES
)

# How far a fraction may sit from the flux it claims to be a fraction of, and how
# far the fractions may sit from summing to 1. Both are rounding room for values
# quoted to three decimal places, fixed here before any declaration was checked.
PARTITION_TOLERANCE = 1e-3


class Refusal:
    """One named reason the claim is refused.

    `code` is what the message says and what a later report is grepped for. It is
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


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_declaration(path: Path = OUTGASSING) -> dict:
    if not path.is_file():
        raise SystemExit(
            f"OUTGASSING-DECLARATION-MISSING: {rel(path)} does not exist. "
            f"There is no default declaration.")
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def read_requirement(declaration: dict) -> tuple[dict | None, Path]:
    """The computed requirement, from the artifact and never from prose."""
    rel_path = ((declaration.get("requirement") or {}).get("artifact")
                or "pedology/analysis/weathering_fluxes.json")
    path = PROJECT_ROOT / rel_path
    if not path.is_file():
        return None, path
    try:
        return json.loads(path.read_text(encoding="utf-8")), path
    except (json.JSONDecodeError, UnicodeDecodeError):
        return {}, path


def _get(declaration: dict, section: str, field: str):
    return (declaration.get(section) or {}).get(field, UNDECLARED)


def _declared(value) -> bool:
    return value is not None and value != UNDECLARED


def margin(requirement: dict | None, declaration: dict, planet: dict) -> dict:
    """Supply over requirement, both as multiples of Earth's.

    The supply half is DERIVED from the declared exponent rather than assumed to
    be the mass ratio. Where the exponent is undeclared this returns no supply at
    all, which is the honest state and is what makes the margin refusable rather
    than merely small.
    """
    out: dict = {"required_over_earth": None, "supply_over_earth": None, "margin": None}
    if requirement:
        check = (((requirement.get("carbon_balance") or {}).get("plausibility_check")) or {})
        value = check.get("required_outgassing_over_earth")
        if isinstance(value, (int, float)) and value > 0:
            out["required_over_earth"] = float(value)
    exponent = _get(declaration, "supply", "mass_scaling_exponent")
    mass = planet.get("planet", {}).get("mass_earth")
    if _declared(exponent) and isinstance(mass, (int, float)):
        out["supply_over_earth"] = float(mass) ** float(exponent)
    if out["required_over_earth"] and out["supply_over_earth"]:
        out["margin"] = out["supply_over_earth"] / out["required_over_earth"]
    return out


def _dig(declaration: dict, path: tuple[str, ...]):
    """Follow a key path, returning the sentinel rather than raising."""
    node = declaration
    for key in path:
        if not isinstance(node, dict) or key not in node:
            return UNDECLARED
        node = node[key]
    return node


def check_supply_arithmetic(declaration: dict) -> list[Refusal]:
    """Is the declared supply side self-consistent?

    Declaring a bracket and declaring a partition are the two places this file
    can be wrong WITHOUT being undeclared, which is the failure the undeclared
    checks above cannot see. A bracket whose ends are the wrong way round, a
    value outside the bracket that is supposed to contain it, and fractions that
    do not sum to 1 or do not match the fluxes they are quoted from are all
    arithmetic, so each has a right answer and can be refused rather than
    reviewed.

    Nothing here fires on an undeclared field: that is already refused by name,
    and refusing it twice would say the same thing in two voices.
    """
    refusals: list[Refusal] = []

    for label, value_path, bracket_path in BRACKETS:
        bracket = _dig(declaration, bracket_path)
        if not _declared(bracket):
            continue
        if (not isinstance(bracket, (list, tuple)) or len(bracket) != 2
                or not all(isinstance(x, (int, float)) for x in bracket)):
            refusals.append(Refusal(
                "OUTGASSING-BRACKET-NOT-A-PAIR",
                f"{label} is bracketed by {bracket!r}, which is not two numbers",
                "volc-8"))
            continue
        low, high = float(bracket[0]), float(bracket[1])
        if not low <= high:
            refusals.append(Refusal(
                "OUTGASSING-BRACKET-DISORDERED",
                f"{label} is bracketed [{low}, {high}], whose low end is above "
                f"its high end", "volc-8"))
            continue
        if value_path is None:
            continue
        value = _dig(declaration, value_path)
        if not isinstance(value, (int, float)):
            continue
        if not low <= float(value) <= high:
            refusals.append(Refusal(
                "OUTGASSING-VALUE-OUTSIDE-ITS-BRACKET",
                f"{label} is {value} and its own bracket is [{low}, {high}]. A "
                f"value its bracket does not contain is one of the two, not both",
                "volc-8"))

    partition = _dig(declaration, ("supply", "source_partition"))
    if not _declared(partition):
        return refusals
    if not isinstance(partition, dict):
        refusals.append(Refusal(
            "OUTGASSING-PARTITION-NOT-A-PARTITION",
            f"supply.source_partition is {type(partition).__name__}, not a "
            f"mapping of source class to fraction", "volc-8"))
        return refusals

    missing = [c for c in PARTITION_CLASSES
               if not isinstance(partition.get(c), dict)
               or not isinstance(partition[c].get("fraction"), (int, float))]
    if missing:
        refusals.append(Refusal(
            "OUTGASSING-PARTITION-CLASS-MISSING",
            f"supply.source_partition names no fraction for "
            f"{', '.join(missing)}, so it partitions something other than the "
            f"four declared source classes", "volc-8"))
        return refusals

    fractions = {c: float(partition[c]["fraction"]) for c in PARTITION_CLASSES}
    total = sum(fractions.values())
    if abs(total - 1.0) > PARTITION_TOLERANCE:
        refusals.append(Refusal(
            "OUTGASSING-PARTITION-DOES-NOT-SUM",
            f"the four source fractions sum to {total:.4f}, not 1", "volc-8"))

    fluxes = {c: partition[c].get("flux_mol_per_year") for c in PARTITION_CLASSES}
    if all(isinstance(v, (int, float)) and v > 0 for v in fluxes.values()):
        flux_total = sum(float(v) for v in fluxes.values())
        for c in PARTITION_CLASSES:
            implied = float(fluxes[c]) / flux_total
            if abs(implied - fractions[c]) > PARTITION_TOLERANCE:
                refusals.append(Refusal(
                    "OUTGASSING-PARTITION-DISAGREES-WITH-ITS-FLUXES",
                    f"source class {c} is declared at fraction {fractions[c]} "
                    f"while its own flux is {implied:.4f} of the declared total. "
                    f"A partition assembled from fluxes has to be the fluxes",
                    "volc-8"))

    return refusals


def evaluate(declaration: dict, planet: dict, requirement: dict | None,
             requirement_path: Path) -> list[Refusal]:
    refusals: list[Refusal] = []

    for section, field, what, issue in PRECONDITIONS:
        if not _declared(_get(declaration, section, field)):
            refusals.append(Refusal(
                f"OUTGASSING-UNDECLARED-{section.upper()}-{field.upper().replace('_', '-')}",
                f"{what} is undeclared", issue))

    refusals.extend(check_supply_arithmetic(declaration))

    acceptance = declaration.get("acceptance") or {}
    floor = acceptance.get("minimum_margin")
    if not isinstance(floor, (int, float)):
        refusals.append(Refusal(
            "OUTGASSING-MARGIN-FLOOR-MISSING",
            "acceptance.minimum_margin is not a number, so nothing decides "
            "whether a margin is enough", "volc-8"))
    elif floor <= 1.0:
        refusals.append(Refusal(
            "OUTGASSING-MARGIN-FLOOR-NOT-A-MARGIN",
            f"acceptance.minimum_margin is {floor}, which is not a margin: at or "
            f"below 1 the supply does not exceed the requirement", "volc-8"))

    if requirement is None:
        refusals.append(Refusal(
            "OUTGASSING-REQUIREMENT-ARTIFACT-MISSING",
            f"{rel(requirement_path)} does not exist, so the required outgassing "
            f"has no artifact behind it. config/planet.yaml already records that "
            f"the check backing the prescribed CO2 is unbacked; this is that "
            f"statement made refusable", "volc-8"))
    else:
        computed = margin(requirement, declaration, planet)
        if computed["required_over_earth"] is None:
            refusals.append(Refusal(
                "OUTGASSING-REQUIREMENT-NOT-IN-ARTIFACT",
                f"{rel(requirement_path)} carries no positive "
                f"carbon_balance.plausibility_check.required_outgassing_over_earth",
                "volc-8"))
        if (acceptance.get("requirement_must_match_configured_build")
                and requirement.get("source_build")
                and requirement["source_build"] != planet.get("source_build")):
            refusals.append(Refusal(
                "OUTGASSING-REQUIREMENT-WRONG-BUILD",
                f"the requirement was computed on build "
                f"{requirement['source_build']!r} and the configured build is "
                f"{planet.get('source_build')!r}. A requirement from another "
                f"terrain backs nothing, whatever its margin", "volc-8"))
        if (computed["margin"] is not None and isinstance(floor, (int, float))
                and computed["margin"] < floor):
            refusals.append(Refusal(
                "OUTGASSING-MARGIN-BELOW-FLOOR",
                f"margin {computed['margin']:.3f} is below the declared floor "
                f"{floor}", "volc-8"))

    premises = declaration.get("premises_the_terrain_cannot_supply") or {}
    if premises.get("claims_an_absolute_rate"):
        refusals.append(Refusal(
            "OUTGASSING-ABSOLUTE-RATE-CLAIMED",
            "an absolute outgassing rate is claimed, and the terrain cannot "
            "support one: Orogen's plate velocities are in arbitrary units and "
            "the model has no time axis, so relative throughput is all that is "
            "available", "volc-8"))

    return refusals


def verdict(declaration: dict, planet: dict, requirement: dict | None,
            requirement_path: Path) -> dict:
    """The conclusion, DERIVED. This is what replaces a constant string.

    Callers that want the sentence rather than the machinery use this; it is what
    `weathering_fluxes.py` writes into its report instead of a verdict it did not
    reach.
    """
    refusals = evaluate(declaration, planet, requirement, requirement_path)
    computed = margin(requirement, declaration, planet)
    floor = (declaration.get("acceptance") or {}).get("minimum_margin")
    if refusals:
        reason = "; ".join(r.code for r in refusals)
        sentence = (
            "NOT SHOWN. The prescribed CO2 is not backed by a demonstrated "
            f"outgassing rate, for {len(refusals)} named reasons: {reason}. It "
            "remains an assumption.")
    else:
        sentence = (
            f"The prescribed CO2 is attainable with margin {computed['margin']:.2f}, "
            f"against a floor of {floor} fixed before the number was computed.")
    return {"granted": not refusals, "verdict": sentence,
            "refusals": [r.as_dict() for r in refusals], **computed}


def granted(declaration: dict | None = None, planet: dict | None = None
            ) -> tuple[bool, list[Refusal]]:
    declaration = declaration if declaration is not None else read_declaration()
    planet = planet if planet is not None else yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    requirement, path = read_requirement(declaration)
    refusals = evaluate(declaration, planet, requirement, path)
    if not declaration.get("requested"):
        # Granted nothing and refuses nothing: no claim is being made.
        return False, refusals
    return not refusals, refusals


def require(declaration: dict | None = None, planet: dict | None = None) -> bool:
    declaration = declaration if declaration is not None else read_declaration()
    ok, refusals = granted(declaration, planet)
    if declaration.get("requested") and refusals:
        raise SystemExit(
            "the prescribed pCO2 is claimed to be outgassing-backed, and it is "
            "refused:\n  " + "\n  ".join(str(r) for r in refusals))
    return ok


def _satisfied(declaration: dict) -> dict:
    """Every precondition met. Placeholders, not proposals: nothing here is a
    value anyone should copy into the declaration."""
    d = copy.deepcopy(declaration)
    for section, field, _what, _issue in PRECONDITIONS:
        d.setdefault(section, {})[field] = "declared, for the fixture"
    d.setdefault("supply", {})["mass_scaling_exponent"] = 1.0
    # The arithmetic checks need shapes rather than sentinels, so the fixture
    # supplies the plainest self-consistent ones there are: an equal four-way
    # split of four equal fluxes, and brackets that are one point wide. These are
    # placeholders in the same sense as the strings above and are not estimates.
    d["supply"]["source_partition"] = {
        cls: {"fraction": 0.25, "flux_mol_per_year": 1.0e+12,
              "bracket_mol_per_year": [1.0e+12, 1.0e+12]}
        for cls in PARTITION_CLASSES
    }
    d["supply"]["mass_scaling_exponent_bracket"] = [1.0, 1.0]
    d["supply"]["per_km_segment_bracket"] = [1.0, 1.0]
    d["supply"]["atmosphere_versus_ocean_split"] = {
        "ocean_fraction": 0.25, "ocean_fraction_bracket": [0.25, 0.25]}
    d["supply"]["degassing_efficiency"] = {
        "arc_recycled": {"subducted_carbon_returned_to_surface": [0.0, 1.0]}}
    d.setdefault("acceptance", {})["minimum_margin"] = 1.5
    d.setdefault("premises_the_terrain_cannot_supply", {})["claims_an_absolute_rate"] = False
    return d


def _repaired_requirement(planet: dict) -> dict:
    """An artifact that would back the claim, for the must-be-granted fixture."""
    return {"source_build": planet.get("source_build"),
            "carbon_balance": {"plausibility_check": {
                "required_outgassing_over_earth": 1.0}}}


def _fixtures(declaration: dict, planet: dict) -> list[dict]:
    """Each is a declaration that MUST be refused by a named code, and each
    asserts that code appears. The first is the declaration as it stands, which
    has to be refused for the reason the contract says it is refused for; the
    rest are mutations built to be wrong in one way each. A fixture that does not
    get its verdict is a defect in this checker, not in the declaration."""
    def mutate(base, section, field, value):
        candidate = copy.deepcopy(base)
        candidate.setdefault(section, {})[field] = value
        return candidate

    good = _satisfied(declaration)
    ok_req = _repaired_requirement(planet)
    # Any path that is not the artifact, so "absent" is what the gate sees.
    missing = Path(__file__).resolve()

    cases = [
        # The live declaration. Its specific code moves as fields are declared
        # and repairs land, and a fixture that fails because someone made
        # progress is a fixture nobody will keep, so this one only has to be
        # refused for SOME reason.
        ("the declaration as it stands", declaration, None, None),
        ("a supply side that will not say what power of mass it claims",
         mutate(good, "supply", "mass_scaling_exponent", UNDECLARED), ok_req,
         "OUTGASSING-UNDECLARED-SUPPLY-MASS-SCALING-EXPONENT"),
        ("a supply side with no per-kilometre segment bracket",
         mutate(good, "supply", "per_km_segment_bracket", UNDECLARED), ok_req,
         "OUTGASSING-UNDECLARED-SUPPLY-PER-KM-SEGMENT-BRACKET"),
        ("an acceptance floor that is not a margin",
         mutate(good, "acceptance", "minimum_margin", 1.0), ok_req,
         "OUTGASSING-MARGIN-FLOOR-NOT-A-MARGIN"),
        ("an acceptance floor that is not a number",
         mutate(good, "acceptance", "minimum_margin", "generous"), ok_req,
         "OUTGASSING-MARGIN-FLOOR-MISSING"),
        ("a claim with no requirement artifact behind it",
         good, None, "OUTGASSING-REQUIREMENT-ARTIFACT-MISSING"),
        ("a requirement computed on a different terrain",
         good, {**ok_req, "source_build": "some-other-build"},
         "OUTGASSING-REQUIREMENT-WRONG-BUILD"),
        ("a requirement whose margin is below the declared floor",
         good, {"source_build": planet.get("source_build"),
                "carbon_balance": {"plausibility_check": {
                    "required_outgassing_over_earth": 100.0}}},
         "OUTGASSING-MARGIN-BELOW-FLOOR"),
        ("an artifact that carries no requirement at all",
         good, {"source_build": planet.get("source_build")},
         "OUTGASSING-REQUIREMENT-NOT-IN-ARTIFACT"),
        # The arithmetic the declaration can get wrong while still being fully
        # declared. Each of these passes every undeclared check.
        ("a bracket whose low end is above its high end",
         mutate(good, "supply", "per_km_segment_bracket", [4.74e+08, 1.52e+06]), ok_req,
         "OUTGASSING-BRACKET-DISORDERED"),
        ("a bracket that is not two numbers",
         mutate(good, "supply", "per_km_segment_bracket", "a hundredfold"), ok_req,
         "OUTGASSING-BRACKET-NOT-A-PAIR"),
        ("an exponent its own bracket does not contain",
         mutate(good, "supply", "mass_scaling_exponent_bracket", [1.1, 1.34]), ok_req,
         "OUTGASSING-VALUE-OUTSIDE-ITS-BRACKET"),
        ("a partition that does not sum to 1",
         mutate(good, "supply", "source_partition",
                {**good["supply"]["source_partition"],
                 "ridge": {**good["supply"]["source_partition"]["ridge"],
                           "fraction": 0.40}}), ok_req,
         "OUTGASSING-PARTITION-DOES-NOT-SUM"),
        ("a partition whose fractions are not the fluxes it quotes them from",
         mutate(good, "supply", "source_partition",
                {**good["supply"]["source_partition"],
                 "ridge": {**good["supply"]["source_partition"]["ridge"],
                           "flux_mol_per_year": 4.0e+12}}), ok_req,
         "OUTGASSING-PARTITION-DISAGREES-WITH-ITS-FLUXES"),
        ("a partition that names three source classes and not four",
         mutate(good, "supply", "source_partition",
                {k: v for k, v in good["supply"]["source_partition"].items()
                 if k != "metamorphic"}), ok_req,
         "OUTGASSING-PARTITION-CLASS-MISSING"),
        ("a supply side claiming an absolute outgassing rate the terrain cannot give",
         mutate(good, "premises_the_terrain_cannot_supply", "claims_an_absolute_rate", True),
         ok_req, "OUTGASSING-ABSOLUTE-RATE-CLAIMED"),
        # The other direction, and the one that says this is a gate rather than a
        # wall: a declaration with every precondition met, against an artifact
        # that backs it, must come back with NO refusals. A gate nothing can
        # satisfy refuses for a reason that is never written down.
        ("every precondition met, against a backing artifact", good, ok_req, "GRANTED"),
    ]

    results = []
    for label, candidate, req, expect in cases:
        codes = {r.code for r in evaluate(candidate, planet, req, missing)}
        if expect == "GRANTED":
            ok = not codes
        elif expect is None:
            ok = bool(codes)
        else:
            ok = expect in codes
        results.append({"fixture": label,
                        "expected": expect or "refused for some reason",
                        "refusals": sorted(codes), "pass": ok})
    return results


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--declaration", type=Path, default=OUTGASSING)
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero on any refusal, even when not requested")
    ap.add_argument("--json", action="store_true", help="print the report instead of prose")
    args = ap.parse_args()

    declaration = read_declaration(args.declaration)
    planet = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    requirement, requirement_path = read_requirement(declaration)
    result = verdict(declaration, planet, requirement, requirement_path)
    fixtures = _fixtures(declaration, planet)
    broken = [f for f in fixtures if not f["pass"]]

    report = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "generator": "pedology/scripts/outgassing_gate.py",
        "declaration": str(rel(args.declaration)),
        "declaration_sha256": sha256(args.declaration),
        "config_sha256": sha256(CONFIG),
        "requested": bool(declaration.get("requested")),
        "requirement_artifact": str(rel(requirement_path)),
        "requirement_artifact_present": requirement is not None,
        **result,
        "fixtures": fixtures,
        "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT,
                                     capture_output=True, text=True).stdout.strip() or None,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        state = "REQUESTED" if report["requested"] else "not requested"
        print(f"outgassing claim: {state}")
        print(f"  requirement artifact  {report['requirement_artifact']} "
              f"{'present' if report['requirement_artifact_present'] else 'ABSENT'}")
        print(f"  verdict  {result['verdict']}")
        if result["refusals"]:
            print(f"\n{len(result['refusals'])} refusals:")
            for r in result["refusals"]:
                tail = f"  [{r['issue']}]" if r["issue"] else ""
                print(f"  {r['code']}: {r['detail']}{tail}")
        print(f"\nfixtures: {len(fixtures) - len(broken)}/{len(fixtures)} pass")
        print(f"wrote {rel(REPORT)}")

    if broken:
        print("\nA fixture did not get the verdict it was built for. That is a "
              "defect in this checker, not in the declaration.")
        for f in broken:
            print(f"  {f['fixture']}: expected {f['expected']}, got {f['refusals']}")
        raise SystemExit(1)
    if result["refusals"] and (args.strict or report["requested"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
