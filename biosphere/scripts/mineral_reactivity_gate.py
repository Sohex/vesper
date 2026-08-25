"""The mineral-reactivity gate: which arm the simulated soil protects carbon and
sorbs phosphorus on, and whether that arm is allowed to be running.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a vegetation model's soil organic matter transfers, the
pedology soil map that feeds them, and the declaration that says which of two
equation sets is active.

The active arm is texture-only: linear functions of clay, or clay plus silt,
inherited from CENTURY through `modules/somdynam.cpp`. The alternative is
mineral-aware, and it refuses, because the Fe-Al oxide, allophane, aggregate and
cation proxies it needs are not produced by anything in this pipeline. What
pedology DOES produce -- an andic areal fraction and an andic phosphate-fixation
share -- now crosses the interface into `Soiltype` and is read by no equation,
which is the state this gate exists to keep honest.

It can fail:

  drift      a coefficient in biosphere/config/mineral_reactivity.yaml whose
             literal no longer appears in modules/somdynam.cpp
  range      a declared equation that leaves its range over the texture domain
  closure    the four microbial-pool fractions summing past 1, which makes the
             remainder transfer to slow SOM run backwards
  column     a soil map that does not carry a declared column, or carries one
             outside its declared range
  reader     a line of the model reading a carried field while the arm that
             would use it refuses
Six reduced fixtures run on every invocation, five built to be wrong in a named
way. A fixture that does not get its verdict is a defect in this checker.

    python biosphere/scripts/mineral_reactivity_gate.py            # status
    python biosphere/scripts/mineral_reactivity_gate.py --strict   # refuses
                                                                   # while the
                                                                   # mineral-aware
                                                                   # arm cannot run

Fail-closed in one direction only, on the same terms as `bvoc_gate.py`. A run on
the declared texture-only arm is a correct run of a declared model boundary, so
the default arm reports and exits 0.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import CONFIG, GENERATED, GUESS_SOURCE, PROJECT_ROOT
from paths import rel  # noqa: E402

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = COMPONENT_ROOT / "config" / "mineral_reactivity.yaml"
REPORT = GENERATED / "mineral_reactivity_gate_report.json"

UNDECLARED = "undeclared"

# The texture simplex is sampled on this grid in each of sand and clay, silt
# taken as the remainder. The declared equations are linear in texture, so the
# extremes are at the corners and a grid this coarse finds them exactly.
TEXTURE_STEPS = 101

EVAL_NAMESPACE = {
    "__builtins__": {},
    "min": min, "max": max, "abs": abs, "pow": pow,
    "exp": math.exp, "log": math.log, "sqrt": math.sqrt,
}


def _simplex(steps: int = TEXTURE_STEPS):
    """Every (sand, clay, silt) with the three summing to 1."""
    for i in range(steps + 1):
        sand = i / steps
        for j in range(steps + 1 - i):
            clay = j / steps
            yield sand, clay, 1.0 - sand - clay


def _extrema(expression: str, points) -> tuple[float, float, tuple]:
    code = compile(expression, "<texture>", "eval")
    lo, hi = math.inf, -math.inf
    worst = None
    for sand, clay, silt in points:
        value = eval(code, EVAL_NAMESPACE,  # noqa: S307 - in-repo declaration
                     {"sand": sand, "clay": clay, "silt": silt})
        if value < lo:
            lo, worst = value, (sand, clay, silt)
        hi = max(hi, value)
    return lo, hi, worst


def _literal_tokens(text: str) -> set[str]:
    tokens = set()
    for raw in re.findall(r"\d+\.?\d*(?:[eE][-+]?\d+)?", text):
        tokens.add(raw)
        try:
            tokens.add(f"{float(raw):g}")
        except ValueError:
            pass
    return tokens


def read_soilmap(path: Path) -> tuple[list[str], list[list[float]]]:
    lines = path.read_text().splitlines()
    header = lines[0].split()
    rows = [[float(v) for v in line.split()] for line in lines[1:] if line.strip()]
    return header, rows


def check(declaration: dict, source_text: str, model_texts: dict[str, str],
          soilmap: tuple[list[str], list[list[float]]] | None) -> list[dict]:
    findings: list[dict] = []
    arm = declaration["texture_only"]
    tokens = _literal_tokens(source_text)
    simplex = list(_simplex())

    for equation in arm["equations"]:
        name = equation["name"]
        for literal in equation.get("literals", []):
            text = f"{literal:g}" if isinstance(literal, float) else str(literal)
            if text not in tokens and str(literal) not in tokens:
                findings.append({"kind": "drift", "what": name,
                                 "detail": f"literal {literal} is not in {arm['source_file']}"})
        try:
            lo, hi, _ = _extrema(equation["expression"], simplex)
        except Exception as error:
            findings.append({"kind": "range", "what": name,
                             "detail": f"could not be evaluated: {error}"})
            continue
        allowed_lo, allowed_hi = equation["range"]
        if lo < allowed_lo - 1e-12 or hi > allowed_hi + 1e-12:
            findings.append({
                "kind": "range", "what": name,
                "detail": (f"reaches [{lo:.6g}, {hi:.6g}] over the texture simplex, "
                           f"outside its declared [{allowed_lo}, {allowed_hi}]")})

    # The closure has two arms. Over the whole texture simplex it is a property
    # of the equations, and it is negative at pure sand, which is registered with
    # its value so a coefficient change that moves it fails here. Over the soil
    # map's own textures it is a property of this world, and it has to stay
    # forwards.
    closure = arm["closure"]
    lo, _, worst = _extrema(closure["expression"], simplex)
    declared = closure["simplex_minimum"]
    if abs(lo - declared) > closure["simplex_tolerance"]:
        findings.append({
            "kind": "closure", "what": closure["name"],
            "detail": (f"over the texture simplex it falls to {lo:.6g} at sand/clay/silt "
                       f"{worst[0]:.3f}/{worst[1]:.3f}/{worst[2]:.3f}, and the declared "
                       f"minimum is {declared}")})

    carried = declaration["carried_state"]
    if soilmap is not None:
        header, rows = soilmap
        index = {name: i for i, name in enumerate(header)}
        if all(c in index for c in ("sand", "clay", "silt")):
            points = [(row[index["sand"]], row[index["clay"]], row[index["silt"]])
                      for row in rows]
            map_lo, _, map_worst = _extrema(closure["expression"], points)
            if map_lo < closure["soilmap_minimum"] - 1e-12:
                findings.append({
                    "kind": "closure", "what": closure["name"],
                    "detail": (f"on the soil map's own textures it falls to {map_lo:.6g} at "
                               f"sand/clay/silt {map_worst[0]:.3f}/{map_worst[1]:.3f}/"
                               f"{map_worst[2]:.3f}, so the transfer runs backwards on this "
                               f"world and not only in principle")})
        lo_allowed, hi_allowed = carried["range"]
        for column in carried["soil_map_columns"]:
            if column not in index:
                findings.append({"kind": "column", "what": column,
                                 "detail": "the soil map does not carry this column"})
                continue
            values = [row[index[column]] for row in rows]
            if min(values) < lo_allowed - 1e-12 or max(values) > hi_allowed + 1e-12:
                findings.append({
                    "kind": "column", "what": column,
                    "detail": (f"spans [{min(values):.6g}, {max(values):.6g}], outside "
                               f"its declared [{lo_allowed}, {hi_allowed}]")})

    # No equation may read a carried field while the arm that would use it
    # refuses. `readers_permitted` is empty, so any use outside the input module
    # that put the value there is a finding.
    permitted = set(carried.get("readers_permitted", []))
    for member in ("andic_frac", "p_fixation_frac"):
        for filename, text in model_texts.items():
            if member in text and filename not in permitted:
                findings.append({
                    "kind": "reader", "what": member,
                    "detail": (f"{filename} reads it while the mineral-aware arm "
                               f"refuses; it is carried, not used")})

    return findings


def _fixtures(declaration, source_text, model_texts, soilmap) -> list[dict]:
    import copy

    def mutate(fn):
        d = copy.deepcopy(declaration)
        fn(d)
        return d

    def add_literal(index, literal):
        return lambda d: d["texture_only"]["equations"][index]["literals"].append(literal)

    def narrow_range(index, span):
        return lambda d: d["texture_only"]["equations"][index].__setitem__("range", span)

    def break_closure(d):
        d["texture_only"]["closure"]["expression"] = "-0.5 + 0.0*sand + 0.0*clay + 0.0*silt"

    def rename_column(d):
        d["carried_state"]["soil_map_columns"] = ["not_a_column"]

    def narrow_column(d):
        d["carried_state"]["range"] = [0.0, 0.001]

    cases = [
        ("the declaration as it stands", declaration, None),
        ("a coefficient the model source does not contain",
         mutate(add_literal(0, 4321.0)), "drift"),
        ("an equation declared narrower than it is",
         mutate(narrow_range(2, [0.0, 0.1])), "range"),
        ("a microbial partition whose remainder runs backwards",
         mutate(break_closure), "closure"),
        ("a carried column the soil map does not have",
         mutate(rename_column), "column"),
        ("a carried column declared narrower than the artifact",
         mutate(narrow_column), "column"),
    ]

    results = []
    for label, candidate, expect in cases:
        findings = check(candidate, source_text, model_texts, soilmap)
        kinds = {f["kind"] for f in findings}
        ok = (not findings) if expect is None else (expect in kinds)
        results.append({"fixture": label, "expected": expect or "clean",
                        "found": sorted(kinds), "pass": ok})
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict", action="store_true",
                        help="refuse while the mineral-aware arm cannot run")
    parser.add_argument("--soilmap", type=Path, default=None,
                        help="soil map to check; defaults to the configured build's")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    declaration = yaml.safe_load(DECLARATION.read_text())
    source_path = PROJECT_ROOT / declaration["texture_only"]["source_file"]
    source_text = source_path.read_text()

    # Every model translation unit that could read a carried field. The input
    # module that writes it is excluded by name, not by accident.
    model_texts = {}
    for path in sorted((GUESS_SOURCE / "modules").glob("*.cpp")):
        if path.name in ("soilinput.cpp",):
            continue
        model_texts[f"modules/{path.name}"] = path.read_text(errors="replace")

    soilmap = None
    soilmap_path = args.soilmap
    if soilmap_path is None:
        import yaml as _yaml
        import builds as _b
        try:
            soilmap_path = _b.soilmap(_yaml.safe_load(CONFIG.read_text()))
        except Exception:
            soilmap_path = None
    if soilmap_path is not None and Path(soilmap_path).is_file():
        soilmap = read_soilmap(Path(soilmap_path))

    fixtures = _fixtures(declaration, source_text, model_texts, soilmap)
    findings = check(declaration, source_text, model_texts, soilmap)

    undeclared = [
        name for name, spec in declaration["mineral_aware"]["required_proxies"].items()
        if spec.get("declared", UNDECLARED) == UNDECLARED
    ]

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": rel(DECLARATION),
        "active_arm": declaration["active_arm"],
        "source": rel(source_path),
        "soilmap": rel(Path(soilmap_path)) if soilmap is not None else None,
        "findings": findings,
        "undeclared_proxies": undeclared,
        "fixtures": fixtures,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print("Mineral reactivity: which arm protects carbon and sorbs phosphorus.\n")
        print(f"  active arm: {declaration['active_arm']}")
        broken = [f for f in fixtures if not f["pass"]]
        print(f"  fixtures: {len(fixtures) - len(broken)} of {len(fixtures)} got their verdict")
        for case in broken:
            print(f"    BROKEN: {case['fixture']}: expected {case['expected']}, found {case['found']}")
        if soilmap is None:
            print("  no soil map resolved, so the carried columns were not checked")
        else:
            print(f"  soil map: {rel(Path(soilmap_path))}, {len(soilmap[1])} land cells")
        if findings:
            print(f"\n  {len(findings)} finding(s):")
            for finding in findings:
                print(f"    [{finding['kind']}] {finding['what']}: {finding['detail']}")
        else:
            print("  every declared equation matches the source and holds its range")
        if undeclared:
            print(f"\n  the mineral-aware arm refuses: {len(undeclared)} proxy(ies) undeclared")
            for name in undeclared:
                spec = declaration["mineral_aware"]["required_proxies"][name]
                print(f"    {name}  (produced by: {spec['produced_by']})")
        print(f"\n  report: {rel(REPORT)}")

    if [f for f in fixtures if not f["pass"]]:
        print("\nA fixture did not get the verdict it was built for. That is a defect in",
              file=sys.stderr)
        print("this checker, not in the declaration.", file=sys.stderr)
        return 2
    if findings:
        return 1
    if args.strict and undeclared:
        print("\n--strict: refused. The mineral-aware arm needs every proxy above, and",
              file=sys.stderr)
        print("nothing in this pipeline produces them.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
