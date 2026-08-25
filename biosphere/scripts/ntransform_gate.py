"""The soil nitrogen transformation gate: what the operator is declared to be,
and whether the model still agrees.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a vegetation model's soil nitrogen chemistry, the file that
declares it, and the source the model actually reads.

`ifntransform 1` is on in the Vesper baseline, so nitrification, denitrification
and volatilisation set how much mineral nitrogen the simulated plants can reach.
The operator is Xu-Ri and Prentice (2008) as ported into LPJ-GUESS: every one of
its response functions is a function of soil temperature, upper-soil water,
water-filled pore space or pH and of nothing else. It sees no pressure, no oxygen partial pressure, no gas
diffusivity and no water table, which is why this world's declared atmosphere
cannot reach it. The argument is in
`biosphere/notes/soil-nitrogen-transformation-parameterisation.md`.

This module is the enforcement, and it can fail:

  drift        a constant in biosphere/config/ntransform.yaml that no longer
               matches the instruction file the model reads it from, or a
               literal in a declared response function that no longer appears
               in modules/ntransform.cpp
  range        a response function that leaves the range its role allows,
               sampled over the whole declared domain. A pool multiplier
               outside [0,1] creates or destroys nitrogen
  product      a chain of factors multiplying one pool mass whose product can
               exceed 1, which takes more nitrogen out of a pool than is in it,
               or a chain declared to be held inside its pool by an explicit
               min() in the operator that the operator no longer contains
  bound        a declared constant outside the bound plib parses it against
  calibration  an entry in the calibration block whose declared verdict is not
               what the arithmetic says: an `agrees` whose value is not the
               paper's value or is outside the paper's range, an `outside` that
               has moved inside it, or an entry naming a constant or a response
               function that does not exist

A dozen reduced fixtures run on every invocation, all but one built to be wrong
in a named way. A fixture that does not get the verdict it was built for is a
defect in this checker rather than in the declaration.

    python biosphere/scripts/ntransform_gate.py            # status, exit 0
    python biosphere/scripts/ntransform_gate.py --strict   # refuses on the
                                                           # named residual

`--strict` names exactly what is still undeclared and refuses on that and
nothing else: the Vesper preconditions that carry the sentinel, plus every
calibration entry whose verdict is `outside` or `unsourced`. An entry the
sources settle is no longer part of the refusal.

It is fail-closed in one direction only, on the same terms as `bvoc_gate.py`. A
run on the Earth-calibrated operator, declared as such, is a correct run of a
declared model boundary, so the default arm reports and exits 0. `--strict` is
the arm that refuses.

Nothing here is verified by execution. LPJ-GUESS does not build on this tree,
so every statement this module makes is against the source and the declaration.
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

from _paths import GENERATED, GUESS_SOURCE, PROJECT_ROOT
from paths import rel  # noqa: E402

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = COMPONENT_ROOT / "config" / "ntransform.yaml"
REPORT = GENERATED / "ntransform_gate_report.json"
OPERATOR_SOURCE = GUESS_SOURCE / "modules" / "ntransform.cpp"

UNDECLARED = "undeclared"

# How finely each declared domain is sampled. A response function here is
# smooth, so a maximum found on a grid this dense is the maximum to well inside
# the slack the role bounds allow.
SAMPLES = 2001

# The role each function's range is judged by. A pool multiplier and a partition
# fraction are both held to [0,1] because both multiply a mass; a factor is held
# to its own declared range and to the products it appears in.
ROLES = ("pool_multiplier", "partition_fraction", "factor")


def richards(a: float, b: float, c: float, d: float, x: float) -> float:
    """guessmath.h's richards_curve, which four of the declared forms call."""
    return a + (b - a) / (1.0 + math.exp(-c * (x - d)))


EVAL_NAMESPACE = {
    "__builtins__": {},
    "min": min, "max": max, "abs": abs, "pow": pow,
    "exp": math.exp, "log": math.log, "sqrt": math.sqrt,
    "richards": richards,
}


def _linspace(lo: float, hi: float, n: int) -> list[float]:
    if n < 2:
        return [lo]
    step = (hi - lo) / (n - 1)
    return [lo + step * i for i in range(n)]


def evaluate(spec: dict, domains: dict, samples: int = SAMPLES
             ) -> tuple[float, float]:
    """The minimum and maximum of one declared form over its whole domain.

    One variable is swept at `samples` points and any others at a coarser grid,
    so a two-variable form is a full product grid rather than a diagonal.
    """
    variables = spec["variables"]
    names = list(variables)
    axes = []
    for index, name in enumerate(names):
        n = samples if index == 0 and len(names) == 1 else 41
        source = variables[name]
        if source == "unit_interval":
            axes.append(_linspace(0.0, 1.0, n))
        else:
            lo, hi = domains[source]
            axes.append(_linspace(float(lo), float(hi), n))

    code = compile(spec["form"], f"<{spec['name']}>", "eval")
    lo = math.inf
    hi = -math.inf

    def walk(index: int, env: dict) -> None:
        nonlocal lo, hi
        if index == len(names):
            value = eval(code, EVAL_NAMESPACE, env)  # noqa: S307 - in-repo form
            lo = min(lo, value)
            hi = max(hi, value)
            return
        for point in axes[index]:
            env[names[index]] = point
            walk(index + 1, env)

    walk(0, {})
    return lo, hi


def _literal_tokens(text: str) -> set[str]:
    """Every numeric token in the operator source, normalised."""
    tokens = set()
    for raw in re.findall(r"\d+\.?\d*(?:[eE][-+]?\d+)?", text):
        tokens.add(raw)
        try:
            tokens.add(f"{float(raw):g}")
        except ValueError:
            pass
    return tokens


# The verdicts the calibration block may carry. `agrees` and `outside` are
# claims about a number and are re-derived here; `unsourced` is a claim that
# no source states the quantity, which no arithmetic can check.
VERDICTS = ("agrees", "outside", "unsourced")

# What `--strict` refuses on, beyond the undeclared preconditions.
REFUSING_VERDICTS = ("outside", "unsourced")

# How an entry is held to its source. `value` is a constant checked against a
# stated value or range; `form` is an equation compared with the paper's, with
# no number to check; `derived` is a bracket on a quantity computed from the
# model rather than on any one of its constants.
COMPARISONS = ("value", "form", "derived")


def _check_calibration(declaration: dict, declared_values: dict,
                       functions: dict) -> list[dict]:
    """Every calibration entry, against the source and the number it names.

    An entry addresses one of three things. `instruction:<key>` is a constant
    the model reads from the instruction file, and its number is the one the
    drift check has already tied to that file, so a `paper_value` or a `bracket`
    on it is checked directly. `function:<name>` is a declared response
    function; a bracket on one is a bracket on a single constant inside it,
    which the entry names as `model_value` and which has to be among that
    function's declared literals, so neither the model nor this file can move
    without the other. `form:<name>` is an expression in the operator that is
    not a declared function, argued in prose and carrying no number.
    """
    findings: list[dict] = []
    calibration = declaration.get("calibration")
    if not calibration:
        return findings
    sources = calibration.get("sources", {})

    def bad(what: str, detail: str) -> None:
        findings.append({"kind": "calibration", "what": what, "detail": detail})

    for entry in calibration.get("entries", []):
        what = entry.get("what", "?")
        verdict = entry.get("verdict")
        if verdict not in VERDICTS:
            bad(what, f"unknown verdict {verdict!r}")
            continue
        if not entry.get("source"):
            bad(what, "carries no source line")
        elif verdict != "unsourced" and not any(k in entry["source"] for k in sources):
            bad(what, "its source line names no declared source key")
        if not entry.get("why"):
            bad(what, "carries no argument")

        kind, _, name = what.partition(":")
        value = None
        if kind == "instruction":
            if name not in declared_values:
                bad(what, f"{name} is not a declared constant")
                continue
            value = declared_values[name]
        elif kind == "function":
            if name not in functions:
                bad(what, f"{name} is not a declared response function")
                continue
            if "model_value" in entry:
                value = float(entry["model_value"])
                literals = [float(x) for x in functions[name].get("literals", [])]
                if not any(abs(value - x) <= 1e-12 for x in literals):
                    bad(what, (f"model_value {value:g} is not among {name}'s declared "
                               "literals, so it is not a constant of that function"))
                    continue
        elif kind == "form":
            if not entry.get("reads"):
                bad(what, "a form entry has to say what it reads")
            continue
        else:
            bad(what, f"unknown target kind {kind!r}")
            continue

        paper_value = entry.get("paper_value")
        bracket = entry.get("bracket")
        compared = entry.get("compared", "value")
        if compared not in COMPARISONS:
            bad(what, f"unknown comparison {compared!r}")
            continue
        if verdict == "unsourced":
            if paper_value is not None:
                bad(what, "declared unsourced but carries a paper value")
            continue
        if compared == "form":
            # The paper states an equation and the model's is the same one, or
            # is not. There is no number in the declaration to re-derive that
            # from, so the argument in `why` is the whole of it.
            if paper_value is not None or bracket is not None:
                bad(what, "compared as a form but carries a number")
            continue
        if paper_value is None and bracket is None:
            bad(what, f"verdict {verdict!r} with neither a paper value nor a bracket")
            continue
        if compared == "derived":
            # The bracket is on a quantity computed from the model rather than
            # on any one constant in it, so `why` carries the derivation.
            if value is not None and "model_value" in entry:
                bad(what, "compared as derived but names a model constant")
            continue
        if value is None:
            bad(what, "compared as a value but names no number to compare")
            continue
        inside = True
        if paper_value is not None:
            inside = abs(value - float(paper_value)) <= 1e-12
        if bracket is not None:
            lo, hi = float(bracket[0]), float(bracket[1])
            inside = inside and (lo - 1e-12 <= value <= hi + 1e-12)
        if verdict == "agrees" and not inside:
            bad(what, (f"declared to agree, but the model's {value:g} is neither the "
                       "paper's value nor inside its range"))
        if verdict == "outside" and inside:
            bad(what, (f"declared to be outside its source's range, but the model's "
                       f"{value:g} is inside it"))
    return findings


def check(declaration: dict, source_text: str, instruction_text: str
          ) -> list[dict]:
    """Every check, as a list of findings. An empty list is a clean gate."""
    findings: list[dict] = []
    domains = declaration["domains"]
    instruction = declaration["instruction_parameters"]
    declared_values = instruction["values"]
    bounds = instruction.get("parser_bounds", {})

    parsed = {}
    for line in instruction_text.splitlines():
        line = line.split("!", 1)[0].strip()
        parts = line.split()
        if len(parts) >= 2 and parts[0] in declared_values:
            try:
                parsed[parts[0]] = float(parts[1])
            except ValueError:
                pass

    for key, value in declared_values.items():
        if key not in parsed:
            findings.append({"kind": "drift", "what": key,
                             "detail": "declared here but not set in the instruction file"})
        elif parsed[key] != value:
            findings.append({"kind": "drift", "what": key,
                             "detail": f"declared {value}, instruction file says {parsed[key]}"})
        if key in bounds:
            lo, hi = bounds[key]
            if not (lo <= value <= hi):
                findings.append({"kind": "bound", "what": key,
                                 "detail": f"{value} is outside the parser bound [{lo}, {hi}]"})

    tokens = _literal_tokens(source_text)
    extrema: dict[str, tuple[float, float]] = {}
    functions = {spec["name"]: spec for spec in declaration["response_functions"]}

    for spec in declaration["response_functions"]:
        name = spec["name"]
        if spec["role"] not in ROLES:
            findings.append({"kind": "role", "what": name,
                             "detail": f"unknown role {spec['role']!r}"})
            continue
        for literal in spec.get("literals", []):
            text = f"{literal:g}" if isinstance(literal, float) else str(literal)
            if text not in tokens and str(literal) not in tokens:
                findings.append({"kind": "drift", "what": name,
                                 "detail": f"literal {literal} is not in modules/ntransform.cpp"})
        try:
            lo, hi = evaluate(spec, domains)
        except Exception as error:  # a form that cannot be evaluated is a finding
            findings.append({"kind": "range", "what": name,
                             "detail": f"could not be evaluated: {error}"})
            continue
        extrema[name] = (lo, hi)
        allowed_lo, allowed_hi = spec["range"]
        if spec["role"] in ("pool_multiplier", "partition_fraction"):
            allowed_lo, allowed_hi = max(allowed_lo, 0.0), min(allowed_hi, 1.0)
        if lo < allowed_lo - 1e-12 or hi > allowed_hi + 1e-12:
            findings.append({
                "kind": "range", "what": name,
                "detail": (f"reaches [{lo:.6g}, {hi:.6g}] over its domain, "
                           f"outside the [{allowed_lo}, {allowed_hi}] its role allows")})

    findings.extend(_check_calibration(declaration, declared_values, functions))

    for product in declaration.get("products", []):
        worst = 1.0
        missing = False
        for factor in product["factors"]:
            if factor.startswith("instruction:"):
                key = factor.split(":", 1)[1]
                if key not in declared_values:
                    findings.append({"kind": "product", "what": product["name"],
                                     "detail": f"factor {key} is not a declared constant"})
                    missing = True
                    break
                worst *= declared_values[key]
            elif factor in extrema:
                worst *= extrema[factor][1]
            else:
                findings.append({"kind": "product", "what": product["name"],
                                 "detail": f"factor {factor} is not a declared response function"})
                missing = True
                break
        if missing:
            continue
        clamp = product.get("clamped_by")
        if clamp is not None and clamp not in source_text:
            findings.append({
                "kind": "product", "what": product["name"],
                "detail": (f"declares that {clamp!r} holds it inside the pool, and "
                           "modules/ntransform.cpp does not contain that clamp")})
        elif worst > product["max"] + 1e-12 and clamp is None:
            findings.append({
                "kind": "product", "what": product["name"],
                "detail": (f"the product of its factors reaches {worst:.6g}, "
                           f"above the {product['max']} that keeps it inside the pool")})

    return findings


# The fixtures. Six are built to be wrong in a named way; the seventh is the
# declaration as it stands, which has to come back clean. A fixture that does
# not get its verdict is a defect in the checker.
def _fixtures(declaration: dict, source_text: str, instruction_text: str
              ) -> list[dict]:
    import copy

    def mutate(fn):
        d = copy.deepcopy(declaration)
        fn(d)
        return d

    def set_value(key, value):
        return lambda d: d["instruction_parameters"]["values"].__setitem__(key, value)

    def widen_domain(name, span):
        return lambda d: d["domains"].__setitem__(name, span)

    def break_literal(index, literal):
        return lambda d: d["response_functions"][index]["literals"].append(literal)

    def unbounded_form(index, form):
        def apply(d):
            d["response_functions"][index]["form"] = form
        return apply

    def bad_factor(name):
        def apply(d):
            d["products"][0]["factors"] = ["instruction:" + name]
        return apply

    def drop_clamp(name):
        def apply(d):
            for product in d["products"]:
                if product["name"] == name:
                    product["clamped_by"] = "min(no_such_pool,"
                    return
            raise KeyError(name)
        return apply

    def _entry(d, what):
        for entry in d["calibration"]["entries"]:
            if entry["what"] == what:
                return entry
        raise KeyError(what)

    def calibration_claim(what, field, value):
        return lambda d: _entry(d, what).__setitem__(field, value)

    def calibration_target(what, target):
        return lambda d: _entry(d, what).__setitem__("what", target)

    cases = [
        ("the declaration as it stands", declaration, None),
        ("a constant that drifted from the instruction file",
         mutate(set_value("k_N", 0.084)), "drift"),
        ("a constant outside the bound plib parses it against",
         mutate(set_value("f_nitri_max", 2.0)), "bound"),
        ("a literal the operator source does not contain",
         mutate(break_literal(0, 1234.5)), "drift"),
        ("a pool multiplier that leaves [0,1] on its domain",
         mutate(unbounded_form(1, "2.0 + 0.0*pH + 0.0*w + 0.0*fT")), "range"),
        ("a pH domain wide enough to make volatilisation unbounded",
         mutate(widen_domain("soil_ph", [3.8, 20.0])), "range"),
        ("a product chain that takes more than the pool holds",
         mutate(bad_factor("nonexistent_constant")), "product"),
        ("a chain declared clamped by a min() the operator does not contain",
         mutate(drop_clamp("denitrification_no3_to_no2")), "product"),
        ("a constant declared to agree with a bracket that excludes it",
         mutate(calibration_claim("instruction:f_nitri_gas_max", "bracket", [0.5, 0.9])),
         "calibration"),
        ("a disagreement declared where the paper's range in fact contains the value",
         mutate(calibration_claim("instruction:k_N", "verdict", "outside")),
         "calibration"),
        ("a bracket held to a number that is not a constant of the function it names",
         mutate(calibration_claim("function:f_denitri_water", "model_value", 9.0)),
         "calibration"),
        ("a calibration entry naming a response function that does not exist",
         mutate(calibration_target("function:f_denitri_water", "function:no_such_function")),
         "calibration"),
    ]

    results = []
    for label, candidate, expect in cases:
        findings = check(candidate, source_text, instruction_text)
        kinds = {f["kind"] for f in findings}
        if expect is None:
            ok = not findings
        else:
            ok = expect in kinds
        results.append({"fixture": label, "expected": expect or "clean",
                        "found": sorted(kinds), "pass": ok})
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict", action="store_true",
                        help="refuse while any Vesper precondition is undeclared")
    parser.add_argument("--json", action="store_true", help="report only")
    args = parser.parse_args()

    declaration = yaml.safe_load(DECLARATION.read_text())
    source_text = OPERATOR_SOURCE.read_text()
    instruction_path = PROJECT_ROOT / declaration["instruction_parameters"]["file"]
    instruction_text = instruction_path.read_text()

    fixtures = _fixtures(declaration, source_text, instruction_text)
    findings = check(declaration, source_text, instruction_text)

    undeclared = [
        name for name, spec in declaration["vesper_preconditions"].items()
        if spec.get("declared", UNDECLARED) == UNDECLARED
    ]
    unsettled = [
        {"what": entry["what"], "verdict": entry["verdict"],
         "owner": entry.get("owner", "?")}
        for entry in declaration.get("calibration", {}).get("entries", [])
        if entry.get("verdict") in REFUSING_VERDICTS
    ]
    settled = [
        entry["what"]
        for entry in declaration.get("calibration", {}).get("entries", [])
        if entry.get("verdict") == "agrees"
    ]

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": rel(DECLARATION),
        "operator_source": rel(OPERATOR_SOURCE),
        "instruction_file": rel(instruction_path),
        "findings": findings,
        "undeclared_preconditions": undeclared,
        "calibration_unsettled": unsettled,
        "calibration_settled": settled,
        "fixtures": fixtures,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print("The soil nitrogen transformation operator, checked against the source.\n")
        broken = [f for f in fixtures if not f["pass"]]
        print(f"  fixtures: {len(fixtures) - len(broken)} of {len(fixtures)} got their verdict")
        for case in broken:
            print(f"    BROKEN: {case['fixture']}: expected {case['expected']}, found {case['found']}")
        if findings:
            print(f"\n  {len(findings)} finding(s):")
            for finding in findings:
                print(f"    [{finding['kind']}] {finding['what']}: {finding['detail']}")
        else:
            print("  the declaration matches the source and every bound holds")
        total = len(declaration.get("calibration", {}).get("entries", []))
        print(f"\n  calibration: {len(settled)} of {total} entries agree with the source they name")
        if unsettled:
            print("  the rest are what remains undeclared, and are what --strict refuses on:")
            for item in unsettled:
                print(f"    [{item['verdict']}] {item['what']}  [{item['owner']}]")
        if undeclared:
            print(f"\n  {len(undeclared)} Vesper precondition(s) undeclared, so this is an")
            print("  Earth-calibrated operator run on this world's water and pH:")
            for name in undeclared:
                owner = declaration["vesper_preconditions"][name].get("owner", "?")
                print(f"    {name}  [{owner}]")
        print(f"\n  report: {rel(REPORT)}")

    if [f for f in fixtures if not f["pass"]]:
        print("\nA fixture did not get the verdict it was built for. That is a defect in",
              file=sys.stderr)
        print("this checker, not in the declaration.", file=sys.stderr)
        return 2
    if findings:
        return 1
    if args.strict and (undeclared or unsettled):
        print("\n--strict: refused, on exactly what is still undeclared and nothing else.",
              file=sys.stderr)
        if undeclared:
            print(f"  {len(undeclared)} Vesper precondition(s): "
                  + ", ".join(undeclared), file=sys.stderr)
        if unsettled:
            print(f"  {len(unsettled)} calibration entr(y/ies) outside or unsourced: "
                  + ", ".join(item["what"] for item in unsettled), file=sys.stderr)
        print("Each has to be declared before this operator's output is a Vesper result.",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
