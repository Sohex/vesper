"""The soil phosphorus gate: what this project's phosphorus path is declared to
be, and whether the model still agrees.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a vegetation model's soil organic phosphorus chemistry, the
file that declares it, and the source the model actually reads.

`biosphere/config/somdynam.yaml` is the declaration and this module is the
enforcement. The subject is the C:P ramps of `modules/somdynam.cpp`: the two
saturation thresholds, the five `setptoc` calls that read them, the lines the
phosphorus argument rests on, and every place this project's source departs from
the vendored LPJ-GUESS-CNP fork.

THE REFERENCE POINT IS A SUBTREE COMMIT AND NOT A RELEASE. `ntransform_gate.py`
checks divergences from LPJ-GUESS 4.1.1, which can be named by Zenodo record and
SVN revision. LPJ-GUESS 4.1.1 has no phosphorus at all, so nothing here can be
compared against a release: the whole path arrived with the CNP fork, and the
only fixed point a divergence can be measured against is the commit that subtree
was imported at. The declaration carries that commit and this module refuses a
register that names none.

It can fail:

  drift        a constant in the declaration whose line is no longer in
               modules/somdynam.cpp, or whose value is not what that line's own
               initialiser evaluates to; a ramp whose call the source no longer
               contains, or contains a different number of times
  range        a ramp that leaves the P:C range its two endpoints allow, swept
               over its whole driver domain and past both ends of it. A
               receiving pool's P:C outside [1/ctop_max, 1/ctop_min] is a
               stoichiometry neither endpoint of the figure's line contains
  invariant    a line the phosphorus argument rests on that the model no longer
               runs. The threshold conversion is one definition propagating only
               while the phosphorus-limitation-off pin reads the same symbol, so
               the pin going is the argument being reverted without the constant
               moving
  divergence   a declared divergence from the vendored CNP fork whose form the
               model no longer records beside the changed one, which it has gone
               back to running, or whose changed line it no longer contains. All
               three are checked, because a divergence that is not recorded is a
               silent fork, one that is recorded but live is a declaration that
               has drifted from the model, and one whose changed line is gone has
               been reverted with the record left behind

Reduced fixtures run on every invocation, all but one built to be wrong in a
named way. A fixture that does not get the verdict it was built for is a defect
in this checker rather than in the declaration.

    python biosphere/scripts/somdynam_gate.py            # status, exit 0
    python biosphere/scripts/somdynam_gate.py --strict   # refuses on the
                                                         # named residual

`--strict` refuses on exactly two things: a constant declared `sourced: false`,
and a divergence whose verdict is `gate`. The first is the refusal
`framework/parameters.cpp` already makes on `ifplim 1`, restated where the
constant is declared rather than left in a C++ error string. There is no gated
divergence; the verdict exists so the next one has somewhere to go.

The default arm reports and exits 0. This project runs `ifplim 0`, where two of
the three divergences are inert in the flows and the third moves one pool's C:P,
so a run on the declared phosphorus path is a correct run of a declared model
boundary.

NOTHING IN THE DIVERGENCE REGISTER IS VERIFIED BY A CONTROLLED EXECUTION, and
that is stated in both arms rather than left to be inferred. LPJ-GUESS now
builds and the ifplim-0 Vesper arm has run, but no matched arm has restored the
fork forms, and two entries are live only under the configuration
`framework/parameters.cpp` refuses. The declaration names that comparison arm.
"""

from __future__ import annotations

import argparse
import ast
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import GENERATED, GUESS_SOURCE, PROJECT_ROOT
from paths import rel  # noqa: E402

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = COMPONENT_ROOT / "config" / "somdynam.yaml"
REPORT = GENERATED / "somdynam_gate_report.json"

# How finely each ramp's driver is swept. setptoc is piecewise linear, so a
# maximum found on a grid this dense is the maximum exactly; the density is
# there so the sweep lands on both breakpoints and inside both branches.
SAMPLES = 2001

# How far past each end of [fmin, fmax] the sweep runs. The clamped branches are
# where a mis-declared pair shows up, and the driver reaches them: pmin_mass is
# pinned at fmax under ifplim 0 and the litter concentration sits three orders
# below fmin.
OVERRUN = 0.5

# What a declared divergence from the fork may say about itself. `keep` is one
# this project reviewed and stood behind, declared at its own site in the
# source; `gate` is one that stands but is not settled, and `--strict` refuses
# while it does. A divergence that does not survive review is reverted and
# leaves the declaration, so there is no `revert`.
DIVERGENCE_VERDICTS = ("keep", "gate")

# What shape the fork's own version has. See the declaration's own comment: the
# two are guarded differently, and `commented_out` is the weaker guard.
MAINLINE_FORMS = ("line", "commented_out")

# Which configuration a divergence is live in. `ifplim_1` is inert in the arm
# this project runs and waiting in the one it does not; `neither` is an entry
# whose value the model overwrites before reading.
BITES = ("both", "ifplim_0", "ifplim_1", "neither")


def _strip_comments(text: str) -> str:
    """The source with C and C++ comments removed.

    The divergence check needs to tell a form RECORDED in a comment from the
    same form LIVE in the code, which is the whole point of recording it: the
    fork's line has to stand beside the changed one without being what the model
    runs.
    """
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", " ", text)


def setptoc(fac: float, ctop_max: float, ctop_min: float,
            fmin: float, fmax: float) -> float:
    """modules/somdynam.cpp's setptoc, as the P:C it writes.

    Reproduced rather than approximated: this is the function whose range the
    declaration bounds, and a paraphrase of it would bound something else.
    """
    if fac <= fmin:
        return 1.0 / ctop_max
    if fac >= fmax:
        return 1.0 / ctop_min
    return 1.0 / (ctop_min + (ctop_max - ctop_min) * (fmax - fac) / (fmax - fmin))


def _linspace(lo: float, hi: float, n: int) -> list[float]:
    if n < 2:
        return [lo]
    step = (hi - lo) / (n - 1)
    return [lo + step * i for i in range(n)]


def _number(value: float) -> str:
    """A C++ double literal as somdynam.cpp writes it: one decimal place."""
    return f"{value:.1f}"


def ramp_call(entry: dict) -> str:
    """The call string the source has to contain, built from the declaration.

    Built and not declared, so a pair changed in the source and a pair changed
    in the declaration fail the same check. Nothing here can move alone.
    """
    return (f"setptoc(soil, {entry['driver']}, {entry['pool']}, "
            f"{_number(entry['ctop_max'])}, {_number(entry['ctop_min'])}, "
            f"{_number(entry['fmin'])}, {entry['saturates_at']});")


def _initialiser(line: str) -> float | None:
    """The value a declaration line's own initialiser evaluates to.

    `static const double PMASS_SAT = 0.002 * 6.6;` has to be 0.0132 in the
    declaration, and the arithmetic that says so is the source's own.
    """
    match = re.search(r"=\s*(.+?);\s*$", line.strip())
    if not match:
        return None
    expression = match.group(1)
    if not re.fullmatch(r"[0-9eE+\-*/.() ]+", expression):
        return None
    try:
        return float(eval(compile(ast.Expression(ast.parse(expression, mode="eval").body),
                                  "<initialiser>", "eval"), {"__builtins__": {}}, {}))
    except Exception:
        return None


def _check_divergences(declaration: dict, sources: dict) -> list[dict]:
    """Every declared divergence from the vendored fork, against its source.

    Three checks on each entry, and a divergence needs all three. Each
    `mainline` string has to be IN the source it names, so the record stands
    beside the changed code and cannot go missing; NOT in that source once
    comments are stripped, so a silent revert to the fork fails; and `live` has
    to be in the stripped source, so deleting the changed line fails as well.
    Without the third, reverting by deleting both would pass.

    A `commented_out` entry is one whose fork version is an absent line. Its
    recorded form carries its own `//` and can never survive comment stripping,
    so the second check is vacuous there and is declared vacuous rather than
    left to look like a check: the gate requires the recorded form to in fact
    begin with a comment marker, so the weaker guard cannot be claimed for a
    line that would support the stronger one, and the revert it guards against
    -- the live line commented out again -- is caught by the third check.
    """
    findings: list[dict] = []
    register = declaration.get("mainline_divergences")
    if not register:
        findings.append({"kind": "divergence", "what": "mainline_divergences",
                         "detail": "the declaration carries no divergence register"})
        return findings

    def bad(what: str, detail: str) -> None:
        findings.append({"kind": "divergence", "what": what, "detail": detail})

    for field in ("reference", "not_a_release", "why_not_verified", "comparison_arm"):
        if not register.get(field):
            bad("mainline_divergences", f"the register carries no {field}")
    reference = register.get("reference") or ""
    if not re.search(r"\b[0-9a-f]{40}\b", reference):
        bad("mainline_divergences",
            "the register's reference names no subtree commit. The phosphorus "
            "path is not in any LPJ-GUESS release, so a divergence here can be "
            "measured against nothing else")
    if register.get("execution_verified") is not False:
        bad("mainline_divergences",
            "execution_verified is not false, but no matched fork-form arm has "
            "run, so no divergence is verified in both directions")

    seen = set()
    for entry in register.get("entries", []):
        what = entry.get("id", "?")
        if what in seen:
            bad(what, "two register entries share one id")
        seen.add(what)
        if entry.get("verdict") not in DIVERGENCE_VERDICTS:
            bad(what, f"unknown divergence verdict {entry.get('verdict')!r}")
        if not entry.get("owner"):
            bad(what, "a divergence naming no owner")
        if entry.get("bites_under") not in BITES:
            bad(what, (f"unknown bites_under {entry.get('bites_under')!r}. An entry "
                       "inert under the configuration this project runs is waiting, "
                       "not harmless, and has to say which one it is live in"))
        for field in ("settles", "worth"):
            if not entry.get(field):
                bad(what, f"a divergence saying nothing about what it {field}")

        if entry.get("where") != "operator":
            bad(what, f"unknown divergence site {entry.get('where')!r}")
            continue
        source_file = entry.get("source_file")
        if source_file not in sources:
            bad(what, f"names {source_file!r}, which this gate does not read")
            continue
        source_text, code_only = sources[source_file]

        form_kind = entry.get("mainline_form")
        if form_kind not in MAINLINE_FORMS:
            bad(what, f"unknown mainline_form {entry.get('mainline_form')!r}")
            continue

        mainline = entry.get("mainline") or []
        if not mainline:
            bad(what, "a divergence naming no form the fork runs")
        for form in mainline:
            if form not in source_text:
                bad(what, (f"declares a divergence from {form!r}, and "
                           f"{source_file} does not record that form beside the "
                           "changed one"))
                continue
            if form_kind == "commented_out":
                if not form.lstrip().startswith("//"):
                    bad(what, (f"declares {form!r} as a line the fork leaves "
                               "commented out, and it carries no comment marker. "
                               "The weaker guard cannot be claimed for a line "
                               "that would support the stronger one"))
            elif form in code_only:
                bad(what, (f"declares a divergence from {form!r}, and "
                           f"{source_file} still runs it"))
        live = entry.get("live")
        if not live:
            bad(what, "a divergence naming no line the model runs instead")
        elif live not in code_only:
            bad(what, (f"declares that {source_file} runs {live!r} instead, "
                       "and it does not"))
    return findings


def check(declaration: dict, sources: dict) -> list[dict]:
    """Every check, as a list of findings. An empty list is a clean gate."""
    findings: list[dict] = []
    findings.extend(_check_divergences(declaration, sources))

    def source_for(block: dict, what: str) -> tuple[str, str] | None:
        name = block.get("source_file")
        if name not in sources:
            findings.append({"kind": "drift", "what": what,
                             "detail": f"names {name!r}, which this gate does not read"})
            return None
        return sources[name]

    constants: dict[str, float] = {}
    block = declaration.get("constants", {})
    pair = source_for(block, "constants")
    if pair:
        _, code_only = pair
        for item in block.get("values", []):
            name = item["name"]
            line = item["declaration"]
            if line not in code_only:
                findings.append({"kind": "drift", "what": name,
                                 "detail": (f"declared as {line!r}, and "
                                            f"{block['source_file']} does not run "
                                            "that line")})
            derived = _initialiser(line)
            if derived is None:
                findings.append({"kind": "drift", "what": name,
                                 "detail": f"the declared line {line!r} has no readable initialiser"})
            elif abs(derived - float(item["value"])) > 1e-12 * max(1.0, abs(derived)):
                findings.append({"kind": "drift", "what": name,
                                 "detail": (f"declared {item['value']}, and its own "
                                            f"line evaluates to {derived:.10g}")})
            if item.get("sourced") and not item.get("source"):
                findings.append({"kind": "drift", "what": name,
                                 "detail": "declared sourced and naming no source"})
            constants[name] = float(item["value"])

    block = declaration.get("ramps", {})
    pair = source_for(block, "ramps")
    if pair:
        _, code_only = pair
        for entry in block.get("entries", []):
            name = entry["name"]
            call = ramp_call(entry)
            found = code_only.count(call)
            wanted = int(entry["occurrences"])
            if found != wanted:
                findings.append({"kind": "drift", "what": name,
                                 "detail": (f"declared as {call!r} appearing {wanted} "
                                            f"time(s), and {block['source_file']} runs "
                                            f"it {found}")})
            symbol = entry["saturates_at"]
            if symbol not in constants:
                findings.append({"kind": "drift", "what": name,
                                 "detail": f"saturates at {symbol}, which is not a declared constant"})
                continue
            fmax = constants[symbol]
            fmin = float(entry["fmin"])
            span = fmax - fmin
            lo_edge = fmin - OVERRUN * abs(span)
            hi_edge = fmax + OVERRUN * abs(span)
            values = [setptoc(fac, float(entry["ctop_max"]), float(entry["ctop_min"]),
                              fmin, fmax)
                      for fac in _linspace(lo_edge, hi_edge, SAMPLES)]
            lo, hi = min(values), max(values)
            allowed_lo, allowed_hi = (float(v) for v in entry["ptoc_range"])
            if not (0.0 < allowed_lo <= allowed_hi <= 1.0):
                findings.append({"kind": "range", "what": name,
                                 "detail": (f"declares a P:C range [{allowed_lo}, "
                                            f"{allowed_hi}], which is not a positive "
                                            "ratio below unity")})
            if lo < allowed_lo - 1e-12 or hi > allowed_hi + 1e-12:
                findings.append({"kind": "range", "what": name,
                                 "detail": (f"reaches [{lo:.6g}, {hi:.6g}] over its "
                                            f"driver domain, outside the "
                                            f"[{allowed_lo:.6g}, {allowed_hi:.6g}] its "
                                            "two endpoints allow")})

    block = declaration.get("invariants", {})
    pair = source_for(block, "invariants")
    if pair:
        _, code_only = pair
        for entry in block.get("entries", []):
            if not entry.get("why"):
                findings.append({"kind": "invariant", "what": entry["name"],
                                 "detail": "an invariant saying nothing about what rests on it"})
            if entry["line"] not in code_only:
                findings.append({"kind": "invariant", "what": entry["name"],
                                 "detail": (f"rests on {entry['line']!r}, and "
                                            f"{block['source_file']} does not run it")})
    return findings


# The fixtures. The first is the declaration as it stands and has to come back
# clean; every other is built to be wrong in a named way. A fixture that does
# not get its verdict is a defect in the checker.
def _fixtures(declaration: dict, sources: dict) -> list[dict]:
    import copy

    def mutate(fn):
        d = copy.deepcopy(declaration)
        fn(d)
        return d

    def _divergence(d, entry_id):
        for entry in d["mainline_divergences"]["entries"]:
            if entry["id"] == entry_id:
                return entry
        raise KeyError(entry_id)

    def divergence_claim(entry_id, field, value):
        return lambda d: _divergence(d, entry_id).__setitem__(field, value)

    def _constant(d, name):
        for item in d["constants"]["values"]:
            if item["name"] == name:
                return item
        raise KeyError(name)

    def constant_claim(name, field, value):
        return lambda d: _constant(d, name).__setitem__(field, value)

    def _ramp(d, name):
        for entry in d["ramps"]["entries"]:
            if entry["name"] == name:
                return entry
        raise KeyError(name)

    def ramp_claim(name, field, value):
        return lambda d: _ramp(d, name).__setitem__(field, value)

    def invariant_claim(name, field, value):
        def apply(d):
            for entry in d["invariants"]["entries"]:
                if entry["name"] == name:
                    entry[field] = value
                    return
            raise KeyError(name)
        return apply

    cases = [
        ("the declaration as it stands", declaration, None),
        ("a divergence whose form the source does not record",
         mutate(divergence_claim("pmass_sat_labile_currency", "mainline",
                                 ["static const double NO_SUCH_LINE = 1.0;"])),
         "divergence"),
        ("a divergence from a form the source in fact still runs",
         mutate(divergence_claim("pmass_sat_labile_currency", "mainline",
                                 ["static const double PCONC_SAT = 0.02;"])),
         "divergence"),
        ("a divergence whose changed line the source no longer contains",
         mutate(divergence_claim("pmass_sat_labile_currency", "live",
                                 "static const double PMASS_SAT = no_such_call();")),
         "divergence"),
        ("a ramp reverted by commenting its call out again",
         mutate(divergence_claim("surfhumus_ptoc_ramp", "live",
                                 "setptoc(soil, pmin_mass, SURFHUMUS, 1.0, 1.0, 0.0, PMASS_SAT);")),
         "divergence"),
        ("the weaker guard claimed for a line that would support the stronger one",
         mutate(divergence_claim("surfhumus_ptoc_init", "mainline_form", "commented_out")),
         "divergence"),
        ("a register that claims a divergence has been run",
         mutate(lambda d: d["mainline_divergences"].__setitem__(
             "execution_verified", True)),
         "divergence"),
        ("a register whose reference names no subtree commit",
         mutate(lambda d: d["mainline_divergences"].__setitem__(
             "reference", "the CNP fork")),
         "divergence"),
        ("a divergence not saying which configuration it is live in",
         mutate(divergence_claim("pmass_sat_labile_currency", "bites_under", "sometimes")),
         "divergence"),
        ("a constant whose declaration line the source does not run",
         mutate(constant_claim("PMASS_SAT", "declaration",
                               "static const double PMASS_SAT = 0.002 * 11.3;")),
         "drift"),
        ("a constant whose value is not what its own line evaluates to",
         mutate(constant_claim("PMASS_SAT", "value", 0.002)),
         "drift"),
        ("a ramp whose call the source does not contain",
         mutate(ramp_claim("slowsom_ptoc", "ctop_min", 91.0)),
         "drift"),
        ("a ramp called from fewer places than declared",
         mutate(ramp_claim("surfmicro_ptoc", "occurrences", 3)),
         "drift"),
        ("a ramp whose declared P:C range its own endpoints exceed",
         mutate(ramp_claim("slowsom_ptoc", "ptoc_range", [0.005, 0.006])),
         "range"),
        ("a ramp declared over a P:C range that is not a ratio below unity",
         mutate(ramp_claim("slowsom_ptoc", "ptoc_range", [0.005, 4.0])),
         "range"),
        ("an invariant the source no longer runs",
         mutate(invariant_claim("labile_p_pin_reads_the_threshold", "line",
                                "soil.pmass_labile = 0.0132;")),
         "invariant"),
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict", action="store_true",
                        help="refuse while a phosphorus constant has no source")
    parser.add_argument("--json", action="store_true", help="report only")
    args = parser.parse_args()

    declaration = yaml.safe_load(DECLARATION.read_text())

    wanted = {declaration[block]["source_file"]
              for block in ("constants", "ramps", "invariants")}
    wanted |= {entry["source_file"]
               for entry in declaration["mainline_divergences"]["entries"]
               if entry.get("source_file")}
    sources: dict[str, tuple[str, str]] = {}
    for name in sorted(wanted):
        text = (PROJECT_ROOT / name).read_text()
        sources[name] = (text, _strip_comments(text))

    fixtures = _fixtures(declaration, sources)
    findings = check(declaration, sources)

    register = declaration["mainline_divergences"]
    divergences = [
        {"id": entry.get("id", "?"), "verdict": entry.get("verdict", "?"),
         "owner": entry.get("owner", "?"),
         "source_file": entry.get("source_file", "?"),
         "bites_under": entry.get("bites_under", "?")}
        for entry in register.get("entries", [])
    ]
    gated = [d for d in divergences if d["verdict"] == "gate"]
    unsourced = [
        {"name": item["name"], "owner": item.get("owner", "?")}
        for item in declaration["constants"]["values"] if not item.get("sourced")
    ]

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": rel(DECLARATION),
        "sources": sorted(sources),
        "findings": findings,
        "fork_reference": register.get("reference"),
        "divergences": divergences,
        "unsourced_constants": unsourced,
        "execution_verified": bool(register.get("execution_verified")),
        "fixtures": fixtures,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print("The simulated soil's phosphorus path, checked against the source.\n")
        broken = [f for f in fixtures if not f["pass"]]
        print(f"  fixtures: {len(fixtures) - len(broken)} of {len(fixtures)} got their verdict")
        for case in broken:
            print(f"    BROKEN: {case['fixture']}: expected {case['expected']}, found {case['found']}")
        if findings:
            print(f"\n  {len(findings)} finding(s):")
            for finding in findings:
                print(f"    [{finding['kind']}] {finding['what']}: {finding['detail']}")
        else:
            print("  the declaration matches the source and every ramp holds its range")
        kept = [d for d in divergences if d["verdict"] == "keep"]
        print(f"\n  {len(divergences)} declared divergence(s) from the vendored CNP fork,")
        print(f"  {len(kept)} kept and {len(gated)} gated. The reference point is a")
        print("  SUBTREE COMMIT and not a release: LPJ-GUESS 4.1.1 has no phosphorus,")
        print(f"  so there is nothing else to measure these against.\n    {register.get('reference')}")
        print("  NONE IS VERIFIED BY A MATCHED FORK-FORM ARM. LPJ-GUESS now builds,")
        print("  but no divergence has been run in both directions, and two of these are")
        print("  live only under the ifplim 1 that parameters.cpp refuses.")
        for item in divergences:
            print(f"    [{item['verdict']}] {item['id']}  bites under "
                  f"{item['bites_under']}  ({item['source_file']})  [{item['owner']}]")
        if unsourced:
            print(f"\n  {len(unsourced)} saturation constant(s) with no phosphorus source,")
            print("  which is what --strict refuses on and what parameters.cpp already")
            print("  refuses ifplim 1 for:")
            for item in unsourced:
                print(f"    {item['name']}  [{item['owner']}]")
        print(f"\n  report: {rel(REPORT)}")

    if broken := [f for f in fixtures if not f["pass"]]:
        print(f"\n{len(broken)} fixture(s) did not get the verdict they were built for.",
              file=sys.stderr)
        print("That is a defect in this checker, not in the declaration.", file=sys.stderr)
        return 2
    if findings:
        return 1
    if args.strict and (unsourced or gated):
        print("\n--strict: refused, on exactly what is still unsourced and nothing else.",
              file=sys.stderr)
        if unsourced:
            print("  saturation constant(s) with no phosphorus source: "
                  + ", ".join(item["name"] for item in unsourced), file=sys.stderr)
        if gated:
            print("  divergence(s) gated: "
                  + ", ".join(item["id"] for item in gated), file=sys.stderr)
        print("Each has to be settled before ifplim 1 is a Vesper result.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
