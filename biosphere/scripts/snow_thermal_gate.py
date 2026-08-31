"""The snow thermal gate: what the vegetation model's snowpack conducts and
stores, and whether its conductivity is still the one the climate model runs.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a vegetation model's snowpack, the two thermal properties it
derives from that pack's density, and the source the model actually reads.

The two are a CONDUCTIVITY, which is a material relation with a source and a
published bracket that both columns have to share, and a volumetric HEAT
CAPACITY, which is a density times a specific heat and is settled by dimensional
argument alone. Only the conductivity is a cross-column declaration, so only it
is held to `lib/snow.py`; the capacity is a divergence from the release like any
other and is checked as one.

`biosphere/config/snow_thermal.yaml` is the declaration and this module is the
enforcement. There are two halves and they fail for different reasons.

THE RELATION IS ONE DECLARATION. `lib/snow.py` states Fourteau et al. (2021)
Eq. (18) once. The climate model's `landmod.f90` and this model's `soil.cpp`
each carry the adopted row as a literal, because a Fortran model and a C++ model
cannot import a Python module at runtime, and `lib/snow.py`'s
`check_restatements()` holds both literals to the one table. This gate runs it,
and so does `scripts/smoke_test.py`. Without that, the register could certify a
divergence into a relation the climate column had since moved off, which is the
two-agreeing-copies state the declaration exists to remove.

THE DIVERGENCES ARE FROM A RELEASE. `modules/soil.cpp`'s
`update_snow_properties` arrived byte-identical to `guess_4.1/modules/soil.cpp`'s
apart from a stripped licence header, and what it ran is default-on:
`data/ins/global.ins` sets `iftwolayersoil 0`, which selects the multilayer soil
temperature scheme, and `Ksnow` and `Csnow` become the conductivity and the
volumetric heat capacity of every active snow layer in that scheme's numerical
solve. So every change here is to default-on behaviour in a widely used
community model, and each is recorded on the same terms `ntransform_gate.py`
records its own: mainline's lines stand verbatim beside the changed one.

It can fail:

  relation     a literal restatement of Fourteau Eq. (18) that no longer carries
               the adopted row's coefficients, or normalises the ice volume
               fraction by something other than the density the fit was made
               against, or has stopped being a quadratic this check can read
  divergence   a declared divergence from the release whose mainline form the
               model no longer records beside the changed one, which it has gone
               back to running, or whose changed line it no longer contains. All
               three are checked on every entry, because a divergence that is not
               recorded is a silent fork, one that is recorded but live is a
               declaration that has drifted from the model, and one whose changed
               line is gone has been reverted with the record left behind
  drift        a density span the declaration names and `modules/soil.h` does
               not declare, or a coefficient of the superseded relation live
               anywhere in the vendored vegetation model

Reduced fixtures run on every invocation, all but one built to be wrong in a
named way. A fixture that does not get the verdict it was built for is a defect
in this checker rather than in the declaration.

    python biosphere/scripts/snow_thermal_gate.py            # status, exit 0
    python biosphere/scripts/snow_thermal_gate.py --strict   # refuses on the
                                                             # named residual

`--strict` refuses on exactly one thing: a divergence whose verdict is `gate`.
There is none; the verdict exists so the next one has somewhere to go that is
not `keep`. It does NOT refuse on the kinetics bracket, because that is not a
residual this project can close: Fourteau says which limit snow is in is
unresolved, the adopted row is the bracket's upper endpoint, and what follows is
that the modelled snow may conduct LESS than the model says and cannot conduct
more. A one-signed open question that the literature holds open is reported and
not gated.

EXECUTION IS A PROPERTY OF EACH DIVERGENCE AND NEVER OF THE REGISTER. The
declaration carries the matched arms that have run, by name, and every entry
claims one of them or says `none` and why. No LPJ-GUESS run this project has
made has passed acceptance, so that block is empty and every statement here is
against the source, the declaration and the papers; an entry appended later
therefore has to state its own standing rather than inherit an answer written
before it existed. The gate refuses an entry with no execution claim, an entry
claiming an arm the register does not carry, an unexecuted entry that says
nothing about why none covers it, and an arm no entry claims. The declaration
names the comparison arm that would price this one in a simulated soil
temperature.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import GENERATED, GUESS_SOURCE, PROJECT_ROOT
from paths import rel  # noqa: E402

import snow  # noqa: E402  -- lib/, put on the path by _paths

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = COMPONENT_ROOT / "config" / "snow_thermal.yaml"
REPORT = GENERATED / "snow_thermal_gate_report.json"

# What a declared divergence may say about itself. `keep` is one this project
# reviewed and stood behind, declared at its own site in the source; `gate` is
# one that stands but is not settled, and `--strict` refuses while it does. A
# divergence that does not survive review is reverted and leaves the
# declaration, so there is no `revert`.
DIVERGENCE_VERDICTS = ("keep", "gate")

# What shape the release's own version has. `line` is the only one used here:
# the release runs a DIFFERENT line, so all three checks bite. The
# `commented_out` kind is carried so the register keeps the same vocabulary as
# `biosphere/config/somdynam.yaml`, where the weaker guard is argued.
MAINLINE_FORMS = ("line", "commented_out")

# Which configuration the divergence is live in. `iftwolayersoil_0` is the
# default and what this project runs, and is the only one in which the
# multilayer scheme's snow layers exist at all.
LIVE_UNDER = ("iftwolayersoil_0", "both", "neither")

# Every C++ translation unit and header of the vendored vegetation model, which
# is where the superseded relation must not reappear. Scanned whole rather than
# one function at a time: a second snow conductivity growing somewhere else is
# exactly what a one-function check would miss.
GUESS_SUFFIXES = (".cpp", ".h")


def _strip_comments(text: str) -> str:
    """The source with C and C++ comments removed.

    The divergence check needs to tell a form RECORDED in a comment from the
    same form LIVE in the code, which is the whole point of recording it: the
    release's lines have to stand beside the changed one without being what the
    model runs.
    """
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", " ", text)


def _guess_sources(root: Path) -> dict[str, str]:
    """Every vendored vegetation-model source, stripped of comments.

    Keyed by repository-relative path so a finding names a file someone can
    open. Build directories are excluded because a generated copy is not a
    second statement of anything.
    """
    out: dict[str, str] = {}
    for path in sorted(root.rglob("*")):
        if path.suffix not in GUESS_SUFFIXES or not path.is_file():
            continue
        if "build" in path.relative_to(root).parts:
            continue
        out[str(path.relative_to(PROJECT_ROOT))] = _strip_comments(
            path.read_text(encoding="utf-8", errors="replace"))
    return out


def _declared_constant(text: str, name: str) -> float | None:
    """One LPJ-GUESS `const double`, read from the source that declares it."""
    match = re.search(rf"^\s*const\s+double\s+{re.escape(name)}\s*=\s*([-\d.eE+]+)",
                      text, re.M)
    return None if match is None else float(match.group(1))


def _check_divergences(declaration: dict, sources: dict) -> list[dict]:
    """The declared divergence from the release, against its source.

    Three checks, and a divergence needs all three. Each `mainline` string has
    to be IN the source it names, so the record stands beside the changed code
    and cannot go missing; NOT in that source once comments are stripped, so a
    silent revert to the release fails; and `live` has to be in the stripped
    source, so deleting the changed line fails as well. Without the third,
    reverting by deleting both would pass.
    """
    findings: list[dict] = []
    register = declaration.get("mainline_divergences")
    if not register:
        findings.append({"kind": "divergence", "what": "mainline_divergences",
                         "detail": "the declaration carries no divergence register"})
        return findings

    def bad(what: str, detail: str) -> None:
        findings.append({"kind": "divergence", "what": what, "detail": detail})

    for field in ("release", "fork_reference", "comparison_arm"):
        if not register.get(field):
            bad("mainline_divergences", f"the register carries no {field}")
    release = register.get("release") or ""
    if not re.search(r"\b\d{6,}\b", release):
        bad("mainline_divergences",
            "the register's release names no Zenodo record. The operator arrived "
            "byte-identical to a release, so that release is what a divergence "
            "here is measured against and it has to be nameable")
    arms = register.get("execution_arms")
    if arms is None or not isinstance(arms, dict):
        bad("mainline_divergences",
            "execution_arms is not a mapping of arm name to its evidence. No "
            "matched arm has run, so it is empty rather than absent: an absent "
            "block cannot be told from one that was never written")
        arms = {}
    for name, evidence in arms.items():
        where = f"execution_arms.{name}"
        if not isinstance(evidence, dict):
            bad(where, "an execution arm that is not a mapping of evidence")
            continue
        for field in ("report", "standing"):
            if not evidence.get(field):
                bad(where, f"an execution arm carrying no {field}")
        report = evidence.get("report")
        if report and not (PROJECT_ROOT / report).is_file():
            bad(where, f"execution report {report} is missing")

    seen = set()
    claimed: set[str] = set()
    for entry in register.get("entries", []):
        what = entry.get("id", "?")
        if what in seen:
            bad(what, "two register entries share one id")
        seen.add(what)
        if entry.get("verdict") not in DIVERGENCE_VERDICTS:
            bad(what, f"unknown divergence verdict {entry.get('verdict')!r}")
        if not entry.get("owner"):
            bad(what, "a divergence naming no owner")
        if entry.get("live_under") not in LIVE_UNDER:
            bad(what, (f"unknown live_under {entry.get('live_under')!r}. A "
                       "divergence has to say which configuration it is live in, "
                       "because one that is inert in the arm this project runs is "
                       "waiting rather than harmless"))
        for field in ("settles", "worth"):
            if not entry.get(field):
                bad(what, f"a divergence saying nothing about what it {field}")
        claim = entry.get("execution")
        if not claim:
            bad(what, "a divergence saying nothing about whether it has been "
                      "executed. Each entry names the arm that covers it or "
                      "says `none`, because one answer over the whole register "
                      "stops being true the moment an entry is appended")
        elif claim == "none":
            if not entry.get("why_no_execution"):
                bad(what, "a divergence claiming no execution arm and saying "
                          "nothing about why none covers it")
        elif claim not in arms:
            bad(what, f"claims execution arm {claim!r}, which the register does "
                      "not carry")
        else:
            claimed.add(claim)

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
            bad(what, "a divergence naming no form the release runs")
        for form in mainline:
            if form not in source_text:
                bad(what, (f"declares a divergence from {form!r}, and "
                           f"{source_file} does not record that form beside the "
                           "changed one"))
                continue
            if form_kind == "commented_out":
                if not form.lstrip().startswith("//"):
                    bad(what, (f"declares {form!r} as a line the release leaves "
                               "commented out, and it carries no comment marker"))
            elif form in code_only:
                bad(what, (f"declares a divergence from {form!r}, and "
                           f"{source_file} still runs it"))
        live = entry.get("live")
        if not live:
            bad(what, "a divergence naming no line the model runs instead")
        elif live not in code_only:
            bad(what, (f"declares that {source_file} runs {live!r} instead, "
                       "and it does not"))

    for name in sorted(set(arms) - claimed):
        bad("mainline_divergences",
            f"execution arm {name!r} is claimed by no divergence, so it is "
            "evidence that has outlived the entries it was taken for")
    return findings


def _check_relation(declaration: dict, sources: dict, guess_sources: dict,
                    root: Path) -> list[dict]:
    """The relation, against `lib/snow.py` and against the source's own constants.

    Three things. Every literal restatement in the tree is held to the one table
    by the module that owns it. The density span the register names is read out
    of `modules/soil.h` rather than trusted, because a compaction ramp that
    moves and a register that does not is the same drift in a different place.
    And the superseded relation's own coefficients are looked for over the whole
    vendored model with comments stripped, so a second snow conductivity growing
    outside the one function fails here too.
    """
    findings: list[dict] = []
    block = declaration.get("relation")
    if not block:
        return [{"kind": "drift", "what": "relation",
                 "detail": "the declaration says nothing about the relation"}]

    def bad(kind: str, what: str, detail: str) -> None:
        findings.append({"kind": kind, "what": what, "detail": detail})

    declared_in = block.get("declared_in")
    if declared_in != "lib/snow.py":
        bad("drift", "relation",
            f"names {declared_in!r} as where the relation is declared, and this "
            "gate checks lib/snow.py")
    for field in ("citation", "bracket"):
        if not block.get(field):
            bad("drift", "relation", f"the relation block carries no {field}")

    for problem in snow.check_restatements(root):
        bad("relation", "restatement", problem)

    densities = block.get("densities") or {}
    source_file = densities.get("source_file")
    if source_file not in sources:
        bad("drift", "densities",
            f"names {source_file!r}, which this gate does not read")
    else:
        text = sources[source_file][1]
        span = {}
        for end in ("start", "end"):
            name = densities.get(end)
            value = None if name is None else _declared_constant(text, name)
            if value is None:
                bad("drift", "densities",
                    f"names {name!r} as the {end} of the compaction ramp, and "
                    f"{source_file} does not declare it")
            span[end] = value
        if span.get("start") is not None and span.get("end") is not None:
            if not span["start"] < span["end"]:
                bad("drift", "densities",
                    f"the compaction ramp runs from {span['start']} to "
                    f"{span['end']}, which is not a ramp")

    forbidden = (block.get("no_other_relation") or {}).get("forbidden") or []
    if not forbidden:
        bad("drift", "no_other_relation",
            "the declaration forbids no coefficient of the superseded relation, "
            "so nothing stops a second snow conductivity growing elsewhere")
    for fragment in forbidden:
        for name, code_only in guess_sources.items():
            if fragment in code_only:
                bad("drift", "no_other_relation",
                    f"{name} runs {fragment!r}, which is the superseded "
                    "relation's own coefficient")
    return findings


def check(declaration: dict, sources: dict, guess_sources: dict,
          root: Path) -> list[dict]:
    """Every check, as a list of findings. An empty list is a clean gate."""
    return (_check_divergences(declaration, sources)
            + _check_relation(declaration, sources, guess_sources, root))


# The fixtures. The first is the declaration as it stands and has to come back
# clean; every other is built to be wrong in a named way. A fixture that does
# not get its verdict is a defect in the checker.
def _fixtures(declaration: dict, sources: dict, guess_sources: dict,
              root: Path) -> list[dict]:
    import copy

    def mutate(fn):
        d = copy.deepcopy(declaration)
        fn(d)
        return d

    def entry(d, entry_id):
        for item in d["mainline_divergences"]["entries"]:
            if item["id"] == entry_id:
                return item
        raise KeyError(entry_id)

    def claim(entry_id, field, value):
        return lambda d: entry(d, entry_id).__setitem__(field, value)

    cases = [
        ("the declaration as it stands", declaration, None),
        ("a divergence whose form the source does not record",
         mutate(claim("snow_conductivity_relation", "mainline",
                      ["Ksnow = no_such_relation(snowdens);"])),
         "divergence"),
        ("a divergence from a form the source in fact still runs",
         mutate(claim("snow_conductivity_relation", "mainline",
                      ["Dsnow = Ksnow / Csnow * SECS_PER_DAY * MM2_PER_M2;"])),
         "divergence"),
        ("a divergence whose changed line the source no longer contains",
         mutate(claim("snow_conductivity_relation", "live",
                      "Ksnow = 1.0 * snowdens_vf;")),
         "divergence"),
        ("a divergence not saying which configuration it is live in",
         mutate(claim("snow_conductivity_relation", "live_under", "sometimes")),
         "divergence"),
        ("a second divergence appended, saying nothing about whether an arm "
         "has executed it",
         mutate(lambda d: d["mainline_divergences"]["entries"].append(
             {k: v for k, v in copy.deepcopy(
                 entry(d, "snow_conductivity_relation")).items()
              if k != "execution"} | {"id": "a_second_divergence"})),
         "divergence"),
        ("a divergence claiming an execution arm the register does not carry",
         mutate(claim("snow_conductivity_relation", "execution", "no_such_arm")),
         "divergence"),
        ("an unexecuted divergence saying nothing about why no arm covers it",
         mutate(lambda d: entry(d, "snow_conductivity_relation")
                .pop("why_no_execution")),
         "divergence"),
        ("an execution arm no divergence claims",
         mutate(lambda d: d["mainline_divergences"]["execution_arms"]
                .__setitem__("an_arm_nobody_claims",
                             {"report": "biosphere/config/snow_thermal.yaml",
                              "standing": "a fixture"})),
         "divergence"),
        ("an execution arm whose report is missing",
         mutate(lambda d: d["mainline_divergences"]["execution_arms"]
                .__setitem__("an_arm_with_no_report",
                             {"report": "biosphere/analysis/no_such.json",
                              "standing": "a fixture"})),
         "divergence"),
        ("an execution_arms block absent rather than empty",
         mutate(lambda d: d["mainline_divergences"].pop("execution_arms")),
         "divergence"),
        ("a register whose release names no Zenodo record",
         mutate(lambda d: d["mainline_divergences"].__setitem__(
             "release", "mainline LPJ-GUESS")),
         "divergence"),
        ("the SECOND divergence's changed line, which the source no longer "
         "contains. The entries are checked one by one, so an appended one is "
         "covered by the loop and not by the first entry's verdict",
         mutate(claim("snow_heat_capacity_density", "live",
                      "Csnow *= ice_density;")),
         "divergence"),
        ("a relation declared somewhere other than the module that owns it",
         mutate(lambda d: d["relation"].__setitem__("declared_in",
                                                    "biosphere/config/snow_thermal.yaml")),
         "drift"),
        ("a compaction ramp whose ends the source does not declare",
         mutate(lambda d: d["relation"]["densities"].__setitem__(
             "start", "snowdens_first")),
         "drift"),
        ("a declaration that forbids no coefficient of the superseded relation",
         mutate(lambda d: d["relation"]["no_other_relation"].__setitem__(
             "forbidden", [])),
         "drift"),
        ("a coefficient the model does in fact still run, declared forbidden",
         mutate(lambda d: d["relation"]["no_other_relation"].__setitem__(
             "forbidden", ["Ksnow = 1.985 * snowdens_vf"])),
         "drift"),
    ]

    results = []
    for label, candidate, expect in cases:
        findings = check(candidate, sources, guess_sources, root)
        kinds = {f["kind"] for f in findings}
        ok = (not findings) if expect is None else (expect in kinds)
        results.append({"fixture": label, "expected": expect or "clean",
                        "found": sorted(kinds), "pass": ok})
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict", action="store_true",
                        help="refuse while a divergence is gated rather than kept")
    parser.add_argument("--json", action="store_true", help="report only")
    args = parser.parse_args()

    declaration = yaml.safe_load(DECLARATION.read_text())

    wanted = {declaration["relation"]["densities"]["source_file"]}
    wanted |= {item["source_file"]
               for item in declaration["mainline_divergences"]["entries"]
               if item.get("source_file")}
    sources: dict[str, tuple[str, str]] = {}
    for name in sorted(wanted):
        text = (PROJECT_ROOT / name).read_text()
        sources[name] = (text, _strip_comments(text))
    guess_sources = _guess_sources(GUESS_SOURCE)

    fixtures = _fixtures(declaration, sources, guess_sources, PROJECT_ROOT)
    findings = check(declaration, sources, guess_sources, PROJECT_ROOT)

    register = declaration["mainline_divergences"]
    divergences = [
        {"id": item.get("id", "?"), "verdict": item.get("verdict", "?"),
         "owner": item.get("owner", "?"),
         "source_file": item.get("source_file", "?"),
         "live_under": item.get("live_under", "?"),
         "execution": item.get("execution", "?")}
        for item in register.get("entries", [])
    ]
    arms = register.get("execution_arms") or {}
    executed = [d for d in divergences if d["execution"] in arms]
    gated = [d for d in divergences if d["verdict"] == "gate"]

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": rel(DECLARATION),
        "relation_declared_in": declaration["relation"]["declared_in"],
        "sources": sorted(sources),
        "vendored_sources_scanned": len(guess_sources),
        "findings": findings,
        "release": register.get("release"),
        "divergences": divergences,
        "execution_arms": sorted(arms),
        "divergences_executed": [d["id"] for d in executed],
        "divergences_unexecuted": [d["id"] for d in divergences
                                   if d["execution"] not in arms],
        "fixtures": fixtures,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print("The simulated snowpack's thermal properties, checked against the source.\n")
        broken = [f for f in fixtures if not f["pass"]]
        print(f"  fixtures: {len(fixtures) - len(broken)} of {len(fixtures)} got their verdict")
        for case in broken:
            print(f"    BROKEN: {case['fixture']}: expected {case['expected']}, found {case['found']}")
        if findings:
            print(f"\n  {len(findings)} finding(s):")
            for finding in findings:
                print(f"    [{finding['kind']}] {finding['what']}: {finding['detail']}")
        else:
            print("  the relation is one declaration and both models restate it")
        print(f"\n  the relation is declared in {declaration['relation']['declared_in']}")
        print(f"  and restated in {len(snow.RESTATEMENTS)} compiled model(s), each held to it:")
        for item in snow.RESTATEMENTS:
            print(f"    {item['path']}  ({item['what']})")
        print("\n  IT IS THE BRACKET'S UPPER ENDPOINT. The modelled snow may conduct")
        print("  LESS than this and it cannot conduct more, which is not a residual")
        print("  --strict can refuse: the literature holds it open.")
        kept = [d for d in divergences if d["verdict"] == "keep"]
        print(f"\n  {len(divergences)} declared divergence(s) from the release,")
        print(f"  {len(kept)} kept and {len(gated)} gated. The operator arrived")
        print("  byte-identical to it and is ON BY DEFAULT:")
        print(f"    {register.get('release')}")
        print(f"  {len(executed)} of {len(divergences)} claim a matched execution "
              "arm, per entry and not per")
        print("  register. Each says for itself which arm covers it, or why none does.")
        for item in divergences:
            print(f"    [{item['verdict']}] {item['id']}  live under "
                  f"{item['live_under']}  ({item['source_file']})  [{item['owner']}]"
                  f"  execution: {item['execution']}")
        print(f"\n  report: {rel(REPORT)}")

    if broken := [f for f in fixtures if not f["pass"]]:
        print(f"\n{len(broken)} fixture(s) did not get the verdict they were built for.",
              file=sys.stderr)
        print("That is a defect in this checker, not in the declaration.", file=sys.stderr)
        return 2
    if findings:
        return 1
    if args.strict and gated:
        print("\n--strict: refused, on exactly what is still unsettled.", file=sys.stderr)
        print("  divergence(s) gated: " + ", ".join(item["id"] for item in gated),
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
