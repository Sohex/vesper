"""The fire gate: what this project's fire operators are declared to depart from,
and whether the model still departs from it.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a vegetation model's fire operators, the file that declares
their departures from mainline, and the source the model actually reads.

`biosphere/config/fire.yaml` is the declaration and this module is the
enforcement. THE SCOPE IS THE FIRE OPERATORS AND NOT ONE SOURCE FILE, which is
the one structural difference from `ntransform_gate.py` and `somdynam_gate.py`:
GLOBFIRM's burned fraction, its mortality and the BLAZE effects path live in
different translation units, so each entry names its own `source_file` and this
gate reads the union of them.

THE DIVERGENCE MAY BE A DELETION. `mainline_form: deleted_line` is the shape the
other three gates have no case for: a line the release runs and this project's
source does not run at all. There is no changed line to find, so the third check
those gates make -- `live` is in the stripped source -- has nothing to look at,
and its work is done by `absent_from_stripped` instead. That block carries the
fragments a REFORMATTED revert would reintroduce, which is what makes the
absence a check rather than a restatement of the second one: the recorded
mainline form carries its own `//` here, so searching the stripped source for it
verbatim can never fail and would be a check in name only. Both halves are run:
the form with its own comment removed, and every declared fragment.

It can fail:

  divergence   a declared divergence whose mainline form the model no longer
               records beside the change, which it has gone back to running, or
               -- for a deletion -- whose forbidden fragments have reappeared in
               code. All are checked, because a divergence that is not recorded
               is a silent fork, one that is recorded but live is a declaration
               that has drifted from the model, and a deletion whose fragments
               are back has been reverted with the record left behind
  execution    an entry that says nothing about whether an arm has run it, an
               entry claiming an arm the register does not carry, an unexecuted
               entry that says nothing about why none covers it, or an arm no
               entry claims

Reduced fixtures run on every invocation, all but one built to be wrong in a
named way. A fixture that does not get the verdict it was built for is a defect
in this checker rather than in the declaration.

    python biosphere/scripts/fire_gate.py            # status, exit 0
    python biosphere/scripts/fire_gate.py --strict   # refuses on the named
                                                     # residual

`--strict` refuses on exactly one thing: a divergence whose verdict is `gate`.
There is none; the verdict exists so the next one has somewhere to go that is
not `keep`.

EXECUTION IS A PROPERTY OF EACH DIVERGENCE AND NEVER OF THE REGISTER. The
declaration carries the arms that have run, by name, and every entry claims one
of them or says `none` and why. No LPJ-GUESS run this project has made has
passed acceptance, so that block is empty here and every entry says so for
itself; an entry appended later has to state its own standing rather than
inherit an answer written before it existed.
"""

from __future__ import annotations

import argparse
import copy
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import GENERATED, PROJECT_ROOT
from paths import rel  # noqa: E402

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = COMPONENT_ROOT / "config" / "fire.yaml"
REPORT = GENERATED / "fire_gate_report.json"

# What a declared divergence may say about itself. `keep` is one this project
# reviewed and stood behind, declared at its own site in the source; `gate` is
# one that stands but is not settled, and `--strict` refuses while it does. A
# divergence that does not survive review is reverted and leaves the
# declaration, so there is no `revert`.
DIVERGENCE_VERDICTS = ("keep", "gate")

# What shape the release's own version has, and the whole reason this gate is
# not `snow_thermal_gate.py` with a different declaration path:
#   line            the release runs a different line, and this project's source
#                   runs its own instead. All three checks bite
#   commented_out   the release's version is a line it leaves commented out, so
#                   the second check is vacuous and is declared vacuous
#   deleted_line    the release runs a line this project's source does not run
#                   at all. There is no `live` line, and the absence is what
#                   `absent_from_stripped` carries
#   added_line      this project's source runs a line the release has NOTHING
#                   in place of. The inverse of a deletion, and its checks are
#                   the inverse too: there is no released form to find absent
#                   from code, so what can fail is `live` going missing and the
#                   recorded statement of what the release does instead going
#                   missing. Both are asserted; the third check the `line` form
#                   makes is declared vacuous here rather than run, because a
#                   release that runs nothing leaves nothing to forbid
MAINLINE_FORMS = ("line", "commented_out", "deleted_line", "added_line")

# Where a divergence sits. `operator` is inside a named function and is the
# ordinary case; `file_scope` is a declaration or a whole function outside any
# other, which the function check cannot be run on -- a DELETED function has no
# enclosing function to name, and a file-scope constant never had one.
DIVERGENCE_SITES = ("operator", "file_scope")

# Which configuration a divergence is live in. Two switches select a fire path
# -- `firemodel` and `vegmode` -- so an entry names the pair it bites under
# rather than a boolean. `neither` is an entry the running configuration never
# reaches, which is waiting rather than harmless and has to say so.
LIVE_UNDER = ("globfirm_cohort", "globfirm_individual", "globfirm_any",
              "blaze", "neither")

# How the constant a divergence disposes of lacked a usable derivation, in
# CLAUDE.md's vocabulary. Three defects with three repairs, and the whole point
# of naming one is that they are not interchangeable: IMPLICIT-EARTH is a
# diagnosis that a sound derivation is for the wrong planet, never an
# endorsement of the number. `not_a_constant` is a divergence that is about an
# operator's form rather than about a number.
CLASSIFICATIONS = ("tuned", "opaque", "implicit-earth", "not_a_constant")


def _strip_comments(text: str) -> str:
    """The source with C and C++ comments removed.

    The divergence check needs to tell a form RECORDED in a comment from the
    same form LIVE in the code, which is the whole point of recording it: the
    release's line has to stand beside the change without being what the model
    runs.
    """
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", " ", text)


def _check_arms(register: dict, bad) -> dict:
    """The register's execution arms, each against the evidence it names.

    An arm is a matched pair of runs and a report. It is refused when it carries
    no report, when the report it names is not on disk, or when it says nothing
    about the standing of that report -- a matched arm between two runs that
    failed acceptance bounds the entries naming it against each other and is not
    a statement that either arm is a Vesper result, and an arm that does not say
    so reads as one that is.
    """
    arms = register.get("execution_arms")
    if arms is None or not isinstance(arms, dict):
        bad("mainline_divergences",
            "execution_arms is not a mapping of arm name to its evidence. No "
            "arm has run, so it is empty rather than absent: an absent block "
            "cannot be told from one that was never written")
        return {}
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
    return arms


def _check_entry_execution(entry: dict, what: str, arms: dict, bad,
                           claimed: set[str]) -> None:
    """One entry's execution claim, against the register's arms.

    Per entry and never per register: one boolean over a set that is appended to
    goes on answering for entries added after it was written, which is what
    `world-h3qg` removed from the other registers.
    """
    claim = entry.get("execution")
    if not claim:
        bad(what, "a divergence saying nothing about whether it has been "
                  "executed. Each entry names the arm that covers it or says "
                  "`none`, because one answer over the whole register stops "
                  "being true the moment an entry is appended")
    elif claim == "none":
        if not entry.get("why_no_execution"):
            bad(what, "a divergence claiming no execution arm and saying "
                      "nothing about why none covers it")
    elif claim not in arms:
        bad(what, f"claims execution arm {claim!r}, which the register does "
                  "not carry")
    else:
        claimed.add(claim)


def _check_forms(entry: dict, what: str, form_kind: str, source_file: str,
                 source_text: str, code_only: str, bad) -> None:
    """The mainline record and the change, against the source.

    Two assertions on every entry and a third that depends on the shape. Each
    `mainline` string has to be IN the source it names, so the record stands
    beside the change and cannot go missing; and NOT in that source once
    comments are stripped, so a silent revert fails.

    The second assertion is run on the form with its OWN comment removed as well
    as verbatim. A recorded mainline line that carries a trailing `// ...` can
    never survive comment stripping, so searching for it verbatim is a check
    that cannot fail; the bare code half of it can, and is what a revert would
    put back.
    """
    mainline = entry.get("mainline") or []
    if not mainline:
        bad(what, "a divergence naming no form the release runs")
    for form in mainline:
        if form not in source_text:
            bad(what, (f"declares a divergence from {form!r}, and "
                       f"{source_file} does not record that form beside the "
                       "change"))
            continue
        if form_kind == "commented_out":
            if not form.lstrip().startswith("//"):
                bad(what, (f"declares {form!r} as a line the release leaves "
                           "commented out, and it carries no comment marker. "
                           "The weaker guard cannot be claimed for a line that "
                           "would support the stronger one"))
            continue
        if form_kind == "added_line":
            # The release runs NOTHING here, so `mainline` is a recorded
            # STATEMENT of that and not a line of code. Asserting it is absent
            # from the stripped source would be a check that cannot fail, which
            # is the shape this module already refuses elsewhere. What can fail
            # is the statement going missing, which the loop above has just
            # checked, and `live` going missing, which the caller checks.
            continue
        for candidate in {form, _strip_comments(form).strip()}:
            if candidate and candidate in code_only:
                bad(what, (f"declares a divergence from {candidate!r}, and "
                           f"{source_file} still runs it"))


def _check_absence(entry: dict, what: str, source_file: str, code_only: str,
                   bad) -> None:
    """A deletion's forbidden fragments, against the stripped source.

    The check the other three gates have no case for. A deletion has no changed
    line, so nothing in the source points at it and a revert can be written in
    any spacing the compiler accepts. `absent_from_stripped` carries the
    fragments that would have to reappear for the release's behaviour to be
    back, and every one of them has to be absent from code.
    """
    fragments = entry.get("absent_from_stripped") or []
    if not fragments:
        bad(what, "a deletion carrying no absent_from_stripped, so a revert "
                  "written in different spacing would pass. A deleted line has "
                  "no changed line to point at and the absence is the whole "
                  "check")
    for fragment in fragments:
        if fragment in code_only:
            bad(what, (f"declares a deletion, and {source_file} runs "
                       f"{fragment!r}"))
    if entry.get("live"):
        bad(what, (f"declares a deletion and names {entry['live']!r} as a line "
                   "the model runs instead. A deletion has no changed line; an "
                   "entry with one is a `line` divergence mis-declared"))


def _check_divergences(declaration: dict, sources: dict) -> list[dict]:
    """Every declared divergence from the release, against its source."""
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
            "the register's release names no Zenodo record. The fire operator "
            "arrived byte-identical to a release apart from the fork's currency "
            "rename, so that release is what a divergence here is measured "
            "against and it has to be nameable")
    arms = _check_arms(register, bad)

    seen: set[str] = set()
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
                       "divergence has to say which fire configuration it is "
                       "live in, because one the running configuration never "
                       "reaches is waiting rather than harmless"))
        if entry.get("classification") not in CLASSIFICATIONS:
            bad(what, (f"unknown classification {entry.get('classification')!r}. "
                       "tuned, opaque and implicit-earth are three defects with "
                       "three repairs and are never synonyms; a divergence about "
                       "an operator's form rather than a number says so"))
        for field in ("settles", "worth", "breaks_if_removed"):
            if not entry.get(field):
                bad(what, f"a divergence saying nothing about what it {field}")
        _check_entry_execution(entry, what, arms, bad, claimed)

        if entry.get("where") not in DIVERGENCE_SITES:
            bad(what, f"unknown divergence site {entry.get('where')!r}")
            continue
        source_file = entry.get("source_file")
        if source_file not in sources:
            bad(what, f"names {source_file!r}, which this gate does not read")
            continue
        source_text, code_only = sources[source_file]

        function = entry.get("function")
        if entry.get("where") == "file_scope":
            if function:
                bad(what, (f"is declared at file scope and names the function "
                           f"{function!r}. A file-scope divergence sits in no "
                           f"function; one that does is an `operator` "
                           f"divergence mis-declared"))
        elif not function:
            bad(what, "a divergence naming no function it sits in")
        elif not re.search(rf"\b{re.escape(function)}\s*\(", code_only):
            bad(what, (f"sits in {function!r}, and {source_file} no longer "
                       "runs a function of that name"))

        form_kind = entry.get("mainline_form")
        if form_kind not in MAINLINE_FORMS:
            bad(what, f"unknown mainline_form {entry.get('mainline_form')!r}")
            continue
        _check_forms(entry, what, form_kind, source_file, source_text,
                     code_only, bad)

        if form_kind == "deleted_line":
            _check_absence(entry, what, source_file, code_only, bad)
        else:
            if entry.get("absent_from_stripped"):
                bad(what, (f"is a {form_kind!r} divergence carrying "
                           "absent_from_stripped, which only a deletion has. "
                           "Nothing was removed, so there is no fragment whose "
                           "reappearance would be a revert"))
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


def check(declaration: dict, sources: dict) -> list[dict]:
    """Every check, as a list of findings. An empty list is a clean gate."""
    return _check_divergences(declaration, sources)


def source_files(declaration: dict) -> set[str]:
    """Every source file the declaration names. The scope is the operators."""
    return {entry["source_file"]
            for entry in declaration["mainline_divergences"].get("entries", [])
            if entry.get("source_file")}


def read_sources(declaration: dict, root: Path) -> dict[str, tuple[str, str]]:
    """Each named source, raw and with comments stripped."""
    out: dict[str, tuple[str, str]] = {}
    for name in sorted(source_files(declaration)):
        text = (root / name).read_text(encoding="utf-8", errors="replace")
        out[name] = (text, _strip_comments(text))
    return out


# The fixtures. The first is the declaration as it stands and has to come back
# clean; every other is built to be wrong in a named way. A fixture that does
# not get its verdict is a defect in the checker.
def _fixtures(declaration: dict, sources: dict) -> list[dict]:

    def mutate(fn):
        d = copy.deepcopy(declaration)
        fn(d)
        return d

    def _entry(d, entry_id):
        for item in d["mainline_divergences"]["entries"]:
            if item["id"] == entry_id:
                return item
        raise KeyError(entry_id)

    def claim(entry_id, field, value):
        return lambda d: _entry(d, entry_id).__setitem__(field, value)

    # NAMED, not indexed and not matched by shape. The fixtures below mutate one
    # entry and assert what the checker then says, and most of them use
    # fragments of `vegdynam.cpp` -- so they need THE FLOOR entry specifically,
    # not merely any deletion. Indexing entries[0] stopped testing them the
    # moment an entry was declared above it, and matching on
    # `mainline_form == "deleted_line"` stopped the moment a second deletion in
    # another file was declared above it. Both failures were silent: the
    # fixtures went on passing while asserting nothing.
    FIXTURE_ENTRY = "globfirm_fireprob_floor"
    ids = {e["id"] for e in declaration["mainline_divergences"]["entries"]}
    if FIXTURE_ENTRY not in ids:
        return [{"fixture": f"the fixtures mutate {FIXTURE_ENTRY!r}",
                 "expected": "that entry to exist", "found": sorted(ids),
                 "pass": False}]
    first = FIXTURE_ENTRY

    # A line the source in fact still runs, and a fragment of one. Both are read
    # out of `fire()` rather than invented, so a fixture that stops being wrong
    # is a fixture whose source moved.
    live_line = "fireprob=s*exp(sm/(0.45*sm*sm*sm+2.83*sm*sm+2.96*sm+1.04));"
    live_fragment = "fireprob=0.0"

    cases = [
        ("the declaration as it stands", declaration, None),
        ("a divergence whose form the source does not record",
         mutate(claim(first, "mainline", ["fireprob=no_such_floor();"])),
         "divergence"),
        ("a divergence from a form the source in fact still runs",
         mutate(claim(first, "mainline", [live_line])),
         "divergence"),
        ("a deletion whose forbidden fragment the source runs",
         mutate(claim(first, "absent_from_stripped", [live_fragment])),
         "divergence"),
        ("a deletion carrying no forbidden fragments, so a reformatted revert "
         "would pass",
         mutate(claim(first, "absent_from_stripped", [])),
         "divergence"),
        # The added_line form and the file_scope site, each able to fail in a
        # named way. They are built by MUTATING the floor entry rather than by
        # writing a synthetic one, so a fixture cannot drift away from the
        # schema the real entries use.
        ("an added line the source does not run",
         mutate(lambda d: (_entry(d, first).update(
             {"mainline_form": "added_line", "absent_from_stripped": None,
              "live": "fireprob = no_such_conversion();"}))),
         "divergence"),
        ("an added line naming no line the model runs",
         mutate(lambda d: (_entry(d, first).update(
             {"mainline_form": "added_line", "absent_from_stripped": None,
              "live": None}))),
         "divergence"),
        ("an added line carrying absent_from_stripped, which only a deletion has",
         mutate(lambda d: (_entry(d, first).update(
             {"mainline_form": "added_line",
              "absent_from_stripped": ["fireprob=0.001"],
              "live": live_line}))),
         "divergence"),
        ("a file-scope divergence that also names a function",
         mutate(claim(first, "where", "file_scope")),
         "divergence"),
        ("a deletion that also names a line the model runs instead",
         mutate(claim(first, "live", live_fragment + ";")),
         "divergence"),
        ("a divergence declared in a shape this gate has no case for",
         mutate(claim(first, "mainline_form", "reworded")),
         "divergence"),
        ("a divergence not saying which fire configuration it is live in",
         mutate(claim(first, "live_under", "sometimes")),
         "divergence"),
        ("a divergence classifying its constant as none of the three defects",
         mutate(claim(first, "classification", "earthlike")),
         "divergence"),
        ("a divergence sitting in a function the source does not run",
         mutate(claim(first, "function", "no_such_function")),
         "divergence"),
        ("a divergence saying nothing about what it breaks_if_removed",
         mutate(lambda d: _entry(d, first).pop("breaks_if_removed")),
         "divergence"),
        ("a second divergence appended, saying nothing about whether an arm "
         "has executed it",
         mutate(lambda d: d["mainline_divergences"]["entries"].append(
             {k: v for k, v in copy.deepcopy(_entry(d, first)).items()
              if k != "execution"} | {"id": "a_second_divergence"})),
         "divergence"),
        ("a divergence claiming an execution arm the register does not carry",
         mutate(claim(first, "execution", "no_such_arm")),
         "divergence"),
        ("an unexecuted divergence saying nothing about why no arm covers it",
         mutate(lambda d: _entry(d, first).pop("why_no_execution")),
         "divergence"),
        ("an execution arm no divergence claims",
         mutate(lambda d: d["mainline_divergences"]["execution_arms"]
                .__setitem__("an_arm_nobody_claims",
                             {"report": "biosphere/config/fire.yaml",
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
                        help="refuse while a divergence is gated rather than kept")
    parser.add_argument("--json", action="store_true", help="report only")
    args = parser.parse_args()

    declaration = yaml.safe_load(DECLARATION.read_text(encoding="utf-8"))
    sources = read_sources(declaration, PROJECT_ROOT)

    fixtures = _fixtures(declaration, sources)
    findings = check(declaration, sources)

    register = declaration["mainline_divergences"]
    divergences = [
        {"id": item.get("id", "?"), "verdict": item.get("verdict", "?"),
         "owner": item.get("owner", "?"),
         "source_file": item.get("source_file", "?"),
         "function": item.get("function", "?"),
         "mainline_form": item.get("mainline_form", "?"),
         "live_under": item.get("live_under", "?"),
         "classification": item.get("classification", "?"),
         "execution": item.get("execution", "?")}
        for item in register.get("entries", [])
    ]
    arms = register.get("execution_arms") or {}
    executed = [d for d in divergences if d["execution"] in arms]
    gated = [d for d in divergences if d["verdict"] == "gate"]

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": rel(DECLARATION),
        "sources": sorted(sources),
        "findings": findings,
        "release": register.get("release"),
        "fork_reference": register.get("fork_reference"),
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
        print("The simulated vegetation's fire operators, checked against the "
              "source.\n")
        broken = [f for f in fixtures if not f["pass"]]
        print(f"  fixtures: {len(fixtures) - len(broken)} of {len(fixtures)} got their verdict")
        for case in broken:
            print(f"    BROKEN: {case['fixture']}: expected {case['expected']}, found {case['found']}")
        if findings:
            print(f"\n  {len(findings)} finding(s):")
            for finding in findings:
                print(f"    [{finding['kind']}] {finding['what']}: {finding['detail']}")
        else:
            print("  every declared divergence is recorded, absent from the code, "
                  "and claims its own execution standing")
        kept = [d for d in divergences if d["verdict"] == "keep"]
        print(f"\n  {len(divergences)} declared divergence(s) from the release,")
        print(f"  {len(kept)} kept and {len(gated)} gated, over "
              f"{len(sources)} operator source(s):")
        print(f"    {register.get('release')}")
        for item in divergences:
            print(f"    [{item['verdict']}] {item['id']}  {item['mainline_form']}  "
                  f"live under {item['live_under']}  ({item['source_file']}: "
                  f"{item['function']})  [{item['owner']}]")
            print(f"        classification: {item['classification']}   "
                  f"execution: {item['execution']}")
        print(f"\n  {len(executed)} of {len(divergences)} claim a matched execution "
              "arm, per entry and not per")
        print("  register. Each says for itself which arm covers it, or why none does.")
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
