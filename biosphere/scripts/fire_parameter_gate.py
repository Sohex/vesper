"""The fire parameter gate: every Earth constant the replacement fire model
needs, held to its declared disposition and to the source it is attributed to.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a vegetation model's fire operators, the file that declares
the constants they run on, and the source the model actually compiles.

`biosphere/config/fire_parameters.yaml` is the declaration and this module is
the enforcement. IT IS NOT `fire_gate.py`, which enforces a different
declaration about the same area: that gate reads `fire.yaml` and asks whether
this project's DIVERGENCES from mainline are still real. This one asks whether
every constant the port will inherit still carries a disposition, still sits
inside its own bracket, and is still on the line it is attributed to.

WHY A REGISTER NEEDS A GATE AT ALL. CLAUDE.md's second working agreement on
field updates: a derived quantity written down somewhere is an update that
cannot propagate, and the question to ask of any number is what re-runs when it
changes. A register of nineteen Earth constants with nothing reading it is the
frozen state in its purest form -- and `world-nfer` already found eleven
brackets in four components that nothing read, one of which had drifted outside
itself. So the declaration is written with a checker in the same commit.

It can fail:

  schema        an entry missing a key every entry needs, or naming a
                disposition, defect or seam that is not one of the declared set
  bracket       a `bracketed` entry with no bracket, a bracket that is not two
                ordered numbers, a central value outside its own bracket, or a
                bracket on an entry whose disposition says it is not swept --
                the last because an unswept bracket is a number wearing an
                interval, which reads as disposed of and is not
  restatement   an entry that says a value is a literal in the compiled source
                and names a line the source no longer carries. THIS IS THE ONE
                THAT FAILS ON SOMEBODY ELSE'S EDIT, and it is why the register
                names `source_file` and `literal` rather than describing them
  loop          an entry whose value reaches no run (`route: none`) and whose
                `consumed_by` row is closed. While the row is open the entry is
                waiting for a consumer that does not exist yet, which is honest;
                once the row closes, a constant nothing reads is frozen and the
                register has drifted from the model
  anchor        a constant declared to follow from `config/planet.yaml` that no
                longer equals what that file declares. `average_fire_duration`
                is the case: its whole finding is that the Earth source anchors
                a fire's duration to one rotation, so a port that wrote 24 has
                reproduced the number instead of the mechanism
  ordering      a threshold pair whose lower end is not below its upper end,
                for the pairs that are one observation reported as an interval
  partition     the fire nitrogen emission ratios not summing to one, or a
                declared share disagreeing with the literal in the source. They
                are a partition of a flux, so a sum that is not one means the
                model creates or destroys nitrogen at the fire
  derivation    a central value that stopped equalling the expression it is
                declared to follow from, or -- the open case -- a source that
                compiles a value the derivation contradicts. The second fires
                every run by design: blaze.cpp's rate-of-spread coefficient is
                Noble (1980)'s conversion with a lost decimal, and the finding
                is meant to report until the source is repaired
  enforcement   a fail-closed entry whose refusal the model no longer makes,
                or one that names no source making it. A declaration the model
                does not share is a preference
  absence       a constant declared DELETED that the source carries again, or
                one claiming both `literal` and `absent_from_source`. Comments
                are stripped first, so the comment left at the deletion site
                naming the constant does not fire it
  unresolved    an entry that names a source as not held while that source is
                in fact in `references/pdf/`. A bracket standing on "the paper
                was not read" is void the moment the paper arrives

Reduced fixtures run on every invocation, each built to be wrong in a named
way. A fixture that does not get the verdict it was built for is a defect in
this checker rather than in the declaration.

    python biosphere/scripts/fire_parameter_gate.py            # status, exit 0
    python biosphere/scripts/fire_parameter_gate.py --strict   # refuses on the
                                                               # named residual

`--strict` refuses on exactly one thing: an entry marked `fail_closed` whose
`central` is still null. There is one -- `phosphorus_volatilised_fraction` --
and the refusal is fire-7's: enabling either fire path under phosphorus
limitation would update carbon and nitrogen while leaving the corresponding
phosphorus in place, which breaks tissue stoichiometry silently rather than
loudly. A null there is a refusal and never a default of zero, because zero is
the claim that fire volatilises no phosphorus and nobody has made it.

IT READS AND RUNS NOTHING. Every input is a file: the declaration, three C++
translation units, `config/planet.yaml`, the tracked issue export and a
directory listing. It belongs in the per-commit tier for the same reason
`mineral_reactivity_gate.py` does.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import GENERATED, PROJECT_ROOT
from paths import rel  # noqa: E402

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = COMPONENT_ROOT / "config" / "fire_parameters.yaml"
REPORT = GENERATED / "fire_parameter_gate_report.json"
PLANET = PROJECT_ROOT / "config" / "planet.yaml"
ISSUE_EXPORT = PROJECT_ROOT / ".beads" / "issues.jsonl"
REFERENCE_PDFS = PROJECT_ROOT / "references" / "pdf"

DISPOSITIONS = ("sourced", "derived", "bracketed", "irreducible")
DEFECTS = ("tuned", "opaque", "implicit_earth", "none")
SEAMS = ("ignition", "occurrence_to_area", "effects", "cnp_destinations")
ROUTES = ("instruction_file", "compiled_literal", "none")

# Keys every entry carries, whatever its disposition. `central` is required and
# may be null: a null is a declared ABSENCE, which is a different statement from
# a missing key and is the shape `phosphorus_volatilised_fraction` needs.
REQUIRED = ("seam", "what", "unit", "central", "disposition", "defect",
            "route", "consumed_by", "source")

# The pairs that are one observation reported as an interval, so a sweep that
# crossed them would report a band running backwards.
ORDERED_PAIRS = (("fuel_threshold_lower", "fuel_threshold_upper"),
                 ("relative_humidity_lower", "relative_humidity_upper"))

# An entry whose central value is not free: it follows from `config/planet.yaml`
# by the named key, and a port that hardcoded the Earth number fails here.
PLANET_ANCHORED = {"average_fire_duration": ("rotation_hours", 1e-9)}


def _strip_c_comments(text: str) -> str:
    """The source with C and C++ comments removed.

    A literal that survives only inside a comment is a record of the constant
    and not the constant, which is the distinction `fire_gate.py` already turns
    on. Strings are not preserved: no declared literal here is inside one.
    """
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", " ", text)


def _numbers(value) -> list[float]:
    """The numbers a declaration states, refusing a bool and a numeric string.

    A YAML 1.1 float wants a sign in its exponent, so `5.0e4` resolves to a
    STRING that reads as a number everywhere until something compares against
    it. `smoke_test.py:check_declared_numerics_resolve_to_numbers` refuses that
    form in the file, so there is nothing here to coerce and coercing would
    absorb the defect.
    """
    if isinstance(value, bool) or value is None:
        return []
    if isinstance(value, (int, float)):
        return [float(value)]
    if isinstance(value, (list, tuple)):
        out = []
        for item in value:
            out.extend(_numbers(item))
        return out
    if isinstance(value, dict):
        out = []
        for item in value.values():
            out.extend(_numbers(item))
        return out
    return []


def _finding(kind: str, what: str, detail: str) -> dict:
    return {"kind": kind, "what": what, "detail": detail}


def check_schema(params: dict) -> list[dict]:
    bad = []
    for name, entry in params.items():
        if not isinstance(entry, dict):
            bad.append(_finding("schema", name, "is not a mapping"))
            continue
        for key in REQUIRED:
            if key not in entry:
                bad.append(_finding("schema", name, f"has no `{key}`"))
        for key, allowed in (("disposition", DISPOSITIONS), ("defect", DEFECTS),
                             ("seam", SEAMS), ("route", ROUTES)):
            got = entry.get(key)
            if key in entry and got not in allowed:
                bad.append(_finding(
                    "schema", name,
                    f"`{key}` is {got!r}, which is not one of {allowed}"))
    return bad


def check_brackets(params: dict) -> list[dict]:
    """A bracketed value inside its own bracket, and no bracket anywhere else.

    Both halves. A central outside its bracket is the failure `world-nfer`
    found in the tree at large; a bracket on an entry that is not swept is the
    quieter one, because an interval printed beside a number reads as a
    disposition that has been taken and there is nothing behind it.
    """
    bad = []
    for name, entry in params.items():
        if not isinstance(entry, dict):
            continue
        bracket = entry.get("bracket")
        swept = entry.get("disposition") == "bracketed"
        if swept and bracket is None:
            bad.append(_finding("bracket", name,
                                "is `bracketed` and declares no bracket, so "
                                "nothing is swept"))
            continue
        if not swept:
            if bracket is not None:
                bad.append(_finding(
                    "bracket", name,
                    f"declares a bracket but its disposition is "
                    f"{entry.get('disposition')!r}, so nothing sweeps it; a "
                    f"bracket nothing sweeps reads as a disposition that was "
                    f"taken"))
            continue
        ends = _numbers(bracket)
        if len(ends) != 2:
            bad.append(_finding("bracket", name,
                                f"bracket {bracket!r} is not two numbers"))
            continue
        low, high = ends
        if not low < high:
            bad.append(_finding("bracket", name,
                                f"bracket [{low}, {high}] is not ordered"))
            continue
        for value in _numbers(entry.get("central")):
            if not low <= value <= high:
                bad.append(_finding(
                    "bracket", name,
                    f"central {value} is outside its own bracket "
                    f"[{low}, {high}]"))
    return bad


def check_restatements(params: dict, root: Path) -> list[dict]:
    """Every value attributed to a line of compiled source is still on it."""
    bad = []
    for name, entry in params.items():
        if not isinstance(entry, dict) or entry.get("route") != "compiled_literal":
            continue
        where = entry.get("source_file")
        literal = entry.get("literal")
        if not where or literal is None:
            bad.append(_finding(
                "restatement", name,
                "is a compiled literal and names no `source_file`/`literal`"))
            continue
        path = root / where
        if not path.is_file():
            bad.append(_finding("restatement", name,
                                f"{where} is not a file"))
            continue
        code = _strip_c_comments(path.read_text(encoding="utf-8",
                                                errors="replace"))
        if literal not in code:
            bad.append(_finding(
                "restatement", name,
                f"{where} no longer carries {literal!r} in code, so the value "
                f"is attributed to a line that is not there"))
    return bad


def _issue_status(root: Path) -> dict[str, str]:
    """Status by issue id, from the TRACKED export rather than from `bd`.

    A static read, which is what keeps this gate in the per-commit tier. The
    export is the copy every other checkout sees, so a check against it is a
    check against what a fresh reader would find.
    """
    path = root / ".beads" / "issues.jsonl"
    if not path.is_file():
        return {}
    out = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            row = json.loads(line)
        except ValueError:
            continue
        if isinstance(row, dict) and row.get("id"):
            out[row["id"]] = str(row.get("status", ""))
    return out


def check_loops(params: dict, statuses: dict[str, str]) -> list[dict]:
    """A constant that reaches no run is waiting, and only while its row is open.

    `route: none` is the honest state for a constant whose consumer has not
    been written: the replacement fire model does not exist yet, so nine of
    these entries reach nothing. What makes it a defect rather than a state is
    the row CLOSING -- at that point a consumer exists, the constant still
    reaches no run, and the register has drifted from the model with nothing
    objecting.
    """
    bad = []
    if not statuses:
        bad.append(_finding("loop", "issue export",
                            ".beads/issues.jsonl is absent or empty, so no "
                            "entry's consumer could be checked; this gate "
                            "cannot tell a waiting constant from a frozen one"))
        return bad
    for name, entry in params.items():
        if not isinstance(entry, dict) or entry.get("route") != "none":
            continue
        if entry.get("absent_from_source"):
            # Deleted, not waiting. `check_absences` holds it to still being
            # gone, which is a stronger statement than naming an open row.
            continue
        row = entry.get("consumed_by")
        status = statuses.get(row)
        if status is None:
            bad.append(_finding(
                "loop", name,
                f"`consumed_by` names {row!r}, which the tracked export does "
                f"not carry"))
        elif status.lower() in ("closed", "done", "resolved"):
            bad.append(_finding(
                "loop", name,
                f"reaches no run (`route: none`) and its consumer {row} is "
                f"{status}, so a constant nothing reads now sits behind a "
                f"model that was written"))
    return bad


def check_anchors(params: dict, planet: dict) -> list[dict]:
    """A constant declared to follow from the planet still equals it."""
    bad = []
    declared = (planet or {}).get("planet") or {}
    for name, (key, tol) in PLANET_ANCHORED.items():
        entry = params.get(name)
        if not isinstance(entry, dict):
            continue  # presence is `check_completeness`'s question, not this one
        want = declared.get(key)
        got = entry.get("central")
        if not isinstance(want, (int, float)) or isinstance(want, bool):
            bad.append(_finding("anchor", name,
                                f"config/planet.yaml declares no numeric "
                                f"`planet.{key}`"))
            continue
        if not isinstance(got, (int, float)) or isinstance(got, bool):
            bad.append(_finding("anchor", name,
                                f"central {got!r} is not a number"))
            continue
        if abs(float(got) - float(want)) > tol:
            bad.append(_finding(
                "anchor", name,
                f"central {got} does not equal config/planet.yaml's "
                f"planet.{key} = {want}. The Earth source anchors this to one "
                f"rotation, so a value that stopped tracking the declaration "
                f"has reproduced a number instead of a mechanism"))
    return bad


def check_ordering(params: dict) -> list[dict]:
    bad = []
    for lower, upper in ORDERED_PAIRS:
        a, b = params.get(lower), params.get(upper)
        if not isinstance(a, dict) or not isinstance(b, dict):
            continue  # presence is `check_completeness`'s question
        av, bv = _numbers(a.get("central")), _numbers(b.get("central"))
        if len(av) != 1 or len(bv) != 1:
            bad.append(_finding("ordering", f"{lower}/{upper}",
                                "the pair does not state one number each"))
            continue
        if not av[0] < bv[0]:
            bad.append(_finding(
                "ordering", f"{lower}/{upper}",
                f"{av[0]} is not below {bv[0]}; the two are one observation "
                f"reported as an interval and a sweep that crossed them would "
                f"report a band running backwards"))
    return bad


def check_partition(params: dict, root: Path) -> list[dict]:
    """The fire nitrogen shares sum to one and match the compiled source.

    A partition that does not sum to one is a model that creates or destroys
    nitrogen at the fire, which the acceptance run's closure test would then
    attribute somewhere else. The second half is the restatement check applied
    per share rather than to the block, because the entry's `literal` names one
    of the four and a drift in the other three would pass it.
    """
    bad = []
    name = "fire_nitrogen_emission_ratios"
    entry = params.get(name)
    if not isinstance(entry, dict):
        return []  # presence is `check_completeness`'s question
    shares = entry.get("central")
    if not isinstance(shares, dict):
        return [_finding("partition", name, "central is not a mapping")]
    total = sum(_numbers(shares))
    if abs(total - 1.0) > 1e-9:
        bad.append(_finding(
            "partition", name,
            f"the four shares sum to {total!r}, not 1. They are a partition "
            f"of a flux, so a sum that is not one is nitrogen created or "
            f"destroyed at the fire"))
    where = entry.get("source_file")
    path = root / where if where else None
    if path and path.is_file():
        code = _strip_c_comments(path.read_text(encoding="utf-8",
                                                errors="replace"))
        for species, value in shares.items():
            pattern = rf"{re.escape(species)}_FIRERATIO\s*=\s*([0-9.eE+-]+)"
            found = re.search(pattern, code)
            if not found:
                bad.append(_finding(
                    "partition", name,
                    f"{where} declares no {species}_FIRERATIO in code"))
                continue
            if abs(float(found.group(1)) - float(value)) > 1e-12:
                bad.append(_finding(
                    "partition", name,
                    f"{species}: the register says {value} and {where} "
                    f"compiles {found.group(1)}"))
    return bad


def check_unresolved(params: dict, pdf_dir: Path) -> list[dict]:
    """A bracket standing on an unread paper, when the paper is now held.

    `unresolved` is a promise that a disposition is provisional because a
    source could not be reached. The promise expires when the source arrives,
    and it expires SILENTLY unless something looks -- so this is the check that
    turns "acquire references proactively" into a thing the tree enforces
    rather than a thing a session remembers.

    Matched PER FILE rather than against one joined string, because
    `references/pdf/` holds nothing but papers and its filenames carry surname
    and year in the same name. Requiring both in ONE filename is what stops a
    Thonicke (2001) already on disk from satisfying a note about Noble (1980)
    just because some other file carries the year.
    """
    bad = []
    if not pdf_dir.is_dir():
        return [_finding("unresolved", "references/pdf",
                         "is not a directory, so no unresolved source could "
                         "be checked")]
    names = [p.name.lower() for p in pdf_dir.iterdir() if p.is_file()]
    # Words that appear in every citation and name nobody.
    stop = {"not", "and", "the", "for", "held", "references", "pdf", "doi",
            "which", "reading", "this", "entry", "from", "sourced", "bracket",
            "stands", "until", "that", "what", "would", "say", "whether",
            "convert", "moves", "tuned", "them", "with", "alone", "comment",
            "identified", "single", "number", "has", "stated", "spread", "are"}
    for name, entry in params.items():
        if not isinstance(entry, dict) or not entry.get("unresolved"):
            continue
        text = str(entry["unresolved"]).lower()
        years = set(re.findall(r"\b(1[89]\d\d|20\d\d)\b", text))
        words = {w for w in re.findall(r"[a-z][a-z-]{2,}", text)
                 if w not in stop}
        if not years or not words:
            continue
        for filename in names:
            if any(y in filename for y in years) and any(
                    w in filename for w in words):
                bad.append(_finding(
                    "unresolved", name,
                    f"says a source is not held, and references/pdf/ now has "
                    f"{filename!r}, which names both. The bracket stands on "
                    f"the source being unreachable, so it is void until the "
                    f"entry is reread against it"))
                break
    return bad


def check_derivations(params: dict, root: Path) -> list[dict]:
    """A value declared to FOLLOW from an expression still follows from it, and
    what the source compiles is still what the entry says it compiles.

    Two halves, and they fail in opposite directions on purpose.

    `derived_from` is arithmetic written out -- a unit conversion, not a fit --
    so it is EVALUATED here rather than trusted. An entry whose `central` stops
    equalling its own expression has had one of the two edited alone.

    `implemented_value` is what the compiled source actually carries where that
    DISAGREES with the derivation. It is the finding, held open deliberately:
    while the two differ, this reports the gap on every run, and the way to
    silence it is to repair the source and delete the field. A gate that let
    the disagreement sit silently would be the frozen state with a defect in
    it.
    """
    bad = []
    for name, entry in params.items():
        if not isinstance(entry, dict):
            continue
        expression = entry.get("derived_from")
        if expression:
            try:
                # Arithmetic only: the expression is a unit conversion written
                # out, and anything that is not a number and an operator is a
                # declaration this gate should refuse rather than execute.
                if not re.fullmatch(r"[0-9eE.+\-*/() ]+", str(expression)):
                    raise ValueError("not arithmetic")
                got = float(eval(str(expression), {"__builtins__": {}}, {}))
            except Exception as exc:                              # noqa: BLE001
                bad.append(_finding(
                    "derivation", name,
                    f"`derived_from` {expression!r} did not evaluate as "
                    f"arithmetic: {exc!r}"))
                continue
            central = entry.get("central")
            if not isinstance(central, (int, float)) or isinstance(central, bool):
                bad.append(_finding("derivation", name,
                                    "declares `derived_from` and its central "
                                    "is not a number"))
            elif abs(float(central) - got) > abs(got) * 1e-6:
                bad.append(_finding(
                    "derivation", name,
                    f"central {central} is not what `derived_from` evaluates "
                    f"to ({got!r}); one of the two was edited alone"))
        implemented = entry.get("implemented_value")
        if implemented is None:
            continue
        where, literal = entry.get("source_file"), entry.get("literal")
        if not where or literal is None:
            bad.append(_finding("derivation", name,
                                "declares `implemented_value` and names no "
                                "source line to read it from"))
            continue
        code = _strip_c_comments((root / where).read_text(
            encoding="utf-8", errors="replace")) if (root / where).is_file() else ""
        found = re.search(r"=\s*([0-9.]+[eE][+-]?[0-9]+|[0-9]*\.?[0-9]+)",
                          literal)
        if not found or found.group(1) not in code:
            bad.append(_finding(
                "derivation", name,
                f"{where} no longer compiles {literal!r}, so the recorded "
                f"disagreement cannot be confirmed against the source"))
            continue
        if not entry.get("row"):
            bad.append(_finding(
                "derivation", name,
                f"{where} compiles {implemented!r} against a derivation of "
                f"{entry.get('central')!r} and the entry names no `row`. A "
                f"disagreement between a source and its own derivation is a "
                f"defect, and an unowned one is a defect nobody is going to "
                f"fix"))
    return bad


def declared_disagreements(params: dict) -> list[dict]:
    """The source-versus-derivation gaps that are DECLARED and OWNED.

    These are reported on every run and do not fail, which is the same split
    `fire_gate.py` makes between a divergence whose verdict is `keep` and one
    whose verdict is `gate`. A gate that failed on a defect already written
    down and already tracked would fail forever, and a check that always fails
    is one people learn to run with `|| true`.

    What DOES fail is the pair drifting: `check_derivations` refuses an entry
    whose recorded `implemented_value` the source no longer compiles, and one
    that names no owning row. So the disagreement is allowed to sit only for
    exactly as long as it is described correctly and owned.
    """
    out = []
    for name, entry in params.items():
        if not isinstance(entry, dict):
            continue
        implemented, central = entry.get("implemented_value"), entry.get("central")
        if implemented is None or not isinstance(central, (int, float)):
            continue
        out.append({
            "parameter": name,
            "row": entry.get("row"),
            "source_file": entry.get("source_file"),
            "compiles": implemented,
            "derivation": central,
            "factor": float(implemented) / float(central) if central else None,
        })
    return out


def check_absences(params: dict, root: Path) -> list[dict]:
    """A constant declared DELETED is still absent from the code.

    The shape `fire_gate.py` already uses for a deleted line. An entry whose
    constant was removed from the source keeps its row -- a register that
    forgets a deleted constant cannot notice it coming back -- and
    `absent_from_source` is what turns that row into a check instead of a
    memorial. Comments are stripped first, because the deletion leaves a comment
    NAMING the constant at the site, and a search that matched it would fire on
    the record of the repair rather than on the repair being undone.
    """
    bad = []
    for name, entry in params.items():
        if not isinstance(entry, dict):
            continue
        token = entry.get("absent_from_source")
        if not token:
            continue
        where = entry.get("source_file")
        if not where:
            bad.append(_finding("absence", name,
                                "declares `absent_from_source` and names no "
                                "`source_file` to check it against"))
            continue
        path = root / where
        if not path.is_file():
            bad.append(_finding("absence", name, f"{where} is not a file"))
            continue
        code = _strip_c_comments(path.read_text(encoding="utf-8",
                                                errors="replace"))
        if re.search(rf"\b{re.escape(str(token))}\b", code):
            bad.append(_finding(
                "absence", name,
                f"is declared deleted and {where} carries {token!r} in code "
                f"again. The repair has been undone and the register still "
                f"says it stands."))
        if entry.get("literal"):
            bad.append(_finding(
                "absence", name,
                "declares both `literal` and `absent_from_source`, which are "
                "opposite claims about the same source"))
    return bad


def check_enforcement(params: dict, root: Path) -> list[dict]:
    """A fail-closed entry's refusal is still in the model that has to make it.

    A declaration the model does not share is a preference. An entry marked
    `fail_closed` may name `enforced_in` and `enforced_by`, and this holds that
    source to still carrying that condition in CODE -- comments stripped, so a
    refusal deleted but described in a comment fails rather than passes.

    It is deliberately a string match on the condition rather than a parse. The
    failure it is built for is the refusal being DELETED or its condition
    narrowed, not its whitespace changing, and a check that tried to understand
    C++ would be a second model of the language that could disagree with the
    compiler.
    """
    bad = []
    for name, entry in params.items():
        if not isinstance(entry, dict) or not entry.get("fail_closed"):
            continue
        where, condition = entry.get("enforced_in"), entry.get("enforced_by")
        if not where or not condition:
            bad.append(_finding(
                "enforcement", name,
                "is fail-closed and names no `enforced_in`/`enforced_by`, so "
                "the refusal lives only in this file and the model does not "
                "make it"))
            continue
        path = root / where
        if not path.is_file():
            bad.append(_finding("enforcement", name, f"{where} is not a file"))
            continue
        code = _strip_c_comments(path.read_text(encoding="utf-8",
                                                errors="replace"))
        if str(condition) not in code:
            bad.append(_finding(
                "enforcement", name,
                f"{where} no longer carries the condition {condition!r} in "
                f"code, so the fail-closed refusal this entry claims is not "
                f"the one the model makes"))
    return bad


# Entries the anchor, ordering and partition checks are ABOUT. Those three ask
# whether a named constant is still right, which is a different question from
# whether it is still there, and a check that answered both would report a
# reduced fixture's missing entry as a defect in the declaration. So presence is
# asked once, here, and only of a declaration claiming to be the whole register.
REQUIRED_IDS = tuple(sorted(
    set(PLANET_ANCHORED)
    | {name for pair in ORDERED_PAIRS for name in pair}
    | {"fire_nitrogen_emission_ratios", "phosphorus_volatilised_fraction"}))


def check_completeness(params: dict) -> list[dict]:
    """The register still carries every entry another check is about.

    Deleting an entry is how a constant escapes its disposition without any
    check firing: the bracket goes with it, the restatement goes with it, and
    the value carries on reaching the model from the source it was registered
    against. Naming the ids here is what makes a deletion louder than an edit.
    """
    return [_finding("missing", name,
                     "is registered as a constant this gate checks and is no "
                     "longer in the declaration; deleting an entry is how a "
                     "value escapes its disposition silently")
            for name in REQUIRED_IDS if not isinstance(params.get(name), dict)]


def check(declaration: dict, root: Path, planet: dict,
          statuses: dict[str, str], pdf_dir: Path,
          complete: bool = False) -> list[dict]:
    """Every check over one declaration.

    `complete` says whether this declaration claims to be the WHOLE register.
    A reduced fixture does not, so it is not asked whether it carries the
    nineteen entries; the real one is.
    """
    params = declaration.get("parameters") or {}
    findings = []
    findings.extend(check_schema(params))
    findings.extend(check_brackets(params))
    findings.extend(check_restatements(params, root))
    findings.extend(check_loops(params, statuses))
    findings.extend(check_anchors(params, planet))
    findings.extend(check_ordering(params))
    findings.extend(check_partition(params, root))
    findings.extend(check_derivations(params, root))
    findings.extend(check_absences(params, root))
    findings.extend(check_enforcement(params, root))
    findings.extend(check_unresolved(params, pdf_dir))
    if complete:
        findings.extend(check_completeness(params))
    return findings


def _fixtures(root: Path) -> list[dict]:
    """Reduced declarations, all but one built to be wrong in a named way.

    A checker that reports nothing is indistinguishable from a checker that
    finds nothing, and this one runs mostly over a declaration its own author
    wrote, so the negative result is the expected one and carries no
    information on its own. Each case below is a positive control for one
    `kind`: if it stops firing, this module has stopped looking rather than the
    declaration having become clean. CLAUDE.md's "a probe carries a control
    that must fire".
    """
    ok = {
        "seam": "ignition", "what": "x", "unit": "dimensionless",
        "central": 0.2, "bracket": [0.1, 0.3], "disposition": "bracketed",
        "defect": "implicit_earth", "route": "none",
        "consumed_by": "fire-3", "source": "x",
    }
    statuses = {"fire-3": "open", "fire-9": "closed"}
    planet = {"planet": {"rotation_hours": 30.0}}

    def one(entry, **over):
        merged = dict(ok)
        merged.update(entry)
        params = {"probe": merged}
        params.update(over.pop("extra", {}))
        return check({"parameters": params}, root, planet, statuses, pdf_dir=REFERENCE_PDFS)

    cases = []

    def case(label, findings, expected):
        found = sorted({f["kind"] for f in findings})
        cases.append({"fixture": label, "expected": expected, "found": found,
                      "pass": expected in found})

    case("a clean entry reports nothing",
         one({}), "nothing")
    cases[-1]["pass"] = cases[-1]["found"] == []
    cases[-1]["expected"] = "no findings"
    cases[-1]["found"] = cases[-1]["found"] or ["nothing"]

    case("a disposition outside the four",
         one({"disposition": "calibrated"}), "schema")
    case("a central below its own bracket",
         one({"central": 0.05}), "bracket")
    case("a bracket whose ends are swapped",
         one({"bracket": [0.3, 0.1]}), "bracket")
    case("a bracket on an entry nothing sweeps",
         one({"disposition": "sourced"}), "bracket")
    case("a literal the source does not carry",
         one({"route": "compiled_literal",
              "source_file": "vendor/lpj-guess/modules/blaze.cpp",
              "literal": "K_LITTER_VESPERIAN = 0.5"}), "restatement")
    case("a value reaching no run behind a closed consumer",
         one({"consumed_by": "fire-9"}), "loop")
    case("a planet-anchored value that stopped tracking the planet",
         check({"parameters": {"average_fire_duration": dict(
             ok, central=24.0, bracket=[15.0, 60.0])}},
             root, planet, statuses, REFERENCE_PDFS), "anchor")
    case("a threshold pair running backwards",
         check({"parameters": {
             "fuel_threshold_lower": dict(ok, central=900.0, bracket=[100.0, 1000.0]),
             "fuel_threshold_upper": dict(ok, central=200.0, bracket=[100.0, 1000.0])}},
             root, planet, statuses, REFERENCE_PDFS), "ordering")
    case("nitrogen shares that do not sum to one",
         check({"parameters": {"fire_nitrogen_emission_ratios": dict(
             ok, disposition="sourced", bracket=None,
             central={"NH3": 0.005, "NOx": 0.237, "N2O": 0.036, "N2": 0.9},
             route="compiled_literal",
             source_file="vendor/lpj-guess/framework/guess.cpp",
             literal="NH3_FIRERATIO = 0.005")}},
             root, planet, statuses, REFERENCE_PDFS), "partition")
    case("an unresolved source that is in fact held",
         one({"unresolved": "Thonicke et al. (2010). Not in references/pdf/."}),
         "unresolved")
    case("a central that stopped following its own derivation",
         one({"central": 9.9, "derived_from": "0.0012 / 100.0 / 3.6"}),
         "derivation")
    case("a source/derivation gap that names no owning row",
         one({"central": 3.333333e-06, "derived_from": "0.0012 / 100.0 / 3.6",
              "implemented_value": 3.3333e-05, "route": "compiled_literal",
              "source_file": "vendor/lpj-guess/modules/blaze.cpp",
              "literal": "A = 3.3333e-05"}), "derivation")
    case("a constant declared deleted that the source carries again",
         one({"absent_from_source": "MIN_FUEL",
              "source_file": "vendor/lpj-guess/modules/blaze.cpp",
              "route": "none"}), "absence")
    case("an entry claiming both a literal and an absence",
         one({"absent_from_source": "NOT_IN_THE_SOURCE_AT_ALL",
              "literal": "HEAT_YIELD = 20.",
              "source_file": "vendor/lpj-guess/modules/blaze.cpp",
              "route": "compiled_literal"}), "absence")
    case("a fail-closed entry whose model refusal is gone",
         one({"fail_closed": "yes",
              "enforced_in": "vendor/lpj-guess/framework/parameters.cpp",
              "enforced_by": "ifplim && firemodel == NOT_A_REAL_CONDITION"}),
         "enforcement")
    case("a fail-closed entry naming no enforcing source",
         one({"fail_closed": "yes"}), "enforcement")
    case("a register missing an entry another check is about",
         check({"parameters": {"probe": dict(ok)}}, root, planet, statuses,
               REFERENCE_PDFS, complete=True), "missing")
    return cases


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict", action="store_true",
                        help="refuse while a fail-closed entry is undeclared")
    parser.add_argument("--json", action="store_true", help="report only")
    args = parser.parse_args()

    declaration = yaml.safe_load(DECLARATION.read_text(encoding="utf-8"))
    planet = yaml.safe_load(PLANET.read_text(encoding="utf-8"))
    statuses = _issue_status(PROJECT_ROOT)
    params = declaration.get("parameters") or {}

    fixtures = _fixtures(PROJECT_ROOT)
    findings = check(declaration, PROJECT_ROOT, planet, statuses,
                     REFERENCE_PDFS, complete=True)

    by_seam: dict[str, list[str]] = {}
    by_defect: dict[str, list[str]] = {}
    for name, entry in params.items():
        by_seam.setdefault(entry.get("seam", "?"), []).append(name)
        by_defect.setdefault(entry.get("defect", "?"), []).append(name)

    blocked = [name for name, entry in params.items()
               if entry.get("fail_closed") and entry.get("central") is None]
    disagreements = declared_disagreements(params)

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": rel(DECLARATION),
        "parameters": len(params),
        "by_seam": {k: sorted(v) for k, v in sorted(by_seam.items())},
        "by_defect": {k: sorted(v) for k, v in sorted(by_defect.items())},
        # A deleted constant also has `route: none` and is NOT waiting for a
        # consumer, so the two are reported apart. Conflating them would say a
        # constant that no longer exists is behind an open row.
        "reaching_no_run": sorted(name for name, e in params.items()
                                  if e.get("route") == "none"
                                  and not e.get("absent_from_source")),
        "deleted_from_source": sorted(name for name, e in params.items()
                                      if e.get("absent_from_source")),
        "fail_closed": sorted(blocked),
        "declared_disagreements": disagreements,
        "findings": findings,
        "fixtures": fixtures,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print("The Earth constants the replacement fire model inherits, held "
              "to their dispositions.\n")
        broken = [f for f in fixtures if not f["pass"]]
        print(f"  fixtures: {len(fixtures) - len(broken)} of {len(fixtures)} "
              f"got their verdict")
        for c in broken:
            print(f"    BROKEN: {c['fixture']}: expected {c['expected']}, "
                  f"found {c['found']}")
        print(f"\n  {len(params)} constants registered")
        for seam in SEAMS:
            names = report["by_seam"].get(seam) or []
            print(f"    {seam:<20} {len(names)}")
        print("\n  by defect class:")
        for defect in DEFECTS:
            names = report["by_defect"].get(defect) or []
            if names:
                print(f"    {defect:<16} {len(names)}  {', '.join(names)}")
        if report["deleted_from_source"]:
            print(f"\n  {len(report['deleted_from_source'])} deleted from the "
                  f"source and held to staying deleted: "
                  f"{', '.join(report['deleted_from_source'])}")
        if report["reaching_no_run"]:
            print(f"\n  {len(report['reaching_no_run'])} reach no run yet, "
                  f"each behind an open row: "
                  f"{', '.join(report['reaching_no_run'])}")
        if findings:
            print(f"\n  {len(findings)} finding(s):")
            for f in findings:
                print(f"    [{f['kind']}] {f['what']}: {f['detail']}")
        else:
            print("\n  every constant carries a disposition, sits inside its "
                  "own bracket, and is on the line it names")
        for d in disagreements:
            print(f"\n  DECLARED DISAGREEMENT ({d['row']}): {d['parameter']}\n"
                  f"    {d['source_file']} compiles {d['compiles']!r}; the "
                  f"derivation gives {d['derivation']!r}, a factor of "
                  f"{d['factor']:.4f}.\n"
                  f"    Owned and reported, not failed. It stops reporting "
                  f"when the source is repaired and `implemented_value` is "
                  f"deleted.")
        if blocked:
            print(f"\n  FAIL CLOSED: {', '.join(blocked)} is undeclared. No "
                  f"fire path may run under phosphorus limitation while it is.")

    if broken:
        return 1
    if findings:
        return 1
    if args.strict and blocked:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
