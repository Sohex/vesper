#!/usr/bin/env python3
"""Is every member of the Soil class accounted for as restart state?

    python biosphere/scripts/soil_restart_state_gate.py
    python biosphere/scripts/soil_restart_state_gate.py --list-unclassified

Worldbuilding frame: this checks the Vesper project's vegetation model source
against a classification of it. Nothing here is a claim about the simulated
planet.

WHAT IT TESTS. `Soil::serialize` in `vendor/lpj-guess/modules/soil.cpp` is the
whole of what a restarted soil column inherits, and the Soil class in
`vendor/lpj-guess/framework/guess.h` is what such a column has. A member of the
second that is absent from the first is either state the restart silently drops
or something the model rebuilds before it needs it, and the two look identical
in a diff. `biosphere/config/soil_restart_state.yaml` says which each one is,
with the file:line that settles it, and this parses both sides and refuses when
they disagree.

It fails on four things, each of which has a right answer rather than merely a
different one:

  unclassified   a member the class declares that neither Soil::serialize
                 streams nor the classification names. A member added to the
                 class reaches this state and stays there until someone decides
                 what it is, which is the point.
  contradicted   a member the classification calls recomputed or diagnostic
                 that Soil::serialize now streams. Either the serializer gained
                 bytes that change no result, or the classification is stale.
  lost           a non-empty `lost` block. That block exists so a finding has
                 somewhere to sit; while anything is in it, a resumed run is
                 not the run it continues.
  unknown        a name in the classification that the Soil class no longer
                 declares, or a class key that is not one of the three.

WHAT IT IS NOT. It reads source, so it can say a member is written before it is
read; it cannot say the resumed run reproduces the uninterrupted one. That is
`biosphere/scripts/verify_lpj_restart_continuity.py`, which needs a compiled
model and a built forcing. This is the half that runs on a tree with neither.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _paths import COMPONENT_ROOT, GENERATED, GUESS_SOURCE, PROJECT_ROOT  # noqa: E402

DECLARATION = COMPONENT_ROOT / "config" / "soil_restart_state.yaml"
CLASSES = ("lost", "recomputed", "diagnostic")

# Declared in the class body but not state: the two references the constructor
# wires to the owning objects, and the file-scope constants. Naming them here
# rather than filtering on a pattern, because a pattern would also swallow a
# member someone adds next to them.
NOT_STATE = ("patch", "soiltype", "ifallowphasechanges", "snowdensityconstant")


def class_body(header: str) -> str:
    """The Soil class declaration, down to its member functions."""
    start = header.find("class Soil : public Serializable {")
    if start < 0:
        raise SystemExit(
            f"{GUESS_SOURCE / 'framework' / 'guess.h'} declares no "
            "`class Soil : public Serializable`, so this gate has nothing to "
            "read. Has the class been renamed?")
    depth = 0
    for index in range(header.find("{", start), len(header)):
        if header[index] == "{":
            depth += 1
        elif header[index] == "}":
            depth -= 1
            if depth == 0:
                body = header[start:index]
                break
    else:
        raise SystemExit("the Soil class declaration does not close")
    functions = body.find("// MEMBER FUNCTIONS")
    return body[:functions] if functions > 0 else body


def declared_members(header: str) -> dict[str, str]:
    """Every data member the Soil class declares, to its declared type.

    A declaration is a line ending in `;` whose type is one token and whose
    remainder is names and array bounds. Anything with a parenthesis is a
    function; anything the parse cannot resolve is left out, and the sweep is
    therefore a floor on the class rather than a claim to have read all of it.
    """
    members: dict[str, str] = {}
    for line in class_body(header).splitlines():
        text = line.strip()
        if not text or text.startswith(("//", "/*", "*")) or "(" in text:
            continue
        match = re.match(
            r"^(?:static\s+)?(?:const\s+)?"
            r"([A-Za-z_]\w*(?:\s*<[^>]*>)?)\s+(.+?);\s*(?://.*)?$", text)
        if not match:
            continue
        kind, rest = match.group(1), match.group(2)
        if kind in ("public", "private", "protected", "return", "class",
                    "struct", "typedef"):
            continue
        for part in rest.split(","):
            name = re.match(r"^([A-Za-z_]\w*)((?:\s*\[[^\]]*\])*)$",
                            part.strip())
            if name:
                members[name.group(1)] = kind + name.group(2).strip()
    for name in NOT_STATE:
        members.pop(name, None)
    return members


def serialized_members(source: str) -> set[str]:
    """Every member name `Soil::serialize` streams into the archive.

    Parsed rather than grepped over the whole file, because what matters is what
    the RESTART carries: a member named anywhere else in soil.cpp is not restart
    state. The same parse as `wetland_gate.py:_serialized_members`.
    """
    start = source.find("void Soil::serialize(")
    if start < 0:
        raise SystemExit(
            f"{GUESS_SOURCE / 'modules' / 'soil.cpp'} defines no "
            "Soil::serialize, so this gate has nothing to read")
    depth = 0
    for index in range(source.find("{", start), len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                body = source[start:index]
                break
    else:
        raise SystemExit("Soil::serialize does not close")
    return set(re.findall(r"&\s*([A-Za-z_]\w*)", body)) - {"arch"}


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT_ROOT))
    except ValueError:
        return str(path)


def check(declaration: dict, members: dict[str, str],
          carried: set[str]) -> tuple[list[str], dict]:
    refusals: list[str] = []

    unknown_keys = sorted(set(declaration) - set(CLASSES))
    if unknown_keys:
        refusals.append(
            f"{rel(DECLARATION)} carries {', '.join(unknown_keys)}, which is "
            f"not one of {', '.join(CLASSES)}. A member has to be in a class "
            "this gate knows how to check.")

    classified: dict[str, str] = {}
    for name in CLASSES:
        for member in (declaration.get(name) or {}):
            classified[member] = name

    absent = sorted(set(members) - carried - set(classified))
    if absent:
        refusals.append(
            f"Soil::serialize does not carry {', '.join(absent)}, and "
            f"{rel(DECLARATION)} does not say why. Each is either restart state "
            "the resumed column loses, or something the model rebuilds before "
            "it reads it; decide which, with the file:line that settles it.")

    contradicted = sorted(m for m, k in classified.items()
                          if k != "lost" and m in carried)
    if contradicted:
        refusals.append(
            f"{rel(DECLARATION)} calls {', '.join(contradicted)} rebuilt or "
            "diagnostic while Soil::serialize now streams it. Either the "
            "serializer gained bytes that change no result, or the "
            "classification is stale.")

    vanished = sorted(m for m in classified if m not in members)
    if vanished:
        refusals.append(
            f"{rel(DECLARATION)} classifies {', '.join(vanished)}, which the "
            "Soil class no longer declares.")

    lost = sorted(declaration.get("lost") or {})
    if lost:
        refusals.append(
            f"{', '.join(lost)} is read before it is written on the first "
            "resumed day, so a run continued from a state file is not the run "
            "it continues.")

    unevidenced = sorted(
        member for name in CLASSES
        for member, why in (declaration.get(name) or {}).items()
        if not str(why or "").strip())
    if unevidenced:
        refusals.append(
            f"{', '.join(unevidenced)} is classified with no evidence. A "
            "classification without a file:line is an assertion, and the next "
            "sweep has to re-derive it.")

    report = {
        "declared_members": len(members),
        "serialized": sorted(carried & set(members)),
        "classified": {name: sorted(declaration.get(name) or {})
                       for name in CLASSES},
        "unclassified": absent,
        "contradicted": contradicted,
        "no_longer_declared": vanished,
    }
    return refusals, report


def fixtures(declaration: dict, members: dict[str, str],
             carried: set[str]) -> list[str]:
    """Reduced cases whose verdict is known before the gate is run.

    Six are built to be wrong in a named way and one is the live tree, which
    has to be GRANTED: a gate nothing can satisfy is a wall, and a gate nothing
    can fail is decoration. Each case is the smallest edit to the live inputs
    that should flip exactly one refusal, so a fixture that stops firing says
    the check it names stopped working rather than that the tree changed.
    """
    failures: list[str] = []

    def expect(label: str, code: str, want: bool,
               decl: dict, mem: dict, car: set) -> None:
        refusals, _ = check(decl, mem, car)
        fired = any(code in refusal for refusal in refusals)
        if fired != want:
            failures.append(
                f"{label}: expected the gate {'to refuse' if want else 'to accept'} "
                f"on {code!r} and it did not")

    a_member = sorted(declaration.get("recomputed") or {})[0]
    a_carried = sorted(carried & set(members))[0]

    # 1. A member the class declares that nothing accounts for.
    stripped = {k: dict(v or {}) for k, v in declaration.items()}
    stripped["recomputed"].pop(a_member)
    expect("an unclassified member", "does not say why", True,
           stripped, members, carried)

    # 2. A classified member the serializer now streams.
    expect("a classified member that is serialized", "rebuilt or", True,
           declaration, members, carried | {a_member})

    # 3. A non-empty lost block.
    lost = {k: dict(v or {}) for k, v in declaration.items()}
    lost["lost"] = {a_member: "a finding"}
    expect("a member recorded as lost", "read before it is written", True,
           lost, members, carried)

    # 4. A classification naming a member the class dropped.
    gone = dict(members)
    gone.pop(a_member)
    expect("a member the class no longer declares", "no longer declares", True,
           declaration, gone, carried)

    # 5. A classification with no evidence.
    bare = {k: dict(v or {}) for k, v in declaration.items()}
    bare["recomputed"][a_member] = ""
    expect("a classification with no evidence", "no evidence", True,
           bare, members, carried)

    # 6. A class key the gate does not know how to check.
    unknown = {k: dict(v or {}) for k, v in declaration.items()}
    unknown["probably_fine"] = {a_member: "a hunch"}
    expect("an unrecognised class", "not one of", True,
           unknown, members, carried)

    # 7. And the live tree, which must be granted.
    refusals, _ = check(declaration, members, carried)
    if refusals:
        failures.append("the live declaration: expected the gate to accept it "
                        f"and it refused: {refusals[0]}")

    # The serialized set has to be non-empty, or every check above passes
    # vacuously on a parse that found nothing.
    if not a_carried:
        failures.append("Soil::serialize parsed to no members at all")

    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list-unclassified", action="store_true",
                        help="print the unclassified members and their declared "
                             "types, ready to paste into the declaration")
    args = parser.parse_args()

    header = (GUESS_SOURCE / "framework" / "guess.h").read_text()
    source = (GUESS_SOURCE / "modules" / "soil.cpp").read_text()
    declaration = yaml.safe_load(DECLARATION.read_text()) or {}

    members = declared_members(header)
    carried = serialized_members(source)
    refusals, report = check(declaration, members, carried)
    broken = fixtures(declaration, members, carried)

    if args.list_unclassified:
        for member in report["unclassified"]:
            print(f"  {member}: \"\"   # {members[member]}")
        return 1 if report["unclassified"] else 0

    report.update({
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "biosphere/scripts/soil_restart_state_gate.py",
        "declaration": rel(DECLARATION),
        "granted": not refusals and not broken,
        "refusals": refusals,
        "fixtures_failed": broken,
    })
    GENERATED.mkdir(parents=True, exist_ok=True)
    out = GENERATED / "soil_restart_state_report.json"
    out.write_text(json.dumps(report, indent=2) + "\n")

    print(f"Soil declares {report['declared_members']} members: "
          f"{len(report['serialized'])} serialized, "
          f"{len(report['classified']['recomputed'])} rebuilt before first read, "
          f"{len(report['classified']['diagnostic'])} diagnostic, "
          f"{len(report['classified']['lost'])} lost")
    for refusal in refusals:
        print(f"\n  REFUSED: {refusal}")
    for broke in broken:
        print(f"\n  FIXTURE: {broke}")
    print(f"\nwrote {rel(out)}")
    return 1 if (refusals or broken) else 0


if __name__ == "__main__":
    sys.exit(main())
