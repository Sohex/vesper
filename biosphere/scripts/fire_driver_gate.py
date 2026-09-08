"""The fire driver gate: which fire arm may be armed, decided by capability
rather than by the name of the thing supplying the data.

Worldbuilding. Vesper is an invented planet; everything below is about the
vegetation model that simulates it -- which of its fire operators may run, and
what has to exist before one may.

`biosphere/config/fire_driver_contract.yaml` is the declaration and this module
is the enforcement. It is the third gate in this area and the three ask three
different questions, which is why they are not one:

  fire_gate.py            where this project's fire source DEPARTS from
                          mainline, and whether it still does
  fire_parameter_gate.py  what disposition every Earth CONSTANT the port
                          inherits carries
  fire_driver_gate.py     what each fire ARM requires, what this world can
                          supply, and which arm the model is actually taking

WHY A CAPABILITY TEST AT ALL. `biosphere/notes/fire-model-audit.md` finding 9:
`modules/driver.cpp:734` runs SIMFIRE's whole annual accounting whenever
`firemodel == BLAZE`, and `getsimfiredata` then opens an archive of human
population history and observed burned area, or aborts. The archive's lookup
keys on `floor(lon*2)/2 + 0.25`, so a Vesper gridcell's coordinates are valid
keys into an Earth half-degree archive: it would find a record, return the Earth
location's population and fire season, and the run would COMPLETE. THE FAILURE
MODE IS SILENT SUCCESS, so a test of the form "is the file there" passes exactly
when it should refuse. What is testable instead is whether the QUANTITY is
something this world has, which is what `external_archive` names and always
refuses.

It can fail:

  reference   a capability naming a forcing process the ecological forcing
              contract does not carry, or a parameter that
              `fire_parameters.yaml` does not register. The contract must not
              restate a field list or a constant it does not own, so what it
              names has to resolve
  disposition a `parameter` capability marked available whose constants do not
              all permit a run. A registered constant whose `central` is null
              is a DECLARED ABSENCE, and an absence does not satisfy a
              requirement
  refusal     an `external_archive` capability not marked refused. The kind
              exists to be always refused; a status that says otherwise is the
              declaration contradicting its own vocabulary
  arm         an arm whose armability, recomputed from its capabilities,
              disagrees with what the contract says about it. The verdict is
              derived and never declared, so this block cannot drift into an
              opinion
  effective   THE ONE THAT GUARDS AN ORDERING RATHER THAN A VALUE. The shipped
              instruction file sets `firemodel "BLAZE"` and `run_lpj_guess.py`
              writes `firemodel "GLOBFIRM"` after importing it. LPJ-GUESS's
              parser takes the LAST assignment, so that write order is the only
              thing disarming a refused arm. Moving the import below the
              override would arm BLAZE with no diagnostic and no change to any
              value. This reads both files and refuses if the override stops
              coming last, if the effective arm is not the expected one, or if
              the effective arm requires a capability this world refuses

Reduced fixtures run on every invocation, all but one built to be wrong in a
named way.

    python biosphere/scripts/fire_driver_gate.py            # status, exit 0
    python biosphere/scripts/fire_driver_gate.py --strict   # also refuses while
                                                            # the replacement arm
                                                            # is unbuildable

`--strict` refuses on one thing the default tolerates: VESPER_REPLACEMENT not
being armable. That is the project's actual state -- fire-1 through fire-7 are
open -- so it is a residual to report rather than a defect to fail on, and
`--strict` is for a caller that wants the arm to exist.

IT READS AND RUNS NOTHING. Three YAML declarations, one generated instruction
file and one Python source, all read as text.
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
DECLARATION = COMPONENT_ROOT / "config" / "fire_driver_contract.yaml"
FORCING = COMPONENT_ROOT / "config" / "ecological_forcing_contract.yaml"
PARAMETERS = COMPONENT_ROOT / "config" / "fire_parameters.yaml"
SHIPPED_INS = COMPONENT_ROOT / "generated" / "vesper_pfts.ins"
RUNNER = COMPONENT_ROOT / "scripts" / "run_lpj_guess.py"
REPORT = GENERATED / "fire_driver_gate_report.json"

STATUSES = ("available", "blocked", "refused")
ALWAYS_REFUSED = ("external_archive",)


def _finding(kind: str, what: str, detail: str) -> dict:
    return {"kind": kind, "what": what, "detail": detail}


def _forcing_processes(forcing: dict) -> dict[str, list]:
    """The forcing contract's process registry, keyed by issue id.

    Read rather than restated: the contract owns the field rows and states each
    edge once, which is the rule `config/pipeline.yaml` applies to `needs`.
    """
    out = {}
    for row in (forcing or {}).get("processes") or []:
        if isinstance(row, dict):
            name = row.get("name") or row.get("id")
            if name:
                out[str(name)] = row.get("needs") or row.get("fields") or []
    return out


def check_references(contract: dict, forcing: dict, params: dict) -> list[dict]:
    """Everything the contract names in another declaration resolves there."""
    bad = []
    processes = _forcing_processes(forcing)
    registered = set((params or {}).get("parameters") or {})
    declared_fields = {row.get("name") for row in (forcing or {}).get("fields") or []
                       if isinstance(row, dict)}
    for name, cap in (contract.get("capabilities") or {}).items():
        if not isinstance(cap, dict):
            bad.append(_finding("reference", name, "is not a mapping"))
            continue
        if cap.get("status") not in STATUSES:
            bad.append(_finding("reference", name,
                                f"status {cap.get('status')!r} is not one of "
                                f"{STATUSES}"))
        process = cap.get("forcing_process")
        if process is not None:
            if process not in processes:
                bad.append(_finding(
                    "reference", name,
                    f"names forcing process {process!r}, which "
                    f"ecological_forcing_contract.yaml does not carry, so the "
                    f"field list it stands for cannot be read"))
            else:
                missing = [f for f in processes[process]
                           if f not in declared_fields]
                if missing:
                    bad.append(_finding(
                        "reference", name,
                        f"forcing process {process!r} needs {missing}, which "
                        f"the contract's own field rows do not declare"))
        for parameter in cap.get("parameters") or []:
            if parameter not in registered:
                bad.append(_finding(
                    "reference", name,
                    f"names parameter {parameter!r}, which "
                    f"fire_parameters.yaml does not register"))
    return bad


def check_dispositions(contract: dict, params: dict) -> list[dict]:
    """A parameter capability is available only if its constants permit a run.

    A registered constant whose `central` is null is a DECLARED ABSENCE -- the
    shape `phosphorus_volatilised_fraction` and `fuel_currency` use -- and an
    absence cannot satisfy a requirement. This is what stops the contract
    marking a capability available because the register merely MENTIONS its
    constants.
    """
    bad = []
    registered = (params or {}).get("parameters") or {}
    for name, cap in (contract.get("capabilities") or {}).items():
        if not isinstance(cap, dict) or cap.get("kind") != "parameter":
            continue
        if cap.get("status") != "available":
            continue
        absent = [p for p in cap.get("parameters") or []
                  if isinstance(registered.get(p), dict)
                  and registered[p].get("central") is None]
        if absent:
            bad.append(_finding(
                "disposition", name,
                f"is marked available and its constants {absent} declare a "
                f"null central, which is an absence rather than a value"))
    return bad


def check_refusals(contract: dict) -> list[dict]:
    """A kind that exists to be refused is still refused everywhere."""
    bad = []
    for name, cap in (contract.get("capabilities") or {}).items():
        if not isinstance(cap, dict):
            continue
        if cap.get("kind") in ALWAYS_REFUSED and cap.get("status") != "refused":
            bad.append(_finding(
                "refusal", name,
                f"is an {cap['kind']} and its status is "
                f"{cap.get('status')!r}. That kind is always refused: this "
                f"world has no observational record, and a file that answers "
                f"anyway is the silent success the contract exists to stop"))
        if cap.get("status") == "refused" and not cap.get("refused_because"):
            bad.append(_finding("refusal", name,
                                "is refused and says nothing about why"))
        if cap.get("status") == "blocked" and not cap.get("owner"):
            bad.append(_finding(
                "blocked", name,
                "is blocked and names no owner. Blocked means there is a path "
                "and somebody on it; without an owner it is refused wearing a "
                "softer word"))
    return bad


def armability(contract: dict) -> dict[str, dict]:
    """Each arm's verdict, DERIVED from its capabilities and never declared."""
    caps = contract.get("capabilities") or {}
    out = {}
    for arm, spec in (contract.get("arms") or {}).items():
        required = (spec or {}).get("requires") or []
        unknown = [c for c in required if c not in caps]
        refused = [c for c in required
                   if caps.get(c, {}).get("status") == "refused"]
        blocked = [c for c in required
                   if caps.get(c, {}).get("status") == "blocked"]
        if unknown:
            verdict = "undeclared"
        elif refused:
            verdict = "refused"
        elif blocked:
            verdict = "blocked"
        else:
            verdict = "armable"
        out[arm] = {"verdict": verdict, "refused_by": refused,
                    "blocked_by": blocked, "undeclared": unknown}
    return out


def check_arms(contract: dict, verdicts: dict) -> list[dict]:
    bad = []
    for arm, got in verdicts.items():
        if got["undeclared"]:
            bad.append(_finding(
                "arm", arm,
                f"requires {got['undeclared']}, which the capabilities block "
                f"does not declare"))
    return bad


def effective_arm(contract: dict, root: Path) -> tuple[str | None, list[dict]]:
    """Which arm the model actually takes, from the two files that decide it.

    LPJ-GUESS's parser takes the LAST assignment of a key, so this is a
    property of write order across two files rather than of either one. The
    check compares positions rather than trusting the comment beside either,
    because a comment is what survives the edit that breaks the ordering.
    """
    bad = []
    spec = contract.get("effective_arm") or {}
    shipped_path = root / (spec.get("shipped_in") or "")
    override_path = root / (spec.get("override_in") or "")
    for label, path in (("shipped_in", shipped_path),
                        ("override_in", override_path)):
        if not path.is_file():
            bad.append(_finding("effective", label,
                                f"{path} is not a file, so the effective arm "
                                f"cannot be determined"))
            return None, bad

    shipped_text = shipped_path.read_text(encoding="utf-8", errors="replace")
    # An LPJ-GUESS instruction comment runs from `!` to end of line.
    shipped_code = re.sub(r"![^\n]*", "", shipped_text)
    shipped = re.findall(r'firemodel\s+"([A-Za-z]+)"', shipped_code)
    if not shipped:
        bad.append(_finding("effective", "shipped_in",
                            f"{spec['shipped_in']} assigns no firemodel, so "
                            f"the override has nothing to override and the "
                            f"contract's account of the ordering is stale"))
    elif shipped[-1] != spec.get("shipped_default"):
        bad.append(_finding(
            "effective", "shipped_in",
            f"{spec['shipped_in']} now ships firemodel {shipped[-1]!r} where "
            f"the contract records {spec.get('shipped_default')!r}"))

    runner_text = override_path.read_text(encoding="utf-8", errors="replace")
    override = list(re.finditer(r'^firemodel\s+"([A-Za-z]+)"', runner_text,
                                re.M))
    imports = list(re.finditer(r'^import\s+"', runner_text, re.M))
    if not override:
        bad.append(_finding(
            "effective", "override_in",
            f"{spec['override_in']} no longer writes a firemodel line, so the "
            f"shipped {spec.get('shipped_default')!r} is what the model takes"))
        return (shipped[-1] if shipped else None), bad
    if not imports:
        bad.append(_finding("effective", "override_in",
                            f"{spec['override_in']} imports no instruction "
                            f"file, so the ordering this rests on is gone"))
        return override[-1].group(1), bad

    last_override, last_import = override[-1], imports[-1]
    if last_override.start() < last_import.start():
        bad.append(_finding(
            "effective", "ordering",
            f"{spec['override_in']} writes firemodel "
            f"{last_override.group(1)!r} BEFORE its import, and the parser "
            f"takes the last assignment, so the shipped "
            f"{spec.get('shipped_default')!r} wins. This is the failure the "
            f"check exists for: an ordering moved, no value changed, and a "
            f"refused arm is armed"))
        return (shipped[-1] if shipped else None), bad

    taken = last_override.group(1)
    if taken != spec.get("expected"):
        bad.append(_finding("effective", "expected",
                            f"the effective arm is {taken!r} where the "
                            f"contract expects {spec.get('expected')!r}"))
    return taken, bad


def check(contract: dict, forcing: dict, params: dict,
          root: Path) -> tuple[list[dict], dict, str | None]:
    findings = []
    findings.extend(check_references(contract, forcing, params))
    findings.extend(check_dispositions(contract, params))
    findings.extend(check_refusals(contract))
    verdicts = armability(contract)
    findings.extend(check_arms(contract, verdicts))
    taken, bad = effective_arm(contract, root)
    findings.extend(bad)
    if taken and verdicts.get(taken, {}).get("verdict") == "refused":
        findings.append(_finding(
            "effective", taken,
            f"is the arm the model takes and it is REFUSED by "
            f"{verdicts[taken]['refused_by']}. A refused arm running is the "
            f"whole failure this contract exists to prevent"))
    return findings, verdicts, taken


def _fixtures(root: Path) -> list[dict]:
    """Reduced declarations, all but one built to be wrong in a named way.

    The ordering fixtures are the ones that matter. `effective_arm` guards a
    property of WRITE ORDER, and the way that check rots is by silently
    matching nothing -- a renamed key, a changed quote style -- after which it
    reports clean forever. So the fixtures write temporary files whose orders
    are known and require the verdict to follow the order.
    """
    import tempfile

    forcing = {"fields": [{"name": "convective_precipitation"},
                          {"name": "air_temperature"}],
               "processes": [{"name": "fire-1",
                              "needs": ["convective_precipitation"]}]}
    params = {"parameters": {"good": {"central": 1.0},
                             "absent": {"central": None}}}

    def contract(**over):
        base = {
            "capabilities": {
                "ok": {"kind": "derived_state", "status": "available"},
                "arch": {"kind": "external_archive", "status": "refused",
                         "refused_because": "no observational record"},
            },
            "arms": {"A": {"requires": ["ok"]}},
            "effective_arm": {"expected": "GLOBFIRM", "shipped_default": "BLAZE",
                              "shipped_in": "shipped.ins",
                              "override_in": "runner.py"},
        }
        base.update(over)
        return base

    def tree(tmp, shipped_text, runner_text):
        d = Path(tmp)
        (d / "shipped.ins").write_text(shipped_text)
        (d / "runner.py").write_text(runner_text)
        return d

    GOOD_SHIPPED = 'vegmode "cohort"\nfiremodel "BLAZE"\n'
    GOOD_RUNNER = 'import "{pfts}"\n\nfiremodel "GLOBFIRM"\n'
    BAD_ORDER = 'firemodel "GLOBFIRM"\n\nimport "{pfts}"\n'

    cases = []

    def case(label, findings, expected):
        found = sorted({f["kind"] for f in findings})
        ok = (found == []) if expected == "no findings" else (expected in found)
        cases.append({"fixture": label, "expected": expected,
                      "found": found or ["nothing"], "pass": ok})

    with tempfile.TemporaryDirectory() as tmp:
        d = tree(tmp, GOOD_SHIPPED, GOOD_RUNNER)
        case("a clean contract reports nothing",
             check(contract(), forcing, params, d)[0], "no findings")

        case("a capability naming a forcing process the contract lacks",
             check(contract(capabilities={
                 "x": {"kind": "forcing_field", "status": "blocked",
                       "owner": "fire-1", "forcing_process": "fire-99"}},
                 arms={}), forcing, params, d)[0], "reference")

        case("a capability naming an unregistered parameter",
             check(contract(capabilities={
                 "x": {"kind": "parameter", "status": "available",
                       "parameters": ["not_registered"]}},
                 arms={}), forcing, params, d)[0], "reference")

        case("a parameter capability available on a null central",
             check(contract(capabilities={
                 "x": {"kind": "parameter", "status": "available",
                       "parameters": ["absent"]}},
                 arms={}), forcing, params, d)[0], "disposition")

        case("an external archive not marked refused",
             check(contract(capabilities={
                 "x": {"kind": "external_archive", "status": "available"}},
                 arms={}), forcing, params, d)[0], "refusal")

        case("a blocked capability with no owner",
             check(contract(capabilities={
                 "x": {"kind": "derived_state", "status": "blocked"}},
                 arms={}), forcing, params, d)[0], "blocked")

        case("an arm requiring an undeclared capability",
             check(contract(arms={"A": {"requires": ["nonexistent"]}}),
                   forcing, params, d)[0], "arm")

    with tempfile.TemporaryDirectory() as tmp:
        d = tree(tmp, GOOD_SHIPPED, BAD_ORDER)
        case("the override moved above the import, arming the shipped default",
             check(contract(), forcing, params, d)[0], "effective")

    with tempfile.TemporaryDirectory() as tmp:
        d = tree(tmp, GOOD_SHIPPED, 'import "{pfts}"\n')
        case("the override deleted entirely",
             check(contract(), forcing, params, d)[0], "effective")

    with tempfile.TemporaryDirectory() as tmp:
        d = tree(tmp, GOOD_SHIPPED, 'import "{pfts}"\n\nfiremodel "BLAZE"\n')
        case("the effective arm is a refused arm",
             check(contract(capabilities={
                 "ok": {"kind": "derived_state", "status": "available"},
                 "arch": {"kind": "external_archive", "status": "refused",
                          "refused_because": "no observational record"}},
                 arms={"BLAZE": {"requires": ["arch"]}}),
                 forcing, params, d)[0], "effective")

    with tempfile.TemporaryDirectory() as tmp:
        d = tree(tmp, 'vegmode "cohort"\n! firemodel "BLAZE"\n', GOOD_RUNNER)
        case("the shipped default is only in a comment",
             check(contract(), forcing, params, d)[0], "effective")

    return cases


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--strict", action="store_true",
                        help="also refuse while the replacement arm is unbuildable")
    parser.add_argument("--json", action="store_true", help="report only")
    args = parser.parse_args()

    contract = yaml.safe_load(DECLARATION.read_text(encoding="utf-8"))
    forcing = yaml.safe_load(FORCING.read_text(encoding="utf-8"))
    params = yaml.safe_load(PARAMETERS.read_text(encoding="utf-8"))

    fixtures = _fixtures(PROJECT_ROOT)
    findings, verdicts, taken = check(contract, forcing, params, PROJECT_ROOT)

    caps = contract.get("capabilities") or {}
    by_status: dict[str, list[str]] = {}
    for name, cap in caps.items():
        by_status.setdefault(cap.get("status", "?"), []).append(name)

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declaration": rel(DECLARATION),
        "capabilities": {k: sorted(v) for k, v in sorted(by_status.items())},
        "arms": verdicts,
        "effective_arm": taken,
        "findings": findings,
        "fixtures": fixtures,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print("Which fire arm may be armed, decided by capability.\n")
        broken = [f for f in fixtures if not f["pass"]]
        print(f"  fixtures: {len(fixtures) - len(broken)} of {len(fixtures)} "
              f"got their verdict")
        for c in broken:
            print(f"    BROKEN: {c['fixture']}: expected {c['expected']}, "
                  f"found {c['found']}")
        print("\n  capabilities:")
        for status in STATUSES:
            names = report["capabilities"].get(status) or []
            if names:
                print(f"    {status:<10} {len(names):>2}  {', '.join(names)}")
        print("\n  arms:")
        for arm, v in verdicts.items():
            reason = ""
            if v["refused_by"]:
                reason = f"  refused by {', '.join(v['refused_by'])}"
            elif v["blocked_by"]:
                reason = f"  blocked by {', '.join(v['blocked_by'])}"
            mark = "<- taken" if arm == taken else ""
            print(f"    {arm:<20} {v['verdict']:<10}{reason} {mark}")
        if findings:
            print(f"\n  {len(findings)} finding(s):")
            for f in findings:
                print(f"    [{f['kind']}] {f['what']}: {f['detail']}")
        else:
            print(f"\n  every capability resolves where it is named, every "
                  f"archive is refused, and the arm the model takes "
                  f"({taken}) is armable")

    replacement = verdicts.get("VESPER_REPLACEMENT", {})
    if broken or findings:
        return 1
    if args.strict and replacement.get("verdict") != "armable":
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
