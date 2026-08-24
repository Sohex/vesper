"""The BVOC activation gate: what has to be declared before the simulated
biosphere's volatile organic source may be switched on, and what an activated
run has to emit before its output is accepted.

Worldbuilding. Vesper is an invented planet; everything below is about the
simulation of it -- a vegetation model's dormant emission module, the
declaration file that governs it, and the output tables it would write.

`ifbvoc 1` reads like one switch and is four models: the simulated plants'
production of speciated volatile carbon, its oxidation in an atmosphere whose
oxidants this project has not derived, the aerosol that oxidation would make,
and what that aerosol would do to radiation and to the model's clouds. Only the
first exists in the vendored source, and it is Earth-calibrated. The argument
is in `biosphere/notes/bvoc-activation-contract.md`.

This module is the enforcement. It reads `biosphere/config/bvoc.yaml` and
refuses activation, by name, for every precondition still carrying the
`undeclared` sentinel -- and it refuses for a declared bracket that is narrower
than the declared floor, and for a cloud branch carried on the prescribed
aerosol optical-depth field.

It is fail-closed in one direction only. A biosphere run with the volatile
source OFF is a correct run, so with `requested: false` the gate reports what is
undeclared and exits 0. Only a request to activate can be refused.

    python biosphere/scripts/bvoc_gate.py              # status, exit 0
    python biosphere/scripts/bvoc_gate.py --check-run biosphere/runs/<run_id>

`run_lpj_guess.py` calls `granted()` and writes `ifbvoc 0` whenever activation
has not been granted, which is every case that exists today.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import CONFIG, GENERATED, GUESS_SOURCE, PROJECT_ROOT
from paths import rel  # noqa: E402

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
DECLARATION = COMPONENT_ROOT / "config" / "bvoc.yaml"
REPORT = GENERATED / "bvoc_gate_report.json"

UNDECLARED = "undeclared"

# The nine monoterpene compounds the vendored module carries internally, from
# framework/parameters.h. The output tables collapse them to an endocyclic and
# an other class, which is the collapse finding 1 of the audit is about.
MONOTERPENE_COMPOUNDS = ("APIN", "BPIN", "LIMO", "MYRC", "SABI", "CAMP",
                         "TRIC", "TBOC", "OTHR")

# What the vendored module reports its volatile fluxes in, against the kgC/m2
# the plant carbon pools use. Six orders of magnitude between the emission and
# the ledger it has to be subtracted from is not a rounding difference, so the
# unit is carried explicitly through every acceptance check rather than assumed.
FLUX_UNIT = "mgC/m2"
POOL_UNIT = "kgC/m2"

# What an ACTIVATED run has to retain, and what `check_run` accepts or rejects.
# One monthly table per compound, because the contract's whole point is that a
# lumped terpene mass cannot reconstruct a chemistry branch. None of these
# output parameters exists in the vendored fork yet, which is BVOC-11 and is
# also why nothing can request them: the gate refuses first.
ACTIVATED_OUTPUTS = ("miso.out",) + tuple(
    f"mmon_{compound.lower()}.out" for compound in MONOTERPENE_COMPOUNDS)

# The declaration fields whose sentinel refuses, and the bd issue that supplies
# each. The order is the order the refusals print in: source first, because a
# source that cannot emit compound-resolved flux makes every later declaration
# unusable, then traits, then the oxidant environment, then clouds.
PRECONDITIONS = (
    ("source", "compound_resolved_output",
     "monthly flux per monoterpene compound", "BVOC-11"),
    ("source", "patch_resolved_output",
     "flux carrying the patch axis, not only the PFT axis", "BVOC-11"),
    ("source", "storage_restart_fixture",
     "a fixture comparing the monoterpene storage pool across a restart",
     "BVOC-3"),
    ("source", "carbon_closure",
     "emitted and stored volatile carbon removed from a named plant pool",
     "BVOC-3"),
    ("source", "environmental_response_port",
     "the rotation, photoperiod, year-length and leaf-energy dependencies "
     "ported off Earth's constants", "BVOC-4"),
    ("source", "missing_scope_bracket",
     "a bracket for the volatile classes the retained module does not "
     "represent", "BVOC-5"),
    ("traits", "registry",
     "the covarying plant trait registry the emission capacities live in",
     "PCAR-5"),
    ("traits", "emitter_fraction", "the emitting fraction per simulated PFT",
     "BVOC-2"),
    ("traits", "speciation", "the per-compound split per simulated PFT",
     "BVOC-2"),
    ("traits", "seasonality",
     "seasonal onset and decay in this world's photoperiod and year",
     "BVOC-2"),
    ("traits", "storage_strategy",
     "emission from synthesis, from storage, or both, per simulated PFT",
     "BVOC-2"),
    ("traits", "declared_bracket_factor",
     "the declared uncertainty brackets on the emission capacities",
     "BVOC-2"),
    ("oxidants", "interface",
     "the reduced chemistry interface, with its inputs and its returns",
     "BVOC-6"),
    ("oxidants", "instrument",
     "what derived the bracket ends, and over what domain", "BVOC-6"),
    ("oxidants", "bracket",
     "both ends of the oxidant bracket, and they must differ", "BVOC-6"),
)


class Refusal:
    """One named reason activation cannot be granted.

    `code` is what the error says and what a later report is grepped for; it is
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


def read_declaration(path: Path = DECLARATION) -> dict:
    if not path.is_file():
        # The declaration's absence is itself a refusal, not an empty default.
        raise SystemExit(f"BVOC-DECLARATION-MISSING: {rel(path)} does not "
                         f"exist. There is no default declaration.")
    return yaml.safe_load(path.read_text())


def _get(declaration: dict, section: str, field: str):
    return (declaration.get(section) or {}).get(field, UNDECLARED)


def _declared(value) -> bool:
    return value is not None and value != UNDECLARED


def probes() -> dict:
    """Read-only evidence attached to the refusals.

    These do not decide anything: the declaration decides, and a probe that
    stops matching because someone reworded a line must not be able to grant
    what a declaration still refuses. What they are for is putting the current
    state of the vendored source next to the refusal, so that the reason is
    checkable against the artifact rather than against this file.
    """
    found: dict = {}

    common = GUESS_SOURCE / "modules" / "commonoutput.cpp"
    if common.is_file():
        text = common.read_text(errors="replace")
        covered = [c for c in MONOTERPENE_COMPOUNDS
                   if re.search(rf'file_\w*{c.lower()}\w*"', text)]
        found["compound_resolved_output_tables"] = {
            "compounds_internal": len(MONOTERPENE_COMPOUNDS),
            "compounds_with_an_output_table": len(covered),
            "collapsed_to": sorted(set(re.findall(r'"(file_[ma]mon_mt\d)"', text))),
        }

    bvoc = GUESS_SOURCE / "modules" / "bvoc.cpp"
    if bvoc.is_file():
        text = bvoc.read_text(errors="replace")
        found["earth_constants_in_the_emission_module"] = {
            "daylength_over_24h": bool(re.search(r"daylength\s*\*\s*PI\s*/\s*24", text)),
            "storage_time_constant_365_days": bool(
                re.search(r"tcstor_max\s*=\s*365", text)),
            "fixed_air_density": bool(re.search(r"rhoair\s*=\s*1\.204", text)),
        }

    core = GUESS_SOURCE / "framework" / "guess.cpp"
    if core.is_file():
        found["monoterpene_storage_in_serializer"] = bool(
            re.search(r"&\s*monstor", core.read_text(errors="replace")))

    radmod = (PROJECT_ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim"
              / "src" / "radmod.f90")
    if radmod.is_file():
        text = radmod.read_text(errors="replace")
        naersp = re.search(r"parameter\(NAERSP\s*=\s*(\d+)\)", text)
        found["climate_model_aerosol_carrier"] = {
            "radiative_species_slots": int(naersp.group(1)) if naersp else None,
            "aerosol_to_droplet_pathway": bool(
                re.search(r"droplet|\bccn\b", text, re.IGNORECASE)),
            "shortwave_aerosol_in_cloudy_fraction": not bool(
                re.search(r"zaerr1\(:,jlev\)\*\(1\.-dcc", text)),
        }
    return found


def evaluate(declaration: dict, planet: dict) -> list[Refusal]:
    """Every reason activation cannot be granted, in declaration order."""
    refusals: list[Refusal] = []

    for section, field, what, issue in PRECONDITIONS:
        if not _declared(_get(declaration, section, field)):
            refusals.append(Refusal(
                f"BVOC-UNDECLARED-{section.upper()}-{field.upper().replace('_', '-')}",
                f"{what} is undeclared", issue))

    # A declared bracket narrower than the declared floor. This is the one
    # check that can fail on a value someone has supplied rather than on a
    # sentinel, and it is the direction that matters: a narrow bracket is a
    # claim to know this world's plants better than the read Earth studies know
    # Earth's, and it arrives looking like progress.
    floors = (declaration.get("traits") or {}).get("minimum_bracket_factor") or {}
    declared = (declaration.get("traits") or {}).get("declared_bracket_factor")
    if isinstance(declared, dict):
        for name, floor in floors.items():
            value = declared.get(name)
            if value is None:
                refusals.append(Refusal(
                    "BVOC-BRACKET-MISSING",
                    f"declared_bracket_factor has no entry for {name}", "BVOC-2"))
            elif float(value) < float(floor):
                refusals.append(Refusal(
                    "BVOC-BRACKET-TOO-NARROW",
                    f"{name} declared at factor {value}, below the declared "
                    f"floor of {floor}", "BVOC-2"))

    # The simulated atmosphere's ozone and methane are prescribed from a
    # photochemical calculation that holds modern Earth's biogenic surface
    # fluxes fixed. Driving a volatile source from this world's own simulated
    # vegetation while those stay fixed is not one atmospheric state, and the
    # inconsistency is the thing that must not be combined silently.
    if not (declaration.get("oxidants") or {}).get("trace_gas_state_closed"):
        prescribed = [k for k in ("pCH4_bar", "ozone")
                      if k in (planet.get("atmosphere") or {})]
        if prescribed:
            refusals.append(Refusal(
                "BVOC-TRACE-GAS-STATE-NOT-CLOSED",
                "config/planet.yaml prescribes " + ", ".join(prescribed) +
                " from a fixed-Earth-flux photochemical calculation, so a "
                "vegetation-derived volatile source cannot be combined with it "
                "as one atmospheric state", "BVOC-6 / WET-10"))

    # The cloud arm. `none` is the registered default and needs nothing.
    branch = (declaration.get("clouds") or {}).get("branch", "none")
    carrier = (declaration.get("clouds") or {}).get("carrier", UNDECLARED)
    if branch not in ("none", "offline_ccn_bound", "activation_port"):
        refusals.append(Refusal(
            "BVOC-CLOUD-BRANCH-UNREGISTERED",
            f"clouds.branch is {branch!r}, which is not one of the three arms "
            "the pre-registration fixed", "BVOC-9"))
    elif branch != "none":
        if not _declared(carrier):
            refusals.append(Refusal(
                "BVOC-CLOUD-CARRIER-UNDECLARED",
                f"clouds.branch is {branch!r} with no declared carrier",
                "BVOC-9"))
        elif re.search(r"optical.?depth|\baod\b|ddustcol|dustsc", str(carrier),
                       re.IGNORECASE):
            refusals.append(Refusal(
                "BVOC-CLOUD-OPTICAL-PROXY",
                "a cloud arm carried on the prescribed aerosol optical-depth "
                "field: a column optical depth is mass and optics and does not "
                "determine particle number, so scaling it is a fourth arm the "
                "pre-registration did not register", "BVOC-9"))
        if branch == "activation_port" and not _declared(
                (declaration.get("clouds") or {}).get("materiality_verdict")):
            refusals.append(Refusal(
                "BVOC-CLOUD-VERDICT-MISSING",
                "the activation port is the arm that is built only on the "
                "registered materiality verdict, and no verdict is declared",
                "BVOC-9"))
    return refusals


def granted(declaration: dict | None = None,
            planet: dict | None = None) -> tuple[bool, list[Refusal]]:
    """Whether an activated run may proceed, and why not.

    `run_lpj_guess.py` calls this. A declaration that does not request
    activation is granted nothing and refuses nothing: the run writes
    `ifbvoc 0` and is a correct run.
    """
    declaration = declaration if declaration is not None else read_declaration()
    planet = planet if planet is not None else yaml.safe_load(CONFIG.read_text())
    refusals = evaluate(declaration, planet)
    if not declaration.get("requested"):
        return False, refusals
    return not refusals, refusals


def require(declaration: dict | None = None,
            planet: dict | None = None) -> bool:
    """`granted`, but a refused REQUEST exits rather than returning False."""
    declaration = declaration if declaration is not None else read_declaration()
    ok, refusals = granted(declaration, planet)
    if declaration.get("requested") and not ok:
        lines = "\n  ".join(str(r) for r in refusals)
        raise SystemExit(
            f"{rel(DECLARATION)} requests BVOC activation and the gate refuses "
            f"it for {len(refusals)} reasons:\n  {lines}\n\n"
            "Each is a declaration to make in that file, with the bd issue that "
            "supplies it. Nothing here has a default value.")
    return ok


# --------------------------------------------------------------------------
# Acceptance, for output that does not exist yet.

def check_run(run_dir: Path) -> list[Refusal]:
    """Reject a run's volatile output rather than accept it by default.

    Every check has an answer that is wrong rather than merely different: a
    cell present in the productivity output and absent here is a missing cell,
    a negative emission is impossible, a monthly table without one column per
    simulated month is not this world's year, and a run whose manifest carries
    no forcing hash cannot be attributed to a climate.
    """
    refusals: list[Refusal] = []
    run_dir = Path(run_dir)

    expected = list(ACTIVATED_OUTPUTS)
    missing = [name for name in expected if not (run_dir / name).is_file()]
    if missing:
        refusals.append(Refusal(
            "BVOC-OUTPUT-MISSING",
            f"{len(missing)} of {len(expected)} compound-resolved monthly "
            f"tables absent from {rel(run_dir)}: " + ", ".join(missing[:4]) +
            ("..." if len(missing) > 4 else ""), "BVOC-11"))

    reference = run_dir / "anpp.out"
    reference_cells: set[tuple[str, str]] = set()
    if reference.is_file():
        for line in reference.read_text().splitlines()[1:]:
            parts = line.split()
            if len(parts) >= 2:
                reference_cells.add((parts[0], parts[1]))
    else:
        refusals.append(Refusal(
            "BVOC-REFERENCE-MISSING",
            f"{rel(run_dir)}/anpp.out is absent, so the simulated cell set the "
            "volatile output has to cover is unknown"))

    months = _simulated_months()
    for name in expected:
        path = run_dir / name
        if not path.is_file():
            continue
        lines = path.read_text().splitlines()
        if not lines:
            refusals.append(Refusal("BVOC-OUTPUT-EMPTY", f"{name} is empty"))
            continue
        header = lines[0].split()
        cells: set[tuple[str, str]] = set()
        for line in lines[1:]:
            parts = line.split()
            if len(parts) < 3:
                continue
            cells.add((parts[0], parts[1]))
            for token in parts[2:]:
                try:
                    value = float(token)
                except ValueError:
                    refusals.append(Refusal(
                        "BVOC-OUTPUT-UNREADABLE",
                        f"{name} carries a non-numeric emission {token!r}"))
                    break
                if math.isnan(value) or math.isinf(value):
                    refusals.append(Refusal(
                        "BVOC-OUTPUT-NOT-FINITE",
                        f"{name} carries a non-finite emission"))
                    break
                if value < 0.0:
                    refusals.append(Refusal(
                        "BVOC-OUTPUT-IMPOSSIBLE",
                        f"{name} carries a negative emission in {FLUX_UNIT}"))
                    break
        if reference_cells and cells != reference_cells:
            refusals.append(Refusal(
                "BVOC-OUTPUT-INCOMPLETE",
                f"{name} covers {len(cells)} simulated cells against "
                f"{len(reference_cells)} in anpp.out"))
        if months is None:
            refusals.append(Refusal(
                "BVOC-CALENDAR-UNKNOWN",
                "biosphere/generated/vesper.h is absent, so how many months "
                "the simulated year has is unknown and a monthly emission "
                "table cannot be accepted"))
        elif len(header) >= 2 and len(header) - 2 < months:
            refusals.append(Refusal(
                "BVOC-OUTPUT-CALENDAR",
                f"{name} has {len(header) - 2} data columns against {months} "
                "months in the simulated year"))

    manifest = run_dir / "run_manifest.json"
    if not manifest.is_file():
        refusals.append(Refusal(
            "BVOC-PROVENANCE-MISSING",
            f"{rel(run_dir)}/run_manifest.json is absent, so the forcing this "
            "volatile flux came from is unrecorded"))
    else:
        record = json.loads(manifest.read_text())
        if not record.get("source_build"):
            refusals.append(Refusal(
                "BVOC-PROVENANCE-NO-BUILD",
                "the run manifest carries no source_build"))
        if not (record.get("inputs") or {}).get("driver", {}).get("sha256"):
            refusals.append(Refusal(
                "BVOC-PROVENANCE-NO-FORCING-HASH",
                "the run manifest carries no hash of the climate forcing"))
    return refusals


def _simulated_months() -> int | None:
    """Months in the simulated year, from the generated header.

    Twelve of them, and each is about half an Earth month, because this world's
    year is shorter. That is why the count is read rather than assumed: a
    monthly emission column here is not an Earth month and the two must not be
    compared without the conversion the run manifest records.
    """
    header = GENERATED / "vesper.h"
    if not header.is_file():
        return None
    match = re.search(r"VESPER_MONTH_LENGTHS\s*\{([^}]*)\}", header.read_text())
    return len([p for p in match.group(1).split(",") if p.strip()]) if match else None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-run", type=Path, default=None,
                        help="accept or reject one run's volatile output")
    parser.add_argument("--declaration", type=Path, default=DECLARATION)
    args = parser.parse_args()

    declaration = read_declaration(args.declaration)
    planet = yaml.safe_load(CONFIG.read_text())
    refusals = evaluate(declaration, planet)
    requested = bool(declaration.get("requested"))
    evidence = probes()

    run_refusals: list[Refusal] = []
    if args.check_run is not None:
        run_refusals = check_run(args.check_run)

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "generator": "biosphere/scripts/bvoc_gate.py",
        "declaration": str(rel(args.declaration)),
        "declaration_sha256": hashlib.sha256(
            args.declaration.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "requested": requested,
        "granted": requested and not refusals,
        "flux_unit": FLUX_UNIT,
        "plant_pool_unit": POOL_UNIT,
        "months_in_the_simulated_year": _simulated_months(),
        "refusals": [r.as_dict() for r in refusals],
        "run_checked": str(args.check_run) if args.check_run else None,
        "run_refusals": [r.as_dict() for r in run_refusals],
        "evidence": evidence,
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")

    state = "REQUESTED" if requested else "not requested"
    print(f"BVOC activation {state}; {len(refusals)} preconditions unmet")
    for refusal in refusals:
        print(f"  {refusal}")
    if args.check_run is not None:
        print(f"\n{rel(args.check_run)}: {len(run_refusals)} acceptance failures")
        for refusal in run_refusals:
            print(f"  {refusal}")
    print(f"\nwrote {rel(REPORT)}")

    if requested and refusals:
        sys.exit(1)
    if run_refusals:
        sys.exit(1)


if __name__ == "__main__":
    main()
