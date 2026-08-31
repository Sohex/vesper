"""Run LPJ-GUESS over Vesper, reproducibly and with its provenance recorded.

Every run before this one was hand-assembled in a scratch directory with a
hand-written instruction file, which is fine for a smoke test and not fine for
anything whose numbers get quoted. This is the equivalent of
`exoplasim/scripts/run_exoplasim.py`: it generates the instruction file, sets up
the rank directories, runs, merges the per-rank output and writes a manifest
pinning every input by hash.

    python biosphere/scripts/run_lpj_guess.py --nyear 50
    python biosphere/scripts/run_lpj_guess.py --ranks 16 --label carved-k25v

Bulk output lands in `runs/<run_id>/`, which is not tracked, in the same way
`exoplasim/runs/` is not. The manifest and a summary go to
`analysis/<run_id>/`, which is.

Three things this handles that catch people out:

**Absolute paths everywhere.** In parallel mode each rank chdirs into its own
`runN/` before reading anything, so a relative instruction path or a relative
`import` silently fails with "could not open ... for input".

**The rank split lives in the input module.** LPJ-GUESS normally splits work by
pre-splitting the gridlist file. `vesperinput` has no gridlist, so it strides the
cells itself; nothing here needs to divide anything.

**Fire and the weather generator.** `global.ins` selects BLAZE, which wants a
SimFIRE input file built from Earth observations, and GWGEN, which wants
sub-daily statistics we do not have. GLOBFIRM and INTERP are what the shipped
demo uses for the same reason.
"""

from __future__ import annotations

import argparse
import hashlib
import uuid
import json
import platform
import shutil
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import (CONFIG, GENERATED, GUESS_BINARY, GUESS_SOURCE, PROJECT_ROOT,
                    RUNS)

import orbit

COMPONENT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = COMPONENT_ROOT / "analysis"
SEED_CONFIG = COMPONENT_ROOT / "config" / "stochastic_seeds.yaml"
SOIL_REPORT = PROJECT_ROOT / "pedology" / "analysis" / "soil_report.json"

# Output files worth keeping. LPJ-GUESS writes one per quantity, per rank.
#
# cflux.out carries the fire carbon flux. GLOBFIRM also exposes its inferred
# return time and burned fraction through firert.out; retain both so fire-driven
# mortality in cmass and dens has an occurrence diagnostic as well as a flux.
#
# mevap.out and mintercep.out are the only route to two of the three ways water
# leaves a gridcell. aaet.out Total is TRANSPIRATION -- commonoutput.cpp builds
# it from indiv.aaet -- while the model's own annual water record is
# patch.aaet + patch.aevap + patch.aintercep. Bare soil evaporation and canopy
# interception evaporation reach no annual table at all, so without these two
# monthly tables the water closure subtracts a third of the losses and calls
# the rest storage. notes/audits/closure-stocks-are-incomplete.md.
OUTPUTS = ("anpp.out", "lai.out", "fpc.out", "cmass.out", "aaet.out",
           "cpool.out", "dens.out", "agpp.out", "nsources.out", "cflux.out",
           "firert.out", "tot_runoff.out", "mevap.out", "mintercep.out",
           "nmass.out", "nuptake.out", "npool.out",
           "nflux.out", "ngases.out", "soil_npool.out", "soil_nflux.out")

NTRANSFORM_PROFILES = ("vesper", "stock-4.1.1")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def short(path: Path) -> str:
    return sha256(path)[:8]


def require_soil_driver_climate(soilmap: Path, driver: Path) -> dict:
    """BIO-19: prove the soil and ecological forcing use one climate state."""
    driver_report = driver.with_name(driver.stem + "_provenance.json")
    for path, role in ((SOIL_REPORT, "soil report"),
                       (driver_report, "LPJ driver provenance")):
        if not path.is_file():
            raise SystemExit(f"{path} is missing; cannot establish {role}")
    try:
        soil = json.loads(SOIL_REPORT.read_text(encoding="utf-8"))
        forcing = json.loads(driver_report.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise SystemExit(f"cannot read BIO-19 provenance: {exc}") from exc
    if soil.get("soilmap_sha256") != sha256(soilmap):
        raise SystemExit(
            f"{soilmap} changed after {SOIL_REPORT}; rebuild iteration soil")
    if soil.get("soil_state") != "baseline":
        raise SystemExit(
            f"{SOIL_REPORT} records soil_state={soil.get('soil_state')!r}; "
            "LPJ-GUESS must use soil weathered under the named baseline")
    soil_climate = soil.get("climatology_sha256")
    driver_climates = forcing.get("climatology_sha256")
    if driver_climates != [soil_climate]:
        raise SystemExit(
            "BIO-19 soil/driver climate mismatch: soil "
            f"{str(soil_climate)[:12]}, driver "
            f"{[str(value)[:12] for value in (driver_climates or [])]}")
    return {
        "soil_iteration": soil.get("iteration"),
        "soil_state": soil.get("soil_state"),
        "climatology_stage": soil.get("climatology_stage"),
        "climatology_sha256": soil_climate,
        "soil_report": str(SOIL_REPORT.relative_to(PROJECT_ROOT)),
        "soil_report_sha256": sha256(SOIL_REPORT),
        "driver_provenance": str(driver_report.relative_to(PROJECT_ROOT)),
        "driver_provenance_sha256": sha256(driver_report),
    }


def build_instruction(paths: dict, settings: dict) -> str:
    """The instruction file, with every path absolute. See the module docstring."""
    output_parameters = {
        # LPJ-GUESS's parameter predates the stock filename used by global.ins.
        "tot_runoff.out": "file_runoff",
    }
    return f"""! GENERATED by biosphere/scripts/run_lpj_guess.py. Do not edit.
!
! Absolute paths throughout: in parallel mode each rank chdirs into its own
! runN/ directory before reading anything.

import "{paths['pfts']}"

! The vendored source is the CNP fork, but Vesper's gridded soil phosphorus
! inputs and replacement productivity prediction are not ready. Keep the
! established C-N model until those land together; the later declaration wins.
ifplim 0
ifwalkernplim 0

! The CNP fork reads this parameter unconditionally in somdynam.cpp even when
! phosphorus limitation is disabled.  Empty selects the declared texture-soil
! weathering route; a non-empty value would activate the gridded route this
! world does not supply (biosphere/notes/cnp-fork-scoping.md).
param "file_pwtr" (str "")

! The volatile organic source. Written here rather than inherited, because an
! inherited zero cannot be told apart from nobody having decided. What requests
! it is biosphere/config/bvoc.yaml and what grants it is
! biosphere/scripts/bvoc_gate.py, which refuses until the preconditions it names
! are declared. plib takes the later declaration, so this overrides the imported
! PFT file whichever way it reads.
ifbvoc {settings['ifbvoc']}

! Acclimated respiration replaces each PFT's respcoeff with a function of the
! GROWTH temperature the tissue has adjusted to, and that growth temperature has
! a memory this project cannot derive but can bracket. The model carries it as
! Climate::tacc_air and Soil::tacc_root and refuses the option outright unless
! acclim_resp_tau declares the length, so nothing runs on an unstated one. The
! bracket, both its sources and the one-factor sensitivity registered over it
! are biosphere/config/respiration_acclimation.yaml; what enforces them is
! biosphere/scripts/acclimation_gate.py.
!
! THE BASELINE TAKES THE STANDARD RESPIRATION PATH, and that is a decision
! rather than an inheritance: it is declared as path.runs in
! biosphere/config/respiration_acclimation.yaml, the gate holds this line to
! it, and the CNP fork's own global_p.ins selects the same. The memory length
! is not what settles it. respiration_acclimated()
! takes no respcoeff argument at all, so switching paths replaces a coefficient
! Pft::init_cton_limits has normalised by the tissue C:N windows with Sprugel
! et al. (1996)'s two fixed reference rates, and that changes three things at
! once: the level of sapwood and fine-root maintenance respiration, what it
! depends on (leaf longevity on one path, growth temperature on the other), and
! whether it is invariant under a tissue C:N window rescaling, which is what
! kept respiration() unmoved by WORLD-XMS4. Acclimation is real and is not what
! is refused; a switch that carries two unsourced changes with it is. The
! arithmetic and what would reopen it are in that same declaration. WORLD-UHFH.
acclimated_respiration 0
! The simulated wetlands, their peat and their methane. Four switches, written
! rather than inherited for the same reason: an inherited zero cannot be told
! apart from nobody having decided. What requests them is
! biosphere/config/wetlands.yaml and what grants them is
! biosphere/scripts/wetland_gate.py, which refuses on two grounds -- a
! precondition still undeclared, and a hydrology defect still present in the
! vendored source, where the low-latitude wetland path adds water the simulated
! world never received. plib takes the later declaration, so these override the
! imported files whichever way they read.
run_peatland {settings['run_peatland']}
ifmethane {settings['ifmethane']}
ifsaturatewetlands {settings['ifsaturatewetlands']}
wetland_runon {settings['wetland_runon']}

title "{settings['title']}"
nyear {settings['nyear']}
vesper_root_seed {settings['root_seed']}

file_driver "{paths['driver']}"
file_soilmap "{paths['soilmap']}"

! BLAZE needs a SimFIRE input file built from Earth observations and GWGEN needs
! sub-daily statistics this world does not have. The shipped demo makes the same
! two substitutions for the same reason.
firemodel "GLOBFIRM"
weathergenerator "INTERP"

npatch {settings['npatch']}
nfix_a {settings['nfix_a']}
nfix_b {settings['nfix_b']}

! Without this, soil carbon never reaches the hydrology and the soil-biosphere
! loop closes only inside pedology/. Stock LPJ-GUESS derives water-holding
! capacity from texture alone; iforganicsoilproperties makes it blend the
! mineral and organic water retention by the organic fraction, which is the
! path pedology's cpool.out feedback has to travel.
!
! It requires a SoilC column in the soil map, which build_soil.py writes, and it
! requires iftwolayersoil 0, which global.ins already sets.
iforganicsoilproperties 1

outputdirectory "./"
""" + "".join(
        f'{output_parameters.get(name, "file_" + name.split(".")[0])} "{name}"\n'
        for name in settings["outputs"])


def merge_outputs(run_dir: Path, ranks: int,
                  outputs: tuple[str, ...] = OUTPUTS) -> dict[str, int]:
    """Concatenate per-rank output, keeping one header.

    Rank directories each hold a complete set of files covering their own cells,
    so merging is a concatenation. The row count is returned per file so a rank
    that died quietly shows up as a short file rather than as a silently smaller
    planet.
    """
    counts = {}
    for name in outputs:
        pieces = [run_dir / f"run{rank}" / name for rank in range(1, ranks + 1)]
        pieces = [p for p in pieces if p.is_file()]
        if not pieces:
            continue
        target = run_dir / name
        with target.open("w") as out:
            for index, piece in enumerate(pieces):
                lines = piece.read_text().splitlines(keepends=True)
                if not lines:
                    continue
                out.writelines(lines if index == 0 else lines[1:])
        counts[name] = sum(1 for _ in target.read_text().splitlines()) - 1
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nyear", type=int, default=50,
                        help="simulation years after spin-up")
    parser.add_argument("--npatch", type=int, default=5,
                        help="replicate patches per gridcell")
    parser.add_argument("--ranks", type=int, default=16,
                        help="MPI ranks. 16 physical cores on this machine.")
    parser.add_argument("--nfix-a", type=float, default=0.234,
                        help="Cleveland fixation slope. LPJ-GUESS brackets this "
                             "0.102 to 0.367 and it is worth 18%% on NPP; sweep it "
                             "rather than quoting the central value alone.")
    parser.add_argument("--nfix-b", type=float, default=-0.172)
    parser.add_argument("--driver", type=Path, default=GENERATED / "vesper_driver.bin")
    parser.add_argument("--soilmap", type=Path,
                        default=None)
    parser.add_argument("--pfts", type=Path, default=GENERATED / "vesper_pfts.ins")
    parser.add_argument("--binary", type=Path, default=GUESS_BINARY,
                        help="LPJ-GUESS binary to execute; its hash is recorded")
    parser.add_argument("--ntransform-profile", choices=NTRANSFORM_PROFILES,
                        default="vesper",
                        help="soil-N operator/instruction profile for a matched arm")
    parser.add_argument("--label", default=None,
                        help="human tag recorded in the run manifest")
    parser.add_argument("--root-seed", type=int, default=None,
                        help="root of all stable ecological random substreams; "
                             "defaults to stochastic_seeds.yaml")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    import sys as _sys
    _sys.path.insert(0, str(PROJECT_ROOT / "lib"))
    if args.soilmap is None:
        import builds as _b
        args.soilmap = _b.soilmap()

    for path in (args.driver, args.soilmap, args.pfts, args.binary):
        if not Path(path).is_file():
            raise SystemExit(f"{path} is missing")

    # THE BINARY MUST CONTAIN THE MODEL IN THIS TREE, and a run is where that is
    # decided. LPJ-GUESS reads none of its source at run time, so a binary built
    # before an edit integrates the code it was built from and reports nothing:
    # the 2026-08-30 case was caught only because the missing parameter reached
    # the model through an instruction file, so the PARSER refused it. A change
    # confined to C++ has no parser in front of it. This refuses instead.
    # world-w62x.
    #
    # It is checked HERE and not only in `scripts/check_consistency.py` because
    # the gate is not run before every run, and the binary can move between the
    # two -- ExoPlaSim split the same question the same way after world-anl,
    # where `--verify` truthfully reported every installed binary current while
    # a probe ran a copy taken from a run directory. It costs a JSON read and a
    # sha over the source set, so it can afford to be here.
    provenance_path = args.binary.with_suffix(".provenance.json")
    binary_provenance = None
    if args.ntransform_profile == "stock-4.1.1":
        # The comparison arm is a build of a COPY of the tree with one
        # compilation unit substituted, so its record is a different contract
        # and `prepare_stock_ntransform_arm.py` is what writes it.
        if not provenance_path.is_file():
            raise SystemExit(
                f"{provenance_path} is missing; prepare the stock arm with "
                "prepare_stock_ntransform_arm.py")
        binary_provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        if binary_provenance.get("profile") != args.ntransform_profile:
            raise SystemExit(f"{provenance_path} is not a stock-4.1.1 arm")
        if binary_provenance.get("binary_sha256") != sha256(args.binary):
            raise SystemExit(f"{args.binary} changed after {provenance_path}")
        # The arm was built from a snapshot of the vendored tree. If that tree
        # has moved since, the arm no longer differs from the active model in
        # `ntransform.cpp` alone, and the comparison it exists for is void.
        from prepare_stock_ntransform_arm import source_fingerprint
        live = source_fingerprint(GUESS_SOURCE)
        if binary_provenance.get("base_source_fingerprint") != live:
            raise SystemExit(
                f"{args.binary} was built from a snapshot of vendor/lpj-guess "
                f"that no longer matches the tree, so the stock arm and the "
                f"active arm differ by more than modules/ntransform.cpp. "
                f"Re-prepare it with prepare_stock_ntransform_arm.py.")
    else:
        from build_lpj_guess import verify as verify_binary
        problems = verify_binary(args.binary)
        if problems:
            raise SystemExit("\n".join(problems) + "\n\nRefusing to run: the "
                             "numbers would be attributed to source this "
                             "binary does not contain.")
        # The full source map lives beside the binary, which is untracked and
        # gets overwritten by the next build. The run manifest is tracked and
        # there is one per run, so it carries the digest and not the 124 shas.
        record = json.loads(provenance_path.read_text(encoding="utf-8"))
        binary_provenance = {k: v for k, v in record.items() if k != "sources"}
        binary_provenance["source_count"] = len(record.get("sources") or {})

    config = yaml.safe_load(CONFIG.read_text())
    from stochastic_seeds import read_declaration
    from lpj_output import driver_header
    seed_contract = read_declaration(SEED_CONFIG)
    forcing_header = driver_header(Path(args.driver))
    root_seed = (seed_contract["root_seed"] if args.root_seed is None
                 else args.root_seed)
    if not 1 <= root_seed <= 2147483646:
        raise SystemExit("--root-seed must be in [1, 2147483646]")

    # The driver is built from a climatology and lives in a FLAT directory, so
    # nothing about its path says which build it belongs to. The soil map is
    # per build and resolves strictly; pairing the two without checking is how
    # one terrain's forcing gets run against another's soil.
    from provenance import require_build
    require_build(Path(args.driver), "LPJ-GUESS driver", config,
                  allow_unstamped=False)
    soil_climate_contract = require_soil_driver_climate(
        Path(args.soilmap), Path(args.driver))

    # The run id pins what the numbers depend on: the forcing, the soil and the
    # compiled planetary constants. Two runs sharing a name share those.
    header = GENERATED / "vesper.h"
    # A UUID, for the reason ExoPlaSim run ids are UUIDs: a derived name
    # separates runs only along the dimensions it happens to encode, and the
    # encoded set is just a list of everything someone has thought of so far.
    # This name carried years, patches, ranks, the fixation slope and a digest
    # of three inputs -- and nothing about CO2, ndep, or the spectrum, so two
    # runs differing only in those computed the same name and would have shared
    # a directory. What a run WAS lives in run_manifest.json, which records all
    # of it including the input hashes.
    run_id = f"lpj_{uuid.uuid4().hex}"

    run_dir = RUNS / run_id
    if run_dir.exists():
        raise SystemExit(f"{run_dir} exists; remove it before re-preparing")

    paths = {
        "driver": Path(args.driver).resolve(),
        "soilmap": Path(args.soilmap).resolve(),
        "pfts": (run_dir / "vesper_pfts.ins").resolve(),
    }
    # The volatile organic source, asked rather than assumed. `require` exits
    # with every unmet precondition named if biosphere/config/bvoc.yaml requests
    # activation; with no request it grants nothing and the run is a correct run
    # with the source off. Retaining the compound-resolved tables is part of
    # activation and not a separate decision: a run that emits volatile carbon
    # and discards the speciation is the collapse the contract exists to stop.
    import bvoc_gate
    bvoc_active = bvoc_gate.require(planet=config)
    outputs = OUTPUTS + (bvoc_gate.ACTIVATED_OUTPUTS if bvoc_active else ())

    # The wetlands, asked on the same terms. `require` exits with every unmet
    # precondition named if biosphere/config/wetlands.yaml requests activation;
    # with no request it grants nothing and the run is a correct run with peat
    # and methane off. Retaining the area, water-table, pathway and residual
    # tables is part of activation and not a separate decision: a run that emits
    # methane and discards the pathway partition is exactly the aggregate the
    # contract exists to stop being accepted.
    import wetland_gate
    wetland_active = wetland_gate.require(planet=config)
    if wetland_active:
        outputs = outputs + tuple(
            wetland_gate.read_declaration()["acceptance"]["retained_outputs"])

    settings = {
        "title": run_id, "nyear": args.nyear, "npatch": args.npatch,
        "root_seed": root_seed,
        "nfix_a": args.nfix_a, "nfix_b": args.nfix_b,
        "ifbvoc": 1 if bvoc_active else 0,
        "outputs": outputs,
        **wetland_gate.switches(wetland_active),
    }

    if args.dry_run:
        print(build_instruction(paths, settings))
        return

    run_dir.mkdir(parents=True)
    # The PFT file and everything it imports have to sit together, because plib
    # resolves imports relative to the importing file's directory.
    shutil.copyfile(args.pfts, run_dir / "vesper_pfts.ins")
    for extra in (GUESS_SOURCE / "data" / "ins").glob("*.ins"):
        if not (run_dir / extra.name).exists():
            shutil.copyfile(extra, run_dir / extra.name)

    soiln_instruction = run_dir / "global_soiln.ins"
    if args.ntransform_profile == "stock-4.1.1":
        soiln_text = soiln_instruction.read_text(encoding="utf-8")
        old = "f_nitri_gas_max \t0.022"
        new = "f_nitri_gas_max \t0.25"
        if soiln_text.count(old) != 1:
            raise SystemExit(
                f"stock profile expected one live {old!r} in {soiln_instruction}")
        soiln_instruction.write_text(soiln_text.replace(old, new), encoding="utf-8")

    instruction = run_dir / "run.ins"
    instruction.write_text(build_instruction(paths, settings))
    for rank in range(1, args.ranks + 1):
        (run_dir / f"run{rank}").mkdir(exist_ok=True)

    command = ["mpirun", "-np", str(args.ranks), "--bind-to", "core",
               str(args.binary.resolve()), "-parallel", "-input", "vesper",
               str(instruction.resolve())]
    print(" ".join(command))
    started = time.time()
    result = subprocess.run(command, cwd=run_dir, capture_output=True, text=True)
    elapsed = time.time() - started
    (run_dir / "mpirun.log").write_text(result.stdout + result.stderr)
    if result.returncode != 0:
        raise SystemExit(
            f"LPJ-GUESS exited {result.returncode} after {elapsed:.0f} s. "
            f"See {run_dir / 'mpirun.log'} and {run_dir}/run*/guess.log")

    counts = merge_outputs(run_dir, args.ranks, outputs)
    cells = 0
    if (run_dir / "anpp.out").is_file():
        rows = (run_dir / "anpp.out").read_text().splitlines()[1:]
        cells = len({(r.split()[0], r.split()[1]) for r in rows if r.split()})

    manifest = {
        "run_id": run_id,
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        # What the run id used to spell out. Ids are UUIDs now, so this block is
        # the only place a run's physical description survives.
        "physical": {
            "nyear": args.nyear,
            "npatch": args.npatch,
            "root_seed": root_seed,
            "ranks": args.ranks,
            "nfix_a": args.nfix_a,
            "nfix_b": args.nfix_b,
            "label": args.label,
            "ntransform_profile": args.ntransform_profile,
            "ifbvoc": settings["ifbvoc"],
            "run_peatland": settings["run_peatland"],
            "ifmethane": settings["ifmethane"],
        },
        "wall_seconds": round(elapsed, 1),
        "ranks": args.ranks,
        "settings": settings,
        "cells_simulated": cells,
        "output_rows": counts,
        "forcing": forcing_header,
        "stochastic_randomness": {
            "contract_version": seed_contract["contract_version"],
            "algorithm": seed_contract["algorithm"],
            "root_seed": root_seed,
            "keys": seed_contract["keys"],
            "processes": seed_contract["processes"],
            "declaration": str(SEED_CONFIG.relative_to(PROJECT_ROOT)),
            "declaration_sha256": sha256(SEED_CONFIG),
            "invariance": "independent of MPI rank and grid-cell traversal order",
        },
        "inputs": {
            "driver": {"path": str(paths["driver"]), "sha256": sha256(args.driver)},
            "soilmap": {"path": str(paths["soilmap"]), "sha256": sha256(args.soilmap)},
            "pfts": {"path": str(args.pfts), "sha256": sha256(args.pfts)},
            "vesper_h": {"path": str(header), "sha256": sha256(header)},
            "binary": {"path": str(args.binary.resolve()),
                       "sha256": sha256(args.binary)},
            "global_soiln": {"path": str(soiln_instruction),
                             "sha256": sha256(soiln_instruction)},
        },
        "ntransform_comparison": {
            "profile": args.ntransform_profile,
            "binary_provenance": binary_provenance,
        },
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "source_build": config.get("source_build"),
        "soil_driver_climate_contract": soil_climate_contract,
        "stellar_spectrum": config.get("radiation", {}).get("stellar_spectrum"),
        "orbital_year_earth_days": orbit.orbital_year_days(config),
        "model_year_days": orbit.model_year_days(config),
        "earth_year_days": orbit.EARTH_SIDEREAL_YEAR_DAYS,
        "earth_years_per_simulation_year": orbit.earth_years_per_model_year(config),
        "annual_flux_per_orbit_to_per_earth_year": (
            1.0 / orbit.earth_years_per_model_year(config)),
        "reporting_interval": (
            "one simulation year, which is one orbit of "
            f"{orbit.model_year_days(config)} absolute days"),
        "unit_note": (
            "Every annual column in the .out files is a sum over one simulation "
            "year, which is one orbit and not one Earth year. Multiply such a "
            "flux by annual_flux_per_orbit_to_per_earth_year before comparing "
            "with anything quoted per Earth year, including "
            "notes/productivity-prediction.md. The name says which direction "
            "the factor goes, because the two directions differ by four in the "
            "answer and look identical in a script. The registry of which "
            "quantity is in which unit is "
            "biosphere/notes/time-base-unit-contract.md."),
        "software": {
            "platform": platform.platform(),
            "lpj_guess": (
                "LPJ-GUESS-CNP v1.0 (b368b893c4324840b43c56866e901c6916858afb) "
                "+ in-tree Vesper port; phosphorus limitation disabled"
            ),
        },
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    (run_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    analysis_dir = ANALYSIS / run_id
    analysis_dir.mkdir(parents=True, exist_ok=True)
    (analysis_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n")

    # BIO-14 makes completion an assessed artifact. No coupled consumer is
    # allowed to infer success from mpirun's exit status or a few present files.
    import assess_lpj_run
    try:
        acceptance = assess_lpj_run.assess(run_dir)
    except Exception as exc:
        refusal = assess_lpj_run.write_failure(run_dir, exc)
        raise SystemExit(
            f"LPJ-GUESS output refused: {refusal['refusal']}. "
            f"See {run_dir / 'acceptance.json'}") from exc

    print(f"\n{run_id}")
    print(f"  {cells} cells, {elapsed:.0f} s wall on {args.ranks} ranks")
    print(f"  accepted {acceptance['coverage']['years'][0]}-"
          f"{acceptance['coverage']['years'][-1]} ({len(acceptance['coverage']['years'])} years)")
    print(f"  output   {run_dir}")
    print(f"  manifest {analysis_dir / 'run_manifest.json'}")


if __name__ == "__main__":
    main()
