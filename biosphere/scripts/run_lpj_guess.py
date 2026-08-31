"""Run LPJ-GUESS over Vesper, reproducibly and with its provenance recorded.

Every run before this one was hand-assembled in a scratch directory with a
hand-written instruction file, which is fine for a smoke test and not fine for
anything whose numbers get quoted. This is the equivalent of
`exoplasim/scripts/run_exoplasim.py`: it generates the instruction file, sets up
the rank directories, runs, merges the per-rank output and writes a manifest
pinning every input by hash.

    python biosphere/scripts/run_lpj_guess.py --nyear 50
    python biosphere/scripts/run_lpj_guess.py --ranks 16 --label carved-k25v
    python biosphere/scripts/run_lpj_guess.py --nyear 8600 --save-state
    python biosphere/scripts/run_lpj_guess.py --nyear 9853 --continue-from lpj_...

Bulk output lands in `runs/<run_id>/`, which is not tracked, in the same way
`exoplasim/runs/` is not. The manifest and a summary go to
`analysis/<run_id>/`, which is.

**The spin-up is bought once.** The ecological spin-up in front of a retained
record is several times the record itself, so a run refused by the acceptance
contract used to be answerable only by a second run from bare ground that
re-integrated the whole of it. `--save-state` leaves the simulated state of the
last year on disk and `--continue-from` resumes from it, which makes more
retained record, another seed or another patch count cost the record alone. The
continuation is a NEW run with a new id whose manifest names its parent, and
`--nyear` is CUMULATIVE: it is the simulated year the run ends at, not the years
it adds. What may be continued is guarded rather than trusted -- see
`continuation_refusals` -- and that the resumed run reproduces the run it
continues is `biosphere/scripts/verify_lpj_restart_continuity.py`.

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
# maet.out, mevap.out and mintercep.out are the only route to the three
# evaporative ways water leaves a gridcell. aaet.out Total is not one of them:
# commonoutput.cpp builds it from indiv.aaet over the individuals still in
# patch.vegetation, and framework.cpp runs mortality, establishment and fire
# disturbance before outannual, so it reports the transpiration of the
# survivors. Bare soil evaporation and canopy interception reach no annual
# table at all. notes/audits/closure-stocks-are-incomplete.md.
OUTPUTS = ("anpp.out", "lai.out", "fpc.out", "cmass.out", "aaet.out",
           "cpool.out", "dens.out", "agpp.out", "nsources.out", "cflux.out",
           "firert.out", "tot_runoff.out",
           "maet.out", "mevap.out", "mintercep.out",
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


def soiln_instruction_text(profile: str) -> str:
    """`global_soiln.ins` as the run will read it, for the profile asked for.

    The stock comparison arm differs from the active model in the vendored
    `ntransform.cpp` and in one instruction-file rate, so the file a run reads
    is a function of the profile and not of the vendored tree alone. Built here
    rather than patched in place, so `--continue-from` can hash what the parent
    read without first laying out a run directory.
    """
    text = (GUESS_SOURCE / "data" / "ins" / "global_soiln.ins").read_text(
        encoding="utf-8")
    if profile != "stock-4.1.1":
        return text
    old, new = "f_nitri_gas_max \t0.022", "f_nitri_gas_max \t0.25"
    if text.count(old) != 1:
        raise SystemExit(
            f"stock profile expected one live {old!r} in "
            f"{GUESS_SOURCE / 'data' / 'ins' / 'global_soiln.ins'}")
    return text.replace(old, new)


# What two runs have to share for the second to be a CONTINUATION of the first
# rather than a chimera. A state file is a simulated state, not a result: read
# it under a different forcing, a different soil, a different binary or a
# different ecological parameter and the run integrates from a state its own
# inputs never produced, while every number out of it is attributed to a
# spin-up that did not happen. That is the failure a saved state creates if it
# is not guarded, and it is worse than having no saved state at all, because
# nothing downstream can see it.
#
# Ranks are deliberately NOT here. The serializer partitions state files by the
# writing rank, but `PartitionedMapDeserializer` searches every file in the
# directory for a cell's coordinates and the stochastic substreams are derived
# from the cell rather than the rank (lib/stochastic_seeds.py), so a
# continuation at a different rank count reads the same state. It is recorded
# rather than refused.
CONTINUATION_INPUTS = ("driver", "soilmap", "pfts", "binary", "vesper_h",
                       "global_soiln")
CONTINUATION_PHYSICAL = ("npatch", "root_seed", "nfix_a", "nfix_b",
                         "ntransform_profile", "ifbvoc", "run_peatland",
                         "ifmethane")


def continuation_refusals(parent: dict, inputs: dict, physical: dict,
                          nyear: int) -> list[str]:
    """Every way this run is not a continuation of the run it names as parent."""
    refusals: list[str] = []
    parent_inputs = parent.get("inputs") or {}
    for name in CONTINUATION_INPUTS:
        was = (parent_inputs.get(name) or {}).get("sha256")
        now = (inputs.get(name) or {}).get("sha256")
        if was is None:
            refusals.append(
                f"the parent manifest pins no {name}, so this run cannot show "
                "it is continuing the same experiment")
        elif was != now:
            refusals.append(
                f"{name} differs: the parent integrated {was[:12]}, this run "
                f"would read {now[:12]}")

    parent_physical = parent.get("physical") or {}
    for name in CONTINUATION_PHYSICAL:
        if name not in parent_physical:
            refusals.append(f"the parent manifest records no {name}")
        elif parent_physical[name] != physical[name]:
            refusals.append(
                f"{name} differs: the parent ran {parent_physical[name]!r}, "
                f"this run would run {physical[name]!r}")

    parent_nyear = parent_physical.get("nyear")
    if not isinstance(parent_nyear, int):
        refusals.append("the parent manifest records no nyear, so the "
                        "simulated year the state covers is not known")
    elif nyear <= parent_nyear:
        refusals.append(
            f"--nyear is CUMULATIVE: the parent already reached year "
            f"{parent_nyear}, so a continuation needs --nyear greater than "
            f"{parent_nyear} and was given {nyear}")
    return refusals


def serialization_block(state: dict | None) -> str:
    """The state-file declarations, or nothing when the run neither saves nor reads one.

    THE ONE PLACE THESE ARE SPELLED. `verify_lpj_restart_continuity.py` builds
    three different save/restart arrangements and this runner builds two more,
    and a state file written by one arrangement is read by another; two copies
    of the block would be two chances for the fixture's meaning of `state_year`
    to drift from the runner's while both keep parsing.

    `state_year`/`state_day` name the last simulated day the state READ covers;
    `save_year`/`save_day` name the last day the state WRITTEN covers. A day of
    -1 is the year boundary, which is the whole of the preceding year and the
    only save point this model had before WORLD-FUJ4. `state_path` is absolute
    for the reason every path in a generated instruction file is: each rank
    chdirs into its own directory before reading anything.
    """
    if not state:
        return ""
    return f"""
! The state file. What a run saves, another run continues from: the spin-up in
! front of a retained record is then bought once rather than once per refusal.
! Written by biosphere/scripts/run_lpj_guess.py; the two instants and the -1
! convention are in its serialization_block.
restart {1 if state['restart'] else 0}
save_state {1 if state['save_state'] else 0}
state_year {state['state_year']}
state_day {state['state_day']}
save_year {state['save_year']}
save_day {state['save_day']}
state_path "{state['state_path']}"
save_path "{state['save_path']}"
"""


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
        for name in settings["outputs"]) + serialization_block(settings.get("state"))


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
    parser.add_argument("--save-state", action="store_true",
                        help="write a state file covering the last simulated "
                             "year, so a later run can continue from it "
                             "instead of re-buying the spin-up")
    parser.add_argument("--continue-from", default=None, metavar="RUN_ID",
                        help="resume from that run's saved state. --nyear is "
                             "CUMULATIVE: it is the year this run ends at, not "
                             "the years it adds. The continuation is a NEW run "
                             "with a new id whose manifest names its parent.")
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

    physical = {
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
    }

    soiln_text = soiln_instruction_text(args.ntransform_profile)
    inputs = {
        "driver": {"path": str(paths["driver"]), "sha256": sha256(args.driver)},
        "soilmap": {"path": str(paths["soilmap"]), "sha256": sha256(args.soilmap)},
        "pfts": {"path": str(args.pfts), "sha256": sha256(args.pfts)},
        "vesper_h": {"path": str(header), "sha256": sha256(header)},
        "binary": {"path": str(args.binary.resolve()),
                   "sha256": sha256(args.binary)},
        "global_soiln": {"path": str(run_dir / "global_soiln.ins"),
                         "sha256": hashlib.sha256(
                             soiln_text.encode("utf-8")).hexdigest()},
    }

    # THE STATE FILE. Two arrangements, and the difference between them is which
    # of the two instants the caller names.
    #
    # `--save-state` writes the state covering the END of the last simulated
    # year, which is where every annual accumulator has flushed and where the
    # save point sat before WORLD-FUJ4 made an arbitrary day expressible.
    #
    # `--continue-from` reads that state and resumes on the day after it. The
    # continued run is a NEW run with a new id (rule 6) whose manifest names its
    # parent, and `--nyear` is the year it ends at rather than the years it
    # adds, because `vesperinput.cpp:getclimate` stops on `date.year` reaching
    # `nyear_spinup + nyear` and counts from year zero whether the state was
    # read or integrated.
    #
    # A run that both reads and writes needs the two directories distinct: the
    # serializer truncates what it opens and it opens before the deserializer
    # reads (`framework/parameters.cpp` refuses the collision).
    continuation = None
    state = None
    if args.continue_from:
        parent_dir = RUNS / args.continue_from
        parent_manifest_path = parent_dir / "run_manifest.json"
        parent_state = parent_dir / "state"
        if not parent_manifest_path.is_file():
            raise SystemExit(
                f"{parent_manifest_path} is absent, so there is no record of "
                f"what {args.continue_from} was and nothing to continue.")
        if not (parent_state / "meta.bin").is_file():
            raise SystemExit(
                f"{parent_state} holds no state file, so {args.continue_from} "
                "cannot be continued. A run saves one only when it is asked "
                "with --save-state.")
        parent = json.loads(parent_manifest_path.read_text(encoding="utf-8"))
        refusals = continuation_refusals(parent, inputs, physical, args.nyear)
        if refusals:
            raise SystemExit(
                f"Refusing to continue {args.continue_from}:\n  "
                + "\n  ".join(refusals)
                + "\n\nA state file is a simulated state, not a result. Read "
                "under different inputs it integrates from a state those "
                "inputs never produced, and every number out of it would be "
                "attributed to a spin-up that did not happen.")
        # THE STATE FILES HAVE TO BE THE ONES THE PARENT WROTE. `runs/` is
        # untracked and nothing stops a directory being cleaned, half-copied
        # from elsewhere or rewritten by a later save, and the archive is an
        # untagged stream of raw object bytes: a truncated or substituted file
        # is read as a simulated state without complaint. The parent recorded
        # what it wrote, so this is an identity with a right answer rather than
        # a plausibility check.
        state_now = {p.name: sha256(p) for p in sorted(parent_state.iterdir())
                     if p.is_file()}
        state_then = (parent.get("saved_state") or {}).get("sha256")
        if state_then is None:
            raise SystemExit(
                f"{parent_manifest_path} records no hashes for the state files "
                f"it wrote, so there is nothing to show {parent_state} still "
                "holds them. That manifest predates --save-state; re-run the "
                "parent rather than continuing a state file of unknown "
                "provenance.")
        if state_then != state_now:
            gone = sorted(set(state_then) - set(state_now))
            extra = sorted(set(state_now) - set(state_then))
            moved = sorted(name for name in set(state_then) & set(state_now)
                           if state_then[name] != state_now[name])
            raise SystemExit(
                f"{parent_state} no longer holds the state files "
                f"{args.continue_from} wrote"
                + (f"; absent: {', '.join(gone)}" if gone else "")
                + (f"; changed: {', '.join(moved)}" if moved else "")
                + (f"; unexpected: {', '.join(extra)}" if extra else "")
                + ". The archive is an untagged byte stream, so a substituted "
                "file would be read as a simulated state without complaint.")

        parent_nyear = parent["physical"]["nyear"]
        continuation = {
            "parent_run_id": parent.get("run_id", args.continue_from),
            "parent_manifest_sha256": sha256(parent_manifest_path),
            "parent_state_dir": str(parent_state.resolve()),
            "parent_state_sha256": state_now,
            "parent_nyear": parent_nyear,
            "parent_ranks": (parent.get("physical") or {}).get("ranks"),
            "resumed_at_year": parent_nyear,
            "years_simulated_here": args.nyear - parent_nyear,
            "chain": (parent.get("continuation") or {}).get("chain", [])
                     + [parent.get("run_id", args.continue_from)],
        }
        state = {
            "restart": True, "save_state": bool(args.save_state),
            "state_year": parent_nyear, "state_day": -1,
            "save_year": args.nyear, "save_day": -1,
            "state_path": str(parent_state.resolve()),
            "save_path": str((run_dir / "state").resolve()),
        }
    elif args.save_state:
        state = {
            "restart": False, "save_state": True,
            "state_year": args.nyear, "state_day": -1,
            "save_year": args.nyear, "save_day": -1,
            "state_path": str((run_dir / "state").resolve()),
            "save_path": str((run_dir / "state").resolve()),
        }
    settings["state"] = state

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
    soiln_instruction.write_text(soiln_text, encoding="utf-8")

    # The serializer opens its files in its constructor and fails if the
    # directory is absent, so a saving run makes it before the model starts.
    if state and state["save_state"]:
        Path(state["save_path"]).mkdir(parents=True, exist_ok=True)

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

    # A run asked to save a state has to have saved one. The serializer writes
    # at a named simulated instant and the model exits zero whether it reached
    # that instant or not, so an empty state directory is how "the save point
    # was never reached" arrives -- and it arrives as a run that looks complete
    # and cannot be continued, discovered only when someone tries.
    if state and state["save_state"]:
        written = [p for p in Path(state["save_path"]).iterdir() if p.is_file()]
        if not any(p.name == "meta.bin" for p in written):
            raise SystemExit(
                f"the run was asked to save a state covering year "
                f"{state['save_year'] - 1} and {state['save_path']} holds "
                f"{len(written)} files with no meta.bin. The save point was "
                "not reached, so nothing can continue this run.")

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
        "physical": physical,
        # Which state this run started from and which one it left behind. A
        # continuation's `nyear` is the year it ended at, so `resumed_at_year`
        # and `years_simulated_here` are what say how much of the record this
        # run actually integrated; `chain` is the whole lineage back to bare
        # ground. Null here means bare ground.
        "continuation": continuation,
        "saved_state": ({
            "dir": str(Path(state["save_path"])),
            "covers_year": state["save_year"] - 1,
            "covers_day": "year boundary",
            "sha256": {p.name: sha256(p)
                       for p in sorted(Path(state["save_path"]).iterdir())
                       if p.is_file()},
        } if state and state["save_state"] else None),
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
        "inputs": inputs,
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
