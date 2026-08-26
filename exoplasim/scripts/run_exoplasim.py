#!/usr/bin/env python3
"""Prepare, validate, and optionally run a reproducible ExoPlaSim experiment."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import shutil
import subprocess
import uuid

import exoplasim as exo
from netCDF4 import Dataset
import numpy as np
import yaml

from _paths import CONFIG, INPUTS, MODEL_SRC, PROJECT_ROOT, RUNS
import build_model
import rebuild_binaries
import reset_restart_accumulators
import rungs


EARTH_STANDARD_GRAVITY = 9.80665
EARTH_SIDEREAL_YEAR_DAYS = 365.2568983
SOLAR_EFFECTIVE_TEMPERATURE_K = 5772.0

# Climate and balance fields retained in compact annual products.
#
# CODE 151, psl, IS DELIBERATELY ABSENT from both lists. world-ld1. It is a
# reduction of surface pressure to sea level, and a reduction is only as good
# as the lapse rate it extrapolates with: pyburn does it at RLAPSE = 0.0065,
# Earth's International Standard Atmosphere rate, against 7.8 K/km measured
# from this model's own sigma profile, and then applies ECMWF's 255 K
# cold-surface guard, which is calibrated on Earth's surface temperature
# distribution and on this world fires over ordinary land for much of the
# orbit. An honest rate has to be MEASURED from a climatology, and
# config/planet.yaml declares `baseline_climatology: null`, so on this build
# no honest rate exists. Nothing in this project reads psl -- checked across
# every component -- so the choice is between not writing it and writing a
# CF-named field that looks authoritative and is not.
#
# Restore 151 when lib/lapse.environmental_lapse_k_per_km has a climatology to
# measure, and pass the rate and the guard through rather than restoring the
# defaults; the exact threading is in world-ld1.
#
# 231 is the land column's drainage, WORLD-P9QQ, and it is here unconditionally
# even though it is identically zero under the default impermeable lower
# boundary. A code the model writes and no product carries is the state
# world-j0az found the air extrema in; the cost of carrying it is one 2D field
# per output bin, and the ledger's `drainage` crossing becomes readable the
# moment a draining lower boundary is selected rather than one change later.
#
# 201 and 202 are the extrema of `dtsa`, the near-surface AIR temperature, and
# outmod has written them unconditionally all along. They were the state the
# comment above describes: a code the model writes and no product carries, and
# the reason `climate.dtr` -- a diurnal AIR temperature range -- had no source.
# 320 and 321 are extrema of `dt(:,NLEP)`, the SURFACE temperature; both pairs
# travel, and pyburn names them apart as `tasmax`/`tasmin` and `maxt`/`mint`.
# world-j0az; biosphere/notes/ecological-forcing-field-contract.md.
#
# 185 and 186 are the surface downward shortwave in the model's two bands,
# WORLD-3QFZ. Their sum is the INCIDENT surface shortwave, which 176 is not:
# 176 is the net flux. Both readers are downstream and both exist today -- the
# biosphere's PAR fraction is a surface quantity anchored on a top-of-atmosphere
# band ratio for want of these, and its driver wants incident where it currently
# passes net.
#
# 163, 168, 171 AND 238 ARE GONE, and they were never produced. `pyburn`
# answers a code either from the raw stream or from a derivation branch; these
# four have neither, so `dataset` marked each derived, matched no branch and
# dropped it without a word. None of the four is in any climatology.
# What they asked for is elsewhere or nowhere: 163 `tcc` is total cloud cover,
# which 164 `clt` already carries from `dcc(:,NLEP)`; 171 `dsw` is a second soil
# wetness the model does not write, its soil water being 140 `mrso` from
# `dwatc`; 238 `tsn` is a snow temperature the model has no prognostic for; and
# 168 `td2m` is a dew point, which the forcing contract does not want because it
# takes near-surface humidity as specific humidity, the linear quantity the
# model integrates. Deriving any of them is a new field with a saturation branch
# to declare, not a restoration -- world-dy7a. `smoke_test.py` now refuses a
# code in either list that nothing writes or derives.
REGULAR_CODES = [
    50, 51, 52, 53, 54, 110, 129, 130, 131, 132, 133, 134, 135, 139,
    140, 141, 142, 143, 144, 146, 147, 157, 160, 164, 167,
    170, 172, 174, 175, 176, 177, 178, 179, 180,
    181, 182, 184, 185, 186, 201, 202, 203, 204, 205, 207, 208, 209, 210,
    211, 218, 221,
    230, 231, 232, 259, 260, 261, 262, 263, 264, 267, 318, 320, 321,
]

# Instantaneous seasonal output omits redundant flux details and expensive
# three-dimensional humidity and vertical-motion fields.
#
# NO EXTREMA IN THIS LIST, and that is the shape of the stream rather than an
# omission. 320/321 and 201/202 are running extrema `outaccu` extends every
# timestep and only `outreset` clears, and `outreset` runs after `outgp` on the
# REGULAR cadence: in a snapshot record they would cover (nstep mod nafter)
# timesteps, a window that varies from record to record and is empty on the
# record after a regular write. Every other field here is an instantaneous
# sample, so a consumer taking a variance or a distribution over such a field
# would be reading a sawtooth of window lengths. 201/202 were kept out on this
# argument under world-j0az; world-adxx took 320/321 out to match and
# `outmod.f90:snapshotgp` no longer writes them at all. The instantaneous
# surface temperature this stream carries is 139; the extrema that mean
# something are the regular stream's.
SNAPSHOT_CODES = [
    50, 51, 52, 53, 54, 129, 130, 131, 132, 134, 139, 140, 141, 142,
    143, 144, 160, 164, 167, 170, 172,
    175, 180, 181, 182, 210, 211, 218, 230, 232, 259, 260, 261, 263,
    318,
]

# THE THREAD STACK, and it belongs here because this is where a run is launched.
#
# `-fopenmp` implies `-frecursive`, so the local arrays gfortran used to place in
# static storage become stack-allocated and therefore per-thread. That is what
# makes the threaded build work without declaring them, and it is also what makes
# the default stack too small: `exoplasim/notes/threads-instead-of-ranks.md`
# measured 463.1 MB of such locals at T170, the largest single one 14 MB, against
# an 8 MB OMP_STACKSIZE default.
#
# TWO SETTINGS, because they size different stacks. OMP_STACKSIZE sizes the
# NON-MASTER threads. The master runs on the process stack, which is the stack
# rlimit, and `bench_ab.py` measured a T127 master overrunning a 16 MB limit and
# segfaulting. Nine gate and bench scripts set both; the production launcher set
# neither, so the one setting that decides whether a threaded run at the top of
# the ladder survives at all existed only in the instruments.
#
# 512M IS BOUNDED BELOW RATHER THAN DERIVED, and the two are different claims.
# `stack_floor.py` walks the model's call graph out of the parallel region and
# sums the declared local arrays along the heaviest chain, which is the deepest
# CONCURRENT set of declared locals: the floor at the top of the ladder is a
# fraction of 512M, and the floor at every rung is checked here rather than
# assumed. What the floor does NOT count is the temporaries gfortran
# materialises for the whole-array and `where` expressions the physics is
# written in; nothing in the source fixes their size, so the margin between the
# floor and 512M is not derived and `world-p4b` is still what closes that.
# `exoplasim/notes/thread-stack-floor.md` carries the numbers and the check
# against a built binary.
#
# Applied to THIS PROCESS before the model is launched. ExoPlaSim builds its
# command as a shell string and runs it with shell=True, so the child inherits
# both and neither the model package nor the namelist has to know.
OMP_STACKSIZE = "512M"


def _stack_bytes(value: str) -> int:
    """An OMP_STACKSIZE string in bytes. Kilobytes with no unit, per the spec."""
    text = str(value).strip()
    unit = {"b": 1, "k": 1024, "m": 1024**2, "g": 1024**3}
    factor = unit.get(text[-1:].lower(), 1024)
    digits = text[:-1] if text[-1:].lower() in unit else text
    return int(float(digits) * factor)


def run_executable_name(model_cfg: dict) -> str:
    """The registry's executable name for a run's configuration.

    ExoPlaSim copies its whole run directory in, so every previously built
    executable is present and a glob picks the wrong one -- `p8` sorts after
    `p16`. Composed the way `build_model.executable_name` composes it, with no
    parallel-mode suffix, and stated once here rather than in each of the three
    places that used to spell it out.
    """
    return (f"most_plasim_t{int(str(model_cfg['resolution']).lstrip('Tt'))}"
            f"_l{int(model_cfg['layers'])}_p{int(model_cfg['ncpus'])}.x")


def check_registered(exe: Path) -> dict:
    """Refuse an executable `binary_manifest.json` does not register.

    AT PREPARE, and refusing rather than warning. A run has cost nothing at this
    point, and the file it is about to be handed is the one every orbit after it
    is integrated by; an orbit is the most expensive thing in this project to
    have to integrate again. The registry gate exists -- `check_consistency.py`
    calls `rebuild_binaries.verify()` and rule 8 asks for it before an expensive
    run -- but it is a gate a caller has to remember, and this is the same
    question asked where it cannot be skipped. `unregistered()` names BOTH shas,
    because "stale" without the two numbers cannot be told from "the manifest
    was never written". world-bdb5.
    """
    reason = rebuild_binaries.unregistered(exe)
    if reason is not None:
        raise RuntimeError(
            f"{reason}.\n"
            f"An unregistered binary has unknown provenance: nothing says what "
            f"model source or what flag line produced it, so every orbit this "
            f"run integrates would be unattributable. Rebuild with "
            f"`python exoplasim/scripts/rebuild_binaries.py` and check with "
            f"`--verify`; to integrate a deliberate ARM instead, build it with "
            f"`build_model.py --no-publish` and pass its path to --binary.")
    entry = rebuild_binaries.registered(exe.name)
    return {"arm": None,
            "manifest_entry": {
                "name": exe.name,
                "sha256": entry.get("sha256"),
                "profile": entry.get("profile"),
                "effective_f90_opts": entry.get("effective_f90_opts"),
                "manifest": str(rebuild_binaries.MANIFEST),
            }}


def install_run_executable(run_dir: Path, model_cfg: dict,
                           arm: Path | None) -> tuple[Path, dict]:
    """Settle which binary this run will be integrated by, and say so.

    Two routes, and the run manifest records which one was taken.

    THE REGISTRY'S, by name, which is every ordinary run. `exo.Model` has
    already copied it into the run directory; this checks that copy against
    `binary_manifest.json` and refuses if it is not an entry.

    AN ARM, named by the caller. `build_model.py` can build an arm that departs
    from the shipped model in precision, in flags, in frame pointers or in a
    patched model source, and `--no-publish` leaves it in its own build
    directory -- which is the point, because the registry's naming carries none
    of those and an arm published under it would be indistinguishable from the
    shipped binary (archive CLIM-22, world-v3d). That left an arm buildable and
    unrunnable. The route is opened here without weakening the refusal: the arm
    is copied OVER the run directory's own copy, under the registry's name,
    inside this run directory and nowhere else. The published binary is not
    touched, the registry is not written, and the arm's build tag goes on the
    manifest so the run cannot later be read as the shipped model. world-u5pf.
    """
    exe_path = run_dir / run_executable_name(model_cfg)
    if not exe_path.is_file():
        raise RuntimeError(
            f"expected executable {exe_path} is not in the run directory")
    if arm is None:
        provenance = check_registered(exe_path)
        entry = provenance["manifest_entry"]
        print(f"executable {exe_path.name} {str(entry['sha256'])[:16]}: "
              f"binary_manifest.json entry, profile {entry['profile']}")
        return exe_path, provenance
    identity = build_model.arm_identity(arm)
    shutil.copy2(identity["path"], exe_path)
    print(f"ARM: {identity['build_tag']} ({file_sha256(exe_path)[:16]}) copied "
          f"over {exe_path.name} in this run directory. It is not the "
          f"registry's binary, it is absent from binary_manifest.json by "
          f"construction, and this run is stamped "
          f"canonical_lineage_eligible = false.")
    return exe_path, {"arm": identity, "manifest_entry": None}


def restore_run_executable(run_dir: Path, model_cfg: dict,
                           manifest: dict) -> tuple[Path, dict]:
    """Put back the binary a prepared run was given, before a resume launches it.

    `exo.Model.__init__` copies the registry's executable into the working
    directory on EVERY construction, and a continuation constructs one. So a run
    prepared on an arm has its arm overwritten by the shipped model at the start
    of its second segment, and nothing in the run would say so: the segment
    record would carry a sha, that sha would be the registry's, and the run
    would read as an arm experiment that had quietly been half a production run.
    Copying the arm back in here is what makes an arm run resumable at all.

    For a run that is NOT an arm this is the prepare's registry check asked
    again, which is the same question at the same boundary: a continuation is
    another chance to hand the model a binary of unknown provenance, and rule 4
    makes that likely rather than exotic, because any change under
    `vendor/exoplasim` is meant to be followed by a rebuild. A rebuilt binary is
    a registered one and passes; an unregistered one is refused before the
    segment costs anything.
    """
    exe_path = run_dir / run_executable_name(model_cfg)
    arm = (manifest.get("executable") or {}).get("arm")
    if arm is None:
        return exe_path, check_registered(exe_path)
    src = Path(arm["path"])
    if not src.is_file():
        raise RuntimeError(
            f"{identifier_of(manifest)} was prepared on the arm "
            f"{arm['build_tag']}, whose executable {src} is gone. An arm lives "
            f"in its build directory and nothing copies it anywhere durable, so "
            f"rebuild it with the same build_model.py arguments -- the "
            f"directory name states them -- or this run cannot be continued as "
            f"the experiment it is.")
    identity = build_model.arm_identity(src)
    if identity["build_tag"] != arm.get("build_tag"):
        raise RuntimeError(
            f"--binary provenance moved: this run was prepared on arm "
            f"{arm.get('build_tag')} and {src} is now {identity['build_tag']}.")
    got, want = file_sha256(src), (manifest.get("executable") or {}).get("sha256")
    if want and got != want:
        raise RuntimeError(
            f"the arm {arm['build_tag']} is {got[:16]} and this run was "
            f"prepared on {want[:16]}: it has been rebuilt under the same tag. "
            f"The tag names every input that changes a byte, so a moved sha "
            f"under an unchanged tag is the model source having moved. Prepare "
            f"a new run rather than splicing two models into one trajectory.")
    shutil.copy2(src, exe_path)
    return exe_path, {"arm": identity, "manifest_entry": None}


def identifier_of(manifest: dict) -> str:
    return str(manifest.get("run_id") or "this run")


def prepare_thread_stack(config: dict) -> dict:
    """Give the threaded model the stack its locals need, and check it is enough.

    Returns what was set, for the run manifest. Unconditional: world-38b left
    one build and it is the threaded one, so there is no longer a parallel mode
    that keeps these locals in static storage and never touches a thread stack.

    TWO STACKS, SET TWO DIFFERENT WAYS. The threaded build compiles with
    `-frecursive`, so the model's large locals are stack-allocated rather than
    static -- which the deleted MPI build never needed, because without it those
    arrays lived in static storage. The OpenMP team's threads take
    `OMP_STACKSIZE`; the MASTER thread runs on the process stack, which takes
    `ulimit -s` and which no OMP variable can change. Both are raised here.

    Neither was set at all until recently, and the symptom was not a crash
    anyone read as one: every rung above T42 died with SIGSEGV before writing a
    record, which read as the model REFUSING that configuration and is what kept
    T85, T127 and T170 out of the stability grid. It was a 16 MB process stack.

    Refuses when the stack is below the floor `stack_floor.py` reads off the
    model source at this rung. That is a check with a right answer rather than a
    number carried forward: a run that cannot hold its own declared locals dies
    inside the model, which reads as the model crashing.
    """
    import os
    import resource

    import stack_floor
    m = config["model"]
    threads = int(m["ncpus"])
    nlat, _, _ = rungs.geometry(str(m["resolution"]))
    params = stack_floor.parameters(nlat, int(m["layers"]), threads)
    frames, _, calls, defined = stack_floor.parse(int(m["precision_bytes"]), params)
    (floor, chain), _ = stack_floor.heaviest_chain(frames, calls, defined)

    os.environ["OMP_STACKSIZE"] = os.environ.get("OMP_STACKSIZE", OMP_STACKSIZE)
    record = {"applied": True, "omp_stacksize": os.environ["OMP_STACKSIZE"],
              "declared_local_floor_bytes": floor,
              "declared_local_floor_chain": chain}
    if _stack_bytes(record["omp_stacksize"]) < floor:
        raise RuntimeError(
            f"OMP_STACKSIZE is {record['omp_stacksize']} and the declared local "
            f"arrays on the heaviest chain out of the parallel region come to "
            f"{floor / 1e6:.1f} MB at {m['resolution']} on {threads} threads "
            f"({' -> '.join(chain)}). That floor counts no compiler temporaries, "
            "so a stack at or near it is already too small. "
            "exoplasim/scripts/stack_floor.py, world-p4b.")
    soft, hard = resource.getrlimit(resource.RLIMIT_STACK)
    try:
        resource.setrlimit(resource.RLIMIT_STACK, (hard, hard))
        record["stack_rlimit_soft"] = ("unlimited" if hard == resource.RLIM_INFINITY
                                       else hard)
    except (ValueError, OSError) as exc:
        # Reported rather than swallowed: a master thread that cannot get its
        # stack segfaults inside the model, which reads as the model crashing.
        record["stack_rlimit_soft"] = soft
        record["stack_rlimit_error"] = str(exc)
    return record


# High-cadence output is a gust distribution and nothing else, so this is the
# near-surface wind alone: `spd` is what the Weibull is fitted to, and `ua`/`va`
# are kept because a direction is what distinguishes a real gust from a reversing
# mean. DUST-5.
#
# This trims the POSTPROCESSED file. The raw one the model writes is 15 GB for a
# single T42 orbit at one sample every fourth timestep regardless of this list,
# because the model's output path takes no field list; it is cleaned up after
# pyburn runs, but the disk has to be free first.
HIGH_CADENCE_CODES = [131, 132, 259]

# The shortwave gas band weights, and the ONE list of them. Lacis and Hansen's
# absorptances are fractions of SOLAR flux, so a non-solar host needs every one
# re-weighted; radmod.f90 carries a key per term and this is where the config
# reaches them.
#
# Each is written only when it differs from the model's OWN default, so a run's
# namelist says what departs from the scheme as shipped -- and the defaults are
# not all 1.0. co2sww's is 0.0, because upstream has no shortwave CO2 term at
# all and zero is what reproduces upstream; comparing it against 1.0 would
# silently drop a weight of 1.0 and silently write one of 0.0, which is the
# inversion of what is meant.
#
# IT LIVES HERE BECAUSE continue_exoplasim.py IMPORTS IT. configure() rewrites
# the namelist on every continuation, so each key has to be reapplied per
# segment, and for a while that was a second copy of this tuple. PHYS-9 then
# added h2o_sw_level to this copy and not to that one, so H2OSWL applied for the
# orbits run_exoplasim.py prepared and silently reverted to the model default on
# every continuation after -- a physics change partway through a run, which is
# the failure mode both copies' comments called the hardest to notice in a long
# spin-up. One list cannot drift from itself.
SHORTWAVE_GAS_KEYS = (
    ("ozone_uv_weight", "O3UVW", 1.0),
    ("ozone_visible_weight", "O3VISW", 1.0),
    ("h2o_sw_weight", "H2OSWW", 1.0),
    # A LEVEL, not a weight, and separate from H2OSWW on purpose: that one is a
    # star-over-Sun ratio and an error in Eq. 21's absolute level divides out of
    # it. PHYS-9.
    ("h2o_sw_level", "H2OSWL", 1.0),
    ("co2_sw_weight", "CO2SWW", 0.0),
)

# The LONGWAVE trace gases, CLIM-42. Same one-list rule and the same reason:
# `continue_exoplasim.py` imports this, because a key applied at prepare and not
# reapplied per segment is a physics change partway through a run.
#
# These come from the `atmosphere` block rather than `model`, because they are
# composition and not a scheme weight -- pCO2_bar's neighbours. radmod_nl wants
# a VOLUME MIXING RATIO in ppmv, so the conversion is the partial pressure over
# the total, and the model's own default is 0.0, meaning absent.
TRACE_GAS_KEYS = (
    ("pCH4_bar", "CH4", 0.0),
    ("pN2O_bar", "N2O", 0.0),
)


def trace_gas_ppmv(config: dict) -> dict:
    """{namelist key: ppmv} for the longwave trace gases the config declares."""
    atmosphere = config.get("atmosphere", {})
    total = sum(float(v) for k, v in atmosphere.items()
                if k.startswith("p") and k.endswith("_bar"))
    if total <= 0.0:
        raise SystemExit("atmosphere block declares no partial pressures")
    out = {}
    for key, name, default in TRACE_GAS_KEYS:
        bar = atmosphere.get(key)
        if bar is None:
            continue
        ppmv = 1e6 * float(bar) / total
        if ppmv != default:
            out[name] = ppmv
    return out


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_version(command: list[str]) -> str | None:
    try:
        result = subprocess.run(command, check=True, text=True, capture_output=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    return (result.stdout or result.stderr).splitlines()[0]


def configure_otherargs(derived: dict) -> dict:
    """Namelist keys passed through `configure()`, for a prepare OR a resume.

    ONE definition, because constructing the model over an existing run
    directory re-copies the shipped namelists over the configured ones, so
    every one of these has to be reapplied per segment. That was a second copy
    in `continue_exoplasim.py` carrying only N_DAYS_PER_YEAR, so CLIM-17's
    TFREEZE applied to the orbits this script prepared and reverted to
    icemod.f90's compiled Earth value on every continuation after.

    At the declared salinity that particular revert is worth 5 mK and nothing
    else, but the salinity BRACKET arms are worth -0.04 to -0.36 K and are run
    as continuation segments: the key would have silently not applied, and an
    A/B designed to have a right answer would have measured zero and looked
    like it had confirmed one.
    """
    lwc = derived["land_water_column"]
    return {
        "N_DAYS_PER_YEAR@plasim_namelist": str(
            derived["rotations_per_orbit_namelist"]
        ),
        # Sea water's freezing point, from the declared salinity rather than
        # from icemod.f90's compiled-in Earth value. `configure()` has no
        # parameter for it and it is an ordinary `icemod_nl` key, so it goes
        # the same route N_DAYS_PER_YEAR does. CLIM-17.
        "TFREEZE@icemod_namelist": f"{derived['sea_water_freezing_point_k']:.4f}",
        # THE OROGRAPHY IS THE STAGED EXPORT AND NOTHING RESCALES IT, declared
        # rather than inherited from p_earth.f90's compiled 1.0. This project's
        # terrain comes from World Orogen and reaches the model as code 129, so
        # a scale factor on it is never a setting this world wants -- but it is
        # a planet_nl key nothing was writing, so the value the model integrated
        # was whatever the vendored default happened to be. Written here rather
        # than through `configure(orography=...)`, whose branch also forces
        # NGLACIER and would fight the glacier settings above. world-6qee.
        "OROSCALE@planet_namelist": "1.0",
        # THE OTHER THREE WAYS SALINITY REACHES THE MODEL, and they travel with
        # the freezing point rather than behind it: a bracket that moves only
        # TFREEZE moves one of four. Density is the mixed-layer heat capacity
        # with the specific heat, and it is the snow-ice flooding threshold as
        # the difference from ice density, where a one per cent density error
        # is a ten per cent threshold error. world-9hb.
        "CRHOS@icemod_namelist": f"{derived['sea_water_density_kg_m3']:.2f}",
        "CPS@icemod_namelist": f"{derived['sea_water_heat_capacity_j_kg_k']:.2f}",
        "CLFI@icemod_namelist": f"{derived['sea_ice_fusion_j_kg']:.1f}",
        # THE DECLARED COLD START, icemod_nl. world-6fh. icemod refuses a cold
        # start with no SST climatology and no declared profile, so these are
        # not optional; they are written on a resume too because constructing
        # the model over an existing run directory re-copies the shipped
        # namelists, and a segment that lost them would refuse to start.
        "TSST_EQ@icemod_namelist": f"{derived['cold_start_sst_equator_k']:.4f}",
        "TSST_POL@icemod_namelist": f"{derived['cold_start_sst_pole_k']:.4f}",
        "HICE_INI@icemod_namelist": f"{derived['cold_start_ice_thickness_m']:.4f}",
        # THE TWO SEA-ICE LENGTHS, icemod_nl. Written unconditionally,
        # defaults included, so the namelist in the run directory records what
        # the arm ran with rather than leaving it to the binary. world-12c.
        "XMAXD@icemod_namelist": f"{derived['sea_ice_max_thickness_m']:.4f}",
        "HLEAD@icemod_namelist": f"{derived['sea_ice_lead_closing_m']:.4f}",
        # WHAT THE MODELLED SEA ICE AND ITS SNOW ARE MADE OF, icemod_nl. All
        # four were compile-time parameters that no configuration could reach.
        # Written unconditionally so the namelist in the run directory records
        # the material the arm integrated rather than leaving it to the binary,
        # and so a bracket arm can move them. world-04ok.
        "CRHOI@icemod_namelist": f"{derived['sea_ice_density_kg_m3']:.2f}",
        "CPI@icemod_namelist": f"{derived['sea_ice_specific_heat_j_kg_k']:.2f}",
        "CKAPI@icemod_namelist": f"{derived['sea_ice_conductivity_w_m_k']:.4f}",
        "CLFSN@icemod_namelist": f"{derived['snow_fusion_j_kg']:.2f}",
        # CLIM-16, oceanmod_nl. Written unconditionally, defaults included, so
        # the namelist in the run directory records what the arm actually ran
        # with instead of leaving it to the binary's compiled value.
        "NHDIFF@oceanmod_namelist": str(derived["ocean_horizontal_diffusion"]),
        "HDIFFK@oceanmod_namelist":
            f"{derived['ocean_horizontal_diffusivity_m2_s']:.6g}",
        # PHYS-11, radmod_nl. It lands on CLOUDABS, which is the co-albedo
        # RATIO itself and not a coefficient carrying it: radmod's compiled
        # default is 1.0, the Sun, so scale 1.0 reproduces the tables as
        # printed.
        #
        # `swr` carries two shortwave cloud schemes and NSWRCL picks between
        # them. The computed branch, NSWRCL = 1, solves a two-stream layer
        # whose single-scattering co-albedo is interpolated from Stephens et
        # al. (1984) Table 1(a), and CLOUDABS multiplies that interpolated
        # co-albedo in both streams. The tables are flux-weighted over range 2
        # against the Sun, which is exactly the denominator
        # `cloud_band_weight.py` divides by, so the weight it derives is
        # unchanged by the table adoption: what changed is the key it arrives
        # on. It used to arrive as TSWR3, the coefficient of the analytic fit
        # the tables replaced, and that coefficient no longer exists.
        # world-f9ig.
        #
        # The prescribed branch, NSWRCL = 0, is the only reader of `acl2`, so
        # scaling `acl2` under NSWRCL = 1 wrote a number the modelled radiation
        # cannot consume. It is no longer written. Dropping it is bit-identical
        # for every run this project makes, and the +2.0 K the PHYS-11 arms
        # measured was always the co-albedo's: those arms ran NSWRCL = 1 too.
        # clim-68.
        "CLOUDABS@radmod_namelist":
            f"{derived['cloud_absorption_scale']:.6g}",
        # NSWRCL, written unconditionally at radmod's compiled default so the
        # run's namelist records WHICH cloud scheme the shortwave used. Nothing
        # in this project had ever named it, so the branch selection was
        # readable only from the model's own diagnostic echo, and a
        # continuation that flipped it would have replaced the two-stream
        # optics with a three-level lookup silently. clim-68.
        "NSWRCL@radmod_namelist": "1",
        # NCLOUDS, likewise at radmod's compiled default. It sits ABOVE NSWRCL:
        # at 0 the shortwave takes no cloud branch at all and NSWRCL selects
        # between two schemes neither of which runs, so a namelist that records
        # NSWRCL and not NCLOUDS does not say what the shortwave did. Nothing
        # in this project sets it -- `columnmode='clear'` is the only route and
        # no config here names one -- which is exactly why the run should say
        # so rather than leaving it to be inferred. world-35en.
        "NCLOUDS@radmod_namelist": "1",
        # world-nfh, landmod_nl. Snow seen through a canopy is a mixture of the
        # band's exposed snow with the band's canopy albedo, so the canopy
        # albedo has to arrive per band. landmod's compiled default is the
        # Earth-Sun broadband endmember in both bands.
        "ALBFOREST@landmod_namelist": ", ".join(
            f"{v:.6g}" for v in derived["vegetation_albedo_bands"]),
        # world-qpe, glacier_nl. See derive().
        "GLACPERSIST@glacier_namelist":
            f"{derived['glacier_persistence_orbits']:.6g}",
        # world-cwc, landmod_nl. configure(maxsnow=...) sets the same key and
        # run_exoplasim passes it, so that the run's exported cfg records it;
        # it is repeated HERE because continue_exoplasim.py's configure() call
        # has no maxsnow argument, and without this the snow cap would revert
        # to landmod's compiled 5 m on every continuation. That is CLIM-17's
        # failure exactly, on a key whose whole purpose is to be lifted.
        "DSMAX@landmod_namelist": f"{derived['max_snow_depth_m']:.6g}",
        # world-qvu, radmod_nl. Written unconditionally, defaults included, for
        # the same reason NHDIFF is: these were a bare literal in lwr until they
        # were named, so the namelist in the run directory is the only place a
        # reader can find out what the modelled surface emitted with.
        "ELWLAND@radmod_namelist":
            f"{derived['land_longwave_emissivity']:.6g}",
        "ELWSEA@radmod_namelist":
            f"{derived['sea_longwave_emissivity']:.6g}",
        # world-ayx, radmod_nl. Where the synthetic ozone profile sits. Written
        # here rather than through configure(ozone=dict) deliberately: that path
        # also rewrites A0O3, A1O3, ACO3 and TOFFO3 from the dict, which would
        # make four more Earth constants into transcriptions in this file. These
        # two are the only ones the gravity argument reaches.
        "BO3@radmod_namelist": f"{derived['ozone_height_m']:.6g}",
        "CO3@radmod_namelist": f"{derived['ozone_spread_m']:.6g}",
        # world-n1nu, rainmod_nl. Written unconditionally, the compiled value
        # included, for the reason NHDIFF and XMAXD are: `clwref` sets the
        # cloud water path that both the shortwave optical depth and the
        # longwave cloud emissivity are built from, and the arms of
        # exoplasim/notes/cloud-water-reference.md are read against a control
        # that has to be bit-identical to the compiled default. A key the run
        # directory does not record is a key no artifact can attribute an arm
        # to.
        "CLWREF@rainmod_namelist":
            f"{derived['cloud_water_reference_kg_m3']:.6g}",
        # world-trs3 and world-o12h, rainmod_nl. Written unconditionally, the
        # compiled sentinels included, on the same argument CLWREF above is
        # written on: both select between a DERIVED form and a literal, and a
        # run directory that does not record which one it used cannot be
        # attributed to either. Negative means derive.
        "GAMMA@rainmod_namelist":
            f"{derived['precip_reevaporation_gamma']:.6g}",
        "RCRITWIDTH@rainmod_namelist":
            f"{derived['cloud_fraction_subgrid_width']:.6g}",
        # world-80ia, fluxmod_nl. Written unconditionally, the compiled
        # sentinel included, on the same argument GAMMA above is written on: it
        # selects between a DERIVED length and a literal, and a run directory
        # that does not record which one it used cannot be attributed to
        # either. Negative means derive.
        "VDIFF_LAMM@fluxmod_namelist":
            f"{derived['asymptotic_mixing_length_m']:.6g}",
        # world-py6p, landmod_nl. THE LAND COLUMN, all eight keys and
        # unconditionally, so the namelist in the run directory says which
        # scheme the segment integrated. They travel together because they are
        # one selection: a continuation that reapplied NLANDWCOL and dropped
        # NLSOILW would run a layered column at landmod's compiled layer count
        # and shape, which is the bucket wearing the other scheme's name.
        "NLANDWCOL@landmod_namelist": str(lwc["NLANDWCOL"]),
        "NLSOILW@landmod_namelist": str(lwc["NLSOILW"]),
        "NLANDWDRAIN@landmod_namelist": str(lwc["NLANDWDRAIN"]),
        # PER-LAYER, and written at the declared layer count rather than at
        # NLSOILWX: landini reads the first `nlsoilw` entries and renormalises
        # the capacity shape over exactly those, so stating more would state a
        # shape the model does not read.
        "DSOILWF@landmod_namelist": ", ".join(
            f"{v:.6g}" for v in lwc["DSOILWF"]),
        "DSOILWZ@landmod_namelist": ", ".join(
            f"{v:.6g}" for v in lwc["DSOILWZ"]),
        "NLANDWPHASE@landmod_namelist": str(lwc["NLANDWPHASE"]),
        # THE EVAPORATION LIMITER'S TWO FREE AXES. The third, `drhsfull`, is
        # landmod's own and is not routed here; these two are the pair that
        # separates this model's limiter from cGENIE's ENTS, which is what
        # makes the disagreement a runtime bracket instead of a code fork.
        "DRHSLOW@landmod_namelist": f"{lwc['DRHSLOW']:.6g}",
        "NRHSEXP@landmod_namelist": str(lwc["NRHSEXP"]),
    }


def freezing_point_k(salinity_psu: float) -> float:
    """Freezing point of sea water at the surface, from its salinity.

    UNESCO (1983) / Millero's polynomial, the standard one, at zero gauge
    pressure. It exists because `icemod.f90` hardcodes `TFREEZE = 271.25` with
    the comment "at S=34.7", so the model carries Earth's ocean unless something
    tells it otherwise, and nothing in this project did until CLIM-17.

    The formula is checked against that constant rather than trusted: at
    S = 34.7 it must reproduce 271.25 K, which is a right answer the model
    already knows and therefore a test that can fail.
    """
    s = float(salinity_psu)
    if not 0.0 <= s <= 42.0:
        raise ValueError(f"salinity {s} psu is outside the range the UNESCO "
                         "polynomial is fitted over (0 to 42)")
    celsius = -0.0575 * s + 1.710523e-3 * s**1.5 - 2.154996e-4 * s**2
    return 273.15 + celsius


def seawater_density(salinity_psu: float, celsius: float) -> float:
    """Density of sea water at the surface, kg/m3.

    EOS-80, the one-atmosphere international equation of state of sea water
    (Millero and Poisson 1981), as printed in UNESCO technical paper 44. It
    exists because `icemod.f90` carried `CRHOS = 1030.` labelled "at S=34.7",
    and that constant is both half the mixed-layer heat capacity and, as the
    difference `CRHOS - CRHOI`, the whole of the snow-ice flooding threshold.

    Checked against EOS-80's own published check values below, which are right
    answers this implementation either reproduces or fails.
    """
    s = float(salinity_psu)
    if not 0.0 <= s <= 42.0:
        raise ValueError(f"salinity {s} psu is outside the range EOS-80 is "
                         "fitted over (0 to 42)")
    t = float(celsius)
    rho_w = (999.842594 + 6.793952e-2 * t - 9.095290e-3 * t**2
             + 1.001685e-4 * t**3 - 1.120083e-6 * t**4 + 6.536332e-9 * t**5)
    a = (8.24493e-1 - 4.0899e-3 * t + 7.6438e-5 * t**2 - 8.2467e-7 * t**3
         + 5.3875e-9 * t**4)
    b = -5.72466e-3 + 1.0227e-4 * t - 1.6546e-6 * t**2
    c = 4.8314e-4
    return rho_w + a * s + b * s**1.5 + c * s**2


def seawater_heat_capacity(salinity_psu: float, celsius: float) -> float:
    """Specific heat of sea water at the surface, J/(kg K).

    UNESCO (1983) / Millero et al. (1973), at zero gauge pressure. It exists
    because `oceanmod.f90` carried `CPS = 4180.` labelled "specific heat of sea
    water", and 4180 is FRESH water at about 25 C: the check below is that this
    polynomial reproduces the model's own constant at exactly the conditions
    that constant actually describes, which is a right answer that can fail.
    """
    s = float(salinity_psu)
    if not 0.0 <= s <= 42.0:
        raise ValueError(f"salinity {s} psu is outside the range the UNESCO "
                         "specific-heat polynomial is fitted over (0 to 42)")
    t = float(celsius)
    cp_pure = (4217.4 - 3.720283 * t + 0.1412855 * t**2
               - 2.654387e-3 * t**3 + 2.093236e-5 * t**4)
    a = -7.643575 + 0.1072763 * t - 1.38385e-3 * t**2
    b = 0.1770383 - 4.07718e-3 * t + 5.148e-5 * t**2
    return cp_pure + a * s + b * s**1.5


# EOS-80's published check values. Tolerance is the precision they are printed
# to; a mistyped coefficient moves the result far further than this.
for _s, _t, _want in ((35.0, 25.0, 1023.343), (35.0, 5.0, 1027.675),
                      (0.0, 5.0, 999.967)):
    if abs(seawater_density(_s, _t) - _want) > 0.001:
        raise RuntimeError(
            f"the EOS-80 density polynomial gives "
            f"{seawater_density(_s, _t):.4f} kg/m3 at S = {_s}, t = {_t} C, "
            f"where EOS-80's own check value is {_want}. One of the two is "
            "wrong and this is not a difference to average over.")

# The specific-heat polynomial against the constant it replaces, at the
# conditions that constant describes: fresh water near room temperature. If
# this passes, `CPS = 4180.` was fresh water's value and not sea water's.
_FRESH_WATER_CPS = 4180.0
if abs(seawater_heat_capacity(0.0, 25.0) - _FRESH_WATER_CPS) > 1.0:
    raise RuntimeError(
        f"the UNESCO specific-heat polynomial gives "
        f"{seawater_heat_capacity(0.0, 25.0):.2f} J/(kg K) for fresh water at "
        f"25 C, where oceanmod.f90's CPS said {_FRESH_WATER_CPS}. The claim "
        "that the compiled constant was fresh water's rests on this agreeing.")


_REFERENCE_SALINITY = 34.7
_REFERENCE_TFREEZE = 271.25
# Tolerance is the precision the model states the constant to, two decimals, and
# not tighter: the polynomial gives 271.2449 and 271.25 is that rounded. A
# tighter bound fails on the rounding rather than on a disagreement, which it
# did when this was first written at 0.005.
if abs(freezing_point_k(_REFERENCE_SALINITY) - _REFERENCE_TFREEZE) > 0.01:
    raise RuntimeError(
        f"the freezing-point polynomial gives "
        f"{freezing_point_k(_REFERENCE_SALINITY):.4f} K at S = "
        f"{_REFERENCE_SALINITY}, where icemod.f90's TFREEZE says "
        f"{_REFERENCE_TFREEZE}. One of the two is wrong and this is not a "
        "difference to average over.")


# COMPILED NAMELIST DEFAULTS, READ FROM THE MODEL RATHER THAN COPIED.
#
# `lib/sea_water.py`'s argument, applied to the keys this file routes: a
# constant copied out of the model into a script stops tracking the model the
# moment either moves, and nothing fails when it does. Reading it means a
# configuration that omits a key reproduces the compiled value by construction
# rather than by a transcription somebody has to re-check.
LANDMOD_SOURCE = MODEL_SRC / "plasim" / "src" / "landmod.f90"
LANDCOLUMN_SOURCE = MODEL_SRC / "plasim" / "src" / "landcolumn.f90"
RAINMOD_SOURCE = MODEL_SRC / "plasim" / "src" / "rainmod.f90"
FLUXMOD_SOURCE = MODEL_SRC / "plasim" / "src" / "fluxmod.f90"


def _fortran_default(source: Path, name: str) -> float:
    """One scalar module default out of a model source file.

    Raises if the model no longer declares the key: falling back to a literal
    here would put the transcription back, silently and at the compiled value's
    old number.
    """
    text = source.read_text(encoding="utf-8", errors="replace")
    m = re.search(rf"^\s*(?:integer|real)\s*::\s*{name}\s*=\s*([-+0-9.eEdD]+)",
                  text, re.IGNORECASE | re.MULTILINE)
    if not m:
        raise SystemExit(
            f"{source.name} no longer declares {name}. It is a namelist key "
            "this script writes, and its default must be read from the model, "
            "never copied into the script that writes it.")
    return float(m.group(1).replace("d", "e").replace("D", "e"))


def landmod_default(name: str) -> float:
    """A scalar `landmod_nl` default."""
    return _fortran_default(LANDMOD_SOURCE, name)


def rainmod_default(name: str) -> float:
    """A scalar `rainmod_nl` default."""
    return _fortran_default(RAINMOD_SOURCE, name)


def fluxmod_default(name: str) -> float:
    """A scalar `fluxmod_nl` default."""
    return _fortran_default(FLUXMOD_SOURCE, name)


# THE LAND LIQUID WATER COLUMN, `landmod_nl`. LSHY-3 and LSHY-5.
#
# Eight keys select the land column's liquid water scheme, its evaporation
# limiter and its soil phase, and until they had a route from
# `config/planet.yaml` the selection was reachable only by hand-editing a
# namelist in a run directory, which the next continuation would overwrite.
# Every registered hypothesis except the compiled default was unrunnable and
# the bracket LSHY-3 exists to make a runtime selection was still a code fork.
# world-py6p.


def landmod_default_array(name: str) -> list[float]:
    """An array `landmod_nl` default, read the same way, across continuations.

    `dsoilwf` and `dsoilwz` are declared as `(/ ... /)` constructors broken over
    two Fortran continuation lines, so the whole constructor is taken and the
    continuation markers stripped before the elements are split.
    """
    text = LANDMOD_SOURCE.read_text(encoding="utf-8", errors="replace")
    m = re.search(rf"^\s*real\s*::\s*{name}\s*\([^)]*\)\s*=\s*\(/(.*?)/\)",
                  text, re.IGNORECASE | re.MULTILINE | re.DOTALL)
    if not m:
        raise SystemExit(
            f"{LANDMOD_SOURCE.name} no longer declares the array {name}. It is "
            "a landmod_nl key and its default must be read from the model, "
            "never copied into the script that writes it.")
    body = re.sub(r"&\s*\n\s*&?", " ", m.group(1))
    return [float(v.replace("d", "e").replace("D", "e"))
            for v in body.replace("\n", " ").split(",") if v.strip()]


def land_water_layer_limit() -> int:
    """`NLSOILWX`, the compiled ceiling on the water column's layer count."""
    text = LANDCOLUMN_SOURCE.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"^\s*integer\s*,\s*parameter\s*::\s*NLSOILWX\s*=\s*(\d+)",
                  text, re.IGNORECASE | re.MULTILINE)
    if not m:
        raise SystemExit(
            f"{LANDCOLUMN_SOURCE.name} no longer declares NLSOILWX, which is "
            "the bound a declared layer count has to be checked against.")
    return int(m.group(1))


# The word each selection is DECLARED by in `config/planet.yaml`, against the
# integer `landmod_nl` takes. Words rather than integers because these are two
# named schemes and not a magnitude: `scheme: bucket` says in the config which
# of LSHY-3's registered hypotheses a run integrated, and 0 does not.
LAND_WATER_SCHEMES = {"bucket": 0, "layered": 1}
LAND_WATER_LOWER_BOUNDARIES = {"impermeable": 0, "free_drainage": 1}


def land_water_column(config: dict) -> dict:
    """The eight `landmod_nl` keys of the land column, from the configuration.

    Returns {NAMELIST KEY: value}, integers and floats and per-layer lists,
    with every key present: they are written unconditionally, defaults
    included, so the namelist in the run directory records which scheme the
    segment integrated rather than leaving it to the compiled binary. That is
    the same argument XMAXD, NHDIFF and GLACPERSIST are written on.

    Read by `derive()`, which passes it to `configure_otherargs`, and by
    `expected_namelist_keys`, which is a check on the STAGED FILE and shares
    this mapping the way it already shares `freezing_point_k`: what must not be
    shared with the staging code is the list of keys, and that list is stated
    separately in both places.

    The refusals mirror `landini`'s, one of them exactly and one it cannot
    make. `nlandwphase = 1` needs `nlandwcol = 1` because phase is a property
    of a layer; landini refuses that too, and refusing here costs a config read
    instead of a launched model. The agreement between the declared layer count
    and the length of the per-layer lists is the one landini CANNOT make: it
    reads the first `nlsoilw` entries of an array whose tail is compiled zeros,
    so a config that declares two layers and one capacity fraction would give
    layer two a zero capacity and run.
    """
    block = config.get("surface", {}).get("land_water_column", {}) or {}
    limit = land_water_layer_limit()

    scheme = str(block.get("scheme", "bucket"))
    if scheme not in LAND_WATER_SCHEMES:
        raise ValueError(
            f"surface.land_water_column.scheme {scheme!r} is not one of "
            f"{sorted(LAND_WATER_SCHEMES)}")
    boundary = str(block.get("lower_boundary", "impermeable"))
    if boundary not in LAND_WATER_LOWER_BOUNDARIES:
        raise ValueError(
            f"surface.land_water_column.lower_boundary {boundary!r} is not one "
            f"of {sorted(LAND_WATER_LOWER_BOUNDARIES)}")
    phase = bool(block.get("soil_phase", False))
    if phase and scheme != "layered":
        raise ValueError(
            "surface.land_water_column.soil_phase needs scheme: layered. "
            "Phase is a property of a LAYER and the scalar bucket has no layer "
            "to freeze; landmod.f90's landini refuses the same combination.")

    layers = int(block.get("layers", landmod_default("nlsoilw")))
    if not 1 <= layers <= limit:
        raise ValueError(
            f"surface.land_water_column.layers {layers} is outside 1 to "
            f"{limit}, which is landcolumn.f90's NLSOILWX")

    fractions = [float(v) for v in block.get(
        "layer_capacity_fraction", landmod_default_array("dsoilwf")[:layers])]
    thicknesses = [float(v) for v in block.get(
        "layer_thickness_m", landmod_default_array("dsoilwz")[:layers])]
    for name, values in (("layer_capacity_fraction", fractions),
                         ("layer_thickness_m", thicknesses)):
        if len(values) != layers:
            raise ValueError(
                f"surface.land_water_column.{name} has {len(values)} entries "
                f"and layers is {layers}. The model reads the first `nlsoilw` "
                "entries of an array whose tail is compiled zeros, so a short "
                "list runs with a zero rather than failing.")
    if any(v < 0.0 for v in fractions) or sum(fractions) <= 0.0:
        raise ValueError(
            "surface.land_water_column.layer_capacity_fraction must be "
            "non-negative and must not sum to zero; landini renormalises it to "
            "sum to one, so it is a SHAPE and its scale carries nothing.")
    if any(v <= 0.0 for v in thicknesses):
        raise ValueError(
            "surface.land_water_column.layer_thickness_m must be positive: it "
            "is a depth in metres, and the mapping onto the soil TEMPERATURE "
            "layers that decide phase is by the water layer's midpoint.")

    limiter = block.get("evaporation_limiter", {}) or {}
    return {
        "NLANDWCOL": LAND_WATER_SCHEMES[scheme],
        "NLSOILW": layers,
        "NLANDWDRAIN": LAND_WATER_LOWER_BOUNDARIES[boundary],
        "DSOILWF": fractions,
        "DSOILWZ": thicknesses,
        "NLANDWPHASE": int(phase),
        "DRHSLOW": float(limiter.get("theta_low", landmod_default("drhslow"))),
        "NRHSEXP": int(limiter.get("exponent", landmod_default("nrhsexp"))),
    }


# THE WORDS AGAINST THE MODEL'S OWN DEFAULT. `landmod.f90` states that the
# default land liquid water scheme is "the scheme that has always run", the
# scalar bucket with an impermeable base, and the two dictionaries above are
# the only place in this project where that claim is turned into an integer. An
# inverted mapping would put every run on the other scheme while its config and
# its namelist both read as the default, which is a difference no gate
# downstream could name. This is a right answer the model already knows, so it
# is a check that can fail.
for _word, _table, _key in (("bucket", LAND_WATER_SCHEMES, "nlandwcol"),
                            ("impermeable", LAND_WATER_LOWER_BOUNDARIES,
                             "nlandwdrain")):
    if _table[_word] != int(landmod_default(_key)):
        raise RuntimeError(
            f"this file maps {_word!r} to {_table[_word]} and landmod.f90 "
            f"declares {_key} = {landmod_default(_key):g}. Either the word is "
            "mapped to the wrong scheme or the model no longer defaults to the "
            "one that has always run; both change what a configuration that "
            "declares the default selects.")


def derive(config: dict, flux_ratio: float) -> dict:
    planet = config["planet"]
    star = config["star"]
    orbit = config["orbit"]
    model = config["model"]
    rotation_days = float(planet["rotation_hours"]) / 24.0
    # Gravity is taken from the configuration rather than derived from mass and
    # radius: World Orogen scaled this terrain's relief as 1/g, so the value it
    # ran with is canonical. Mass is checked against it so the two cannot drift.
    gravity = float(planet["gravity_m_s2"])
    implied_mass = gravity * float(planet["radius_earth"]) ** 2 / EARTH_STANDARD_GRAVITY
    if abs(implied_mass - float(planet["mass_earth"])) > 5e-4:
        raise ValueError(
            f"planet.mass_earth {planet['mass_earth']} disagrees with "
            f"planet.gravity_m_s2 {gravity} at radius {planet['radius_earth']} "
            f"(implied mass {implied_mass:.6f})"
        )
    semimajor_axis_au = math.sqrt(float(star["luminosity_solar"]) / flux_ratio)
    orbital_year_days = (
        EARTH_SIDEREAL_YEAR_DAYS
        * math.sqrt(semimajor_axis_au**3 / float(star["mass_solar"]))
    )
    timestep_minutes = float(model["timestep_minutes"])
    runsteps = round(orbital_year_days * 1440.0 / timestep_minutes)
    rotations_per_orbit = round(orbital_year_days / rotation_days)
    stellar_radius_solar = (
        math.sqrt(float(star["luminosity_solar"]))
        / (float(star["effective_temperature_k"]) / SOLAR_EFFECTIVE_TEMPERATURE_K) ** 2
    )
    return {
        "gravity_m_s2": gravity,
        "gravity_earth": gravity / EARTH_STANDARD_GRAVITY,
        "rotation_days": rotation_days,
        "stellar_radius_solar": stellar_radius_solar,
        "stellar_flux_ratio_earth": flux_ratio,
        "stellar_flux_w_m2": flux_ratio * float(orbit["earth_solar_constant_w_m2"]),
        "semimajor_axis_au": semimajor_axis_au,
        "orbital_year_earth_days": orbital_year_days,
        "orbital_year_seconds": orbital_year_days * 86400.0,
        "rotations_per_orbit_namelist": rotations_per_orbit,
        "runsteps_per_orbit": runsteps,
        "actual_run_duration_earth_days": runsteps * timestep_minutes / 1440.0,
        "snapshot_interval_steps": max(
            1, runsteps // int(model["seasonal_samples_per_orbit"])
        ),
        "ocean_salinity_psu": float(config["ocean"]["salinity_psu"]),
        "sea_water_freezing_point_k": freezing_point_k(
            config["ocean"]["salinity_psu"]),
        # world-9hb. Density and specific heat are evaluated AT the freezing point
        # of the declared salinity, not at some mean sea surface temperature:
        # the flooding threshold is a property of the water the modelled ice
        # floats in, and over the mixed layer's range the specific heat moves
        # by a tenth of a per cent. One temperature for both keeps the pair
        # consistent with each other and with TFREEZE.
        "sea_water_density_kg_m3": seawater_density(
            config["ocean"]["salinity_psu"],
            freezing_point_k(config["ocean"]["salinity_psu"]) - 273.15),
        "sea_water_heat_capacity_j_kg_k": seawater_heat_capacity(
            config["ocean"]["salinity_psu"],
            freezing_point_k(config["ocean"]["salinity_psu"]) - 273.15),
        # DECLARED, not derived: the heat of fusion of sea ice is a function of
        # the ice's own salinity and temperature, and icemod carries neither as
        # a variable.
        "sea_ice_fusion_j_kg": float(config["ocean"]["sea_ice_fusion_j_kg"]),
        # world-6fh. The declared cold start. DECLARED and not derived: an SST
        # field is what this model produces, so there is nothing here to
        # compute one from, and the alternative was a constructed Earth field.
        "cold_start_sst_equator_k": float(
            config["ocean"]["cold_start"]["sst_equator_k"]),
        "cold_start_sst_pole_k": float(
            config["ocean"]["cold_start"]["sst_pole_k"]),
        "cold_start_ice_thickness_m": float(
            config["ocean"]["cold_start"]["sea_ice_thickness_m"]),
        # world-12c. Negative means no maximum thickness. The lead-closing
        # scale is DECLARED at ExoPlaSim's value rather than recalibrated;
        # nothing here can recalibrate it.
        "sea_ice_max_thickness_m": float(
            config["surface"]["sea_ice_max_thickness_m"]),
        "sea_ice_lead_closing_m": float(
            config["surface"]["sea_ice_lead_closing_m"]),
        # world-04ok. The material the modelled sea ice and its snow are made
        # of. The three sea-ice properties are DECLARED, because each is a
        # function of a brine volume the model does not carry; the snow's
        # melting enthalpy is DERIVED, because snow has no brine in it and
        # IAPWS-06 therefore gives it exactly. Indexed rather than `.get`, so a
        # config that has lost the block fails here instead of silently
        # reverting the model to its compiled Earth values.
        "sea_ice_density_kg_m3": float(
            config["surface"]["cryosphere"]["sea_ice_density_kg_m3"]),
        "sea_ice_specific_heat_j_kg_k": float(
            config["surface"]["cryosphere"]["sea_ice_specific_heat_j_kg_k"]),
        "sea_ice_conductivity_w_m_k": float(
            config["surface"]["cryosphere"]["sea_ice_conductivity_w_m_k"]),
        "snow_fusion_j_kg": float(
            config["surface"]["cryosphere"]["snow_fusion_j_kg"]),
        # CLIM-16. Ocean horizontal heat transport EXISTS; a constant
        # diffusivity is a BOUND on the missing transport rather than the
        # transport, and it is bracketed rather than tuned. NLEV_OCE is 1
        # (oceanmod.f90:15) so hdiffk is a one-element array and a scalar
        # assignment fills it.
        "ocean_horizontal_diffusion": int(bool(
            config["ocean"].get("horizontal_diffusion", False))),
        "ocean_horizontal_diffusivity_m2_s": float(
            config["ocean"].get("horizontal_diffusivity_m2_s", 1.0e3)),
        # PHYS-11. Scales the range-2 cloud CO-ALBEDO and nothing else. The
        # backscatter fractions are tabulated alongside it and are scattering
        # rather than absorption, so they carry no stellar re-weighting, and
        # acllwr is a thermal-band constant with no stellar dependence; neither
        # gets an arm on this argument. 1.0 must reproduce radmod's compiled
        # value exactly, which is the arm the prediction says has to come out
        # bit-identical.
        "cloud_absorption_scale": float(
            config["model"].get("cloud_absorption_scale", 1.0)),
        # CLIM/world-nfh. The canopy albedo the two forested-snow endmembers
        # are mixed against, per band. landmod's default is the Earth-Sun
        # broadband endmember in BOTH bands, which is the thing the finding is
        # about, so this has to be passed for the fix to mean anything.
        "vegetation_albedo_bands": [
            float(v) for v in config["model"]["vegetation_albedo_bands"]],
        # world-qpe. Orbits of continuous snow cover at or above GLACELIM
        # before a gridcell becomes glacier. Written unconditionally so the
        # namelist in the run directory records the criterion the segment
        # actually ran with, rather than leaving it to the compiled default.
        "glacier_persistence_orbits": float(
            config["surface"].get("glaciers", {}).get(
                "persistence_orbits", 1.0)),
        # world-cwc. Metres of water equivalent, -1 for no limit.
        "max_snow_depth_m": float(
            config["surface"].get("glaciers", {}).get(
                "max_snow_depth_m", -1.0)),
        # world-qvu. The defaults here are radmod.f90's own, so a config that
        # says nothing reproduces the compiled literal exactly.
        "land_longwave_emissivity": float(
            config["model"].get("land_longwave_emissivity", 1.0)),
        "sea_longwave_emissivity": float(
            config["model"].get("sea_longwave_emissivity", 0.98)),
        # world-ayx. Defaults are radmod.f90's own, so a config that says
        # nothing keeps upstream's fixed geometric placement.
        "ozone_height_m": float(config["model"].get("ozone_height_m", 20000.0)),
        "ozone_spread_m": float(config["model"].get("ozone_spread_m", 5000.0)),
        # world-n1nu, rainmod_nl. rho_l0, the CCM3 reference in-cloud liquid
        # water density `mkclouds` anchors its exponential cloud water profile
        # on. DECLARED at CCM3's own value and read out of rainmod.f90 rather
        # than copied, so a config that says nothing reproduces the compiled
        # value exactly and the bracket arms move only what they declare.
        "cloud_water_reference_kg_m3": float(
            config["model"].get("cloud_water_reference_kg_m3",
                                rainmod_default("clwref"))),
        # world-trs3, rainmod_nl. The precipitation re-evaporation fraction.
        # NEGATIVE is the sentinel and the default: `rainini` then derives it
        # per cell and per level from Kessler (1969), and a POSITIVE value
        # restores the pre-derivation literal at all four re-evaporation sites.
        # Read out of rainmod.f90 for the same reason `clwref` is, so a config
        # that says nothing reproduces the compiled sentinel exactly. The route
        # exists because the derived form is a factor of 3.6 to 20 above the
        # constant it replaced over the fluxes that matter, which makes it the
        # largest single suspect in the bundle drift
        # `exoplasim/notes/forcing-bundle-predictions.md` measures, and until
        # now it was the one term in that bundle with a control and no way to
        # reach it.
        "precip_reevaporation_gamma": float(
            config["model"].get("precip_reevaporation_gamma",
                                rainmod_default("gamma"))),
        # world-o12h, rainmod_nl. The subgrid humidity width. NEGATIVE is the
        # sentinel and the default, deriving `(RCNLATREF/NLAT)^(1/3)` from the
        # grid; a POSITIVE value is used as the factor directly, and 1.0 is the
        # identity. At T21 the derived value IS 1.0, so the two agree there by
        # construction and the path is untestable at this rung rather than
        # verified at it.
        "cloud_fraction_subgrid_width": float(
            config["model"].get("cloud_fraction_subgrid_width",
                                rainmod_default("rcritwidth"))),
        # world-80ia, fluxmod_nl. The asymptotic mixing length of the
        # boundary-layer scheme, in metres. NEGATIVE is the sentinel and the
        # default: `fluxini` then derives it as `160 * (OMEGA_EARTH/ww)` from
        # Blackadar (1962) eq. 25, which at this world's rotation is 200.5 m
        # against the 160 m ECHAM anchors at Earth's; a POSITIVE value is the
        # length as declared and 160 restores the pre-derivation constant.
        # Read out of fluxmod.f90 for the same reason `gamma` is, so a config
        # that says nothing reproduces the compiled sentinel exactly. The route
        # exists because this is a SIXTH term in the batch-2 forcing bundle
        # that no registered prediction covers: it sets turbulent exchange over
        # land and ocean alike, it moved by 25 per cent, and without a key the
        # only control was a code fork.
        "asymptotic_mixing_length_m": float(
            config["model"].get("asymptotic_mixing_length_m",
                                fluxmod_default("vdiff_lamm"))),
        # world-py6p, landmod_nl. The eight keys of the land liquid water
        # column: which of LSHY-3's registered hypotheses runs, its evaporation
        # limiter, and LSHY-5's soil phase. See `land_water_column`.
        "land_water_column": land_water_column(config),
    }


def surface_sra(config: dict, code: int) -> Path:
    """Path to a generated surface field, keyed by the configured resolution."""
    resolution = str(config["model"]["resolution"]).upper()
    return INPUTS / resolution.lower() / f"orogen_{resolution}_surf_{code:04d}.sra"


def stellar_spectrum_path(config: dict) -> str | None:
    """Resolve the configured stellar spectrum to a file ExoPlaSim will accept.

    `configure(starspec=...)` needs a real path and both the file and its
    `_hr.dat` companion; it does not search the package. Returning None falls
    back to the blackbody at `startemp`.

    `inputs/stellarspectra/` is searched before the ExoPlaSim package, because
    spectra we generate are tracked here and the package lives in an untracked
    `.venv` that any reinstall resets. It also lets a project spectrum shadow a
    package one of the same name rather than silently losing to it.
    """
    name = config.get("radiation", {}).get("stellar_spectrum")
    if not name:
        return None
    roots = (
        INPUTS / "stellarspectra",
        Path(exo.__file__).resolve().parent / "stellarspectra",
    )
    for base in roots:
        path = base / f"{name}.dat"
        companion = base / f"{name}_hr.dat"
        if path.is_file() and companion.is_file():
            return str(path)
    available = sorted(
        {f.stem for base in roots for f in base.glob("*.dat")
         if not f.stem.endswith("_hr")}
    )
    raise RuntimeError(
        f"stellar spectrum {name!r} needs {name}.dat and {name}_hr.dat in one of "
        f"{[str(r) for r in roots]}, which do not both exist. Available: {available}"
    )


def stage_stellar_spectrum(model, run_dir: Path, spectrum: str | None) -> str | None:
    """Copy the spectrum into the run directory and shorten the namelist entry.

    `radmod.f90` declares `starfile` and `starfilehr` as `character(len=80)`,
    and `configure()` absolutises whatever path it is given. A venv path is
    easily longer than that: ours is 90 characters, and Fortran truncates it at
    80, emits a namelist *warning* rather than an error, then dies with an
    end-of-file inside `readdat`. The failure names a path that does not exist,
    which is not obviously a length problem.

    The model runs with the run directory as its working directory, so staging
    the two files there and referring to them by bare name sidesteps the limit
    entirely.
    """
    if not spectrum:
        return None
    src = Path(spectrum)
    hires = src.with_name(f"{src.stem}_hr.dat")
    for f in (src, hires):
        shutil.copyfile(f, run_dir / f.name)
    model._edit_namelist("radmod_namelist", "STARFILE", f"'{src.name}'")
    model._edit_namelist("radmod_namelist", "STARFILEHR", f"'{hires.name}'")
    return src.name


def verify_stellar_spectrum(model, config: dict) -> None:
    """Refuse to run unless the CONFIGURED spectrum is the one in the namelist.

    `radmod.f90:813` takes the spectrum branch only when `NSTARFILE > 0`;
    otherwise `solarini` builds a Planck curve at `STARBBTEMP` and says so
    nowhere except one line of `MOST_DIAG` that nothing reads. Constructing an
    `exo.Earthlike` on an existing run directory re-copies the shipped
    namelists over the configured ones, so any driver that rebuilds the model
    and calls `configure()` without `starspec` reverts the star in silence.

    `continue_exoplasim.py` and `run_stellar_cycle.py` both did. Every run on
    this build integrated its first orbit against `k25v` and every orbit after
    it against a 4965 K blackbody: the model's own log reports
    `Energy fraction below 0.75 microns` as 0.384383 in `MOST_DIAG.00000` and
    0.418350 from `MOST_DIAG.00001` onward, and that fraction weights the snow
    and sea-ice albedo, which moved the broadband snow albedo the model uses
    from 0.5382 to 0.5623.
    """
    name = config.get("radiation", {}).get("stellar_spectrum")
    if not name:
        return
    namelist = Path(model.workdir) / "radmod_namelist"
    entries = {}
    for line in namelist.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, _, value = line.partition("=")
            entries[key.strip().upper()] = value.strip().strip(" ,").strip("'\"")
    if entries.get("NSTARFILE") != "1" or not entries.get("STARFILEHR"):
        raise RuntimeError(
            f"config names stellar spectrum {name!r}, but {namelist} carries "
            f"NSTARFILE={entries.get('NSTARFILE')!r} and "
            f"STARFILEHR={entries.get('STARFILEHR')!r}. The model would run on a "
            "blackbody at STARBBTEMP instead. Pass starspec= to configure() and "
            "call stage_stellar_spectrum() after it.")
    # THE NAME IS COMPARED, not merely read. `name` was bound from the config
    # and then used only inside the error messages above, so the gate asserted
    # that A spectrum was configured rather than THE configured one, and a run
    # directory carrying another config's STARFILEHR passed the check written
    # against exactly this class of silent revert. Both entries are checked,
    # because `stage_stellar_spectrum` writes both and `solarini` reads both:
    # the low-resolution file sets the band split and the high-resolution one
    # the within-band weighting. world-60x0.
    want = {"STARFILE": f"{name}.dat", "STARFILEHR": f"{name}_hr.dat"}
    for key, expect in want.items():
        if entries.get(key) != expect:
            raise RuntimeError(
                f"config names stellar spectrum {name!r}, so {namelist} must "
                f"carry {key}={expect!r}; it carries {entries.get(key)!r}. The "
                "model would integrate against a different star from the one "
                "this run is configured for.")
        staged = Path(model.workdir) / expect
        if not staged.is_file():
            raise RuntimeError(
                f"radmod_namelist names {expect} but it is not in "
                f"{model.workdir}; readdat would die at end of file.")
def stellar_spectrum_digest(config: dict) -> dict | None:
    """Content digests of the two spectrum FILES the model will actually read.

    The config names a spectrum; the model reads a file. `build_stellar_spectrum.py`
    writes `<name>.dat` and `<name>_hr.dat` from `star.spectral_type` and
    `star.effective_temperature_k`, IN PLACE and under the same name, so the
    whole radiative input can change while every recorded name stays `k25v`.
    That difference is worth +0.024 on broadband snow albedo at the size it has
    already been measured at, and it moves every derived surface albedo with it.
    A name is a label here; the sha256 is the identity.

    `stellar_spectrum_path` has already proved both files exist, so this does
    not have to handle a half-present pair. Returns None only when the config
    names no spectrum, which means a blackbody at `startemp` and nothing to
    hash.

    This is the run-manifest twin of `spectrum_sha256` in
    `biosphere/generated/vesper_provenance.json`, deliberately: the biosphere
    half of this was already safe, and one mechanism recorded in two places
    beats a second mechanism.
    """
    path = stellar_spectrum_path(config)
    if not path:
        return None
    src = Path(path)
    hires = src.with_name(f"{src.stem}_hr.dat")
    return {
        "name": config.get("radiation", {}).get("stellar_spectrum"),
        "file": src.name,
        "sha256": file_sha256(src),
        "hr_file": hires.name,
        "hr_sha256": file_sha256(hires),
    }


def require_stellar_spectrum(manifest: dict, config: dict) -> bool:
    """Refuse to continue a run whose spectrum FILE has changed under it.

    The resume guard compares parsed config values, and `star.spectral_type` is
    on its inert list. That entry is CORRECT about the config key -- no script
    passes the spectral type to the model -- and WRONG about the artifact,
    because the value reaches the radiation through `k25v.dat`, which the key
    only names. Comparing the file is what closes that route; loosening the
    config comparison would not, since the spectrum can also be regenerated with
    no config edit at all.

    Returns True when the manifest was BACKFILLED and the caller must write it
    out. A run prepared before this existed carries no digest, and refusing
    those outright would make every existing run unresumable for a defect they
    predate. So an unstamped run is stamped with what it is about to run on,
    said out loud, and guarded from the next resume onward -- the same direction
    `lib/provenance.py:require_build` takes with `allow_unstamped`. Segments
    record their own digest either way, so which orbits are covered stays
    readable rather than being inferred from the top-level value.
    """
    current = stellar_spectrum_digest(config)
    recorded = manifest.get("stellar_spectrum_digest")
    if recorded is None:
        manifest["stellar_spectrum_digest"] = current
        print("  warning: this run was prepared before its spectrum was recorded "
              "by content, so the orbits already in it cannot be checked against "
              f"{(current or {}).get('file', 'a blackbody')}. Stamping it now; "
              "the next resume is guarded.")
        return True
    if recorded == current:
        return False
    differing = sorted(
        key for key in ("name", "file", "sha256", "hr_file", "hr_sha256")
        if (recorded or {}).get(key) != (current or {}).get(key))
    raise RuntimeError(
        "The stellar spectrum this run was prepared on is not the one on disk; "
        "refusing to resume:\n  "
        + "\n  ".join(f"{k}: {(recorded or {}).get(k)!r} -> {(current or {}).get(k)!r}"
                      for k in differing)
        + "\nThe spectrum weights every snow, ice and glacier albedo, so orbits "
          "either side of this are different climates. Start a new run rather "
          "than extending this one."
    )


def surface_input_paths(config: dict) -> list[Path]:
    """Every SRA file that defines this run's surface, in a stable order."""
    return [surface_sra(config, code)
            for code in sorted(intended_surface_codes(config))]


def geography_tag(config: dict) -> str:
    """Short digest of every surface input, carried in the run's identity.

    NOT the directory name any more: names are UUIDs, per CLAUDE.md rule 6,
    because a parameter-built name separates runs only along the dimensions it
    encodes and both ExoPlaSim and LPJ-GUESS collided that way. What this value
    is now is `physical_fingerprint["geography"]`, and it is load-bearing there
    -- it is the field that lets `INDEX.json` tell two runs on different surfaces
    apart, which is how the albedo bracket's cases are distinguishable at all.

    What it used to prevent is now prevented by the prepare guard. Two runs
    sharing a config landing in one directory made ExoPlaSim's finalize() select
    output as the last match of sorted(glob("MOST*")), so a shorter new run in a
    directory holding a longer old one copied out the previous run's final year
    under the new name, silently. The guard refuses that case directly.

    It covers albedo as well as mask and topography, because the planned
    experiment is exactly two runs differing only in albedo. Digesting the
    geography alone would put both endmembers in one directory.
    """
    digest = hashlib.sha256()
    for path in surface_input_paths(config):
        if not path.is_file():
            raise RuntimeError(
                f"{path} is missing; run build_boundary_conditions.py and "
                "build_surface_albedo.py before preparing a run"
            )
        digest.update(file_sha256(path).encode("ascii"))
    return digest.hexdigest()[:8]


# Surface fields we supply ourselves. Everything else falls back to a uniform
# namelist default, which for a world that is not Earth is mostly the right
# answer: Earth's roughness and vegetation maps are tied to Earth's continents.
# Albedo is the exception, because the substrate genuinely varies and we know how
# from the lithology, so it is supplied when land_albedo_source is not uniform.
BASE_SURFACE_CODES = {129, 172}
ALBEDO_SURFACE_CODES = {174, 175, 176, 212}
# dwmax, the soil water bucket whose overflow *is* ExoPlaSim's runoff. Supplied
# from pedology when asked for; otherwise the uniform namelist default stands.
SOIL_WATER_SURFACE_CODES = {229}
# dsoilwfc, the capacity split of that bucket over the land column's water
# layers. WORLD-VJBZ. Staged only when the column has more than one layer,
# because at one layer the split is identically one and `landmod_nl`'s dsoilwf
# says so exactly -- so requiring the file there would be a staging cost for a
# field that carries no information. Above one layer landmod REFUSES to run
# without it, which is the same position surfmod takes on an absent land mask.
SOIL_WATER_SPLIT_SURFACE_CODES = {2290}
ROUGHNESS_SURFACE_CODES = {173}
# ddustcol, the prescribed band-1 column dust optical depth. DUST-11. Only read
# by a binary carrying patches/exoplasim-3.4.2-prescribed-dust.patch; an
# unpatched one cannot parse NDUSTRAD and aborts in radini_, which is the loud
# failure this arrangement is designed to produce rather than avoid.
# CLIM-39 made ndustrad a COUNT of prescribed species rather than a switch, and
# the model now reads one surface field per species: species s is code 1810+s,
# up to NAERSP in radmod. Only 1811 is listed here because only 1811 has a
# generator -- `build_surface_dust.py`. Routing a second species in is CLIM-40,
# and it adds the code here in the same commit as the generator that writes it,
# per the rule that an artifact no step produces does not exist.
DUST_SURFACE_CODES = {1811}
# dsrcw, ddrage and dwpr, the source map for the INTERACTIVE emission scheme.
# DUST-3. A different thing from 1811 and not an alternative to it: 1811 is a
# prescribed column optical depth that nothing responds to, while these three
# are a source map that the model's own winds drive. Only read by a binary
# carrying patches/exoplasim-3.4.2-dust-emission.patch; an unpatched one cannot
# parse LDUSTEMIT and aborts in aero_ini, which is the loud failure this
# arrangement is designed to produce rather than avoid.
DUST_EMISSION_SURFACE_CODES = {1801, 1802, 1803}

# PlaSim's own 28-term energy decomposition, denergy(NHOR,28), written to these
# codes when nenergy > 0. Turned on to name the gap between the top of the
# atmosphere and the surface; it closed to 0.02 W/m2 and so ruled itself out, and
# is kept as the audit on the gridpoint physics. See
# notes/water-and-energy-closure.md.
#
# Term 15 needs patches/exoplasim-3.4.2-energy-diagnostics.patch applied and
# ExoPlaSim rebuilt, or it reports the latent heating with sublimation and
# vaporisation swapped.
ENERGY_DIAGNOSTIC_CODES = list(range(360, 388))

# The same 28 terms resolved per level, under nener3d. A column total says the
# atmosphere is losing energy; these say at which level, which is the difference
# between a number and a lead.
ENERGY_3D_CODES = list(range(460, 488))

# Which module each denergy term is accumulated in, from the assignment sites in
# the PlaSim source. Attribution only: the individual terms are not separately
# documented upstream, and naming them by physics would be inventing detail that
# has not been traced. Module plus index is enough to say where a residual is
# being introduced, which is the question being asked.
ENERGY_TERM_MODULE = {
    1: "plasim", 2: "plasim", 3: "plasim", 4: "plasim", 5: "plasim",
    6: "miscmod", 7: "fluxmod", 8: "fluxmod",
    9: "radmod", 10: "radmod",
    11: "rainmod", 12: "rainmod", 13: "rainmod", 14: "rainmod",
    15: "rainmod", 16: "rainmod",
    17: "radmod", 18: "radmod", 19: "radmod", 20: "radmod",
    21: "fluxmod", 22: "fluxmod",
    23: "plasim", 24: "plasim", 25: "plasim",
    26: "plasim", 27: "plasim", 28: "radmod",
}


def register_energy_diagnostic_codes() -> int:
    """Teach the postprocessor about codes 360-387 before it runs.

    pyburn resolves every requested code against its own `ilibrary`, which ships
    119 entries and none of these. Asking for them produces a run that completes,
    writes a restart, and then dies in postprocessing with "Going to stop here
    just in case", naming neither the code nor the reason. So enabling the
    diagnostics without this crashes at the last step of an otherwise good run.

    `ilibrary` is consulted at call time rather than captured at import, so
    inserting here is enough and the vendored tree stays untouched. That matters:
    `.venv` is untracked and any reinstall would silently drop a patch, whereas
    this travels with the script that depends on it.
    """
    from exoplasim import pyburn
    for i, code in enumerate(ENERGY_DIAGNOSTIC_CODES, start=1):
        module = ENERGY_TERM_MODULE.get(i, "unknown")
        pyburn.ilibrary.setdefault(str(code), [
            f"denergy{i:02d}",
            f"plasim_energy_budget_term_{i:02d}_{module}",
            "W m-2",
        ])
    for i, code in enumerate(ENERGY_3D_CODES, start=1):
        module = ENERGY_TERM_MODULE.get(i, "unknown")
        pyburn.ilibrary.setdefault(str(code), [
            f"dener3d{i:02d}",
            f"plasim_energy_budget_term_{i:02d}_{module}_per_level",
            "W m-2",
        ])
    return len(ENERGY_DIAGNOSTIC_CODES) + len(ENERGY_3D_CODES)


def load_conversion_report(restart: Path, resolved: Path) -> dict | None:
    """The report `convert_restart.py` wrote beside a converted restart.

    Looked for beside the file rather than passed as a flag, so a converted
    state carries its own provenance wherever it is moved to and cannot be
    seeded from without it. The report's own output hash must match the file:
    a report beside a DIFFERENT restart is worse than none, because it would
    vouch for a state it does not describe.
    """
    for candidate in (Path(str(restart) + ".conversion.json"),
                      restart.with_suffix(restart.suffix + ".conversion.json"),
                      restart.with_suffix(".conversion.json")):
        if candidate.is_file():
            report = json.loads(candidate.read_text(encoding="utf-8"))
            break
    else:
        return None
    got = file_sha256(resolved)
    want = (report.get("output") or {}).get("sha256")
    if want != got:
        raise RuntimeError(
            f"{candidate} describes a restart with sha256 {want} and "
            f"{resolved} has {got}. A conversion report beside the wrong "
            "restart vouches for a state it does not describe; regenerate it "
            "with convert_restart.py --report.")
    report["__path__"] = str(candidate)
    return report


def conversion_surface_reason(report: dict, config: dict) -> str | None:
    """Why this converted state's surface does not match what the run staged.

    A converted restart's static surface records are the TEMPLATE's, and the
    template was cut from a run that staged a particular set of `.sra` files.
    If this run stages different ones, landmod will read the template's from
    the restart and the staged files will be silently discarded -- which is the
    same failure the donor-run guard below exists for, reached by a different
    route.
    """
    template = report.get("target_template") or {}
    hashes = template.get("surface_field_sha256")
    if not hashes:
        return ("the conversion report carries no surface hashes for its "
                "target template, so nothing can say which surface the "
                "converted state froze")
    now = intended_surface_codes(config)
    transparent = DUST_SURFACE_CODES
    was = {int(c) for c in hashes}
    if (was - transparent) != (now - transparent):
        return (f"the template staged codes {sorted(was)} and this run stages "
                f"{sorted(now)}")
    changed = [c for c in sorted(now)
               if surface_sra(config, c).is_file()
               and hashes.get(str(c)) != file_sha256(surface_sra(config, c))]
    if changed:
        return f"content changed for code(s) {changed} since the template"
    return None


def intended_surface_codes(config: dict) -> set[int]:
    """Which surface fields this run supplies rather than leaving at defaults.

    `model.soil_water_source` defaults to `uniform` when the key is absent, so a
    config predating this option behaves exactly as it did.

    Setting it to `pedology` makes 229 mandatory here, and 229 comes from a soil
    weathered under a climatology. On a terrain that has none yet, the bootstrap
    run has to go out with the key at `uniform`; see `docs/src/pipeline/sequencing.md` loop A, which
    also says why the flip cannot happen mid-run.
    """
    codes = set(BASE_SURFACE_CODES)
    if str(config["model"].get("land_albedo_source", "uniform")) != "uniform":
        codes |= ALBEDO_SURFACE_CODES
    if str(config["model"].get("soil_water_source", "uniform")) != "uniform":
        codes |= SOIL_WATER_SURFACE_CODES
        if int((config.get("surface", {}).get("land_water_column", {})
                or {}).get("layers", 1)) > 1:
            codes |= SOIL_WATER_SPLIT_SURFACE_CODES
    if str(config["model"].get("roughness_source", "uniform")) != "uniform":
        codes |= ROUGHNESS_SURFACE_CODES
    if str(config["model"].get("dust_source", "none")) != "none":
        codes |= DUST_SURFACE_CODES
    if str(config["model"].get("dust_emission", "none")) != "none":
        codes |= DUST_EMISSION_SURFACE_CODES
    return codes


def energy_diagnostics_enabled(config: dict) -> bool:
    """Whether to ask PlaSim for its 28-term energy decomposition.

    Defaults False when the key is absent, so a config predating this option
    behaves exactly as it did and `config_sha256` does not move. Turning it on
    adds 28 output fields and is meant for short diagnostic segments, not for
    production runs.
    """
    return bool(config["model"].get("energy_diagnostics", False))


def energy_diagnostics_level(config: dict) -> int:
    """NENERGY: 1 for the 28-term decomposition, 2 to add the conversion control.

    `energy_diagnostics: 2` asks the model for the split of its adiabatic
    conversion across the semi-implicit scheme, which is a CONTROL for world-0ov
    and costs three extra spectral transforms a timestep. `true` and `1` both
    mean 1, so a config predating this reads exactly as it did.
    """
    level = config["model"].get("energy_diagnostics", False)
    if isinstance(level, bool):
        return 1 if level else 0
    return int(level)


def set_low_io(model, low_io: bool) -> None:
    """Turn off PlaSim's low-I/O accumulation, so output is instantaneous.

    `NLOWIO = 1` is PlaSim's default (`plasimmod.f90:151`) and accumulates fields
    over each output interval, dividing in place at write time. `NLOWIO = 0`
    writes instantaneous records instead, 182 an orbit rather than 12 bins, and
    `pyburn` averages them to the same 12 for the .nc, so nothing downstream
    changes shape.

    WHY THIS IS THE DEFAULT. An accumulation cannot be undone and a sample set
    can always be averaged, so samples are strictly more information. The
    accumulation is also not the mean you would compute yourself: binned `spd`
    under `NLOWIO = 1` sits between the speed of the time-mean vector and the
    mean of instantaneous speeds. Anything reading variance, extremes or single
    records -- DUST-5's gust distribution above all -- needs the samples, and a
    spin-up reading scalars does not, which is what `--low-io` in
    `continue_exoplasim.py` is for. See `exoplasim/README.md`.

    A SEPARATE DEFECT, NOW FIXED. `NLOWIO = 1` used to corrupt the first output
    record of every model call, because `naccuout` survived a restart while the
    accumulators did not: bottom-level wind read 7.5x the other bins and
    humidity 27% low, while scalars stayed within 2%.
    `exoplasim-3.4.2-lowio-first-record.patch` repairs that and was verified on
    2026-08-18 against a binary with it reversed, which cut the artifact by about
    four orders of magnitude. It does NOT change what is accumulated, so it does
    not make the two regimes interchangeable and it is not why this function
    exists. Runs written before it still carry the defect; see
    `exoplasim/notes/first-output-bin.md`.

    WHAT IT COSTS. Raw output is several times larger, and with the postprocessor
    fixed the clean regime is about 1.27x an orbit rather than the 3.7x once on
    record -- almost all of that gap was a quadratic reader in `pyburn` and not
    the I/O mode. Model time is identical either way. See
    `notes/audits/pyburn-postprocessing-cost.md`. Delete run directories once
    their climatologies are extracted.

    It has to be reapplied on every continuation, because `configure()` rewrites
    the namelist each time.

    WRITTEN EXPLICITLY IN BOTH DIRECTIONS, never left to the compiled default.
    `nlowio` defaults to 1 in `plasimmod.f90` and is read on NROOT only; the
    broadcast that makes it agree across ranks is ours, and a namelist that
    states the regime is the artifact that says which one ran. Leaving it unsaid
    is what made the deadlock in `notes/audits/nlowio-collective-deadlock.md`
    invisible for a day.
    """
    model._edit_namelist("plasim_namelist", "NLOWIO", "1" if low_io else "0")


def declare_cold_start_seed(model, config: dict, is_cold: bool) -> None:
    """Give a cold start a DECLARED initial kick instead of the clock's.

    `initrandom` at `plasim.f90:2068` reads the namelist `SEED` when `seed(1)`
    is non-zero and calls `system_clock` otherwise. Nothing here wrote `SEED`,
    so every cold-started run drew a different initial condition and no cold
    result could be checked by re-running it. The resolution-ladder stability
    boundary was measured entirely on cold starts, which makes each of its cells
    one draw rather than a verdict.

    Only for a cold start. A run seeded from a restart reads its state, so the
    kick never happens and writing a seed there would suggest a control this has
    no part in.
    """
    if not is_cold:
        return
    seed = int(config["model"].get("cold_start_seed", 0))
    if seed == 0:
        raise RuntimeError(
            "model.cold_start_seed is 0 or absent, which is the value that makes "
            "initrandom fall back to the system clock. Declare it; a cold run "
            "nobody can reproduce is a measurement nobody can check.")
    model._edit_namelist("plasim_namelist", "SEED", str(seed))
    print(f"cold start: SEED = {seed} (declared; without it initrandom takes "
          f"the system clock and the run is unreproducible)")


def read_applied_energy_fix(run_dir: Path) -> dict | None:
    """What the energy fixer actually had to put back, from the run's own diag.

    The fixer restores the energy the adiabatic conversion loses (`world-0ov`),
    and by doing so it HIDES that loss: with it on, `denergy26 - denergy27`
    reports a residual near zero and the size of the defect is the correction
    instead. So the correction has to be readable, or a known 0.85 W/m2 quietly
    becomes an unknown one that can grow. `world-mzy`.

    Taken from the model's own print rather than recomputed, so this cannot
    disagree with what was applied. The first samples are dropped: the first
    step out of a restart shows an imbalance of order 250 W/m2 and the rate
    limiter spends its first steps climbing out of it, which is startup and not
    the defect.
    """
    diags = sorted(run_dir.glob("MOST_DIAG.*"))
    if not diags:
        return None
    # Anchored on the header rather than counted from the end, because the line
    # carries a varying number of trailing diagnostics: after the label come the
    # step and then the applied W/m2, and that order is what the model writes.
    values = []
    for line in diags[-1].read_text(errors="replace").splitlines():
        if "ENERGY FIXER applied" not in line:
            continue
        tail = line.split("K/day", 1)[-1].split()
        if len(tail) < 2:
            continue
        try:
            values.append(float(tail[1]))
        except ValueError:
            continue
    if len(values) < 8:
        return None
    settled = values[len(values) // 2:]
    return {"applied_w_m2_mean": sum(settled) / len(settled),
            "applied_w_m2_min": min(settled),
            "applied_w_m2_max": max(settled),
            "samples": len(settled),
            "note": "the energy the fixer had to restore; a CORRECTION and not "
                    "physics. It is the size of the defect on world-0ov, which "
                    "the fixer otherwise hides. Second half of the run only."}


def declare_energy_fixer(model, config: dict) -> bool:
    """Switch on the global energy fixer, which is a CORRECTION and not physics.

    The adiabatic step must conserve total energy: the column enthalpy it gives
    up has to equal the kinetic energy it takes on. It does not. Measured on a
    dry adiabatic run it loses about 0.8 W/m2, and the whole of that is the
    reference conversion's two halves being taken at different time levels of
    the semi-implicit scheme: the advective half is explicit at t, the
    divergence half is applied on `sdt`, the centred mean of t-dt and t+dt, and
    the displacement between them is -0.96 W/m2. That defect is `world-0ov`.
    This does not fix it.

    What this does is put the missing energy back as a uniform warming, the way
    ECHAM, the IFS and CAM all do. Without it the surface silently supplies the
    shortfall -- across three filter settings the adiabatic term moved 2.08 W/m2
    while the closed budget moved 0.013, with sensible heat flux up 2.17 and
    latent up 1.48, about a tenth of the sensible flux on this planet -- and
    downstream components read that as physics.

    IT MASKS WHAT IT COMPENSATES, so the applied increment is printed by the
    model at `ndiag` cadence and belongs on the manifest as a health metric: a
    change in the underlying defect has to be able to show rather than be
    absorbed. `world-mzy`.

    The fixer reads `denergy26` and `denergy27`, so it needs the energy
    diagnostics. Rather than switch them on quietly underneath the caller, which
    would make a declared setting mean something it does not say, this raises.
    """
    fix = config["model"].get("energy_fixer")
    if fix is None:
        raise RuntimeError(
            "model.energy_fixer is absent. It is a correction for a known "
            "defect (world-0ov) and whether a run carries it changes what its "
            "surface fluxes mean, so it is declared and never defaulted.")
    if not fix:
        model._edit_namelist("plasim_namelist", "NENERGYFIX", "0")
        print("energy fixer: OFF (declared). The core's conversion loses about "
              "0.8 W/m2 and the surface supplies it; see world-0ov.")
        return False
    if not config["model"].get("energy_diagnostics"):
        raise RuntimeError(
            "model.energy_fixer needs model.energy_diagnostics, because the "
            "fixer is driven by denergy26 and denergy27. Declare both rather "
            "than having one switch the other on underneath you.")
    model._edit_namelist("plasim_namelist", "NENERGYFIX", "1")
    print("energy fixer: ON (declared). A CORRECTION, not physics: it restores "
          "the energy the adiabatic conversion loses without restoring where it "
          "went. world-mzy; the defect is world-0ov.")
    return True


def declare_robert_filter(model, config: dict) -> float | None:
    """PNU, the leapfrog time filter's coefficient. world-0ov.

    The leapfrog scheme carries a computational mode that alternates sign every
    step, and the Robert-Asselin filter is what damps it. It matters here beyond
    stability: the adiabatic sink is the reference conversion's two halves being
    taken at different time levels, the displacement between them is the second
    time difference of the divergence, and a mode that alternates sign every
    step contributes to that difference whatever the timestep is.

    IT IS IN `planet_nl` AND NOT `plasim_nl`, which is where a reader looks
    first: `pnu` is declared in `plasimmod.f90` at 0.0 and set to 0.1 by
    `p_earth.f90`'s `planet_ini`, so the value that runs comes from the planet
    module and the key that overrides it sits in the planet namelist. Written to
    the wrong one the model aborts with "Cannot match namelist object name pnu",
    which is the loud failure and not the silent one.

    Absent means the value the planet module sets, so a config predating this
    reads exactly as it did.
    """
    value = config["model"].get("robert_filter")
    if value is None:
        return None
    model._edit_namelist("planet_namelist", "PNU", repr(float(value)))
    print(f"Robert time filter: PNU = {float(value)}")
    return float(value)


def declare_dealias_conversion(model, config: dict) -> bool:
    """Truncate V.grad(ln ps) to the retained modes before the products. world-ly5.

    Hoskins and Simmons (1975) section 2, the paper `spectrala` cites: the
    transform grid "is insufficient for removing aliased interactions" for "the
    triple correlation involved in the energy conversion term and in the
    vertical advection terms", which need `M_g >= 4M + 1`. This model runs
    `NLON = 3*NTRU + 1`, the grid that dealiases a product of two band-limited
    fields and not of three, so the conversion is aliased by construction.

    ON, `calcgp` projects `zvgpg` onto the retained modes before anything
    multiplies it, which makes every product downstream quadratic in band
    limited fields again. That is HS75's own option (ii) and the cheaper of the
    two they name; the other is a `4M + 1` transform grid.

    IT CHANGES WHAT THE MODEL INTEGRATES, and it costs two extra transform pairs
    a timestep. It also has no SHTns path -- the projection uses the Legendre
    decomposition's partials -- and the model refuses rather than transforming
    through a path whose partials mean something else.

    Absent means off, so a config predating this reads exactly as it did.
    """
    if not config["model"].get("dealias_conversion", False):
        return False
    model._edit_namelist("plasim_namelist", "NDEALIAS", "1")
    print("dealiasing: V.grad(ln ps) is truncated to the retained modes before "
          "the conversion and the vertical advection use it. A CHANGE TO THE "
          "DYNAMICS; world-ly5.")
    return True


def declare_storm_diagnostics(model) -> dict:
    """Both halves of hurricanemod, declared OFF rather than left at a default.

    CLIM-54 asked whether the baseline carries the cyclogenesis index
    diagnostics and CLIM-55 asked that the storm-capture half stay a declared
    gap. This is where both answers live, because this is where a run meets the
    switch. Neither is a config key: making them settings would put a broken
    diagnostic back on the menu, and what is wrong with them is not a value a
    namelist can carry.

    NSTORMDIAG = 0. The index half is an ENVIRONMENTAL diagnostic and its
    resolution is not the problem -- indices of this kind exist because a model
    of this class cannot resolve a storm, so coarse fields are their intended
    support. What is wrong is that the indices are fitted to Earth in places a
    namelist cannot reach. `gpot` normalises the genesis potential index by
    compiled literals -- absolute vorticity against 1e-5 per second, relative
    humidity against 0.5, potential intensity against 70 m/s, shear against
    10 s -- none of which is in `hurricane_nl`; `vreducedvmax` folds `VITHRESH`,
    an Earth cyclogenesis cutoff, into the CONTINUOUS field on code 328 rather
    than into a mask; and `CPD` is hardcoded at 1005.7 instead of following
    `acpd`, while `CL` is an admitted fudge. Eight fields would be written and
    two of them would be Earth's calibration wearing this world's units. The
    thresholds that mask them are worse: `LAVTHRESH` is 25 per cent too strict
    against a planetary vorticity 0.8 of Earth's at a 30-hour day. Turning this
    on is a re-derivation of the indices, not a namelist key. world-9d1,
    world-khn, and the argument is `notes/audits/dormant-exoplasim-modules.md`
    finding 3.

    HC_CAPTURE = 0, and this one is a DECLARED GAP on DUST-16's precedent
    rather than a decision that could go the other way. The capture half hunts a
    RESOLVED vortex, and the bar is higher on this world than on Earth by two
    compounding factors: a given truncation buys 20 per cent coarser spacing in
    kilometres because the planet is 1.2 radii, and the storms are not
    correspondingly larger. The support Earth GCMs want before cyclone-like
    vortices reach realistic intensity is off the top of the resolution ladder
    this project has, and converting up does not put it on. A storm statistic
    produced by lowering `SIZETHRESH` and `WINDTHRESH` until something triggered
    would be a property of the grid rather than of the world. CLIM-55 is why no
    storm count exists and why that absence is not a null result.

    WRITTEN, NOT ASSUMED. Both keys already read 0 from the namelist ExoPlaSim
    ships, so this changes nothing the model integrates. What it changes is that
    the value is DECLARED: `expected_namelist_keys` checks both on every prepare
    and every continuation, so a `stormclim=True` reaching `configure()` from
    anywhere -- a loaded `.cfg`, a `modify()`, a future caller -- is refused
    instead of quietly adding eight Earth-calibrated fields to a segment.
    failure-modes class 22.
    """
    model._edit_namelist("hurricane_namelist", "NSTORMDIAG", "0")
    model._edit_namelist("hurricane_namelist", "HC_CAPTURE", "0")
    return {"nstormdiag": 0, "hc_capture": 0}


def declare_timestep(config: dict) -> float:
    """The step this run integrates at, resolved through the ladder registry.

    THE DEFECT THIS CLOSES. `model.timestep_minutes` was a bare scalar with no
    dependence on `model.resolution`, so changing the rung left the step where
    it was: the first T42 arm ran at T21's 45 minutes and took a SIGFPE in its
    forty-seventh orbit. WORLD-TD3, WORLD-J37.

    THE VALUE STILL COMES FROM THE CONFIG, and that is deliberate rather than a
    half-measure. The step selects nothing the registry owns -- it is an
    operational choice, and the arms that establish what a rung can take are
    exactly the ones that move it: `stability_probe.py` sweeps six steps per
    rung and `filter_timestep_matrix.py` runs T42 at dt 90 on purpose, to find
    a trap boundary. What moved into `lib/rungs.py` is the AUTHORITY: the
    coarsest step each rung is measured to start clean at, and the step the
    escalation route runs it at.

    SO THIS REPORTS AND DOES NOT REFUSE. The place a step is refused is where
    the DECLARED configuration is judged -- `scripts/check_consistency.py` and
    `scripts/smoke_test.py`, both of which read `config/planet.yaml` itself. A
    refusal here would refuse the measurements the registry is made of.

    IT ALSO PRINTS CAVEATS, which are the other half of what the registry knows.
    A pair can carry a blow-up that does not BAR it -- the run and its
    reproducer gone, so the claim has no artifact left to check it against --
    and that is not the same as a pair nothing has ever tried. Printing it here
    is what keeps an attempt at the first from looking like an attempt at the
    second. `lib/rungs.py:commissioning_caveats`.
    """
    timestep, rung = rungs.configured_timestep(config)
    problems = rungs.timestep_problems(rung, timestep)
    if not problems:
        print(f"timestep: {rung} at {timestep} min, on the escalation route "
              f"and at or below the measured ceiling "
              f"{rungs.stability_ceiling(rung)}")
    else:
        print(f"timestep: {rung} at {timestep} min, OFF THE DECLARED LADDER. "
              "A diagnostic arm is what this is for; a commissioning run at "
              "this step is not the escalation the project declared.")
        for problem in problems:
            print(f"  - {problem}")
    for caveat in rungs.commissioning_caveats(rung, timestep):
        print(f"timestep CAVEAT: {caveat}")
    return timestep


def declare_conversion_time_level(model, config: dict) -> bool:
    """Put the reference conversion's two halves on one time level. world-0ov.

    `calcgp` carries the advective half of the adiabatic reference conversion at
    time t and `spectrala` applies its divergence half on `sdt`, the centred
    mean of t-dt and t+dt. The identity that makes the two cancel the momentum
    equations' reference pressure-gradient work is the global integral of a
    mass-flux divergence, which holds at one time level and not across two, and
    the reference geopotential weights the split by up to 3.8 times `t0`. On the
    dry adiabatic arm the displacement is 0.96 W/m2 against a sink of 0.79.

    ON, this puts the divergence half back at time t. IT CHANGES WHAT THE MODEL
    INTEGRATES, and it takes that half out of the semi-implicit treatment in the
    temperature equation while the divergence solve still treats the temperature
    implicitly, so the timestep it is stable at is its own question and is NOT
    the one `lib/rungs.py`'s measured stability ceilings answer.

    Absent means off, so a config predating this reads exactly as it did.
    """
    if not config["model"].get("conversion_time_level", False):
        return False
    model._edit_namelist("plasim_namelist", "NCONVTIME", "1")
    print("conversion time level: the reference conversion's divergence half is "
          "taken at t, not on sdt. A CHANGE TO THE DYNAMICS; world-0ov.")
    return True


def declare_hyperdiffusion(model, config: dict) -> dict:
    """Write the derived horizontal diffusion, overriding the compiled branch.

    ExoPlaSim hard-codes these for T21 and T42 only (`plasim.f90:1544`), and
    every finer rung silently falls through to T21's -- which at T170 is four
    times too weak at the truncation and, worse, damps at 41% of the LOCAL
    CASCADE RATE at half the truncation, because `nhdiff` is an absolute
    wavenumber worth 71% of T21's spectrum and 8.8% of T170's.

    Setting them here is enough and needs no source change: `prolog` calls
    `readnl` then `initpm`, `readnl` applies the branch and THEN reads the
    namelist, and `initpm` builds the operator from whatever survived. Rule 4
    does not apply.

    `nhdiff` is an integer wavenumber, so the fraction is rounded; the rounding
    is reported rather than hidden, because at coarse truncations one wavenumber
    is a several-percent change in where the damping starts.

    exoplasim/notes/resolution-tuned-parameters.md has the derivation.
    """
    model_cfg = config["model"]
    hd = model_cfg.get("hyperdiffusion")
    if not hd:
        raise RuntimeError(
            "model.hyperdiffusion is absent. Without it the run inherits T21's "
            "damping at every rung above T42, which is the defect that block "
            "exists to fix; see exoplasim/notes/resolution-tuned-parameters.md.")
    rung = str(model_cfg["resolution"]).upper()
    table = hd["timescales_days"]
    if rung not in table:
        raise RuntimeError(
            f"model.hyperdiffusion.timescales_days has no entry for {rung}. The "
            f"damping for a rung is derived from the rule, so a rung with no "
            f"entry has no derived damping and must not fall back to T21's.")
    tau = table[rung]
    ntru = int(rung.lstrip("Tt"))
    nhdiff = int(round(float(hd["cutoff_fraction"]) * ntru))
    # EVERY LEVEL, NOT JUST THE FIRST. ndel and the four tdiss are ndel(NLEV)
    # and tdiss*(NLEV) in plasimmod, and a Fortran namelist scalar assigns
    # ELEMENT ONE and leaves the rest. Written as scalars, these reached the
    # model top and nine levels in ten kept whatever readnl had preset -- at
    # T21 `NDEL=4, 9*2`, so world-1nz's grad^8 landed on one level and the
    # other nine stayed on grad^4 with humidity damped 7.4 times too hard.
    #
    # At T42 it was worse than wrong, it was in the wrong UNIT. readnl's
    # NTRU==42 branch fills the arrays in SECONDS, and `dayseccheck` decides
    # days-against-seconds from MAXVAL over the whole array, so the preset's
    # 65664 suppressed the conversion element one needed and tdisst(1) stayed
    # 2.8224 SECONDS where 2.8224 days was meant: the top level damped 86400
    # times too hard, and only the top level. That is where the T42 blow-up of
    # world-td3 starts, and it is why T21, which has no such preset, survives.
    #
    # The repeat count is Fortran namelist syntax and the model echoes what it
    # read; smoke_test checks the echo rather than this line. world-720 is the
    # same defect found from the model side, and world-1nz is the continuation
    # that dropped these keys entirely; both are closed by this and by the
    # call in continue_exoplasim.py.
    nlev = int(model_cfg["layers"])
    def _every_level(value) -> str:
        return f"{nlev}*{value}"
    keys = {"NDEL": _every_level(int(hd['order_alpha'])),
            "NHDIFF": f"{nhdiff}",
            "TDISSD": _every_level(float(tau['divergence'])),
            "TDISSZ": _every_level(float(tau['vorticity'])),
            "TDISST": _every_level(float(tau['temperature'])),
            "TDISSQ": _every_level(float(tau['humidity']))}
    for key, value in keys.items():
        model._edit_namelist("plasim_namelist", key, value)
    print(f"hyperdiffusion: {rung} alpha={hd['order_alpha']} "
          f"nhdiff={nhdiff} (n*/N={nhdiff/ntru:.3f}, asked {hd['cutoff_fraction']}) "
          f"tau_vorticity={tau['vorticity']} d")
    return {"rung": rung, "nhdiff": nhdiff, "cutoff_fraction_actual": nhdiff / ntru,
            "order_alpha": int(hd["order_alpha"]), "timescales_days": dict(tau),
            "eddy_wind_m_s": float(hd["eddy_wind_m_s"])}


def declare_dry_constants(model, config: dict) -> dict:
    """The dry thermodynamic constants, the cold-start profile and the sponge.

    Four things ExoPlaSim otherwise takes from `p_earth.f90` and its compiled
    defaults, grouped because they are all read by the same dry column and all
    were inherited rather than chosen.

    AKAP. `p_earth.f90:41` declares R/cp as 0.286 and `readnl` derives the
    specific heat of dry air from it as `acpd = gascon/akap`. GASCON reaches the
    model derived from the declared composition and this did not, so the model
    ran a cp inconsistent with its own gas constant. Both come off
    `lib/lapse.py:gas_properties` here, which is the one place in this tree that
    turns the composition into (R, cp), and the R it returns is checked against
    the GASCON the model was actually configured with rather than assumed to
    agree. world-cwu.

    T0. The semi-implicit reference temperature, written as a per-level list
    because it is an array and a namelist scalar would set element 1 only.
    world-bmf.

    TGR, ALR, DTROP. The cold start's initial profile, and TGR additionally sets
    the orographic reduction of surface pressure on EVERY start. world-wmw.

    TFRC. Rayleigh drag on the top levels, in seconds, from a declaration in
    rotations. Written for every level so the value is the config's and not
    `readnl`'s unconditional NLEV==10 branch. world-aee.
    """
    import lapse

    model_cfg = config["model"]
    layers = int(model_cfg["layers"])
    gas_constant, cp = lapse.gas_properties(config)
    akap = gas_constant / cp
    # A CHECK THAT CAN FAIL: the same composition through two implementations.
    # If ExoPlaSim's mean molecular weight and this module's molar masses ever
    # disagree, akap and gascon stop being two views of one gas and the model
    # gets a cp that belongs to neither.
    configured = float(model.gascon)
    if abs(configured - gas_constant) > 0.01:
        raise RuntimeError(
            f"lib/lapse.py derives R = {gas_constant:.4f} J/kg/K from the "
            f"declared composition and ExoPlaSim was configured with GASCON = "
            f"{configured:.4f}. akap = R/cp is only meaningful beside the "
            "gascon it was derived with; reconcile the two before running.")

    profile = model_cfg["cold_start_profile"]
    t0 = float(model_cfg["semi_implicit_reference_temperature_k"])
    sponge = [float(x) for x in model_cfg["rayleigh_sponge_rotations"]]
    if len(sponge) != layers:
        raise ValueError(
            f"model.rayleigh_sponge_rotations has {len(sponge)} entries and the "
            f"model has {layers} levels. TFRC is per level and a short list is "
            "a level whose drag nobody declared.")
    # TFRC is entered in SECONDS. One rotation is the SIDEREAL day, which is
    # also the model's own unit of time after world-rt1, so a sponge declared in
    # rotations is the same number of nondimensional units at any rotation rate.
    rotation_s = float(config["planet"]["rotation_hours"]) * 3600.0
    # Fixed point rather than %g: a namelist real is parsed by the Fortran
    # runtime and there is no reason to hand it an exponent it has to read back.
    tfrc = ",".join(f"{x * rotation_s:.4f}" for x in sponge)

    model._edit_namelist("planet_namelist", "AKAP", f"{akap:.8g}")
    model._edit_namelist("planet_namelist", "ALR",
                         f"{float(profile['lapse_rate_k_per_m']):.8g}")
    model._edit_namelist("plasim_namelist", "TGR",
                         f"{float(profile['surface_temperature_k']):.8g}")
    model._edit_namelist("plasim_namelist", "DTROP",
                         f"{float(profile['tropopause_height_m']):.8g}")
    model._edit_namelist("plasim_namelist", "T0", f"{layers}*{t0:.8g}")
    model._edit_namelist("plasim_namelist", "TFRC", tfrc)
    print(f"dry constants: akap={akap:.6f} (cp={cp:.2f} J/kg/K beside gascon "
          f"{configured:.4f}), t0={t0:g} K, cold start "
          f"{profile['surface_temperature_k']:g} K / "
          f"{float(profile['lapse_rate_k_per_m']) * 1000:g} K/km / "
          f"{profile['tropopause_height_m']:g} m, "
          f"sponge {sponge[:2]} rotations")
    return {"akap": akap, "cp_j_kg_k": cp, "gascon_j_kg_k": configured,
            "t0_k": t0, "cold_start_profile": dict(profile),
            "rayleigh_sponge_rotations": sponge,
            "rayleigh_sponge_seconds": [x * rotation_s for x in sponge]}


def declare_dynamics_only(model, config: dict) -> bool:
    """Strip every diabatic source, leaving the dynamical core alone.

    A DIAGNOSTIC CONFIGURATION, not a climate. With radiation, vertical
    diffusion and surface fluxes off and precipitation, convection and dry
    convective adjustment switched out of `rainstep`, nothing heats or cools the
    atmosphere. Anything that then appears in the energy budget is the numerics,
    and -- the point of it -- there is no surface exchange left to absorb it.

    That matters because the adiabatic residual is otherwise INVISIBLE: measured
    across three filters, `26 - 27` moved 2.08 W/m2 while the closed budget moved
    0.013, because the surface silently supplied whatever the dynamics destroyed.
    Take the surface away and the same sink has nowhere to hide.

    THE DAMPING STAYS ON, deliberately. `26 - 27` is measured across `spectrala`
    and the hyperdiffusion is applied in `spectrald`, so the damping does not
    enter the identity; removing it would only make the run blow up sooner. The
    physics filter is left to the config, so it can be varied as its own arm.
    """
    if not bool(config["model"].get("dynamics_only", False)):
        return False
    # RADIATION STAYS SWITCHED ON AND IS EMPTIED INSTEAD. `nrad = 0` skips
    # `radstep` entirely, which leaves `eccf` at its initialised zero, and
    # `outsc_` then writes `1.0/sqrt(eccf)` -- a divide by zero that traps under
    # the project's floating-point traps at the first output record. Setting
    # `nswr` and `nlwr` to zero leaves the orbital geometry computed and removes
    # every radiative tendency, which is what this mode actually wants.
    # NOT `nflux`. It is accepted by `plasim_nl` and read by nothing -- the only
    # occurrence in the source is its own namelist declaration -- so setting it
    # is silent and does nothing. Measured: with `nflux = 0` the latent flux was
    # still 6.5 W/m2. The live switches are fluxmod's.
    for key in ("NEVAP", "NSHFL", "NSTRESS", "NVDIFF"):
        model._edit_namelist("fluxmod_namelist", key, "0")
    for key in ("NSWR", "NLWR"):
        model._edit_namelist("radmod_namelist", key, "0")
    for key in ("NPRL", "NPRC", "NDCA", "NSHALLOW"):
        model._edit_namelist("rainmod_namelist", key, "0")
    print("dynamics only: nswr=0 nlwr=0, no evaporation, sensible flux, "
          "surface stress or vertical diffusion, and large-scale, convective, "
          "shallow and dry-adjustment heating all off. Not a climate.")
    return True


def enable_energy_diagnostics(model, config: dict) -> bool:
    """Set nenergy in plasim_nl, which the Python API does not expose.

    Same mechanism `stage_stellar_spectrum` uses for STARFILE: the namelist is
    edited directly because `configure()` has no parameter for it.
    """
    if not energy_diagnostics_enabled(config):
        return False
    model._edit_namelist("plasim_namelist", "NENERGY",
                         str(energy_diagnostics_level(config)))
    if config["model"].get("energy_diagnostics_3d", False):
        model._edit_namelist("plasim_namelist", "NENER3D", "1")
    return True


# The ecological stream's own codes, `outmod.f90:ecogp`. 600 to 603 are the
# interval record -- the three bounds and the orbital position -- and 610 to 628
# the fields; `compare_eco_streams.FIELD_NAMES` is the one table of what each
# field IS. They are listed here only so a postprocessor code list can be
# refused for carrying one, and the refusal is below.
ECO_STREAM_CODES = frozenset(range(600, 604)) | frozenset(range(610, 629))


def refuse_eco_codes(codes: list[int], which: str) -> None:
    """The ecological stream does not go through pyburn, and this says so.

    NOT AN OVERSIGHT AND NOT A GAP TO FILL LATER. `pyburn.readallvariables`
    builds its time axis by counting records of code 139 and has no other
    notion of time, so a stream is postprocessable only if a consumer is
    willing to infer each record's interval from its position in the file. The
    ecological stream exists in part to refuse exactly that inference: it
    writes the interval start, end and duration with every block, in absolute
    seconds from the model's own step counter, because a consumer that infers
    an interval from a record number cannot tell a missing block from a short
    one. Routing it through pyburn would drop the bounds and restore the
    inference.

    So the stream is read directly. `exoplasim/scripts/compare_eco_streams.py`
    reads it today and EFOR-3 reads it next, both from the raw records.

    What this stops is the failure mode that would follow from asking anyway: a
    code pyburn's `ilibrary` does not carry makes it stop with "Going to stop
    here just in case", naming neither the code nor the reason, at the end of an
    otherwise good run.
    """
    bad = sorted(set(codes) & ECO_STREAM_CODES)
    if bad:
        raise ValueError(
            f"{which} carries ecological stream codes {bad}. That stream is "
            f"not postprocessed: it is written to its own unit with its own "
            f"interval bounds and is read directly, by "
            f"exoplasim/scripts/compare_eco_streams.py today and by EFOR-3 "
            f"next. See refuse_eco_codes.")


def declare_ecological_stream(model, enabled: bool, interval_steps: int | None) -> dict:
    """Set NECO and NECOSTEP in plasim_nl, which the Python API does not expose.

    THE STREAM. A second output stream written for the biosphere rather than
    for a climate diagnostic: nineteen surface fields plus its own interval
    bounds, reduced at the producer with the operator each field's meaning
    calls for. `exoplasim/README.md` has what it is for and
    `biosphere/notes/ecological-forcing-field-contract.md` what it carries.

    ON THE COMMAND LINE AND NOT IN `config/planet.yaml`, for the reason the
    benchmark overrides are: it cannot change a result. The restarts written
    with `NECO = 1` and `NECO = 0` are byte-identical, because the stream reads
    state that already exists and writes to a unit of its own. A config edit
    would move `config_sha256` and make every existing run unresumable, for a
    switch that changes only what is WRITTEN.

    THE INTERVAL IS DERIVED AND NOT TYPED. `NECOSTEP = 0` is the sentinel
    `plasim.f90` reads as `mtspd`, the number of timesteps in one absolute
    24-hour day EXACTLY: prolog derives `mtspd` from `day_24hr` and then
    recomputes `mpstep` so that `mtspd * mpstep * 60 = day_24hr`, so the
    interval is 86400 s with no rounding and no drift and the bounds `ecogp`
    writes are exact. Computing the same number here in Python would duplicate
    arithmetic that depends on the timestep and therefore on the resolution
    rung, and a duplicate that disagrees produces intervals whose declared
    bounds are right and whose contents are not.

    24 h is the ecological step LPJ-GUESS integrates on, not a claim about this
    world's rotation. `biosphere/notes/time-base-unit-contract.md` settles it.

    WRITTEN IN BOTH DIRECTIONS, never left to the compiled default, so the
    namelist on disk records which regime ran -- and reapplied on every
    continuation, because `configure()` rewrites the namelist each time.
    """
    if interval_steps is not None and interval_steps < 1:
        raise ValueError("--eco-interval-steps must be positive")
    step = 0 if interval_steps is None else int(interval_steps)
    model._edit_namelist("plasim_namelist", "NECO", "1" if enabled else "0")
    model._edit_namelist("plasim_namelist", "NECOSTEP", str(step))
    return {"enabled": bool(enabled), "necostep": step,
            "interval": ("mtspd, one absolute 24-hour day, derived by the model"
                         if step == 0 else f"{step} timesteps, declared")}


def enable_prescribed_dust(model, run_dir: Path, config: dict) -> dict | None:
    """Switch on the prescribed dust the patched radiation reads. DUST-11.

    Every value comes from the provenance file `build_surface_dust.py` wrote
    beside the field, so the run cannot be given a longwave ratio, a scale height
    or a burden that disagrees with the field it is applied to. There is no
    default for `DUSTQLW` here and none in the Fortran either: a shortwave-only
    dust is a larger error than no dust, so `radini` aborts rather than assume
    one, and this raises rather than supply one.

    `AEROFILE` goes in `radmod_namelist` rather than `aero_namelist`, because the
    prescribed path never runs `aero_ini` -- that is only called when the
    semi-Lagrangian tracer transport is on, and nothing here is transported. The
    file is staged into the run directory and named bare, for the same
    `character(len=80)` reason `stage_stellar_spectrum` does it.
    """
    if str(config["model"].get("dust_source", "none")) == "none":
        return None
    field = surface_sra(config, sorted(DUST_SURFACE_CODES)[0])
    prov_path = field.with_name(field.stem + "_provenance.json")
    if not prov_path.is_file():
        raise RuntimeError(
            f"{prov_path} is missing. Run exoplasim/scripts/build_surface_dust.py; "
            "the namelist values live with the field, not in the config.")
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    values = prov["namelist_values"]
    if float(values.get("DUSTQLW", 0.0)) <= 0.0:
        raise RuntimeError(
            f"{prov_path} carries no thermal-infrared absorption ratio. "
            "Shortwave-only dust is worse than no dust; see notes/dust.md.")
    aerofile = PROJECT_ROOT / prov["aerofile"]
    if not aerofile.is_file():
        raise RuntimeError(
            f"{aerofile} is missing. Run exoplasim/scripts/dust_aerofile.py.")
    shutil.copyfile(aerofile, run_dir / aerofile.name)
    model._edit_namelist("radmod_namelist", "AEROFILE", f"'{aerofile.name}'")
    scale = float(config["model"].get("dust_scale", values["DUSTSC"]))
    for key, value in (("NDUSTRAD", "1"),
                       ("DUSTSC", f"{scale}"),
                       ("DUSTHSC", f"{float(values['DUSTHSC'])}"),
                       ("DUSTQLW", f"{float(values['DUSTQLW'])}")):
        model._edit_namelist("radmod_namelist", key, value)
    return {"variant": prov["variant"], "dustsc": scale,
            "dusthsc": float(values["DUSTHSC"]),
            "dustqlw": float(values["DUSTQLW"]),
            "aerofile": aerofile.name,
            "field_sha256": prov["output_sha256"],
            "land_mean_band1_optical_depth":
                prov["field_statistics"]["land_mean_band1_optical_depth"]}


def enable_dust_emission(model, run_dir: Path, config: dict) -> dict | None:
    """Write the `aero_nl` group the in-model emission scheme reads. DUST-3.

    This is the half of the item that is easiest to skip and it is load-bearing:
    the emission patch and the deposition patch both add namelist keys with NO
    Fortran defaults, so until something writes this group none of them is
    reachable and `aero_ini` aborts rather than run. `aero_namelist` is written
    by ExoPlaSim's own Python API and this project had never touched it.

    Every value comes from the provenance file
    `aeolian/scripts/build_dust_source_fields.py` wrote beside the fields, so a
    run cannot be given a threshold, a Weibull shape or a roughness that
    disagrees with the source map it is applied to. Same arrangement as
    `enable_prescribed_dust` and for the same reason.

    THE REMOVAL SWITCHES RIDE HERE TOO. `ldepvel` and `lwetdep` come from
    exoplasim-3.4.2-aerosol-deposition.patch, they default off, and their
    coefficients live in the same `aeolian/config/dust.yaml`. They are read
    directly from that config rather than from the field provenance because they
    describe the atmosphere rather than the surface, and nothing about them
    depends on which terrain the fields were built on.

    Turning emission on also turns on the aerosol machinery it lives in:
    `L_AERO` gates `aero_main` and `l_source = 2` selects the dust case. Every
    run this project has made so far set `L_AERO = 0`, which is why `aerocore`
    had never executed and why seven defects survived in it.
    """
    if str(config["model"].get("dust_emission", "none")) == "none":
        return None
    field = surface_sra(config, min(DUST_EMISSION_SURFACE_CODES))
    prov_path = field.with_name(field.stem + "_provenance.json")
    if not prov_path.is_file():
        raise RuntimeError(
            f"{prov_path} is missing. Run "
            "aeolian/scripts/build_dust_source_fields.py; the namelist values "
            "live with the fields, not in the config.")
    prov = json.loads(prov_path.read_text(encoding="utf-8"))
    values = dict(prov["namelist_values"])

    dust_cfg_path = PROJECT_ROOT / "aeolian" / "config" / "dust.yaml"
    dust_cfg = yaml.safe_load(dust_cfg_path.read_text(encoding="utf-8"))
    removal = dust_cfg["removal"]

    # The removal terms are opt-in per run, because each one is a separate A/B
    # arm and A3 wants them switchable from the namelist alone.
    if bool(config["model"].get("dust_dry_deposition", False)):
        vd = config["model"].get("dust_deposition_velocity_m_s")
        if vd is None:
            raise RuntimeError(
                "model.dust_dry_deposition is on but "
                "model.dust_deposition_velocity_m_s is unset. vdaero is the "
                "NON-gravitational part of dry deposition and there is no "
                "defensible default; see the removal block of "
                "aeolian/config/dust.yaml and the deposition patch header.")
        values["LDEPVEL"] = 1
        values["VDAERO"] = float(vd)
    if bool(config["model"].get("dust_wet_scavenging", False)):
        values["LWETDEP"] = 1
        values["SCAVA"] = float(removal["scavenging_a"])
        values["SCAVB"] = float(removal["scavenging_b"])

    model._edit_namelist("plasim_namelist", "L_AERO", "1")
    model._edit_namelist("aero_namelist", "l_source", "2")

    # THE SHIPPED aero_namelist CANNOT BE READ, and nothing has noticed because
    # `aero_ini` runs only when `L_AERO > 0` and every run this project has made
    # set it to 0. It carries `aerofile = 0`, an unquoted integer for a
    # `character(len=80)`, which gfortran rejects with iostat 5010 -- and
    # `aero_ini` reads without an iostat, so the first run to turn the aerosol
    # on aborts in the namelist read before any of this is reached. Quoting it
    # is the fix and it has to happen here, because this is the first thing that
    # sets L_AERO at all.
    aerofile = PROJECT_ROOT / "exoplasim" / "data" / "dust" / "vesper_dust_aerosol.dat"
    if not aerofile.is_file():
        raise RuntimeError(
            f"{aerofile} is missing. Run exoplasim/scripts/dust_aerofile.py.")
    shutil.copyfile(aerofile, run_dir / aerofile.name)
    model._edit_namelist("aero_namelist", "aerofile", f"'{aerofile.name}'")

    # APART AND RHOP COME FROM THE FILE THE OPTICS WERE BUILT FOR, and they were
    # written by nothing at all. `aeromod.f90` declares apart = 50 nm and
    # rhop = 1000 kg/m3, a photochemical haze grain at water density, and
    # `aero_ini` validates the fourteen emission constants and not these two, so
    # an unwritten pair sits silently at the compiled values.
    #
    # WHAT THAT IS WORTH. `vels` is Stokes, v = 2 beta apart^2 ga (rhop - rhog) /
    # (9 mu), so the radius enters squared: (2.2068e-6/5e-8)^2 = 1948, density
    # 2.6, Cunningham 0.354 the other way, about 1790 net. The compiled default
    # settles the grain at roughly a micrometre a second where the intended one
    # falls at 0.21 cm/s, and sedimentation is the only removal term active at
    # ldepvel = 0 and lwetdep = 0, so the burden would instead be set by the
    # timestep-dependent bottom-layer scrub. `mmr2n` is off by the cube.
    #
    # The sidecar is the aerofile's own, beside the .dat: it is the file
    # `dust_aerofile.py` writes the burden-matched radius into, derived for THIS
    # planet's gravity, and it already carried a `namelist_values` block that
    # nothing opened. world-906.
    aero_prov = aerofile.with_name(aerofile.stem + ".provenance.json")
    if not aero_prov.is_file():
        raise RuntimeError(
            f"{aero_prov} is missing. APART and RHOP live with the aerofile "
            "they were derived for, not in the config; run "
            "exoplasim/scripts/dust_aerofile.py.")
    aero_values = dict(json.loads(aero_prov.read_text(encoding="utf-8"))
                       ["namelist_values"])
    for key in ("APART", "RHOP"):
        if float(aero_values.get(key, 0.0)) <= 0.0:
            raise RuntimeError(
                f"{aero_prov} carries no positive {key}. The tracer's radius "
                "and density set the settling velocity and the number density; "
                "aeromod's compiled defaults are a 50 nm haze grain at water "
                "density and are wrong by about three orders of magnitude in "
                "the settling velocity.")
        model._edit_namelist("aero_namelist", key, f"{float(aero_values[key])!r}")
    print(f"  dust grain: APART = {float(aero_values['APART']):.4e} m, "
          f"RHOP = {float(aero_values['RHOP']):g} kg/m3, from "
          f"{aero_prov.name}")

    # WHETHER THE EMITTED DUST IS RADIATIVELY ACTIVE, and it is a switch now
    # rather than a pin. Both reasons it was pinned off have been repaired:
    #
    #   * the shortwave. `l_aerorad = 1` put the emitted dust into the shortwave
    #     through `radmod`'s OWN `apart`, which `aero_ini` did not populate from
    #     the namelist, so the optical depth came out (50e-9/apart)**2 of intent
    #     -- 1/385 at this world's optical effective radius. `aeromod.f90`'s
    #     `aero_ini` now hands the value across as `rad_apart = apart`, inside
    #     the NROOT block because `radini` broadcasts it.
    #   * the longwave. `aeroqlw` had no writer anywhere in the project, so the
    #     interactive path had no thermal term structurally available, and
    #     shortwave-only dust is worth several kelvin of spurious cooling.
    #     `build_surface_dust.py` now derives it beside `DUSTQLW`; the two are
    #     the same ratio and, on this world, the same number, for the reason
    #     recorded there.
    #
    # It stays OFF by default, which is what every existing run has, and dust-13
    # is the decision. Turning it on needs the prescribed field's provenance for
    # the ratio, which is where the ratio is derived; there is no default for it
    # here and none in the Fortran either, because `radini` aborts on an
    # interactive aerosol with `aeroqlw` at zero rather than assume one.
    radiative = bool(config["model"].get("dust_emission_radiative", False))
    model._edit_namelist("aero_namelist", "l_aerorad", "1" if radiative else "0")
    if radiative:
        dust_field = surface_sra(config, sorted(DUST_SURFACE_CODES)[0])
        dust_prov = dust_field.with_name(dust_field.stem + "_provenance.json")
        if not dust_prov.is_file():
            raise RuntimeError(
                "model.dust_emission_radiative is on but "
                f"{dust_prov} is missing. AEROQLW is derived there, beside "
                "DUSTQLW and from the same optics; run "
                "exoplasim/scripts/build_surface_dust.py.")
        aeroqlw = float(json.loads(dust_prov.read_text(encoding="utf-8"))
                        ["namelist_values"].get("AEROQLW", 0.0))
        if aeroqlw <= 0.0:
            raise RuntimeError(
                f"{dust_prov} carries no AEROQLW. Shortwave-only dust is worse "
                "than no dust; see notes/dust.md and radini's own refusal.")
        # radmod_nl, NOT aero_nl, so it stays out of `values`: everything in
        # that dict is written into aero_namelist by the loop below.
        model._edit_namelist("radmod_namelist", "AEROQLW", f"{aeroqlw}")

    for key, value in values.items():
        model._edit_namelist("aero_namelist", key, f"{value}")
    return {
        "z0_bracket_end": prov["z0_bracket_end"],
        "aeolian_z0_m": prov["aeolian_z0_m"],
        "terrain_hash": prov["terrain_hash"],
        "field_sha256": prov["output_sha256"],
        "namelist_values": values,
        "aerosol_namelist_values": aero_values,
        "radiatively_active": radiative,
        "aeroqlw": aeroqlw if radiative else None,
        "source_cells": prov["field_statistics"]["srcw_cells_nonzero"],
    }


def stage_surface_extras(run_dir: Path, config: dict) -> list[int]:
    """Copy the surface fields we generate into the run directory.

    Has to run after `configure()`, which does `rm workdir/*.sra` whenever a
    landmap or topomap is given and then writes back only those two
    (`__init__.py:2938`). Anything staged before that call is silently deleted.
    """
    staged = []
    for code in sorted(intended_surface_codes(config) - BASE_SURFACE_CODES):
        src = surface_sra(config, code)
        if not src.is_file():
            if code in SOIL_WATER_SURFACE_CODES | SOIL_WATER_SPLIT_SURFACE_CODES:
                builder, setting = ("build_surface_soil_water.py",
                                    "model.soil_water_source")
            elif code in DUST_SURFACE_CODES:
                builder, setting = ("build_surface_dust.py", "model.dust_source")
            elif code in DUST_EMISSION_SURFACE_CODES:
                builder, setting = ("build_dust_source_fields.py",
                                    "model.dust_emission")
            elif code in ROUGHNESS_SURFACE_CODES:
                builder, setting = ("build_surface_roughness.py",
                                    "model.roughness_source")
            else:
                builder, setting = ("build_surface_albedo.py",
                                    "model.land_albedo_source")
            fallback = ("none"
                        if code in DUST_SURFACE_CODES | DUST_EMISSION_SURFACE_CODES
                        else "uniform")
            raise RuntimeError(
                f"{src} is missing. Run {builder}, or set {setting}: {fallback} "
                f"to accept ExoPlaSim's namelist default."
            )
        shutil.copyfile(src, run_dir / f"N{int(config['model']['latitudes']):03d}"
                                       f"_surf_{code:04d}.sra")
        staged.append(code)
    return staged

# The uniform values those fallbacks take, from plasim/src/landmod.f90 preset
# block and array declarations. Recorded in the manifest so a run says what its
# land surface actually was rather than leaving it to be inferred.
UNIFORM_LAND_DEFAULTS = {
    "dz0clim_roughness_m": 2.0,
    "dwmax_field_capacity_m": "wsmax (Earth value)",
    "dalbcl_background_albedo": 0.22,
    "dforest_fraction": 0.5,
    "dglac_glacier_mask": 0.0,
}


def surface_field_report(run_dir: Path, config: dict) -> dict:
    """Record which surface fields are set from file and which are uniform.

    `configure()` runs `rm workdir/*.sra` whenever a landmap or topomap is
    given (__init__.py:2938) and then writes only those two. So every other
    surface field is absent, and surfmod reads each one with
    `inquire(exist=)` and skips it silently when it is (surfmod.f90:125).

    That is not a failure. landmod presets the fields to uniform namelist
    values before reading, so the fallback is 2.0 m roughness, 0.22 albedo,
    Earth field capacity, 0.5 forest fraction and no glaciers, not zeros. For a
    custom world that is a defensible baseline and better than imprinting
    Earth's geography. It just has to be a declared choice rather than an
    accident, which is what `model.uniform_land_surface` is for.

    Because the wipe happens regardless of resolution, the shipped-data limit
    at T42 does not constrain us: a T85 run would default identically.
    """
    nlat = int(config["model"]["latitudes"])
    present = set()
    for f in run_dir.glob(f"N{nlat:03d}_surf_*.sra"):
        try:
            present.add(int(f.stem.rsplit("_", 1)[1]))
        except (IndexError, ValueError):
            continue

    missing = sorted(intended_surface_codes(config) - present)
    if missing:
        raise RuntimeError(
            f"{run_dir} is missing surface fields we supply ourselves: "
            f"{missing}. surfmod skips these silently, so the run would use "
            "flat terrain or an all-ocean mask without saying so."
        )

    declared = bool(config["model"].get("uniform_land_surface", False))
    defaulted = sorted(present - intended_surface_codes(config))
    if not declared:
        raise RuntimeError(
            "Every surface field other than topography and the land mask will "
            "fall back to a uniform namelist default, because configure() "
            "clears them. Set model.uniform_land_surface: true in "
            "config/planet.yaml to declare that this is intended."
        )
    return {
        "from_file": sorted(present),
        "unexpected_extra_files": defaulted,
        "uniform_defaults_declared": declared,
        "uniform_values": UNIFORM_LAND_DEFAULTS,
        "note": ("Fields absent from the run directory are preset to uniform "
                 "namelist values by landmod and never read from file."),
    }


def physical_fingerprint(config: dict, flux_ratio: float) -> dict:
    """Everything about a run that makes it a different climate.

    Recorded in the manifest and indexed by `index_runs.py` so a run is findable
    by what it is. It is deliberately NOT the run's identity -- see `run_id`.
    """
    p = config["planet"]
    a = config["atmosphere"]
    m = config["model"]
    return {
        "resolution": str(m["resolution"]),
        "layers": int(m["layers"]),
        "ranks": int(m["ncpus"]),
        "precision_bytes": int(m["precision_bytes"]),
        "flux_ratio": round(float(flux_ratio), 6),
        "co2_ppm": round(1e6 * float(a["pCO2_bar"]), 3),
        "ch4_ppmv": round(1e6 * float(a.get("pCH4_bar", 0.0)), 4),
        "n2o_ppmv": round(1e6 * float(a.get("pN2O_bar", 0.0)), 4),
        "rotation_hours": float(p["rotation_hours"]),
        "obliquity_degrees": float(p["obliquity_degrees"]),
        "eccentricity": float(p["eccentricity"]),
        "glaciers": bool(config["surface"].get("glaciers", {}).get("enabled")),
        "stellar_spectrum": config.get("radiation", {}).get("stellar_spectrum"),
        "geography": geography_tag(config),
    }


def run_id(config: dict, flux_ratio: float) -> str:
    """A UUID. Not derived from anything.

    This used to spell out the physical parameters in the directory name, so
    that two different climates could not share a directory, which matters
    because ExoPlaSim's `finalize()` picks output as the last glob match.

    **A derived identifier can only separate runs along the dimensions it
    encodes, and that is not a property you can maintain.** The ozone band-weight
    patch changed the physics and moved nothing in the encoded set, so the
    pre-patch and post-patch runs at the same flux computed the same identifier
    and landed in the same directory. Adding the patch stack to the name would
    have fixed that instance and not the next one; the encoded set is a list of
    everything someone has thought of so far.

    A UUID is unique unconditionally, so no run can ever collide with another
    regardless of what changed between them, including things nothing here
    models. What the run *was* belongs in the manifest, which is written anyway,
    and in `exoplasim/runs/INDEX.json`, which is generated from the manifests.

    The cost is that a continuation cannot recompute the directory name and must
    be given it. That is the safer direction: recomputation could silently
    resolve to a different run, and being handed an id cannot.
    """
    return "run_" + uuid.uuid4().hex[:12]


def validate_outputs(run_dir: Path) -> dict:
    regular = run_dir / "MOST.00000.nc"
    snapshot = run_dir / "snapshots" / "MOST_SNAP.00000.nc"
    required_regular = {
        "ts", "tas", "pr", "ua", "va", "snd", "sic",
        "sit", "mrso", "mrro", "evap", "ntr", "hfns", "lsm", "sg", "nu",
    }
    required_snapshot = {
        "ts", "tas", "pr", "ua", "va", "snd", "sic",
        "sit", "mrso", "mrro", "evap", "ntr", "hfns", "lsm", "sg", "nu",
    }
    result = {}
    for label, path, required in [
        ("regular", regular, required_regular),
        ("snapshot", snapshot, required_snapshot),
    ]:
        if not path.is_file():
            raise RuntimeError(f"Missing {label} output: {path}")
        with Dataset(path) as dataset:
            present = set(dataset.variables)
            missing = sorted(required - present)
            if missing:
                raise RuntimeError(f"{path} is missing required fields: {missing}")
            nonfinite = [
                name for name in sorted(required)
                if not np.isfinite(dataset.variables[name][:]).all()
            ]
            if nonfinite:
                raise RuntimeError(f"{path} has non-finite required fields: {nonfinite}")
            result[label] = {
                "path": str(path),
                "sha256": file_sha256(path),
                "time_samples": len(dataset.dimensions["time"]),
                "required_fields_present": sorted(required),
            }
    return result


def read_sra(path: Path, expected_code: int, nlat: int, nlon: int) -> np.ndarray:
    with path.open("r", encoding="ascii") as handle:
        header = np.fromstring(handle.readline(), sep=" ", dtype=np.int64)
        values = np.loadtxt(handle).reshape(-1)
    expected = np.array([expected_code, 0, 20260811, 0, nlon, nlat, 0, 0])
    if not np.array_equal(header, expected):
        raise ValueError(f"Unexpected SRA header in {path}: {header.tolist()}")
    if values.size != nlat * nlon or not np.isfinite(values).all():
        raise ValueError(f"Invalid SRA payload in {path}")
    return values.reshape(nlat, nlon)


def expected_namelist_keys(config: dict) -> dict:
    """What the CONFIG says the run's namelists must contain, {file: {KEY: value}}.

    Derived from `config/planet.yaml` and from nothing else -- deliberately not
    from the staging code, because a check built out of the staging lists tests
    only that they agree with themselves. CONS-9.

    A value is either a number or a SEQUENCE of numbers, and both go through one
    comparison: the staged text is expanded into its elements, Fortran's
    `n*value` replication included, and compared element by element at 1e-9.
    LENGTH is half the check and is what tells a replicated array apart from a
    scalar -- a namelist scalar assigned to an array sets element 1 only, so a
    one-element ARRAY where the config wants NLEV of them is the world-720
    failure and must not pass.
    """
    m = config["model"]
    want: dict = {"radmod_namelist": {}, "icemod_namelist": {}, "plasim_namelist": {},
                  "planet_namelist": {}, "landmod_namelist": {},
                  "glacier_namelist": {}, "oceanmod_namelist": {},
                  "rainmod_namelist": {}, "hurricane_namelist": {},
                  "fluxmod_namelist": {}}
    # CLIM-54 and CLIM-55, unconditional and not derived from a config key.
    # `declare_storm_diagnostics` says why both halves of hurricanemod are off;
    # they are checked here so that off is a state this project asserts rather
    # than the shipped namelist's default surviving by luck. Every other entry
    # in this function answers to `config/planet.yaml`, and these two do not
    # because a diagnostic whose indices are fitted to Earth is not a setting.
    want["hurricane_namelist"]["NSTORMDIAG"] = 0.0
    want["hurricane_namelist"]["HC_CAPTURE"] = 0.0
    for key, name, default in SHORTWAVE_GAS_KEYS:
        v = m.get(key)
        if v is not None and float(v) != default:
            want["radmod_namelist"][name] = float(v)
    o3 = m.get("ozone_scale")
    if o3 is not None and float(o3) != 1.0:
        want["radmod_namelist"]["O3SCALE"] = float(o3)
    # world-qvu. Written unconditionally by configure_otherargs, so they are
    # checked unconditionally: an absent key here means a continuation dropped
    # them and the modelled surface reverted to radmod.f90's compiled literal.
    want["radmod_namelist"]["ELWLAND"] = float(
        m.get("land_longwave_emissivity", 1.0))
    want["radmod_namelist"]["ELWSEA"] = float(
        m.get("sea_longwave_emissivity", 0.98))
    # world-ayx, and unconditional for the same reason.
    want["radmod_namelist"]["BO3"] = float(m.get("ozone_height_m", 20000.0))
    want["radmod_namelist"]["CO3"] = float(m.get("ozone_spread_m", 5000.0))
    # world-nfh, the per-band canopy albedo, and PHYS-11's cloud absorption
    # scale. Both are written unconditionally by `configure_otherargs`, so both
    # are checked unconditionally: ALBFOREST is per band and dropping it reverts
    # the modelled canopy to landmod's compiled Earth-Sun broadband endmember in
    # both bands, and dropping CLOUDABS reverts the modelled cloud co-albedo to
    # the solar weighting Stephens's tables are printed at. ALBFOREST is the one
    # this gate could not hold until the comparison carried sequences;
    # world-2wd.
    #
    # Rounded to the six significant digits `configure_otherargs` writes, for
    # AKAP's reason below: the comparison is against the value the namelist can
    # hold and not against one more digit than it carries.
    want["landmod_namelist"]["ALBFOREST"] = [
        float(f"{float(v):.6g}") for v in m["vegetation_albedo_bands"]]
    cloud_scale = float(m.get("cloud_absorption_scale", 1.0))
    want["radmod_namelist"]["CLOUDABS"] = float(f"{cloud_scale:.6g}")
    # clim-68. NSWRCL selects the shortwave cloud scheme, and CLOUDABS acts
    # only in the computed branch it selects at 1. Checking it is what makes
    # the line above mean anything: a segment that reached the prescribed
    # branch would carry a scale no code reads, and `acl2`, which that branch
    # does read, is left at radmod's Earth tuning deliberately.
    want["radmod_namelist"]["NSWRCL"] = 1.0
    # world-35en. NCLOUDS gates the branch NSWRCL selects within, so checking
    # NSWRCL alone leaves the shortwave's cloud treatment half-stated.
    want["radmod_namelist"]["NCLOUDS"] = 1.0
    # CLIM-16, oceanmod_nl, and unconditional for the same reason. NLEV_OCE is
    # 1 so HDIFFK is a one-element array, which is why one element is what this
    # asks for.
    ocean = config.get("ocean", {})
    want["oceanmod_namelist"]["NHDIFF"] = float(
        int(bool(ocean.get("horizontal_diffusion", False))))
    want["oceanmod_namelist"]["HDIFFK"] = [
        float(f"{float(ocean.get('horizontal_diffusivity_m2_s', 1.0e3)):.6g}")]
    # world-cwc and world-qpe, written unconditionally by configure_otherargs
    # for the same reason. DSMAX is the one CLIM-17 already lost once: a key
    # whose whole purpose is to be lifted reverts to the compiled 5 m on any
    # continuation that drops it, and the glacier module then cannot build an
    # ice sheet at all.
    glaciers = config.get("surface", {}).get("glaciers", {})
    want["landmod_namelist"]["DSMAX"] = float(
        glaciers.get("max_snow_depth_m", -1.0))
    want["glacier_namelist"]["GLACPERSIST"] = float(
        glaciers.get("persistence_orbits", 1.0))
    # world-n1nu, rainmod_nl, and unconditional for the same reason: the
    # control arm of the cloud water bracket is the compiled value, so a
    # continuation that dropped CLWREF would produce a control that agrees with
    # its own arm by accident rather than by construction.
    clwref = float(m.get("cloud_water_reference_kg_m3",
                         rainmod_default("clwref")))
    want["rainmod_namelist"]["CLWREF"] = float(f"{clwref:.6g}")
    # world-trs3 and world-o12h, rainmod_nl, unconditional on the same
    # argument: each selects between a derived form and a literal, so a
    # continuation that dropped one would silently return the segment to the
    # derived branch and no artifact would say it had.
    gamma = float(m.get("precip_reevaporation_gamma", rainmod_default("gamma")))
    want["rainmod_namelist"]["GAMMA"] = float(f"{gamma:.6g}")
    rcritwidth = float(m.get("cloud_fraction_subgrid_width",
                             rainmod_default("rcritwidth")))
    want["rainmod_namelist"]["RCRITWIDTH"] = float(f"{rcritwidth:.6g}")
    # world-80ia, fluxmod_nl, unconditional on the same argument: the
    # asymptotic mixing length selects between a rotation-derived value and a
    # declared literal, so a continuation that dropped it would return the
    # segment to the derived branch with nothing in the run directory saying so.
    lamm = float(m.get("asymptotic_mixing_length_m",
                       fluxmod_default("vdiff_lamm")))
    want["fluxmod_namelist"]["VDIFF_LAMM"] = float(f"{lamm:.6g}")
    # world-py6p, landmod_nl. THE LAND COLUMN, all eight, unconditionally.
    # `land_water_column` is shared with the staging side the way
    # `freezing_point_k` is: it maps the CONFIG to a value and knows nothing
    # about what gets written, and the list of keys below is stated here
    # independently, which is the half that must not come from the staging
    # code. The two per-layer keys are asked for at the DECLARED layer count,
    # so a config that declares three layers and a namelist carrying one fails
    # on its length -- world-720's failure, on the land column.
    lwc = land_water_column(config)
    for _key in ("NLANDWCOL", "NLSOILW", "NLANDWDRAIN", "NLANDWPHASE",
                 "NRHSEXP"):
        want["landmod_namelist"][_key] = float(lwc[_key])
    want["landmod_namelist"]["DRHSLOW"] = float(f"{lwc['DRHSLOW']:.6g}")
    for _key in ("DSOILWF", "DSOILWZ"):
        want["landmod_namelist"][_key] = [
            float(f"{v:.6g}") for v in lwc[_key]]
    salinity = config.get("ocean", {}).get("salinity_psu")
    if salinity is not None:
        celsius = freezing_point_k(salinity) - 273.15
        want["icemod_namelist"]["TFREEZE"] = round(freezing_point_k(salinity), 4)
        want["icemod_namelist"]["CRHOS"] = round(
            seawater_density(salinity, celsius), 2)
        want["icemod_namelist"]["CPS"] = round(
            seawater_heat_capacity(salinity, celsius), 2)
    fusion = config.get("ocean", {}).get("sea_ice_fusion_j_kg")
    if fusion is not None:
        want["icemod_namelist"]["CLFI"] = round(float(fusion), 1)
    cold = config.get("ocean", {}).get("cold_start")
    if cold is not None:
        want["icemod_namelist"]["TSST_EQ"] = round(
            float(cold["sst_equator_k"]), 4)
        want["icemod_namelist"]["TSST_POL"] = round(
            float(cold["sst_pole_k"]), 4)
        want["icemod_namelist"]["HICE_INI"] = round(
            float(cold["sea_ice_thickness_m"]), 4)
    ice_max = config.get("surface", {}).get("sea_ice_max_thickness_m")
    if ice_max is not None:
        want["icemod_namelist"]["XMAXD"] = round(float(ice_max), 4)
    lead = config.get("surface", {}).get("sea_ice_lead_closing_m")
    if lead is not None:
        want["icemod_namelist"]["HLEAD"] = round(float(lead), 4)
    cryo = config.get("surface", {}).get("cryosphere")
    if cryo is not None:
        for _key, _name, _dp in (("sea_ice_density_kg_m3", "CRHOI", 2),
                                 ("sea_ice_specific_heat_j_kg_k", "CPI", 2),
                                 ("sea_ice_conductivity_w_m_k", "CKAPI", 4),
                                 ("snow_fusion_j_kg", "CLFSN", 2)):
            want["icemod_namelist"][_name] = round(float(cryo[_key]), _dp)
    if energy_diagnostics_enabled(config):
        want["plasim_namelist"]["NENERGY"] = float(energy_diagnostics_level(config))
    if m.get("conversion_time_level", False):
        want["plasim_namelist"]["NCONVTIME"] = 1.0
    if m.get("dealias_conversion", False):
        want["plasim_namelist"]["NDEALIAS"] = 1.0
    # THE FIXER IS A DECLARED SWITCH AND BELONGS IN THE CHECK. It is written on
    # both branches -- 1 when on and 0 when off -- so its absence from a staged
    # namelist means a continuation dropped it, which is exactly what this check
    # exists to catch. failure-modes class 22.
    if m.get("energy_fixer") is not None:
        want["plasim_namelist"]["NENERGYFIX"] = 1.0 if m["energy_fixer"] else 0.0
    # world-6qee, and unconditional for NSTORMDIAG's reason rather than
    # derived from a config key: no orography scaling is an assertion about
    # this world, not a knob it sets. `configure_otherargs` writes it on every
    # prepare and every resume, so an absent key means a continuation dropped
    # it and the model reverted to p_earth.f90's compiled literal.
    want["planet_namelist"]["OROSCALE"] = 1.0
    if m.get("robert_filter") is not None:
        want["planet_namelist"]["PNU"] = float(m["robert_filter"])
        if m.get("energy_diagnostics_3d", False):
            want["plasim_namelist"]["NENER3D"] = 1.0
    # The scalar half of `declare_dry_constants`. Each of these reverts to a
    # p_earth.f90 or plasimmod.f90 Earth default if a continuation drops it, and
    # TGR reverting moves the realised mean surface pressure, so the class-22
    # failure here is a change in the mean state rather than a lost switch. T0
    # and TFRC are per-level lists and cannot go through this float comparison.
    profile = m.get("cold_start_profile")
    if profile is not None:
        import lapse
        gas_constant, cp = lapse.gas_properties(config)
        # Formatted the way `declare_dry_constants` writes it, so the comparison
        # is against the value the namelist can hold and not against one more
        # digit than it carries.
        want["planet_namelist"]["AKAP"] = float(f"{gas_constant / cp:.8g}")
        want["planet_namelist"]["ALR"] = float(profile["lapse_rate_k_per_m"])
        want["plasim_namelist"]["TGR"] = float(profile["surface_temperature_k"])
        want["plasim_namelist"]["DTROP"] = float(profile["tropopause_height_m"])
    # THE FILTER AND THE HYPERDIFFUSION, the two sets world-8bs found missing
    # from every continuation. They fail in opposite ways and both are silent.
    # `configure()` writes FILTERKAPPA and NFILTEREXP unconditionally from its
    # own defaults, so a dropped key is a WRONG value; it does not write NDEL,
    # NHDIFF or TDISS* at all, so a dropped key is an ABSENT one and the model
    # falls back to `readnl`'s compiled T21/T42 branch. Checking only the config
    # keys catches the first; checking presence catches the second.
    # world-a05: NEQSIG and PTOP came from a library subclass default. Now they
    # are config, and a continuation that dropped them would put the model top
    # somewhere else without saying so. PTOP is written in Pa from hPa.
    if m.get("vertical_grid") is not None:
        want["plasim_namelist"]["NEQSIG"] = float(m["vertical_grid"])
    if m.get("model_top_hpa") is not None:
        want["plasim_namelist"]["PTOP"] = float(m["model_top_hpa"]) * 100.0
    if m.get("filter_kappa") is not None:
        want["plasim_namelist"]["FILTERKAPPA"] = float(m["filter_kappa"])
    if m.get("filter_power") is not None:
        want["plasim_namelist"]["NFILTEREXP"] = float(m["filter_power"])
    hd = m.get("hyperdiffusion")
    if hd:
        rung = str(m["resolution"]).upper()
        tau = hd["timescales_days"][rung]
        ntru = int(rung.lstrip("Tt"))
        want["plasim_namelist"]["NHDIFF"] = float(
            round(float(hd["cutoff_fraction"]) * ntru))
        # PER-LEVEL ARRAYS, asked for at their full length. A scalar here would
        # pass a value comparison while levels 2..NLEV kept `readnl`'s values,
        # which is world-720's failure; it fails a LENGTH comparison, which is
        # exactly what this must not miss.
        layers = int(m["layers"])
        want["plasim_namelist"]["NDEL"] = [float(int(hd["order_alpha"]))] * layers
        for key, field in (("TDISSD", "divergence"), ("TDISSZ", "vorticity"),
                           ("TDISST", "temperature"), ("TDISSQ", "humidity")):
            want["plasim_namelist"][key] = [float(tau[field])] * layers
    return {f: keys for f, keys in want.items() if keys}


def verify_staged_namelists(run_dir: Path, config: dict) -> dict:
    """Every config key that departs from the model default reached the namelist.

    CONS-9. Reads the RUN, which is the only artifact that records what was
    actually integrated -- `docs/src/practice/failure-modes.md` class 12 said so and class 22
    is what happens when nobody checks: `h2o_sw_level` was declared at 1.127,
    staged by the prepare script, dropped by every continuation, and no guard
    noticed for a whole run. `config_drift` compared the config against the
    manifest and both said 1.127; this compares the config against the FILE.

    Raises rather than warns. A run whose radiation does not match its own
    configuration is not a cheaper run, it is a different world.
    """
    checked, wrong = {}, []
    for fname, keys in expected_namelist_keys(config).items():
        path = run_dir / fname
        for key, want in keys.items():
            try:
                raw = namelist_value(path, key)
            except (KeyError, FileNotFoundError):
                wrong.append(f"{key} absent from {fname}, config wants {want}")
                continue
            wanted = [float(v) for v in want] if isinstance(want, (list, tuple)) \
                else [float(want)]
            try:
                got = namelist_elements(raw)
            except ValueError:
                wrong.append(f"{key} in {fname} is not a number")
                continue
            if len(got) != len(wanted):
                # LENGTH FIRST. A scalar staged where the config wants a
                # per-level array sets element 1 and leaves levels 2..NLEV at
                # `readnl`'s compiled values, which reads as correct at every
                # element the comparison would otherwise look at. world-720.
                wrong.append(
                    f"{key} in {fname} has {len(got)} element(s), config wants "
                    f"{len(wanted)}: {raw!r}")
            else:
                for i, (g, w) in enumerate(zip(got, wanted), start=1):
                    if abs(g - w) > 1e-9:
                        where = f"{key} in {fname}" if len(wanted) == 1 \
                            else f"element {i} of {key} in {fname}"
                        wrong.append(f"{where} is {g:g}, config wants {w:g}")
            checked[f"{key}@{fname}"] = got[0] if len(got) == 1 else got
    if wrong:
        raise SystemExit(
            "the staged namelists do not match config/planet.yaml:\n  "
            + "\n  ".join(wrong)
            + "\nThe namelist in the run directory is what the model integrates. "
              "Fix the staging rather than the check; see docs/src/practice/failure-modes.md "
              "class 22.")
    return checked


def namelist_elements(raw: str) -> list[float]:
    """A staged namelist value as the list of numbers the model will read.

    Fortran's `n*value` replication is expanded, because a replicated array and
    the same array written out element by element are the same value to the
    model and must be the same value here. Raises ValueError on anything that
    is not a number, which is what puts a logical or a string beyond this gate.
    """
    elements: list[float] = []
    for token in raw.replace(",", " ").split():
        count, star, value = token.partition("*")
        if star:
            elements.extend([float(value)] * int(count))
        else:
            elements.append(float(count))
    if not elements:
        raise ValueError(f"no numbers in {raw!r}")
    return elements


def namelist_value(path: Path, key: str) -> str:
    for line in path.read_text(encoding="ascii").splitlines():
        stripped = line.strip().upper()
        if stripped.startswith(key.upper() + " ") or stripped.startswith(key.upper() + "="):
            return line.split("=", 1)[1].strip()
    raise KeyError(f"{key} not found in {path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--flux-ratio", type=float, default=None)
    # THE DEFAULT IS THE CHEAP REGIME, because the orbits this script writes are
    # a spin-up: a run is prepared cold or off a seed and integrates toward
    # equilibrium, and nothing should be building a climatology from the first
    # block. The clean regime is bought per segment, at the end, for the orbits
    # an analysis actually reads -- `continue_exoplasim.py --purpose
    # post_equilibrium_climatology`, which refuses low I/O outright.
    #
    # This was `disable_low_io()` with no way to say otherwise, so every orbit
    # of every run paid for instantaneous samples and the spin-up ones threw
    # them away. Safe to flip because the bin-0 defect is fixed and verified
    # (`exoplasim/notes/first-output-bin.md`) and because `nlowio` is now
    # broadcast, without which NLOWIO = 0 deadlocks a patched model
    # (`notes/audits/nlowio-collective-deadlock.md`).
    #
    # What low I/O costs is not accuracy but INFORMATION: it writes interval
    # accumulations where the clean regime writes instantaneous samples, and an
    # accumulation cannot be undone. Anything reading variance, extremes or
    # single records needs --clean-io.
    # BENCHMARK OVERRIDES. Both change the machine's layout and neither changes
    # the physics, so they belong on the command line rather than in
    # config/planet.yaml -- a config edit would move config_sha256 and make
    # every existing run unresumable, for a question about core placement.
    # `exoplasim/notes/rank-layout-benchmark.md` is what they are for.
    parser.add_argument(
        "--ncpus", type=int, default=None,
        help="override model.ncpus for this run. Selects a different compiled "
             "binary, since ExoPlaSim builds one per (resolution, layers, "
             "threads) configuration; NLAT must divide by it")
    # DIAGNOSTIC. The model writes an output record every `nafter` timesteps and
    # defaults to one per EARTH day, 32 steps of 45 minutes. This planet's day is
    # 30 hours, which is exactly 40 steps, so the default samples the diurnal
    # cycle at a 32/40 = 4/5 ratio: consecutive records advance the phase by
    # four fifths of a day and land on FIVE PHASES, forever. The regular output's
    # annual mean therefore carries a fixed aliasing bias, while the ocean and
    # ice streams accumulate every step and do not. CLIM-11.
    #
    # Choose a value coprime with the 40-step day to spread the sampling over
    # all forty phases -- 37 or 39 -- and the bias should collapse if that is
    # what it is. Diagnostic only: it changes what is written, not what is
    # integrated, but it moves the record count and so the climatology.
    parser.add_argument(
        "--writes-per-day", type=int, default=None,
        help="NWPD, output records per day. The model derives its write "
             "interval as nafter = mtspd / nwpd, so this is the control and "
             "NAFTER is not a namelist key -- setting that one aborts the run "
             "with `Cannot match namelist object name nafter`. Raising it "
             "samples more diurnal phases")
    # EFOR-2. On the command line for the same reason as the overrides above,
    # and it is a stronger reason here: the stream cannot change a result at
    # all. `declare_ecological_stream` has the argument.
    parser.add_argument(
        "--ecological-stream", action="store_true",
        help="NECO = 1, writing the biosphere's own output stream to "
             "plasim_eco alongside the climate one. Off by default. It reads "
             "state that already exists and writes to a unit of its own, so "
             "the restart is byte-identical either way")
    parser.add_argument(
        "--eco-interval-steps", type=int, default=None,
        help="NECOSTEP, timesteps per ecological interval. Left unset the "
             "model uses mtspd, one absolute 24-hour day exactly, which is "
             "the step LPJ-GUESS integrates on. Set it only for a run whose "
             "consumer integrates on something else")
    parser.add_argument(
        "--clean-io", dest="low_io", action="store_false", default=True,
        help="run this block at NLOWIO = 0, writing instantaneous samples "
             "rather than interval accumulations. The default is the cheap "
             "regime, because a prepared run is a spin-up; pass this when the "
             "block itself has to be read as data")
    parser.add_argument(
        "--run-years", type=int, default=0,
        help="Run this many model orbits after preparation",
    )
    parser.add_argument(
        "--force-prepare", action="store_true",
        help="Allow re-preparing an existing run with no climate outputs",
    )
    parser.add_argument(
        "--superseded-surface-ok", action="store_true",
        help="allow --restart-from when the staged surface fields have moved "
             "since the donor run. The new fields are DISCARDED: a restart "
             "reads albedo, roughness and soil water from itself. Only valid "
             "for a paired A/B, where every arm inherits the same superseded "
             "surface and the difference is what is measured. Stamped on the "
             "run manifest so the run is self-labelling.")
    parser.add_argument(
        "--binary", type=Path, default=None,
        help="integrate this run with a named executable instead of the one "
             "the registry's naming resolves to. It must be an ARM: a build "
             "made by exoplasim/scripts/build_model.py with --no-publish, "
             "still in its own build directory, whose directory name states "
             "every input that makes it differ from the shipped model. The arm "
             "is copied into THIS run directory only; the published binary and "
             "binary_manifest.json are untouched. The run is stamped with the "
             "arm's build tag and with canonical_lineage_eligible = false, and "
             "continue_exoplasim.py will not take a "
             "--purpose post_equilibrium_climatology segment on it")
    parser.add_argument(
        "--restart-from", type=Path, default=None,
        help="Seed the initial state from an existing MOST_REST file instead of "
             "cold-starting. Only the spin-up path changes, not the equilibrium.",
    )
    args = parser.parse_args()
    config_path = args.config.resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    # Applied to the loaded config rather than threaded through, so the binary
    # name, the physical fingerprint and the model construction cannot disagree
    # about how many ranks this run has. `config_sha256` hashes the FILE, so
    # provenance still records the config as written, and the fingerprint
    # records the ranks that actually ran.
    if args.ncpus is not None:
        if args.ncpus < 1:
            raise ValueError("--ncpus must be positive")
        config["model"]["ncpus"] = int(args.ncpus)
    flux_ratio = float(
        config["orbit"]["baseline_flux_earth"]
        if args.flux_ratio is None else args.flux_ratio
    )
    if flux_ratio <= 0:
        raise ValueError("Flux ratio must be positive")
    derived = derive(config, flux_ratio)
    identifier = run_id(config, flux_ratio)
    run_dir = (RUNS / identifier).resolve()
    # Announce the id on its own line, in a machine-readable form. It is a UUID,
    # so a caller cannot reconstruct it, and the alternative -- picking the
    # newest or first run_* directory afterwards -- is selection by sort order,
    # which is the pattern this project has now been bitten by four times.
    print(f"RUN_ID={identifier}", flush=True)
    run_dir.parent.mkdir(parents=True, exist_ok=True)
    # Any of these left over from an earlier run makes finalize() and the
    # crash-tolerant rewind pick up the wrong world's state, so refuse the lot
    # rather than only the primary output, and do it regardless of --force-prepare.
    stale = []
    if run_dir.exists():
        for pattern in ("MOST.*", "MOST_REST.*", "MOST*DIAG*",
                        "snapshots/*", "highcadence/*"):
            stale.extend(sorted(run_dir.glob(pattern)))
    if stale:
        raise RuntimeError(
            f"{run_dir} already holds model artifacts ({len(stale)}, e.g. "
            f"{stale[0].name}); this command will not overwrite or implicitly "
            "resume them. Use a fresh directory or move the old run aside."
        )
    if run_dir.exists() and any(run_dir.iterdir()) and not args.force_prepare:
        raise RuntimeError(
            f"{run_dir} already exists and is nonempty; inspect it or pass --force-prepare"
        )

    nlat = int(config["model"]["latitudes"])
    nlon = int(config["model"]["longitudes"])
    landmap = surface_sra(config, 172).resolve()
    topomap = surface_sra(config, 129).resolve()
    land = read_sra(landmap, 172, nlat, nlon)
    topo = read_sra(topomap, 129, nlat, nlon)
    if not np.all(np.isin(land, [0.0, 1.0])):
        raise ValueError("Land SRA is not binary")
    if np.any(topo[land == 0.0] != 0.0):
        raise ValueError("Ocean geopotential must be zero")

    # Seeding the initial state from a nearby equilibrium is worth about 20 of
    # the 70 orbits a cold start spends climbing from its 269 K initial state.
    # It is an initial condition and nothing else: the equilibrium a run settles
    # to is set by its forcing, so this changes the path and not the answer. The
    # source is recorded in the manifest because the path is no longer a function
    # of the configuration alone.
    restart_seed = None
    superseded_surface = None
    conversion = None
    if args.restart_from is not None:
        restart_seed = args.restart_from.resolve()
        if not restart_seed.is_file():
            raise RuntimeError(f"--restart-from {restart_seed} does not exist")
        # A restart FREEZES every land surface boundary condition.
        #
        # landmod.f90's landini reads dwmax, dz0clim, dz0climo, dalbcl, dalbcl1
        # and dalbcl2 from the restart file when restart > 0, and from the
        # surface .sra files only on a cold start. So seeding from a restart
        # silently discards any surface field that has changed since -- soil
        # water, roughness and all three albedo bands.
        #
        # This was found by running exactly that: a baseline seeded from the
        # bootstrap, with a new soil-water map and lake-composited albedo, would
        # have reproduced the bootstrap and shown that lakes and soil water do
        # nothing. It only surfaced because it also happened to trap a SIGFPE in
        # landini. A silent null result is the failure mode this guard exists to
        # prevent.
        # CLIM-31. A restart carries the donor's PARTIAL ACCUMULATION and the
        # counters that normalise it, so a seeded run opens mid-window: it adds
        # its own steps to somebody else's sum and divides by its own count. The
        # damage is one output record and it is visible on a field that cannot
        # vary -- under low I/O the land mask reads 0.95341 in orbit 0 bin 0,
        # exactly 1.0 everywhere else, and a cold-started low-I/O run reads 1.0
        # throughout.
        #
        # Zeroed in a COPY, so the donor run is never modified: a restart is the
        # only record of where a run was, and editing one in place would make a
        # completed run unreproducible to fix a defect in the run seeded from
        # it. Every record it rewrites is one the model itself puts back to a
        # clean value at an interval boundary -- which is zero for most and a
        # large sentinel for the two running minima -- so the copy describes a
        # run sitting exactly at the start of an accumulation window.
        seeded_from = restart_seed

        # BEFORE the copy, and before the run directory exists. A guard that
        # fires after doing work leaves the work behind: this one used to run
        # after the seed was written and left a run directory holding nothing
        # but a MOST_REST.seed on every refusal.
        #
        # A CONVERTED restart is not a donor run's file and has no run
        # directory to carry a manifest. Its provenance is the conversion
        # report `convert_restart.py` wrote beside it, which names the template
        # its static surface records came from and the hash of every staged
        # surface file that template was cut against. That is exactly what the
        # guard below needs, so a converted state is checked rather than
        # exempted.
        conversion = load_conversion_report(seeded_from, restart_seed)
        src_manifest = seeded_from.parent / "run_manifest.json"
        if conversion is not None:
            reason = conversion_surface_reason(conversion, config)
            if reason and args.superseded_surface_ok:
                print(f"  OVERRIDDEN (--superseded-surface-ok): {reason}")
                superseded_surface = reason
                reason = None
            if reason:
                raise RuntimeError(
                    "--restart-from refused a converted restart: " + reason
                    + ". The converted state's static surface records come "
                    "from its target template, not from this run's staged "
                    ".sra files, so a mismatch means the run would integrate "
                    "a surface it did not stage. Rebuild the template against "
                    "the current surface with build_restart_template.py.")
            print(f"  conversion report: "
                  f"{conversion['source']['truncation']} -> "
                  f"{conversion['target_template']['truncation']}, donor "
                  f"{(conversion.get('source_manifest') or {}).get('run_id')}, "
                  "target surface matches this run's staged files")
        elif not src_manifest.is_file():
            raise RuntimeError(
                f"--restart-from refused: {seeded_from} has neither a "
                "run_manifest.json beside it nor a conversion report, so "
                "nothing says which surface it froze. A restart FREEZES soil "
                "water, roughness and albedo -- landmod reads them from the "
                "restart, not from the .sra -- and with no provenance that "
                "cannot be ruled out. Seed from a run directory, or convert "
                "with convert_restart.py, which writes the report this needs.")
        if src_manifest.is_file():
            src = json.loads(src_manifest.read_text(encoding="utf-8"))
            was = set((src.get("surface_fields") or {}).get("from_file") or [])
            now = intended_surface_codes(config)
            reason = None
            # The guard above is about fields landmod reads only on a cold
            # start. It does not apply to every surface field, and refusing a
            # restart for one it does not apply to costs a cold start for no
            # reason. `radini` reads the dust field with `mpsurfgp` outside any
            # `nrestart` test, where landmod's block is inside `if (nrestart ==
            # 0)`, so a restarted run picks it up. Adding or removing dust is
            # therefore a forcing change and an ordinary perturbation
            # experiment, which is exactly what seeding from an equilibrium is
            # for. The exemption is stated per code rather than assumed for the
            # class: anything else added here must be checked the same way,
            # in the Fortran and not from the name.
            transparent = DUST_SURFACE_CODES
            if (was - transparent) != (now - transparent):
                reason = (f"codes differ: {sorted(was)} then, {sorted(now)} now")
            else:
                # Same codes, possibly different content.
                old_h = src.get("surface_field_sha256") or {}
                changed = [c for c in sorted(now)
                           if surface_sra(config, c).is_file()
                           and old_h.get(str(c)) not in
                           (None, file_sha256(surface_sra(config, c)))]
                if changed:
                    reason = f"content changed for code(s) {changed}"
                elif not old_h:
                    reason = ("the source run recorded no surface field hashes, "
                              "so content changes cannot be ruled out")
            if reason and args.superseded_surface_ok:
                # Deliberately overridden. The guard protects the CANONICAL
                # chain: a spin-up seeded this way would silently discard the
                # newer surface and reproduce its parent. An A3 forcing A/B is
                # the one case where that does not matter, because every arm
                # inherits the SAME superseded surface and what is being
                # measured is the DIFFERENCE between arms, not the absolute
                # climate. The override is a flag rather than a deletion so it
                # cannot be taken by accident, and it is stamped on the manifest
                # below so the resulting run says what it is.
                print(f"  OVERRIDDEN (--superseded-surface-ok): {reason}. "
                      "The staged surface is DISCARDED and this run inherits "
                      "its parent's. Valid for a paired A/B, invalid for "
                      "anything whose absolute climate is read.")
                superseded_surface = reason
                reason = None
            if reason:
                raise RuntimeError(
                    "--restart-from refused: " + reason + ". A restart reads soil "
                    "water, roughness and albedo from ITSELF, not from the .sra "
                    "files -- landmod's landini takes dwmax, dz0clim and all "
                    "three dalbcl bands from the restart when restart > 0 -- so "
                    "those changes would be silently discarded and the run would "
                    "reproduce its parent. Cold-start instead.")

        # Every refusal above has now been passed, so the copy is made and the
        # run directory created. Zeroed in a COPY so the donor is never
        # modified: a restart is the only record of where a run was, and
        # editing one in place would make a completed run unreproducible to fix
        # a defect in the run seeded from it.
        restart_seed = run_dir / "MOST_REST.seed"
        run_dir.mkdir(parents=True, exist_ok=True)
        reset_info = reset_restart_accumulators.reset(seeded_from, restart_seed)
        print(f"seed {seeded_from.name} from {seeded_from.parent.name}: "
              f"zeroed {len(reset_info['zeroed'])} accumulator records, "
              f"{len(reset_info['sentinels'])} put back to a nonzero clean "
              f"value, {len(reset_info['never_reset_by_the_model'])} the model "
              "resets nowhere left as they were (CLIM-31)")

    atmosphere = config["atmosphere"]
    planet = config["planet"]
    star = config["star"]
    model_cfg = config["model"]
    surface = config["surface"]
    thread_stack = prepare_thread_stack(config)
    model = exo.Earthlike(
        resolution=model_cfg["resolution"],
        layers=int(model_cfg["layers"]),
        ncpus=int(model_cfg["ncpus"]),
        precision=int(model_cfg["precision_bytes"]),
        workdir=str(run_dir),
        modelname=identifier,
        outputtype=model_cfg["output_type"],
    )
    print(f"threads: {int(model_cfg['ncpus'])}, SHTns transform")
    timestep_minutes = declare_timestep(config)
    model.configure(
        restartfile=None if restart_seed is None else str(restart_seed),
        flux=derived["stellar_flux_w_m2"],
        startemp=float(star["effective_temperature_k"]),
        starspec=stellar_spectrum_path(config),
        starradius=derived["stellar_radius_solar"],
        pN2=float(atmosphere["pN2_bar"]),
        pO2=float(atmosphere["pO2_bar"]),
        pAr=float(atmosphere["pAr_bar"]),
        pCO2=float(atmosphere["pCO2_bar"]),
        rotationperiod=derived["rotation_days"],
        year=derived["orbital_year_earth_days"],
        gravity=derived["gravity_m_s2"],
        radius=float(planet["radius_earth"]),
        eccentricity=float(planet["eccentricity"]),
        obliquity=float(planet["obliquity_degrees"]),
        lonvernaleq=float(config["orbit"]["longitude_vernal_equinox_degrees"]),
        fixedorbit=True,
        keplerian=True,
        meananomaly0=0.0,
        seaice=bool(surface["sea_ice"]),
        glaciers={
            "toggle": bool(surface.get("glaciers", {}).get("enabled", False)),
            "mindepth": float(surface.get("glaciers", {}).get("min_snow_depth_m", 2.0)),
            "initialh": float(surface.get("glaciers", {}).get("initial_height_m", -1.0)),
        },
        # world-cwc. landmod's compiled DSMAX is 5 m water equivalent, which is
        # PlaSim's Earth seasonal-snow ceiling and caps the same accumulator
        # glaciermod builds ice orography from. With the glacier module on, a
        # cell that should glaciate instead manufactures meltwater at the
        # snowfall rate forever.
        maxsnow=derived["max_snow_depth_m"],
        ozone=bool(atmosphere["ozone"]),
        mldepth=float(surface["mixed_layer_depth_m"]),
        twobandalbedo=bool(config["radiation"]["two_band_albedo"]),
        # DECLARED, not inherited. Earthlike.configure supplies vtype=4 and
        # modeltop=50.0 from its own signature, so NEQSIG and PTOP reached the
        # model from a library subclass default that no document here named.
        # Passed explicitly so the vertical grid and the model top are config,
        # and verified below like every other config-set key. world-a05.
        vtype=int(model_cfg["vertical_grid"]),
        modeltop=float(model_cfg["model_top_hpa"]),
        timestep=timestep_minutes,
        physicsfilter=model_cfg["physics_filter"],
        filterkappa=float(model_cfg["filter_kappa"]),
        filterpower=int(model_cfg["filter_power"]),
        landmap=str(landmap),
        topomap=str(topomap),
        runsteps=int(derived["runsteps_per_orbit"]),
        snapshots=int(derived["snapshot_interval_steps"]),
        otherargs=configure_otherargs(derived),
    )
    # Shipped output lists omit orbital phase and several hydrology fields.
    # This private helper is stable in the pinned release and edits those
    # postprocessor-code lists only.
    regular_codes = list(REGULAR_CODES)
    if energy_diagnostics_enabled(config):
        regular_codes = regular_codes + ENERGY_DIAGNOSTIC_CODES
        if config["model"].get("energy_diagnostics_3d", False):
            regular_codes = regular_codes + ENERGY_3D_CODES
    refuse_eco_codes(regular_codes, "REGULAR_CODES")
    refuse_eco_codes(SNAPSHOT_CODES, "SNAPSHOT_CODES")
    model._add_postcodes("example.nl", regular_codes)
    model._add_postcodes("snapshot.nl", SNAPSHOT_CODES)
    model.cfgpostprocessor(
        ftype="regular",
        extension=model_cfg["output_type"],
        # ExoPlaSim 3.4.2 documents integer codes but silently drops derived
        # variables when integers are used. String codes take the correct path.
        variables=[str(code) for code in regular_codes],
        mode="grid",
        times=int(model_cfg["regular_output_bins_per_orbit"]),
        timeaverage=True,
        interpolatetimes=False,
    )
    model.cfgpostprocessor(
        ftype="snapshot",
        extension=model_cfg["output_type"],
        variables=[str(code) for code in SNAPSHOT_CODES],
        mode="grid",
        times=None,
        timeaverage=False,
        interpolatetimes=False,
    )
    staged = stage_surface_extras(run_dir, config)
    spectrum = stage_stellar_spectrum(model, run_dir, stellar_spectrum_path(config))
    verify_stellar_spectrum(model, config)
    # Ozone column scaling, set through the namelist because the Python API does
    # not expose o3scale. The model prescribes an Earth column and derives
    # nothing about it from the host star; Segura et al. (2003) measure 0.794 of
    # Earth's for a K2V host at 1 PAL O2.
    o3 = config["model"].get("ozone_scale")
    if o3 is not None and float(o3) != 1.0:
        model._edit_namelist("radmod_namelist", "O3SCALE", f"{float(o3)}")
        print(f"ozone column scaled to {float(o3)} of Earth's (Segura et al. 2003)")
    # SHORTWAVE_GAS_KEYS is the one list; see its definition for why it is not
    # written out here.
    for key, name, default in SHORTWAVE_GAS_KEYS:
        w = config["model"].get(key)
        if w is not None and float(w) != default:
            model._edit_namelist("radmod_namelist", name, f"{float(w)}")
            print(f"{name} = {float(w)} (radmod.f90 default {default})")

    for name, ppmv in trace_gas_ppmv(config).items():
        model._edit_namelist("radmod_namelist", name, f"{ppmv:.6g}")
        print(f"{name} = {ppmv:.6g} ppmv (radmod.f90 default 0.0, meaning absent)")

    dust = enable_prescribed_dust(model, run_dir, config)
    if dust is not None:
        print(f"prescribed dust ON: {dust['variant']} burden, land-mean band-1 "
              f"optical depth {dust['land_mean_band1_optical_depth']:.4f} x "
              f"{dust['dustsc']}, DUSTQLW {dust['dustqlw']:.5f}. Needs "
              "patches/exoplasim-3.4.2-prescribed-dust.patch and a rebuild.")

    emission = enable_dust_emission(model, run_dir, config)
    if emission is not None:
        print(f"in-model dust emission ON: {emission['source_cells']} source "
              f"cells, aeolian z0 {emission['aeolian_z0_m']:g} m "
              f"({emission['z0_bracket_end']} end). Needs "
              "patches/exoplasim-3.4.2-dust-emission.patch, the two aerosol "
              "patches under it, and a rebuild.")
        if emission["radiatively_active"]:
            print(f"  and RADIATIVELY ACTIVE: l_aerorad = 1, AEROQLW "
                  f"{emission['aeroqlw']:.6f}")
        else:
            print("  and radiatively inert: l_aerorad = 0. dust-13 is the "
                  "decision; model.dust_emission_radiative turns it on.")

    if args.writes_per_day is not None:
        if args.writes_per_day < 1:
            raise ValueError("--writes-per-day must be positive")
        model._edit_namelist("plasim_namelist", "NWPD",
                             str(int(args.writes_per_day)))
        print(f"NWPD = {args.writes_per_day} writes/day "
              f"(default 1 samples 5 diurnal phases; see CLIM-11)")
    declare_cold_start_seed(model, config, args.restart_from is None)
    hyperdiffusion = declare_hyperdiffusion(model, config)
    dry_constants = declare_dry_constants(model, config)
    dynamics_only = declare_dynamics_only(model, config)
    robert_filter = declare_robert_filter(model, config)
    conversion_time_level = declare_conversion_time_level(model, config)
    dealias_conversion = declare_dealias_conversion(model, config)
    declare_storm_diagnostics(model)
    set_low_io(model, args.low_io)
    eco_stream = declare_ecological_stream(
        model, args.ecological_stream, args.eco_interval_steps)
    if eco_stream["enabled"]:
        print(f"ecological stream ON: NECO = 1, NECOSTEP = "
              f"{eco_stream['necostep']} ({eco_stream['interval']}). Written "
              f"to plasim_eco and moved aside per orbit as MOST_ECO.NNNNN; "
              f"read directly, not through pyburn")
    energy_fixer = declare_energy_fixer(model, config)
    if enable_energy_diagnostics(model, config):
        n = register_energy_diagnostic_codes()
        print(f"energy diagnostics on: nenergy=1, {n} codes 360-387 registered "
              "with the postprocessor. Term 15 needs "
              "patches/exoplasim-3.4.2-energy-diagnostics.patch and a rebuild.")
    model.exportcfg(str(run_dir / f"{identifier}.cfg"))
    surface_report = surface_field_report(run_dir, config)

    # CONS-9. AFTER every staging call and BEFORE the run, because the namelist
    # on disk is the only artifact that records what will actually be
    # integrated. Raises, so a run whose radiation does not match its own
    # configuration never starts.
    staged_namelists = verify_staged_namelists(run_dir, config)
    print(f"namelists verified {len(staged_namelists)} config-set keys present "
          f"with the declared values: {', '.join(sorted(staged_namelists))}")

    checks = {
        "N_RUN_STEPS": namelist_value(run_dir / "plasim_namelist", "N_RUN_STEPS"),
        "N_DAYS_PER_YEAR": namelist_value(run_dir / "plasim_namelist", "N_DAYS_PER_YEAR"),
        "NSTPS": namelist_value(run_dir / "plasim_namelist", "NSTPS"),
        "SIDEREAL_YEAR": namelist_value(run_dir / "planet_namelist", "SIDEREAL_YEAR"),
        "ROTSPD": namelist_value(run_dir / "planet_namelist", "ROTSPD"),
        "GSOL0": namelist_value(run_dir / "planet_namelist", "GSOL0"),
        "ECCEN": namelist_value(run_dir / "planet_namelist", "ECCEN"),
        "OBLIQ": namelist_value(run_dir / "planet_namelist", "OBLIQ"),
        "STARBBTEMP": namelist_value(run_dir / "radmod_namelist", "STARBBTEMP"),
        "NSIMPLEALBEDO": namelist_value(run_dir / "radmod_namelist", "NSIMPLEALBEDO"),
    }
    exe_path, exe_provenance = install_run_executable(run_dir, model_cfg, args.binary)

    manifest = {
        "schema_version": 1,
        "run_id": identifier,
        # What this run IS, as opposed to what it is called. The id is a UUID and
        # carries no meaning, so this is the only place the physics is written
        # down in a machine-readable form; exoplasim/scripts/index_runs.py builds
        # exoplasim/runs/INDEX.json from it.
        "physical": physical_fingerprint(config, flux_ratio),
        # Which terrain this run is of, by name, at the top level rather than
        # only inside the config copy below. Everything downstream verifies build
        # identity through lib/provenance.py, and build_climatology was looking
        # for it here and finding nothing, so every climatology went out
        # unstamped and every check on it went quiet.
        "source_build": config.get("source_build"),
        "status": "prepared",
        "config_path": str(config_path),
        "config_sha256": file_sha256(config_path),
        "source_config": config,
        "derived_parameters": derived,
        # The DONOR is the provenance, not the copy this run reads. Recording
        # the copy would name a file that exists only inside this run directory
        # and whose sha nothing else can match, which is the opposite of what a
        # provenance field is for.
        "initial_state": ({"cold_start": True} if restart_seed is None else {
            "cold_start": False,
            "superseded_surface_override": superseded_surface,
            "restart_from": str(seeded_from),
            "restart_from_sha256": file_sha256(seeded_from),
            "restart_from_run": seeded_from.parent.name,
            "accumulators_reset": [n for n, _ in reset_info["zeroed"]],
            # A converted state BEGINS a lineage. The donor is provenance, not
            # continuation: resolution conversion changes the represented state
            # on purpose and the target model owns the equilibrium.
            "conversion": None if conversion is None else {
                "report": str(conversion["__path__"]),
                "converter_version": conversion.get("converter_version"),
                "from": conversion["source"]["truncation"],
                "to": conversion["target_template"]["truncation"],
                "real_bytes": [conversion["source"]["real_bytes"],
                               conversion["target_template"]["real_bytes"]],
                "donor_run": (conversion.get("source_manifest") or {}).get("run_id"),
                "template": conversion["target_template"].get("path"),
                "template_cut_from_run":
                    conversion["target_template"].get("cut_from_run"),
                "records_the_model_must_rebuild":
                    conversion.get("expected_to_change_in_model_fixup"),
                "whole_run_accumulators_restarted":
                    conversion.get("whole_run_accumulators_restarted"),
            },
            "seed_copy_sha256": file_sha256(restart_seed),
            "note": "Initial condition only; equilibrium is set by the forcing. "
                    "The accumulator records listed were zeroed in a copy before "
                    "the run read it, so this run does not inherit the donor's "
                    "partial accumulation window (CLIM-31); every other record "
                    "is byte-identical to the donor.",
        }),
        "geography": {
            "landmap": str(landmap),
            "landmap_sha256": file_sha256(landmap),
            "topomap": str(topomap),
            "topomap_sha256": file_sha256(topomap),
        },
        # ExoPlaSim names the binary most_plasim_t<res>_l<layers>_p<ncpus>.x, with
        # no precision in the path, and reuses whatever is already there rather
        # than recompiling. So the compiled precision is invisible shared state.
        # Recording the binary's digest is what makes a silent swap auditable.
        "executable": {
            "path": str(exe_path),
            "sha256": file_sha256(exe_path),
            # WHICH MANIFEST ENTRY, or which ARM. The sha above says what file
            # ran; on its own it cannot say whether that file was one the
            # registry knows, and a run whose whole life was integrated by an
            # unregistered binary is unattributable after the fact. The prepare
            # refuses that case outright, so this block records the verdict it
            # passed. world-bdb5, world-u5pf.
            **exe_provenance,
            "note": ("Built precision is not encoded in the filename. The "
                     "registry's naming carries no precision at all, which is "
                     "why `arm` below is the only statement of an arm's."),
        },
        # An ARM RUN IS NOT PRODUCTION. Set false when the caller named the
        # binary rather than letting the registry's naming resolve it: the
        # executable is then one build_model.py refused to publish, and no
        # climatology built on it belongs to the canonical lineage.
        # `continue_exoplasim.py` refuses a post_equilibrium_climatology
        # segment on a run carrying false here.
        "canonical_lineage_eligible": exe_provenance.get("arm") is None,
        "software": {
            "python": platform.python_version(),
            "exoplasim": getattr(exo, "__version__", "3.4.2"),
            "numpy": np.__version__,
            "gfortran": command_version(["gfortran", "--version"]),
        },
        "namelist_checks": {**checks, **staged_namelists},
        "surface_fields": surface_report,
        "surface_fields_staged": staged,
        # Per-code hashes, so a surface field that CHANGED CONTENT under the same
        # code is detectable. The set of codes alone is not enough: the albedo
        # rebuilt with lakes composited is still code 174, and a restart would
        # have discarded it silently.
        "surface_field_sha256": {
            str(c): file_sha256(surface_sra(config, c))
            for c in sorted(intended_surface_codes(config))
            if surface_sra(config, c).is_file()},
        "stellar_spectrum": spectrum,
        # The spectrum by CONTENT, because the line above is a filename and the
        # file is regenerated in place. `continue_exoplasim.py` compares this on
        # every resume; without it a regenerated `k25v.dat` changes the star
        # mid-run and nothing says so. CONS-3.
        "stellar_spectrum_digest": stellar_spectrum_digest(config),
        # Null when the run has no prescribed dust, which is most of them. The
        # values are copied from the field's own provenance rather than the
        # config, so a run says what burden and what longwave ratio it actually
        # carried instead of what a setting asked for.
        "prescribed_dust": dust,
        # What damping this run actually got, including the rounded nhdiff, so a
        # result can be read against the operator that produced it rather than
        # against the config that was meant to.
        "hyperdiffusion": hyperdiffusion,
        # The dry thermodynamic constants, the cold-start profile and the top
        # sponge, all of which used to be ExoPlaSim's Earth defaults and none of
        # which appeared in a manifest.
        "dry_constants": dry_constants,
        "dynamics_only": dynamics_only,
        "conversion_time_level": conversion_time_level,
        "dealias_conversion": dealias_conversion,
        "robert_filter": robert_filter,
        # A CORRECTION and not physics; what it is correcting is world-0ov and
        # the fixer itself is world-mzy. On the manifest because whether a run
        # carries it changes what its surface fluxes mean.
        "energy_fixer": energy_fixer,
        # Null when the run has no interactive emission, which is most of them.
        # Like prescribed_dust, the values are copied from the fields' own
        # provenance rather than from the config, so the manifest says what the
        # run carried and not what a setting asked for.
        "dust_emission": emission,
        # What stack the threaded model was actually launched with. On the
        # manifest because a run that dies for want of it dies inside the model,
        # which reads as the model crashing rather than as a launch setting.
        "thread_stack": thread_stack,
        "postprocessor": {
            "regular_codes": regular_codes,
            "energy_diagnostics": energy_diagnostics_enabled(config),
            "snapshot_codes": SNAPSHOT_CODES,
        },
    }
    manifest_path = run_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {"run_dir": str(run_dir), "derived": derived, "namelists": checks},
            indent=2,
        )
    )

    if args.run_years:
        try:
            model.run(years=args.run_years, crashifbroken=True, clean=True)
            manifest["output_validation"] = validate_outputs(run_dir)
            manifest["status"] = (
                "smoke_complete" if args.run_years == 1 else "run_complete"
            )
            manifest["completed_orbits"] = args.run_years
            if energy_fixer:
                applied = read_applied_energy_fix(run_dir)
                if applied is None:
                    raise RuntimeError(
                        "the energy fixer was declared on and the run's diag "
                        "carries no readable correction. The fixer hides the "
                        "defect it compensates, so a run that cannot report "
                        "what it applied is a run whose fluxes nobody can "
                        "check. world-mzy.")
                manifest["energy_fixer_applied"] = applied
                print(f"energy fixer applied {applied['applied_w_m2_mean']:+.4f} "
                      f"W/m2 on average over the second half "
                      f"({applied['applied_w_m2_min']:+.3f} to "
                      f"{applied['applied_w_m2_max']:+.3f}); that is the size of "
                      f"the defect on world-0ov, not a physical flux")
            # REGISTER THE BLOCK AS A SEGMENT, like a continuation does.
            # Without this the orbits this script writes carry no purpose and
            # no `low_io`, and CLIM-9's two defaults then disagree about them:
            # unlabelled counts as production for convergence and as TAINTED
            # for I/O regime. That was harmless while the block was always
            # NLOWIO = 0 and merely conservative; with the cheap regime as the
            # default it would be a run whose first block really is
            # accumulations and says so nowhere.
            #
            # The purpose is `spinup` and is not a parameter: this script
            # prepares a run and integrates it from a cold start or a seed,
            # which is the definition of one. Orbits meant to be read as data
            # come from `continue_exoplasim.py --purpose
            # post_equilibrium_climatology`, on a run something has assessed.
            manifest.setdefault("segments", []).append({
                "start_year_index": 0,
                "end_year_index": args.run_years - 1,
                "seasonal_output": True,
                "low_io": bool(args.low_io),
                "ecological_stream": eco_stream,
                "high_cadence": False,
                "purpose": "spinup",
                # WHY A FIRST RECORD IS REFUSABLE, and it is no longer CLIM-31.
                #
                # CLIM-31 was the donor's partial accumulation window arriving
                # in the seed: under low I/O the land mask, which cannot vary,
                # read 0.95341 in the first record. That is FIXED. The
                # accumulator reset takes its record set and its clean values
                # from `restart_schema.py` rather than from a hand list that
                # covered 49 of 114, and a converted seeded run now reads
                # exactly 1.0 there.
                #
                # What remains is narrower and applies only to a CONVERTED
                # state. Its derived surface records -- albedo, roughness --
                # come from the target template and are not consistent with the
                # prognostics that were remapped onto it, so the model's first
                # timestep runs on the template's albedo before rebuilding it.
                # One timestep, and then it is gone: measured at T21 to T42
                # against a template 0.104 out in area-weighted albedo, the
                # first record carries 0.213 W/m2 of reflected shortwave, which
                # is 33.6 W/m2 spread over the 160 steps the record averages.
                # Under clean I/O the first record is an instantaneous sample
                # rather than an accumulation, so it does not pick this up.
                # world-eyb.
                "first_record_tainted": bool(args.low_io and conversion is not None),
                "stellar_spectrum_digest": stellar_spectrum_digest(config),
                "finished_utc": datetime.now(timezone.utc).isoformat(),
            })
        except Exception:
            manifest["status"] = "failed"
            manifest_path.write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
            )
            raise
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
