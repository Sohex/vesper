#!/usr/bin/env python3
"""Build cGENIE through a config, and measure what one EMIC ocean-year costs.

    python analysis/cgenie_cost.py --years 20 100 --repeats 2 --probe72 --verify
    python analysis/cgenie_cost.py --case worbe2_36x36x8 --years 20

Worldbuilding. Vesper is an invented planet; this measures a candidate ocean
component on this host, not any real ocean.

## Why the build has to run through a config

`vendor/cgenie/genie-main/genie_control.f90` fixes the atmosphere as
`ilon1_atm = GENIENX, ilat1_atm = GENIENY` and the ocean and sea ice as
`GOLDSTEINNLONS/NLATS/NLEVS`, all cpp macros with `#ifndef` fallbacks. The
fallbacks do not agree with each other: the atmosphere falls back to the
64 x 32 IGCM grid and the ocean to 36 x 36, so a bare `make` fails in `genie.F`
where a sea-ice array accumulates an atmosphere field. That is a missing
configuration, not a source defect, and the macros cannot be supplied on the
make line either, because `makefile.arc` reads them from `GENIE_FPPFLAGS`,
which `genie.job` composes from the `<build>` block of a config. So one
executable exists per grid, exactly as ExoPlaSim compiles one per
(resolution, layers, ranks), and CLAUDE.md rule 4 applies here too.

## What this measures

Per case: a full `make cleanall` and rebuild through a generated XML config,
then the model run at two lengths. Two lengths are the point, because cGENIE's
initialisation reads the topography, the wind stress fields and the restoring
climatologies before it integrates anything, and a single run cannot tell that
apart from the integration. Fitting wall clock against simulated years separates
the per-year slope, which is what a spin-up costs, from the intercept, which is
what every pass pays once. Whether the intercept is large is the measurement,
not the assumption: `--stability` reports the same split against `nyear` at one
grid.

The component is serial, and more completely than the build flags suggest.
`makefile.arc`'s gfortran section carries `-fopenmp` commented out, but turning
it on would parallelise nothing: every `$omp` string in `genie-goldstein`,
`genie-embm`, `genie-biogem` and `genie-main` sits behind a `c` or a `!` in
column one, so none of them is a sentinel and none is a directive. What is left
live is three `use omp_lib` statements with no call against them. MPI appears
only in `genie-plasim/src/fortran/plasimmod.f90`, the vendored PlaSim this
project does not use. So there is no thread-count axis to sweep here: cost is
one core against simulated years, cores cannot absorb a grid refinement, and
parallelism is across EXPERIMENTS rather than within a run. That is why an
offline ocean can run beside a commissioning rather than queue behind it.

## What the horizontal arm is a test of, stated before it ran

`ubarsolv.f` eliminates over `do i=1,n*m-1` with an inner width of at most
`n+1`, where `n = imax` and `m = jmax+1`. Work per call is therefore about
`imax^2 * jmax`: QUADRATIC in the number of longitudes and only LINEAR in the
number of latitudes. Tracer work in `tstepo` is `imax*jmax*kmax*lmax`,
symmetric in the two. Two cases here have nearly the same cell count and very
different barotropic work, which makes the two hypotheses separable:

    36 x 60 :  imax^2*jmax =  77760
    64 x 32 :  imax^2*jmax = 131072   ratio 1.69

If the barotropic solve dominates, the 64 x 32 case costs about 1.69 times the
36 x 60 one per model year; if tracer work dominates it costs about 0.95 times.
The decision rule, fixed before the runs: above 1.30 reads barotropic-dominated,
below 1.10 reads tracer-dominated, and between the two leaves the r^4 estimate
unverified. The ratio does not depend on the level count, so the arm is run at
the eight levels all three of those topographies are built for.

## The timestep is a namelist consequence, and nothing checks it

GOLDSTEIN's tracer timestep is `dt = sodaylen*yearlen/(nyear*tsc)`, set by the
namelist `nyear` alone. Neither `genie-goldstein` nor `genie-embm` contains a
Courant, CFL or stability test of any kind, so the criterion has to be applied
from outside. This script states it explicitly per case and reports the margin;
see `cfl_margin` in the output. The criterion is fixed here before any run:
the advective Courant number on the smallest zonal cell width must stay below
one, evaluated at a nominal surface current speed, and the vertical and
horizontal diffusive numbers must also stay below one half.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CGENIE = PROJECT_ROOT / "vendor" / "cgenie"
GENIE_MAIN = CGENIE / "genie-main"
CONFIG_DIR = GENIE_MAIN / "configs"
OUT_ROOT = Path(os.environ.get("CGENIE_OUT", Path.home() / "cgenie_output"))
NETCDF_DIR = os.environ.get("NETCDF_DIR", "/usr")
PROBE_DIR = Path(os.environ.get("CGENIE_PROBE", "/tmp/cgenie_probe72"))

# Non-dimensionalising scales, hardcoded in
# genie-goldstein/src/fortran/initialise_goldstein.F and exposed nowhere.
RSC_M = 6.37e6
# Nominal surface current for the advective criterion, m/s. GOLDSTEIN is
# frictional-geostrophic, so there is no gravity-wave constraint; advection of
# tracers by the diagnosed flow is what can go unstable.
U_NOMINAL_MS = 0.5
DIFF_H_M2S = 2000.0
DIFF_V_M2S = 1.0e-4
DSC_M = 5000.0

# Every case is a topography that ships complete at the stated grid: .k1,
# .paths, .psiles, the four wind-stress components, the two wind-speed fields
# and the restoring climatologies. Grids are (nlons, nlats, nlevs).
#
# LEVEL COUNT IS A PROPERTY OF THE TOPOGRAPHY, NOT A FREE PARAMETER. A .k1
# entry is the index of the level a column's floor sits on, so the highest wet
# value in the file IS the grid's level count and land is coded above it. Give
# a 36 x 36 x 8 topography sixteen levels and every column becomes bottomless;
# give it thirty-two and initialisation stops with "wind stress not defined
# outside domain" -- a message about the wrong thing, returning exit status
# zero, which is why run_ok below reads the log rather than the status.
# The vertical arm of this sweep therefore changes topography with the levels.
CASES = {
    "worbe2_36x36x8": dict(world="worbe2", nlons=36, nlats=36, nlevs=8, prefix=""),
    "worjh2_36x36x16": dict(world="worjh2", nlons=36, nlats=36, nlevs=16, prefix=""),
    "worri4_36x36x32": dict(world="worri4", nlons=36, nlats=36, nlevs=32, prefix=""),
    "g3660l_36x60x8": dict(world="g3660l", nlons=36, nlats=60, nlevs=8, prefix="g3660l_"),
    "igcmv3_64x32x8": dict(world="igcmv3", nlons=64, nlats=32, nlevs=8, prefix="igcmv3_"),
    # The control for the code model. Same topography, same grid and same run
    # lengths as worjh2_36x36x16, compiled and linked -mcmodel=medium instead of
    # the shipped small. The 72 x 72 probe below has to be built that way, so
    # without this pair the flag's own cost would sit inside the grid's.
    "worjh2_36x36x16_medium": dict(world="worjh2", nlons=36, nlats=36, nlevs=16,
                                   prefix="", mcmodel="medium"),
}

# A doubling of the shipped horizontal grid. It is separate from CASES because
# its forcing is not complete in the tree: see build_probe72_inputs. Selected
# with --probe72, never by default.
PROBE72 = {
    "dan_72_72x72x16_probe": dict(
        world="dan_72", nlons=72, nlats=72, nlevs=16, prefix="", synthetic_forcing=True,
        mcmodel="medium",
        # At the shipped nyear = 100 this grid needs ndta >= 8: below that
        # EMBM's surface solve stops the model, and it is not close. 10 is the
        # smallest round value inside the boundary. It does not distort the
        # comparison with the 36 x 36 cases at ndta = 5, because the --stability
        # sweep shows the wall clock here is flat in ndta from 8 to 20: at this
        # grid EMBM sets the timestep and the ocean sets the cost.
        ndta=10,
        files=dict(
            xu="dan_72_NCEP-DOE_Reanalysis_2_average_taux_u.dat",
            yu="dan_72_NCEP-DOE_Reanalysis_2_average_tauy_u.dat",
            xv="dan_72_NCEP-DOE_Reanalysis_2_average_taux_v.dat",
            yv="dan_72_NCEP-DOE_Reanalysis_2_average_tauy_v.dat",
            uw="probe72_uncep.silo", vw="probe72_vncep.silo",
            # Named but never opened: every reader of these is behind
            # ctrl_diagend or behind debug_loop, and both are off here.
            ta="probe72_uncep.silo", qa="probe72_uncep.silo",
            gold_t="probe72_uncep.silo", gold_s="probe72_uncep.silo",
            indir_embm=PROBE_DIR / "embm",
            indir_gold=PROBE_DIR / "goldstein",
            indir_sic=PROBE_DIR / "goldsteinseaice",
            diagend=".false.",
        ),
    ),
}

CONFIG_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<job author="analysis/cgenie_cost.py">
\t<vars>
\t\t<var name="EXPID">{expid}</var>
\t</vars>
\t<config>
\t\t<model name="goldstein"/>
\t\t<model name="goldsteinseaice"/>
\t\t<model name="embm"/>
\t</config>
\t<parameters>
\t\t<control>
\t\t\t<param name="write_flag_atm">.false.</param>
\t\t\t<param name="write_flag_ocn">.false.</param>
\t\t\t<param name="write_flag_sic">.false.</param>
\t\t\t<param name="ksic_loop">{kloop}</param>
\t\t\t<param name="kocn_loop">{kloop}</param>
\t\t\t<param name="koverall_total">{koverall}</param>
\t\t\t<param name="dt_write">{big}</param>
\t\t\t<param name="genie_timestep">{genie_timestep}</param>
\t\t\t<param name="lgraphics">.false.</param>
\t\t</control>
\t\t<model name="goldstein">
\t\t\t<param name="world">{world}</param>
\t\t\t<param name="nyear">{nyear}</param>
\t\t\t<param name="debug_loop">{debug_loop}</param>
\t\t\t<param name="npstp">{npstp}</param>
\t\t\t<param name="ianav">{big}</param>
\t\t\t<param name="itstp">{big}</param>
\t\t\t<param name="iwstp">{big}</param>
\t\t\t<param name="temp0">10.</param>
\t\t\t<param name="temp1">10.</param>
\t\t\t<param name="indir_name">{indir_gold}</param>
\t\t\t<param name="ctrl_diagend">{diagend}</param>
\t\t\t<param name="tdatafile">{gold_t}</param>
\t\t\t<param name="sdatafile">{gold_s}</param>
\t\t\t<paramArray name="diff">
\t\t\t\t<param index="1">{diff_h}</param>
\t\t\t\t<param index="2">{diff_v}</param>
\t\t\t</paramArray>
\t\t</model>
\t\t<model name="goldsteinseaice">
\t\t\t<param name="world">{world}</param>
\t\t\t<param name="indir_name">{indir_sic}</param>
\t\t\t<param name="nyear">{nyear}</param>
\t\t\t<param name="npstp">{big}</param>
\t\t\t<param name="ianav">{big}</param>
\t\t\t<param name="itstp">{big}</param>
\t\t\t<param name="iwstp">{big}</param>
\t\t</model>
\t\t<model name="embm">
\t\t\t<param name="world">{world}</param>
\t\t\t<param name="nyear">{nyear}</param>
\t\t\t<param name="ndta">{ndta}</param>
\t\t\t<param name="npstp">{big}</param>
\t\t\t<param name="ianav">{big}</param>
\t\t\t<param name="itstp">{big}</param>
\t\t\t<param name="iwstp">{big}</param>
\t\t\t<param name="tatm">0.0</param>
\t\t\t<param name="indir_name">{indir_embm}</param>
\t\t\t<param name="xu_wstress">{xu}</param>
\t\t\t<param name="yu_wstress">{yu}</param>
\t\t\t<param name="xv_wstress">{xv}</param>
\t\t\t<param name="yv_wstress">{yv}</param>
\t\t\t<param name="u_wspeed">{uw}</param>
\t\t\t<param name="v_wspeed">{vw}</param>
\t\t\t<param name="tdatafile">{ta}</param>
\t\t\t<param name="qdatafile">{qa}</param>
\t\t</model>
\t</parameters>
\t<build>
\t\t<make-arg name="IGCMATMOSDP">TRUE</make-arg>
\t\t<make-arg name="GENIEDP">TRUE</make-arg>
\t\t<macro handle="GENIENXOPTS" status="defined">
\t\t\t<identifier>GENIENX</identifier>
\t\t\t<replacement>{nlons}</replacement>
\t\t</macro>
\t\t<macro handle="GENIENYOPTS" status="defined">
\t\t\t<identifier>GENIENY</identifier>
\t\t\t<replacement>{nlats}</replacement>
\t\t</macro>
\t\t<macro handle="GOLDSTEINNLONSOPTS" status="defined">
\t\t\t<identifier>GOLDSTEINNLONS</identifier>
\t\t\t<replacement>{nlons}</replacement>
\t\t</macro>
\t\t<macro handle="GOLDSTEINNLATSOPTS" status="defined">
\t\t\t<identifier>GOLDSTEINNLATS</identifier>
\t\t\t<replacement>{nlats}</replacement>
\t\t</macro>
\t\t<macro handle="GOLDSTEINNLEVSOPTS" status="defined">
\t\t\t<identifier>GOLDSTEINNLEVS</identifier>
\t\t\t<replacement>{nlevs}</replacement>
\t\t</macro>
\t\t<macro handle="GOLDSTEINMAXISLESOPTS" status="defined">
\t\t\t<identifier>GOLDSTEINMAXISLES</identifier>
\t\t\t<replacement>{maxisles}</replacement>
\t\t</macro>
\t</build>
</job>
"""


# Set by --mcmodel. Every field cGENIE holds is in a named COMMON block sized
# from the grid macros, and user.mak compiles with -fno-automatic, so all of it
# is static. genie-embm/src/fortran/embm.cmn dimensions its seasonal arrays
# (maxi, maxj, maxnyr) with maxnyr = 400, which alone is about 286 kB of static
# data per surface cell. Past roughly a gigabyte of COMMON the default
# -mcmodel=small cannot reach it and the LINK fails with "relocation truncated
# to fit: R_X86_64_PC32 ... defined in COMMON section". That is the ceiling this
# option exists to move, and moving it is a build flag rather than a source
# change.
MCMODEL = ""


def make_args(target: str = "") -> str:
    """Overrides this host needs. user.mak expects the tree at ~/cgenie.muffin
    and netCDF under /usr/local; neither is true here, and the vendored tree is
    read through the repository rather than a home-directory checkout."""
    parts = [
        target,
        f"GENIE_ROOT={CGENIE}",
        f"RUNTIME_ROOT={CGENIE}",
        f"NETCDF_DIR={NETCDF_DIR}",
        f"OUT_DIR={OUT_ROOT}",
        "-j 1",
    ]
    if MCMODEL:
        parts += [f"GENIE_FFLAGS=-mcmodel={MCMODEL}",
                  f"GENIE_LDFLAGS=-mcmodel={MCMODEL}"]
    return " ".join(p for p in parts if p)


def island_count(world: str) -> int:
    """Highest island index in the topography's .psiles file. This is a
    compile-time bound, GOLDSTEINMAXISLES, so a geography with more islands
    than the built executable allows needs a rebuild before it runs."""
    path = CGENIE / "genie-goldstein" / "data" / "input" / f"{world}.psiles"
    highest = 0
    for line in path.read_text().split("\n"):
        for tok in line.split():
            try:
                highest = max(highest, int(tok))
            except ValueError:
                pass
    return highest


def cfl(nlons: int, nlats: int, nlevs: int, nyear: int) -> dict:
    """The stability criterion, stated before any run.

    The grid is equal-area in sin(latitude), so the zonal cell width shrinks
    towards the poles as cos(lat) and the narrowest cell sets the bound."""
    dt = 86400.0 * 365.25 / nyear
    dx_eq = 2.0 * math.pi * RSC_M / nlons
    # centre of the polemost row on an equal-area grid
    sin_top = 1.0 - 1.0 / nlats
    dx_min = dx_eq * math.sqrt(max(1.0 - sin_top * sin_top, 1e-12))
    dy = math.pi * RSC_M / nlats
    dz = DSC_M / nlevs
    return {
        "dt_s": dt,
        "dx_min_m": dx_min,
        "dz_m": dz,
        "courant_advective": U_NOMINAL_MS * dt / dx_min,
        "diffusive_horizontal": DIFF_H_M2S * dt / (dx_min * dx_min),
        "diffusive_vertical": DIFF_V_M2S * dt / (dz * dz),
        "diffusive_meridional": DIFF_H_M2S * dt / (dy * dy),
    }


def write_config(path: Path, case: dict, years: int, nyear: int, maxisles: int,
                 npstp: int | None = None, debug_loop: bool = False,
                 ndta: int = 5) -> None:
    # THE ATMOSPHERE'S TIMESTEP AND ITS CALL FREQUENCY ARE TWO SETTINGS AND
    # THEY MUST AGREE. `initialise_embm.F:578` sets EMBM's step as
    # `dtatm = dt_ocean/ndta`, but how often EMBM is CALLED is `katm_loop` in
    # `genie.F:300`, which defaults to 1, so EMBM steps once per genie step.
    # With `kocn_loop` genie steps per ocean step, EMBM therefore advances
    # `kocn_loop/ndta` times the ocean's time in the same interval. That is one
    # only when `kocn_loop == ndta`, which is why the shipped default pairs
    # ndta = 5 with kocn_loop = 5. Raising ndta on its own does not sub-step the
    # atmosphere: it slows the atmosphere's clock relative to the ocean's, and
    # the run is then not a model of anything even where it is stable. So the
    # ocean loop count is tied to ndta here rather than fixed at 5.
    kloop = ndta
    koverall = years * kloop * nyear
    pre = case["prefix"]
    files = dict(
        xu=f"{pre}taux_u.interp", yu=f"{pre}tauy_u.interp",
        xv=f"{pre}taux_v.interp", yv=f"{pre}tauy_v.interp",
        uw=f"{pre}uncep.silo", vw=f"{pre}vncep.silo",
        ta=f"{pre}ta_ncep.silo", qa=f"{pre}qa_ncep.silo",
        gold_t=f"{pre}tempann.silo", gold_s=f"{pre}saliann.silo",
        indir_embm=CGENIE / "genie-embm" / "data" / "input",
        indir_gold=CGENIE / "genie-goldstein" / "data" / "input",
        indir_sic=CGENIE / "genie-goldsteinseaice" / "data" / "input",
        diagend=".true.",
    )
    files.update(case.get("files", {}))
    path.write_text(
        CONFIG_TEMPLATE.format(
            expid=path.stem,
            world=case["world"],
            **files,
            nlons=case["nlons"],
            nlats=case["nlats"],
            nlevs=case["nlevs"],
            nyear=nyear,
            kloop=kloop,
            koverall=koverall,
            big=koverall + 1,
            genie_timestep=f"{86400.0 * 365.25 / (kloop * nyear):.4f}",
            ndta=ndta,
            npstp=npstp if npstp is not None else koverall + 1,
            debug_loop=".true." if debug_loop else ".false.",
            diff_h=DIFF_H_M2S,
            diff_v=DIFF_V_M2S,
            maxisles=maxisles,
        )
    )


# This host runs other work. cGENIE is one serial process on a 32-thread part,
# so a run is not competing for a core, but it does compete for last-level cache
# and memory bandwidth, and that is a real effect on a model whose state is one
# large block of static COMMON. Two things follow, and both are in the numbers
# below rather than in a caveat: every timed quantity is the MINIMUM over
# repeats, which is the least-contended sample rather than an average of a
# contaminated population, and the load average is recorded beside each run so a
# reader can see what the machine was doing.
def run(cmd: list[str], log: Path, cwd: Path = GENIE_MAIN) -> tuple[int, float]:
    start = time.perf_counter()
    with log.open("w") as fh:
        rc = subprocess.call(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)
    return rc, time.perf_counter() - start


def loadavg() -> float:
    return round(os.getloadavg()[0], 2)


def run_ok(log: Path) -> tuple[bool, str]:
    """The executable's exit status does not distinguish a completed run from a
    Fortran STOP, so completion is read from the model's own closing banner."""
    text = log.read_text(errors="replace")
    if "Shutdown complete; home time" not in text:
        for line in text.split("\n"):
            if line.startswith("STOP") or "ERROR" in line:
                return False, line.strip()
        return False, "no shutdown banner and no message"
    return True, ""


def model_seconds(log: Path) -> float | None:
    """genie.job wraps the executable in `time`, so the log carries the model's
    own wall clock separately from the xsltproc configuration around it."""
    for line in log.read_text(errors="replace").split("\n"):
        if line.startswith("real"):
            token = line.split()[-1]
            minutes, _, seconds = token.partition("m")
            return int(minutes) * 60 + float(seconds.rstrip("s"))
    return None


def job(config_rel: str, extra: list[str], target: str = "") -> list[str]:
    return [
        "./genie.job",
        *extra,
        "-f", config_rel,
        "-o", str(OUT_ROOT),
        "-c", str(CGENIE),
        "-g", str(CGENIE),
        "-h", ".",
        "-m", make_args(target),
    ]


def provenance() -> dict:
    def sh(*cmd: str) -> str:
        try:
            return subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL).strip()
        except Exception:
            return "unavailable"

    cpu = "unknown"
    try:
        for line in Path("/proc/cpuinfo").read_text().split("\n"):
            if line.startswith("model name"):
                cpu = line.split(":", 1)[1].strip()
                break
    except OSError:
        pass
    return {
        "repo_head": sh("git", "-C", str(PROJECT_ROOT), "rev-parse", "HEAD"),
        "repo_dirty": bool(sh("git", "-C", str(PROJECT_ROOT), "status", "--porcelain", "vendor/cgenie")),
        "fortran": sh("gfortran", "--version").split("\n")[0],
        "netcdf_c": sh("nc-config", "--version"),
        "netcdf_fortran": sh("nf-config", "--version"),
        "netcdf_flibs": sh("nf-config", "--flibs"),
        "host": platform.node(),
        "cores": os.cpu_count(),
        "load_average_at_start": loadavg(),
        "timing_estimator": (
            "minimum over --repeats, of the model's own `time real` rather than the"
            " wrapper's, on a host that runs other work. cGENIE is one serial process"
            " and cores were free throughout; cache and memory bandwidth were not"
            " exclusively its own. Each run records the load average it finished under."
        ),
        "cpu": cpu,
        "kernel": platform.release(),
        "build_type": "SHIP (user.mak default), gfortran, serial",
        "parallelism": (
            "None available. -fopenmp is commented out in makefile.arc, and every $omp"
            " string in goldstein, embm, biogem and genie-main is commented out in column"
            " one, so it is not a sentinel: enabling the flag would parallelise nothing."
            " Three live `use omp_lib` statements have no call against them. MPI exists only"
            " in genie-plasim, which this project does not use. Cost is one core against"
            " simulated years; parallelism is across experiments, not within a run."
        ),
        "code_model": MCMODEL or "small (compiler default, what user.mak leaves in place)",
        "make_overrides": make_args(),
        "build_command": (
            "cd vendor/cgenie/genie-main && ./genie.job -x -f <config> "
            "-o <outdir> -c <cgenie root> -g <cgenie root> -h . "
            "-m 'genie.exe GENIE_ROOT=<cgenie root> RUNTIME_ROOT=<cgenie root> "
            "NETCDF_DIR=/usr OUT_DIR=<outdir> -j 1'"
        ),
        "notes": [
            "The make target genie.exe is passed through -m because the default target"
            " also builds nccompare, whose src/c/compare.cpp needs the legacy netCDF C++"
            " header netcdf.hh, absent on this host.",
            "-j 1 is required: the dependency generator genie-main/finc.py is Python 2 and"
            " fails under this host's python3, so no .d files are produced and a parallel"
            " make races on .mod files.",
        ],
    }


def build_probe72_inputs() -> Path:
    """Assemble a 72 x 72 input set, and say plainly what is real in it.

    `genie-goldstein/data/input/dan_72.k1`, `.paths` and `.psiles` are a real
    72 x 72 x 16 topography, and `genie-embm/data/input` carries dan_72's four
    real wind-stress components. What does NOT ship for it is the pair of
    advective wind-speed fields EMBM reads unconditionally -- `flag_wind` is
    hardcoded false at `initialise_embm.F:303` and reachable from no namelist,
    so the externally-supplied-wind branch is dead and the files must exist.
    Those two are made here by replicating the 36 x 36 fields into 2 x 2 blocks.

    That makes this a COST probe and not a simulation. Doubling a field's
    resolution by replication changes no arithmetic the timestep performs, so
    the wall clock is the wall clock of a 72 x 72 x 16 ocean; the circulation it
    produces means nothing and is not read."""
    src_embm = CGENIE / "genie-embm" / "data" / "input"
    src_gold = CGENIE / "genie-goldstein" / "data" / "input"
    src_sic = CGENIE / "genie-goldsteinseaice" / "data" / "input"
    shutil.rmtree(PROBE_DIR, ignore_errors=True)
    for name, src in (("embm", src_embm), ("goldstein", src_gold), ("goldsteinseaice", src_sic)):
        dst = PROBE_DIR / name
        dst.mkdir(parents=True)
        for f in src.iterdir():
            (dst / f.name).symlink_to(f)

    def upsample(src: Path, dst: Path, n: int = 36) -> None:
        vals = [float(t) for t in src.read_text().split()]
        assert len(vals) == n * n, (src, len(vals))
        rows = [vals[j * n:(j + 1) * n] for j in range(n)]
        with dst.open("w") as fh:
            for row in rows:
                doubled = [v for v in row for _ in (0, 1)]
                for _ in (0, 1):
                    fh.writelines(f"{v:14.6E}\n" for v in doubled)

    for base in ("uncep.silo", "vncep.silo"):
        upsample(src_embm / base, PROBE_DIR / "embm" / f"probe72_{base}")
    # EMBM and the sea ice each read <world>.k1 from their own input directory,
    # and dan_72's only copy is GOLDSTEIN's.
    for name in ("embm", "goldsteinseaice"):
        link = PROBE_DIR / name / "dan_72.k1"
        link.unlink(missing_ok=True)
        link.symlink_to(src_gold / "dan_72.k1")
    return PROBE_DIR


def biogem_cost(logdir: Path, years: list[int], repeats: int, reuse: bool = False,
                perf: bool = False) -> dict:
    """What the biogeochemistry costs, from the shipped BIOGEM regression case.

    The physics-only sweep above is EMBM plus GOLDSTEIN plus sea ice with two
    tracers, and the published EMIC costs this project has been quoting are for
    a configuration carrying ocean biogeochemistry. They are not the same
    measurement, so this runs `configs/eb_go_gs_ac_bg_test.xml` unaltered except
    for its length: fourteen GOLDSTEIN tracers, ATCHEM and BIOGEM on, the same
    36 x 36 x 8 grid as the cheapest physics case, so the ratio between them is
    the price of the biogeochemistry and nothing else."""
    src = (CONFIG_DIR / "eb_go_gs_ac_bg_test.xml").read_text()
    cfg = CONFIG_DIR / "bench_biogem.xml"
    out = {"from": "configs/eb_go_gs_ac_bg_test.xml", "nlons": 36, "nlats": 36,
           "nlevs": 8, "tracers": 14, "runs": []}

    def write(n: int) -> None:
        koverall = n * 5 * 100
        t = re.sub(r'(<param name="koverall_total">)\d+', rf"\g<1>{koverall}", src)
        t = re.sub(r'(<param name="dt_write">)\d+', rf"\g<1>{koverall}", t)
        t = re.sub(r'(<param name="par_misc_t_runtime">)\d+', rf"\g<1>{n}", t)
        t = t.replace("genie_eb_go_gs_ac_bg", cfg.stem, 1)
        cfg.write_text(t)

    per_rep: dict[int, dict[int, float]] = {}
    counts: dict[int, float] = {}
    if reuse:
        for n in years:
            best = None
            for rep in range(repeats):
                log = logdir / f"biogem.run{n}.{rep}.log"
                if not log.exists() or not run_ok(log)[0]:
                    continue
                inner = model_seconds(log)
                if inner is None:
                    continue
                per_rep.setdefault(rep, {})[n] = inner
                best = inner if best is None else min(best, inner)
            if best is not None:
                out["runs"].append({"years": n, "model_seconds": best})
        out["build_ok"] = bool(per_rep)
        out.update(slope_and_intercept(per_rep, years))
        return out
    try:
        write(years[0])
        run(["/usr/bin/make", *make_args().split(), "cleanall"], logdir / "biogem.clean.log")
        rc, secs = run(job(f"configs/{cfg.name}", ["-x"], target="genie.exe"),
                       logdir / "biogem.build.log")
        out["build_seconds"] = round(secs, 2)
        out["build_ok"] = rc == 0
        if rc != 0:
            return out
        for n in years:
            write(n)
            best = None
            for rep in range(repeats):
                log = logdir / f"biogem.run{n}.{rep}.log"
                rc, secs = run(job(f"configs/{cfg.name}", ["-z"]), log)
                ok, why = run_ok(log)
                if rc != 0 or not ok:
                    out["error"] = f"run {n}y failed ({why or rc}), see {log}"
                    return out
                inner = model_seconds(log)
                if inner is not None:
                    per_rep.setdefault(rep, {})[n] = inner
                best = inner if best is None else min(best, inner)
            outdir = OUT_ROOT / cfg.stem
            size = sum(f.stat().st_size for f in outdir.rglob("*")
                       if f.is_file() and f.name != "genie.exe")
            row = {"years": n, "model_seconds": best, "output_bytes": size}
            if perf:
                row["perf"] = perf_run(outdir, logdir / f"biogem.perf{n}.log")
                counts[n] = row["perf"].get("instructions")
            out["runs"].append(row)
            shutil.rmtree(outdir, ignore_errors=True)
        out.update(slope_and_intercept(per_rep, years))
        if len(counts) >= 2:
            y0, y1 = years[0], years[1]
            slope = (counts[y1] - counts[y0]) / (y1 - y0)
            out["perf"] = {"instructions_per_model_year": slope,
                           "fixed_instructions": counts[y0] - slope * y0}
    finally:
        cfg.unlink(missing_ok=True)
    return out


def verify_knowngood(logdir: Path) -> dict:
    """Run the shipped eb_go_gs regression case and compare the result against
    the reference cGENIE ships for it.

    `make testebgogs` cannot answer this on its own for two reasons. Its
    comparison tool is `src/c/nccompare.exe`, which needs the legacy netCDF C++
    header, and its reference file is a GOLDSTEIN annual-average netCDF that
    `goldstein.F:732` only writes when the namelist `debug_loop` is true. That
    flag defaults false, so the shipped test runs, writes no reference-shaped
    output, and reports a failure that is about its own configuration. Setting
    `debug_loop` changes no physics: every use of it in `goldstein.F` guards a
    print, a diagnostic dump or the averaging call, never the state."""
    src = (CONFIG_DIR / "eb_go_gs_test.xml").read_text()
    marker = '\t\t<model name="goldstein">\n'
    assert marker in src
    cfg = CONFIG_DIR / "verify_eb_go_gs.xml"
    cfg.write_text(src.replace(
        marker, marker + '\t\t\t<param name="debug_loop">.true.</param>\n', 1
    ).replace("<var name=\"EXPID\">genie_eb_go_gs</var>",
              f"<var name=\"EXPID\">{cfg.stem}</var>"))
    out = {"reference": "genie-knowngood/genie_eb_go_gs_knowngood/goldstein/"
                        "gold_spn_av_0000000020_00.nc"}
    try:
        run(["/usr/bin/make", *make_args().split(), "cleanall"], logdir / "verify.clean.log")
        rc, secs = run(job(f"configs/{cfg.name}", ["-x"], target="genie.exe"),
                       logdir / "verify.build.log")
        out["build_ok"] = rc == 0
        if rc != 0:
            return out
        log = logdir / "verify.run.log"
        rc, secs = run(job(f"configs/{cfg.name}", ["-z"]), log)
        ok, why = run_ok(log)
        out["run_ok"] = bool(rc == 0 and ok)
        out["run_seconds"] = model_seconds(log)
        if not out["run_ok"]:
            out["error"] = why or f"exit {rc}"
            return out

        import netCDF4  # noqa: PLC0415

        ref = CGENIE / out["reference"]
        got = OUT_ROOT / cfg.stem / "goldstein" / "gold_spn_av_0000000020_00.nc"
        out["produced"] = got.exists()
        if not got.exists():
            return out
        fields = {}
        with netCDF4.Dataset(ref) as a, netCDF4.Dataset(got) as b:
            for name in sorted(set(a.variables) & set(b.variables)):
                x, y = a[name][:], b[name][:]
                if x.shape != y.shape or x.dtype.kind != "f":
                    continue
                d = abs(y - x)
                scale = float(abs(x).max()) or 1.0
                fields[name] = {
                    "max_abs_diff": float(d.max()),
                    "max_rel_to_field_range": float(d.max()) / scale,
                }
            out["variables_only_in_reference"] = sorted(set(a.variables) - set(b.variables))
        out["fields"] = fields
        out["worst_relative"] = max((v["max_rel_to_field_range"] for v in fields.values()),
                                    default=None)
        shutil.rmtree(OUT_ROOT / cfg.stem, ignore_errors=True)
    finally:
        cfg.unlink(missing_ok=True)
    return out


# The stability criterion is applied with the model's OWN instrument.
# `genie-goldstein/src/fortran/diag.f` forms `cnmax` as the largest of
# |u|dt/dx, |v|dt/dy and |w|dt/dz over every wet cell and prints it as `Cn`,
# from the flow GOLDSTEIN diagnosed rather than from a nominal speed. Nothing in
# the component acts on it: there is no CFL or stability test anywhere in
# `genie-goldstein` or `genie-embm`, and `diag` is called only behind
# `debug_loop`, which defaults false. So a run past the limit prints nothing,
# completes, and reports success -- which is why a grid that RUNS is not
# evidence of a grid that is STABLE, and why this sweep turns the diagnostic on.
#
# The criterion, fixed before any of these runs and identical in form to the
# a-priori one in this file's docstring:
#
#   STABLE = the model's own Cn stays below 1 at every diagnostic step, AND the
#            ocean tracer field stays finite and bounded.
#
# The bound on T is a blow-up detector rather than a physical range, so that a
# marginal but real circulation is never called unstable on it alone.
STABILITY_CN_MAX = 1.0
STABILITY_T_ABS_MAX_C = 1.0e3

_CN = re.compile(r"^\s*Cn\s+(\S+)\s*$")
_TK = re.compile(r"^\s*max and min T at kmax\s+(\S+)\s+(\S+)\s*$")
_SST = re.compile(r"^\s*average SST\s+(\S+)\s*$")


def _f(token: str) -> float:
    """Fortran list-directed output of a broken field. gfortran writes NaN and
    Infinity as words; a field too wide for its edit descriptor comes out as
    asterisks. All three mean the same thing here."""
    try:
        return float(token)
    except ValueError:
        return math.nan


def read_diag(log: Path) -> dict:
    sic_fail: dict | None = None
    cn: list[float] = []
    tmax: list[float] = []
    tmin: list[float] = []
    sst: list[float] = []
    for line in log.read_text(errors="replace").split("\n"):
        m = _CN.match(line)
        if m:
            cn.append(_f(m.group(1)))
            continue
        m = _TK.match(line)
        if m:
            tmax.append(_f(m.group(1)))
            tmin.append(_f(m.group(2)))
            continue
        m = _SST.match(line)
        if m:
            sst.append(_f(m.group(1)))
            continue
        m = _SIC_FAIL.search(line)
        if m and sic_fail is None:
            sic_fail = {"step": int(m.group(1)), "i": int(m.group(2)), "j": int(m.group(3))}
    return {"cn": cn, "t_max": tmax, "t_min": tmin, "sst": sst, "sic_fail": sic_fail}


# The one stability guard anywhere in the component, and it is not in the ocean.
# `genie-embm/src/fortran/surflux.F:900-954` solves the sea-ice surface
# temperature by Newton iteration and, when it does not converge, prints this
# and executes a bare `stop` -- unless EMBM's own `debug_loop` is set, in which
# case it clamps `tice` and carries on. So turning EMBM diagnostics on converts
# the only hard guard in the model into a silent clamp. This sweep therefore
# sets `debug_loop` on GOLDSTEIN only: the ocean's `diag` runs and EMBM's guard
# stays live.
_SIC_FAIL = re.compile(r"warning sea-ice iteration failed at\s+(\d+)\s+(\d+)\s+(\d+)")


def stability_verdict(diag: dict) -> tuple[str, str]:
    cn, tmax, tmin = diag["cn"], diag["t_max"], diag["t_min"]
    if diag["sic_fail"]:
        f = diag["sic_fail"]
        return "unstable", ("EMBM's sea-ice surface temperature solve did not converge at"
                            f" step {f['step']}, cell ({f['i']},{f['j']}), and surflux.F"
                            " stopped the model")
    if not cn:
        return "no-diagnostics", "diag never ran; debug_loop or npstp is wrong"
    values = cn + tmax + tmin + diag["sst"]
    if any(not math.isfinite(v) for v in values):
        return "unstable", "the field stopped being finite"
    worst = max(cn)
    if worst >= STABILITY_CN_MAX:
        return "unstable", f"the model's own Cn reached {worst:.3f}"
    excursion = max(max(abs(v) for v in tmax), max(abs(v) for v in tmin))
    if excursion > STABILITY_T_ABS_MAX_C:
        return "unstable", f"surface T reached {excursion:.3g} C"
    return "stable", f"the model's own Cn peaked at {worst:.3f}"


def stability_sweep(name: str, case: dict, nyears: list[int], years: int,
                    maxisles: int, logdir: Path, mcmodel: str = "",
                    ndtas: list[int] | None = None) -> dict:
    """One grid, one executable, several timesteps.

    `nyear` is the only control on the tracer timestep, dt = sodaylen * yearlen
    / nyear, and it is a namelist value, so the whole sweep runs on one build.
    That is the point: it separates what the GRID costs from what the timestep
    the grid needs costs, and those are different halves of a ceiling.

    `ndta` is the second knob and it is not the same knob. EMBM's own step is
    `dtatm = dt_ocean / ndta` (`initialise_embm.F:578`), namelist-reachable in
    `ini_embm_nml`, so the atmosphere can be sub-stepped without shortening the
    ocean's tracer step. Which of the two a grid actually needs decides whether
    a refinement costs the whole model or only EMBM."""
    global MCMODEL  # noqa: PLW0603
    MCMODEL = mcmodel or case.get("mcmodel", "")
    out = {"case": name, "nlons": case["nlons"], "nlats": case["nlats"],
           "nlevs": case["nlevs"], "years": years, "code_model": MCMODEL or "small",
           "synthetic_forcing": case.get("synthetic_forcing", False), "points": []}
    cfg = CONFIG_DIR / f"stab_{name}.xml"
    cfg_rel = f"configs/{cfg.name}"
    try:
        write_config(cfg, case, years, nyears[0], maxisles, npstp=nyears[0],
                     debug_loop=True)
        run(["/usr/bin/make", *make_args().split(), "cleanall"],
            logdir / f"stab_{name}.clean.log")
        rc, secs = run(job(cfg_rel, ["-x"], target="genie.exe"),
                       logdir / f"stab_{name}.build.log")
        out["build_ok"] = rc == 0
        out["build_seconds"] = round(secs, 2)
        if rc != 0:
            out["error"] = f"build failed, see {logdir / ('stab_' + name + '.build.log')}"
            return out
        for nyear, ndta in [(n, d) for n in nyears for d in (ndtas or [5])]:
            # print once per simulated year, whatever the timestep
            write_config(cfg, case, years, nyear, maxisles, npstp=nyear,
                         debug_loop=True, ndta=ndta)
            log = logdir / f"stab_{name}.nyear{nyear}.ndta{ndta}.log"
            rc, secs = run(job(cfg_rel, ["-z"]), log)
            ok, why = run_ok(log)
            diag = read_diag(log)
            state, reason = stability_verdict(diag)
            if not ok and state != "unstable":
                # only when the log says nothing about WHY; a named instability
                # is a better answer than the absence of a shutdown banner
                state, reason = "did-not-complete", why or f"exit {rc}"
            point = {
                "nyear": nyear,
                "ndta": ndta,
                "atmosphere_dt_s": 86400.0 * 365.25 / (nyear * ndta),
                "predicted": cfl(case["nlons"], case["nlats"], case["nlevs"], nyear),
                "completed": ok,
                "model_seconds": model_seconds(log),
                "load_average_after": loadavg(),
                "sea_ice_solve_failed_at": diag["sic_fail"],
                "cn_peak": max(diag["cn"]) if diag["cn"] else None,
                "cn_final": diag["cn"][-1] if diag["cn"] else None,
                "diagnostic_steps": len(diag["cn"]),
                "t_max_c": max(diag["t_max"]) if diag["t_max"] else None,
                "t_min_c": min(diag["t_min"]) if diag["t_min"] else None,
                "verdict": state,
                "because": reason,
            }
            out["points"].append(point)
            print(json.dumps(point, indent=2), flush=True)
            shutil.rmtree(OUT_ROOT / cfg.stem, ignore_errors=True)
    finally:
        cfg.unlink(missing_ok=True)
    return out


# Why the timing arm builds every case FIRST and only then runs them.
#
# `genie.job -z` sets REMAKE=FALSE: no make runs, and the job copies
# `genie-main/genie.exe` into the run directory and executes it there. So an
# executable can be built, stashed, and restored later, which makes it possible
# to run the cases ROUND ROBIN rather than one case at a time.
#
# It matters because the cases are compared with each other. A sweep that
# finishes one grid before starting the next gives each grid a different slice
# of whatever else the host is doing, and a ratio between two of them is then a
# ratio of two machine states as much as of two grids. Round robin gives every
# case the same distribution of conditions, and the minimum over repeats takes
# the least contended sample of each. The absolute numbers remain upper bounds
# under contention; the ratios are what this buys.
def build_all(cases: list[str], args, results: dict) -> Path:
    """Build one executable per case and stash it. Returns the stash directory."""
    global MCMODEL  # noqa: PLW0603
    stash = args.logdir / "exe"
    stash.mkdir(parents=True, exist_ok=True)
    for name in cases:
        case = CASES[name]
        MCMODEL = args.mcmodel or case.get("mcmodel", "")
        record = results[name]
        record["code_model"] = MCMODEL or "small (the compiler default user.mak leaves in place)"
        if (stash / f"{name}.exe").exists() and args.reuse_builds:
            # The executable is a pure function of the grid macros and the code
            # model, so re-timing does not need a rebuild. Said explicitly here
            # because rule 4 is about exactly this: an executable that no longer
            # matches its configuration is silently wrong.
            record["build_ok"] = True
            record["build_reused"] = True
            record["exe_bytes"] = (stash / f"{name}.exe").stat().st_size
            continue
        cfg = CONFIG_DIR / f"bench_{name}.xml"
        try:
            write_config(cfg, case, args.years[0], args.nyear, args.maxisles,
                         ndta=case.get("ndta", 5))
            run(["/usr/bin/make", *make_args().split(), "cleanall"],
                args.logdir / f"{name}.clean.log")
            rc, secs = run(job(f"configs/{cfg.name}", ["-x"], target="genie.exe"),
                           args.logdir / f"{name}.build.log")
            record["build_seconds"] = round(secs, 2)
            record["build_ok"] = rc == 0
            exe = GENIE_MAIN / "genie.exe"
            record["exe_bytes"] = exe.stat().st_size if exe.exists() else None
            if rc != 0 or not exe.exists():
                record["error"] = f"build failed, see {args.logdir / (name + '.build.log')}"
                continue
            shutil.copy2(exe, stash / f"{name}.exe")
        finally:
            cfg.unlink(missing_ok=True)
    return stash


# Above this one-minute load average a wall-clock number is a measurement of the
# host rather than of the model, and is labelled as such rather than dropped: a
# number recorded as unreliable is evidence about the instrument, a deleted one
# is nothing. scripts/machine.py uses the same threshold.
QUIET_LOAD = 4.0


def slope_and_intercept(per_rep: dict[int, dict[int, float]], years: list[int],
                        loads: dict[int, float] | None = None) -> dict:
    """Fit each REPEAT separately and keep the cheapest fit.

    The obvious estimator -- minimum over repeats at each length, then one line
    through the two minima -- is wrong on a busy host, and wrong in a way that
    hides. The two minima can come from different repeats, so their difference
    is not a difference of two runs that saw the same machine, and the slope it
    gives can be anything. A first pass in that form returned a 36 x 36 x 32
    ocean as CHEAPER per model year than the 36 x 36 x 16 one, which cannot be
    true and is the reason this function exists.

    Fitting within a repeat keeps the two lengths adjacent in time, so they saw
    nearly the same machine, and taking the minimum slope over repeats picks the
    least contended of those fits. Every repeat's fit is reported so the scatter
    is visible rather than asserted."""
    y0, y1 = years[0], years[1]
    fits = []
    for rep in sorted(per_rep):
        t = per_rep[rep]
        if y0 not in t or y1 not in t:
            continue
        slope = (t[y1] - t[y0]) / (y1 - y0)
        load = loads.get(rep) if loads else None
        fits.append({"repeat": rep, "seconds_per_model_year": round(slope, 4),
                     "fixed_seconds": round(t[y0] - slope * y0, 3),
                     "model_seconds": {str(y): t[y] for y in years if y in t},
                     "load_average": load,
                     "contaminated": None if load is None else load > QUIET_LOAD})
    if not fits:
        return {}
    clean = [f for f in fits if f["contaminated"] is False]
    best = min(clean or fits, key=lambda f: f["seconds_per_model_year"])
    return {
        "seconds_per_model_year": best["seconds_per_model_year"],
        "fixed_seconds": best["fixed_seconds"],
        "from_repeat": best["repeat"],
        "from_a_quiet_host": bool(clean),
        "load_average_of_that_fit": best["load_average"],
        "per_repeat_fits": fits,
    }


# The cost figure that does not have to be re-taken.
#
# Wall clock on this host is a measurement of a machine state as much as of the
# model: it runs many agents at once, and a run that shares the last-level cache
# and the memory controller with a sixteen-thread integration is slower without
# anything about cGENIE having changed. RETIRED INSTRUCTIONS are not. They are a
# property of the binary and its input, they are the same on a busy host and an
# idle one, and for a strictly serial process there is no spin-wait to correct
# for -- there are no barriers to spin at. So the durable price of an
# ocean-year is quoted here in instructions, and seconds are quoted beside it
# for whoever has to wait.
_PERF_EVENTS = "instructions,task-clock"


def perf_run(outdir: Path, log: Path) -> dict:
    """Re-run an already configured experiment under perf, in its own directory.

    `genie.job` has by this point written the namelists, staged the inputs and
    copied the executable into `outdir`, so running it again there repeats
    exactly the same integration."""
    cmd = ["perf", "stat", "-x,", "-e", _PERF_EVENTS, "./genie.exe"]
    proc = subprocess.run(cmd, cwd=outdir, capture_output=True, text=True, check=False)
    log.write_text(proc.stdout + "\n----- perf -----\n" + proc.stderr)
    out: dict = {}
    for line in proc.stderr.split("\n"):
        parts = line.split(",")
        if len(parts) < 3:
            continue
        value, _, event = parts[0], parts[1], parts[2]
        try:
            out[event.strip()] = float(value)
        except ValueError:
            out[event.strip()] = None
    if "Shutdown complete; home time" not in proc.stdout:
        out["error"] = "the run under perf did not reach the shutdown banner"
    return out


def perf_all(cases: list[str], args, results: dict, stash: Path) -> None:
    """Instructions per ocean-year, per case, from two lengths as before."""
    global MCMODEL  # noqa: PLW0603
    for name in cases:
        case = CASES[name]
        record = results[name]
        if not record.get("build_ok") or record.get("error"):
            continue
        MCMODEL = args.mcmodel or case.get("mcmodel", "")
        cfg = CONFIG_DIR / f"bench_{name}.xml"
        counts: dict[int, float] = {}
        arm: list[dict] = []
        for years in args.years:
            shutil.copy2(stash / f"{name}.exe", GENIE_MAIN / "genie.exe")
            try:
                write_config(cfg, case, years, args.nyear, args.maxisles,
                             ndta=case.get("ndta", 5))
                setup = args.logdir / f"{name}.perfsetup{years}.log"
                rc, _ = run(job(f"configs/{cfg.name}", ["-z"]), setup)
                ok, why = run_ok(setup)
                if rc != 0 or not ok:
                    record["perf_error"] = f"setup {years}y failed ({why or rc})"
                    break
                outdir = OUT_ROOT / cfg.stem
                got = perf_run(outdir, args.logdir / f"{name}.perf{years}.log")
                got["years"] = years
                got["load_average"] = loadavg()
                arm.append(got)
                if got.get("instructions"):
                    counts[years] = got["instructions"]
                shutil.rmtree(outdir, ignore_errors=True)
            finally:
                cfg.unlink(missing_ok=True)
        record["perf"] = {"runs": arm}
        if len(counts) >= 2:
            y0, y1 = args.years[0], args.years[1]
            slope = (counts[y1] - counts[y0]) / (y1 - y0)
            record["perf"]["instructions_per_model_year"] = slope
            record["perf"]["fixed_instructions"] = counts[y0] - slope * y0
            record["perf"]["note"] = (
                "retired instructions are a property of the binary and its input,"
                " not of what else the host was doing")


def time_from_logs(cases: list[str], args, results: dict) -> None:
    """Re-derive the fits from runs that already happened.

    The measurement is the log: `genie.job` wraps the executable in `time`, so
    every run's own wall clock is in the file it wrote. Re-reading them costs no
    machine time and changes nothing about what was measured, which is what
    makes it the right way to correct an ESTIMATOR after the fact rather than
    re-running and quietly getting different numbers."""
    for name in cases:
        record = results[name]
        per_rep: dict[int, dict[int, float]] = {}
        for rep in range(args.repeats):
            for years in args.years:
                log = args.logdir / f"{name}.run{years}.{rep}.log"
                if not log.exists():
                    continue
                ok, _ = run_ok(log)
                inner = model_seconds(log)
                if ok and inner is not None:
                    per_rep.setdefault(rep, {})[years] = inner
        if not per_rep:
            record["error"] = f"no usable run logs under {args.logdir}"
            continue
        record["build_ok"] = True
        record["code_model"] = (args.mcmodel or CASES[name].get("mcmodel", "")
                                or "small (the compiler default user.mak leaves in place)")
        for years in args.years:
            times = [t[years] for t in per_rep.values() if years in t]
            if times:
                record["runs"].append({"years": years, "model_seconds": min(times),
                                       "repeats": len(times)})
        record.update(slope_and_intercept(per_rep, args.years))
        record["loads_unknown_from_logs"] = True


def time_all(cases: list[str], args, results: dict, stash: Path) -> None:
    """Run every built case at every length, round robin, `repeats` times."""
    global MCMODEL  # noqa: PLW0603
    per_rep: dict[str, dict[int, dict[int, float]]] = {}
    rep_load: dict[str, dict[int, float]] = {}
    best: dict[tuple[str, int], float] = {}
    size: dict[tuple[str, int], int] = {}
    load: dict[tuple[str, int], list[float]] = {}
    for rep in range(args.repeats):
        for name in cases:
            case = CASES[name]
            record = results[name]
            if not record.get("build_ok") or record.get("error"):
                continue
            MCMODEL = args.mcmodel or case.get("mcmodel", "")
            cfg = CONFIG_DIR / f"bench_{name}.xml"
            for years in args.years:
                shutil.copy2(stash / f"{name}.exe", GENIE_MAIN / "genie.exe")
                try:
                    write_config(cfg, case, years, args.nyear, args.maxisles,
                                 ndta=case.get("ndta", 5))
                    log = args.logdir / f"{name}.run{years}.{rep}.log"
                    rc, _ = run(job(f"configs/{cfg.name}", ["-z"]), log)
                    ok, why = run_ok(log)
                    if rc != 0 or not ok:
                        record["error"] = f"run {years}y failed ({why or rc}), see {log}"
                        continue
                    inner = model_seconds(log)
                    key = (name, years)
                    if inner is not None:
                        best[key] = inner if key not in best else min(best[key], inner)
                        per_rep.setdefault(name, {}).setdefault(rep, {})[years] = inner
                        rep_load.setdefault(name, {})[rep] = max(
                            rep_load.get(name, {}).get(rep, 0.0), loadavg())
                    load.setdefault(key, []).append(loadavg())
                    outdir = OUT_ROOT / cfg.stem
                    # genie.job copies the executable into the run directory and
                    # it is larger than everything the run writes. Storage means
                    # what the run produced.
                    size[key] = sum(f.stat().st_size for f in outdir.rglob("*")
                                    if f.is_file() and f.name != "genie.exe")
                    shutil.rmtree(outdir, ignore_errors=True)
                finally:
                    cfg.unlink(missing_ok=True)

    for name in cases:
        record = results[name]
        for years in args.years:
            key = (name, years)
            if key in best:
                record["runs"].append({
                    "years": years,
                    "model_seconds": best[key],
                    "repeats": args.repeats,
                    "load_average_after": load.get(key, []),
                    "output_bytes": size.get(key),
                })
        record.update(slope_and_intercept(per_rep.get(name, {}), args.years,
                                          rep_load.get(name, {})))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--case", action="append", choices=sorted(CASES) + sorted(PROBE72),
                    default=None)
    ap.add_argument("--years", type=int, nargs="+", default=[20, 100])
    ap.add_argument("--repeats", type=int, default=2)
    ap.add_argument("--nyear", type=int, default=100)
    ap.add_argument("--maxisles", type=int, default=20)
    ap.add_argument("--out", type=Path, default=PROJECT_ROOT / "analysis" / "cgenie_cost.json")
    ap.add_argument("--logdir", type=Path, default=Path("/tmp/cgenie_cost_logs"))
    ap.add_argument("--biogem", action="store_true",
                    help="also price the shipped BIOGEM regression case at the same grid")
    ap.add_argument("--biogem-years", type=int, nargs="+", default=[5, 25])
    ap.add_argument("--mcmodel", default="", choices=["", "medium", "large"],
                    help="code model to compile and link with; empty is the shipped small")
    ap.add_argument("--probe72", action="store_true",
                    help="add the 72 x 72 x 16 cost probe, whose wind speeds are synthesised")
    ap.add_argument("--verify", action="store_true",
                    help="also run the shipped eb_go_gs case and compare against genie-knowngood")
    ap.add_argument("--stability", action="store_true",
                    help="sweep nyear on one build and read the model's own Cn, instead of timing")
    ap.add_argument("--nyear-sweep", type=int, nargs="+", default=[50, 100, 200, 400])
    ap.add_argument("--stability-years", type=int, default=10)
    ap.add_argument("--perf", action="store_true",
                    help="also count retired instructions per ocean-year, which no host load changes")
    ap.add_argument("--reuse-builds", action="store_true",
                    help="use the executables already stashed under --logdir/exe instead of rebuilding")
    ap.add_argument("--reuse-runs", action="store_true",
                    help="re-derive the fits from run logs already in --logdir; builds nothing")
    ap.add_argument("--ndta-sweep", type=int, nargs="+", default=[5],
                    help="EMBM sub-steps per ocean step; the atmosphere's own timestep knob")
    args = ap.parse_args()

    global MCMODEL  # noqa: PLW0603
    MCMODEL = args.mcmodel

    if args.probe72:
        CASES.update(PROBE72)
        build_probe72_inputs()
    cases = args.case or [c for c in CASES if args.probe72 or c not in PROBE72]
    args.logdir.mkdir(parents=True, exist_ok=True)
    results = []

    if args.stability:
        sweeps = [stability_sweep(n, CASES[n], args.nyear_sweep, args.stability_years,
                                  args.maxisles, args.logdir, args.mcmodel,
                                  args.ndta_sweep)
                  for n in cases]
        payload = {"provenance": provenance(), "criterion": {
            "instrument": "the model's own cnmax, printed as Cn by genie-goldstein/src/fortran/diag.f",
            "cn_max": STABILITY_CN_MAX,
            "t_abs_max_c": STABILITY_T_ABS_MAX_C,
            "stated": "before any run, beside stability_sweep in this file",
        }, "sweeps": sweeps}
        args.out.write_text(json.dumps(payload, indent=2) + "\n")
        print(f"wrote {args.out}")
        return 0

    records = {}
    for name in cases:
        case = CASES[name]
        isles = island_count(case["world"])
        records[name] = {
            "case": name,
            "world": case["world"],
            "nlons": case["nlons"],
            "nlats": case["nlats"],
            "nlevs": case["nlevs"],
            "ocean_cells": case["nlons"] * case["nlats"] * case["nlevs"],
            "surface_cells": case["nlons"] * case["nlats"],
            "ubarsolv_work": case["nlons"] ** 2 * case["nlats"],
            "synthetic_forcing": case.get("synthetic_forcing", False),
            "islands_in_psiles": isles,
            "maxisles_compiled": args.maxisles,
            "nyear": args.nyear,
            "ndta": case.get("ndta", 5),
            "stability": cfl(case["nlons"], case["nlats"], case["nlevs"], args.nyear),
            "runs": [],
        }

    stash = args.logdir / "exe"
    if args.reuse_runs:
        time_from_logs(cases, args, records)
    else:
        stash = build_all(cases, args, records)
        time_all(cases, args, records, stash)
    if args.perf:
        if args.reuse_runs:
            for n in cases:
                records[n]["build_ok"] = (stash / f"{n}.exe").exists()
        perf_all(cases, args, records, stash)
    results = [records[n] for n in cases]
    for record in results:
        print(json.dumps(record, indent=2), flush=True)

    payload = {"provenance": provenance(), "criterion": {
        "advective_courant_max": 1.0,
        "diffusive_number_max": 0.5,
        "u_nominal_ms": U_NOMINAL_MS,
        "stated": "before any run in this file's docstring",
    }, "cases": results,
        "knowngood": verify_knowngood(args.logdir) if args.verify else None,
        "biogem": biogem_cost(args.logdir, args.biogem_years, args.repeats,
                              args.reuse_runs, args.perf)
        if args.biogem else None}
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
