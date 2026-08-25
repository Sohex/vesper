#!/usr/bin/env python3
"""Build cGENIE through a config, and measure what one EMIC ocean-year costs.

    python analysis/cgenie_cost.py --years 20 100 --repeats 2
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
then the model run at two lengths. Two lengths are the point. cGENIE's
initialisation reads the topography, the wind stress fields and the restoring
climatologies before it integrates anything, and at the shipped grid that fixed
cost is a large fraction of a short run. Fitting wall clock against simulated
years separates the per-year slope, which is what a spin-up costs, from the
intercept, which is what every pass pays once.

The component is serial: `makefile.arc`'s gfortran section carries `-fopenmp`
commented out, so the OpenMP directives in `tstepo.F`, `biogem.f90` and
`genie.F` are inert in this build and cores cannot absorb a grid refinement.

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

# Every case is a topography that ships complete: .k1, .paths, .psiles, the
# four wind-stress components, the two wind-speed fields and the restoring
# climatologies, all at the stated grid. Grids are (nlons, nlats, nlevs).
CASES = {
    "worbe2_36x36x8": dict(world="worbe2", nlons=36, nlats=36, nlevs=8, prefix=""),
    "worbe2_36x36x16": dict(world="worbe2", nlons=36, nlats=36, nlevs=16, prefix=""),
    "worbe2_36x36x32": dict(world="worbe2", nlons=36, nlats=36, nlevs=32, prefix=""),
    "g3660l_36x60x16": dict(world="g3660l", nlons=36, nlats=60, nlevs=16, prefix="g3660l_"),
    "igcmv3_64x32x16": dict(world="igcmv3", nlons=64, nlats=32, nlevs=16, prefix="igcmv3_"),
    "igcmv3_64x32x32": dict(world="igcmv3", nlons=64, nlats=32, nlevs=32, prefix="igcmv3_"),
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
\t\t\t<param name="npstp">{big}</param>
\t\t\t<param name="ianav">{big}</param>
\t\t\t<param name="itstp">{big}</param>
\t\t\t<param name="iwstp">{big}</param>
\t\t\t<param name="temp0">10.</param>
\t\t\t<param name="temp1">10.</param>
\t\t\t<param name="tdatafile">{prefix}tempann.silo</param>
\t\t\t<param name="sdatafile">{prefix}saliann.silo</param>
\t\t\t<paramArray name="diff">
\t\t\t\t<param index="1">{diff_h}</param>
\t\t\t\t<param index="2">{diff_v}</param>
\t\t\t</paramArray>
\t\t</model>
\t\t<model name="goldsteinseaice">
\t\t\t<param name="world">{world}</param>
\t\t\t<param name="nyear">{nyear}</param>
\t\t\t<param name="npstp">{big}</param>
\t\t\t<param name="ianav">{big}</param>
\t\t\t<param name="itstp">{big}</param>
\t\t\t<param name="iwstp">{big}</param>
\t\t</model>
\t\t<model name="embm">
\t\t\t<param name="world">{world}</param>
\t\t\t<param name="nyear">{nyear}</param>
\t\t\t<param name="npstp">{big}</param>
\t\t\t<param name="ianav">{big}</param>
\t\t\t<param name="itstp">{big}</param>
\t\t\t<param name="iwstp">{big}</param>
\t\t\t<param name="tatm">0.0</param>
\t\t\t<param name="xu_wstress">{prefix}taux_u.interp</param>
\t\t\t<param name="yu_wstress">{prefix}tauy_u.interp</param>
\t\t\t<param name="xv_wstress">{prefix}taux_v.interp</param>
\t\t\t<param name="yv_wstress">{prefix}tauy_v.interp</param>
\t\t\t<param name="u_wspeed">{prefix}uncep.silo</param>
\t\t\t<param name="v_wspeed">{prefix}vncep.silo</param>
\t\t\t<param name="tdatafile">{prefix}ta_ncep.silo</param>
\t\t\t<param name="qdatafile">{prefix}qa_ncep.silo</param>
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


def write_config(path: Path, case: dict, years: int, nyear: int, maxisles: int) -> None:
    kloop = 5
    koverall = years * kloop * nyear
    path.write_text(
        CONFIG_TEMPLATE.format(
            expid=path.stem,
            world=case["world"],
            prefix=case["prefix"],
            nlons=case["nlons"],
            nlats=case["nlats"],
            nlevs=case["nlevs"],
            nyear=nyear,
            kloop=kloop,
            koverall=koverall,
            big=koverall + 1,
            genie_timestep=f"{86400.0 * 365.25 / (kloop * nyear):.4f}",
            diff_h=DIFF_H_M2S,
            diff_v=DIFF_V_M2S,
            maxisles=maxisles,
        )
    )


def run(cmd: list[str], log: Path, cwd: Path = GENIE_MAIN) -> tuple[int, float]:
    start = time.perf_counter()
    with log.open("w") as fh:
        rc = subprocess.call(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)
    return rc, time.perf_counter() - start


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
        "cpu": cpu,
        "kernel": platform.release(),
        "build_type": "SHIP (user.mak default), gfortran, serial: -fopenmp is commented out in makefile.arc",
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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--case", action="append", choices=sorted(CASES), default=None)
    ap.add_argument("--years", type=int, nargs="+", default=[20, 100])
    ap.add_argument("--repeats", type=int, default=2)
    ap.add_argument("--nyear", type=int, default=100)
    ap.add_argument("--maxisles", type=int, default=20)
    ap.add_argument("--out", type=Path, default=PROJECT_ROOT / "analysis" / "cgenie_cost.json")
    ap.add_argument("--logdir", type=Path, default=Path("/tmp/cgenie_cost_logs"))
    args = ap.parse_args()

    cases = args.case or list(CASES)
    args.logdir.mkdir(parents=True, exist_ok=True)
    results = []

    for name in cases:
        case = CASES[name]
        isles = island_count(case["world"])
        record = {
            "case": name,
            "world": case["world"],
            "nlons": case["nlons"],
            "nlats": case["nlats"],
            "nlevs": case["nlevs"],
            "ocean_cells": case["nlons"] * case["nlats"] * case["nlevs"],
            "surface_cells": case["nlons"] * case["nlats"],
            "islands_in_psiles": isles,
            "maxisles_compiled": args.maxisles,
            "nyear": args.nyear,
            "stability": cfl(case["nlons"], case["nlats"], case["nlevs"], args.nyear),
            "runs": [],
        }

        cfg = CONFIG_DIR / f"bench_{name}.xml"
        cfg_rel = f"configs/{cfg.name}"
        try:
            write_config(cfg, case, args.years[0], args.nyear, args.maxisles)
            rc, _ = run(["/usr/bin/make", *make_args().split(), "cleanall"],
                        args.logdir / f"{name}.clean.log")
            rc, secs = run(job(cfg_rel, ["-x"], target="genie.exe"),
                           args.logdir / f"{name}.build.log")
            record["build_seconds"] = round(secs, 2)
            record["build_ok"] = rc == 0
            exe = GENIE_MAIN / "genie.exe"
            record["exe_bytes"] = exe.stat().st_size if exe.exists() else None
            if rc != 0:
                record["error"] = f"build failed, see {args.logdir / (name + '.build.log')}"
                results.append(record)
                continue

            for years in args.years:
                write_config(cfg, case, years, args.nyear, args.maxisles)
                best = None
                best_model = None
                for rep in range(args.repeats):
                    log = args.logdir / f"{name}.run{years}.{rep}.log"
                    rc, secs = run(job(cfg_rel, ["-z"]), log)
                    if rc != 0:
                        record["error"] = f"run {years}y failed, see {log}"
                        best = None
                        break
                    best = secs if best is None else min(best, secs)
                    inner = model_seconds(log)
                    if inner is not None:
                        best_model = inner if best_model is None else min(best_model, inner)
                if best is None:
                    break
                outdir = OUT_ROOT / cfg.stem
                size = sum(f.stat().st_size for f in outdir.rglob("*") if f.is_file())
                record["runs"].append({
                    "years": years,
                    "model_seconds": best_model,
                    "wall_seconds": round(best, 3),
                    "repeats": args.repeats,
                    "output_bytes": size,
                })
                shutil.rmtree(outdir, ignore_errors=True)
        finally:
            cfg.unlink(missing_ok=True)

        if len(record["runs"]) >= 2:
            pairs = [(r["years"], r["model_seconds"] or r["wall_seconds"]) for r in record["runs"][:2]]
            (y0, t0), (y1, t1) = pairs
            slope = (t1 - t0) / (y1 - y0)
            record["seconds_per_model_year"] = round(slope, 4)
            record["fixed_seconds"] = round(t0 - slope * y0, 3)
        results.append(record)
        print(json.dumps(record, indent=2), flush=True)

    payload = {"provenance": provenance(), "criterion": {
        "advective_courant_max": 1.0,
        "diffusive_number_max": 0.5,
        "u_nominal_ms": U_NOMINAL_MS,
        "stated": "before any run in this file's docstring",
    }, "cases": results}
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
