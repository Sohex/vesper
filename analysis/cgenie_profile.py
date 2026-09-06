#!/usr/bin/env python3
"""Profile the EMIC ocean: where the retired instructions actually go.

    python analysis/cgenie_profile.py --years 10 --period 20000000

Writes `analysis/cgenie_profile.json`;
`notes/audits/cgenie-parallelism-and-coupling-support.md` is what it found.

Worldbuilding. Vesper is an invented planet; this profiles a candidate ocean
component on this host, not any real ocean.

## Why this exists and what it is FOR

`notes/audits/cgenie-build-cost-and-grid-ceiling.md` prices a configuration and
says nothing about which routine spends it. Its section 4a settles one binary
question without a profiler -- the barotropic solve is invisible at these grids,
measured 0.956 against a decision rule fixed beforehand -- and explicitly leaves
open where the time goes WITHIN the tracer work.

That open question is not academic. cGENIE is strictly serial: no `!$omp`
directive anywhere in `vendor/cgenie` is live, `-fopenmp` is commented out in
`makefile.arc`, and MPI exists only in `genie-plasim`, which this project does
not use. So the cost of every spin-up is one core against simulated years, and
the coupling resolution has to be argued against that budget. Whether the budget
COULD be different is a function of what fraction of the work sits in loops over
cells and tracers -- which parallelise -- against what sits in a global
recurrence, which does not. This measures that fraction.

## The instrument

Retired instructions, sampled. `perf record -e instructions -c <period>` takes
one sample every `period` retired instructions rather than at a wall-clock
frequency, so the per-symbol shares are a property of the binary and its input
and do not move with what else the host is doing. That matters here for the same
reason it mattered in the cost note: this machine runs several agents at once.

Both grids are compiled `-mcmodel=medium`. The 72 x 72 x 16 grid has no choice
-- its static COMMON is past what `-mcmodel=small` can address -- and building
36 x 36 x 16 the same way means the two profiles differ in the grid and in
nothing else. The cost note prices that flag separately at 1.030 times the
instructions, so it is not free, but it is a uniform 3 per cent and it does not
move between routines.

No source file is edited and no compiler flag beyond the code model is changed.
The tree is exported out of the repository with `git archive` and built there
because the code model IS the change: `-mcmodel=medium` is not what
`vendor/cgenie` is otherwise built at, and `make` decides what to recompile from
timestamps, so medium-model objects left in the tree are what the next build
there silently links against. The export gives this profile a tree it is the
only writer of.

## What a "hot routine" means here

gfortran emits one symbol per external procedure with a trailing underscore, and
these are separate source files, so there is no cross-file inlining to blur the
attribution. `-Ofast` will inline within a file; every routine named below is
the whole of its own file.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "analysis"))

import cgenie_cost as cc  # noqa: E402

# The vendored tree is READ through the repository and BUILT outside it. Both
# grids here compile at -mcmodel=medium, which is not the tree's code model, and
# make recompiles by timestamp, so medium-model objects left in vendor/cgenie
# would be what the next build there links against without saying so.
WORK_ROOT = Path(os.environ.get("CGENIE_Z1", Path.home() / "cgenie_z1"))
CGENIE = WORK_ROOT / "vendor" / "cgenie"
OUT_ROOT = Path(os.environ.get("CGENIE_PROFILE_OUT", Path.home() / "cgenie_profile_out"))
# cgenie_cost builds the 72 x 72 forcing set under its own PROBE_DIR and bakes
# that path into the case's namelist entries at import, so this must be the same
# directory rather than a private one.
PROBE_DIR = cc.PROBE_DIR

CASES = {
    "worjh2_36x36x16": dict(cc.CASES["worjh2_36x36x16"], mcmodel="medium"),
    "dan_72_72x72x16_probe": dict(cc.PROBE72["dan_72_72x72x16_probe"]),
}


def retarget() -> None:
    """Point every path cgenie_cost resolves at the exported tree."""
    cc.CGENIE = CGENIE
    cc.GENIE_MAIN = CGENIE / "genie-main"
    cc.CONFIG_DIR = cc.GENIE_MAIN / "configs"
    cc.OUT_ROOT = OUT_ROOT
    cc.PROBE_DIR = PROBE_DIR
    cc.MCMODEL = "medium"


def sh(cmd: list[str], log: Path, cwd: Path) -> int:
    with log.open("w") as fh:
        return subprocess.call(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT,
                               env=cc.serial_env())


def build(name: str, case: dict, years: int, nyear: int, maxisles: int,
          logdir: Path, stash: Path, reuse: bool) -> dict:
    exe = stash / f"{name}.exe"
    if reuse and exe.exists():
        return {"build_ok": True, "build_reused": True,
                "exe_bytes": exe.stat().st_size, "code_model": "medium"}
    cfg = cc.CONFIG_DIR / f"prof_{name}.xml"
    try:
        cc.write_config(cfg, case, years, nyear, maxisles,
                        ndta=case.get("ndta", 5))
        sh(["/usr/bin/make", *cc.make_args().split(), "cleanall"],
           logdir / f"{name}.clean.log", cwd=cc.GENIE_MAIN)
        t0 = time.perf_counter()
        rc = sh(cc.job(f"configs/{cfg.name}", ["-x"], target="genie.exe"),
                logdir / f"{name}.build.log", cwd=cc.GENIE_MAIN)
        built = cc.GENIE_MAIN / "genie.exe"
        if rc != 0 or not built.exists():
            return {"build_ok": False,
                    "error": f"build failed, see {logdir / (name + '.build.log')}"}
        stash.mkdir(parents=True, exist_ok=True)
        shutil.copy2(built, exe)
        return {"build_ok": True, "build_seconds": round(time.perf_counter() - t0, 1),
                "exe_bytes": exe.stat().st_size, "code_model": "medium"}
    finally:
        cfg.unlink(missing_ok=True)


# Which component a symbol belongs to, resolved from the source tree rather than
# guessed from the name: gfortran's `tstepo_` could be anyone's.
COMPONENTS = {
    "genie-goldstein": "ocean",
    "genie-embm": "atmosphere (EMBM)",
    "genie-goldsteinseaice": "sea ice",
    "genie-main": "driver",
    "genie-lib": "library",
    "genie-biogem": "biogeochemistry",
    "genie-atchem": "biogeochemistry",
    "genie-ents": "land",
}

_DECL = re.compile(
    r"^\s{0,10}(?:recursive\s+)?(?:(?:real|integer|logical|double\s+precision|character)"
    r"[^!\n]*?\s+)?(subroutine|function)\s+([A-Za-z_][A-Za-z_0-9]*)",
    re.IGNORECASE)


def symbol_index() -> dict[str, dict]:
    """Map a lowercase procedure name to the file and component that defines it."""
    index: dict[str, dict] = {}
    for comp_dir, comp in COMPONENTS.items():
        root = CGENIE / comp_dir / "src" / "fortran"
        if not root.is_dir():
            continue
        for path in sorted(root.rglob("*")):
            if path.suffix.lower() not in (".f", ".f90", ".f95"):
                continue
            for line in path.read_text(errors="replace").split("\n"):
                if line[:1] in ("c", "C", "*", "!"):
                    continue
                m = _DECL.match(line)
                if m:
                    index.setdefault(m.group(2).lower(),
                                     {"file": str(path.relative_to(CGENIE)),
                                      "component": comp})
    return index


# THE TRAP THIS GUARD EXISTS FOR, and it cost a whole pass to find.
#
# `kernel.perf_event_max_sample_rate` is 3000 samples per second on this host.
# Above it the kernel THROTTLES: it stops delivering samples for the rest of the
# window and resumes afterwards, and `perf report` says "Total Lost Samples: 0"
# because no RECORD was lost. What is lost is the assumption the profile rests
# on, that every retired instruction is equally likely to be sampled -- the
# surviving samples are the ones that fell in the unthrottled part of each
# window, which is a systematic bias and not noise.
#
# Taken at a period of 2e7 the shipped grid delivers exactly 1409 samples for
# 28.18e9 instructions, so nothing was dropped. Taken at 2e6 the same run
# delivered 7133 of the 14090 the period demands, and the shares moved by up to
# a factor of 1.6 -- eight standard errors on the counts involved. So the tighter
# sampling was not a better measurement of the same thing; it was a measurement
# of a different, throttled process. More samples come from a LONGER run at a
# period the kernel will honour, never from a shorter period.
def sample_rate_ok(period: int, instr_per_second: float = 1.5e10) -> tuple[bool, float, int]:
    """Would this period ask the kernel for more samples than it will deliver?"""
    try:
        cap = int(Path("/proc/sys/kernel/perf_event_max_sample_rate").read_text().strip())
    except OSError:
        return True, 0.0, 0
    wanted = instr_per_second / period
    return wanted <= cap, wanted, cap


def perf_profile(outdir: Path, log: Path, period: int) -> dict:
    """`perf record` the already-configured experiment, then read the flat
    profile back per symbol.

    Sampling on retired instructions rather than on cycles or on a wall-clock
    frequency is what makes the shares independent of host contention: the
    period is counted in the process's own retired instructions."""
    data = outdir / "perf.data"
    rec = subprocess.run(
        ["perf", "record", "-e", "instructions", "-c", str(period),
         "--no-buildid-cache", "-o", str(data), "./genie.exe"],
        cwd=outdir, capture_output=True, text=True, check=False,
        env=cc.serial_env())
    log.write_text(rec.stdout + "\n----- perf record -----\n" + rec.stderr)
    out: dict = {}
    if "Shutdown complete; home time" not in rec.stdout:
        out["error"] = "the run under perf did not reach the shutdown banner"
        return out
    rep = subprocess.run(
        ["perf", "report", "-i", str(data), "--stdio", "--no-children",
         "-g", "none", "-F", "overhead,sample,symbol,dso", "--percent-limit", "0.0"],
        cwd=outdir, capture_output=True, text=True, check=False)
    (log.parent / (log.stem + ".report.txt")).write_text(rep.stdout + rep.stderr)
    rows = []
    total = 0
    for line in rep.stdout.split("\n"):
        if not line or line.startswith("#"):
            continue
        m = re.match(r"\s*([0-9.]+)%\s+(\d+)\s+\[[.k]\]\s+(\S+)\s+(\S+)\s*$", line)
        if not m:
            continue
        rows.append({"percent": float(m.group(1)), "samples": int(m.group(2)),
                     "symbol": m.group(3), "dso": m.group(4)})
        total += int(m.group(2))
    out["samples"] = total
    out["period"] = period
    out["instructions_sampled"] = total * period
    ok, wanted, cap = sample_rate_ok(period)
    out["sample_rate_wanted_per_second"] = round(wanted)
    out["sample_rate_cap_per_second"] = cap
    if not ok:
        out["throttled"] = (
            f"the period asks for about {wanted:.0f} samples per second against a"
            f" kernel cap of {cap}: the kernel dropped samples in bursts and the"
            " shares below are biased, not merely noisy. Lengthen the run instead")
    out["symbols"] = rows
    shutil.rmtree(outdir / "perf.data", ignore_errors=True)
    data.unlink(missing_ok=True)
    return out


# The classification the note's section 0 fixed before any profile was taken.
# A routine is PARALLEL when its loop nest runs over the grid indices and every
# loop-carried dependence in it is a face flux the previous cell already
# computed, recoverable by recomputing one row, column or level per thread, with
# no global reduction other than a sum or a maximum. It is SERIAL when it carries
# a dependence over the whole domain that no such recomputation removes. LEAF is
# a pure routine called from inside one of the parallel nests, which inherits the
# caller's classification and needs nothing of its own. Anything absent from this
# table is UNCLASSIFIED and counts against the parallel fraction, never for it.
CLASSIFICATION = {
    # ocean
    "tstepo_flux_": ("parallel", "k,j,i,l nest; fw/fs/fb carry a face flux forward, recoverable by recomputing one row per thread; limps and dmax are a count and a max"),
    "tstepo_": ("parallel", "copies and boundary loops over the grid; the six ts_t1/ts1_t1/rho_t1/... copies serve commented-out code"),
    "co_": ("parallel", "convective adjustment, a vertical algorithm inside each (i,j)"),
    "krausturner_": ("parallel", "mixed-layer deepening, one column at a time"),
    "velc_": ("parallel", "vertical integration inside each (i,j)"),
    "jbar_": ("parallel", "pressure integral down each column"),
    "wind_": ("parallel", "stress interpolation per cell"),
    "get_hosing_": ("parallel", "per-cell freshwater forcing"),
    "eos_": ("leaf", "pure, no locals, called from the tracer and convection nests"),
    "eosd_": ("leaf", "pure, called once per wet cell from tstepo_flux"),
    "ediff_": ("leaf", "per-cell diffusivity"),
    "goldstein_": ("parallel", "the driver's own (i,j,k) loops: flux assembly, island superposition, boundary copies"),
    "ubarsolv_": ("serial", "banded triangular solve; the forward elimination and back substitution both carry a dependence along the band and no recomputation removes it"),
    "invert_": ("serial", "the one-off LU factorisation behind ubarsolv, at initialisation only"),
    "island_": ("serial", "a path integral accumulated around each island boundary"),
    "matinv_gold_": ("serial", "the isles x isles island system"),
    # biogeochemistry. BIOGEM's own work sits in two loops in biogem.f90: the
    # column sweep `do n=1,n_vocn`, whose every argument is indexed by n with the
    # reduction onto the grid after the loop, and an (i,j) loop beside it. The
    # static checks that would break the first were already done under OCN-19: no
    # SAVE or DATA, no writes to module-scope arrays since the state arrives
    # through dummy arguments, and no initialised locals anywhere in
    # biogem_box.f90.
    "biogem_": ("parallel", "the driver's column sweep and its (i,j) loop"),
    "biogem_tracercoupling_": ("parallel", "do n=1,n_vocn over the vectorised columns"),
    "biogem_climate_": ("parallel", "per-cell climate handoff"),
    "sub_box_remin_part": ("parallel", "one water column, reached from the n_vocn sweep"),
    "sub_box_remin_dom": ("parallel", "one water column, reached from the n_vocn sweep"),
    "sub_box_remin_redfield": ("parallel", "per-column stoichiometry"),
    "sub_calc_bio_uptake": ("parallel", "per-column production"),
    "sub_box_misc_geochem": ("parallel", "per-column geochemistry"),
    "sub_calc_carb": ("leaf", "carbonate system solve, per cell"),
    "sub_calc_carbconst": ("leaf", "carbonate constants, per cell"),
    "fun_lib_conv_vsedtosed": ("parallel", "vector-to-grid conversion over cells"),
    "fun_lib_conv_vocntoocn": ("parallel", "vector-to-grid conversion over cells"),
    "fun_lib_conv_sedtovsed": ("parallel", "grid-to-vector conversion over cells"),
    "fun_lib_conv_ocntovocn": ("parallel", "grid-to-vector conversion over cells"),
    "diag_biogem_timeseries_": ("parallel", "diagnostic accumulation, a sum reduction over cells"),
    "diag_biogem_": ("parallel", "diagnostic accumulation over cells"),
    "cpl_comp_atmocn_": ("parallel", "per-cell composition exchange"),
    "cpl_flux_ocnatm_": ("parallel", "per-cell flux exchange"),
    "atchem_": ("parallel", "per-cell atmospheric chemistry"),
    # atmosphere (EMBM)
    "tstipa_": ("parallel", "Jacobi iteration: tq2 holds the whole previous iterate, so each sweep is data-parallel over cells with a barrier between sweeps"),
    "tstepa_": ("parallel", "the two-dimensional form of tstepo_flux's flux recycling"),
    "surflux_": ("parallel", "one (i,j) loop carrying the sea-ice surface Newton iteration and the land column"),
    "embm_": ("parallel", "the driver's own per-cell loops"),
    "radfor_": ("parallel", "per-latitude insolation"),
    "ocean_alb_": ("leaf", "per-cell albedo"),
    # sea ice
    "tstepsic_": ("parallel", "the same flux recycling on the ice fields"),
    "tstipsic_": ("parallel", "the same Jacobi iteration on the ice fields"),
    "gold_seaice_": ("parallel", "the driver's own per-cell loops"),
}

# libm leaves called from inside the per-cell nests. They are pure and reentrant.
_LIBM_LEAF = {"pow", "log", "exp", "sqrt", "log10", "atan2", "__ieee754_pow_fma",
              "pow@plt", "log@plt", "exp@plt"}
# libgfortran leaves. BIOGEM dispatches tracer behaviour on string names, so these
# are called from inside the per-cell loops; they are pure and reentrant.
_LIBGFORTRAN_LEAF = {"_gfortran_select_string", "_gfortran_compare_string",
                     "_gfortran_string_trim", "_gfortran_concat_string"}
# The heap traffic tstepo_flux's array-section argument creates, one allocation
# per wet cell per timestep. It is not a decomposition question: it is work that
# should not exist. Counted apart so it never flatters either fraction.
_HEAP = {"malloc", "cfree", "free", "malloc@plt", "free@plt", "memcpy@plt",
         "memset@plt", "memmove", "__libc_malloc"}


_MODMANGLE = re.compile(r"^__[A-Za-z_0-9]+_MOD_(.+)$")


def base_symbol(sym: str) -> str:
    """gfortran mangles a module procedure as __<module>_MOD_<name>. Strip that,
    and the trailing underscore an external procedure carries, so one table
    covers both shapes."""
    m = _MODMANGLE.match(sym)
    if m:
        return m.group(1).lower()
    return sym.lower()


def classify(rows: list[dict]) -> dict:
    """Split the profile into what a thread team could divide and what it could not."""
    buckets = {"parallel": 0.0, "serial": 0.0, "heap": 0.0, "unclassified": 0.0}
    for row in rows:
        sym = row["symbol"]
        kind, why = CLASSIFICATION.get(sym, CLASSIFICATION.get(base_symbol(sym), (None, None)))
        if kind in ("parallel", "leaf"):
            bucket = "parallel"
        elif kind == "serial":
            bucket = "serial"
        elif sym in _LIBM_LEAF:
            bucket, why = "parallel", "pure libm leaf of a per-cell nest"
        elif sym in _LIBGFORTRAN_LEAF:
            bucket, why = "parallel", "pure libgfortran leaf of a per-cell nest"
        elif sym in _HEAP or row.get("dso") == "libc.so.6":
            bucket, why = "heap", "heap and byte-moving traffic, mostly the array-section temporaries tstepo_flux creates for eosd"
        else:
            bucket, why = "unclassified", "not in the table; counted against the parallel fraction"
        row["decomposition"] = bucket
        if why:
            row["decomposition_reason"] = why
        buckets[bucket] += row["percent"]
    out = {"shares_percent": {k: round(v, 3) for k, v in buckets.items()}}
    # Amdahl, with the heap and the unclassified counted on the SERIAL side, which
    # is the pessimistic reading and the one a budget should be built on.
    p = buckets["parallel"] / 100.0
    s = 1.0 - p
    out["parallel_fraction"] = round(p, 4)
    out["amdahl_bound"] = {str(n): round(1.0 / (s + p / n), 2)
                           for n in (2, 4, 8, 16, 32)}
    # And the optimistic reading, where only the genuine recurrence is serial.
    p2 = (buckets["parallel"] + buckets["heap"]) / 100.0
    s2 = 1.0 - p2
    out["parallel_fraction_if_heap_removed"] = round(p2, 4)
    out["amdahl_bound_if_heap_removed"] = {str(n): round(1.0 / (s2 + p2 / n), 2)
                                           for n in (2, 4, 8, 16, 32)}
    return out


def annotate(rows: list[dict], index: dict[str, dict]) -> list[dict]:
    for row in rows:
        key = base_symbol(row["symbol"]).rstrip("_")
        meta = index.get(key)
        if meta:
            row.update(meta)
        elif row["dso"] not in ("genie.exe",):
            row["component"] = f"external ({row['dso']})"
        else:
            row["component"] = "unattributed"
    return rows


def profile_case(name: str, case: dict, args, stash: Path) -> dict:
    record: dict = {}
    cfg = cc.CONFIG_DIR / f"prof_{name}.xml"
    try:
        shutil.copy2(stash / f"{name}.exe", cc.GENIE_MAIN / "genie.exe")
        cc.write_config(cfg, case, args.years, args.nyear, args.maxisles,
                        ndta=case.get("ndta", 5))
        setup = args.logdir / f"{name}.setup.log"
        rc = sh(cc.job(f"configs/{cfg.name}", ["-z"]), setup, cwd=cc.GENIE_MAIN)
        ok, why = cc.run_ok(setup)
        if rc != 0 or not ok:
            return {"error": f"setup failed ({why or rc})"}
        outdir = OUT_ROOT / cfg.stem
        record = perf_profile(outdir, args.logdir / f"{name}.perf.log", args.period)
        record["years"] = args.years
        record["load_average"] = cc.loadavg()
        shutil.rmtree(outdir, ignore_errors=True)
    finally:
        cfg.unlink(missing_ok=True)
    return record


def profile_biogem(args, logdir: Path) -> dict:
    """The same profile against the shipped BIOGEM regression case.

    The physics-only arms are EMBM plus GOLDSTEIN plus sea ice with two tracers.
    A spin-up this project would actually run carries the biogeochemistry, which
    the cost note prices at 2.68 times the physics in instructions and 5.4 times
    in seconds at the same grid. So the fraction that matters to a real budget is
    this one's, not the physics-only one's, and it is measured here from
    `configs/eb_go_gs_ac_bg_test.xml` unaltered except for its length."""
    src = (cc.CONFIG_DIR / "eb_go_gs_ac_bg_test.xml").read_text()
    cfg = cc.CONFIG_DIR / "prof_biogem.xml"
    koverall = args.years * 5 * 100
    text = re.sub(r'(<param name="koverall_total">)\d+', rf"\g<1>{koverall}", src)
    text = re.sub(r'(<param name="dt_write">)\d+', rf"\g<1>{koverall}", text)
    text = re.sub(r'(<param name="par_misc_t_runtime">)\d+', rf"\g<1>{args.years}", text)
    text = text.replace("genie_eb_go_gs_ac_bg", cfg.stem, 1)
    out: dict = {"from": "configs/eb_go_gs_ac_bg_test.xml", "grid": "36x36x8",
                 "tracers": 14, "code_model": "medium"}
    try:
        cfg.write_text(text)
        sh(["/usr/bin/make", *cc.make_args().split(), "cleanall"],
           logdir / "biogem.clean.log", cwd=cc.GENIE_MAIN)
        rc = sh(cc.job(f"configs/{cfg.name}", ["-x"], target="genie.exe"),
                logdir / "biogem.build.log", cwd=cc.GENIE_MAIN)
        if rc != 0 or not (cc.GENIE_MAIN / "genie.exe").exists():
            out["error"] = "build failed"
            return out
        out["build_ok"] = True
        setup = logdir / "biogem.setup.log"
        rc = sh(cc.job(f"configs/{cfg.name}", ["-z"]), setup, cwd=cc.GENIE_MAIN)
        ok, why = cc.run_ok(setup)
        if rc != 0 or not ok:
            out["error"] = f"setup failed ({why or rc})"
            return out
        outdir = cc.OUT_ROOT / cfg.stem
        out.update(perf_profile(outdir, logdir / "biogem.perf.log", args.period))
        out["years"] = args.years
        out["load_average"] = cc.loadavg()
        shutil.rmtree(outdir, ignore_errors=True)
    finally:
        cfg.unlink(missing_ok=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--years", type=int, default=10,
                    help="model years per profiled run")
    ap.add_argument("--nyear", type=int, default=100,
                    help="ocean timesteps per model year (the shipped default)")
    ap.add_argument("--maxisles", type=int, default=20)
    ap.add_argument("--period", type=int, default=20_000_000,
                    help="retired instructions between samples")
    ap.add_argument("--case", action="append", default=None)
    ap.add_argument("--biogem", action="store_true",
                    help="profile the shipped BIOGEM regression case as well")
    ap.add_argument("--only-biogem", action="store_true")
    ap.add_argument("--tag", default="",
                    help="suffix for the case key, so an arm at a different run"
                         " length lands beside the others instead of replacing them")
    ap.add_argument("--reclassify", action="store_true",
                    help="re-apply the decomposition table to an existing artifact"
                         " without running anything")
    ap.add_argument("--logdir", type=Path,
                    default=Path.home() / "cgenie_log" / "profile")
    ap.add_argument("--reuse-builds", action="store_true")
    ap.add_argument("--out", type=Path,
                    default=PROJECT_ROOT / "analysis" / "cgenie_profile.json")
    args = ap.parse_args()

    retarget()
    if args.reclassify:
        results = json.loads(args.out.read_text())
        index = symbol_index()
        for name, rec in results.get("cases", {}).items():
            if rec.get("symbols"):
                annotate(rec["symbols"], index)
                rec["decomposition"] = classify(rec["symbols"])
        args.out.write_text(json.dumps(results, indent=2) + "\n")
        for name, rec in results.get("cases", {}).items():
            dec = rec.get("decomposition")
            if dec:
                print(f"== {name} ({rec.get('grid')}) samples={rec.get('samples')}")
                print("   shares:", dec["shares_percent"])
                print("   Amdahl:", dec["amdahl_bound"])
                print("   Amdahl if heap removed:", dec["amdahl_bound_if_heap_removed"])
        return 0
    args.logdir.mkdir(parents=True, exist_ok=True)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    cc.build_probe72_inputs()

    names = [] if args.only_biogem else (args.case or list(CASES))
    stash = args.logdir / "exe"
    index = symbol_index()
    # Arms are run separately on a contended host, so a later arm merges into
    # what an earlier one wrote rather than replacing it.
    results: dict = {"provenance": cc.provenance(), "cases": {}}
    if args.out.exists():
        try:
            results = json.loads(args.out.read_text())
            results.setdefault("cases", {})
        except json.JSONDecodeError:
            pass
    results["provenance"]["profile_instrument"] = (
        "perf record -e instructions -c PERIOD: one sample per PERIOD retired"
        " instructions, so per-symbol shares do not move with host load")
    results["provenance"]["cgenie_tree"] = str(CGENIE)
    results["provenance"]["code_model"] = "medium for both grids, so the two profiles differ only in the grid"

    for name in names:
        case = CASES[name]
        rec = build(name, case, args.years, args.nyear, args.maxisles,
                    args.logdir, stash, args.reuse_builds)
        if rec.get("build_ok"):
            rec.update(profile_case(name, case, args, stash))
            if rec.get("symbols"):
                annotate(rec["symbols"], index)
                by_comp: dict[str, float] = {}
                for row in rec["symbols"]:
                    by_comp[row["component"]] = by_comp.get(row["component"], 0.0) + row["percent"]
                rec["by_component"] = dict(sorted(by_comp.items(),
                                                  key=lambda kv: -kv[1]))
                rec["decomposition"] = classify(rec["symbols"])
        rec["grid"] = f"{case['nlons']}x{case['nlats']}x{case['nlevs']}"
        results["cases"][name + args.tag] = rec

    if args.biogem or args.only_biogem:
        rec = profile_biogem(args, args.logdir)
        if rec.get("symbols"):
            annotate(rec["symbols"], index)
            by_comp: dict[str, float] = {}
            for row in rec["symbols"]:
                by_comp[row["component"]] = by_comp.get(row["component"], 0.0) + row["percent"]
            rec["by_component"] = dict(sorted(by_comp.items(), key=lambda kv: -kv[1]))
            rec["decomposition"] = classify(rec["symbols"])
        results["cases"]["eb_go_gs_ac_bg_36x36x8" + args.tag] = rec

    args.out.write_text(json.dumps(results, indent=2, sort_keys=False) + "\n")
    print(f"wrote {args.out}")
    for name, rec in results["cases"].items():
        print(f"\n== {name} ({rec.get('grid')})  samples={rec.get('samples')}")
        if rec.get("error"):
            print("   ERROR:", rec["error"])
        if rec.get("throttled"):
            print("   THROTTLED:", rec["throttled"])
        for row in (rec.get("symbols") or [])[:15]:
            print(f"   {row['percent']:6.2f}%  {row['symbol']:<24} "
                  f"{row.get('decomposition',''):<14} {row.get('component','')}")
        dec = rec.get("decomposition")
        if dec:
            print("   shares:", dec["shares_percent"])
            print("   Amdahl:", dec["amdahl_bound"],
                  " if heap removed:", dec["amdahl_bound_if_heap_removed"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
