#!/usr/bin/env python3
"""Build, verify and time the EMIC ocean's optimisation arms.

    python analysis/cgenie_omp.py --export
    python analysis/cgenie_omp.py --arm base --arm eosd --verify --instructions

Worldbuilding. Vesper is an invented planet; this drives a candidate ocean
component on this host and measures its cost. Ocean, sea ice and tracer
transport name modelled quantities, not observations of anything.

## What this is FOR

`notes/audits/cgenie-parallelism-and-coupling-support.md` measured where the
ocean's instructions go and left three changes proposed and none made. This
driver is what makes them, one ARM at a time, and holds every arm to the same
acceptance test.

An ARM is a named set of edits to the exported tree plus the compiler flags it
is built with. `base` is the tree as it stands. Every other arm has to
reproduce `base`'s answer.

## The acceptance test, fixed before any arm was built

`genie-knowngood/` ships a GOLDSTEIN annual average for the eb_go_gs regression
case. Each of these changes is correctness-preserving by construction -- an
array temporary removed, a storage class changed, a loop divided over threads
that carries no dependence -- so the bar is not a tolerance:

  EVERY float variable in the regression case's `gold_spn_av` output must be
  BIT-FOR-BIT identical to what the `base` arm writes.

A difference is a defect to explain, not a tolerance to widen. The one place
this driver will accept less is a threaded arm whose only difference is in a
routine holding a floating-point REDUCTION, because summation order is a
property of the thread count; no such reduction has been threaded here, so the
allowance is declared and unused.

The reference the arms are compared against is `base`'s own output rather than
the shipped netCDF, because the shipped file was written by a different
compiler on a different host and the difference between it and `base` is not a
property of any arm. `--verify` reports both.

## The instrument

Retired instructions, counted rather than sampled: `perf stat -e instructions`.
On a single-threaded process that is exact and needs no correction. On a
THREADED process it is not a cost: a thread spinning at an OpenMP barrier
retires instructions in proportion to how long it waits, so a threaded arm is
reported with `OMP_WAIT_POLICY=passive` and with the wall clock beside it, and
the instruction count is labelled as what it is.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "analysis"))

import cgenie_cost as cc  # noqa: E402

# A worktree's ignored build products are SYMLINKS into the shared checkout, so
# building in place would write there. Everything happens in an exported copy.
WORK_ROOT = Path(os.environ.get("CGENIE_Z3", Path.home() / "cgenie_z3"))
CGENIE = WORK_ROOT / "vendor" / "cgenie"
OUT_ROOT = Path(os.environ.get("CGENIE_OMP_OUT", Path.home() / "cgenie_omp_out"))
LOGDIR = Path(os.environ.get("CGENIE_OMP_LOG", Path.home() / "cgenie_log" / "omp"))


def retarget() -> None:
    cc.CGENIE = CGENIE
    cc.GENIE_MAIN = CGENIE / "genie-main"
    cc.CONFIG_DIR = cc.GENIE_MAIN / "configs"
    cc.OUT_ROOT = OUT_ROOT
    cc.MCMODEL = "medium"


def export_tree(rev: str | None = None) -> None:
    """Copy the tracked vendored tree out of the repository and into WORK_ROOT.

    With no revision it takes the WORKING-TREE state, so uncommitted edits are
    what gets built. With one it takes that revision, which is how the arm a
    change is measured against is built from the tree as it stood before it.

    Building in place is not an option: a worktree's ignored build products are
    symlinks into the shared checkout, so `make` there writes through them."""
    if WORK_ROOT.exists():
        shutil.rmtree(WORK_ROOT)
    WORK_ROOT.mkdir(parents=True)
    if rev:
        src = subprocess.Popen(["git", "archive", rev, "vendor/cgenie"],
                               cwd=PROJECT_ROOT, stdout=subprocess.PIPE)
        dst = subprocess.Popen(["tar", "-xf", "-", "-C", str(WORK_ROOT)],
                               stdin=src.stdout)
        src.stdout.close()
        src.wait()
        dst.wait()
        print(f"exported vendor/cgenie at {rev} to {WORK_ROOT}")
        return
    names = subprocess.check_output(
        ["git", "ls-files", "-z", "vendor/cgenie"], cwd=PROJECT_ROOT).split(b"\0")
    src = subprocess.Popen(["tar", "--null", "-T", "-", "-cf", "-"],
                           cwd=PROJECT_ROOT, stdin=subprocess.PIPE,
                           stdout=subprocess.PIPE)
    dst = subprocess.Popen(["tar", "-xf", "-", "-C", str(WORK_ROOT)],
                           stdin=src.stdout)
    src.stdout.close()
    src.stdin.write(b"\0".join(names))
    src.stdin.close()
    src.wait()
    dst.wait()
    print(f"exported {len(names)} working-tree paths to {WORK_ROOT}")


def sh(cmd: list[str], log: Path, cwd: Path, env: dict | None = None) -> int:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as fh:
        return subprocess.call(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT,
                               env=env)


# ---------------------------------------------------------------- the arms
#
# Each arm is (extra compiler flags, whether the tree is expected to be edited).
# The EDITS are made in the repository and travel through --export; an arm that
# names a flag only is a build-flag arm and needs no edit.
ARMS = {
    # The tree as it stands, with whatever `makefile.arc` chooses.
    "base": dict(fflags=[]),
    # THE READ-BEFORE-WRITE DETECTOR, and it is here because dropping
    # -fno-automatic changes more than where a local LIVES.
    #
    # A static local sits in .bss and the loader zeroes it, so a routine that
    # reads one of its own locals before writing it reads zero on the first
    # call. An automatic local is whatever was on the stack. That is a
    # behaviour change no regression test is guaranteed to reach, because the
    # value it changes is one nothing was supposed to read.
    #
    # -finit-real=snan -finit-integer=-2147483647 makes every uninitialised
    # local a value that cannot pass for a plausible answer, so a read before
    # write propagates instead of hiding. The test that can FAIL is that this
    # arm reproduces `base` bit-for-bit: if it does, nothing in the exercised
    # paths read a local it had not written.
    "initpoison": dict(fflags=["-finit-real=snan",
                               "-finit-integer=-2147483647"]),
    # The tree compiled SERIALLY now that makefile.arc carries -fopenmp:
    # the !$OMP sentinels go back to being comments. This is the arm that
    # separates what the source changes cost from what the threads buy.
    "serial": dict(fflags=["-fno-openmp"]),
    # The tree compiled serially AND with the storage class it had before,
    # so an instruction count against `serial` prices -frecursive alone.
    "serialstatic": dict(fflags=["-fno-openmp", "-fno-automatic"]),
    # The threaded build. Without -fopenmp the !$OMP sentinels are comments,
    # so `base` is the same source compiled serially -- which is what makes
    # "the restructuring alone changed nothing" a separate, answerable
    # question from "the threads changed nothing".
    "omp": dict(fflags=["-fopenmp"]),
    # -fopenmp and the poison together: an OpenMP private copy is
    # uninitialised on entry to the region exactly as an automatic local is,
    # so this is the arm that catches a variable that needed the value it had
    # outside.
    "omppoison": dict(fflags=["-fopenmp", "-finit-real=snan",
                              "-finit-integer=-2147483647"]),
}

# The regression CASES the acceptance test runs. Each names the shipped
# configuration, the shipped reference tree, and the netCDF inside both that is
# compared variable by variable.
#
# `eb_go_gs` is the physics: EMBM, GOLDSTEIN and the sea ice. It is the arm
# every change here touches directly.
#
# `eb_go_gs_ac_bg` carries ATCHEM and BIOGEM as well, and it is in the list for
# one specific reason rather than for coverage: it is the configuration whose
# per-column arrays in biogem.f90 are what gfortran's own -fmax-stack-var-size
# warning is about, so it is the case that would fail if -frecursive moved more
# onto the stack than the stack holds.
CASES = {
    "eb_go_gs": dict(
        config="eb_go_gs_test.xml",
        expid="genie_eb_go_gs",
        reference="genie-knowngood/genie_eb_go_gs_knowngood",
        nc="goldstein/gold_spn_av_0000000020_00.nc",
        # `goldstein.F:732` writes that annual average only when the namelist
        # `debug_loop` is true, and it defaults false, which is why `make
        # testebgogs` runs, writes no reference-shaped output and reports a
        # failure about its own configuration. Setting it changes no state:
        # every use of it in `goldstein.F` guards a print, a dump or the
        # averaging call.
        debug_loop=True,
    ),
    "eb_go_gs_ac_bg": dict(
        config="eb_go_gs_ac_bg_test.xml",
        expid="genie_eb_go_gs_ac_bg",
        reference="genie-knowngood/genie_eb_go_gs_ac_bg_knowngood",
        nc="biogem/fields_biogem_3d.nc",
        debug_loop=False,
    ),
}


def make_args(target: str = "") -> str:
    """The make overrides this host needs. `user.mak` expects the tree at
    ~/cgenie.muffin and netCDF under /usr/local, and neither is true here."""
    parts = [target, f"GENIE_ROOT={CGENIE}", f"RUNTIME_ROOT={CGENIE}",
             f"NETCDF_DIR={cc.NETCDF_DIR}", f"OUT_DIR={OUT_ROOT}", "-j 1"]
    return " ".join(p for p in parts if p)


def arm_env(arm: str, threads: int = 1) -> dict:
    """An arm's compiler flags travel in the ENVIRONMENT, not on make's command
    line. `genie.job` interpolates its -m argument unquoted, so a make variable
    whose value contains a space is resplit by the shell and `-finit-real=snan`
    arrives as make's own `-f init-real=snan`. Make reads an environment
    variable as a make variable, and `makefile.arc` does not override
    GENIE_FFLAGS, so the environment is where a multi-flag arm belongs.

    -mcmodel=medium is in every arm rather than only where it is needed: the
    doubled grid's static COMMON is past what the small model can address, and
    building every arm the same way means two arms differ in the arm and in
    nothing else."""
    flags = " ".join(["-mcmodel=medium", *ARMS[arm]["fflags"]])
    env = dict(os.environ)
    env["GENIE_FFLAGS"] = flags
    env["GENIE_LDFLAGS"] = flags
    env["OMP_NUM_THREADS"] = str(threads)
    # A thread spinning at a barrier retires instructions in proportion to how
    # long it waits, so an active wait policy puts the idle threads' spin INTO
    # the instruction count. Passive is what makes the count comparable across
    # thread counts.
    env["OMP_WAIT_POLICY"] = "passive"
    return env


def knowngood_config(case: str) -> Path:
    """The shipped regression configuration, renamed so its output lands in a
    directory of this driver's own, and with `debug_loop` turned on where the
    reference the case compares against is a file only `debug_loop` writes."""
    spec = CASES[case]
    src = (cc.CONFIG_DIR / spec["config"]).read_text()
    cfg = cc.CONFIG_DIR / f"omp_{case}.xml"
    if spec["debug_loop"]:
        marker = '\t\t<model name="goldstein">\n'
        assert marker in src, f"{spec['config']} no longer has a goldstein block"
        src = src.replace(
            marker, marker + '\t\t\t<param name="debug_loop">.true.</param>\n', 1)
    src = src.replace(f'<var name="EXPID">{spec["expid"]}</var>',
                      f'<var name="EXPID">{cfg.stem}</var>')
    assert cfg.stem in src, f"{spec['config']} does not set EXPID to {spec['expid']}"
    cfg.write_text(src)
    return cfg


def build(arm: str, case: str) -> dict:
    cfg = knowngood_config(case)
    try:
        env = arm_env(arm)
        sh(["/usr/bin/make", *make_args().split(), "cleanall"],
           LOGDIR / f"{arm}.{case}.clean.log", cwd=cc.GENIE_MAIN, env=env)
        t0 = time.perf_counter()
        rc = sh(["./genie.job", "-x", "-f", f"configs/{cfg.name}",
                 "-o", str(OUT_ROOT), "-c", str(CGENIE), "-g", str(CGENIE),
                 "-h", ".", "-m", make_args(target="genie.exe")],
                LOGDIR / f"{arm}.{case}.build.log", cwd=cc.GENIE_MAIN, env=env)
        built = cc.GENIE_MAIN / "genie.exe"
        out = {"arm": arm, "flags": ARMS[arm]["fflags"],
               "build_seconds": round(time.perf_counter() - t0, 1)}
        if rc != 0 or not built.exists():
            out["build_ok"] = False
            out["error"] = f"see {LOGDIR / (arm + '.' + case + '.build.log')}"
            return out
        stash = LOGDIR / "exe"
        stash.mkdir(parents=True, exist_ok=True)
        shutil.copy2(built, stash / f"{arm}.{case}.exe")
        out["build_ok"] = True
        out["case"] = case
        out["exe"] = str(stash / f"{arm}.{case}.exe")
        return out
    finally:
        cfg.unlink(missing_ok=True)


def run_case(arm: str, case: str, threads: int, tag: str) -> dict:
    """Run the regression case from an already-built, stashed executable.

    `genie.job -z` means REMAKE=FALSE, so the executable sitting in genie-main
    is the one that runs; the arm's stashed copy is put there first, which is
    what keeps the arms from silently sharing a binary."""
    exe = LOGDIR / "exe" / f"{arm}.{case}.exe"
    if not exe.exists():
        return {"run_ok": False, "why": f"no stashed executable for arm {arm}"}
    shutil.copy2(exe, cc.GENIE_MAIN / "genie.exe")
    cfg = knowngood_config(case)
    env = arm_env(arm, threads)
    outdir = OUT_ROOT / cfg.stem
    shutil.rmtree(outdir, ignore_errors=True)
    try:
        rc = sh(["./genie.job", "-z", "-f", f"configs/{cfg.name}",
                 "-o", str(OUT_ROOT), "-c", str(CGENIE), "-g", str(CGENIE),
                 "-h", ".", "-m", make_args()],
                LOGDIR / f"{tag}.run.log", cwd=cc.GENIE_MAIN, env=env)
        log = LOGDIR / f"{tag}.run.log"
        ok, why = cc.run_ok(log)
        rec = {"arm": arm, "case": case, "threads": threads,
               "run_ok": bool(rc == 0 and ok),
               "seconds": cc.model_seconds(log)}
        if why:
            rec["why"] = why
        if rec["run_ok"]:
            keep = LOGDIR / "nc" / tag
            keep.parent.mkdir(parents=True, exist_ok=True)
            shutil.rmtree(keep, ignore_errors=True)
            shutil.copytree(outdir, keep)
            rec["output"] = str(keep)
        return rec
    finally:
        cfg.unlink(missing_ok=True)


# The COST case. The regression cases above answer whether an arm is
# correct; neither answers what it costs, because both are twenty timesteps
# and dominated by startup. This is the audit's own arm: the shipped
# worjh2 topography at 36 x 36 x 16, run the length the audit ran it, so a
# before and an after here are comparable with the shares in
# `notes/audits/cgenie-parallelism-and-coupling-support.md` section 2c.
COST_CASE = "worjh2_36x36x16"


def cost_config(years: int, nyear: int, maxisles: int) -> Path:
    case = dict(cc.CASES[COST_CASE], mcmodel="medium")
    cfg = cc.CONFIG_DIR / f"omp_cost_{COST_CASE}.xml"
    cc.write_config(cfg, case, years, nyear, maxisles, ndta=case.get("ndta", 5))
    return cfg


def cost_run(arm: str, threads: int, tag: str, years: int, nyear: int,
             maxisles: int, build_it: bool) -> dict:
    """Build and run the cost case, counting retired instructions.

    `perf stat -e instructions` COUNTS rather than samples, so on a
    single-threaded process it is exact and needs no correction for host
    load: it is a property of the binary and its input. On a THREADED
    process it is not a cost, because a thread waiting at a barrier retires
    instructions in proportion to how long it waits. OMP_WAIT_POLICY is
    passive in `arm_env` for exactly that reason, and the wall clock and the
    load average are recorded beside the count rather than in place of it,
    so a reader can see what the machine was doing."""
    env = arm_env(arm, threads)
    cfg = cost_config(years, nyear, maxisles)
    rec = {"arm": arm, "threads": threads, "case": COST_CASE,
           "years": years, "nyear": nyear, "tree": str(CGENIE)}
    try:
        if build_it:
            sh(["/usr/bin/make", *make_args().split(), "cleanall"],
               LOGDIR / f"{tag}.clean.log", cwd=cc.GENIE_MAIN, env=env)
            rc = sh(["./genie.job", "-x", "-f", f"configs/{cfg.name}",
                     "-o", str(OUT_ROOT), "-c", str(CGENIE), "-g", str(CGENIE),
                     "-h", ".", "-m", make_args(target="genie.exe")],
                    LOGDIR / f"{tag}.build.log", cwd=cc.GENIE_MAIN, env=env)
            if rc != 0 or not (cc.GENIE_MAIN / "genie.exe").exists():
                rec["error"] = f"build failed, see {LOGDIR / (tag + '.build.log')}"
                return rec
        # Lay the run directory out once with genie.job, then run the
        # executable under perf inside it, so perf sees the model and not
        # the xsltproc configuration around it.
        rc = sh(["./genie.job", "-z", "-f", f"configs/{cfg.name}",
                 "-o", str(OUT_ROOT), "-c", str(CGENIE), "-g", str(CGENIE),
                 "-h", ".", "-m", make_args()],
                LOGDIR / f"{tag}.setup.log", cwd=cc.GENIE_MAIN, env=env)
        ok, why = cc.run_ok(LOGDIR / f"{tag}.setup.log")
        if not ok:
            rec["error"] = why or f"setup run exit {rc}"
            return rec
        rundir = OUT_ROOT / cfg.stem
        rec["load_before"] = cc.loadavg()
        t0 = time.perf_counter()
        res = subprocess.run(
            ["perf", "stat", "-e", "instructions,task-clock", "-x,",
             "./genie.exe"],
            cwd=rundir, capture_output=True, text=True, check=False, env=env)
        rec["wall_seconds"] = round(time.perf_counter() - t0, 2)
        rec["load_after"] = cc.loadavg()
        (LOGDIR / f"{tag}.perf.log").write_text(res.stdout + "\n--- perf ---\n"
                                                + res.stderr)
        if "Shutdown complete; home time" not in res.stdout:
            rec["error"] = "the run under perf did not reach the shutdown banner"
            return rec
        for line in res.stderr.split("\n"):
            parts = line.split(",")
            if len(parts) > 2 and parts[2] == "instructions":
                rec["instructions"] = int(parts[0])
            if len(parts) > 2 and parts[2] == "task-clock":
                rec["task_clock_ms"] = float(parts[0])
        if "instructions" in rec:
            rec["instructions_per_model_year"] = round(
                rec["instructions"] / years / 1e9, 3)
        shutil.rmtree(rundir, ignore_errors=True)
        return rec
    finally:
        cfg.unlink(missing_ok=True)


def compare(ref_dir: Path, got_dir: Path, nc: str) -> dict:
    """Every float variable in the regression case's GOLDSTEIN annual average,
    compared BIT-FOR-BIT. The bar was fixed before any arm was built: these
    changes are correctness-preserving by construction, so a difference is a
    defect to explain rather than a tolerance to widen."""
    import netCDF4  # noqa: PLC0415
    import numpy as np  # noqa: PLC0415

    ref, got = ref_dir / nc, got_dir / nc
    out = {"reference": str(ref), "compared": str(got)}
    if not ref.exists() or not got.exists():
        out["error"] = "one side is missing"
        return out
    differing, worst = {}, 0.0
    with netCDF4.Dataset(ref) as a, netCDF4.Dataset(got) as b:
        names = sorted(set(a.variables) & set(b.variables))
        out["variables"] = len(names)
        out["variables_only_in_reference"] = sorted(set(a.variables) - set(b.variables))
        for name in names:
            x, y = np.asarray(a[name][:]), np.asarray(b[name][:])
            if x.shape != y.shape or x.dtype.kind != "f":
                continue
            if np.array_equal(x, y):
                continue
            d = np.abs(y - x)
            scale = float(np.abs(x).max()) or 1.0
            rel = float(d.max()) / scale
            worst = max(worst, rel)
            differing[name] = {"max_abs_diff": float(d.max()),
                               "max_rel_to_field_range": rel,
                               "cells_differing": int((x != y).sum())}
    out["bit_for_bit"] = not differing
    out["differing"] = differing
    out["worst_relative"] = worst
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--export", action="store_true",
                    help="refresh the exported tree from the working tree")
    ap.add_argument("--rev", default=None,
                    help="export vendor/cgenie at this git revision instead of"
                         " the working tree, so an arm can be measured against"
                         " the tree as it stood before a change")
    ap.add_argument("--arm", action="append", default=None, choices=list(ARMS))
    ap.add_argument("--case", action="append", default=None, choices=list(CASES))
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--threads", action="append", type=int, default=None)
    ap.add_argument("--tag", default="")
    ap.add_argument("--against", default=None,
                    help="the name of a directory under the log directory's nc/"
                         " whose output this run must reproduce")
    ap.add_argument("--shipped", action="store_true",
                    help="compare against genie-knowngood/ rather than an arm")
    ap.add_argument("--cost", action="store_true",
                    help="build and run the cost case under perf stat")
    ap.add_argument("--years", type=int, default=100)
    ap.add_argument("--nyear", type=int, default=100)
    ap.add_argument("--maxisles", type=int, default=20)
    ap.add_argument("--out", type=Path,
                    default=PROJECT_ROOT / "analysis" / "cgenie_omp.json")
    args = ap.parse_args()
    retarget()
    LOGDIR.mkdir(parents=True, exist_ok=True)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    results: dict = {"cases": {}}
    if args.out.exists():
        try:
            results = json.loads(args.out.read_text())
            results.setdefault("cases", {})
        except json.JSONDecodeError:
            pass
    results["provenance"] = cc.provenance()
    results["provenance"]["cgenie_tree"] = str(CGENIE)

    if args.export:
        export_tree(args.rev)
    if args.cost:
        for arm in args.arm or ["base"]:
            first = True
            for threads in (args.threads or [1]):
                tag = f"cost.{arm}{args.tag}.t{threads}"
                rec = cost_run(arm, threads, tag, args.years, args.nyear,
                               args.maxisles, build_it=first)
                first = False
                results["cases"][tag] = rec
                print(json.dumps(rec, indent=2))
        args.out.write_text(json.dumps(results, indent=2) + "\n")
        print(f"wrote {args.out}")
        return 0
    for arm in args.arm or []:
        for case in (args.case or ["eb_go_gs"]):
            spec = CASES[case]
            if args.build:
                rec = build(arm, case)
                results["cases"][f"{arm}.{case}{args.tag}.build"] = rec
                print(json.dumps(rec, indent=2))
                if not rec.get("build_ok"):
                    continue
            if args.run:
                for threads in (args.threads or [1]):
                    tag = f"{arm}.{case}{args.tag}.t{threads}"
                    rec = run_case(arm, case, threads, tag)
                    if rec.get("run_ok"):
                        if args.shipped:
                            rec["vs_shipped"] = compare(
                                CGENIE / spec["reference"], Path(rec["output"]),
                                spec["nc"])
                        if args.against:
                            rec["vs_" + args.against] = compare(
                                LOGDIR / "nc" / args.against,
                                Path(rec["output"]), spec["nc"])
                    results["cases"][tag] = rec
                    print(json.dumps(rec, indent=2))
    args.out.write_text(json.dumps(results, indent=2) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
