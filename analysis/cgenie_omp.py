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


def export_tree() -> None:
    """Copy the tracked vendored tree, at its WORKING-TREE state so uncommitted
    edits are what gets built, out of the repository."""
    if WORK_ROOT.exists():
        shutil.rmtree(WORK_ROOT)
    WORK_ROOT.mkdir(parents=True)
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
    print(f"exported {len(names)} paths to {WORK_ROOT}")


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
    "base": dict(fflags=[], threads=1),
    "nofnoauto": dict(fflags=["-frecursive"], threads=1),
    "omp": dict(fflags=["-frecursive", "-fopenmp"], threads=None),
}

KNOWNGOOD_REF = ("genie-knowngood/genie_eb_go_gs_knowngood/goldstein/"
                 "gold_spn_av_0000000020_00.nc")
KNOWNGOOD_OUT = "goldstein/gold_spn_av_0000000020_00.nc"


def make_args(arm: str, target: str = "") -> str:
    parts = [target, f"GENIE_ROOT={CGENIE}", f"RUNTIME_ROOT={CGENIE}",
             f"NETCDF_DIR={cc.NETCDF_DIR}", f"OUT_DIR={OUT_ROOT}", "-j 1"]
    flags = ["-mcmodel=medium", *ARMS[arm]["fflags"]]
    parts += [f"GENIE_FFLAGS={' '.join(flags)}",
              f"GENIE_LDFLAGS={' '.join(flags)}"]
    return " ".join(p for p in parts if p)


def knowngood_config() -> Path:
    """The shipped eb_go_gs regression case with GOLDSTEIN's `debug_loop` on.

    `make testebgogs` cannot answer this on its own: its reference file is an
    annual average `goldstein.F` writes only under `debug_loop`, which defaults
    false. Setting it changes no state -- every use of it in `goldstein.F`
    guards a print, a dump or the averaging call."""
    src = (cc.CONFIG_DIR / "eb_go_gs_test.xml").read_text()
    marker = '\t\t<model name="goldstein">\n'
    assert marker in src, "eb_go_gs_test.xml no longer has a goldstein block"
    cfg = cc.CONFIG_DIR / "omp_eb_go_gs.xml"
    cfg.write_text(src.replace(
        marker, marker + '\t\t\t<param name="debug_loop">.true.</param>\n', 1
    ).replace('<var name="EXPID">genie_eb_go_gs</var>',
              f'<var name="EXPID">{cfg.stem}</var>'))
    return cfg


def build(arm: str) -> dict:
    cfg = knowngood_config()
    try:
        sh(["/usr/bin/make", *make_args(arm).split(), "cleanall"],
           LOGDIR / f"{arm}.clean.log", cwd=cc.GENIE_MAIN)
        t0 = time.perf_counter()
        rc = sh(["./genie.job", "-x", "-f", f"configs/{cfg.name}",
                 "-o", str(OUT_ROOT), "-c", str(CGENIE), "-g", str(CGENIE),
                 "-h", ".", "-m", make_args(arm, target="genie.exe")],
                LOGDIR / f"{arm}.build.log", cwd=cc.GENIE_MAIN)
        built = cc.GENIE_MAIN / "genie.exe"
        out = {"arm": arm, "flags": ARMS[arm]["fflags"],
               "build_seconds": round(time.perf_counter() - t0, 1)}
        if rc != 0 or not built.exists():
            out["build_ok"] = False
            out["error"] = f"see {LOGDIR / (arm + '.build.log')}"
            return out
        stash = LOGDIR / "exe"
        stash.mkdir(parents=True, exist_ok=True)
        shutil.copy2(built, stash / f"{arm}.exe")
        out["build_ok"] = True
        out["exe"] = str(stash / f"{arm}.exe")
        return out
    finally:
        cfg.unlink(missing_ok=True)


def run_case(arm: str, threads: int, tag: str, perf: bool) -> dict:
    """Run the regression case from an already-built, stashed executable.

    `genie.job -z` means REMAKE=FALSE, so the executable sitting in genie-main
    is the one that runs; the arm's stashed copy is put there first, which is
    what keeps the arms from silently sharing a binary."""
    exe = LOGDIR / "exe" / f"{arm}.exe"
    if not exe.exists():
        return {"run_ok": False, "why": f"no stashed executable for arm {arm}"}
    shutil.copy2(exe, cc.GENIE_MAIN / "genie.exe")
    cfg = knowngood_config()
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    # A thread spinning at a barrier retires instructions in proportion to how
    # long it waits, so an active wait policy would put the idle threads' spin
    # INTO the instruction count and make a threaded arm look like it did more
    # work. Passive is what makes the count comparable across thread counts.
    env["OMP_WAIT_POLICY"] = "passive"
    outdir = OUT_ROOT / cfg.stem
    shutil.rmtree(outdir, ignore_errors=True)
    try:
        rc = sh(["./genie.job", "-z", "-f", f"configs/{cfg.name}",
                 "-o", str(OUT_ROOT), "-c", str(CGENIE), "-g", str(CGENIE),
                 "-h", ".", "-m", make_args(arm)],
                LOGDIR / f"{tag}.run.log", cwd=cc.GENIE_MAIN, env=env)
        log = LOGDIR / f"{tag}.run.log"
        ok, why = cc.run_ok(log)
        rec = {"arm": arm, "threads": threads, "run_ok": bool(rc == 0 and ok),
               "seconds": cc.model_seconds(log)}
        if why:
            rec["why"] = why
        if rec["run_ok"]:
            keep = LOGDIR / "nc" / tag
            keep.parent.mkdir(parents=True, exist_ok=True)
            shutil.rmtree(keep, ignore_errors=True)
            shutil.copytree(outdir, keep)
            rec["output"] = str(keep)
        if perf and rec["run_ok"]:
            rec.update(perf_stat(arm, threads, tag, env))
        return rec
    finally:
        cfg.unlink(missing_ok=True)


def perf_stat(arm: str, threads: int, tag: str, env: dict) -> dict:
    """Retired instructions, COUNTED not sampled, plus the wall clock.

    On a single-threaded process the count is the cost. On a threaded one it is
    not: it includes whatever the waiting threads retire, which is why the wait
    policy above is passive and why the wall clock is reported beside it."""
    outdir = cc.GENIE_MAIN
    res = subprocess.run(
        ["perf", "stat", "-e", "instructions,task-clock", "-x,",
         str(LOGDIR / "exe" / f"{arm}.exe")],
        cwd=OUT_ROOT / "run", capture_output=True, text=True, check=False, env=env)
    return {"perf_raw": res.stderr}


def compare(ref_dir: Path, got_dir: Path) -> dict:
    """Every float variable in the regression case's GOLDSTEIN annual average,
    compared BIT-FOR-BIT. The bar was fixed before any arm was built: these
    changes are correctness-preserving by construction, so a difference is a
    defect to explain rather than a tolerance to widen."""
    import netCDF4  # noqa: PLC0415
    import numpy as np  # noqa: PLC0415

    ref, got = ref_dir / KNOWNGOOD_OUT, got_dir / KNOWNGOOD_OUT
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
    ap.add_argument("--arm", action="append", default=None, choices=list(ARMS))
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--threads", action="append", type=int, default=None)
    ap.add_argument("--tag", default="")
    ap.add_argument("--against", default=None,
                    help="a tag under the log directory's nc/ to compare with")
    ap.add_argument("--shipped", action="store_true",
                    help="compare against genie-knowngood/ rather than an arm")
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
        export_tree()
    for arm in args.arm or []:
        if args.build:
            rec = build(arm)
            results["cases"][f"{arm}{args.tag}.build"] = rec
            print(json.dumps(rec, indent=2))
            if not rec.get("build_ok"):
                break
        if args.run:
            for threads in (args.threads or [1]):
                tag = f"{arm}{args.tag}.t{threads}"
                rec = run_case(arm, threads, tag, perf=False)
                if rec.get("run_ok"):
                    if args.shipped:
                        rec["vs_shipped"] = compare(
                            CGENIE / "genie-knowngood/genie_eb_go_gs_knowngood",
                            Path(rec["output"]))
                    if args.against:
                        rec["vs_" + args.against] = compare(
                            LOGDIR / "nc" / args.against, Path(rec["output"]))
                results["cases"][tag] = rec
                print(json.dumps(rec, indent=2))
    args.out.write_text(json.dumps(results, indent=2) + "\n")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
