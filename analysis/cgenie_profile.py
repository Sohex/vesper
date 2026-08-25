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
The tree is exported out of the repository with `git archive` and built there,
so nothing writes through the worktree's symlinks into the shared checkout.

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

# The vendored tree is read through the repository, but a worktree's ignored
# build products are SYMLINKS into the main checkout, so building in place would
# write there. Build from an exported copy outside the repository instead.
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
        return subprocess.call(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)


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
        cwd=outdir, capture_output=True, text=True, check=False)
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
    out["symbols"] = rows
    shutil.rmtree(outdir / "perf.data", ignore_errors=True)
    data.unlink(missing_ok=True)
    return out


def annotate(rows: list[dict], index: dict[str, dict]) -> list[dict]:
    for row in rows:
        key = row["symbol"].lower().rstrip("_")
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
    ap.add_argument("--logdir", type=Path,
                    default=Path.home() / "cgenie_log" / "profile")
    ap.add_argument("--reuse-builds", action="store_true")
    ap.add_argument("--out", type=Path,
                    default=PROJECT_ROOT / "analysis" / "cgenie_profile.json")
    args = ap.parse_args()

    retarget()
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
        rec["grid"] = f"{case['nlons']}x{case['nlats']}x{case['nlevs']}"
        results["cases"][name] = rec

    if args.biogem or args.only_biogem:
        rec = profile_biogem(args, args.logdir)
        if rec.get("symbols"):
            annotate(rec["symbols"], index)
            by_comp: dict[str, float] = {}
            for row in rec["symbols"]:
                by_comp[row["component"]] = by_comp.get(row["component"], 0.0) + row["percent"]
            rec["by_component"] = dict(sorted(by_comp.items(), key=lambda kv: -kv[1]))
        results["cases"]["eb_go_gs_ac_bg_36x36x8"] = rec

    args.out.write_text(json.dumps(results, indent=2, sort_keys=False) + "\n")
    print(f"wrote {args.out}")
    for name, rec in results["cases"].items():
        print(f"\n== {name} ({rec.get('grid')})  samples={rec.get('samples')}")
        if rec.get("error"):
            print("   ERROR:", rec["error"])
        for row in (rec.get("symbols") or [])[:15]:
            print(f"   {row['percent']:6.2f}%  {row['symbol']:<24} {row.get('component','')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
