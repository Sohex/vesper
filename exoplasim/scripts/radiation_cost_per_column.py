#!/usr/bin/env python
"""The cost of one radiation call, per column, for both schemes on this host.

    scripts/lock_and_run -m "CLIM-61 radiation cost" \
        python exoplasim/scripts/radiation_cost_per_column.py --steps 500 2500

Worldbuilding frame: a COMPUTE measurement of the Vesper climate model and of a
candidate radiation code on this desktop. Nothing here is about the simulated
planet.

`exoplasim/notes/radiation-scheme-price.md` establishes the SHAPE of both
schemes' cost in retired instructions and then refuses to convert either into a
cost per timestep, because "instructions are not cycles" and WORLD-43RK measured
an instruction count of exactly that kind over-stating a wall clock by 5.6. This
script produces the number that note says is still missing, for both schemes, in
the same unit.

## The unit is CPU-seconds per column per radiation call, and that choice is the
## whole design

The model is threaded and the candidate's driver is not, so wall clock is not a
common currency. CPU time is: both schemes are embarrassingly parallel over
columns, so "how much processor time does one column of one radiation call take"
is a property of the scheme rather than of the thread count. `task-clock` is read
straight out of `perf stat`, and the model is run with `OMP_WAIT_POLICY=passive`
so that a thread waiting at a barrier sleeps instead of retiring spin.

## The bed is sized against its own startup, and the sizing is the measurement

`exoplasim/notes/where-the-time-goes.md` records a confident T42 result on a
2.94 s bed against a 1.6 s startup that REVERSED SIGN when the bed was
lengthened. So the model is run at two step counts and the per-step cost is the
DIFFERENCE over the extra steps, which cancels the startup exactly rather than
assuming it small. The two counts and the difference they resolve are reported,
and if the shorter arm is not clearly above the startup the report says so.

The candidate's driver does the same thing by construction: its spectral file
read is timed separately and reported, and one warm-up call is discarded before
the timed loop, so the k-table touch and the allocator are not in the number.

## The attribution, and what it deliberately includes

`perf record` gives the model's cost by symbol. Radiation is `swr_` and `lwr_`
plus the libm they reach: `where-the-time-goes.md` establishes that the call
graph puts all of the `pow` and `log` under `radstep`, and the masked scalar
`pow` is the present scheme's kernel. So the radiation share is `swr_` + `lwr_`
+ `radstep_` self + everything in libm, and the report carries each term
separately so a reader can take a narrower definition if they want one.

## Interleaved, never blocked

`notes/audits/aocl-and-model-build-flags.md` measured a 6 per cent swing between
blocked sessions: whichever arm ran first after an idle stretch took the boost
clock. The arms here alternate, and the host load is sampled throughout and
reported beside every number, because a timing without the machine state it was
taken under cannot be compared with a later one.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

import _paths  # noqa: F401  anchors every path on this file and adds lib/

ROOT = Path(_paths.__file__).resolve().parents[2]
BENCH = ROOT / "exoplasim" / "bench"
WORK = Path(os.environ.get("RADCOST_WORK", BENCH / "climprice"))

# The rung comes from config/planet.yaml through lib/rungs.py, never from a
# literal in a path: SPAT-2's defect was a graph naming one rung while the file
# on disk carried another. The grid dimensions the report divides by come from
# the same call, so the column count and the bed cannot disagree.
import rungs  # noqa: E402
import yaml  # noqa: E402

CONFIG = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text())
RUNG, NLAT, NLON = rungs.model_grid(CONFIG)
SOURCE_BED = BENCH / "43rk" / f"bed_{RUNG.lower()}"

# The candidate. Built outside the tree, because SOCRATES is read-only
# reference under `references/` and a build there would write through a
# worktree's directory link into the main checkout.
SOC_BENCH = Path(os.environ.get("SOCRATES_BENCH", "")) or None


def load() -> tuple[float, float, float]:
    return os.getloadavg()


def run_perf_stat(cmd, cwd, env=None):
    """`perf stat` around one command, returning task-clock and instructions."""
    full = ["perf", "stat", "-x,", "-e", "task-clock,instructions:u", "--"] + cmd
    e = dict(os.environ)
    e.update(env or {})
    t0 = time.perf_counter()
    p = subprocess.run(full, cwd=str(cwd), env=e, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    out = dict(wall_s=wall, returncode=p.returncode)
    for line in p.stderr.splitlines():
        f = line.split(",")
        if len(f) < 3:
            continue
        try:
            v = float(f[0])
        except ValueError:
            continue
        if f[2].startswith("task-clock"):
            out["task_clock_ms"] = v
        elif f[2].startswith("instructions"):
            out["instructions"] = v
    if p.returncode != 0:
        out["stderr_tail"] = p.stderr[-2000:]
    return out


def make_bed(steps: int) -> Path:
    """A private copy of the T21 profiling bed at a given step count.

    Copied rather than edited in place: `exoplasim/bench/43rk` belongs to
    WORLD-43RK's measurement and a gate there deletes its own subdirectory.
    """
    dest = WORK / f"bed_{RUNG.lower()}_{steps}"
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(SOURCE_BED, dest, symlinks=True)
    nl = dest / "plasim_namelist"
    text = nl.read_text()
    text = re.sub(r"N_RUN_STEPS\s*=\s*\d+", f"N_RUN_STEPS = {steps}", text)
    nl.write_text(text)
    return dest


def model_exe(bed: Path) -> Path:
    exes = sorted(bed.glob("most_plasim_*.x"))
    if not exes:
        raise SystemExit(f"no model executable in {bed}")
    return exes[0]


def perf_record_shares(bed: Path, exe: Path, threads: int) -> dict:
    """The model's cost by symbol, as a share of its own samples."""
    data = bed / "perf.data"
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS=str(threads), OMP_WAIT_POLICY="passive")
    subprocess.run(["perf", "record", "-q", "-F", "999", "-o", str(data),
                    "--", str(exe)], cwd=str(bed), env=env,
                   capture_output=True, text=True)
    rep = subprocess.run(["perf", "report", "-q", "-i", str(data),
                          "--no-children", "--percentage", "absolute",
                          "--sort", "dso,symbol", "--stdio"],
                         cwd=str(bed), capture_output=True, text=True)
    shares, dso_total = {}, {}
    for line in rep.stdout.splitlines():
        m = re.match(r"\s*([\d.]+)%\s+(\S+)\s+\[[.k]\]\s+(\S+)", line)
        if not m:
            continue
        pct, dso, sym = float(m.group(1)), m.group(2), m.group(3)
        shares[f"{dso}:{sym}"] = shares.get(f"{dso}:{sym}", 0.0) + pct
        dso_total[dso] = dso_total.get(dso, 0.0) + pct
    return dict(by_symbol=shares, by_object=dso_total)


def radiation_share(shares: dict) -> dict:
    by_sym = shares["by_symbol"]
    by_dso = shares["by_object"]

    def sym(name):
        return sum(v for k, v in by_sym.items() if k.endswith(":" + name))

    libm = sum(v for k, v in by_dso.items() if "libm" in k)
    swr, lwr, radstep = sym("swr_"), sym("lwr_"), sym("radstep_")
    return dict(swr_percent=round(swr, 3), lwr_percent=round(lwr, 3),
                radstep_self_percent=round(radstep, 3),
                libm_percent=round(libm, 3),
                radiation_percent=round(swr + lwr + radstep + libm, 3),
                radiation_percent_model_code_only=round(swr + lwr + radstep, 3))


def socrates_arm(nprofile: int, nrep: int, nlayer: int, cloud: int) -> dict:
    if SOC_BENCH is None or not SOC_BENCH.is_file():
        return dict(unavailable="SOCRATES_BENCH is not set to a built driver")
    env = dict(SB_NPROFILE=str(nprofile), SB_NREP=str(nrep),
               SB_NLAYER=str(nlayer), SB_CLOUD=str(cloud))
    cwd = SOC_BENCH.parent
    e = dict(os.environ)
    e.update(env)
    stat = run_perf_stat([str(SOC_BENCH)], cwd, env)
    p = subprocess.run([str(SOC_BENCH)], cwd=str(cwd), env=e,
                       capture_output=True, text=True)
    parsed = {}
    for line in p.stdout.splitlines():
        f = line.split()
        if len(f) == 2:
            try:
                parsed[f[0]] = float(f[1])
            except ValueError:
                pass
    parsed.update(task_clock_ms=stat.get("task_clock_ms"),
                  instructions=stat.get("instructions"),
                  nprofile=nprofile, nrep=nrep, nlayer=nlayer, cloud=cloud)
    return parsed


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--steps", type=int, nargs=2, default=[500, 2500],
                    help="the short and long model arms; the difference is the "
                         "per-step cost and the startup cancels")
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--nlon", type=int, default=NLON)
    ap.add_argument("--nlat", type=int, default=NLAT)
    ap.add_argument("--socrates-profiles", type=int, default=2048)
    ap.add_argument("--socrates-reps", type=int, default=200)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "exoplasim/analysis/radiation_cost_per_column.json")
    args = ap.parse_args()

    WORK.mkdir(parents=True, exist_ok=True)
    short, long = sorted(args.steps)
    beds = {n: make_bed(n) for n in (short, long)}
    exe = model_exe(beds[short])
    env = dict(OMP_NUM_THREADS=str(args.threads), OMP_WAIT_POLICY="passive")

    loads = [load()[0]]
    model = {short: [], long: []}
    socr = []
    for r in range(args.rounds):
        order = [short, long] if r % 2 == 0 else [long, short]
        for n in order:
            model[n].append(run_perf_stat([str(exe)], beds[n], env))
            loads.append(load()[0])
        socr.append(socrates_arm(args.socrates_profiles, args.socrates_reps,
                                 10, 0))
        loads.append(load()[0])

    def med(rows, key):
        vals = [r[key] for r in rows if r.get(key) is not None]
        return statistics.median(vals) if vals else None

    cpu_short = med(model[short], "task_clock_ms")
    cpu_long = med(model[long], "task_clock_ms")
    ins_short = med(model[short], "instructions")
    ins_long = med(model[long], "instructions")
    wall_short = med(model[short], "wall_s")
    wall_long = med(model[long], "wall_s")

    d_steps = long - short
    cpu_per_step_ms = (cpu_long - cpu_short) / d_steps if cpu_long else None
    startup_cpu_ms = (cpu_short - cpu_per_step_ms * short) if cpu_per_step_ms else None

    shares = perf_record_shares(beds[long], exe, args.threads)
    rad = radiation_share(shares)

    columns = args.nlon * args.nlat
    # `radstep` is called from `gridpointd` on every timestep when nrad > 0,
    # which is what config/planet.yaml configures, so one radiation call per
    # step; the report carries the assumption so it can be checked.
    rad_cpu_per_step_ms = (cpu_per_step_ms * rad["radiation_percent"] / 100.0
                           if cpu_per_step_ms else None)
    present_s_per_col_call = (rad_cpu_per_step_ms / 1000.0 / columns
                              if rad_cpu_per_step_ms else None)

    cand = [s for s in socr if "both_s_per_col_call" in s]
    cand_med = statistics.median(s["both_s_per_col_call"] for s in cand) if cand else None
    ratio = (cand_med / present_s_per_col_call
             if cand_med and present_s_per_col_call else None)
    slowdown = (1.0 + rad["radiation_percent"] / 100.0 * (ratio - 1.0)
                if ratio else None)

    report = dict(
        measured_on=time.strftime("%Y-%m-%d"),
        host_load=dict(min=round(min(loads), 2), max=round(max(loads), 2),
                       median=round(statistics.median(loads), 2),
                       samples=len(loads)),
        rung=RUNG,
        threads=args.threads,
        rounds=args.rounds,
        bed=dict(source=str(SOURCE_BED.relative_to(ROOT)),
                 steps_short=short, steps_long=long,
                 wall_s_short=wall_short, wall_s_long=wall_long,
                 cpu_ms_short=cpu_short, cpu_ms_long=cpu_long,
                 startup_cpu_ms=startup_cpu_ms,
                 short_arm_over_startup=(cpu_short / startup_cpu_ms
                                         if startup_cpu_ms else None)),
        present_scheme=dict(
            cpu_ms_per_step=cpu_per_step_ms,
            instructions_per_step=((ins_long - ins_short) / d_steps
                                   if ins_long else None),
            attribution=rad,
            columns=columns,
            radiation_calls_per_step=1,
            cpu_s_per_column_per_call=present_s_per_col_call),
        candidate=dict(
            spectral_files="ga7",
            arms=socr,
            cpu_s_per_column_per_call=cand_med),
        comparison=dict(
            candidate_over_present=ratio,
            whole_model_slowdown_at_this_rung=slowdown,
            note="the slowdown assumes the radiation share measured here and "
                 "that the rest of the model is untouched; it is stated at the "
                 "rung the bed is, and the note carries what it means at T85"),
        by_object=shares["by_object"],
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["comparison"], indent=2))
    print(json.dumps(report["present_scheme"]["attribution"], indent=2))
    print("candidate s/col/call", cand_med)
    print("present   s/col/call", present_s_per_col_call)


if __name__ == "__main__":
    main()
