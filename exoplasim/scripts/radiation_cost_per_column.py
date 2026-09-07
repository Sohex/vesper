#!/usr/bin/env python
"""The cost of one radiation call, per column, for both schemes on this host.

    qrun -p exoplasim-omp -- \
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

## The binary is built here rather than found

The T21 profiling bed carries the executable WORLD-43RK measured with, and
`radmod.f90` has moved seventeen times since -- WORLD-F9IG alone replaced three
analytic fits with the Stephens (1984) tables, which changes the shortwave
cloud kernel's transcendental count. A cost measurement taken on that binary
would be a measurement of a scheme this project no longer runs, and it would
look exactly like a measurement of the one it does. So a fresh executable is
built from the tree's own source at every invocation, with `--no-publish` so it
stays in its build directory and cannot become the binary a run picks up, and
its sha is reported beside the number.

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


def make_bed(steps: int, from_run: Path, exe: Path) -> Path:
    """A bed at a given step count, built by `make_profile_bed.py`.

    NOT a hand copy of an existing bed, and the difference is the whole
    measurement. Every run directory and every bed on disk was written by an
    older model and carries keys the current `namelist` statements no longer
    declare; `make_profile_bed.py` PRUNES those against `plasim/src` itself and
    FORCES what `config/planet.yaml` declares. Copying a bed instead leaves the
    old keys in place and the model aborts in `readnl` before its first step --
    which it did here, on `tswr3`, a key WORLD-F9IG removed when `swr` stopped
    reading three tuned coefficients. An aborted model is 0.04 s of wall clock
    and a radiation share of zero, and a zero share divides into an
    ordinary-looking slowdown of one.
    """
    dest = WORK / f"bed_{RUNG.lower()}_{steps}"
    if dest.exists():
        shutil.rmtree(dest)
    cmd = [sys.executable, str(ROOT / "exoplasim/scripts/make_profile_bed.py"),
           "--from-run", str(from_run), "--dest", str(dest),
           "--steps", str(steps), "--binary", exe.name,
           "--binary-dir", str(exe.parent), "--allow-unregistered"]
    p = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit("make_profile_bed.py failed:\n" + p.stdout[-3000:]
                         + "\n" + p.stderr[-3000:])
    return dest


def build_current_model(threads: int, verbose: bool = False) -> Path:
    """A fresh executable from THIS tree's source, unpublished.

    `--no-publish` leaves it in its build directory: an arm that publishes
    overwrites the shipped binary and gives the next run unknown provenance,
    which is what that flag exists to prevent.
    """
    cmd = [sys.executable, str(ROOT / "exoplasim/scripts/build_model.py"),
           "--res", RUNG, "--ranks", str(threads),
           "--no-publish", "--print-path"]
    p = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit("build_model.py failed:\n" + p.stdout[-3000:]
                         + "\n" + p.stderr[-3000:])
    path = Path(p.stdout.strip().splitlines()[-1])
    if not path.is_file():
        raise SystemExit(f"build_model.py named {path}, which is not a file")
    return path


def model_exe(bed: Path, name: str) -> Path:
    """The executable inside a bed, BY NAME rather than by glob.

    An unpublished build is `plasim.x` in its build directory, not
    `most_plasim_<res>_l<n>_p<r>.x`, which is the name the registry publishes
    under. Globbing the published pattern found nothing and the run died after
    taking the host lock, which is the expensive way to learn it. The name is
    the one this script handed `make_profile_bed.py`, so the two cannot drift.
    """
    exe = bed / name
    if not exe.is_file():
        raise SystemExit(f"no executable {name} in {bed}; "
                         f"it holds {sorted(q.name for q in bed.glob('*.x'))}")
    return exe


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
    if not data.is_file() or data.stat().st_size == 0:
        raise SystemExit(
            "perf record produced nothing. A profile that returns no samples "
            "returns a radiation share of zero, and a zero share is a ratio of "
            "one and a slowdown of one -- an ordinary-looking answer with "
            "nothing behind it. Check /proc/sys/kernel/perf_event_paranoid.")
    for line in rep.stdout.splitlines():
        m = re.match(r"\s*([\d.]+)%\s+(\S+)\s+\[[.k]\]\s+(\S+)", line)
        if not m:
            continue
        pct, dso, sym = float(m.group(1)), m.group(2), m.group(3)
        shares[f"{dso}:{sym}"] = shares.get(f"{dso}:{sym}", 0.0) + pct
        dso_total[dso] = dso_total.get(dso, 0.0) + pct
    total = sum(dso_total.values())
    if total < 50.0:
        raise SystemExit(
            f"perf report attributed only {total:.1f}% of samples to a named "
            "object. The share this script divides by would be a fraction of a "
            "profile rather than a profile.")
    return dict(by_symbol=shares, by_object=dso_total,
                attributed_percent=round(total, 2))


def radiation_share(shares: dict) -> dict:
    by_sym = shares["by_symbol"]
    by_dso = shares["by_object"]

    def sym(name):
        return sum(v for k, v in by_sym.items() if k.endswith(":" + name))

    # libm.so and libmvec.so are SEPARATE objects here and only the first
    # belongs to radiation. `radiation-scheme-price.md` establishes that every
    # non-integer `**` in the radiation sits inside a `where` block, so none of
    # them reaches libmvec's vector calls: the scalar `pow`, `exp` and `log` are
    # radiation's and the vector ones are the unmasked code elsewhere in the
    # physics. Counting libmvec as radiation would inflate the share with work
    # that is not radiation's, so it is reported beside rather than inside.
    libm = sum(v for k, v in by_dso.items()
               if "libm.so" in k or k.endswith("libm"))
    libmvec = sum(v for k, v in by_dso.items() if "libmvec" in k)
    swr, lwr, radstep = sym("swr_"), sym("lwr_"), sym("radstep_")
    return dict(swr_percent=round(swr, 3), lwr_percent=round(lwr, 3),
                radstep_self_percent=round(radstep, 3),
                libm_percent=round(libm, 3),
                libmvec_percent_not_counted=round(libmvec, 3),
                radiation_percent=round(swr + lwr + radstep + libm, 3),
                radiation_percent_model_code_only=round(swr + lwr + radstep, 3))


def socrates_arm(nprofile: int, nrep: int, nlayer: int, cloud: int,
                 zenith: int = 0) -> dict:
    """One candidate arm.

    `zenith` is an ARM and not a setting. 0 spreads the zenith cosine over a
    full diurnal cycle, so half the columns are dark and whatever the candidate
    does with a dark column is in the number; 1 lights every column at 0.5,
    which is the arm no skipped column can flatter. The broadband scheme it is
    compared against evaluates both branches of a `where` and selects, so it
    does not save on a dark column either, and reporting both is what keeps the
    comparison from resting on that difference.
    """
    if SOC_BENCH is None or not SOC_BENCH.is_file():
        return dict(unavailable="SOCRATES_BENCH is not set to a built driver")
    env = dict(SB_NPROFILE=str(nprofile), SB_NREP=str(nrep),
               SB_NLAYER=str(nlayer), SB_CLOUD=str(cloud),
               SB_ZENITH=str(zenith))
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
                  nprofile=nprofile, nrep=nrep, nlayer=nlayer, cloud=cloud,
                  zenith=zenith)
    return parsed


BLOCK_SWEEP = tuple(25 * 2 ** k for k in range(7))


def socrates_working_set(blocks=BLOCK_SWEEP, nlayer=10, nrep=3) -> dict:
    """The candidate's per-die state, and what blocking the columns costs.

    A thread team's working set on one die targets 32 MB, so what matters is
    not the driver's total but its MARGINAL state per column: the fixed part is
    the spectral tables, one copy a process, and the marginal part is what a
    thread multiplies by its own block. Both are read off a sweep in block size
    rather than counted from declarations, because a count of declarations is
    what missed 4.94 MB of threadprivate copies in the Legendre weight work.

    The same sweep prices blocking, which is the disposition the criterion in
    `radiation-scheme-price.md` names: if the per-column cost is flat in the
    block size then the working set is a caller's choice at no cost in time.
    """
    if SOC_BENCH is None or not SOC_BENCH.is_file():
        return dict(unavailable="SOCRATES_BENCH is not set to a built driver")
    rows = []
    for b in blocks:
        e = dict(os.environ)
        e.update(SB_NPROFILE=str(b), SB_NREP=str(nrep), SB_NLAYER=str(nlayer),
                 SB_ZENITH="1", SB_CLOUD="0")
        p = subprocess.run(["/usr/bin/time", "-f", "%M", str(SOC_BENCH)],
                           cwd=str(SOC_BENCH.parent), env=e,
                           capture_output=True, text=True)
        rss = None
        for line in p.stderr.strip().splitlines()[::-1]:
            try:
                rss = int(line.strip())
                break
            except ValueError:
                continue
        cost = None
        for line in p.stdout.splitlines():
            f = line.split()
            if len(f) == 2 and f[0] == "both_s_per_col_call":
                cost = float(f[1])
        rows.append(dict(columns=b, max_rss_kb=rss,
                         cpu_s_per_column_per_call=cost))
    known = [r for r in rows if r["max_rss_kb"]]
    marginal = fixed = None
    if len(known) >= 2:
        lo, hi = known[0], known[-1]
        marginal = ((hi["max_rss_kb"] - lo["max_rss_kb"])
                    / (hi["columns"] - lo["columns"]))
        fixed = lo["max_rss_kb"] - marginal * lo["columns"]
    # What block a team of eight on one die can carry inside 32 MB.
    budget_kb = 32 * 1024
    per_thread = (int((budget_kb - (fixed or 0)) / (8 * marginal))
                  if marginal else None)
    return dict(rows=rows, marginal_kb_per_column=marginal,
                fixed_kb=fixed, die_budget_kb=budget_kb,
                columns_per_thread_inside_the_budget=per_thread,
                note="eight threads a die, one shared copy of the spectral "
                     "tables; the marginal term is what a thread multiplies")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--steps", type=int, nargs=2, default=[8000, 16000],
                    help="the short and long model arms; the difference is the "
                         "per-step cost and the startup cancels. The default "
                         "pair is WORLD-43RK's own, and it is chosen against the "
                         "startup rather than for convenience: the T21 bed "
                         "starts in about 3.3 s, so a 500-step arm is a quarter "
                         "of its own startup and the difference it would resolve "
                         "sits near this machine's round-to-round scatter. 8,000 "
                         "and 16,000 steps are 16.7 s and 30.1 s there, and the "
                         "13.4 s between them carries 2 per cent scatter rather "
                         "than 6")
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--from-run", type=Path,
                    default=ROOT / "exoplasim/runs/run_0d41aa82c287",
                    help="the run directory the bed is cut from. Its namelists "
                         "are pruned and forced against the current model "
                         "source, so an older run is fine; what it supplies is "
                         "a state and a set of inputs")
    ap.add_argument("--verbose", action="store_true")
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
    exe = build_current_model(args.threads, args.verbose)
    exe_sha = __import__("hashlib").sha256(exe.read_bytes()).hexdigest()
    built = exe
    beds = {n: make_bed(n, args.from_run, built) for n in (short, long)}
    exe = model_exe(beds[short], built.name)
    env = dict(OMP_NUM_THREADS=str(args.threads), OMP_WAIT_POLICY="passive")

    loads = [load()[0]]
    model = {short: [], long: []}
    socr = []
    for r in range(args.rounds):
        order = [short, long] if r % 2 == 0 else [long, short]
        for n in order:
            model[n].append(run_perf_stat([str(exe)], beds[n], env))
            loads.append(load()[0])
        for zen in (0, 1):
            for cloud in (0, 1):
                socr.append(socrates_arm(args.socrates_profiles,
                                         args.socrates_reps, 10, cloud, zen))
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
    bad = [r for arm in model.values() for r in arm if r.get("returncode")]
    if bad:
        raise SystemExit(
            "a model arm exited non-zero, so there is no per-step cost to "
            "report. The first failure's stderr:\n"
            + str(bad[0].get("stderr_tail", ""))[-2000:])
    cpu_per_step_ms = (cpu_long - cpu_short) / d_steps if cpu_long else None
    if not cpu_per_step_ms or cpu_per_step_ms <= 0:
        raise SystemExit(
            f"the long arm cost {cpu_long} ms against the short arm's "
            f"{cpu_short} over {d_steps} more steps, which is a per-step cost "
            "of {cpu_per_step_ms}. A model that does not cost more for running "
            "longer did not run; refusing rather than reporting it.")
    startup_cpu_ms = (cpu_short - cpu_per_step_ms * short) if cpu_per_step_ms else None

    working_set = socrates_working_set()
    shares = perf_record_shares(beds[long], model_exe(beds[long], built.name),
                                args.threads)
    rad = radiation_share(shares)

    columns = args.nlon * args.nlat
    # `radstep` is called from `gridpointd` on every timestep when nrad > 0,
    # which is what config/planet.yaml configures, so one radiation call per
    # step; the report carries the assumption so it can be checked.
    rad_cpu_per_step_ms = (cpu_per_step_ms * rad["radiation_percent"] / 100.0
                           if cpu_per_step_ms else None)
    present_s_per_col_call = (rad_cpu_per_step_ms / 1000.0 / columns
                              if rad_cpu_per_step_ms else None)

    def arm_median(zen, cloud):
        vals = [x["both_s_per_col_call"] for x in socr
                if x.get("zenith") == zen and x.get("cloud") == cloud
                and "both_s_per_col_call" in x]
        return statistics.median(vals) if vals else None

    by_arm = {f"zenith{z}_cloud{c}": arm_median(z, c)
              for z in (0, 1) for c in (0, 1)}
    have = [v for v in by_arm.values() if v]
    cand_lo, cand_hi = (min(have), max(have)) if have else (None, None)

    def slow(v):
        if not (v and present_s_per_col_call):
            return None, None
        r = v / present_s_per_col_call
        return r, 1.0 + rad["radiation_percent"] / 100.0 * (r - 1.0)

    ratio_lo, slowdown_lo = slow(cand_lo)
    ratio_hi, slowdown_hi = slow(cand_hi)

    report = dict(
        measured_on=time.strftime("%Y-%m-%d"),
        host_load=dict(min=round(min(loads), 2), max=round(max(loads), 2),
                       median=round(statistics.median(loads), 2),
                       samples=len(loads)),
        rung=RUNG,
        threads=args.threads,
        rounds=args.rounds,
        executable=dict(path=str(exe), sha256=exe_sha,
                        built_here=True, from_run=str(args.from_run)),
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
            by_arm=by_arm,
            working_set=working_set,
            cpu_s_per_column_per_call=[cand_lo, cand_hi]),
        comparison=dict(
            candidate_over_present=[ratio_lo, ratio_hi],
            whole_model_slowdown_at_this_rung=[slowdown_lo, slowdown_hi],
            note="the slowdown assumes the radiation share measured here and "
                 "that the rest of the model is untouched; it is stated at the "
                 "rung the bed is, and the note carries what it means at T85"),
        by_object=shares["by_object"],
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report["comparison"], indent=2))
    print(json.dumps(report["present_scheme"]["attribution"], indent=2))
    print("candidate s/col/call", by_arm)
    print("present   s/col/call", present_s_per_col_call)


if __name__ == "__main__":
    main()
