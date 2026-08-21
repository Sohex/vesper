#!/usr/bin/env python3
"""Measure rank layouts, including oversubscribed ones, on latency AND throughput.

    python exoplasim/scripts/bench_rank_layout.py --bed <bed> --resolution T127 --rounds 6

Worldbuilding frame: a COMPUTE measurement of the Vesper climate model on this
desktop. Nothing here is about the simulated planet.

An arm is N concurrent jobs, each with its own rank count and its own pinning.
Both numbers are reported and they are DIFFERENT QUESTIONS, which
`exoplasim/notes/rank-layout-benchmark.md` established and this keeps:

  latency     seconds for one job to finish. What "how fast is this run" means.
  throughput  jobs finished per wall second across everything in flight. What
              "how many runs an hour does the machine produce" means.

A layout can lose on one and win on the other, and averaging them into a single
recommendation destroys the only useful thing the pair says.

**Pinning is verified, not assumed.** `rank-layout-benchmark.md`'s central
finding is that every spelling without `:ordered` CONFINES a job to the right
CPUs and silently drops the per-rank binding, so ranks stack on some cores and
leave others idle -- and an arm whose ranks are colliding reports a slow layout
when the launch line was the fault. Before any arm is timed, each rank is asked
for its own `Cpus_allowed_list` and the answer is checked against what the arm
declared. A mismatch aborts rather than measures.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import statistics as st
import shutil
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

# Physical core i owns hardware threads i and i+16 on this host (lscpu -e=CPU,CORE).
# Cores 0-7 are the V-Cache die, 8-15 the frequency die. `rank-layout-benchmark.md`
# measured the dies as indistinguishable for this workload and closed that
# question; these arms are about SMT, which was never measured at all.
SMT_OFFSET = 16
NCORE = 16


def cpus(spec: str) -> list[int]:
    """Parse a comma list, tolerating ranges in what the KERNEL reports back.

    `pe-list` itself refuses ranges, but /proc/self/status writes affinities as
    `0-1,16-17`, so the parser has to read both.
    """
    out: list[int] = []
    for part in spec.split(","):
        if "-" in part:
            a, b = part.split("-")
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def cpu_of(pe: int, bind: str) -> set[int]:
    """The CPUs Open MPI actually gives PE `pe` under this binding unit.

    THE PE INDEX SPACE CHANGES MEANING WITH THE BINDING UNIT, which is the trap
    here and cost an arm that landed on the wrong cores while looking right:

      --bind-to core      PE p is physical core p, and the rank gets BOTH of its
                          threads, {p, p+16}. Range 0-15.
      --bind-to hwthread  PE p is (core p//2, thread p%2), so one cpu,
                          p//2 + 16*(p%2). Range 0-31, enumerated core-major.

    Verified directly: `pe-list=0,1 --bind-to hwthread` puts ranks on cpus 0 and
    16 -- the two threads of core 0 -- while `pe-list=0,16` puts them on cpus 0
    and 8, two different cores. Writing an arm in cpu numbers gets it wrong.

    `pe-list` also refuses RANGES. `0-7` fails with "not enough CPUs in the
    specified PE-LIST"; it has to be spelled out with commas.
    """
    if bind.split(":")[0] == "core":
        return {pe, pe + SMT_OFFSET}
    return {pe // 2 + SMT_OFFSET * (pe % 2)}


def pes_to_list(pes: list[int]) -> str:
    return ",".join(str(p) for p in pes)


# Thread 0 and thread 1 of a given set of cores, in PE-index space.
def t0(cores: range | list[int]) -> list[int]:
    return [2 * c for c in cores]


def t1(cores: range | list[int]) -> list[int]:
    return [2 * c + 1 for c in cores]


# name -> list of jobs, each (ranks, pe-list, bind unit).
ARMS: dict[str, list[tuple[int, str, str]]] = {
    # --- latency ---
    # incumbent: one rank per core, both threads each, as an unpinned launch does
    "L1_1x16":     [(16, pes_to_list(list(range(16))), "core")],
    # SMT: one rank per hardware thread, two per core
    "L2_1x32":     [(32, pes_to_list(list(range(32))), "hwthread")],
    # --- throughput ---
    # control for T2: one job on the eight V-Cache cores
    "T0_1x8":      [(8, pes_to_list(list(range(8))), "core")],
    # incumbent: one job per die
    "T1_2x8_dies": [(8, pes_to_list(list(range(8))), "core"),
                    (8, pes_to_list(list(range(8, 16))), "core")],
    # SMT isolated: two jobs sharing the SAME eight physical cores, one on each
    # core's first thread and one on its second
    "T2_2x8_smt":  [(8, pes_to_list(t0(range(8))), "hwthread"),
                    (8, pes_to_list(t1(range(8))), "hwthread")],
    # SMT, whole machine: two 16-rank jobs, one per thread of every core
    "T3_2x16":     [(16, pes_to_list(t0(range(16))), "hwthread"),
                    (16, pes_to_list(t1(range(16))), "hwthread")],
    # same 32 threads, four jobs
    "T4_4x8":      [(8, pes_to_list(t0(range(8))), "hwthread"),
                    (8, pes_to_list(t1(range(8))), "hwthread"),
                    (8, pes_to_list(t0(range(8, 16))), "hwthread"),
                    (8, pes_to_list(t1(range(8, 16))), "hwthread")],

    # --- the same oversubscription, but MIGRATION-TOLERANT -------------------
    # rank-layout-benchmark.md:135 is explicit that `--bind-to core` giving each
    # rank both SMT siblings is deliberate and should be left alone: "a rank that
    # can move to its sibling can step around a daemon landing on its thread,
    # where a rank pinned to one hardware thread simply waits -- and in a
    # bulk-synchronous code that wait propagates to every other rank through the
    # next collective." It measured 19 to 383 migrations per rank over 80 s
    # during a production run, so this is routine rather than hypothetical.
    #
    # Every arm above binds to hwthread, which is exactly the pinned-to-one-
    # thread case that warning describes. That handicaps all of them against the
    # core-bound controls for a reason unrelated to SMT. These variants put the
    # same two ranks on a core while leaving each free to use either of its
    # threads, which is the honest way to ask whether SMT helps.
    "L2b_1x32_core": [(32, "*", "core:overload-allowed")],
    "T2b_2x8_core":  [(8, pes_to_list(list(range(8))), "core"),
                      (8, pes_to_list(list(range(8))), "core")],
    "T3b_2x16_core": [(16, pes_to_list(list(range(16))), "core"),
                      (16, pes_to_list(list(range(16))), "core")],
    "T4b_4x8_core":  [(8, pes_to_list(list(range(8))), "core"),
                      (8, pes_to_list(list(range(8))), "core"),
                      (8, pes_to_list(list(range(8, 16))), "core"),
                      (8, pes_to_list(list(range(8, 16))), "core")],

    # --- does the L3 difference matter UNDER CONTENTION ----------------------
    # rank-layout-benchmark.md measured the dies as indistinguishable and closed
    # the question. It measured ONE job per die. Two jobs sharing one die's L3 is
    # a different cache regime: CCD0 carries 96 MB and CCD1 32 MB, and doubling
    # the jobs on a die doubles the footprint competing for it. That is a reason
    # to extend the measurement rather than to reopen a settled one, and it is
    # not "the dies look different on paper".
    #
    # A 2x2, so the contention penalty is attributable per die rather than
    # confounded with the die itself:
    #     T0        1 job  on CCD0        T0b_1x8_ccd1   1 job  on CCD1
    #     T2b       2 jobs on CCD0        T2c_2x8_ccd1   2 jobs on CCD1
    # If (T2c/T0b) is worse than (T2b/T0), the smaller L3 is paying for the
    # contention and the die matters once a die is shared.
    #
    # Expect this to be resolution-dependent and to show at T127 if anywhere: the
    # Legendre weights are 0.23 MB per rank at T42, where 16 ranks fit either
    # die, and 6 MB per rank at T127, where eight ranks want 48 MB against
    # CCD1's 32.
    "T0b_1x8_ccd1":  [(8, pes_to_list(list(range(8, 16))), "core")],
    "T2c_2x8_ccd1":  [(8, pes_to_list(list(range(8, 16))), "core"),
                      (8, pes_to_list(list(range(8, 16))), "core")],
}


def mpi_opts(pe: str, bind: str, ranks: int = 0) -> list[str]:
    # `:ordered` is load-bearing. Without it the job is confined to the right
    # CPUs and the per-rank binding is silently dropped.
    # `*` means the whole machine rather than a pe-list. Needed for the single
    # 32-rank core-bound arm: two ranks per core cannot be spelled as a pe-list,
    # since a PE appears once.
    if pe == "*":
        opts = ["--map-by", "core", "--bind-to", bind]
    else:
        opts = ["--map-by", f"pe-list={pe}:ordered", "--bind-to", bind]
    # Open MPI counts a slot per physical CORE, so a single job asking for more
    # ranks than there are cores does not launch at all -- "not enough slots
    # available" -- however the binding is spelled. This raises the slot count
    # and nothing else: verified that 32 ranks still land on 32 distinct cpus,
    # one each. Only the one-job-32-rank arm needs it; the concurrent arms are 8
    # or 16 ranks per mpiexec and each fits.
    if ranks > NCORE:
        opts.append("--oversubscribe")
    return opts


def job_env(yield_when_idle: bool) -> dict:
    """Environment for one job.

    `mpi_yield_when_idle` decides whether a rank waiting in a collective spins
    or calls sched_yield. It defaults to FALSE, and Open MPI turns it on by
    itself only "when oversubscribing nodes" -- that is, when one mpiexec knows
    it has asked for more ranks than slots.

    That auto-detection does not fire for the arms that need it most. T2, T3 and
    T4 are several INDEPENDENT mpiexec invocations of 8 or 16 ranks each: every
    one of them believes it has the machine to itself, so yield stays off and two
    aggressive spinners share every physical core. Measuring SMT that way would
    be measuring it at close to its worst case.

    So it is a declared dimension of the matrix rather than a default, and the
    profile says why it could matter: MPI is 11 to 27 percent of samples at T127
    and most of it is opal_progress, which is spinning rather than moving bytes.
    """
    env = dict(os.environ)
    env["OMPI_MCA_mpi_yield_when_idle"] = "1" if yield_when_idle else "0"
    return env


def check_pinning(bed: Path, ranks: int, pe: str, bind: str) -> list[str]:
    """Ask each rank for its own affinity and compare against cpu_of."""
    r = subprocess.run(
        ["mpiexec", "-np", str(ranks), *mpi_opts(pe, bind, ranks),
         "sh", "-c", "grep Cpus_allowed_list /proc/self/status"],
        cwd=bed, capture_output=True, text=True)
    got = re.findall(r"Cpus_allowed_list:\s*(\S+)", r.stdout)
    if len(got) != ranks:
        tail = (r.stderr.strip().splitlines() or ["(no stderr)"])[-1]
        return [f"expected {ranks} affinity lines, got {len(got)}: {tail}"]

    if pe == "*":
        want = [cpu_of(c, bind) for c in range(NCORE)]
    else:
        want = [cpu_of(p, bind) for p in cpus(pe)]
    # Two ranks on one core is the POINT of the overload arms, so sharing is
    # only a fault when the arm did not ask for it.
    shared_ok = ranks > len(want)
    want_all = set().union(*want) if want else set()
    seen: set[int] = set()
    bad = []
    for line in got:
        c = set(cpus(line))
        if c not in want:
            bad.append(f"rank on cpus {sorted(c)}, which is not any declared PE "
                       f"of pe-list={pe} under --bind-to {bind}")
        if not c <= want_all:
            bad.append(f"rank on cpus {sorted(c)}, outside the arm's cpus "
                       f"{sorted(want_all)}")
        if c & seen and not shared_ok:
            bad.append(f"two ranks share cpus {sorted(c & seen)}; binding did not take")
        seen |= c
    return bad


def run_job(bed: Path, exe: Path, ranks: int, pe: str, bind: str, tag: str,
            env: dict) -> float:
    local = bed / f"probe_{tag}.x"
    shutil.copy2(exe, local)
    (bed / "Abort_Message").unlink(missing_ok=True)
    t0 = time.perf_counter()
    r = subprocess.run(["mpiexec", "-np", str(ranks), *mpi_opts(pe, bind, ranks),
                        f"./{local.name}"],
                       cwd=bed, capture_output=True, text=True, env=env)
    dt = time.perf_counter() - t0
    local.unlink(missing_ok=True)
    # A FAILING job is fast, and an unchecked timing reports that as a win. The
    # 32-rank arm died instantly in get_restart_array -- a 16-rank restart is not
    # readable at 32 ranks -- and came back as 2.18 s against 12 s for the
    # control, which reads as a 5x speedup rather than as a crash. Never time a
    # process without checking whether it ran.
    if r.returncode != 0 or (bed / "Abort_Message").exists():
        tail = "\n".join(l for l in (r.stderr or r.stdout).splitlines()
                         if "Fortran runtime warning" not in l
                         and not l.startswith("At line"))[-500:]
        raise RuntimeError(
            f"job failed after {dt:.2f}s: {ranks} ranks, pe-list={pe}, "
            f"bind={bind}\n{tail}")
    return dt


def run_arm(beds: list[Path], exes: list[Path], jobs, env: dict) -> tuple[list[float], float]:
    """Launch every job at once; return per-job seconds and the wall for all of them."""
    t0 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=len(jobs)) as pool:
        futures = [pool.submit(run_job, beds[i], exes[i], r, pe, bind, str(i), env)
                   for i, (r, pe, bind) in enumerate(jobs)]
        per_job = [f.result() for f in futures]
    return per_job, time.perf_counter() - t0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", type=Path, required=True,
                    help="template bed; one private copy is made per concurrent job")
    ap.add_argument("--resolution", required=True)
    ap.add_argument("--rounds", type=int, default=6)
    ap.add_argument("--arms", nargs="*", default=list(ARMS))
    ap.add_argument("--yield-when-idle", action="store_true",
                    help="set OMPI_MCA_mpi_yield_when_idle=1 for every job. Off "
                         "by default, which is Open MPI's default and is what "
                         "the concurrent arms get on their own, since each "
                         "mpiexec believes it has the machine to itself.")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    template = args.bed.resolve()
    low = args.resolution.lower()
    ref = Path("exoplasim/bench/ref").resolve()
    run = Path("vendor/exoplasim/exoplasim/plasim/run").resolve()

    def exe_for(ranks: int) -> Path:
        name = f"most_plasim_{low}_l10_p{ranks}.x"
        for d in (run, ref):
            if (d / name).is_file():
                return d / name
        raise SystemExit(f"no binary {name} in {run} or {ref}")

    # One private bed per concurrent job. Concurrent jobs writing the same
    # directory would race on plasim_status and on the probe executable.
    nmax = max(len(ARMS[a]) for a in args.arms)
    beds = []
    for i in range(nmax):
        b = template.parent / f"{template.name}_job{i}"
        if not b.is_dir():
            shutil.copytree(template, b)
        for stale in ("plasim_status", "Abort_Message"):
            (b / stale).unlink(missing_ok=True)
        beds.append(b)
    print(f"{nmax} private bed copies under {template.parent}")

    print("\nverifying pinning before any arm is timed")
    for name in args.arms:
        for ranks, pe, bind in ARMS[name]:
            bad = check_pinning(beds[0], ranks, pe, bind)
            if bad:
                print(f"  {name}: FAILED")
                for b in bad:
                    print(f"    {b}")
                raise SystemExit("pinning did not take; refusing to measure a launch line")
        print(f"  {name}: ok")

    env = job_env(args.yield_when_idle)
    print(f"\nmpi_yield_when_idle = {env['OMPI_MCA_mpi_yield_when_idle']}")

    print(f"warm-up ({args.arms[0]})")
    run_arm(beds, [exe_for(r) for r, _, _ in ARMS[args.arms[0]]], ARMS[args.arms[0]], env)

    results: dict[str, dict] = {a: {"per_job": [], "wall": []} for a in args.arms}
    for rnd in range(1, args.rounds + 1):
        # Rotate arm order each round so none is always first after an idle gap.
        order = args.arms[rnd % len(args.arms):] + args.arms[:rnd % len(args.arms)]
        for name in order:
            jobs = ARMS[name]
            per_job, wall = run_arm(beds, [exe_for(r) for r, _, _ in jobs], jobs, env)
            results[name]["per_job"].append(per_job)
            results[name]["wall"].append(wall)
        print(f"  round {rnd} done", flush=True)

    print(f"\n{'arm':16} {'jobs':>4} {'ranks':>6} {'latency s':>10} "
          f"{'throughput':>11} {'spread':>7}")
    table = {}
    for name in args.arms:
        jobs = ARMS[name]
        lat = st.median([st.median(p) for p in results[name]["per_job"]])
        walls = results[name]["wall"]
        thr = len(jobs) / st.median(walls)          # jobs per wall second
        spread = 100 * (max(walls) - min(walls)) / st.median(walls)
        table[name] = {"jobs": len(jobs), "ranks_total": sum(r for r, _, _ in jobs),
                       "latency_s": lat, "throughput_jobs_per_s": thr,
                       "wall_spread_pct": spread,
                       "per_job": results[name]["per_job"], "wall": walls}
        print(f"{name:16} {len(jobs):4d} {sum(r for r,_,_ in jobs):6d} "
              f"{lat:10.2f} {thr*3600:11.1f} {spread:6.1f}%")
    print("  throughput is beds per wall HOUR across all jobs in flight")

    # The declared verdicts. Thresholds are rank-layout-benchmark.md's and are
    # not renegotiated here: under 5% is no difference, and a throughput layout
    # is adopted only by beating the incumbent by more than 10%.
    print()
    if "L1_1x16" in table and "L2_1x32" in table:
        base, smt = table["L1_1x16"]["latency_s"], table["L2_1x32"]["latency_s"]
        d = 100 * (base - smt) / base
        verdict = "ADOPT" if d > 5 else ("no difference" if abs(d) <= 5 else "REFUSE")
        print(f"LATENCY   1x32 vs 1x16: {d:+.1f}%  -> {verdict}")
    if "T1_2x8_dies" in table:
        inc = table["T1_2x8_dies"]["throughput_jobs_per_s"]
        for name in args.arms:
            if name.startswith("T") and name != "T1_2x8_dies":
                d = 100 * (table[name]["throughput_jobs_per_s"] - inc) / inc
                verdict = "ADOPT" if d > 10 else "keep incumbent"
                print(f"THROUGHPUT {name} vs 2x8-dies: {d:+.1f}%  -> {verdict}")

    for name in args.arms:
        if table[name]["wall_spread_pct"] > 5.0:
            print(f"  ** {name} self-scatter {table[name]['wall_spread_pct']:.1f}% "
                  f"exceeds the 5% floor: that arm is NOT a result.")

    if args.out:
        args.out.write_text(json.dumps(
            {"resolution": args.resolution, "rounds": args.rounds,
             "yield_when_idle": args.yield_when_idle,
             "arms": {k: ARMS[k] for k in args.arms}, "results": table}, indent=2) + "\n")
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
