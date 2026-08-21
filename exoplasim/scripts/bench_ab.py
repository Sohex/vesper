#!/usr/bin/env python3
"""Paired, interleaved A/B of two executables on one bed.

    python exoplasim/scripts/bench_ab.py --bed exoplasim/bench/bed_t85 \
        --a exoplasim/bench/ref/unpatched_t85.x \
        --b exoplasim/bench/ref/patched_t85.x --rounds 6

Worldbuilding frame: a COMPUTE measurement of the Vesper climate model on this
desktop. Nothing here is about the simulated planet.

**Interleaved, never blocked.** `notes/audits/aocl-and-model-build-flags.md`
measured a 6% swing in one binary between blocked sessions: whichever arm ran
first after an idle stretch took the boost clock. Alternating A and B within
each round, and flipping which goes first on odd rounds, makes every comparison
paired against a run seconds away and shares any drift equally.

Reports the paired per-round gain, which is the number to quote, alongside the
medians. Also reports each arm's restart sha: if the two differ that is a
numerics change and must be stated separately from the timing, never folded
into it.
"""
from __future__ import annotations

import argparse
import json
import shutil
import statistics as st
import subprocess
import time
from pathlib import Path


def launcher(spec: str, ranks: int) -> tuple[list[str], dict]:
    """How to start one arm, and the environment it needs.

    `mpi`             mpiexec -np <ranks>, which binds rank r to core r.
    `omp[:policy]`    one process; the thread count is compiled in. The two
                      exports are not optional and not tuning: libgomp does NOT
                      bind by default, its placement changes run to run, and
                      with sixteen threads some land on SMT siblings while whole
                      cores sit idle -- which on this processor also randomises
                      which die a thread gets. Bound this way a thread takes
                      core t, exactly as Open MPI gives rank t core t, so the
                      two arms are compared on the same placement rather than on
                      their defaults. `policy` sets OMP_WAIT_POLICY.
    """
    import os
    if spec == "mpi":
        return ["mpiexec", "-np", str(ranks), "./probe_ab.x"], dict(os.environ)
    if spec.split(":")[0] != "omp":
        raise SystemExit(f"unknown launcher {spec!r}; use mpi or omp[:active|passive]")
    env = dict(os.environ)
    env["OMP_STACKSIZE"] = env.get("OMP_STACKSIZE", "512M")
    env["OMP_PROC_BIND"] = "close"
    env["OMP_PLACES"] = "cores"
    if ":" in spec:
        env["OMP_WAIT_POLICY"] = spec.split(":", 1)[1]
    # OMP_STACKSIZE sizes the NON-MASTER threads only; the master runs on the
    # process stack, which is ulimit -s. The model's large locals become
    # stack-allocated under -frecursive, and at T127 the master overruns a
    # 16 MB limit and segfaults. The MPI build never needed this, because
    # without -frecursive those same arrays are static.
    return ["bash", "-c", "ulimit -s unlimited; exec ./probe_ab.x"], env


def run_once(bed: Path, exe: Path, ranks: int, spec: str = "mpi") -> tuple[float, str]:
    local = bed / "probe_ab.x"
    shutil.copy2(exe, local)
    for stale in ("plasim_status", "Abort_Message"):
        (bed / stale).unlink(missing_ok=True)
    cmd, env = launcher(spec, ranks)
    t0 = time.perf_counter()
    r = subprocess.run(cmd, cwd=bed, capture_output=True, text=True, env=env)
    dt = time.perf_counter() - t0
    local.unlink(missing_ok=True)
    if r.returncode != 0 or (bed / "Abort_Message").exists():
        tail = "\n".join(l for l in (r.stderr or r.stdout).splitlines()
                         if "Fortran runtime warning" not in l
                         and not l.startswith("At line"))[-400:]
        raise SystemExit(f"model failed running {exe.name} on {bed.name}:\n{tail}")
    status = bed / "plasim_status"
    sha = subprocess.run(["sha256sum", str(status)], capture_output=True,
                         text=True).stdout[:16] if status.is_file() else "(none)"
    return dt, sha


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", type=Path, required=True)
    ap.add_argument("--a", type=Path, required=True)
    ap.add_argument("--b", type=Path, required=True)
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--ranks", type=int, default=16)
    ap.add_argument("--a-launch", default="mpi",
                    help="mpi, or omp[:active|passive]")
    ap.add_argument("--b-launch", default="mpi",
                    help="mpi, or omp[:active|passive]")
    ap.add_argument("--rounds", type=int, default=6)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    bed, a, b = args.bed.resolve(), args.a.resolve(), args.b.resolve()
    print(f"bed {bed.name}: {args.label_a} vs {args.label_b}, "
          f"{args.rounds} interleaved rounds, {args.ranks} ranks")

    # TWO warm-ups, one per arm. One was not enough: with a single warm-up the
    # first timed round came in 30-40% slow in EVERY arm at T21 and moved the
    # medians, because the machine was still climbing to its boost clock and the
    # bed's inputs were not yet in page cache. This is the same effect
    # `notes/audits/aocl-and-model-build-flags.md` found between blocked
    # sessions, one level down.
    print("warm-up (both arms) ...", flush=True)
    run_once(bed, a, args.ranks, args.a_launch)
    run_once(bed, b, args.ranks, args.b_launch)

    ta: list[float] = []
    tb: list[float] = []
    sa: set[str] = set()
    sb: set[str] = set()
    for rnd in range(1, args.rounds + 1):
        # Flip the order each round so neither arm is always first.
        first_is_a = rnd % 2 == 1
        pair = [(a, ta, sa), (b, tb, sb)] if first_is_a else [(b, tb, sb), (a, ta, sa)]
        for exe, times, shas in pair:
            spec = args.a_launch if exe == a else args.b_launch
            dt, sha = run_once(bed, exe, args.ranks, spec)
            times.append(dt)
            shas.add(sha)
        print(f"  round {rnd:2d}: {args.label_a}={ta[-1]:7.2f}  "
              f"{args.label_b}={tb[-1]:7.2f}  "
              f"({'A first' if first_is_a else 'B first'})", flush=True)

    gains = [100.0 * (x - y) / x for x, y in zip(ta, tb)]
    med_gain = st.median(gains)
    # Reported alongside, never instead: a reader can see whether the answer
    # depends on dropping the first round.
    tail = gains[1:] or gains
    med_tail = st.median(tail)
    result = {
        "bed": bed.name, "ranks": args.ranks, "rounds": args.rounds,
        "a_launch": args.a_launch, "b_launch": args.b_launch,
        "label_a": args.label_a, "label_b": args.label_b,
        "a_median_s": st.median(ta), "b_median_s": st.median(tb),
        "a_times": ta, "b_times": tb,
        "a_sha": sorted(sa), "b_sha": sorted(sb),
        "paired_gain_pct_median": med_gain,
        "paired_gain_pct_min": min(gains), "paired_gain_pct_max": max(gains),
        "b_faster_in_rounds": sum(g > 0 for g in gains),
        "paired_gain_pct_median_excl_round1": med_tail,
        "numerics": "unchanged" if sa == sb else "CHANGED",
    }

    print(f"\n  {args.label_a:12} median {st.median(ta):7.2f} s  "
          f"spread {100 * (max(ta) - min(ta)) / st.median(ta):.1f}%  sha {sorted(sa)[0]}")
    print(f"  {args.label_b:12} median {st.median(tb):7.2f} s  "
          f"spread {100 * (max(tb) - min(tb)) / st.median(tb):.1f}%  sha {sorted(sb)[0]}")
    print(f"  PAIRED GAIN {med_gain:+.2f}%  [{min(gains):+.2f}, {max(gains):+.2f}]  "
          f"faster in {result['b_faster_in_rounds']}/{len(gains)} rounds")
    print(f"  excluding round 1: {med_tail:+.2f}%")
    print(f"  numerics: {result['numerics']}")
    # The guard fires on the GAIN being unreadable against the noise, not on the
    # noise alone. A 19% gain measured on a bed with 7% self-scatter is still a
    # 19% gain; a 1% gain on the same bed is not a result. Firing on scatter
    # alone cried wolf on exactly the arm that mattered most.
    worst = max(100 * (max(t) - min(t)) / st.median(t) for t in (ta, tb))
    if worst > 5.0 and abs(med_gain) < worst:
        print(f"  ** self-scatter {worst:.1f}% exceeds the 5% floor this project "
              f"declares as 'no difference', and the gain ({med_gain:+.2f}%) does "
              f"not clear it. NOT a result; quiet the machine or lengthen the bed.")
    elif worst > 5.0:
        print(f"  (self-scatter {worst:.1f}% is above the 5% floor, but the gain "
              f"{med_gain:+.2f}% is larger than the noise and is readable.)")

    if args.out:
        args.out.write_text(json.dumps(result, indent=2) + "\n")
        print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()
