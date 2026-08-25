#!/usr/bin/env python3
"""Paired, interleaved A/B of two executables on one bed.

    python exoplasim/scripts/bench_ab.py --bed exoplasim/bench/bed_t85 \
        --a exoplasim/bench/ref/unpatched_<rung>.x \
        --b exoplasim/bench/ref/patched_<rung>.x --rounds 6

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
import re
import shutil
import statistics as st
import subprocess
import time
from pathlib import Path


def launcher(spec: str, threads: int) -> tuple[list[str], dict]:
    """How to start one arm, and the environment it needs.

    There is one build configuration to launch. The model is a threaded build,
    the thread count is compiled into the binary, and `threads` here only says
    how long an `omp@` core list has to be.

    `omp[:policy]`    one process. The two exports are not optional and not
                      tuning: libgomp does NOT bind by default, its placement
                      changes run to run, and with sixteen threads some land on
                      SMT siblings while whole cores sit idle -- which on this
                      processor also randomises which die a thread gets. Bound
                      this way a thread takes core t, which is what
                      `exoplasim/__init__.py` launches production under, so the
                      two arms are compared on production's placement rather
                      than on their defaults. `policy` sets OMP_WAIT_POLICY.
    `omp@<cores>`     as `omp`, but with an explicit thread-to-core order:
                      a comma-separated core list, thread t taking the t'th
                      entry. Each core expands to both its SMT siblings, so a
                      thread has exactly the freedom `OMP_PLACES=cores` gives
                      it and only the ORDER differs from the control. The
                      decomposition is unchanged, so an arm that differs only
                      in placement must stay bit identical.
    """
    import os
    head = spec.split(":")[0].split("@")[0]
    if head != "omp":
        raise SystemExit(f"unknown launcher {spec!r}; use omp[:active|passive] or omp@<cores>")
    env = dict(os.environ)
    env["OMP_STACKSIZE"] = env.get("OMP_STACKSIZE", "512M")
    env["OMP_PROC_BIND"] = "close"
    if "@" in spec.split(":")[0]:
        order = [int(c) for c in spec.split(":")[0].split("@", 1)[1].split(",")]
        if len(order) != threads or sorted(order) != sorted(set(order)):
            raise SystemExit(f"omp@ needs {threads} distinct cores, got {len(order)}")
        nsib = os.cpu_count() // 2
        env["OMP_PLACES"] = ",".join(f"{{{c},{c + nsib}}}" for c in order)
    else:
        env["OMP_PLACES"] = "cores"
    if ":" in spec:
        env["OMP_WAIT_POLICY"] = spec.split(":", 1)[1]
    # OMP_STACKSIZE sizes the NON-MASTER threads only; the master runs on the
    # process stack, which is ulimit -s. The model's large locals become
    # stack-allocated under -frecursive, and at T127 the master overruns a
    # 16 MB limit and segfaults.
    return ["bash", "-c", "ulimit -s unlimited; exec ./probe_ab.x"], env


def set_namelist(bed: Path, settings: list[str]) -> None:
    """Force KEY=VALUE into plasim_namelist, in place.

    Some arms differ by a NAMELIST switch rather than by an executable -- NSHTNS
    is the reason this exists, and it is a switch precisely so that the two
    transforms can be timed without a second build in the comparison. Setting it
    per arm keeps that property: the two arms then differ by the one thing under
    test and share every compiler decision.

    Rewritten before every run rather than once at setup, because the arms share
    one bed and whichever ran last would otherwise leave its setting behind.
    """
    nl = bed / "plasim_namelist"
    lines = nl.read_text().splitlines()
    for kv in settings:
        key, _, val = kv.partition("=")
        key, val = key.strip().upper(), val.strip()
        pat = re.compile(rf"^\s*{re.escape(key)}\s*=", re.IGNORECASE)
        lines = [l for l in lines if not pat.match(l)]
        lines.insert(1, f" {key} = {val} ")
    nl.write_text("\n".join(lines) + "\n")


def run_once(bed: Path, exe: Path, threads: int, spec: str = "omp",
             settings: list[str] | None = None) -> tuple[float, str]:
    local = bed / "probe_ab.x"
    shutil.copy2(exe, local)
    if settings:
        set_namelist(bed, settings)
    for stale in ("plasim_status", "Abort_Message"):
        (bed / stale).unlink(missing_ok=True)
    cmd, env = launcher(spec, threads)
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
    ap.add_argument("--threads", type=int, default=16,
                    help="the thread count the binaries were compiled for; "
                         "only an omp@ core list is checked against it")
    ap.add_argument("--a-launch", default="omp",
                    help="omp[:active|passive], or omp@<core order>")
    ap.add_argument("--b-launch", default="omp",
                    help="omp[:active|passive], or omp@<core order>")
    ap.add_argument("--a-nl", action="append", default=[], metavar="KEY=VALUE",
                    help="namelist setting forced for arm A; repeatable")
    ap.add_argument("--b-nl", action="append", default=[], metavar="KEY=VALUE",
                    help="namelist setting forced for arm B; repeatable")
    ap.add_argument("--rounds", type=int, default=6)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    bed, a, b = args.bed.resolve(), args.a.resolve(), args.b.resolve()
    print(f"bed {bed.name}: {args.label_a} vs {args.label_b}, "
          f"{args.rounds} interleaved rounds, {args.threads} threads")

    # TWO warm-ups, one per arm. One was not enough: with a single warm-up the
    # first timed round came in 30-40% slow in EVERY arm at T21 and moved the
    # medians, because the machine was still climbing to its boost clock and the
    # bed's inputs were not yet in page cache. This is the same effect
    # `notes/audits/aocl-and-model-build-flags.md` found between blocked
    # sessions, one level down.
    print("warm-up (both arms) ...", flush=True)
    run_once(bed, a, args.threads, args.a_launch, args.a_nl)
    run_once(bed, b, args.threads, args.b_launch, args.b_nl)

    ta: list[float] = []
    tb: list[float] = []
    sa: set[str] = set()
    sb: set[str] = set()
    # An arm carries its own launcher and namelist rather than being looked up
    # by executable, because the two arms may BE the same executable -- which is
    # the point when what differs is a namelist switch.
    arm_a = (a, args.a_launch, args.a_nl, ta, sa)
    arm_b = (b, args.b_launch, args.b_nl, tb, sb)
    for rnd in range(1, args.rounds + 1):
        # Flip the order each round so neither arm is always first.
        first_is_a = rnd % 2 == 1
        pair = [arm_a, arm_b] if first_is_a else [arm_b, arm_a]
        for exe, spec, settings, times, shas in pair:
            dt, sha = run_once(bed, exe, args.threads, spec, settings)
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
        "bed": bed.name, "threads": args.threads, "rounds": args.rounds,
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
