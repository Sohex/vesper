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

**Two arms is the minimum, not the shape.** `--arm LABEL=PATH`, repeated, runs
any number of arms in one interleaved schedule, and the reason is that a
two-arm run measures one effect with every other held at whatever value the
person building the arms happened to pick. When two changes are believed to
COMPOUND, the arbitrary setting is the answer: measuring change X against a
baseline that already has change Y is a different number from measuring it
against one that does not, and neither on its own says whether X and Y add,
multiply, or interfere.

    python exoplasim/scripts/bench_ab.py --bed <bed> --rounds 8 --factorial \
        --arm "fp64 masked=<x>" --arm "fp64 vector=<x>" \
        --arm "fp32 masked=<x>" --arm "fp32 vector=<x>"

`--factorial` declares that the four arms are, in that order, the base, the
first factor alone, the second factor alone, and both. It then reports each
factor's gain WITHIN each level of the other, and the excess of the joint gain
over what multiplying the two single-factor gains predicts. That excess is the
whole answer to "do they compound", and it is computed per round and reported
with its spread, because a median excess of two points means nothing beside a
round-to-round spread of ten.

The rotation is by round: with n arms, round r starts at arm r mod n, so over
any n consecutive rounds each arm occupies each position once and the drift
that position carries is shared equally. Two arms recover the original A/B/B/A
flip exactly.

**When the host cannot be made quiet, change the instrument.** `--counter
instructions:u` runs the same schedule under `perf stat` and compares retired
instructions instead of seconds. A counter is a property of the code and the
input and does not move with what else the machine is doing; a clock is a
property of the machine as well, and a contended clock is worse than no number
because it looks like a measurement. The counter cannot say what a change
COSTS -- a vector call is more instruction-efficient than eight scalar ones and
need not be eight times faster -- so it answers where the work went, and a clock
arm on a quiet host answers what that was worth. Report both or say which one
is missing.
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


# Set from --counter. Module level because it changes what run_once MEASURES,
# not how one arm is launched, and every arm must be measured the same way.
COUNTER: str | None = None


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
    cmd = ["bash", "-c", "ulimit -s unlimited; exec ./probe_ab.x"]
    if COUNTER:
        # A HARDWARE COUNTER, NOT A CLOCK. A retired-instruction count is a
        # property of the code and the input: it does not move with what else
        # the machine is doing, and on this project it reproduces to about one
        # part in 10^8 across repeats, where a wall clock on a shared desktop
        # does not. `exoplasim/notes/radiation-scheme-price.md` argues that
        # choice; `scripts/machine.py` states the same thing from the other
        # side. WHAT IT CANNOT SAY is cycles: a vector call is more
        # instruction-efficient than eight scalar ones and need not be eight
        # times faster, because the limit may be elsewhere. So a counter arm
        # settles WHERE THE WORK WENT and a clock arm settles what it cost, and
        # neither substitutes for the other.
        #
        # OMP_WAIT_POLICY=passive unless the caller set one: with the default
        # active policy an idle thread SPINS, and the instructions it retires
        # spinning are a measurement of the host's scheduler rather than of the
        # model.
        env.setdefault("OMP_WAIT_POLICY", "passive")
        cmd = ["perf", "stat", "-x,", "-e", COUNTER, "--"] + cmd
    return cmd, env


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
    if COUNTER:
        # perf -x, writes "<value>,<unit>,<event>,..." on stderr, and a counter
        # that could not be read writes "<not counted>" in the value field. That
        # is a failed measurement and never a zero.
        got = None
        for line in (r.stderr or "").splitlines():
            f = line.split(",")
            if len(f) > 2 and f[2].startswith(COUNTER.split(":")[0]):
                if not f[0].replace(".", "").isdigit():
                    raise SystemExit(f"perf could not count {COUNTER}: {line}")
                got = float(f[0])
                break
        if got is None:
            raise SystemExit(f"perf printed no {COUNTER} row:\n{r.stderr[-400:]}")
        dt = got
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
    ap.add_argument("--a", type=Path)
    ap.add_argument("--b", type=Path)
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--arm", action="append", default=[], metavar="LABEL=PATH",
                    help="one arm, repeatable; replaces --a/--b. Every --arm "
                         "shares --launch and --nl, so the arms differ by the "
                         "executable and by nothing else.")
    ap.add_argument("--launch", default="omp",
                    help="launcher for every --arm; see --a-launch")
    ap.add_argument("--nl", action="append", default=[], metavar="KEY=VALUE",
                    help="namelist setting forced for every --arm; repeatable")
    ap.add_argument("--factorial", action="store_true",
                    help="the four --arm entries are base, factor A, factor B, "
                         "both, in that order; report whether the two compound")
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
    ap.add_argument("--counter", metavar="EVENT", default=None,
                    help="measure a perf hardware counter instead of the wall "
                         "clock, e.g. instructions:u. The whole paired and "
                         "interleaved schedule is unchanged; only what is "
                         "compared changes. Use this when the host cannot be "
                         "made quiet: a counter is a property of the code and "
                         "the input, a clock is a property of the machine too.")
    ap.add_argument("--outlier-factor", type=float, default=3.0,
                    help="a round in which any arm exceeds this multiple of "
                         "that arm's own median is a MACHINE EVENT, not a "
                         "measurement, and is dropped whole -- every arm of it "
                         "together, so the pairing survives. The dropped rounds "
                         "and their times are always reported. 0 disables it.")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    global COUNTER
    COUNTER = args.counter
    unit = "" if COUNTER else " s"
    fmt = "12.4g" if COUNTER else "7.2f"
    bed = args.bed.resolve()
    if args.arm:
        if args.a or args.b:
            raise SystemExit("--arm replaces --a/--b; do not mix the two forms")
        arms = []
        for spec in args.arm:
            label, sep, path = spec.partition("=")
            if not sep:
                raise SystemExit(f"--arm {spec!r} is not LABEL=PATH")
            arms.append((label.strip(), Path(path).resolve(),
                         args.launch, list(args.nl)))
    else:
        if not (args.a and args.b):
            raise SystemExit("give --a and --b, or two or more --arm LABEL=PATH")
        arms = [(args.label_a, args.a.resolve(), args.a_launch, args.a_nl),
                (args.label_b, args.b.resolve(), args.b_launch, args.b_nl)]
    if args.factorial and len(arms) != 4:
        raise SystemExit("--factorial takes exactly four arms: base, factor A, "
                         f"factor B, both. Got {len(arms)}.")
    n = len(arms)
    print(f"bed {bed.name}: {', '.join(l for l, *_ in arms)}, "
          f"{args.rounds} interleaved rounds, {args.threads} threads")

    # ONE WARM-UP PER ARM. One was not enough: with a single warm-up the first
    # timed round came in 30-40% slow in EVERY arm at T21 and moved the medians,
    # because the machine was still climbing to its boost clock and the bed's
    # inputs were not yet in page cache. This is the same effect
    # `notes/audits/aocl-and-model-build-flags.md` found between blocked
    # sessions, one level down.
    print("warm-up (every arm) ...", flush=True)
    for _label, exe, spec, settings in arms:
        run_once(bed, exe, args.threads, spec, settings)

    times: list[list[float]] = [[] for _ in arms]
    shas: list[set[str]] = [set() for _ in arms]
    for rnd in range(1, args.rounds + 1):
        # Rotate the starting arm each round so no arm is always first and the
        # drift each position carries is shared equally over any n rounds.
        order = [(rnd - 1 + i) % n for i in range(n)]
        for i in order:
            _label, exe, spec, settings = arms[i]
            dt, sha = run_once(bed, exe, args.threads, spec, settings)
            times[i].append(dt)
            shas[i].add(sha)
        print(f"  round {rnd:2d}: "
              + "  ".join(f"{arms[i][0]}={times[i][-1]:{fmt}}" for i in range(n))
              + f"   (started {arms[order[0]][0]})", flush=True)

    # A ROUND IS DROPPED WHOLE OR NOT AT ALL. Another job landing on the host
    # for one round does not slow the arm that happened to be running by ten per
    # cent, it slows it by a factor; and it lands on ONE arm, so the pairing
    # that defeats drift cannot share it. Dropping the whole round keeps every
    # remaining comparison paired, and dropping on a FACTOR rather than on a
    # rank means the rule is about a different phenomenon rather than about the
    # tail of this one. The rule is stated as a multiple so it can be fixed
    # before the run; what it removed is printed either way, because a rule that
    # silently discards data is a rule that can be tuned after the fact.
    dropped: list[int] = []
    raw = [list(t) for t in times]
    if args.outlier_factor > 0:
        med0 = [st.median(t) for t in times]
        dropped = [r for r in range(args.rounds)
                   if any(times[i][r] > args.outlier_factor * med0[i]
                          for i in range(n))]
        if dropped:
            print(f"\n  DROPPED {len(dropped)} round(s) as machine events, each "
                  f"holding an arm above {args.outlier_factor:g}x its own median:")
            for r in dropped:
                print("    round {:2d}: ".format(r + 1)
                      + "  ".join(f"{arms[i][0]}={times[i][r]:.2f}"
                                  for i in range(n)))
            keep = [r for r in range(args.rounds) if r not in dropped]
            if len(keep) < 2:
                raise SystemExit("fewer than two rounds survive the outlier "
                                 "rule; the host was not quiet enough to "
                                 "measure on. Nothing is reported.")
            times = [[t[r] for r in keep] for t in times]

    def scatter(t: list[float]) -> float:
        return 100.0 * (max(t) - min(t)) / st.median(t)

    ref = 0
    gains = [[100.0 * (r - x) / r for r, x in zip(times[ref], times[i])]
             for i in range(n)]
    result = {
        "bed": bed.name, "threads": args.threads, "rounds": args.rounds,
        "measured": COUNTER or "wall_clock_seconds",
        "outlier_factor": args.outlier_factor,
        "rounds_dropped": [r + 1 for r in dropped],
        "rounds_kept": args.rounds - len(dropped),
        "reference_arm": arms[ref][0],
        "arms": [
            {"label": arms[i][0], "exe": str(arms[i][1]),
             "launch": arms[i][2], "namelist": arms[i][3],
             "times_s": times[i], "times_all_rounds_s": raw[i],
             "median_s": st.median(times[i]),
             "self_scatter_pct": scatter(times[i]), "sha": sorted(shas[i]),
             "paired_gain_pct_median": st.median(gains[i]),
             "paired_gain_pct_min": min(gains[i]),
             "paired_gain_pct_max": max(gains[i]),
             "faster_in_rounds": sum(g > 0 for g in gains[i]),
             "numerics": ("reference" if i == ref else
                          "unchanged" if shas[i] == shas[ref] else "CHANGED")}
            for i in range(n)],
    }
    if len(arms) == 2:
        # The keys the two-arm callers and the audit notes already read.
        result.update({
            "label_a": arms[0][0], "label_b": arms[1][0],
            "a_launch": arms[0][2], "b_launch": arms[1][2],
            "a_median_s": st.median(times[0]), "b_median_s": st.median(times[1]),
            "a_times": times[0], "b_times": times[1],
            "a_sha": sorted(shas[0]), "b_sha": sorted(shas[1]),
            "paired_gain_pct_median": st.median(gains[1]),
            "paired_gain_pct_min": min(gains[1]),
            "paired_gain_pct_max": max(gains[1]),
            "b_faster_in_rounds": sum(g > 0 for g in gains[1]),
            "paired_gain_pct_median_excl_round1":
                st.median(gains[1][1:] or gains[1]),
            "numerics": "unchanged" if shas[0] == shas[1] else "CHANGED",
        })

    print()
    for i in range(n):
        print(f"  {arms[i][0]:16} median {st.median(times[i]):{fmt}}{unit}  "
              f"self-scatter {scatter(times[i]):5.2f}%  sha {sorted(shas[i])[0]}")
    print(f"\n  paired gain against {arms[ref][0]}:")
    for i in range(n):
        if i == ref:
            continue
        g = gains[i]
        print(f"    {arms[i][0]:16} {st.median(g):+7.2f}%  "
              f"[{min(g):+.2f}, {max(g):+.2f}]  "
              f"faster in {sum(x > 0 for x in g)}/{len(g)} rounds  "
              f"numerics {result['arms'][i]['numerics']}")

    if args.factorial:
        # base = 0, factor A alone = 1, factor B alone = 2, both = 3.
        # Each factor's gain WITHIN each level of the other, per round, so the
        # question "is X worth more once Y is in place" is answered by two
        # paired numbers rather than by differencing two medians.
        def paired(i: int, j: int) -> list[float]:
            return [100.0 * (x - y) / x for x, y in zip(times[i], times[j])]
        a_in_base, a_in_b = paired(0, 1), paired(2, 3)
        b_in_base, b_in_a = paired(0, 2), paired(1, 3)
        joint = gains[3]
        # What multiplying the two single-factor gains predicts for the joint
        # arm, per round. Excess above zero is compounding beyond multiplicative;
        # below zero is interference.
        pred = [100.0 * (1.0 - (1.0 - x / 100.0) * (1.0 - y / 100.0))
                for x, y in zip(a_in_base, b_in_base)]
        excess = [j - p for j, p in zip(joint, pred)]
        result["factorial"] = {
            "factor_a_label": arms[1][0], "factor_b_label": arms[2][0],
            "a_in_base_pct": st.median(a_in_base), "a_in_b_pct": st.median(a_in_b),
            "b_in_base_pct": st.median(b_in_base), "b_in_a_pct": st.median(b_in_a),
            "joint_pct": st.median(joint),
            "multiplicative_prediction_pct": st.median(pred),
            "excess_over_multiplicative_pct": st.median(excess),
            "excess_min": min(excess), "excess_max": max(excess),
            "a_in_base_rounds": a_in_base, "a_in_b_rounds": a_in_b,
            "b_in_base_rounds": b_in_base, "b_in_a_rounds": b_in_a,
            "excess_rounds": excess,
        }
        print(f"\n  factorial, {arms[1][0]} x {arms[2][0]}:")
        print(f"    {arms[1][0]:16} alone           {st.median(a_in_base):+7.2f}%  "
              f"[{min(a_in_base):+.2f}, {max(a_in_base):+.2f}]")
        print(f"    {arms[1][0]:16} with {arms[2][0]:10} {st.median(a_in_b):+7.2f}%  "
              f"[{min(a_in_b):+.2f}, {max(a_in_b):+.2f}]")
        print(f"    {arms[2][0]:16} alone           {st.median(b_in_base):+7.2f}%  "
              f"[{min(b_in_base):+.2f}, {max(b_in_base):+.2f}]")
        print(f"    {arms[2][0]:16} with {arms[1][0]:10} {st.median(b_in_a):+7.2f}%  "
              f"[{min(b_in_a):+.2f}, {max(b_in_a):+.2f}]")
        print(f"    both together                    {st.median(joint):+7.2f}%  "
              f"[{min(joint):+.2f}, {max(joint):+.2f}]")
        print(f"    multiplying the two alone predicts {st.median(pred):+.2f}%")
        print(f"    EXCESS over that prediction      {st.median(excess):+7.2f}%  "
              f"[{min(excess):+.2f}, {max(excess):+.2f}]")

    # The guard fires on the GAIN being unreadable against the noise, not on the
    # noise alone. A 19% gain measured on a bed with 7% self-scatter is still a
    # 19% gain; a 1% gain on the same bed is not a result. Firing on scatter
    # alone cried wolf on exactly the arm that mattered most.
    worst = max(scatter(t) for t in times)
    smallest = min(abs(st.median(gains[i])) for i in range(n) if i != ref)
    if worst > 5.0 and smallest < worst:
        print(f"\n  ** self-scatter {worst:.1f}% exceeds the 5% floor this project "
              f"declares as 'no difference', and the smallest gain "
              f"({smallest:.2f}%) does not clear it. NOT a result for that arm; "
              f"quiet the machine or lengthen the bed.")
    elif worst > 5.0:
        print(f"\n  (self-scatter {worst:.1f}% is above the 5% floor, but every "
              f"gain is larger than the noise and is readable.)")
    if args.factorial:
        exc = result["factorial"]
        span = exc["excess_max"] - exc["excess_min"]
        if abs(exc["excess_over_multiplicative_pct"]) < span:
            print(f"  ** the excess ({exc['excess_over_multiplicative_pct']:+.2f}%) is "
                  f"smaller than its own round-to-round spread ({span:.2f}%). "
                  f"The two factors are consistent with multiplying; a departure "
                  f"from that is NOT resolved by this bed.")

    if args.out:
        args.out.write_text(json.dumps(result, indent=2) + "\n")
        print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()
