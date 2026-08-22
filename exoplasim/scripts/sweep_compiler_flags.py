#!/usr/bin/env python3
"""Interleaved A/B of compiler flag sets, timed on a profiling bed.

    python exoplasim/scripts/sweep_compiler_flags.py --bed exoplasim/bench/bed_t127 \
        --resolution T127 --rounds 8

Worldbuilding frame: a COMPUTE measurement of the Vesper climate model on this
desktop. Nothing here is about the simulated planet.

**Arms are interleaved, never blocked**, which is
`notes/audits/aocl-and-model-build-flags.md`'s finding and the reason this is a
script rather than a shell loop. Run in blocks, the stock baseline moved 6%
between sessions while every other arm reproduced to 0.5%: whichever arm ran
first after an idle stretch took the boost clock and the rest ran warm, so
blocking hands the whole artefact to one arm. Round-robin inside each round
makes the comparison paired -- every arm is measured against a baseline that ran
seconds before it -- and a warm-up precedes the first timed run.

**Every arm records the restart sha.** A flag that changes the integration is
reported as a numerics change whatever its time, because that is a different
question from whether it is faster and must not be settled by the same number.

The build is the slow part and is done once per arm, up front. Timing is the
wall time of the `mpiexec` call on the bed, which is model compute: the bed has
its output streams off.
"""
from __future__ import annotations

import argparse
import json
import shutil
import statistics as st
import subprocess
import sys
import time
from pathlib import Path

from _paths import COMPONENT_ROOT, MODEL_RUN, PROJECT_ROOT

PKG = PROJECT_ROOT / "vendor" / "exoplasim" / "exoplasim"

# Arm name -> what it does to the DECLARED flag line in config/planet.yaml.
# `add` and `drop` reach the compiler through build_model.py's --extra-flag and
# --drop-flag, which put both in the build directory's name, so two arms cannot
# land in one directory and reuse each other's objects.
#
# There is no `opt` any more. It existed to feed compile.sh's -O hook, which
# could only APPEND one whitespace-free token to a copy of a generated options
# file; -march=znver4 is in the declared line itself now, so every arm carries it
# by construction rather than by each arm remembering to.
#
# A `drop` naming a flag the declared line does not contain is an ERROR, not a
# no-op: an arm that silently dropped nothing would measure the control and
# report it as the treatment.
ARMS: dict[str, dict] = {
    "stock":        {},
    "no_fcheck":    {"drop": ["-fcheck=all"]},
    "unroll":       {"add": ["-funroll-loops"]},
    "no_fcheck_unroll": {
                         "drop": ["-fcheck=all"], "add": ["-funroll-loops"]},
    # -Ofast implies -ffast-math, which reassociates and would fight
    # -ffpe-trap=invalid; the trap is dropped WITH it so the arm measures
    # codegen rather than dying on the first denormal. It is a numerics change
    # by construction and is here to size the ceiling, not to be adopted.
    "ofast":        {"drop": ["-fcheck=all", "-O3",
                                                     "-ffpe-trap=invalid,zero,overflow"],
                     "add": ["-Ofast"]},
    "prefetch":     {"drop": ["-fcheck=all"],
                     "add": ["-fprefetch-loop-arrays"]},
    "ipa":          {"drop": ["-fcheck=all"],
                     "add": ["-fipa-pta", "-fno-semantic-interposition"]},
}

BASE_LINE_PREFIX = "MOST_F90_OPTS="


def build(arm: str, spec: dict, res: str, layers: int, ranks: int) -> Path | None:
    """Compile one arm and stash the executable under bench/flagsweep/."""
    name = f"most_plasim_t{res.lstrip('Tt').lower()}_l{layers}_p{ranks}.x"
    stash = COMPONENT_ROOT / "bench" / "flagsweep" / arm
    stash.mkdir(parents=True, exist_ok=True)
    target = stash / name
    if target.is_file():
        print(f"  {arm}: already built, reusing")
        return target

    cmd = [sys.executable, str(Path(__file__).resolve().parent / "build_model.py"),
           "--res", res, "--levels", str(layers), "--ranks", str(ranks),
           "--parmode", "mpi"]
    for f in spec.get("drop", []):
        cmd.append(f"--drop-flag={f}")
    for f in spec.get("add", []):
        cmd.append(f"--extra-flag={f}")
    r = subprocess.run(cmd, capture_output=True, text=True)
    built = MODEL_RUN / name
    if r.returncode != 0 or not built.is_file():
        # One arm that will not compile must not take the sweep with it: the
        # remaining arms are still a valid interleaved comparison as long as
        # `stock` survives, and a flag that does not build IS the result for
        # that flag.
        tail = (r.stderr or r.stdout).strip().splitlines()[-3:]
        print(f"  {arm}: BUILD FAILED, arm dropped -- {' / '.join(tail)}")
        return None
    shutil.move(str(built), target)
    print(f"  {arm}: built")
    return target


def time_once(bed: Path, exe: Path, ranks: int) -> tuple[float, str]:
    local = bed / exe.name
    shutil.copy2(exe, local)
    t0 = time.perf_counter()
    subprocess.run(["mpiexec", "-np", str(ranks), f"./{local.name}"],
                   cwd=bed, capture_output=True, text=True)
    dt = time.perf_counter() - t0
    status = bed / "plasim_status"
    sha = subprocess.run(["sha256sum", str(status)], capture_output=True,
                         text=True).stdout[:16] if status.is_file() else "(none)"
    local.unlink()
    return dt, sha


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", type=Path, required=True)
    ap.add_argument("--resolution", required=True)
    ap.add_argument("--layers", type=int, default=10)
    ap.add_argument("--ranks", type=int, default=16)
    ap.add_argument("--rounds", type=int, default=8)
    ap.add_argument("--arms", nargs="*", default=list(ARMS))
    args = ap.parse_args()

    bed = args.bed.resolve()
    arms = {a: ARMS[a] for a in args.arms}

    # No try/finally: nothing to restore. No arm edits a file in the tree any
    # more -- flags reach the compiler as arguments -- so a run that dies leaves
    # nothing behind for whatever builds next, which is what the old route did.
    print("building arms")
    exes = {a: build(a, s, args.resolution, args.layers, args.ranks)
            for a, s in arms.items()}
    exes = {a: e for a, e in exes.items() if e is not None}
    arms = {a: s for a, s in arms.items() if a in exes}
    if "stock" not in arms:
        raise SystemExit("the stock arm failed to build; nothing to compare against")

    print(f"\nwarm-up ({list(arms)[0]})")
    time_once(bed, exes[list(arms)[0]], args.ranks)

    times: dict[str, list[float]] = {a: [] for a in arms}
    shas: dict[str, set] = {a: set() for a in arms}
    order = list(arms)
    for rnd in range(1, args.rounds + 1):
        # Rotate the order each round so no arm is always first in a round.
        rot = order[rnd % len(order):] + order[:rnd % len(order)]
        line = []
        for a in rot:
            dt, sha = time_once(bed, exes[a], args.ranks)
            times[a].append(dt)
            shas[a].add(sha)
            line.append(f"{a}={dt:.2f}")
        print(f"round {rnd:2d}: " + "  ".join(line), flush=True)

    base = st.median(times["stock"]) if "stock" in times else None
    print(f"\n{'arm':18} {'median s':>9} {'vs stock':>9} {'spread':>7}  numerics")
    for a in arms:
        med = st.median(times[a])
        spread = 100 * (max(times[a]) - min(times[a])) / med
        vs = f"{100 * (base - med) / base:+.1f}%" if base else "--"
        same = "unchanged" if shas[a] == shas.get("stock", shas[a]) else "CHANGED"
        print(f"{a:18} {med:9.2f} {vs:>9} {spread:6.1f}%  {same}")

    paired = {}
    if "stock" in times:
        for a in arms:
            if a == "stock":
                continue
            gains = [100 * (s - x) / s for s, x in zip(times["stock"], times[a])]
            paired[a] = {"median_gain_pct": st.median(gains),
                         "wins": sum(g > 0 for g in gains), "rounds": len(gains)}
        print(f"\npaired against the stock run in the SAME round:")
        for a, p in paired.items():
            print(f"  {a:18} {p['median_gain_pct']:+6.2f}%  faster in "
                  f"{p['wins']}/{p['rounds']} rounds")

    out = COMPONENT_ROOT / "bench" / f"flagsweep_{args.resolution.lower()}.json"
    out.write_text(json.dumps({
        "resolution": args.resolution, "layers": args.layers, "ranks": args.ranks,
        "bed": str(bed), "rounds": args.rounds,
        "times": times, "shas": {a: sorted(s) for a, s in shas.items()},
        "paired_vs_stock": paired,
    }, indent=2) + "\n")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
