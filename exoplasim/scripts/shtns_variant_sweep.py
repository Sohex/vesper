#!/usr/bin/env python3
"""Which on-the-fly variant is fastest per transform type, per resolution?

    python exoplasim/scripts/shtns_variant_sweep.py --rounds 4 --reps 2001

Worldbuilding frame: a COMPUTE measurement of the transform library the Vesper
climate model calls. Nothing here is about the simulated planet.

CLIM-74. `exoplasim/notes/shtns-algorithm-selection.md` establishes what this
measures and why it is measured HERE rather than on the model: SHTns 3.7.5 has
no matrix-based algorithms, so the choice is on-the-fly unroll depth, and its
model-level worth is of order 0.1% against a bench floor of 5%. The transform
is the only level at which the effect is larger than the instrument --
`docs/src/practice/failure-modes.md` class 34.

THE PICK IS PER RUNG AND PER TYPE. It is not a search for one best setting: the
candidate set itself differs by resolution, since `init_sht_array_func` caps
`alg_lim` by `nlat_2 / VSIZE2`. A rung disagreeing with its neighbour is an
expected outcome, not an inconsistency to be smoothed away.

THE THRESHOLD IS DECLARED HERE AND WAS FIXED BEFORE THIS SWEEP RAN. A variant
other than the fly2 default is adopted for a (rung, type) only when it wins the
median in EVERY round and its median gain over fly2 exceeds MIN_GAIN. The pilot
in the note measured run-to-run scatter at 2 to 4% of a median, and 5% clears
it. Anything under that stays on fly2, which is what `SHT_QUICK_INIT` pins with
no file and no startup race.

INTERLEAVED, NEVER BLOCKED, for the reason `bench_ab.py` documents: a machine
that drifts hands the whole of one arm to a different clock. Each round runs
every variant once and rotates which goes first.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _paths  # noqa: F401  -- anchors ROOT below the way every script does

ROOT = Path(__file__).resolve().parents[2]
PREFIX = ROOT / "vendor" / "shtns-install"
PROBE_SRC = ROOT / "exoplasim" / "scripts" / "probe_shtns_variant_cost.c"
GEN_MODE = "gauss"

# The ladder, as build_model.py declares it. NLON is 2*NLAT throughout.
RUNGS = [("T21", 21, 32), ("T42", 42, 64), ("T85", 85, 128),
         ("T127", 127, 192), ("T170", 170, 256)]
VARIANTS = ["fly1", "fly2", "fly3", "fly4"]
DEFAULT = "fly2"          # what SHT_QUICK_INIT pins, with no file involved
MIN_GAIN = 5.0            # percent, declared before the sweep; see the docstring

RESULT = re.compile(r"RESULT (\w+) \S+ nlat=(\d+) median_us=([\d.]+) "
                    r"iqr_us=([\d.]+) reps=(\d+)")


def build(work: Path) -> Path:
    exe = work / "probe_cost"
    cmd = ["gcc", "-O2", f"-I{PREFIX}/include", str(PROBE_SRC), "-o", str(exe),
           f"-L{PREFIX}/lib", "-lshtns_omp", "-lfftw3_omp", "-lfftw3", "-lm",
           "-fopenmp"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"probe build failed:\n{r.stderr}")
    return exe


def seed_config(work: Path, ntru: int, nlat: int, nlon: int) -> str:
    """One timed initialisation, purely to obtain a well-formed config line.

    The line's ALGORITHM NAMES are overwritten below, so which variants this
    race happened to pick does not reach the result. What is taken from it is
    the grid parameters, the request flags and the SIMD id, none of which this
    script should be reconstructing by hand.
    """
    src = work / "gencfg.c"
    src.write_text(
        "#include <stdlib.h>\n#include <shtns.h>\n"
        "int main(int c, char**v){ shtns_verbose(0); shtns_use_threads(1);\n"
        "  shtns_cfg s = shtns_create(atoi(v[1]), atoi(v[1]), 1, sht_orthonormal);\n"
        "  shtns_set_grid(s, sht_gauss | SHT_PHI_CONTIGUOUS | SHT_LOAD_SAVE_CFG,"
        " 0.0, atoi(v[2]), atoi(v[3])); return 0; }\n")
    exe = work / "gencfg"
    r = subprocess.run(["gcc", "-O2", f"-I{PREFIX}/include", str(src), "-o", str(exe),
                        f"-L{PREFIX}/lib", "-lshtns_omp", "-lfftw3_omp", "-lfftw3",
                        "-lm", "-fopenmp"], capture_output=True, text=True)
    if r.returncode:
        raise SystemExit(f"config seeder build failed:\n{r.stderr}")
    cfg = work / "shtns_cfg"
    cfg.unlink(missing_ok=True)
    subprocess.run([str(exe), str(ntru), str(nlat), str(nlon)], cwd=work, check=True,
                   capture_output=True)
    return cfg.read_text()


def run_arm(exe: Path, work: Path, template: str, variant: str,
            ntru: int, nlat: int, nlon: int, reps: int, procs: int) -> dict:
    """One variant, forced and then CONFIRMED from what the library reports.

    A forced arm that silently fell back to the default would compare fly2
    against itself and report a clean null, which is the failure this check
    exists to make impossible.
    """
    (work / "shtns_cfg").write_text(re.sub(r"fly[0-9]", variant, template))
    # SEVERAL PROCESSES PER ARM, and this is not belt-and-braces. The first
    # sweep found the T42 `syn` arms BIMODAL between processes -- fly2 landing
    # at 8.1 us with an IQR of 0.16 in one round and 13.2 us with an IQR of
    # 4.61 in another, while fly1 and fly4 held to a fifth of a microsecond in
    # the same rounds. Buffers are allocated once per process, so an arm-level
    # median inherits ONE draw of whatever placement that process got, and a
    # comparison between single draws measures the draw. Taking the median over
    # processes is what makes the arm a measurement of the variant.
    runs, checked = [], None
    for _ in range(procs):
        r = subprocess.run([str(exe), str(ntru), str(nlat), str(nlon), str(reps)],
                           cwd=work, capture_output=True, text=True)
        if r.returncode:
            raise SystemExit(f"probe failed at {nlat}/{variant}:\n{r.stderr}")
        runs.append(r.stdout)
        checked = checked or r.stdout
    loaded = [ln for ln in checked.splitlines() if "std:" in ln]
    if not loaded:
        raise SystemExit(f"{nlat}/{variant}: the library reported no ftable to check")
    # POSITIONAL, because `fprint_ftable` writes the types in `sht_type`
    # order -- syn, ana, vsy, van, gsp, gto, v3s, v3a -- and only the first
    # five are the calls shtnsmod.f90 makes and this probe times. The trailing
    # entries legitimately stay on `s+v` where a fly variant does not apply,
    # so a whole-line set comparison would read that as a failed arm.
    entries = loaded[0].split(":", 1)[1].split()
    timed = entries[:5]
    if any(e != variant for e in timed):
        # NOT AVAILABLE AT THIS RUNG, and that is data rather than a failure.
        # `init_sht_array_func` caps `alg_lim` by `nlat_2 / VSIZE2`, so the
        # deeper unrolls do not exist on the small grids; `config_load` accepts
        # only non-null pointers, so an unavailable name leaves the default in
        # place. Detecting it is what stops the sweep comparing the default
        # against itself and reporting a tidy null, which is why the request is
        # confirmed rather than assumed.
        if all(e == DEFAULT for e in timed):
            return None
        raise SystemExit(
            f"{nlat}/{variant}: asked for {variant} and the library loaded "
            f"{timed} for the timed types, which is neither the request nor "
            f"the default throughout.")
    per_type: dict[str, list[float]] = {}
    iqrs: dict[str, list[float]] = {}
    for text in runs:
        for m in RESULT.finditer(text):
            per_type.setdefault(m.group(1), []).append(float(m.group(3)))
            iqrs.setdefault(m.group(1), []).append(float(m.group(4)))
    out = {}
    for typ, xs in per_type.items():
        xs = sorted(xs)
        out[typ] = {"median_us": xs[len(xs) // 2],
                    # The SPREAD ACROSS PROCESSES, reported rather than
                    # averaged away: it is the scatter any gain has to clear,
                    # and hiding it is what made the first sweep's +35.6% look
                    # like a result.
                    "proc_spread_pct": round(100.0 * (xs[-1] - xs[0]) / xs[len(xs) // 2], 2),
                    "procs": len(xs),
                    "iqr_us": sorted(iqrs[typ])[len(iqrs[typ]) // 2]}
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", type=int, default=4)
    ap.add_argument("--procs", type=int, default=5,
                    help="processes per arm per round; the arm's value is the "
                         "median over them. See run_arm on why one is not enough.")
    ap.add_argument("--reps", type=int, default=2001,
                    help="transform calls timed per type per arm per round")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "exoplasim" / "analysis" / "shtns_variant_sweep.json")
    args = ap.parse_args()

    work = Path(subprocess.run(["mktemp", "-d"], capture_output=True, text=True)
                .stdout.strip())
    exe = build(work)
    rounds: list[dict] = []
    unavailable: set[tuple[str, str]] = set()
    for rnd in range(args.rounds):
        # ROTATED, so no variant keeps the first slot and its warm-up cost.
        order = VARIANTS[rnd % len(VARIANTS):] + VARIANTS[:rnd % len(VARIANTS)]
        for name, ntru, nlat in RUNGS:
            template = seed_config(work, ntru, nlat, 2 * nlat)
            for variant in order:
                got = run_arm(exe, work, template, variant, ntru, nlat, 2 * nlat,
                              args.reps, args.procs)
                if got is None:
                    unavailable.add((name, variant))
                    continue
                for typ, rec in got.items():
                    rounds.append({"round": rnd, "rung": name, "nlat": nlat,
                                   "type": typ, "variant": variant, **rec})
        print(f"round {rnd + 1}/{args.rounds} done", flush=True)

    # THE VERDICT, by the rule declared at the top of this file.
    picks = {}
    for name, _, nlat in RUNGS:
        for typ in sorted({r["type"] for r in rounds if r["rung"] == name}):
            def med(v):
                xs = [r["median_us"] for r in rounds
                      if r["rung"] == name and r["type"] == typ and r["variant"] == v]
                return sorted(xs)[len(xs) // 2] if xs else None
            base = med(DEFAULT)
            best, gain, every = DEFAULT, 0.0, False
            for v in VARIANTS:
                if v == DEFAULT or (name, v) in unavailable:
                    continue
                g = 100.0 * (base - med(v)) / base
                wins = all(
                    min(r["median_us"] for r in rounds
                        if r["rung"] == name and r["type"] == typ
                        and r["variant"] == v and r["round"] == rr)
                    < min(r["median_us"] for r in rounds
                          if r["rung"] == name and r["type"] == typ
                          and r["variant"] == DEFAULT and r["round"] == rr)
                    for rr in range(args.rounds))
                if g > gain:
                    best, gain, every = v, g, wins
            adopt = best != DEFAULT and gain > MIN_GAIN and every
            picks[f"{name}/{typ}"] = {
                "pick": best if adopt else DEFAULT,
                "gain_pct_vs_default": round(gain, 2),
                "wins_every_round": every,
                "adopted": adopt,
                "candidates": [v for v in VARIANTS if (name, v) not in unavailable],
                "medians_us": {v: med(v) for v in VARIANTS
                               if (name, v) not in unavailable}}

    payload = {
        "note": "Which SHTns on-the-fly variant is fastest per transform type "
                "per resolution. CLIM-74; see "
                "exoplasim/notes/shtns-algorithm-selection.md.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "threshold": {"min_gain_pct": MIN_GAIN, "default": DEFAULT,
                      "rule": "adopt only if it beats the default in every "
                              "round AND by more than min_gain_pct",
                      "fixed_before_the_sweep": True},
        "provenance": {
            "shtns_prefix": str(PREFIX),
            "probe_sha256": hashlib.sha256(PROBE_SRC.read_bytes()).hexdigest(),
            "rounds": args.rounds, "reps": args.reps,
            "host": platform.node(), "machine": platform.machine()},
        "unavailable_by_rung": {r: sorted(v for n, v in unavailable if n == r)
                                for r, _, _ in RUNGS},
        "picks": picks,
        "rounds_raw": rounds}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2))
    shutil.rmtree(work, ignore_errors=True)

    print(f"\nwrote {args.out.relative_to(ROOT)}")
    adopted = [k for k, v in picks.items() if v["adopted"]]
    print(f"{len(adopted)} of {len(picks)} (rung, type) pairs beat {DEFAULT} "
          f"by more than {MIN_GAIN}% in every round"
          + (":" if adopted else "; every one stays on the default."))
    for k in adopted:
        print(f"  {k}: {picks[k]['pick']}  +{picks[k]['gain_pct_vs_default']}%")


if __name__ == "__main__":
    main()
