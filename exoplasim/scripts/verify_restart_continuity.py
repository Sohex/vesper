#!/usr/bin/env python3
"""Is a segmented run the same experiment as the same run taken whole?

    python exoplasim/scripts/verify_restart_continuity.py --exe <plasim.x> --bed <dir>
    python exoplasim/scripts/verify_restart_continuity.py --exe <plasim.x> --bed <dir> --round-trip

Worldbuilding frame: this checks the Vesper project's climate model against
itself. Nothing here is a claim about the simulated planet.

WHAT IT TESTS, and why it is not the reproducibility matrix.
`reproducibility_matrix.py` asks whether the model gives the same bytes twice
from ONE restart, and the answer recorded in
`notes/audits/model-reproducibility.md` is yes. That is a different property
from the one every segmented run depends on: that N steps taken as one model
call and N steps taken as two, with the first call's `plasim_status` handed to
the second as its `plasim_restart`, end in the same state. A model can be
perfectly reproducible per call and still integrate a different weather either
side of every boundary, and this one was: `glacierini` scattered the current
surface pressure into the leapfrog minus level on the restart path as well as
the cold one, so the first step after every segment boundary was a step the
uninterrupted run never took. world-8yyh.

TWO MODES, and the second is the one that localises a failure.

  default       N steps whole against `--split` + (N - split), comparing the
                two final restarts record by record. This is the property
                callers actually need, and it is what fails first.
  --round-trip  N steps, then a restart from that state taking NO STEPS. No
                integration happens between the two `plasim_status` files, so
                every record that differs is state the write-and-read of a
                restart does not carry, named directly. A defect the default
                mode reports as forty-odd differing records shows up here as
                the one or two that are actually lost.

THE BED HAS TO PIN THE SEED. A cold start with `KICK > 0` and no `SEED` in
`plasim_namelist` draws its initial perturbation from the system clock, so two
cold starts are two different experiments and any comparison between them is
noise wearing the shape of a result. This refuses such a bed rather than
reporting the difference, because the difference would be real and would mean
nothing.

WHAT A BED IS. Any directory the model can be run in: the namelists, the
surface `.sra` files for the compiled grid, and whatever else those name. A run
directory under `exoplasim/runs/` works, and so does a purpose-built one. The
run length is taken from this script, not from the bed's `N_RUN_STEPS`.
"""

from __future__ import annotations

import argparse
import os
import shutil
import struct
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import restart_format  # noqa: E402

# The leapfrog state. Everything else in a restart is either a boundary field
# that no arithmetic touches or a diagnostic accumulator; these twelve are what
# the next timestep is computed from, so a difference here is a different
# integration rather than a different report of one.
PROGNOSTIC = ("sz", "sd", "st", "sq", "sr", "sp", "so",
              "szm", "sdm", "stm", "sqm", "spm")

# The step-count keys a bed may carry. All of them are overridden, because a
# bed that also declares months or days would run past the step count.
LENGTH_KEYS = ("N_RUN_STEPS", "N_RUN_YEARS", "N_RUN_MONTHS", "N_RUN_DAYS")


def check_seed(bed: Path) -> None:
    """Refuse a bed whose cold start is seeded from the clock."""
    nl = bed / "plasim_namelist"
    if not nl.is_file():
        raise SystemExit(f"{bed} has no plasim_namelist, so it is not a bed")
    keys = {}
    for line in nl.read_text(errors="replace").splitlines():
        if "=" in line:
            keys[line.split("=")[0].strip().upper()] = line.split("=", 1)[1].strip()
    kick = keys.get("KICK", "0").rstrip(",").strip()
    seed = keys.get("SEED")
    try:
        kicking = int(float(kick)) > 0
    except ValueError:
        kicking = True
    if kicking and (seed is None or seed.replace(",", " ").split()[:1] in ([], ["0"])):
        raise SystemExit(
            f"{nl} sets KICK = {kick} and no SEED. The initial perturbation is "
            f"then drawn from the system clock, so two cold starts are two "
            f"different experiments and this comparison would measure that "
            f"instead of the restart. Declare SEED in plasim_namelist.")


def stage(run: Path, bed: Path, steps: int, restart: Path | None) -> None:
    shutil.rmtree(run, ignore_errors=True)
    run.mkdir(parents=True)
    for f in sorted(bed.iterdir()):
        if f.is_file():
            shutil.copy2(f, run)
    nl = run / "plasim_namelist"
    lines = [l for l in nl.read_text().splitlines()
             if l.split("=")[0].strip().upper() not in LENGTH_KEYS]
    lines.insert(1, f" N_RUN_STEPS = {steps}")
    for key in LENGTH_KEYS[1:]:
        lines.insert(1, f" {key} = 0")
    nl.write_text("\n".join(lines) + "\n")
    if restart is not None:
        shutil.copy2(restart, run / "plasim_restart")


def execute(exe: Path, run: Path, threads: int) -> None:
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    with open(run / "run.log", "w") as log:
        r = subprocess.run([str(exe)], cwd=run, stdout=log,
                           stderr=subprocess.STDOUT, env=env)
    status = run / "plasim_status"
    if r.returncode != 0 or not status.is_file():
        tail = "".join((run / "run.log").read_text(
            errors="replace").splitlines(keepends=True)[-25:])
        raise SystemExit(
            f"the model failed in {run} (exit {r.returncode}); its last lines:\n{tail}")


def numbers(body: bytes):
    if len(body) % 8 == 0:
        return struct.unpack(f"<{len(body) // 8}d", body)
    return struct.unpack(f"<{len(body) // 4}i", body)


def report(a: Path, b: Path, label_a: str, label_b: str) -> int:
    A = {r.name: r.payload for r in restart_format.read(a)}
    B = {r.name: r.payload for r in restart_format.read(b)}
    only_a = [n for n in A if n not in B]
    only_b = [n for n in B if n not in A]
    common = [n for n in A if n in B]
    bad = [n for n in common if A[n] != B[n]]

    print(f"  {label_a}: {len(A)} records")
    print(f"  {label_b}: {len(B)} records")
    if only_a:
        print(f"  only in {label_a}: {' '.join(only_a)}")
    if only_b:
        print(f"  only in {label_b}: {' '.join(only_b)}")
    print(f"  differing: {len(bad)} of {len(common)}")
    for n in bad:
        x, y = numbers(A[n]), numbers(B[n])
        if len(x) != len(y):
            print(f"    {n:12s} length {len(x)} against {len(y)}")
            continue
        worst = max(abs(p - q) for p, q in zip(x, y))
        scale = max(abs(p) for p in x) or 1.0
        moved = sum(1 for p, q in zip(x, y) if p != q)
        flag = "  PROGNOSTIC" if n in PROGNOSTIC else ""
        print(f"    {n:12s} max |a-b| {worst:12.5g}  relative {worst / scale:10.4g}"
              f"  {moved}/{len(x)} elements{flag}")
    return 1 if (bad or only_a or only_b) else 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Compare a run taken whole with the same run taken in two segments.")
    ap.add_argument("--exe", required=True, help="the plasim executable to run")
    ap.add_argument("--bed", required=True,
                    help="a directory the model can run in: namelists and surface files")
    ap.add_argument("--workdir", required=True,
                    help="where the arms are staged; must not be under exoplasim/runs/")
    ap.add_argument("--steps", type=int, default=6, help="timesteps in total (default 6)")
    ap.add_argument("--split", type=int, default=3,
                    help="timesteps in the first segment (default 3)")
    ap.add_argument("--threads", type=int, default=4,
                    help="OMP_NUM_THREADS; must match the executable's compiled NPRO")
    ap.add_argument("--round-trip", action="store_true",
                    help="instead: take --steps steps, then restart and take NONE, and "
                         "name every record the restart did not carry")
    a = ap.parse_args()

    exe, bed = Path(a.exe).resolve(), Path(a.bed).resolve()
    work = Path(a.workdir).resolve()
    if not exe.is_file() or not os.access(exe, os.X_OK):
        raise SystemExit(f"{exe} is not an executable")
    check_seed(bed)
    if not 0 < a.split < a.steps:
        raise SystemExit(f"--split must be between 1 and --steps - 1, got {a.split}")

    if a.round_trip:
        first, again = work / "first", work / "again"
        stage(first, bed, a.steps, None)
        execute(exe, first, a.threads)
        stage(again, bed, 0, first / "plasim_status")
        execute(exe, again, a.threads)
        print(f"restart round trip, {a.steps} steps then a zero-step restart:")
        rc = report(first / "plasim_status", again / "plasim_status",
                    "written", "rewritten")
        print("  a difference here is state the restart does not carry: no "
              "integration happened between the two files."
              if rc else "  the restart carries the whole state.")
        return rc

    whole, seg1, seg2 = work / "whole", work / "seg1", work / "seg2"
    stage(whole, bed, a.steps, None)
    execute(exe, whole, a.threads)
    stage(seg1, bed, a.split, None)
    execute(exe, seg1, a.threads)
    stage(seg2, bed, a.steps - a.split, seg1 / "plasim_status")
    execute(exe, seg2, a.threads)

    print(f"{a.steps} steps whole against {a.split} + {a.steps - a.split}:")
    rc = report(whole / "plasim_status", seg2 / "plasim_status", "whole", "segmented")
    print("  the segmented run is a different experiment from the whole one."
          if rc else "  the segmented run is the whole run, record for record.")
    return rc


if __name__ == "__main__":
    sys.exit(main())
