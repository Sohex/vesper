#!/usr/bin/env python3
"""Is the climate model run-to-run reproducible, and what decides it?

    python exoplasim/scripts/reproducibility_matrix.py <bed_dir> [--repeats 3]
    python exoplasim/scripts/reproducibility_matrix.py --plan

Worldbuilding. Vesper is an invented planet and this script is about the
simulation of it: whether the climate model gives the same bytes twice from the
same inputs. Nothing here is a claim about the world.

CLIM-44. This exists because two measurements of this project's own model
disagree, and every A/B this project runs depends on which of them is general:

- CLIM-39, at T42 L10 on **8 ranks** with output ON, found gridpoint output
  bit-identical at 1 timestep and DIFFERENT by 16, and found `plasim_status`
  different between two runs of one binary even at 1 timestep.
- The build-flag benchmark (`notes/audits/aocl-and-model-build-flags.md`), at
  T42 L10 on **16 ranks** with `NOUTPUT = 0`, found `plasim_status`
  BIT-IDENTICAL across about 48 repeats spanning four distinct binaries.

So the model is not non-reproducible in general, and the question is what makes
the difference. Two candidates were already named in the task row -- rank count
and whether the segment writes output -- and reading the source before running
anything added a third that outranks both.

## The third candidate, and why it is first

**The CLIM-39 measurement was taken on a tree that does not contain the
`zsolars` fix.** `git merge-base --is-ancestor 80f11e9 b43f17e` is false: the
worktree branched before the toolchain commit that fixed it.

Before that fix, `radstop` wrote the two solar constants through `mpputgp`,
which gathers `NHOR` elements per rank from an array holding 2 and writes a
record of `NUGP` doubles -- 8190 of them memory past the end of `zsolars(2)`.
`notes/audits/zsolars-restart-overread.md` has the measurement, including the
part that matters here: those bytes are stable only because `-finit-real=zero`
holds them, and the `-flto` arm returned a different restart sha on every one
of five runs from one binary on one input.

That is a complete and sufficient explanation for the `plasim_status` half of
CLIM-39's observation, and it is ALREADY FIXED. It does not explain the
gridpoint-output half, because nothing reads the record back and the audit
establishes the integration is untouched by it.

## What this measures

A full factorial, because the three candidates are not nested and a one-factor
sweep would confound whichever two it holds fixed:

| factor | levels |
| --- | --- |
| ranks | 8, 16 |
| output | `NOUTPUT`/`NSNAPSHOT` 0, or the production 1 |
| steps | derived from the write cadence, see below |

**The segment lengths are derived rather than chosen, and that is a correction
to CLIM-39's control.** A gridpoint record is written when
`mod(nstep, nafter) == 0` on the ABSOLUTE step count, so whether a short segment
writes anything is a property of the restart. On this project's own production
restart `nstep` is 589672 and the model logs `nafter` as 32, leaving 24 steps to
the next write -- so NEITHER of CLIM-39's 1-step and 16-step segments writes a
gridpoint record at all. An output comparison over files that hold no model
records cannot fail, and a cell that cannot fail is not evidence.

So one length is set to stop one step short of the first write and the other to
cross two. The second also straddles the FIRST record after a resume, which is
the record `exoplasim/notes/first-output-bin.md`'s class corrupts: `naccuout`
survives the restart holding a partial window -- 26 here -- so that record is
divided by 50 rather than by 32. That defect is DETERMINISTIC and cannot produce
a run-to-run difference, but it is the one thing that makes the first record
after a resume unlike the records after it, so the matrix covers it explicitly
rather than by accident.

Every cell is run `--repeats` times from the SAME restart, and the comparison
is between repeats of one cell -- one binary, one input, one configuration. A
cell is REPRODUCIBLE if every repeat agrees bit for bit on every artifact it
produced, and NOT otherwise. There is no tolerance and no threshold to choose:
bytes are equal or they are not.

## Declared before running, so the reading cannot be fitted

Each hypothesis owns a pattern, and the matrix can refute all three:

- **H1, writing output is what does it.** Cells at `NOUTPUT = 0` reproduce and
  cells at `NOUTPUT = 1` do not, at BOTH rank counts.
- **H2, MPI reduction ordering at 8 ranks.** Cells at 8 ranks do not
  reproduce and cells at 16 ranks do, at both output settings.
- **H3, the model reproduces now.** Every cell reproduces, CLIM-39's own
  configuration included. Note what this does NOT establish: H3 is a pattern,
  not a mechanism. Running the same matrix on the pre-fix binaries a run
  directory keeps -- which still write `zsolars` as a 65536-byte record --
  gives the SAME pattern, so "it was the overread" is refuted as the cause of
  run-to-run variation even though the overread is real. It is deterministic
  per rank count under `-finit-real=zero`, which makes it a source of
  RANK-dependence rather than of run-to-run difference.

H3 is the one to try hardest to refute, because it is the comfortable answer
and because it is the only one that would let this project go back to
expecting bit-identity. If H3 holds, the standing instruction that every A/B
must establish its own reproducibility horizon is still right -- it just costs
one control run rather than a bound on every claim.

A cell that reproduces at the short length and not at the long one refutes none
of the three on its own: it separates the model's INTEGRATION from its
initialisation, and it is reported as its own column rather than folded into the
verdict.

A cell in the output arm that writes no gridpoint record is reported
INCONCLUSIVE and is excluded from the scoring, rather than counted as
reproducible. That is the failure mode this design exists to avoid.

## Building a bed

There is no step that makes one. Copy a run directory, drop its outputs, and
REPLACE ITS EXECUTABLES:

    cp -a exoplasim/runs/<run> bed && rm -f bed/MOST* bed/plasim_output \
        bed/plasim_snapshot bed/plasim_status && rm -rf bed/snapshots/*
    cp -f vendor/exoplasim/exoplasim/plasim/run/most_plasim_t42_l10_p*.x bed/

The second line is not optional and `check_binaries` refuses without it: a run
directory keeps the executables it was STARTED with, so a bed copied from a
recent run silently pins the model source to whenever that run began. The
outputs are dropped because `snapshots/` alone can be gigabytes and none of it
is read.

## Instrumenting the model, if a cell needs opening up

Two things cost a cycle each here and are not obvious.

**Only NROOT's writes to `nud` survive.** `nud` is unit 6 and every rank writes
to it, but the diagnostics that reach `plasim_diag` are the root's; the others
go nowhere. So `if (mypid == NROOT)` is not a filter, it is the only thing that
works -- and at T42 on 8 ranks NROOT holds the first eight latitude rows, which
are all polar. A diagnostic written that way is sampling one climate zone. To
see any other cell, write to a per-rank unit: `write(70+mypid,...)` gives one
`fort.7N` per rank.

**A temporary namelist key can break the run for reasons that are not about
the key.** Adding one to `radmod_nl` and setting it in `radmod_namelist` failed
with "Cannot match namelist object name" even with the key compiled in and
present. Rather than chase it, read the switch from the ENVIRONMENT --
`get_environment_variable` -- which touches no namelist group and no shared
file. That is what CLIM-42's per-rank probe ended up doing.

## What it does not measure

Whether a run that diverges diverges CHAOTICALLY. CLIM-39 established that
separately -- 2.4e-6 relative in the first differing record growing to order 1 --
and a byte comparison cannot see it. If a cell fails here, that is the next
measurement and not this one.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import ANALYSIS, PROJECT_ROOT  # noqa: E402

OUTPUT = ANALYSIS / "reproducibility_matrix.json"

RANKS = (8, 16)
OUTPUT_LEVELS = (0, 1)

# Artifacts compared between repeats. Missing ones are recorded as missing
# rather than skipped, so a cell that silently stopped producing output cannot
# pass by having nothing to disagree about.
ARTIFACTS = ("plasim_status", "plasim_output", "plasim_snapshot", "plasim_diag")

# `plasim_diag` is an ASCII log that ends with the run's own resource
# accounting, so two identical integrations differ in it by construction. These
# lines are removed before it is hashed. The file is NOT dropped: everything
# above them is model diagnostics -- the initial step number, the mean surface
# pressure, the spin-up trace -- and that is worth comparing.
#
# The list is taken from the block that writes it, `plasim.f90:1011-1042`,
# rather than from what one pair of runs happened to differ on. EVERY line
# there is conditional on its counter being non-zero, so a filter built by
# reading a diff is complete only for the runs that produced it: "Page faults"
# appears in some repeats and not others, which is exactly how a
# read-the-diff filter passes three cells and fails three more.
WALL_CLOCK = re.compile(
    r"^\* (User   time|System time|Total CPU time|Memory usage"
    r"|Page reclaims|Page faults|Page swaps|Disk read|Disk write"
    r"|Seconds per sim year|Minutes per sim year|Days per sim year"
    r"|Sim years per day)")

HYPOTHESES = {
    "H1_output_writing": "NOUTPUT = 0 reproduces, NOUTPUT = 1 does not, at both rank counts",
    "H2_rank_reduction_order": "8 ranks does not reproduce, 16 ranks does, at both output settings",
    "H3_reproducible_now": "every cell reproduces, CLIM-39's own configuration included",
}


def restart_scalars(path: Path) -> dict:
    """`nstep` and `naccuout` out of a restart, by walking its record pairs.

    The write cadence is `mod(nstep, nafter) == 0` on the ABSOLUTE step count,
    so where a short segment falls against it is a property of the restart and
    not of the segment length. That is why this is read rather than assumed.
    """
    import struct
    data = path.read_bytes()
    out, i = {}, 0
    while i < len(data) - 8:
        (n,) = struct.unpack_from("<i", data, i)
        if n <= 0 or i + 8 + n > len(data):
            break
        payload = data[i + 4:i + 4 + n]
        if struct.unpack_from("<i", data, i + 4 + n)[0] != n:
            break
        i += 8 + n
        if n != 16:
            continue
        name = payload.decode("latin-1").strip()
        (m,) = struct.unpack_from("<i", data, i)
        if name in ("nstep", "naccuout") and m == 4:
            out[name] = struct.unpack_from("<i", data, i + 4)[0]
        i += 8 + m
    for key in ("nstep", "naccuout"):
        if key not in out:
            raise SystemExit(f"{path.name} carries no {key} record")
    return out


def segment_lengths(nstep: int, nafter: int) -> dict:
    """Two segment lengths chosen against the write cadence, not round numbers.

    CLIM-39's control was run at 1 and 16 timesteps. On this project's own
    production restart `nstep` is 589672 and `nafter` is 32, so the next write
    is 24 steps away and NEITHER of those segments writes a gridpoint record at
    all. A comparison of output files that contain no model records cannot fail,
    and a cell that cannot fail is not evidence.

    So the lengths are derived: one segment that deliberately crosses no write,
    and one that crosses two. The second also straddles the FIRST record after a
    resume, which is the record the accumulator class corrupts -- `naccuout`
    survives the restart holding a partial window, so that record is divided by
    a count that is neither `nafter` nor the number of samples in it.
    `exoplasim/notes/first-output-bin.md` is the argument.
    """
    to_first = nafter - (nstep % nafter)
    return {"no_write": max(1, to_first - 1), "two_writes": to_first + nafter,
            "steps_to_first_write": to_first}


def sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    if path.name == "plasim_diag":
        for line in path.read_text(encoding="latin-1").splitlines(keepends=True):
            if not WALL_CLOCK.search(line):
                h.update(line.encode("latin-1"))
        return h.hexdigest()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def set_namelist(path: Path, **keys: int) -> None:
    """Rewrite `plasim_namelist` keys in place, preserving everything else."""
    text = path.read_text(encoding="utf-8")
    for key, value in keys.items():
        pattern = re.compile(rf"^(\s*{key}\s*=\s*)\S+\s*$", re.M | re.I)
        if not pattern.search(text):
            raise SystemExit(f"{key} is not in {path.name}; refusing to add it, "
                             "because a key this script invents is a key the bed "
                             "was not measured with")
        text = pattern.sub(rf"\g<1>{value} ", text)
    path.write_text(text, encoding="utf-8")


def check_binaries(bed: Path) -> dict:
    """Refuse a bed whose executables are not the ones the manifest describes.

    A run directory keeps the binaries it was STARTED with, so copying a bed out
    of `exoplasim/runs/` silently pins the model source to whenever that run
    began. That is how the first pass of this matrix came to measure a tree
    three commits behind the one `rebuild_binaries.py --verify` calls current,
    and nothing about the result said so: it ran, it was self-consistent, and it
    was about the wrong code.

    `--verify` answers "are the INSTALLED binaries current". It cannot answer
    "is the bed running them", and that is the question here.
    """
    manifest = json.loads((PROJECT_ROOT / "exoplasim" / "binary_manifest.json")
                          .read_text(encoding="utf-8"))["binaries"]
    seen, stale = {}, []
    for ranks in RANKS:
        name = f"most_plasim_t42_l10_p{ranks}.x"
        path = bed / name
        if not path.exists():
            raise SystemExit(f"{name} is not in the bed")
        got = hashlib.sha256(path.read_bytes()).hexdigest()
        seen[name] = got
        if name not in manifest:
            stale.append(f"{name} is in no manifest entry")
        elif got != manifest[name]["sha256"]:
            stale.append(f"{name} is {got[:12]}, manifest says "
                         f"{manifest[name]['sha256'][:12]}")
    if stale:
        raise SystemExit(
            "the bed's executables are not the ones binary_manifest.json "
            "describes, so this would measure a model source nobody named:\n  "
            + "\n  ".join(stale)
            + "\nCopy them from vendor/exoplasim/exoplasim/plasim/run/, or "
              "rebuild. CLAUDE.md rule 4.")
    return seen


def run_cell(bed: Path, work: Path, ranks: int, noutput: int, steps: int) -> dict:
    if work.exists():
        shutil.rmtree(work)
    shutil.copytree(bed, work, symlinks=True)

    set_namelist(work / "plasim_namelist",
                 N_RUN_STEPS=steps, NOUTPUT=noutput, NSNAPSHOT=noutput)

    binary = f"most_plasim_t42_l10_p{ranks}.x"
    if not (work / binary).exists():
        raise SystemExit(f"{binary} is not in the bed. ExoPlaSim compiles one "
                         "executable per (resolution, layers, ranks) triple; "
                         "CLAUDE.md rule 4.")

    started = datetime.datetime.now(datetime.timezone.utc)
    proc = subprocess.run(["mpiexec", "-np", str(ranks), f"./{binary}"],
                          cwd=work, capture_output=True, text=True)
    elapsed = (datetime.datetime.now(datetime.timezone.utc) - started).total_seconds()

    if (work / "Abort_Message").exists() or proc.returncode != 0:
        raise SystemExit(f"the model aborted at ranks={ranks} noutput={noutput} "
                         f"steps={steps} (rc {proc.returncode})\n"
                         f"{proc.stdout[-2000:]}\n{proc.stderr[-2000:]}")

    return {"seconds": round(elapsed, 2),
            "sha256": {name: sha256(work / name) for name in ARTIFACTS},
            "bytes": {name: ((work / name).stat().st_size
                             if (work / name).exists() else None)
                      for name in ARTIFACTS}}


def score(cells: dict) -> dict:
    """Which hypothesis survives. Patterns are the ones declared above."""
    vacuous = sorted(k for k, c in cells.items() if c["vacuous"])
    repro = {key: cell["reproducible"] for key, cell in cells.items()
             if not cell["vacuous"]}
    if not repro:
        return {"per_hypothesis": {}, "survivors": [], "vacuous_cells": vacuous,
                "reading": "every cell was vacuous; nothing was measured"}

    def all_of(pred, want):
        return all(v is want for k, v in repro.items() if pred(k))

    def part(key, i):
        return int(key.split("_")[i])

    verdict = {
        "H1_output_writing": (all_of(lambda k: part(k, 1) == 0, True)
                              and all_of(lambda k: part(k, 1) == 1, False)),
        "H2_rank_reduction_order": (all_of(lambda k: part(k, 0) == 8, False)
                                    and all_of(lambda k: part(k, 0) == 16, True)),
        "H3_reproducible_now": all(repro.values()),
    }
    survivors = [k for k, v in verdict.items() if v]
    # Reproducible WITHIN a rank count is not the same as rank-independent, and
    # the matrix already holds the evidence for both. This is reported rather
    # than scored: no hypothesis above is about it, and it was not predicted.
    across = {}
    for noutput in OUTPUT_LEVELS:
        for steps in {c["steps"] for c in cells.values()}:
            keys = [k for k, c in cells.items()
                    if c["noutput"] == noutput and c["steps"] == steps]
            shas = {cells[k]["repeats"][0]["sha256"]["plasim_status"] for k in keys}
            if len(keys) > 1:
                across[f"noutput={noutput},steps={steps}"] = (
                    "same across rank counts" if len(shas) == 1
                    else "DIFFERS across rank counts")

    return {"per_hypothesis": verdict,
            "survivors": survivors,
            "vacuous_cells": vacuous,
            "across_rank_counts": across,
            "reading": (survivors[0] if len(survivors) == 1 else
                        "NONE of the declared hypotheses matches the pattern; "
                        "the cause is something the matrix did not vary")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("bed", nargs="?", type=Path,
                    help="a run directory holding the namelists, the .sra "
                         "inputs, plasim_restart and the T42 executables")
    ap.add_argument("--repeats", type=int, default=3,
                    help="repeats per cell; two is the minimum that can "
                         "disagree and three separates one odd run from a "
                         "cell that never repeats")
    ap.add_argument("--work", type=Path, default=None,
                    help="scratch directory for the run copies")
    ap.add_argument("--nafter", type=int, default=32,
                    help="timesteps between writes. The model prints it as "
                         "\"Timesteps / write\" in MOST_DIAG; do not guess it, "
                         "because the whole point of the segment lengths is "
                         "where they fall against it")
    ap.add_argument("--plan", action="store_true",
                    help="print the matrix and the declared readings, run nothing")
    args = ap.parse_args()

    if args.bed is not None:
        scalars = restart_scalars(args.bed.resolve() / "plasim_restart")
        lengths = segment_lengths(scalars["nstep"], args.nafter)
    else:
        scalars = {"nstep": 589672, "naccuout": 26}
        lengths = segment_lengths(scalars["nstep"], args.nafter)
    steps = (lengths["no_write"], lengths["two_writes"])

    cells = list(product(RANKS, OUTPUT_LEVELS, steps))
    if args.plan:
        print(f"{len(cells)} cells x {args.repeats} repeats "
              f"= {len(cells) * args.repeats} model runs\n")
        print(f"  nstep {scalars['nstep']}, naccuout {scalars['naccuout']}, "
              f"nafter {args.nafter}")
        print(f"  first write {lengths['steps_to_first_write']} steps in, so "
              f"{steps[0]} crosses none and {steps[1]} crosses two\n")
        for ranks, noutput, steps in cells:
            print(f"  ranks={ranks:2d}  NOUTPUT={noutput}  steps={steps:2d}")
        print("\ndeclared readings:")
        for name, pattern in HYPOTHESES.items():
            print(f"  {name}: {pattern}")
        return

    if args.bed is None:
        ap.error("a bed directory is required unless --plan is given")
    bed = args.bed.resolve()
    if not (bed / "plasim_restart").exists():
        raise SystemExit(f"{bed} has no plasim_restart; this measures a RESUME, "
                         "because a cold start seeds its RNG from the clock "
                         "(plasim.f90 initrandom) and would be non-reproducible "
                         "by construction rather than by defect")

    binaries = check_binaries(bed)
    work_root = (args.work or bed.parent / "work").resolve()
    work_root.mkdir(parents=True, exist_ok=True)

    results = {}
    for ranks, noutput, steps in cells:
        key = f"{ranks}_{noutput}_{steps}"
        repeats = []
        for i in range(args.repeats):
            print(f"  ranks={ranks:2d} NOUTPUT={noutput} steps={steps:2d} "
                  f"repeat {i + 1}/{args.repeats}", flush=True)
            repeats.append(run_cell(bed, work_root / f"{key}_{i}", ranks, noutput, steps))

        first = repeats[0]["sha256"]
        disagreeing = sorted(name for name in ARTIFACTS
                             if any(r["sha256"][name] != first[name] for r in repeats))
        produced = sorted(name for name in ARTIFACTS if first[name] is not None)

        # An output arm that wrote no gridpoint record cannot disagree, and a
        # cell that cannot fail is not evidence. CLIM-39's control was run at 1
        # and 16 steps against nafter = 32 with nstep leaving 24 steps to the
        # next write, so this is the trap that reading recorded as a result.
        # One gridpoint record at T42 is 8192 four-byte values plus its name
        # record, about 32 kB. A file at that size holds a single constant
        # field and no timestep write, which is the near-vacuous case the
        # segment lengths exist to straddle: it is reported, not hidden.
        out_bytes = repeats[0]["bytes"]["plasim_output"] or 0
        wrote_output = out_bytes > 2 * 32816
        vacuous = noutput == 1 and not wrote_output
        results[key] = {
            "ranks": ranks, "noutput": noutput, "steps": steps,
            "repeats": repeats,
            "artifacts_produced": produced,
            "artifacts_disagreeing": disagreeing,
            "output_bytes": out_bytes,
            "reproducible": None if vacuous else not disagreeing,
            "vacuous": vacuous,
        }
        if vacuous:
            verdict = (f"INCONCLUSIVE: output on, but plasim_output is "
                       f"{out_bytes} bytes -- under two gridpoint records")
        elif disagreeing:
            verdict = "DIFFERS on " + ", ".join(disagreeing)
        else:
            verdict = "reproducible"
        print(f"    -> {verdict}  (produced {', '.join(produced) or 'nothing'})",
              flush=True)

    report = {
        "generated": datetime.datetime.now(datetime.timezone.utc)
                             .replace(microsecond=0).isoformat(),
        "generator": "exoplasim/scripts/reproducibility_matrix.py",
        "task": "CLIM-44",
        "bed": str(bed),
        "repeats": args.repeats,
        "host_cpu_count": os.cpu_count(),
        "binaries": binaries,
        "model_source": subprocess.run(
            ["git", "log", "-1", "--format=%H %s", "--", "vendor/exoplasim"],
            cwd=PROJECT_ROOT, capture_output=True, text=True).stdout.strip(),
        "hypotheses": HYPOTHESES,
        "cells": results,
        "verdict": score(results),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"\nverdict: {report['verdict']['reading']}")
    print(f"wrote {OUTPUT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
