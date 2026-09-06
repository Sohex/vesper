#!/usr/bin/env python3
"""Measure where `conversion_time_level` is stable, against the model itself.

    scripts/lock_and_run -m "world-bt3b arms" \\
      python exoplasim/scripts/conversion_time_arms.py \\
        --from-run exoplasim/runs/run_0d41aa82c287 --rung T21 \\
        --steps 4000 --timesteps 15,12,10,8,6 --controls 15,10

Worldbuilding frame: a numerical experiment on the Vesper project's climate
model. Nothing here is about the simulated planet.

WHAT IT IS FOR. `exoplasim/scripts/conversion_time_stability.py` derives the
timestep the modified scheme is stable at from the model's own linearised step,
and `plasim.f90`'s guard measures the same thing at startup. Both are the same
linearisation. This is the third check and the only one against the FULL
nonlinear model: arms either side of the derived boundary, run long enough for
the predicted growth to reach the range the declared `-ffpe-trap=overflow`
refuses.

TERMINATION IS READ TWO WAYS AND NEVER FROM THE EXIT STATUS ALONE.
`plasim.f90:stability_check` writes `Abort_Message` and executes a bare Fortran
`stop`, which exits 0, so an aborted run and a finished one carry the same
status; only a trap, which arrives as a signal, is distinguishable from either.
Both are recorded, separately, for every arm. `docs/src/practice/failure-modes.md`
class 38.

THE CONTROL ARMS ARE NOT OPTIONAL. An arm that dies says the configuration
died, not that the TERM killed it. `--controls` runs the same timesteps with
`nconvtime = 0`, and a control that dies means the bed or the step is what
failed and the measured arm says nothing.

THE GUARD REFUSES THE ARMS THAT WOULD SHOW THE MOST, and that is the repair
working rather than a limitation. `plasim.f90` measures its own amplification at
startup and stops before integrating a configuration that grows, so an arm above
the boundary does not run on a current binary at all -- it is recorded here as
`refused_by_guard`, which is a THIRD outcome and neither a completion nor a
death. Arms taken above the boundary before the guard existed are dated records
in `exoplasim/notes/convdecomp-reproducibility.md`, named against the binary
that ran them. What this can still measure on a current binary is the boundary
from BELOW, which is what a re-take of a comparison needs.

THE STEP AN ARM DIES AT IS BRACKETED, NOT KNOWN. The model reports at `ndiag`
intervals, so an arm that wrote k reports died after step `k * ndiag` and before
`(k+1) * ndiag`. That bracket is what turns a predicted growth rate into a
predicted amplification, and a rate is confirmed when the arms' amplification
intervals OVERLAP: one exponential, at the predicted rates, accounting for
every death. A single arm dying proves only that something grew.

WHAT IT COSTS. One bed copy per arm and one short run each; at T21 p8 with 4000
steps an arm is a few seconds. It writes no run into `exoplasim/runs/` and
registers nothing: these are beds, not runs.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _paths  # noqa: F401
from _paths import ANALYSIS, CONFIG, MODEL_RUN, PROJECT_ROOT  # noqa: E402
import rungs  # noqa: E402  the one rung-to-dimension mapping

OUT = ANALYSIS / "conversion_time_arms.json"

# The dry adiabatic set, and which namelist each key lives in.
# `run_exoplasim.py:declare_dynamics_only` is the authority for the set; NRAD is
# deliberately NOT in it, because skipping `radstep` leaves `eccf` zero and
# `outsc` then divides by its square root.
DRY = {
    "radmod_namelist": {"NSWR": 0, "NLWR": 0},
    "fluxmod_namelist": {"NVDIFF": 0, "NSHFL": 0, "NEVAP": 0, "NSTRESS": 0},
    "rainmod_namelist": {"NPRL": 0, "NPRC": 0, "NDCA": 0, "NSHALLOW": 0},
}


def set_keys(path: Path, keys: dict) -> None:
    """Set namelist keys, inserting any the file does not already carry."""
    lines = path.read_text(encoding="latin-1").splitlines()
    out, seen = [], set()
    for line in lines:
        match = re.match(r"^\s*([A-Za-z_0-9]+)\s*=", line)
        name = match.group(1).upper() if match else None
        if name in keys:
            out.append(f" {name} = {keys[name]} ")
            seen.add(name)
        else:
            out.append(line)
    missing = [k for k in keys if k not in seen]
    if missing:
        head = next((i for i, l in enumerate(out) if l.strip().startswith("&")), 0) + 1
        for key in reversed(missing):
            out.insert(head, f" {key} = {keys[key]} ")
    path.write_text("\n".join(out) + "\n", encoding="latin-1")


def run_arm(bed: Path, binary: str, threads: int) -> dict:
    """One arm, with both channels read separately."""
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_PROC_BIND="close",
               OMP_PLACES="cores", OMP_STACKSIZE="512M")
    started = time.time()
    proc = subprocess.run(["bash", "-c", f"ulimit -s unlimited; exec ./{binary}"],
                          cwd=bed, env=env, capture_output=True, text=True)
    wall = time.time() - started
    abort = bed / "Abort_Message"
    abort_step = None
    if abort.is_file():
        m = re.search(r"timestep\s*=\s*(\d+)", abort.read_text(errors="replace"))
        abort_step = int(m.group(1)) if m else -1
    diag_path = bed / "plasim_diag"
    diag = diag_path.read_text(errors="replace") if diag_path.is_file() else ""
    # THE GUARD'S REFUSAL IS ITS OWN OUTCOME. It ends with a bare Fortran
    # `stop`, so on the exit status alone it is indistinguishable from a
    # completed run -- and it wrote no timestep at all.
    refused = ("nconvtime grows at this timestep" in diag
               or "nconvtime grows at this timestep" in proc.stderr)
    return {
        "returncode": proc.returncode,
        "signalled": proc.returncode < 0,
        "wall_s": round(wall, 1),
        "abort_message": abort.is_file(),
        "abort_message_step": abort_step,
        "refused_by_guard": refused,
        # The model reports its cross sections once per `ndiag`; counting them
        # is how far the arm got, and it is a BRACKET and not a step.
        "diagnostic_reports": diag.count("Temperature [C]"),
        "guard": guard_report(diag),
    }


def guard_report(diag: str) -> dict | None:
    """What `plasim.f90`'s NCONVTIME guard printed, if it ran."""
    def last(pattern):
        found = None
        for m in re.finditer(pattern, diag, re.S):
            found = m
        return found
    on = last(r"NCONVTIME: amplification\s+([0-9.EeDd+-]+)\s+per step at total "
              r"wavenumber\s+(\d+)")
    off = last(r"NCONVTIME: control arm\s+([0-9.EeDd+-]+)\s+per step at total "
               r"wavenumber\s+(\d+)")
    if on is None or off is None:
        return None
    real = lambda x: float(x.replace("D", "E").replace("d", "e"))
    return {"measured_arm": real(on.group(1)),
            "measured_arm_wavenumber": int(on.group(2)),
            "control_arm": real(off.group(1)),
            "control_arm_wavenumber": int(off.group(2))}


def predicted(sweep: Path, dt: float) -> float | None:
    """The growth per step this timestep is predicted to have, if recorded."""
    if not sweep.is_file():
        return None
    for row in json.loads(sweep.read_text()).get("rows", []):
        if abs(float(row.get("timestep_minutes", -1)) - dt) < 1e-9:
            return row.get("growth_per_step")
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--from-run", type=Path, required=True,
                        help="the run whose restart and inputs the bed is cut "
                             "from; its restart must carry the record set the "
                             "current binary reads")
    parser.add_argument("--rung", default="T21")
    parser.add_argument("--ranks", type=int, default=8)
    parser.add_argument("--steps", type=int, default=4000)
    parser.add_argument("--ndiag", type=int, default=200)
    parser.add_argument("--timesteps", default="15,12,10,8,6",
                        help="minutes, run with the term ON")
    parser.add_argument("--controls", default="15,10",
                        help="minutes, run with the term OFF; a control that "
                             "dies means the bed failed and the measured arm "
                             "says nothing")
    parser.add_argument("--work", type=Path,
                        default=Path("/tmp/vesper-conversion-time-arms"))
    parser.add_argument("--sweep", type=Path,
                        default=ANALYSIS / "conversion_time_stability_sweep.json",
                        help="the derived growth rates each arm is judged "
                             "against")
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()

    rung = args.rung.upper()
    rungs.geometry(rung)          # refuses a rung that is not on the ladder
    binary = f"most_plasim_{rung.lower()}_l10_p{args.ranks}.x"
    on = [float(x) for x in args.timesteps.replace(",", " ").split()]
    off = [float(x) for x in args.controls.replace(",", " ").split()]
    if not off:
        raise SystemExit(
            "--controls is empty. An arm that dies says the configuration "
            "died, not that the term killed it; without a control arm at the "
            "same step this measures nothing.")

    args.work.mkdir(parents=True, exist_ok=True)
    bed0 = args.work / "bed"
    if bed0.exists():
        shutil.rmtree(bed0)
    make_bed = subprocess.run(
        [sys.executable, str(Path(__file__).parent / "make_profile_bed.py"),
         "--from-run", str(args.from_run.resolve()), "--dest", str(bed0),
         "--binary", binary, "--binary-dir", str(MODEL_RUN),
         "--steps", str(args.steps)],
        capture_output=True, text=True)
    print(make_bed.stdout, make_bed.stderr)
    if make_bed.returncode != 0:
        return 1

    arms = []
    for dt, nconv in [(x, 1) for x in on] + [(x, 0) for x in off]:
        bed = args.work / f"dt{dt:g}_term{nconv}"
        if bed.exists():
            shutil.rmtree(bed)
        shutil.copytree(bed0, bed, symlinks=True)
        set_keys(bed / "plasim_namelist",
                 {"MPSTEP": f"{dt}", "N_RUN_STEPS": args.steps,
                  "NDIAG": args.ndiag, "NENERGY": 1, "NENERGYFIX": 0,
                  "NCONVTIME": nconv})
        for name, keys in DRY.items():
            set_keys(bed / name, keys)
        row = {"timestep_minutes": dt, "conversion_time_level": bool(nconv),
               "steps_requested": args.steps, "ndiag": args.ndiag}
        row.update(run_arm(bed, binary, args.ranks))
        refused = row["refused_by_guard"]
        finished = (not refused and row["returncode"] == 0
                    and not row["abort_message"])
        reports = row["diagnostic_reports"]
        row["completed"] = finished
        died = not finished and not refused
        row["died_after_step"] = reports * args.ndiag if died else None
        row["died_before_step"] = (reports + 1) * args.ndiag if died else None
        rate = predicted(args.sweep, dt)
        row["predicted_growth_per_step"] = rate
        # An arm the guard REFUSED never integrated, so it has no death to
        # bracket and no amplification to have reached. `not finished` covers
        # both refusals and deaths and is not the test.
        if rate and rate > 1.0 and row["died_after_step"] is not None:
            ln = math.log(rate)
            row["amplification_reached"] = [
                float(f"{math.exp(ln * row['died_after_step']):.3g}"),
                float(f"{math.exp(ln * row['died_before_step']):.3g}")]
        arms.append(row)
        print(json.dumps(row), flush=True)

    # THE RATE IS CONFIRMED BY AN OVERLAP AND NOT BY A SINGLE ARM. Each death is
    # bracketed by the report cadence, and the predicted rate turns that bracket
    # into the amplification the arm reached. If one exponential at the
    # predicted rates accounts for every death, those intervals share a value.
    reached = [a["amplification_reached"] for a in arms
               if a.get("amplification_reached")]
    overlap, agreement = None, None
    if len(reached) < 2:
        # NOT THE SAME AS A DISAGREEMENT, and a null that meant both would be
        # unreadable. Fewer than two deaths is nothing to compare, which is the
        # ordinary case once the guard refuses the arms that would have died.
        agreement = (f"{len(reached)} arm(s) died, so there is nothing to "
                     "compare. Two deaths are the fewest that can share a rate.")
    else:
        lo = max(x[0] for x in reached)
        hi = min(x[1] for x in reached)
        if lo <= hi:
            overlap = [float(f"{lo:.3g}"), float(f"{hi:.3g}")]
            agreement = ("one exponential, at the predicted rates, accounts "
                         "for every death")
        else:
            agreement = ("the deaths do NOT share a rate: no single "
                         "amplification is consistent with all their brackets, "
                         "so the prediction is wrong about the size of the "
                         "effect even where it is right about its sign")
    controls_held = all(a["completed"] for a in arms
                        if not a["conversion_time_level"])
    refused = [a["timestep_minutes"] for a in arms if a["refused_by_guard"]]

    result = {
        "what": ("where model.conversion_time_level is stable, measured "
                 "against the full model rather than its linearisation"),
        "bed": {"rung": rung, "ranks": args.ranks, "binary": binary,
                "parent_run": str(args.from_run),
                "steps_requested": args.steps, "ndiag": args.ndiag,
                "dry_adiabatic": True, "nenergy": 1, "energy_fixer": False},
        "how_termination_is_read": (
            "Abort_Message and the exit status separately, never the exit "
            "status alone: plasim.f90:stability_check writes that file and "
            "executes a bare Fortran stop, which exits 0. A negative "
            "returncode is a signal, which under the declared "
            "-ffpe-trap=invalid,zero,overflow is the trap."),
        "controls_held": controls_held,
        "refused_by_guard_minutes": refused,
        "refused_by_guard_note": (
            "plasim.f90 measures its own amplification at startup and stops "
            "before integrating a configuration that grows, so an arm above the "
            "boundary does not run on a current binary. That is a third "
            "outcome and neither a completion nor a death; the guard's stop "
            "exits 0, so the exit status alone cannot tell it from a "
            "completed run."),
        "shared_amplification": overlap,
        "shared_amplification_verdict": agreement,
        "shared_amplification_note": (
            "The interval every dying arm's bracket is consistent with, and "
            "the verdict beside it says which of the three things a null "
            "means: too few deaths to compare, or deaths that share no rate. "
            "A non-null interval says one exponential, at the rates "
            "conversion_time_stability.py predicts, accounts for every death."),
        "provenance": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "sweep": str(args.sweep.relative_to(PROJECT_ROOT))
                     if args.sweep.is_relative_to(PROJECT_ROOT) else str(args.sweep),
            "config": str(CONFIG.relative_to(PROJECT_ROOT)),
        },
        "arms": arms,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")
    if not controls_held:
        print("A CONTROL ARM DID NOT COMPLETE, so the bed or the timestep is "
              "what failed and the measured arms say nothing about the term.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
