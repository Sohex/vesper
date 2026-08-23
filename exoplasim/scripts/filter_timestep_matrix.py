#!/usr/bin/env python3
"""The filter-strength x timestep x resolution matrix, on a wall-clock budget.

    python exoplasim/scripts/filter_timestep_matrix.py --hours 6.5

Worldbuilding frame: a numerical experiment on the Vesper project's climate
model. Nothing here is about the simulated planet.

WHAT IT IS FOR. CLIM's physics-filter row and SPAT-11 turned out to be one
question. `exoplasim/notes/physics-filter-stability.md` measures the filter
buying timestep -- the model traps in the shortwave at a 45-minute step without
it and runs at 22.5 -- so the largest acceptable step is a function of filter
strength, and a ladder that fixes one and sweeps the other measures one
arbitrary point on a trade. This sweeps both.

TWO TRACKS, because they answer to different states.

  A. THE ENERGY IDENTITY, at T42 only. `denergy26 - denergy27` is read from a
     SPUN-UP state, and the only restart this project has is T42's. Cold-started
     arms would give the identity on a state that has not settled, which is not
     comparable. So track A fills the kappa x dt grid from
     `run_4182235e9781`'s restart and extends the five arms already measured.

  B. STABILITY AND COST, every rung, cold. Whether a configuration traps, and
     what it costs per orbit, are properties of the resolution and the namelist
     rather than of a settled state -- `make_profile_bed.py` makes the same
     argument for profiling. So track B cold-starts, which is the only thing
     available above T42 and is adequate for both questions. Its cost column IS
     SPAT-11's deliverable.

THE BUDGET IS A DEADLINE AND THE ESTIMATES ARE MEASURED, NOT GUESSED. Every job
is priced before it starts and skipped if it will not fit in what is left. The
price starts from this project's own numbers -- 87 to 94 s per orbit at T42 dt
45, 187 at dt 22.5, so cost is linear in step count -- and each completed job
RE-MEASURES its rung's factor, so the schedule sharpens as it runs. The
resolution factors below are DECLARED BRACKETS from SPAT-11's recorded work
ratios, not measurements, and they are labelled as such in the output until a
job replaces them.

NO SILENT CAPS. Every job the deadline refuses is written to the result with
its estimate and the reason, because a matrix that quietly stopped early reads
as a matrix that covered everything.

FAIL-SOFT. Each job is independently tried; a failure is recorded with its
diagnosis -- trap, crash, or refusal -- and the sweep continues. The result file
is rewritten after every job, so killing this at any point leaves usable output.
"""
from __future__ import annotations

import argparse
import copy
import json
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _paths  # noqa: F401

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config" / "planet.yaml"
RUNS = ROOT / "exoplasim" / "runs"
OUT = ROOT / "exoplasim" / "analysis" / "filter_timestep_matrix.json"
# The sprint writes BESIDE the matrix, never over it. They are different
# experiments and the matrix took five and a half hours; an output path shared
# between a long run and a focused follow-up is one --sprint away from deleting
# the thing the follow-up exists to extend.
OUT_SPRINT = ROOT / "exoplasim" / "analysis" / "filter_timestep_sprint.json"
WORK = Path("/tmp/vesper-filter-matrix")
PARENT_RESTART = RUNS / "run_4182235e9781" / "MOST_REST.00039"

# MEASURED on this host, 2026-08-23: T42 dt 45 gives 87 and 94 s per simulated
# orbit on two arms, dt 22.5 gives 187 on two more. Twice the steps, twice the
# time, so the model cost is linear in step count and the anchor is per-step.
T42_SECONDS_PER_ORBIT_AT_DT45 = 90.0

# DECLARED BRACKETS, not measurements. From SPAT-11's recorded ratios: Legendre
# work 7.9x at T85 and 26.2x at T127 against T42, grid work 4.0x and 9.0x, and
# transform shares of 5.5, 10.4 and 14.9 percent across the rungs. Composed as
# share-weighted work; T170 extrapolates grid work as nlat^2 and Legendre as
# ntru^3. Each is REPLACED by measurement the first time a job at that rung
# completes, and the output says which of the two any number is.
# MEASURED on this host by the 2026-08-23 sweep, normalised to dt 45: T21 37,
# T42 93, T85 296 and T127 760 seconds per orbit. The brackets these replace
# said 0.30, 4.2 and 10.0, so they ran 25 to 30 percent high at the top --
# worth knowing before anything budgets off a bracket. T170 is STILL A BRACKET:
# every T170 arm in that sweep trapped, so nothing has priced it.
RUNG_FACTOR = {"T21": 0.407, "T42": 1.00, "T85": 3.20, "T127": 8.21, "T170": 18.7}

# `model.resolution` does NOT carry the grid on its own. `read_sra` validates
# every staged surface file against `model.latitudes` and `model.longitudes`,
# so a config with the resolution changed and those left behind refuses its own
# inputs with an "Unexpected SRA header". That is SPAT-2's finding -- the
# artifact path still carries resolution literals -- met head on, and it is why
# these three move together here.
RUNG_GRID = {"T21": (32, 64), "T42": (64, 128), "T85": (128, 256),
             "T127": (192, 384), "T170": (256, 512)}

# Prep, staging and postprocessing, which the model's own seconds-per-orbit does
# not count. Scaled by rung because the postprocessor reads the raw file.
OVERHEAD_S = 150.0

TRAP = re.compile(r"SIGFPE|Floating-point exception")
# A trap at the first radiation call and a blow-up twenty minutes in are both
# SIGFPE and are not the same event. The wall clock separates them, and the
# outcome records which: a refusal is a property of the configuration, a late
# failure is the integration going unstable and has to be read as such.
LATE_FAILURE_S = 120.0
CRASH = re.compile(r"crashed or begun producing garbage")


def now() -> float:
    return time.monotonic()


def estimate(rung: str, dt: float, orbits: int) -> float:
    model = RUNG_FACTOR[rung] * T42_SECONDS_PER_ORBIT_AT_DT45 * (45.0 / dt) * orbits
    return model + OVERHEAD_S * (1.0 + RUNG_FACTOR[rung] * 0.5)


def write_config(base: dict, name: str, rung: str, kappa: float | None,
                 dt: float, energy: bool, bootstrap: bool = False) -> Path:
    """One arm's config, derived from planet.yaml rather than hand-copied.

    `bootstrap` is what makes the ladder reachable, and it is the project's own
    vocabulary rather than a trick: a BOOTSTRAP RUN is the first climate run on
    a build, on TERRAIN-ONLY fields. `intended_surface_codes` requires only
    BASE_SURFACE_CODES -- orography 129 and the land mask 172 -- when the three
    `*_source` keys are `uniform`, and every rung on the ladder has exactly
    those two staged. The derived fields that only T42 carries are what a
    BASELINE run needs, and this sweep is not measuring a baseline: it is
    measuring where the model traps and what an orbit costs, both of which are
    properties of the resolution and the namelist.

    The sharp orography that the filter exists to tame is code 129, so it is
    present in this configuration. That is what makes the trap question askable
    here at all.
    """
    cfg = copy.deepcopy(base)
    m = cfg["model"]
    if bootstrap:
        m["land_albedo_source"] = "uniform"
        m["soil_water_source"] = "uniform"
        m["roughness_source"] = "uniform"
        m["dust_source"] = "none"
        m["dust_emission"] = "none"
    m["resolution"] = rung
    m["latitudes"], m["longitudes"] = RUNG_GRID[rung]
    m["timestep_minutes"] = dt
    m["energy_diagnostics"] = energy
    m["energy_diagnostics_3d"] = False
    if kappa is None:
        m["physics_filter"] = ""
    else:
        m["physics_filter"] = "gp|exp|sp"
        m["filter_kappa"] = float(kappa)
    WORK.mkdir(parents=True, exist_ok=True)
    path = WORK / f"planet_{name}.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return path


def seconds_per_orbit(run_dir: Path) -> float | None:
    """The model's own timing line, which counts integration and not staging."""
    for diag in sorted(run_dir.glob("MOST_DIAG.*"), reverse=True):
        try:
            text = diag.read_bytes()[-4000:].decode("latin-1")
        except OSError:
            continue
        # BOTH UNITS. The model prints "Seconds per sim year" when it is quick
        # and "Minutes per sim year" when it is not, so a seconds-only pattern
        # silently loses exactly the expensive rungs whose cost is the reason
        # for measuring -- T85 at dt 15 and T127 at dt 30 both read as no data.
        hit = re.search(r"(Seconds|Minutes) per sim year\s*:?\s*(\d+)", text)
        if hit:
            return float(hit.group(2)) * (60.0 if hit.group(1) == "Minutes" else 1.0)
    return None


def run_job(job: dict, base: dict, log_dir: Path) -> dict:
    """One arm. Never raises: a failure is a result with a diagnosis."""
    name = job["name"]
    from_restart = job["track"] in ("A", "S")
    cfg = write_config(base, name, job["rung"], job["kappa"], job["dt"],
                       job["track"] == "A", bootstrap=not from_restart)
    cmd = [sys.executable, str(ROOT / "exoplasim/scripts/run_exoplasim.py"),
           "--config", str(cfg), "--clean-io", "--run-years", str(job["orbits"])]
    if from_restart:
        # From the spun-up T42 restart. Both arms of every comparison inherit the
        # same superseded surface, which is what that flag is for.
        cmd += ["--restart-from", str(PARENT_RESTART), "--superseded-surface-ok"]
    log = log_dir / f"{name}.log"
    started = now()
    try:
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                              timeout=job.get("timeout_s", 14400))
        text = (proc.stdout or "") + (proc.stderr or "")
        rc = proc.returncode
    except subprocess.TimeoutExpired as exc:
        text = f"TIMEOUT after {exc.timeout}s"
        rc = -9
    elapsed = now() - started
    log.write_text(text)

    run_id = None
    hit = re.search(r"RUN_ID=(run_[0-9a-f]+)", text)
    if hit:
        run_id = hit.group(1)
    trapped = bool(TRAP.search(text))
    crashed = bool(CRASH.search(text))
    result = {**job, "returncode": rc, "wall_s": round(elapsed, 1),
              "run_id": run_id, "trapped": trapped, "crashed": crashed,
              "log": str(log.relative_to(ROOT)) if log.is_relative_to(ROOT) else str(log),
              "outcome": ("late_failure" if trapped and elapsed > LATE_FAILURE_S else
                          "trap" if trapped else
                          "crash" if crashed else
                          "ok" if rc == 0 else "failed")}
    if run_id and (RUNS / run_id).is_dir():
        spo = seconds_per_orbit(RUNS / run_id)
        if spo:
            result["seconds_per_orbit"] = spo
            # REFINE the rung's factor from what actually ran, and say so.
            RUNG_FACTOR[job["rung"]] = (spo * (job["dt"] / 45.0)
                                        / T42_SECONDS_PER_ORBIT_AT_DT45)
            result["rung_factor_measured"] = round(RUNG_FACTOR[job["rung"]], 3)
    return result


def measure_identity(run_id: str, orbits: int) -> dict | None:
    """denergy26 - denergy27 over the arm's own window."""
    cmd = [sys.executable, str(ROOT / "exoplasim/scripts/close_term_energy.py"),
           str(RUNS / run_id), "--first", "0", "--last", str(orbits - 1)]
    try:
        proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                              timeout=3600)
        if proc.returncode != 0:
            return {"error": (proc.stderr or proc.stdout)[-400:]}
        obj = json.loads(proc.stdout[proc.stdout.index("{"):])
        ad = obj["identities"]["adiabatic_energy_conservation"]["residual"]
        ke = obj["identities"]["kinetic_energy_steady_state"]["residual"]
        return {"adiabatic_residual_w_m2": ad["mean_w_m2"],
                "adiabatic_spread_w_m2": ad["spread_w_m2"],
                "kinetic_identity_w_m2": ke["mean_w_m2"]}
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}


def build_jobs() -> list[dict]:
    """The matrix, in priority order: cheapest and most decisive first.

    DEEP AT T42 RATHER THAN WIDE ACROSS THE LADDER, and that is a finding
    rather than a preference. Only T42 carries a full staged surface set --
    every other rung has orography and the land mask alone, and the soil-water
    generator needs a climatology `config/planet.yaml` declares null. Track B
    therefore attempts each rung, fails in seconds, and is abandoned with the
    reason recorded; SPAT-2 is the row that makes those rungs reachable. What
    that frees is spent on the question both rows actually hinge on, which is
    the filter-timestep trade, and which T42 can answer completely.

    Track A extends a grid five arms of which are already measured, so the
    duplicates are omitted rather than re-run.
    """
    jobs = []
    have = {(8.0, 45.0), (2.0, 45.0), (8.0, 22.5), (2.0, 22.5), (None, 22.5)}

    # A. The energy identity across the grid. Three timesteps give the ORDER of
    # the dt scaling rather than a ratio between two points; five filter
    # strengths give the kappa dependence a shape rather than a slope.
    for dt in (45.0, 30.0, 22.5):
        for kappa in (8.0, 4.0, 2.0, 1.0, None):
            if (kappa, dt) in have:
                continue
            tag = "off" if kappa is None else f"k{kappa:g}"
            jobs.append({"track": "A", "rung": "T42", "kappa": kappa, "dt": dt,
                         "orbits": 2, "name": f"A_T42_{tag}_dt{dt:g}",
                         "priority": 10})

    # COST ANCHORS, run early and deliberately out of cost order. SPAT-11's
    # deliverable is what each rung COSTS, and a budget that runs out at T85
    # leaves the top of the ladder unpriced -- which is the one number the whole
    # optimisation workstream rests on and nobody has written down. One orbit at
    # the production step and the production filter, per rung, before the wide
    # sweep gets any of the budget.
    for rung in ("T85", "T127", "T170"):
        jobs.append({"track": "B", "rung": rung, "kappa": 8.0, "dt": 45.0,
                     "orbits": 1, "name": f"C_{rung}_k8_dt45_anchor",
                     "priority": 15})

    # S. THE TRAP BOUNDARY, which is world-qoo's actual deliverable: the
    # weakest filter that still runs at each step. One orbit is enough, because
    # the failure is immediate -- the arm that trapped wrote no output record at
    # all -- so the cells that fail cost seconds and only the survivors cost an
    # orbit. Longer steps than production are included deliberately: 45 is
    # where the trade currently sits and the boundary may be above it.
    for dt in (90.0, 60.0, 45.0, 30.0, 22.5, 15.0):
        for kappa in (8.0, 4.0, 2.0, 1.0, 0.5, None):
            if (kappa, dt) in have or (dt in (45.0, 30.0, 22.5)
                                       and kappa in (8.0, 4.0, 2.0, 1.0, None)):
                continue  # track A already runs that cell and reports its outcome
            tag = "off" if kappa is None else f"k{kappa:g}"
            jobs.append({"track": "S", "rung": "T42", "kappa": kappa, "dt": dt,
                         "orbits": 1, "name": f"S_T42_{tag}_dt{dt:g}",
                         "priority": 20})

    # B. THE LADDER, on bootstrap surface fields and cold. Reachable after all:
    # every rung has orography and the land mask staged, which is exactly what a
    # bootstrap run reads. Cheapest rung first inside each priority, and the
    # deadline decides how far up it gets -- what it does not reach is written
    # out with its estimate rather than dropped quietly. The cost column here is
    # SPAT-11's deliverable.
    ladder = {"T21": (30, (90.0, 60.0, 45.0, 30.0, 22.5, 15.0)),
              "T85": (31, (45.0, 30.0, 22.5, 15.0)),
              "T127": (32, (45.0, 30.0, 22.5)),
              "T170": (33, (45.0, 30.0))}
    for rung, (prio, steps) in ladder.items():
        for dt in steps:
            for kappa in (8.0, None):
                if rung != "T21" and dt == 45.0 and kappa == 8.0:
                    continue  # the cost anchor above already runs this cell
                tag = "off" if kappa is None else f"k{kappa:g}"
                jobs.append({"track": "B", "rung": rung, "kappa": kappa,
                             "dt": dt, "orbits": 1,
                             "name": f"B_{rung}_{tag}_dt{dt:g}",
                             "priority": prio})

    # F. DEEPENING, priority last. The residual's spread is 4 to 13 percent of
    # its value on a two-orbit window; four orbits roughly halves that, and
    # eight halves it again at the two anchor points. This is what the budget
    # buys once the grid is filled, and it is worth more than a wider grid
    # measured badly.
    for kappa in (8.0, 2.0, None):
        for dt in (22.5, 30.0):
            tag = "off" if kappa is None else f"k{kappa:g}"
            jobs.append({"track": "A", "rung": "T42", "kappa": kappa, "dt": dt,
                         "orbits": 4, "name": f"F4_T42_{tag}_dt{dt:g}",
                         "priority": 60})
    for kappa, dt in ((8.0, 45.0), (8.0, 22.5)):
        jobs.append({"track": "A", "rung": "T42", "kappa": kappa, "dt": dt,
                     "orbits": 8, "name": f"F8_T42_k{kappa:g}_dt{dt:g}",
                     "priority": 70})

    jobs.sort(key=lambda j: (j["priority"], estimate(j["rung"], j["dt"], j["orbits"])))
    return jobs


def build_sprint_jobs() -> list[dict]:
    """T170, and the one cell in the matrix that does not fit its pattern.

    The full sweep left T170 unpriced -- every arm trapped, at dt 45 and dt 30,
    filtered and not -- so the ladder's top rung has no cost number and the
    optimisation workstream's central claim still rests on an extrapolation.

    WHERE TO START, and it is arithmetic rather than a guess. T127 runs at dt 30
    and the stable step scales roughly as 1/N for a fixed wind speed, so T170
    wants about 30 * 127/170 = 22.4 minutes. dt 22.5 is therefore the first arm
    and dt 15 the fallback, ordered by PRIORITY rather than by cost so the
    number that is missing arrives first even though it is not the cheapest
    thing here.

    T127 at dt 22.5 is the anomaly: it runs at dt 30 and fails LATE at 22.5,
    after 20.6 minutes, and in `writegp_` rather than `swr_` -- the output
    routine, not the radiation. Non-monotone stability in the timestep should
    not happen, so the first move is to reproduce it and find out whether it
    lands in the same place twice.
    """
    jobs = []
    for prio, rung, kappa, dt in (
            (10, "T170", 8.0, 22.5),
            (11, "T127", 8.0, 22.5),
            (20, "T170", None, 22.5),
            (30, "T170", 8.0, 15.0),
            (40, "T170", None, 15.0),
            (50, "T170", 8.0, 10.0)):
        tag = "off" if kappa is None else f"k{kappa:g}"
        jobs.append({"track": "B", "rung": rung, "kappa": kappa, "dt": dt,
                     "orbits": 1, "name": f"P_{rung}_{tag}_dt{dt:g}",
                     "priority": prio})
    return jobs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hours", type=float, default=6.5,
                    help="wall-clock budget. Jobs are priced before they start "
                         "and skipped when they will not fit what is left.")
    ap.add_argument("--sprint", action="store_true",
                    help="run the focused T170 and T127-anomaly list instead of "
                         "the full matrix. Ordered by priority rather than cost, "
                         "because the missing number is not the cheapest job.")
    ap.add_argument("--min-free-gb", type=float, default=60.0,
                    help="stop scheduling when free disk falls below this")
    args = ap.parse_args()

    out_path = OUT_SPRINT if args.sprint else OUT
    base = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    WORK.mkdir(parents=True, exist_ok=True)
    log_dir = WORK / "logs"
    log_dir.mkdir(exist_ok=True)

    deadline = now() + args.hours * 3600.0
    jobs = build_sprint_jobs() if args.sprint else build_jobs()
    started_at = datetime.now(timezone.utc).isoformat()
    done: list[dict] = []
    skipped: list[dict] = []
    # Declared BEFORE flush() closes over it: a rung that cannot be LAUNCHED is
    # different from one that traps -- the first is this tree's artifact path,
    # the second is the physics being measured -- so the reason is carried in
    # the result rather than inferred from a gap in it.
    abandoned: dict[str, str] = {}

    def flush(note: str = "") -> None:
        payload = {
            "note": "Filter strength x timestep x resolution. CLIM's "
                    "physics-filter row and SPAT-11, which are one question. "
                    "See exoplasim/notes/physics-filter-stability.md.",
            "generated": started_at,
            "finished": datetime.now(timezone.utc).isoformat(),
            "budget_hours": args.hours,
            "status": note,
            "anchor": {"t42_seconds_per_orbit_at_dt45": T42_SECONDS_PER_ORBIT_AT_DT45,
                       "measured_on": "2026-08-23, two arms at each timestep"},
            "rung_factor_current": {k: round(v, 3) for k, v in RUNG_FACTOR.items()},
            "already_measured_elsewhere": [
                {"rung": "T42", "kappa": 8.0, "dt": 45.0, "adiabatic_residual_w_m2": 0.3147},
                {"rung": "T42", "kappa": 2.0, "dt": 45.0, "adiabatic_residual_w_m2": 0.4020},
                {"rung": "T42", "kappa": 8.0, "dt": 22.5, "adiabatic_residual_w_m2": 0.1204},
                {"rung": "T42", "kappa": 2.0, "dt": 22.5, "adiabatic_residual_w_m2": 0.1637},
                {"rung": "T42", "kappa": None, "dt": 22.5, "adiabatic_residual_w_m2": 0.5575},
                {"rung": "T42", "kappa": None, "dt": 45.0, "outcome": "trap"},
            ],
            "abandoned_rungs": abandoned,
            "completed": done,
            "skipped_for_budget": skipped,
        }
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2))

    flush("running")
    for job in jobs:
        if job["rung"] in abandoned:
            skipped.append({**job, "reason": f"rung abandoned: {abandoned[job['rung']]}"})
            continue
        remaining = deadline - now()
        cost = estimate(job["rung"], job["dt"], job["orbits"])
        free_gb = shutil.disk_usage(ROOT).free / 1e9
        if free_gb < args.min_free_gb:
            skipped.append({**job, "reason": f"free disk {free_gb:.0f} GB below floor"})
            continue
        if cost > remaining:
            skipped.append({**job, "estimate_s": round(cost),
                            "remaining_s": round(remaining),
                            "reason": "would not fit the budget"})
            continue
        print(f"[{time.strftime('%H:%M:%S')}] {job['name']} "
              f"est {cost/60:.1f} min, {remaining/60:.0f} min left", flush=True)
        result = run_job(job, base, log_dir)
        if result["outcome"] == "ok" and job["track"] == "A" and result["run_id"]:
            result["identity"] = measure_identity(result["run_id"], job["orbits"])
        if (job["track"] == "B" and result["outcome"] == "failed"
                and not any(d["rung"] == job["rung"] and d["outcome"] != "failed"
                            for d in done)):
            abandoned[job["rung"]] = (
                f"{job['name']} failed to launch (rc {result['returncode']}); "
                f"not a trap, so this is the artifact path rather than the physics")
            print(f"    -> ABANDONING {job['rung']}: {abandoned[job['rung']]}", flush=True)
        done.append(result)
        print(f"    -> {result['outcome']} in {result['wall_s']/60:.1f} min"
              + (f", 26-27 = {result['identity'].get('adiabatic_residual_w_m2')}"
                 if result.get("identity", {}).get("adiabatic_residual_w_m2") is not None else ""),
              flush=True)
        flush("running")

    flush("finished")
    print(f"\n{len(done)} jobs run, {len(skipped)} skipped for budget or disk")
    print(f"wrote {out_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
