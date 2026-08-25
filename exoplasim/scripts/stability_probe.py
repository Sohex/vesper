#!/usr/bin/env python3
"""Where does the model refuse to start, and what does a step cost? In minutes.

    python exoplasim/scripts/stability_probe.py --rung T170 --dt 22.5 --steps 300
    python exoplasim/scripts/stability_probe.py --rung T170 --sweep 45,30,22.5,15,10

Worldbuilding frame: a COMPUTE and stability probe on the Vesper project's
climate model. Nothing here is about the simulated planet.

WHY THIS EXISTS RATHER THAN ANOTHER ORBIT. The matrix answered the same two
questions by integrating whole orbits, which costs 80 minutes at T170 and gave
one cell per 80 minutes. Neither question needs an orbit:

  COST is linear in step count. Measured at four rungs across six timesteps --
  T21 19 s at dt 90 rising to 112 at dt 15, T85 295 at dt 45 to 900 at dt 15 --
  every one within a few percent of proportional. So one short probe gives the
  per-step cost, and every other timestep at that rung follows by arithmetic.
  A rung needs ONE measurement, not one per timestep.

  THE REFUSAL happens at the first radiation call. The trapping arms wrote no
  output record at all and died in under a tenth of a minute, so three hundred
  steps is already two orders of magnitude more than that failure needs.

WHAT IT CANNOT SEE, stated so nothing reads more into it than it holds: the
LATE blow-up. T127 at dt 22.5 integrated twenty minutes and 4.3 GB before a
field exceeded single precision in the output writer. A short probe passes that
cell. So this maps the refusal boundary and prices the rung; whether a
configuration survives a whole orbit is a longer question and this does not
answer it.

THE TEMPLATE is a run directory that already has the rung's surface fields
staged and its binary beside them -- a crashed arm serves, since what failed
there was the integration and not the staging. Output is switched off, which is
what keeps a T170 probe from writing twelve gigabytes to measure a step.
"""
from __future__ import annotations

import argparse
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
from paths import rel  # noqa: E402  from lib/, put on sys.path by _paths

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "exoplasim" / "runs"
WORK = Path("/tmp/vesper-stability-probe")
OUT = ROOT / "exoplasim" / "analysis" / "stability_probe.json"

# Hours in a Vesper orbit: 182.801 d x 30.0 h. Both come from config/planet.yaml
# through derive(); restated here only to turn a per-step cost into a per-orbit
# one, and checked against the matrix's own measured orbits below.
ORBIT_HOURS = 182.801 * 30.0
NLAT = {"T21": 32, "T42": 64, "T85": 128, "T127": 192, "T170": 256}
TRAP = re.compile(r"SIGFPE|Floating-point exception|signal 8")


def steps_per_orbit(dt_minutes: float) -> float:
    return ORBIT_HOURS * 60.0 / dt_minutes


def find_template(rung: str) -> Path:
    """A run directory carrying this rung's staged surface fields and binary."""
    want = f"N{NLAT[rung]:03d}_surf_"
    best = None
    for d in sorted(RUNS.glob("run_*"), key=lambda p: p.stat().st_mtime, reverse=True):
        if not d.is_dir():
            continue
        if not any(d.glob(f"{want}*.sra")):
            continue
        if not any(d.glob(f"most_plasim_t{rung[1:]}_l*_p*.x")):
            continue
        best = d
        break
    if best is None:
        raise SystemExit(
            f"no run directory carries {want}*.sra and a {rung} binary. "
            f"Run one arm at {rung} first, even a failing one: what this needs "
            f"from it is the staging, not the integration.")
    return best


def build_bed(rung: str, template: Path, tag: str) -> tuple[Path, str]:
    bed = WORK / f"bed_{rung}_{tag}"
    if bed.exists():
        shutil.rmtree(bed)
    bed.mkdir(parents=True)
    binary = sorted(template.glob(f"most_plasim_t{rung[1:]}_l*_p*.x"))
    # The MPI binary, not the threaded one: the matrix measured on MPI and a
    # cost compared across parallel modes is not a cost comparison.
    binary = [b for b in binary if not b.name.endswith("_omp.x")] or binary
    exe = binary[0]
    for pattern in ("*_namelist", "*.nl", f"N{NLAT[rung]:03d}_surf_*.sra",
                    "k25v*.dat", "GUI.cfg"):
        for f in template.glob(pattern):
            if f.is_file():
                shutil.copy2(f, bed / f.name)
    shutil.copy2(exe, bed / exe.name)
    return bed, exe.name


def set_keys(bed: Path, keys: dict[str, str]) -> None:
    path = bed / "plasim_namelist"
    text = path.read_text(encoding="latin-1")
    for key, value in keys.items():
        pattern = re.compile(rf"^\s*{key}\s*=.*$", re.MULTILINE | re.IGNORECASE)
        line = f" {key} = {value}"
        if pattern.search(text):
            text = pattern.sub(line, text, count=1)
        else:
            # Namelists are order-free, so an absent key is APPENDED rather than
            # left to a compiled-in default nobody named.
            text = text.rstrip()
            text = text[:text.rfind("/")] + line + "\n/\n"
    path.write_text(text, encoding="latin-1")


def time_run(bed: Path, exe: str, ranks: int) -> tuple[float, bool, str]:
    started = time.monotonic()
    proc = subprocess.run(["mpiexec", "-np", str(ranks), f"./{exe}"], cwd=bed,
                          capture_output=True, text=True, timeout=3600)
    elapsed = time.monotonic() - started
    text = (proc.stdout or "") + (proc.stderr or "")
    return elapsed, bool(TRAP.search(text)) or proc.returncode != 0, text


def probe(rung: str, dt: float, kappa: float | None, steps: int,
          ranks: int, template: Path, gamma: int,
          tau_scale: float | None = None) -> dict:
    tag = ("off" if kappa is None else f"k{kappa:g}") + f"_dt{dt:g}"
    bed, exe = build_bed(rung, template, tag)
    # SEED is declared, not inherited: the template is a run directory made
    # before `model.cold_start_seed` existed, and `initrandom` falls back to the
    # system clock when seed(1) is zero. A probe nobody can re-run is a boundary
    # nobody can check.
    cfg_all = yaml.safe_load(
        (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    seed = int(yaml.safe_load(
        (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8")
    )["model"]["cold_start_seed"])
    keys = {"N_RUN_STEPS": str(steps), "NOUTPUT": "0", "NSNAPSHOT": "0",
            "NDIAG": "0", "MPSTEP": f"{dt}", "NENERGY": "0", "NENER3D": "0",
            "SEED": str(seed)}
    if kappa is None:
        keys |= {"NFILTER": "0", "NGPTFILTER": "0", "NSPVFILTER": "0"}
    else:
        # NFILTEREXP too: kappa sets the damping AT the truncation and gamma
        # sets how far down it reaches, and leaving gamma to whatever the
        # template carried would vary the confinement between arms that differ
        # only in kappa.
        keys |= {"NFILTER": "2", "NGPTFILTER": "1", "NSPVFILTER": "1",
                 "FILTERKAPPA": f"{kappa}", "NFILTEREXP": f"{gamma}"}
    # TWO LENGTHS, AND THE SLOPE BETWEEN THEM. A single short run divides the
    # STARTUP cost -- Legendre setup, FFTW planning, staging -- over its own few
    # hundred steps, while an orbit divides it over fourteen thousand. Measured
    # here: a 200-step probe of T170 at dt 22.5 implied 69.2 minutes an orbit
    # where the full orbit took 51, an overstatement of 36% that is entirely
    # startup. Differencing two lengths cancels it, which is the correction
    # `docs/src/practice/failure-modes.md` class 34 exists to demand -- a bed
    # shorter than its own startup measuring the startup.
    short_steps = max(50, steps // 3)
    # HYPERDIFFUSION TOO, and scaled. With the filter off this becomes the
    # model's only damping, which is the whole point of being able to vary it:
    # the strength the core actually demands cannot be measured while a second
    # mechanism supplies several hundred times more.
    if tau_scale is not None:
        hd = cfg_all["model"]["hyperdiffusion"]["timescales_days"][rung]
        # Every level. These are NLEV arrays and a namelist scalar sets element
        # one only, which is what left nine levels in ten on readnl's presets
        # and, at T42, element one in seconds rather than days. The probe reads
        # a stability boundary, so it has to declare the damping the runs
        # declare or it is measuring a different model. See
        # declare_hyperdiffusion in run_exoplasim.py.
        nlev = int(cfg_all["model"]["layers"])
        keys |= {"TDISSD": f"{nlev}*{hd['divergence'] / tau_scale}",
                 "TDISSZ": f"{nlev}*{hd['vorticity'] / tau_scale}",
                 "TDISST": f"{nlev}*{hd['temperature'] / tau_scale}",
                 "TDISSQ": f"{nlev}*{hd['humidity'] / tau_scale}",
                 "NDEL": f"{nlev}*{int(cfg_all['model']['hyperdiffusion']['order_alpha'])}"}
    set_keys(bed, keys | {"N_RUN_STEPS": str(short_steps)})
    t_short, trapped, text = time_run(bed, exe, ranks)
    result = {"rung": rung, "dt_minutes": dt, "kappa": kappa, "ranks": ranks,
              "steps_short": short_steps, "steps_long": steps,
              "wall_short_s": round(t_short, 2),
              "outcome": "refused" if trapped else "no_refusal_in_steps"}
    if trapped:
        result["failed_after_s"] = round(t_short, 2)
        frames = [ln.strip() for ln in text.splitlines() if re.match(r"^#\d+ ", ln.strip())]
        result["backtrace"] = frames[:6]
        shutil.rmtree(bed, ignore_errors=True)
        return result

    set_keys(bed, keys | {"N_RUN_STEPS": str(steps)})
    t_long, trapped_long, _ = time_run(bed, exe, ranks)
    result["wall_long_s"] = round(t_long, 2)
    if trapped_long:
        # Ran short and refused long: that is a LATE failure inside the probe's
        # own range, and worth more than either number.
        result["outcome"] = "refused_only_at_length"
        shutil.rmtree(bed, ignore_errors=True)
        return result
    per_step = (t_long - t_short) / (steps - short_steps)
    result["seconds_per_step"] = round(per_step, 5)
    result["startup_s"] = round(t_long - per_step * steps, 2)
    result["implied_seconds_per_orbit"] = round(per_step * steps_per_orbit(dt), 1)
    result["naive_single_run_seconds_per_orbit"] = round(
        (t_long / steps) * steps_per_orbit(dt), 1)
    shutil.rmtree(bed, ignore_errors=True)
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rung", required=True, choices=sorted(NLAT))
    ap.add_argument("--dt", type=float, default=None)
    ap.add_argument("--sweep", default=None,
                    help="comma-separated timesteps in minutes, cheapest first")
    ap.add_argument("--kappa", default="8",
                    help="filter strength, or 'off'. Comma-separated to sweep.")
    ap.add_argument("--steps", type=int, default=600,
                    help="timesteps per probe. The refusal fires on the first "
                         "radiation call, so this is already far more than that "
                         "failure needs; it is long enough to price a step.")
    ap.add_argument("--tau-scale", type=float, default=None,
                    help="multiply the derived hyperdiffusion STRENGTH by this "
                         "(so tau is divided by it). 1 is the cascade-absorbing "
                         "value; larger is stronger damping.")
    ap.add_argument("--gamma", type=int, default=None,
                    help="filter power; default is config/planet.yaml's")
    ap.add_argument("--ranks", type=int, default=16)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    if args.gamma is None:
        args.gamma = int(yaml.safe_load(
            (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8")
        )["model"]["filter_power"])
    dts = ([float(x) for x in args.sweep.split(",")] if args.sweep
           else [args.dt if args.dt else 45.0])
    kappas = [None if k.strip().lower() == "off" else float(k)
              for k in args.kappa.split(",")]
    WORK.mkdir(parents=True, exist_ok=True)
    template = find_template(args.rung)
    print(f"template: {rel(template)}")

    results = []
    for dt in dts:
        for kappa in kappas:
            r = probe(args.rung, dt, kappa, args.steps, args.ranks, template,
                      args.gamma, args.tau_scale)
            r["tau_scale"] = args.tau_scale
            results.append(r)
            k = "off" if kappa is None else f"{kappa:g}"
            if r["outcome"] != "refused":
                print(f"  {args.rung} kappa {k:>3s} dt {dt:5.1f}  no refusal in "
                      f"{r['steps_long']} steps, "
                      f"{r['seconds_per_step']:.4f} s/step -> "
                      f"{r['implied_seconds_per_orbit']/60:.1f} min/orbit  "
                      f"(startup {r['startup_s']:.0f} s; naive single run would "
                      f"say {r['naive_single_run_seconds_per_orbit']/60:.1f})", flush=True)
            else:
                print(f"  {args.rung} kappa {k:>3s} dt {dt:5.1f}  REFUSED "
                      f"after {r['failed_after_s']:.1f} s", flush=True)
            prior = json.loads(args.out.read_text()) if args.out.is_file() else {"probes": []}
            prior.setdefault("probes", [])
            prior["probes"] = [p for p in prior["probes"]
                               if not (p["rung"] == r["rung"] and p["dt_minutes"] == r["dt_minutes"]
                                       and p["kappa"] == r["kappa"])] + [r]
            prior["note"] = ("Refusal boundary and per-step cost. Short probes: "
                             "this cannot see a LATE blow-up. "
                             "exoplasim/scripts/stability_probe.py")
            prior["generated"] = datetime.now(timezone.utc).isoformat()
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(json.dumps(prior, indent=2))
    print(f"\nwrote {rel(args.out)}")


if __name__ == "__main__":
    main()
