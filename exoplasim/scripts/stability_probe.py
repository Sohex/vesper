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

THE DAMPING IS DECLARED BY DEFAULT. `--tau-scale` is 1.0 unless a caller asks
otherwise, so every probe writes the derived NHDIFF, NDEL and TDISS* over all
NLEV levels. It used to write them only when `--tau-scale` was given, and the
2026-08-23 grid was taken without it -- so every cell of that grid ran the
model's own compiled defaults rather than the damping the runs use, which makes
a refusal boundary for a different model. `--inherited-damping` is how a caller
asks for the old behaviour, and it is an arm rather than an oversight.

REFUSAL AND COST ARE TWO VERDICTS. `--refusal-only` takes the first and reports
no cost at all, which is what to use on a contended host: a per-step cost
measured against a shared wall clock is worse than no cost, because it looks
like a measurement. The cost path differences two lengths and needs the machine
to itself.

THE TEMPLATE is a run directory that already has the rung's surface fields
staged and its binary beside them -- a crashed arm serves, since what failed
there was the integration and not the staging. Output is switched off, which is
what keeps a T170 probe from writing twelve gigabytes to measure a step. THIS
PROBE COPIES STAGING AND DOES NOT BUILD IT: with no run directory on disk there
is no bed to be had at any rung, and the surface family under
`exoplasim/inputs/<rung>/` has to exist before a run can be staged to copy from.

AND THE STAGING IS RECORDED BY CONTENT. Every cell carries the sha256 of the
namelist and of each surface `.sra` in its bed, because a cell that names its
executable and not its surface fields is half attributable -- the same probe on
the same binary over a different land mask is a different boundary. Run
directories are untracked and their ids are UUIDs, so `template` alone is a path
nobody can resolve later. `--template` names the run directory outright; the
default takes the most recently modified candidate, which is a property of the
filesystem rather than of the model, and says how many it chose from.

WHICH BINARY, AND HOW IT IS ESTABLISHED. A run directory holds the executable
that was copied into it when that run was staged, so after any rebuild the copy
beside the surface fields is an older model while the registry's is the current
one. That is not rule 4 -- the rebuild was complete -- it is a CONSUMER reading
a copy the manifest does not track, and it cost a cycle on 2026-08-24 when the
probe printed pre-repair surface albedos off a superseded executable and they
read as a patch having failed to take. So the executable is checked against
`exoplasim/binary_manifest.json` before anything is integrated, this REFUSES
when the two disagree, and every result records the name and sha of what ran.
`--restage` takes the registry's current executable instead and says so.
world-anl.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
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
from _paths import INPUTS, MODEL_RUN  # noqa: E402
import build_model  # noqa: E402
import rebuild_binaries  # noqa: E402
from paths import rel  # noqa: E402  from lib/, put on sys.path by _paths

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "exoplasim" / "runs"
WORK = Path("/tmp/vesper-stability-probe")
OUT = ROOT / "exoplasim" / "analysis" / "stability_probe.json"

# THE ORBIT, DERIVED, AND IN THE UNIT lib/orbit.py RETURNS IT IN.
# `orbital_year_days` is Earth days of 24 hours -- that is what
# EARTH_SIDEREAL_YEAR_DAYS scales -- and this file used to multiply it by the
# planet's 30-hour ROTATION instead, which made every orbit 25 per cent too
# long. The model settles it: `run_exoplasim.py` asks for 5850 steps at dt 45
# and 8774 at dt 30 for one orbit, both 4387 hours, where the restated constant
# implied 7312 and 10968. So a per-orbit cost read off this file was a quarter
# too high and every `orbits_covered` a fifth too low.
#
# Read rather than restated for that reason. The model's own calendar rounds
# the orbit to a whole number of 30-hour solar days -- 146 of them, which is
# 4380 hours -- and the 0.16 per cent between that and the true period is below
# anything this file reports.
def _orbit_hours() -> float:
    sys.path.insert(0, str(ROOT / "lib"))
    import orbit as orbit_lib
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    return orbit_lib.orbital_year_days(cfg) * 24.0


ORBIT_HOURS = _orbit_hours()
NLAT = {"T21": 32, "T42": 64, "T85": 128, "T127": 192, "T170": 256}
TRAP = re.compile(r"SIGFPE|Floating-point exception|signal 8")


def steps_per_orbit(dt_minutes: float) -> float:
    return ORBIT_HOURS * 60.0 / dt_minutes


def find_template(rung: str, need_binary: bool) -> Path:
    """A run directory carrying this rung's staged surface fields.

    `need_binary` is the default and asks for the executable beside them too,
    because that is the one the probe will run. Under `--restage` the
    executable comes from the registry, so the template is wanted for its
    STAGING alone and a run directory without a binary serves.

    THE CHOICE IS PRINTED AND COUNTED. Selection is by modification time, which
    is a property of the filesystem and not of the model: with more than one
    candidate the probe would take a different staging on a different day and
    the cell would not say so. `--template` names one outright, and every result
    records the staging's own shas either way, so a cell can be checked rather
    than re-derived. world-qnue and world-37tn.
    """
    want = f"N{NLAT[rung]:03d}_surf_"
    candidates = []
    for d in sorted(RUNS.glob("run_*"), key=lambda p: p.stat().st_mtime, reverse=True):
        if not d.is_dir():
            continue
        if not any(d.glob(f"{want}*.sra")):
            continue
        if need_binary and not any(d.glob(f"most_plasim_t{rung[1:]}_l*_p*.x")):
            continue
        candidates.append(d)
    if candidates:
        if len(candidates) > 1:
            print(f"{len(candidates)} run directories could serve as the {rung} "
                  f"template; taking the most recently modified. Pass --template "
                  f"to name one.")
        return candidates[0]
    also = " and a {} binary".format(rung) if need_binary else ""
    staged = INPUTS / rung.lower()
    raise SystemExit(
        f"no run directory carries {want}*.sra{also}, so there is nothing to "
        f"stage a probe bed from.\n"
        f"  A bed needs a namelist and a staged surface family, and both come "
        f"from a run directory: this probe copies staging, it does not build "
        f"it.\n"
        f"  {rel(staged)} {'exists' if staged.is_dir() else 'DOES NOT EXIST'}, "
        f"and the surface family under it is what run_exoplasim.py stages a run "
        f"from. `python scripts/pipeline.py --status` names the steps that write "
        f"it.\n"
        f"  So: build the surface family for {rung}, run one arm at {rung} -- a "
        f"failing one serves, since what this needs from it is the staging and "
        f"not the integration -- and then re-run this.")


def staging_shas(bed: Path) -> dict:
    """What this bed was staged FROM, by content, not by path.

    A cell that names its executable and not its surface fields is half
    attributable: the same probe on the same binary over a different land mask
    is a different boundary. Run directories are untracked and their ids are
    UUIDs, so `template` alone is a path nobody can resolve three months later.
    These shas can be checked against any staging that still exists.
    """
    out = {}
    for f in sorted(bed.iterdir()):
        if f.is_file() and (f.suffix == ".sra" or f.name.endswith("_namelist")
                            or f.suffix == ".nl"):
            out[f.name] = hashlib.sha256(f.read_bytes()).hexdigest()
    return out


def resolve_binary(rung: str, levels: int, threads: int, template: Path,
                   restage: bool) -> tuple[Path, dict]:
    """The executable this probe will run, and its provenance. CHECKED.

    A run directory holds the executable that was copied into it when that run
    was staged, so its copy is as old as the run and nothing about the probe's
    numbers says which model produced them. Refusing is preferred to running:
    the numbers look ordinary either way, which is
    `docs/src/practice/failure-modes.md` class 34, and a boundary measured on a
    superseded executable is worse than no boundary because it will be quoted.

    Restaging from the registry is the friendlier answer and is available, but
    it is `--restage` and never automatic. Doing both -- checking, then quietly
    substituting -- would leave "which binary was this" exactly as unanswerable
    as it is now.
    """
    if restage:
        exe = MODEL_RUN / build_model.executable_name(rung, levels, threads, False)
        if not exe.is_file():
            raise SystemExit(
                f"--restage wants {exe.name} in {rel(MODEL_RUN)} and it is not "
                f"built. Run exoplasim/scripts/rebuild_binaries.py.")
    else:
        # A registry binary, not a profiling arm: `_fp` carries
        # -fno-omit-frame-pointer and costs -1.17%, so a cost measured on one
        # is not the model's cost. There is one parallel mode, so nothing else
        # to exclude.
        candidates = sorted(b for b in template.glob(f"most_plasim_t{rung[1:]}_l*_p*.x")
                            if not b.name.endswith("_fp.x"))
        if not candidates:
            raise SystemExit(f"{rel(template)} carries no {rung} executable")
        exe = candidates[0]

    problem = rebuild_binaries.unregistered(exe)
    if problem is not None:
        remedy = ("Rebuild (exoplasim/scripts/rebuild_binaries.py) so the "
                  "registry and the manifest agree."
                  if restage else
                  "Pass --restage to take the registry's current executable "
                  "and use the template for its staging alone, or rebuild "
                  "(exoplasim/scripts/rebuild_binaries.py) and restage the run.")
        raise SystemExit(
            f"REFUSING: {problem}.\n"
            f"  taken from {rel(exe.parent)}\n"
            f"  This probe would price and bound a model nobody can name, and "
            f"its numbers would read as current. {remedy}")

    record = rebuild_binaries.registered(exe.name)
    return exe, {"executable": exe.name,
                 "executable_sha256": record["sha256"],
                 "executable_from": "registry" if restage else "template",
                 "build_profile": record.get("profile"),
                 "template": rel(template)}


def build_bed(rung: str, template: Path, tag: str, exe: Path) -> Path:
    bed = WORK / f"bed_{rung}_{tag}"
    if bed.exists():
        shutil.rmtree(bed)
    bed.mkdir(parents=True)
    # THE EXECUTABLE COMES FROM THE REGISTRY, verified against
    # binary_manifest.json by the caller, not from whatever sits in the run
    # directory this bed was staged from. What that selection used to do is why
    # no cell of the pre-2026-08-25 grid can be believed: it preferred an MPI
    # binary, fell through to the threaded one when none existed -- which is
    # every rung now, since rebuild_binaries' MATRIX is five omp rows -- and
    # then launched it under `mpiexec -np 16`. That does not run one model on
    # sixteen threads; it runs SIXTEEN INDEPENDENT MODELS in one directory,
    # over each other's output. world-anl and world-37tn.
    for pattern in ("*_namelist", "*.nl", f"N{NLAT[rung]:03d}_surf_*.sra",
                    "k25v*.dat", "GUI.cfg"):
        for f in template.glob(pattern):
            if f.is_file():
                shutil.copy2(f, bed / f.name)
    # The template's own executable is NOT swept in by the patterns above, so a
    # bed carries exactly the one `resolve_binary` checked. Under --restage the
    # template's copy would otherwise sit beside it under a different name.
    shutil.copy2(exe, bed / exe.name)
    return bed


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


def time_run(bed: Path, exe: str, threads: int) -> tuple[float, bool, str]:
    """ONE PROCESS AND `threads` THREADS, launched directly, on a big stack.

    There is no `mpiexec` here and there must not be: the thread count is
    compiled into the executable, so `mpiexec -np N ./most_plasim_..._pN.x`
    would start N copies of an N-thread binary in one directory over one set of
    restart files.

    `ulimit -s unlimited` through a shell, because the master thread runs on the
    PROCESS stack and no OMP variable moves it: the threaded build is compiled
    `-frecursive`, so the model's large locals are stack-allocated, and above
    T42 the master overruns a 16 MB limit and takes SIGSEGV before writing a
    record. A probe reads that as a refusal and reports the rung as impossible
    at every timestep, which is the shape the grid had for T85 and above.

    The three OMP exports are what `exoplasim/__init__.py` launches production
    under; unbound, libgomp lands threads on SMT siblings and across both dies,
    and this probe reports a per-step cost.
    """
    env = dict(os.environ, OMP_NUM_THREADS=str(threads),
               OMP_PLACES="cores", OMP_PROC_BIND="close",
               OMP_STACKSIZE=os.environ.get("OMP_STACKSIZE", "512M"))
    started = time.monotonic()
    proc = subprocess.run(["bash", "-c", f"ulimit -s unlimited; exec ./{exe}"],
                          cwd=bed, env=env, capture_output=True, text=True,
                          timeout=3600)
    elapsed = time.monotonic() - started
    text = (proc.stdout or "") + (proc.stderr or "")
    return elapsed, bool(TRAP.search(text)) or proc.returncode != 0, text


def probe(rung: str, dt: float, kappa: float | None, steps: int,
          threads: int, template: Path, gamma: int, exe: Path,
          provenance: dict, tau_scale: float | None = None,
          refusal_only: bool = False) -> dict:
    tag = ("off" if kappa is None else f"k{kappa:g}") + f"_dt{dt:g}"
    bed = build_bed(rung, template, tag, exe)
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
    #
    # EVERY LEVEL, WRITTEN OUT. These are per-level arrays and a namelist scalar
    # sets element 1 only, leaving levels 2 upward on the compiled branch's
    # values -- which at T42 are already in SECONDS, so the array is mixed-unit
    # and `dayseccheck` converts none of it. An arm that scales tau by a factor
    # and reaches one level of ten is not the arm it reports. world-720.
    if tau_scale is not None:
        hd = cfg_all["model"]["hyperdiffusion"]["timescales_days"][rung]
        # Every level. These are NLEV arrays and a namelist scalar sets element
        # one only, which is what left nine levels in ten on readnl's presets
        # and, at T42, element one in seconds rather than days. The probe reads
        # a stability boundary, so it has to declare the damping the runs
        # declare or it is measuring a different model. See
        # declare_hyperdiffusion in run_exoplasim.py.
        nlev = int(cfg_all["model"]["layers"])
        # NHDIFF TOO, AND IT IS PER-RUNG. The four timescales say how hard the
        # damping is at the truncation; `nhdiff` is the absolute wavenumber it
        # starts from, and it is derived as `cutoff_fraction * ntru` -- 8 at
        # T21, 16 at T42, 32 at T85. A bed staged from a T21 run carries 8, so
        # a probe that wrote the timescales and left this alone damped a T42
        # bed from T21's wavenumber: 8 is 38% of T21's spectrum and 19% of
        # T42's, and the operator that results is not the one any run of that
        # rung uses. That is the same contamination the timescales had, one
        # key over. `declare_hyperdiffusion` in run_exoplasim.py is where the
        # runs get it, and the two derivations are the same line.
        hyper = cfg_all["model"]["hyperdiffusion"]
        ntru = int(rung.lstrip("Tt"))
        nhdiff = int(round(float(hyper["cutoff_fraction"]) * ntru))
        keys |= {"TDISSD": f"{nlev}*{hd['divergence'] / tau_scale}",
                 "TDISSZ": f"{nlev}*{hd['vorticity'] / tau_scale}",
                 "TDISST": f"{nlev}*{hd['temperature'] / tau_scale}",
                 "TDISSQ": f"{nlev}*{hd['humidity'] / tau_scale}",
                 "NHDIFF": f"{nhdiff}",
                 "NDEL": f"{nlev}*{int(hyper['order_alpha'])}"}
    # WHICH BINARY, AND WHAT DAMPING, IN EVERY RESULT. A refusal boundary and a
    # per-step cost are both properties of one executable measured under one
    # damping, and a result naming neither can only be re-derived, never
    # checked. The 2026-08-23 grid recorded no filter power and no
    # hyperdiffusion, and both had moved by the time anyone read it; separately,
    # no entry named a binary at all. world-anl and world-37tn, and the two
    # contaminations are independent.
    declared = {"gamma": gamma, "tau_scale": tau_scale,
                "hyperdiffusion": "inherited" if tau_scale is None else "derived",
                "per_level": tau_scale is not None,
                "nhdiff": None if tau_scale is None else nhdiff,
                "ndel": None if tau_scale is None else int(
                    cfg_all["model"]["hyperdiffusion"]["order_alpha"]),
                "tdiss_days": None if tau_scale is None else {
                    "divergence": hd["divergence"] / tau_scale,
                    "vorticity": hd["vorticity"] / tau_scale,
                    "temperature": hd["temperature"] / tau_scale,
                    "humidity": hd["humidity"] / tau_scale}}
    first_steps = steps if refusal_only else short_steps
    set_keys(bed, keys | {"N_RUN_STEPS": str(first_steps)})
    staging = staging_shas(bed)
    t_short, trapped, text = time_run(bed, exe.name, threads)
    result = {"rung": rung, "dt_minutes": dt, "kappa": kappa, "threads": threads,
              "steps_short": first_steps, "steps_long": steps,
              # HOW MUCH OF AN ORBIT THIS CELL ACTUALLY COVERS, in the units the
              # endurance question is asked in. A refusal verdict is a statement
              # about the first fraction of an orbit and nothing else, and the
              # fraction is small: 600 steps at dt 45 is 0.08 orbit, while the
              # T42 blow-up this grid exists beside happened in orbit 47. Left
              # to be re-derived, that ratio is exactly what a reader does not
              # do -- `docs/src/practice/failure-modes.md` class 34 -- so the
              # cell carries it.
              "orbits_covered": round(steps / steps_per_orbit(dt), 4),
              "declared": declared,
              "outcome": "refused" if trapped else "no_refusal_in_steps",
              "staging_sha256": staging,
              **provenance}
    # A CONTENDED HOST CAN TAKE THE REFUSAL AND NOT THE COST. world-37tn.
    if not refusal_only:
        result["wall_short_s"] = round(t_short, 2)
    if trapped:
        if not refusal_only:
            result["failed_after_s"] = round(t_short, 2)
        frames = [ln.strip() for ln in text.splitlines() if re.match(r"^#\d+ ", ln.strip())]
        result["backtrace"] = frames[:6]
        shutil.rmtree(bed, ignore_errors=True)
        return result

    if refusal_only:
        # NO COST FIELD AT ALL, not a cost with a caveat beside it. A number in
        # this row would be quoted.
        shutil.rmtree(bed, ignore_errors=True)
        return result

    set_keys(bed, keys | {"N_RUN_STEPS": str(steps)})
    t_long, trapped_long, _ = time_run(bed, exe.name, threads)
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
    ap.add_argument("--tau-scale", type=float, default=1.0,
                    help="multiply the derived hyperdiffusion STRENGTH by this "
                         "(so tau is divided by it). 1 is the cascade-absorbing "
                         "value and the DEFAULT; larger is stronger damping.")
    ap.add_argument("--inherited-damping", action="store_true",
                    help="write no hyperdiffusion keys at all, leaving every "
                         "level on plasimmod's compiled defaults. This is what "
                         "contaminated the 2026-08-23 grid, so it is an "
                         "explicit arm and never a default. world-37tn.")
    ap.add_argument("--refusal-only", action="store_true",
                    help="one probe, refusal verdict, NO cost. Use this when "
                         "the host is contended: a per-step cost taken against "
                         "a shared wall clock is worse than none, because it "
                         "looks like a measurement.")
    ap.add_argument("--gamma", type=int, default=None,
                    help="filter power; default is config/planet.yaml's")
    ap.add_argument("--threads", type=int, default=16,
                    help="the thread count the binary was compiled for")
    ap.add_argument("--restage", action="store_true",
                    help="take the executable from the model's run/ directory "
                         "rather than from the template, and use the template "
                         "for its staging alone. The default REFUSES when the "
                         "template's copy is not what binary_manifest.json "
                         "registers; this is the other answer, and it is asked "
                         "for rather than applied silently.")
    ap.add_argument("--template", type=Path, default=None,
                    help="the run directory to stage every bed from. The "
                         "default takes the most recently modified one that "
                         "carries this rung's surface fields, which is a "
                         "property of the filesystem rather than of the model; "
                         "name one here when a grid has to be repeatable.")
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    cfg_model = yaml.safe_load(
        (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))["model"]
    if args.gamma is None:
        args.gamma = int(cfg_model["filter_power"])
    dts = ([float(x) for x in args.sweep.split(",")] if args.sweep
           else [args.dt if args.dt else 45.0])
    kappas = [None if k.strip().lower() == "off" else float(k)
              for k in args.kappa.split(",")]
    WORK.mkdir(parents=True, exist_ok=True)
    if args.template is not None:
        template = args.template.resolve()
        if not template.is_dir():
            raise SystemExit(f"--template {template} is not a directory")
        if not any(template.glob(f"N{NLAT[args.rung]:03d}_surf_*.sra")):
            raise SystemExit(
                f"--template {rel(template)} carries no "
                f"N{NLAT[args.rung]:03d}_surf_*.sra, so it cannot stage a "
                f"{args.rung} bed")
    else:
        template = find_template(args.rung, need_binary=not args.restage)
    exe, provenance = resolve_binary(args.rung, int(cfg_model["layers"]),
                                     args.threads, template, args.restage)
    print(f"template: {rel(template)}")
    print(f"binary:   {exe.name} {provenance['executable_sha256'][:16]} "
          f"(profile {provenance['build_profile']}, "
          f"from the {provenance['executable_from']})")

    results = []
    for dt in dts:
        for kappa in kappas:
            r = probe(args.rung, dt, kappa, args.steps, args.threads, template,
                      args.gamma, exe, provenance,
                      None if args.inherited_damping else args.tau_scale,
                      args.refusal_only)
            r["tau_scale"] = None if args.inherited_damping else args.tau_scale
            results.append(r)
            k = "off" if kappa is None else f"{kappa:g}"
            head = f"  {args.rung} kappa {k:>3s} dt {dt:5.1f}"
            if r["outcome"] == "refused":
                # NO WALL CLOCK UNDER --refusal-only, and the line has to say so
                # rather than read it: `failed_after_s` is a cost field and the
                # refusal path deliberately writes none, so reading it here
                # killed the sweep on its FIRST refusing cell -- which is the
                # one cell a refusal sweep exists to find, and the one whose
                # result was then never written to the grid.
                when = (f" after {r['failed_after_s']:.1f} s"
                        if "failed_after_s" in r else "")
                print(f"{head}  REFUSED{when}", flush=True)
            elif r["outcome"] == "refused_only_at_length":
                # The cell this grid exists to find: it starts, and dies inside
                # the probe's own range. Reported as its own verdict because a
                # step that refuses and a step that fails late are different
                # facts about a rung, and collapsing them is what let T42 at dt
                # 45 be read as usable when it dies in its forty-seventh orbit.
                print(f"{head}  RAN {r['steps_short']} steps, REFUSED at "
                      f"{r['steps_long']} ({r['wall_long_s']:.1f} s)",
                      flush=True)
            elif args.refusal_only:
                # No cost field on a contended host, so nothing to report but
                # the verdict. world-37tn.
                print(f"{head}  no refusal in {r['steps_long']} steps",
                      flush=True)
            else:
                print(f"{head}  no refusal in {r['steps_long']} steps, "
                      f"{r['seconds_per_step']:.4f} s/step -> "
                      f"{r['implied_seconds_per_orbit']/60:.1f} min/orbit  "
                      f"(startup {r['startup_s']:.0f} s; naive single run would "
                      f"say {r['naive_single_run_seconds_per_orbit']/60:.1f})", flush=True)
            prior = json.loads(args.out.read_text()) if args.out.is_file() else {"probes": []}
            prior.setdefault("probes", [])
            prior["probes"] = [p for p in prior["probes"]
                               if not (p["rung"] == r["rung"] and p["dt_minutes"] == r["dt_minutes"]
                                       and p["kappa"] == r["kappa"])] + [r]
            prior["note"] = (
                "Refusal boundary and per-step cost. Short probes: this cannot "
                "see a LATE blow-up. Each probe names the executable that "
                "produced it and the sha binary_manifest.json registers for it; "
                "an entry without `executable_sha256` predates that and cannot "
                "be attributed to a model source. "
                "exoplasim/scripts/stability_probe.py")
            prior["generated"] = datetime.now(timezone.utc).isoformat()
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(json.dumps(prior, indent=2))
    print(f"\nwrote {rel(args.out)}")


if __name__ == "__main__":
    main()
