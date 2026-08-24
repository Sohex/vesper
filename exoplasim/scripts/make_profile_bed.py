#!/usr/bin/env python3
"""Build a short-run working copy of a run directory for profiling.

    python exoplasim/scripts/make_profile_bed.py --from-run <run_dir> --dest <bed>

Worldbuilding frame: this is a COMPUTE instrument. It measures where the Vesper
climate model spends its cycles on this desktop, and nothing it produces is
about the simulated planet.

The bed is `notes/audits/aocl-and-model-build-flags.md`'s: a copy of a run
directory with `N_RUN_STEPS` cut and every output stream switched off, so the
instrument reads model compute and not postprocessing. It copies the inputs and
the executable and leaves the source run untouched; nothing here writes into
`exoplasim/runs/`.

Two modes, chosen by whether a restart of the right shape exists:

  warm  the source run's `plasim_restart` is copied in, and the model resumes a
        spun-up state. The only honest bed for the resolution that run was made
        at, because the transform share depends on nothing about the state but
        the physics branches taken do.
  cold  no restart is copied, so `restart_ini` finds none and the model starts
        from `KICK`. This is what a resolution with no run behind it gets, and
        it is adequate for a profile: which routines run, and how much work each
        carries, is set by the resolution and the namelist, not by the state.

The manifest records the source run, the sha256 of the executable and of the
restart, and every namelist key this script changed, so a profile can be traced
back to the pair it was taken on.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import yaml

from _paths import CONFIG, COMPONENT_ROOT, INPUTS, MODEL_RUN, RUNS  # noqa: F401

BINARY_MANIFEST = COMPONENT_ROOT / "binary_manifest.json"

# Output off, and short. NOUTPUT and NSNAPSHOT are the two streams; NDIAG is the
# diagnostic print interval and is already 0 in production, set here so a bed
# built from some other run cannot inherit one.
BED_NAMELIST = {
    "N_RUN_STEPS": "{steps}",
    "NOUTPUT": "0",
    "NSNAPSHOT": "0",
    "NDIAG": "0",
}

# Keys the CONFIG declares, forced to what it says rather than to whatever the
# source run happened to use.
#
# A bed copies a run directory, so it inherits that run's namelist -- which is
# right for reproducing that run and wrong for the question a bed is asked,
# "what does a FUTURE run cost". run_4182235e9781 was made with the energy
# diagnostics ON; config/planet.yaml turned them off on 2026-08-19 once they had
# answered CLIM-1. Left inherited, every bed accumulated adener3d, a
# (NHOR, NLEV, 28) array, every timestep -- work no current run does. That
# inflates the diagnostics bucket of a profile and understates the transform's
# share against it.
CONFIG_FORCED = {
    "NENERGY": ("energy_diagnostics", lambda v: "1" if v else "0"),
    "NENER3D": ("energy_diagnostics_3d", lambda v: "1" if v else "0"),
}

# Everything a plasim run directory needs that is not output. Globs, because the
# surface field set is resolution-dependent (N064_* at T42, N128_* at T85).
BED_INPUTS = (
    "*_namelist",
    "*.nl",
    "N???_surf_*.sra",
    "k25v*.dat",
    "GUI.cfg",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def config_forced() -> dict[str, str]:
    """The CONFIG_FORCED keys resolved against config/planet.yaml."""
    model = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))["model"]
    out = {}
    for key, (cfg_key, render) in CONFIG_FORCED.items():
        if cfg_key in model:
            out[key] = render(model[cfg_key])
    return out


def edit_namelist(text: str, steps: int, mpstep: float | None,
                  seed: int | None) -> tuple[str, dict]:
    """Rewrite the bed keys in place; report what moved.

    Fortran namelists are order-free, so a key that is absent is appended rather
    than assumed to hold its compiled-in default. NOUTPUT defaults to 1.
    """
    changed = {}
    lines = text.splitlines()
    keys = dict(BED_NAMELIST)
    keys.update(config_forced())
    if mpstep is not None:
        keys["MPSTEP"] = f"{mpstep}"
    if seed is not None:
        keys["SEED"] = f"{seed}"
    for key, value in keys.items():
        want = value.format(steps=steps)
        pattern = re.compile(rf"^(\s*){key}(\s*)=(\s*)(\S+)(\s*)$", re.IGNORECASE)
        for i, line in enumerate(lines):
            m = pattern.match(line)
            if m:
                if m.group(4) != want:
                    changed[key] = {"was": m.group(4), "now": want}
                lines[i] = f"{m.group(1)}{key} = {want} "
                break
        else:
            end = next(i for i, l in enumerate(lines) if l.strip().upper() == "/END")
            lines.insert(end, f" {key} = {want} ")
            changed[key] = {"was": None, "now": want}
    return "\n".join(lines) + "\n", changed


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-run", type=Path, required=True,
                    help="run directory to copy inputs and the executable from")
    ap.add_argument("--dest", type=Path, required=True,
                    help="bed directory to create; refuses to overwrite")
    ap.add_argument("--binary", required=True,
                    help="executable name, e.g. "
                         "most_plasim_t<res>_l<layers>_p<ranks>.x")
    ap.add_argument("--binary-dir", type=Path, default=MODEL_RUN,
                    help="where to take the executable from (default: the vendored "
                         "model's run/ directory, which is what binary_manifest.json "
                         "registers and what a new run would launch). A run directory "
                         "holds the executable ITS orbits were made with, which after "
                         "any rebuild is not the same file.")
    ap.add_argument("--allow-unregistered", action="store_true",
                    help="proceed with an executable absent from binary_manifest.json")
    ap.add_argument("--steps", type=int, default=600,
                    help="N_RUN_STEPS for the bed (default 600, the AOCL bed's)")
    ap.add_argument("--cold", action="store_true",
                    help="do not copy plasim_restart; start from KICK")
    ap.add_argument("--surface-from", type=Path, default=None,
                    help="directory of orogen_T??_surf_NNNN.sra to use instead of "
                         "the run directory's, renamed to the N<NLAT>_surf_NNNN.sra "
                         "the model asks for (surfmod.f90:185). Needed for any "
                         "resolution the source run was not made at.")
    ap.add_argument("--nlat", type=int, default=None,
                    help="latitude count, required with --surface-from; it is what "
                         "the model builds the surface filename from")
    ap.add_argument("--seed", type=int, default=88888888,
                    help="SEED for a cold bed's initial kick. Applied only with "
                         "--cold, and it is what makes a cold bed reproducible: "
                         "with SEED absent, initrandom (plasim.f90:2025) takes the "
                         "system clock, every pass integrates different weather, and "
                         "the restart sha stops being able to tell a different model "
                         "from a different kick.")
    ap.add_argument("--mpstep", type=float, default=None,
                    help="override MPSTEP, minutes per timestep. A bed inherits the "
                         "source run's, which was chosen for the source run's "
                         "resolution; a finer grid needs a shorter step and a bed "
                         "that goes unstable profiles nothing.")
    args = ap.parse_args()

    src = args.from_run.resolve()
    dest = args.dest.resolve()
    if dest.exists():
        raise SystemExit(f"{dest} exists; delete it or pick another --dest")
    binary = args.binary_dir.resolve() / args.binary
    if not binary.is_file():
        raise SystemExit(f"no executable {args.binary} in {args.binary_dir}")

    # An executable absent from the manifest has unknown provenance, which is the
    # manifest's own rule. A profile is a claim about a specific pair of source
    # and flags, so it is worth less than nothing if the pair cannot be named.
    binary_sha = sha256(binary)
    registered = json.loads(BINARY_MANIFEST.read_text())["binaries"].get(args.binary, {})
    if registered.get("sha256") != binary_sha:
        msg = (f"{args.binary} sha {binary_sha[:16]} is not what binary_manifest.json "
               f"registers ({str(registered.get('sha256'))[:16]}); rebuild rather than "
               f"reason about it, or pass --allow-unregistered")
        if not args.allow_unregistered:
            raise SystemExit(msg)
        print(f"WARNING: {msg}")

    if args.surface_from and not args.nlat:
        raise SystemExit("--surface-from needs --nlat")

    dest.mkdir(parents=True)
    copied = []
    for pattern in BED_INPUTS:
        if args.surface_from and pattern == "N???_surf_*.sra":
            continue  # supplied below, from another resolution
        for f in sorted(src.glob(pattern)):
            if f.is_file():
                shutil.copy2(f, dest / f.name)
                copied.append(f.name)

    if args.surface_from:
        for f in sorted(args.surface_from.resolve().glob("orogen_T*_surf_*.sra")):
            code = int(f.stem.rsplit("_", 1)[1])
            name = f"N{args.nlat:03d}_surf_{code:04d}.sra"
            shutil.copy2(f, dest / name)
            copied.append(f"{name} <- {f.name}")

    shutil.copy2(binary, dest / binary.name)

    restart_sha = None
    if not args.cold:
        restart = src / "plasim_restart"
        if not restart.is_file():
            raise SystemExit(f"no plasim_restart in {src}; use --cold for a fresh start")
        shutil.copy2(restart, dest / "plasim_restart")
        restart_sha = sha256(restart)

    namelist = dest / "plasim_namelist"
    text, changed = edit_namelist(namelist.read_text(), args.steps, args.mpstep,
                                  args.seed if args.cold else None)
    namelist.write_text(text)

    manifest = {
        "note": "Profiling bed. A short, output-free copy of a run directory; "
                "not a climate run, and nothing here reaches a climatology.",
        "created": datetime.now(timezone.utc).isoformat(),
        "source_run": str(src),
        "binary": binary.name,
        "binary_from": str(binary.parent),
        "binary_sha256": binary_sha,
        "binary_registered": registered.get("sha256") == binary_sha,
        "restart_sha256": restart_sha,
        "cold_start": args.cold,
        "steps": args.steps,
        "namelist_changes": changed,
        "inputs_copied": copied,
        "host": platform.node(),
        "kernel": platform.release(),
        "gcc": subprocess.run(["gcc", "--version"], capture_output=True, text=True
                              ).stdout.splitlines()[0],
    }
    (dest / "bed_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"bed:      {dest}")
    print(f"binary:   {binary.name}  {manifest['binary_sha256'][:16]}")
    print(f"restart:  {'(cold start)' if args.cold else restart_sha[:16]}")
    print(f"steps:    {args.steps}")
    moves = ", ".join(f"{k} {v['was']}->{v['now']}" for k, v in changed.items())
    print(f"namelist: {moves or '(already a bed)'}")


if __name__ == "__main__":
    main()
