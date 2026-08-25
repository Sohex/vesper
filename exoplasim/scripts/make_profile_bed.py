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

A BED IS A COPY OF A NAMELIST THE MODEL HAS MOVED PAST, and that is the failure
this script exists to survive. Every run directory on disk was written by an
older model, so it carries keys the current `namelist` statements no longer
declare and lacks keys the current model requires. Both kill the bed before it
reaches a timestep: an undeclared key aborts in `readnl` with "Cannot match
namelist object name", and a cold start with `tsst_eq` left at its sentinel
aborts in `ice_cold_start`. So this script does two things to every namelist it
copies, and records both in the manifest:

  PRUNE   any key the model's own `namelist` statements do not declare, read out
          of `plasim/src` rather than listed here, so the list cannot go stale.
  FORCE   the keys `config/planet.yaml` declares, to what it says rather than to
          whatever the source run happened to hold.

The manifest records the source run, the sha256 of the executable and of the
restart, every namelist key this script changed and every key it dropped, so a
profile can be traced back to the pair it was taken on.
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

from _paths import CONFIG, COMPONENT_ROOT, INPUTS, MODEL_SRC, MODEL_RUN, RUNS  # noqa: F401

BINARY_MANIFEST = COMPONENT_ROOT / "binary_manifest.json"
MODEL_NAMELIST_SRC = MODEL_SRC / "plasim" / "src"

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
# source run happened to use. Keyed by the namelist FILE, because the model
# reads one group per file and a key written into the wrong one reaches a reader
# that does not declare it.
#
# A bed copies a run directory, so it inherits that run's namelist -- which is
# right for reproducing that run and wrong for the question a bed is asked,
# "what does a FUTURE run cost". run_4182235e9781 was made with the energy
# diagnostics ON; config/planet.yaml turned them off on 2026-08-19 once they had
# answered CLIM-1. Left inherited, every bed accumulated adener3d, a
# (NHOR, NLEV, 28) array, every timestep -- work no current run does. That
# inflates the diagnostics bucket of a profile and understates the transform's
# share against it.
#
# THE THREE icemod_nl ROWS ARE WHAT MAKES A COLD BED STARTABLE. `ice_cold_start`
# refuses a cold start that was given neither a sea surface temperature
# climatology (surface code 169) nor a declared profile, and no run directory
# and no rung under `exoplasim/inputs/` carries a 169 field, so the declared
# profile is the only route a bed has. `ocean.cold_start` is that declaration
# and `run_exoplasim.py` stages the same three keys from it for a live run; a
# bed that did not would start from a different initial condition than the runs
# it stands in for. The number formats match `run_exoplasim.py`'s so the two
# routes write the same line.
CONFIG_FORCED = {
    "plasim_namelist": {
        "NENERGY": (("model", "energy_diagnostics"), lambda v: "1" if v else "0"),
        "NENER3D": (("model", "energy_diagnostics_3d"), lambda v: "1" if v else "0"),
    },
    "icemod_namelist": {
        "TSST_EQ": (("ocean", "cold_start", "sst_equator_k"),
                    lambda v: f"{float(v):.4f}"),
        "TSST_POL": (("ocean", "cold_start", "sst_pole_k"),
                     lambda v: f"{float(v):.4f}"),
        "HICE_INI": (("ocean", "cold_start", "sea_ice_thickness_m"),
                     lambda v: f"{float(v):.4f}"),
    },
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

# One key per line, optionally subscripted: ` TFRC(1) = 20.0 `. That is the shape
# every namelist in a run directory has, because each was written a key at a time
# by `run_exoplasim.py` or by the model's own `write(nud,<group>)`.
NAMELIST_KEY = re.compile(r"^\s*([A-Za-z_]\w*)\s*(?:\([^)]*\))?\s*=")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def config() -> dict:
    return yaml.safe_load(CONFIG.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# What the model declares
# ---------------------------------------------------------------------------

_FILE_VAR = re.compile(
    r"character\s*\(\s*\d+\s*\)\s*::\s*(\w+)\s*=\s*['\"]([^'\"]+)['\"]", re.I)
_OPEN = re.compile(r"open\s*\(\s*(\w+)\s*,\s*file\s*=\s*(\w+)\s*[,)]", re.I)
_READ = re.compile(r"read\s*\(\s*(\w+)\s*,\s*(\w+_nl)\s*\)", re.I)
_NAMELIST_DECL = re.compile(r"^\s*namelist\s*/\s*(\w+)\s*/(.*)$", re.I)


def _uncomment(line: str) -> str:
    if line.lstrip().startswith("!"):
        return ""
    cut = line.find("!")
    return line[:cut] if cut >= 0 else line


def declared_namelist_keys(src: Path = MODEL_NAMELIST_SRC) -> dict[str, set[str]]:
    """Namelist FILE name -> the lowercase keys the model declares in it.

    Read out of `plasim/src` and not listed here, because a list here is a
    second statement of the model's own `namelist` statements and would go
    stale the next time one of them loses a key -- which is exactly the failure
    this function exists to catch. `nguidbg` and `zeta` were both dropped from
    their groups while every run directory on disk kept writing them.

    The file-to-group binding is DERIVED and not assumed. Each module opens a
    unit on a filename variable and then reads a group from that unit, so the
    two statements together say which group belongs to which file; matching
    `<x>_namelist` to `<x>_nl` by name would agree today and would stop being a
    check the first time a module broke the convention. Every group declared in
    the tree pairs this way, and `--verify-namelist-map` prints the pairing.
    """
    files = sorted(src.glob("*.f90"))
    if not files:
        raise SystemExit(f"no model source at {src}")

    var_to_file: dict[str, str] = {}
    for f in files:
        for m in _FILE_VAR.finditer(f.read_text(errors="replace")):
            var_to_file[m.group(1).lower()] = m.group(2)

    groups: dict[str, set[str]] = {}
    file_to_group: dict[str, str] = {}
    for f in files:
        lines = f.read_text(errors="replace").splitlines()
        unit_var: dict[str, str] = {}
        i = 0
        while i < len(lines):
            line = _uncomment(lines[i])
            decl = _NAMELIST_DECL.match(line)
            if decl:
                body = decl.group(2)
                while body.rstrip().endswith("&") and i + 1 < len(lines):
                    body = body.rstrip()[:-1]
                    i += 1
                    body += re.sub(r"^\s*&", "", _uncomment(lines[i]))
                keys = {k.strip().lower()
                        for k in body.replace("&", "").split(",") if k.strip()}
                groups.setdefault(decl.group(1).lower(), set()).update(keys)
            for om in _OPEN.finditer(line):
                unit_var[om.group(1).lower()] = om.group(2).lower()
            for rm in _READ.finditer(line):
                var = unit_var.get(rm.group(1).lower())
                if var in var_to_file:
                    file_to_group[var_to_file[var]] = rm.group(2).lower()
            i += 1

    unpaired = sorted(set(groups) - set(file_to_group.values()))
    if unpaired:
        raise SystemExit(
            f"{src} declares namelist groups nothing reads from a file: "
            f"{', '.join(unpaired)}. A bed cannot be checked against a group "
            f"whose file is unknown; give the group an open/read pair, or say "
            f"here why it has none.")
    return {fname: groups[grp] for fname, grp in file_to_group.items()}


def prune_namelist(text: str, declared: set[str]) -> tuple[str, list[str]]:
    """Drop every key the model does not declare; report what went.

    An undeclared key is not ignored by the reader. gfortran's namelist input
    aborts on the first name it cannot match, so one stale key costs the whole
    run before the first timestep.
    """
    kept, dropped = [], []
    for line in text.splitlines():
        m = NAMELIST_KEY.match(line)
        if m and m.group(1).lower() not in declared:
            dropped.append(m.group(1))
            continue
        kept.append(line)
    return "\n".join(kept) + "\n", dropped


def config_forced(fname: str, cfg: dict) -> dict[str, str]:
    """The CONFIG_FORCED keys for one namelist file, resolved against the config."""
    out = {}
    for key, (path, render) in CONFIG_FORCED.get(fname, {}).items():
        node = cfg
        for part in path:
            if not isinstance(node, dict) or part not in node:
                node = None
                break
            node = node[part]
        if node is not None:
            out[key] = render(node)
    return out


def edit_namelist(text: str, keys: dict[str, str]) -> tuple[str, dict]:
    """Rewrite the given keys in place; report what moved.

    Fortran namelists are order-free, so a key that is absent is appended rather
    than assumed to hold its compiled-in default. NOUTPUT defaults to 1.
    """
    changed = {}
    lines = text.splitlines()
    for key, want in keys.items():
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


def require_settled_bed(dest: Path, cold: bool) -> None:
    """Refuse a bed whose cold start would be seeded from the clock.

    The precondition `exoplasim/scripts/_bed_guard.sh:require_settled_bed` and
    `verify_restart_continuity.py:check_seed` impose on a bed they are handed,
    imposed here on the bed as it is written, so a bed that cannot support a
    comparison is never produced in the first place. `initrandom` takes the
    system clock when `seed(1)` is zero, so two passes over such a bed integrate
    different weather and every comparison built on it reports a broken model.
    """
    if not cold or (dest / "plasim_restart").is_file():
        return
    keys = {}
    for line in (dest / "plasim_namelist").read_text(encoding="utf-8").splitlines():
        m = NAMELIST_KEY.match(line)
        if m:
            keys[m.group(1).upper()] = line.split("=", 1)[1].strip()
    kick = keys.get("KICK", "0").rstrip(",").strip()
    try:
        kicking = int(float(kick)) != 0
    except ValueError:
        kicking = True
    if not kicking:
        return
    seed = (keys.get("SEED") or "").replace(",", " ").split()[:1]
    if seed and seed != ["0"]:
        return
    raise SystemExit(
        f"refusing to write {dest}: a cold bed with KICK = {kick} and no fixed "
        f"SEED draws its initial perturbation from the clock, so two passes "
        f"over it are two different experiments. Pass a non-zero --seed.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-run", type=Path,
                    help="run directory to copy inputs and the executable from")
    ap.add_argument("--dest", type=Path,
                    help="bed directory to create; refuses to overwrite")
    ap.add_argument("--binary",
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
    ap.add_argument("--seed", type=int, default=None,
                    help="SEED for a cold bed's initial kick. Applied only with "
                         "--cold, and defaults to config/planet.yaml's "
                         "model.cold_start_seed, which is what a live run uses. It "
                         "is what makes a cold bed reproducible: with SEED absent, "
                         "initrandom (plasim.f90) takes the system clock, every "
                         "pass integrates different weather, and the restart sha "
                         "stops being able to tell a different model from a "
                         "different kick.")
    ap.add_argument("--mpstep", type=float, default=None,
                    help="override MPSTEP, minutes per timestep. A bed inherits the "
                         "source run's, which was chosen for the source run's "
                         "resolution; a finer grid needs a shorter step and a bed "
                         "that goes unstable profiles nothing.")
    ap.add_argument("--verify-namelist-map", action="store_true",
                    help="print the namelist file to group pairing this script "
                         "derives from the model source, with the key count of "
                         "each, and exit. Builds nothing.")
    args = ap.parse_args()

    declared = declared_namelist_keys()
    if args.verify_namelist_map:
        for fname in sorted(declared):
            print(f"{fname:22s} {len(declared[fname]):3d} keys")
        return

    for required in ("from_run", "dest", "binary"):
        if getattr(args, required) is None:
            raise SystemExit(f"--{required.replace('_', '-')} is required")

    cfg = config()
    seed = args.seed
    if seed is None:
        seed = int(cfg["model"].get("cold_start_seed", 0))

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

    # PRUNE, THEN FORCE, on every namelist the model reads. A namelist file this
    # cannot account for is an error rather than a pass-through, because a key
    # nothing checked is the failure this whole block exists to stop.
    unchecked, changes, dropped = [], {}, {}
    for path in sorted(dest.glob("*_namelist")):
        if path.name not in declared:
            unchecked.append(path.name)
            continue
        text, gone = prune_namelist(path.read_text(encoding="utf-8"),
                                    declared[path.name])
        keys = dict(config_forced(path.name, cfg))
        if path.name == "plasim_namelist":
            keys.update({k: v.format(steps=args.steps)
                         for k, v in BED_NAMELIST.items()})
            if args.mpstep is not None:
                keys["MPSTEP"] = f"{args.mpstep}"
            if args.cold:
                keys["SEED"] = f"{seed}"
        text, moved = edit_namelist(text, keys)
        path.write_text(text, encoding="utf-8")
        if gone:
            dropped[path.name] = gone
        if moved:
            changes[path.name] = moved
    if unchecked:
        raise SystemExit(
            f"{dest} carries namelist files the model source does not account "
            f"for: {', '.join(unchecked)}. They would reach readnl unchecked.")

    require_settled_bed(dest, args.cold)

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
        "namelist_changes": changes,
        "namelist_keys_dropped": dropped,
        "inputs_copied": copied,
        "host": platform.node(),
        "kernel": platform.release(),
        "gcc": subprocess.run(["gcc", "--version"], capture_output=True, text=True
                              ).stdout.splitlines()[0],
    }
    (dest / "bed_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"bed:      {dest}")
    print(f"binary:   {binary.name}  {manifest['binary_sha256'][:16]}")
    print(f"restart:  {'(cold start)' if not restart_sha else restart_sha[:16]}")
    print(f"steps:    {args.steps}")
    for fname in sorted(set(changes) | set(dropped)):
        moves = ", ".join(f"{k} {v['was']}->{v['now']}"
                          for k, v in changes.get(fname, {}).items())
        gone = ", ".join(dropped.get(fname, ()))
        print(f"{fname}: {moves or '(no changes)'}"
              + (f"; dropped undeclared {gone}" if gone else ""))


if __name__ == "__main__":
    main()
