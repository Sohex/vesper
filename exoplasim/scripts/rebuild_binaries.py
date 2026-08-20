#!/usr/bin/env python3
"""Rebuild every ExoPlaSim executable and record which patches each contains.

    python exoplasim/scripts/rebuild_binaries.py            # rebuild all
    python exoplasim/scripts/rebuild_binaries.py --verify   # check, build nothing

ExoPlaSim compiles a separate executable per (resolution, layers, ranks) triple.
Patching the model source and running rebuilds **only the configuration you are
running**; every other binary keeps the old code, silently, until something asks
for it. That is `docs/src/practice/failure-modes.md` class 11 and it fired three times in one
day: the T42 16-rank binary was three days stale, the star-cycle tree was four
days stale and had no 16-rank binary at all, and of six executables on disk
exactly one was newer than the patched source.

That failure is survivable only by luck. The ozone patch added a *namelist key*,
which an old binary cannot parse, so it aborts loudly in `radini_`. A patch that
changed the value or meaning of an existing quantity would have run to
convergence on unpatched physics and nothing would have said so.

The model source is the `vendor/exoplasim` subtree, so there is no patch stack
to police any more: what the binaries are built from is whatever is in that
directory. This script records the sha of every `plasim/src` file that goes into
a build, and `--verify` reports any executable whose recorded sources no longer
match what is on disk. It records the COMPILER AND THE FLAGS beside them, so
an executable whose sha has moved under unchanged sources is attributed to the
toolchain rather than left as unknown provenance.

So: **after any change under `vendor/exoplasim`, rebuild everything.** This
script is that operation, and it writes `exoplasim/binary_manifest.json` mapping
every executable's sha256 to the sources it was built from and the subtree
commit it came from. `scripts/check_consistency.py` reads that manifest and
fails when a binary on disk is absent from it or was built from a source that
has since changed.

Hashing the files rather than trusting the subtree commit is deliberate: a
commit says what was committed, the shas say what is on disk, and an
uncommitted edit under `vendor/` is caught rather than waved through.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import yaml

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "vendor" / "exoplasim" / "exoplasim"
SRC = PKG / "plasim" / "src"
RUN = PKG / "plasim" / "run"
BIN = PKG / "plasim" / "bin"
MANIFEST = ROOT / "exoplasim" / "binary_manifest.json"

# Precision is declared, never inferred: compile.sh takes it in BYTES and
# silently falls back to 4 for anything it does not recognise.
_MODEL = yaml.safe_load(
    (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8")
)["model"]
PRECISION = int(_MODEL["precision_bytes"])

# The flag line is DECLARED IN FULL by `model.compile_flags.f90_opts`, and this script
# writes it into `most_compiler_mpi` before every build.
#
# Not compile.sh's `-O` hook, which is what this used to use. That hook appends
# (`sed '3s/$/ '$optimization'/'`), so it can add a flag and can never remove
# one -- and removal is the operation that matters here, because `-fcheck=all`
# arrives in configure.sh's generated default line and has to come OUT of
# production builds. It also took exactly one flag, since the sed quotes nothing
# and a value with whitespace would corrupt the expression rather than fail.
#
# Writing the whole line also makes a reconfigure harmless: configure.sh may
# regenerate most_compiler_mpi with whatever defaults it likes, and the next
# build overwrites line 3 from the config regardless.
_BUILD = _MODEL["compile_flags"]
BASE_F90_OPTS = [str(f) for f in _BUILD["f90_opts"]]
PROFILES = {k: [str(f) for f in v] for k, v in _BUILD["profiles"].items()}
DEFAULT_PROFILE = str(_BUILD["profile"])
if DEFAULT_PROFILE not in PROFILES:
    raise SystemExit(
        f"model.compile_flags.profile is {DEFAULT_PROFILE!r}, which is not one of "
        f"{sorted(PROFILES)}")
for _flag in BASE_F90_OPTS + [f for v in PROFILES.values() for f in v]:
    if not _flag.startswith("-") or _flag.split() != [_flag]:
        raise SystemExit(
            f"model.compile_flags entries must each be a single token starting with '-', "
            f"got {_flag!r}")

# MOST_F90_OPTS is line 3 of most_compiler_mpi, which compile.sh's own sed also
# assumes. Asserted rather than searched, so a reshaped file fails loudly here
# instead of producing a build with the wrong flags.
F90_OPTS_LINE = 3


def f90_opts(profile: str) -> list[str]:
    return BASE_F90_OPTS + PROFILES[profile]


def write_flag_line(profile: str) -> str:
    """Put the declared flag line into most_compiler_mpi; return what was written."""
    f = PKG / "most_compiler_mpi"
    lines = f.read_text(encoding="utf-8").splitlines()
    if not lines[F90_OPTS_LINE - 1].startswith("MOST_F90_OPTS="):
        raise SystemExit(
            f"{f} line {F90_OPTS_LINE} is not MOST_F90_OPTS=; compile.sh's own "
            f"sed assumes it is, so refusing to guess")
    line = "MOST_F90_OPTS=" + " ".join(f90_opts(profile))
    lines[F90_OPTS_LINE - 1] = line
    # The trailing newline is load-bearing: compile.sh cats this file into the
    # makefile, and without it the next file's first line joins this one and the
    # build fails on a target named after two concatenated variables.
    f.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return line

# The matrix. NLAT must divide by ranks: 32 at T21, 64 at T42, 128 at T85.
MATRIX = [("T21", 10, 8), ("T21", 10, 16),
          ("T42", 10, 8), ("T42", 10, 16),
          ("T85", 10, 16), ("T127", 10, 16)]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def model_sources() -> dict[str, str]:
    """Every file a binary is compiled FROM, with its current sha.

    This is what makes a stale executable detectable. `docs/src/practice/failure-modes.md`
    class 11 is a binary built from source that has since moved, and it has
    fired repeatedly: a T42 binary three days stale, a cycle tree four days
    stale, and of six executables on disk exactly one newer than its source.
    Vendoring the model as a subtree does not fix that -- a stale binary is
    still a stale binary -- so the check survives the migration even though the
    patch stack it used to police does not.

    Hashing the files rather than reading the subtree's commit is deliberate.
    A commit says what was COMMITTED; these shas say what is on disk right now,
    so an uncommitted edit under vendor/ is caught rather than waved through.
    """
    files = (sorted(SRC.glob("*.f90")) + sorted(SRC.glob("*.c"))
             + sorted(SRC.glob("make_*")))
    return {str(f.relative_to(PKG)): sha256(f) for f in files if f.is_file()}


def toolchain(profile: str) -> dict:
    """The compiler and the flags a build is about to use.

    `model_sources()` answers "what source went in". It cannot answer "what
    turned that source into this executable", and the two are separate inputs:
    the same sources through a different compiler, or through the same compiler
    with different flags, give a different binary. Recording only the sources
    means a sha that moved can be NOTICED and not ATTRIBUTED -- `--verify`
    reports unknown provenance and cannot say whether the model changed or the
    toolchain under it did.

    That is not hypothetical here. `-fcheck=all` is in the shipped flags and
    removing it changes the integration, so the flag line is as much a part of
    what a binary IS as the source it was compiled from. See
    `notes/audits/aocl-and-model-build-flags.md`, which also records that this
    build is bit-reproducible: same sources, same flags, same compiler give a
    byte-identical executable, which is what makes the comparison below worth
    making at all.

    `configure.sh` writes the compiler names and options into
    `most_compiler_mpi` and `most_compiler`; the versions are asked of the
    compilers those files actually name, so this stays right if they change.
    """
    args, versions = {}, {}
    for name in ("most_compiler_mpi", "most_compiler", "most_precision_options",
                 "most_precision_optionsx"):
        f = PKG / name
        args[name] = f.read_text(encoding="utf-8") if f.is_file() else None

    named = {}
    for line in (args.get("most_compiler_mpi") or "").splitlines():
        key, _, value = line.partition("=")
        if key in ("MOST_F90", "MOST_CC"):
            named[key] = value.split("#")[0].strip()
    for key, exe in sorted(named.items()):
        try:
            r = subprocess.run([exe, "--version"], capture_output=True, text=True)
            versions[exe] = r.stdout.splitlines()[0].strip() if r.stdout else None
        except (OSError, IndexError):
            versions[exe] = None

    # What compile.sh ACTUALLY fed the makefile. It copies most_compiler_mpi to
    # plasim/bld/compilerargs and appends the -O flag THERE, so the file this
    # function reads above never carries it and recording only that file
    # under-reports the build. Adopting -march=znver4 is what exposed that: the
    # executable carried 20,371 AVX-512 references while the manifest's flag
    # line showed none. This is read after the build, so it describes the
    # binaries that were just written rather than the previous set.
    effective = None
    bld = PKG / "plasim" / "bld" / "compilerargs"
    if bld.is_file():
        for line in bld.read_text(encoding="utf-8").splitlines():
            if line.startswith("MOST_F90_OPTS="):
                effective = line
                break

    # The DECLARED line is what the config says this profile should compile
    # with; the EFFECTIVE line is what the makefile actually used. Recording
    # both is what makes `--verify` able to say whether a moved sha is the
    # config, the toolchain, or neither -- and verifying under a different
    # profile than the binaries were built with is a real drift, reported as
    # one rather than waved through.
    return {"compiler_args": args, "compiler_versions": versions,
            "precision_bytes": PRECISION,
            "build_profile": profile,
            "declared_f90_opts": "MOST_F90_OPTS=" + " ".join(f90_opts(profile)),
            "effective_f90_opts": effective}


def describe_toolchain_drift(prior: dict, current: dict) -> list[str]:
    """What moved between two `toolchain()` readings, as lines a person can act on."""
    if not prior:
        return ["no toolchain recorded (manifest predates the stamping)"]
    out = []
    for name, want in (prior.get("compiler_args") or {}).items():
        got = (current.get("compiler_args") or {}).get(name)
        if got != want:
            out.append(f"{name} differs")
    for exe, want in (prior.get("compiler_versions") or {}).items():
        got = (current.get("compiler_versions") or {}).get(exe)
        if got != want:
            out.append(f"{exe}: built with {want!r}, now {got!r}")
    if prior.get("precision_bytes") != current.get("precision_bytes"):
        out.append(f"precision_bytes: built at {prior.get('precision_bytes')}, "
                   f"config now declares {current.get('precision_bytes')}")
    for key, label in (("declared_f90_opts", "model.compile_flags.f90_opts, as written"),
                       ("build_profile", "model.compile_flags profile"),
                       ("effective_f90_opts", "the flag line the makefile used")):
        if prior.get(key) != current.get(key):
            out.append(f"{label}: built with {prior.get(key)!r}, "
                       f"now {current.get(key)!r}")
    return out


def expected(res: str, lev: int, ranks: int) -> str:
    return f"most_plasim_{res.lower()}_l{lev}_p{ranks}.x"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default=DEFAULT_PROFILE, choices=sorted(PROFILES),
                    help=f"build profile from model.compile_flags.profiles "
                         f"(default {DEFAULT_PROFILE}). 'production' is what "
                         f"orbits are integrated with; 'checked' adds runtime "
                         f"bounds checking and is for anything whose answer is "
                         f"not yet trusted. The profile is part of a binary's "
                         f"identity and is recorded beside its sha.")
    ap.add_argument("--verify", action="store_true",
                    help="report staleness and exit without building")
    args = ap.parse_args()

    if not SRC.is_dir():
        raise SystemExit(f"no ExoPlaSim source at {SRC}; is the vendor/exoplasim subtree present?")

    # `configure.sh` writes most_compiler_mpi and friends, and ExoPlaSim only runs
    # it from `sysconfigure()`. A pip install used to trigger that; a subtree
    # checkout has never been configured, so do it here when the output is
    # absent. It is idempotent.
    if not (PKG / "most_compiler_mpi").is_file():
        print("configuring the vendored source (most_compiler_mpi absent) ...")
        import exoplasim
        exoplasim.sysconfigure()

    rev = subprocess.run(["git", "-C", str(ROOT), "log", "-1", "--format=%h %s",
                          "--", "vendor/exoplasim"],
                         capture_output=True, text=True)
    print(f"model source: vendor/exoplasim @ {rev.stdout.strip() or 'unknown'}")
    dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain",
                            "--", "vendor/exoplasim"],
                           capture_output=True, text=True).stdout.strip()
    if dirty:
        print(f"  UNCOMMITTED changes under vendor/exoplasim ({len(dirty.splitlines())} paths). "
              "The build will use them; the manifest records their shas.")

    sources = model_sources()
    tools = toolchain(args.profile)

    if args.verify:
        stale, absent = [], []
        prior = (json.loads(MANIFEST.read_text(encoding="utf-8"))
                 if MANIFEST.is_file() else {"binaries": {}})
        drift = describe_toolchain_drift(prior.get("toolchain") or {}, tools)
        for exe in sorted(RUN.glob("most_plasim_*.x")):
            rec = prior["binaries"].get(exe.name)
            if rec is None or rec.get("sha256") != sha256(exe):
                absent.append(exe.name)
            elif rec.get("sources") != sources:
                stale.append(exe.name)
        if absent:
            print("NOT IN MANIFEST (provenance unknown):", ", ".join(absent))
        if stale:
            print("BUILT FROM OLDER SOURCE:", ", ".join(stale))
        if drift:
            # Reported whatever the shas say. A toolchain that has moved means
            # the binaries on disk are no longer what this tree would build,
            # even when every one of them still matches the manifest.
            print("TOOLCHAIN HAS MOVED SINCE THE BUILD:")
            for line in drift:
                print(f"  {line}")
        if not absent and not stale and not drift:
            print(f"all {len(list(RUN.glob('most_plasim_*.x')))} binaries current")
        raise SystemExit(1 if (absent or stale or drift) else 0)

    removed = 0
    for d in (RUN, BIN):
        for exe in sorted(d.glob("most_plasim_*.x")):
            exe.unlink()
            removed += 1
    print(f"removed {removed} existing executables\n")

    flag_line = write_flag_line(args.profile)
    print(f"profile {args.profile}: {flag_line}\n")

    built = {}
    for res, lev, ranks in MATRIX:
        name = expected(res, lev, ranks)
        print(f"building {name} ...", flush=True)
        # `-n` is ncpus and `-p` is PRECISION IN BYTES. This passed the rank
        # count to both, so an 8-rank build got -p 8 and was double precision by
        # accident, while a 16-rank build got -p 16, matched no case in
        # compile.sh, and fell back to its default of 4. Every 16-rank binary
        # was therefore single precision while config/planet.yaml declared 8,
        # and ExoPlaSim quietly recompiled at -r8 on first use -- so the
        # manifest described binaries that never ran a single orbit.
        cmd = ["./compile.sh", "-n", str(ranks), "-p", str(PRECISION),
               "-r", res, "-v", str(lev)]
        # No -O. The whole flag line is written into most_compiler_mpi above,
        # because -O can only append and this has to be able to REMOVE.
        r = subprocess.run(cmd, cwd=PKG, capture_output=True, text=True)
        target = RUN / name
        if r.returncode != 0 or not target.is_file():
            print(r.stdout[-1500:])
            print(r.stderr[-1500:])
            raise SystemExit(f"build failed for {name}")
        built[name] = {"sha256": sha256(target),
                       "resolution": res, "layers": lev, "ranks": ranks,
                       "profile": args.profile,
                       "sources": sources}
        print(f"  {built[name]['sha256'][:16]}")

    payload = {
        "note": "Which model source each executable was compiled from. Generated "
                "by exoplasim/scripts/rebuild_binaries.py; do not edit. A binary "
                "absent from here has unknown provenance and must not be "
                "trusted -- rebuild rather than reason about it.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "exoplasim_version": "3.4.2",
        # Recomputed AFTER the build loop: `effective_f90_opts` is read from
        # plasim/bld/compilerargs, which compile.sh writes, so before the
        # loop it describes the previous build.
        "build_profile": args.profile,
        "toolchain": toolchain(args.profile),
        "model_source": subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD:vendor/exoplasim"],
            capture_output=True, text=True).stdout.strip() or None,
        "binaries": built,
    }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {MANIFEST.relative_to(ROOT)} ({len(built)} executables)")
    print("\nThe stellar cycle is a namelist switch (nsolcycle) on these same "
          "binaries; there is no separate cycle executable to build.")


if __name__ == "__main__":
    main()
