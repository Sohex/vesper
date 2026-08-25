#!/usr/bin/env python3
"""Rebuild every ExoPlaSim executable and record which patches each contains.

    python exoplasim/scripts/rebuild_binaries.py            # rebuild all
    python exoplasim/scripts/rebuild_binaries.py --verify   # check, build nothing

ExoPlaSim compiles a separate executable per (resolution, layers, ranks)
configuration -- the geometry and the thread count are `parameter`s compiled in,
not settings a binary reads. Patching the model source and running rebuilds
**only the configuration you are running**; every other binary keeps the old
code, silently, until something asks for it. That is `docs/src/practice/failure-modes.md` class 11 and it fired three times in one
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

import concurrent.futures as cf

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_model
import restart_schema

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "vendor" / "exoplasim" / "exoplasim"
SRC = PKG / "plasim" / "src"
RUN = PKG / "plasim" / "run"
BIN = PKG / "plasim" / "bin"
MANIFEST = ROOT / "exoplasim" / "binary_manifest.json"

# Precision is declared, never inferred. It reaches the compiler through
# build_model.py, which refuses a value it does not recognise rather than
# falling back to single -- archive CLIM-22 is that fallback having shipped
# twelve single-precision binaries under a declaration of eight.
_MODEL = yaml.safe_load(
    (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8")
)["model"]
PRECISION = int(_MODEL["precision_bytes"])

# The flag line is DECLARED IN FULL by `model.compile_flags.f90_opts` and reaches
# the compiler through `build_model.py`, which reads the same declaration. There
# is no generated compiler-options file to keep in step any more: the three that
# existed disagreed, and only one of them was ever written from this declaration.
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

def f90_opts(profile: str) -> list[str]:
    return BASE_F90_OPTS + PROFILES[profile]


# THREADED_FLAGS restates what `build_model.build` adds on top of the declared
# line, and the restatement is checked rather than trusted:
# `assert_threaded_flags_agree()` below fails the run if the two lists ever
# diverge. Recording a flag line the compiler did not receive is the defect
# `toolchain()` documents, and it reached the manifest once already.
THREADED_FLAGS = ["-fopenmp", "-DOMPSHARED"]


def effective_opts(profile: str) -> list[str]:
    """The flag line the compiler actually receives for this configuration."""
    return (f90_opts(profile) + THREADED_FLAGS
            + (["-fdefault-real-8"] if PRECISION == 8 else []))


def assert_threaded_flags_agree() -> None:
    """Check THREADED_FLAGS against the source that adds them.

    `build_model.build` is where a flag reaches the compiler. This module only
    RECORDS what went in, so the two can drift apart silently and the manifest
    would keep describing a build nobody made. Read the flags out of the one
    that does the work.
    """
    src = Path(build_model.__file__).read_text(encoding="utf-8")
    want = "flags = flags + [" + ", ".join(repr(f) for f in THREADED_FLAGS) + "]"
    if want.replace("'", '"') not in src.replace("'", '"'):
        raise SystemExit(
            f"THREADED_FLAGS is {THREADED_FLAGS}, and build_model.py no longer "
            f"contains {want!r}. One of the two has moved; the manifest must "
            f"record the line the compiler receives.")


# The matrix. NLAT must divide the rank count, and 8, 16 and 32 all divide every
# NLAT on the T21/T42/T85/T127/T170 ladder: 32, 64, 128, 192, 256. Every rung has
# an export under `source/<build>/exoplasim-<T>`, and a rung without a binary is
# terrain nothing can run.
#
# The 32-thread rows exist for the SMT arms of
# `exoplasim/notes/smt-rank-layout.md`, and the 8-rank rows at the high
# resolutions exist because several of those arms are concurrent 8-rank jobs. A
# rank-layout matrix needs every rank count it compares, at every resolution it
# compares them at, built from ONE source through ONE flag line -- so they are
# registered here rather than built off to one side. A binary absent from the
# manifest has unknown provenance, and that note's whole argument is a
# comparison between binaries.
# THE THREADED ROWS. `exoplasim/notes/non-mpi-performance-roadmap.md` puts
# "establish the completed OpenMP/SHTns build as the new baseline" first in its
# measurement order, and CLIM-59, CLIM-74, CLIM-84 and SPAT-11 all measure on
# that build -- so it is registry work, not a side build. Sixteen threads at
# every rung and no other count: `exoplasim/notes/thread-count-by-resolution.md`
# measures sixteen as the winner at T42 with every lower count monotonically
# worse (-36.0% at eight), `smt-rank-layout.md` refuses everything above it at
# every resolution, and NPRO must divide NLAT. The sweep binaries that settled
# that are arms, built by `thread_count_sweep.sh` and not registered; these are
# the configurations something is meant to RUN.
# ONE RUNTIME, EXERCISED. The model has one parallel layer, threads over a
# shared address space, and it runs SHTns: `nshtns` is 1 here, and SHTns is 5.1x
# legmod's transform on this host (`exoplasim/notes/shtns-viability.md`). A
# second registered runtime that nothing exercises is how `writesp`,
# `writescalar` and `writecolumn` kept an unguarded write for as long as they
# did -- harmless under distributed ranks, and corruption under threads.
# EVERY ROW IS A RUNG SOMETHING IS MEANT TO RUN. None of them is an
# executable-only cost artifact, and the question "is this rung runnable" is not
# answered by the presence of its binary. Three things have to line up, and each
# refuses by name when it does not: this matrix builds the executable,
# `config/planet.yaml` declares the rung's hyperdiffusion timescales and its
# timestep -- `declare_hyperdiffusion` refuses a rung with no entry rather than
# letting it inherit T21's damping -- and the surface family for the rung has to
# be staged under `exoplasim/inputs/<rung>/`, which `surface_field_report`
# checks before the model is reached. A rung with a binary and no surface family
# is not a runnable rung; it is an unfinished one, and `scripts/pipeline.py
# --status` with that rung configured is what says which. world-3oj.
MATRIX = [("T21", 10, 16),
          ("T42", 10, 16),
          ("T85", 10, 16),
          ("T127", 10, 16),
          ("T170", 10, 16)]


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

    The set is READ OUT OF `plasim/CMakeLists.txt`, not globbed from
    `plasim/src`. A glob hashes files no configuration compiles, so editing one
    reports every executable as built from changed source when no compiled line
    moved -- and files that compile nowhere are kept on purpose here, the
    offline glacier accelerator among them. The one configurable slot left is
    PLASIM_FFT, and it is expanded to every value the build file admits, because
    this manifest spans the configurations rule 4 counts binaries over and one
    recorded source set has to cover them all.
    """
    names = restart_schema.compiled_modules(PKG / "plasim", (".f90", ".c"))
    files = [SRC / n for n in sorted(names)]
    missing = [n for n, f in zip(sorted(names), files) if not f.is_file()]
    if missing:
        raise SystemExit(f"{SRC} is missing sources CMakeLists.txt lists: "
                         f"{', '.join(missing)}")
    return {str(f.relative_to(PKG)): sha256(f) for f in files}


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

    There is one declaration and it is `config/planet.yaml`. What used to be
    recorded here was the contents of three generated compiler-options files, of
    which only one was ever written from that declaration -- so a threaded binary
    carried a provenance record of flags it was not built with. The flag line
    below is the one `build_model.py` hands the compiler, and the versions are
    asked of the compilers it actually invokes.
    """
    versions = {}
    for exe in ("gfortran", "gcc", "cmake", "ninja"):
        try:
            r = subprocess.run([exe, "--version"], capture_output=True, text=True)
            versions[exe] = r.stdout.splitlines()[0].strip() if r.stdout else None
        except (OSError, IndexError):
            versions[exe] = None

    # DECLARED and EFFECTIVE are now the same list plus the precision flag,
    # and both are recorded because they were once able to disagree: the old
    # driver appended `-O` to a COPY of the options file, so the file this
    # recorded never carried it and the manifest under-reported the build.
    return {"compiler_versions": versions,
            "precision_bytes": PRECISION,
            "build_profile": profile,
            "declared_f90_opts": " ".join(f90_opts(profile)),
            # THE LINE THE COMPILER RECEIVED, declared plus the threading and
            # precision flags build_model.py appends. Recording the declared line
            # alone is how a binary came to carry a provenance record of flags it
            # was not built with, which is the defect this function's docstring is
            # about. Each binary also records its own line, and this is the summary.
            "effective_f90_opts": " ".join(effective_opts(profile)),
            "note": "The registry builds every configuration in MATRIX. A "
                    "profiling build adds -fno-omit-frame-pointer and is an arm "
                    "rather than a registry entry. build_model.py is where a flag "
                    "reaches the compiler."}


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
                       ("effective_f90_opts",
                        "the flag line the compiler received")):
        if prior.get(key) != current.get(key):
            out.append(f"{label}: built with {prior.get(key)!r}, "
                       f"now {current.get(key)!r}")
    return out


def expected(res: str, lev: int, ranks: int) -> str:
    """The registry's name for a configuration.

    Deliberately a SECOND statement of the rule `build_model.executable_name`
    states, rather than a call to it. The two are cross-checked against each
    other after every build below, and a check that asks one function twice
    cannot fail -- `docs/src/practice/failure-modes.md` class 17.
    """
    return f"most_plasim_{res.lower()}_l{lev}_p{ranks}.x"


def verify(profile: str = DEFAULT_PROFILE) -> dict:
    """What `--verify` reports, as data, so a caller can act on it.

    Builds nothing and links nothing: a manifest read and a sha per source.
    `scripts/check_consistency.py` calls this rather than restating it, because
    a second implementation of "is this executable current" is a second thing to
    keep in step and the two had already drifted -- that one globbed `RUN`
    alone, so an executable in `BIN` kept a provenance nothing looked at, and it
    could not see a MATRIX row with no binary at all or a toolchain that had
    moved under an unchanged source.

    Four separate findings, kept apart because they are acted on differently:

      absent   on disk and not in the manifest, or in it under another sha:
               unknown provenance
      stale    in the manifest, built from a source file that has since moved
      missing  named by MATRIX and not built at all -- rule 4's other half, and
               invisible to anything that globs the directory
      drift    the compiler, the flag line or the declared precision is not what
               the manifest records, so the binaries are no longer what this
               tree would build even where every sha still matches
    """
    prior = (json.loads(MANIFEST.read_text(encoding="utf-8"))
             if MANIFEST.is_file() else {"binaries": {}})
    sources = model_sources()
    drift = describe_toolchain_drift(prior.get("toolchain") or {}, toolchain(profile))
    # BOTH DIRECTORIES, because the rebuild clears both. Globbing RUN alone let
    # an executable in BIN keep a provenance nothing looked at, which is rule
    # 4's failure mode one directory over: the verify must cover exactly the set
    # the build replaces. world-60x0.
    built = sorted(list(RUN.glob("most_plasim_*.x"))
                   + list(BIN.glob("most_plasim_*.x")))
    absent, stale = [], []
    for exe in built:
        rec = prior.get("binaries", {}).get(exe.name)
        if rec is None or rec.get("sha256") != sha256(exe):
            absent.append(exe.name)
        elif rec.get("sources") != sources:
            stale.append(exe.name)
    # A binary on disk built from moved source is one failure; a MATRIX entry
    # with NO binary at all is the mirror of it, and globbing the directory
    # cannot see that one. It is the same class 11 shape: silent until something
    # asks for the configuration. The ladder in `MATRIX` is the declaration of
    # what must exist, so check against it.
    missing = [expected(r, l, k) for r, l, k in MATRIX
               if not (RUN / expected(r, l, k)).is_file()]
    return {"manifest_present": MANIFEST.is_file(),
            "built": [e.name for e in built],
            "absent": absent, "stale": stale,
            "missing": missing, "drift": drift}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default=DEFAULT_PROFILE, choices=sorted(PROFILES),
                    help=f"build profile from model.compile_flags.profiles "
                         f"(default {DEFAULT_PROFILE}). 'production' is what "
                         f"orbits are integrated with; 'checked' adds runtime "
                         f"bounds checking and is for anything whose answer is "
                         f"not yet trusted. The profile is part of a binary's "
                         f"identity and is recorded beside its sha.")
    ap.add_argument("--jobs", type=int, default=4,
                    help="configurations to build at once (default 4). Each one "
                         "compiles in its own directory, so they do not collide; "
                         "each also runs its own parallel compile underneath.")
    ap.add_argument("--verify", action="store_true",
                    help="report staleness and exit without building")
    args = ap.parse_args()

    if not SRC.is_dir():
        raise SystemExit(f"no ExoPlaSim source at {SRC}; is the vendor/exoplasim subtree present?")

    assert_threaded_flags_agree()

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
        report = verify(args.profile)
        built = report["built"]
        absent, stale = report["absent"], report["stale"]
        missing, drift = report["missing"], report["drift"]
        if missing:
            print("IN THE MATRIX AND NOT BUILT:", ", ".join(missing))
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
        if not absent and not stale and not drift and not missing:
            print(f"all {len(built)} binaries current")
        # `missing` SETS THE EXIT CODE. It was printed and then dropped from the
        # verdict, so `--verify` announced "IN THE MATRIX AND NOT BUILT" and
        # exited 0 -- and rule 4 is the rule that a rebuild only refreshes the
        # configuration you ran and leaves the rest stale in silence, which is
        # precisely what a configuration with no binary at all is.
        raise SystemExit(1 if (absent or stale or drift or missing) else 0)

    removed = 0
    for d in (RUN, BIN):
        for exe in sorted(d.glob("most_plasim_*.x")):
            exe.unlink()
            removed += 1
    print(f"removed {removed} existing executables\n")

    flag_line = "MOST_F90_OPTS=" + " ".join(f90_opts(args.profile))
    print(f"profile {args.profile}: {flag_line}\n")

    # CONCURRENTLY, which is new and is the whole reason the build directory is
    # per-configuration. The old driver compiled every configuration in one
    # `plasim/bld` that it emptied on entry, so two builds could not run at once
    # without one deleting the other's objects and this loop had to be serial.
    # `notes/audits/aocl-and-model-build-flags.md` measures parallelism inside a
    # single build at 3.5x; this is the other half, across the twelve.
    built = {}
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {
            pool.submit(build_model.build, res, lev, ranks,
                        args.profile, False, None, False): (res, lev, ranks)
            for res, lev, ranks in MATRIX
        }
        for fut in cf.as_completed(futures):
            res, lev, ranks = futures[fut]
            name = expected(res, lev, ranks)
            try:
                target = fut.result()
            except SystemExit as exc:
                raise SystemExit(f"build failed for {name}: {exc}") from exc
            if target.name != name:
                raise SystemExit(
                    f"build_model produced {target.name} where the registry "
                    f"expects {name}; the two naming rules have diverged")
            built[name] = {"sha256": sha256(target),
                           "resolution": res, "layers": lev, "ranks": ranks,
                           "profile": args.profile,
                           # PER BINARY. The top-level toolchain summarises; this
                           # is the line THIS executable was compiled with.
                           "effective_f90_opts": " ".join(
                               effective_opts(args.profile)),
                           "sources": sources}
            print(f"built {name}  {built[name]['sha256'][:16]}", flush=True)
    built = {k: built[k] for k in sorted(built)}

    payload = {
        "note": "Which model source each executable was compiled from. Generated "
                "by exoplasim/scripts/rebuild_binaries.py; do not edit. A binary "
                "absent from here has unknown provenance and must not be "
                "trusted -- rebuild rather than reason about it.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "exoplasim_version": "3.4.2",
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
