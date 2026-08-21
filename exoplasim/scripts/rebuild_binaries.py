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

# compile.sh's -O hook appends this to MOST_F90_OPTS. Declared in the config
# rather than edited into the generated `most_compiler_mpi`, which configure.sh
# rewrites; `toolchain()` records the resulting flag line either way, but only a
# declared flag survives a reconfigure. compile.sh strips no dash and quotes
# nothing (`sed '3s/$/ '$optimization'/'`), so it takes exactly one flag and a
# value with whitespace would corrupt the sed expression rather than fail.
OPTIMIZATION = str(_MODEL.get("optimization_flag") or "").strip()
if OPTIMIZATION.split() != ([OPTIMIZATION] if OPTIMIZATION else []):
    raise SystemExit(
        f"model.optimization_flag must be a single flag, got {OPTIMIZATION!r}")
if OPTIMIZATION and not OPTIMIZATION.startswith("-"):
    raise SystemExit(
        f"model.optimization_flag must start with '-', got {OPTIMIZATION!r}")

# The matrix. NLAT must divide by ranks: 32 at T21, 64 at T42, 128 at T85,
# 192 at T127, 256 at T170. The ladder is T21/T42/T85/T127/T170 and every rung
# has an export under `source/<build>/exoplasim-<T>`; a rung without a binary is
# terrain nothing can run. T127 and T170 take p16 for the same reason T85 does,
# 16 dividing 192 and 256 at 12 and 16 rows a rank; p32 also divides both, at 6
# and 8, and is the option if throughput at the top of the ladder matters more
# than the compile.
MATRIX = [("T21", 10, 8), ("T21", 10, 16),
          ("T42", 10, 8), ("T42", 10, 16),
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
    """
    files = (sorted(SRC.glob("*.f90")) + sorted(SRC.glob("*.c"))
             + sorted(SRC.glob("make_*")))
    return {str(f.relative_to(PKG)): sha256(f) for f in files if f.is_file()}


def toolchain() -> dict:
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

    return {"compiler_args": args, "compiler_versions": versions,
            "precision_bytes": PRECISION,
            "optimization_flag": OPTIMIZATION or None,
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
    for key, label in (("optimization_flag", "model.optimization_flag"),
                       ("effective_f90_opts", "the flag line the makefile used")):
        if prior.get(key) != current.get(key):
            out.append(f"{label}: built with {prior.get(key)!r}, "
                       f"now {current.get(key)!r}")
    return out


def expected(res: str, lev: int, ranks: int) -> str:
    return f"most_plasim_{res.lower()}_l{lev}_p{ranks}.x"


def main() -> None:
    ap = argparse.ArgumentParser()
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
    tools = toolchain()

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
        # A binary on disk built from moved source is one failure; a MATRIX
        # entry with NO binary at all is the mirror of it, and globbing the
        # directory cannot see that one. It is the same class 11 shape: silent
        # until something asks for the configuration. The ladder in `MATRIX` is
        # the declaration of what must exist, so check against it.
        missing = [expected(r, l, k) for r, l, k in MATRIX
                   if not (RUN / expected(r, l, k)).is_file()]
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
            print(f"all {len(list(RUN.glob('most_plasim_*.x')))} binaries current")
        raise SystemExit(1 if (absent or stale or drift) else 0)

    removed = 0
    for d in (RUN, BIN):
        for exe in sorted(d.glob("most_plasim_*.x")):
            exe.unlink()
            removed += 1
    print(f"removed {removed} existing executables\n")

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
        if OPTIMIZATION:
            # -O takes the flag WITHOUT its leading dash and prepends one.
            cmd += ["-O", OPTIMIZATION[1:]]
        r = subprocess.run(cmd, cwd=PKG, capture_output=True, text=True)
        target = RUN / name
        if r.returncode != 0 or not target.is_file():
            print(r.stdout[-1500:])
            print(r.stderr[-1500:])
            raise SystemExit(f"build failed for {name}")
        built[name] = {"sha256": sha256(target),
                       "resolution": res, "layers": lev, "ranks": ranks,
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
        "toolchain": toolchain(),
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
