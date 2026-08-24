#!/usr/bin/env python3
"""Build one PlaSim executable, and refuse anything it does not understand.

    python exoplasim/scripts/build_model.py --res T170 --ranks 16
    python exoplasim/scripts/build_model.py --res T21 --ranks 8 --profile checked

Worldbuilding frame: this builds the Vesper project's climate model, a hard fork
of ExoPlaSim. Nothing here is about the simulated planet.

WHAT THIS REPLACED, and why the replacement is not a tidy-up. `compile.sh` and
`configure.sh` between them defaulted silently in three places that each produced
a DIFFERENT BUILD rather than an error: an unrecognised `-r` built T21, an
unrecognised `-p` built single precision, and the executable was then named for
the resolution the parse arrived at, so a caller copying the name it asked for
got an untouched binary from a previous session. Every one of those is a check
that could only pass. `notes/audits/model-build-driver.md` has the measurements;
archive CLIM-22 is the precision one having already cost a session.

So the rule here is that every argument is validated against a list and an
unrecognised value exits non-zero having written nothing. There is no default
resolution, no default rank count and no default precision.

WHERE THE FLAGS COME FROM. `config/planet.yaml`, `model.compile_flags`, and
nowhere else. The old build had three hand-edited files -- `most_compiler`,
`most_compiler_mpi`, `most_compiler_omp` -- of which only one was written from
the declaration, so the threaded build, which is what this project measures,
took whatever was in a file nothing derived. That was CONS-14.

ONE PARALLEL LAYER. This project builds and runs threads over a shared address
space and nothing else, so there is no parmode argument, no parmode in the build
tag and no parmode in the executable name: an axis with one value separates no
two builds and only invites a caller to ask for a configuration that does not
exist. `-fopenmp -DOMPSHARED` therefore reaches every compile below. world-38b.

WHERE IT BUILDS. One directory per configuration under `vendor/exoplasim/build/`,
so two configurations cannot delete each other's objects. The old build had one
`plasim/bld` that it emptied on entry, which is why `rebuild_binaries.py` had to
be serial across its twelve executables and why the tree needed marker files to
notice it held the wrong configuration. Both go away.

AND NOT THERE WHEN THE SOURCE IS PATCHED. A verification arm that corrupts
`plasim/src` in place builds under `build/patched/` with a hash of what it
patched in the directory name, and cannot publish. `patched_sources` carries
that argument.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

from _paths import CONFIG, MODEL_RUN, PROJECT_ROOT  # noqa: E402
import rungs  # noqa: E402  -- the ladder registry; _paths put lib on the path

PKG = PROJECT_ROOT / "vendor" / "exoplasim" / "exoplasim"
PLASIM = PKG / "plasim"
BUILD_ROOT = PROJECT_ROOT / "vendor" / "exoplasim" / "build"
# WHERE A PATCHED SOURCE BUILDS. Never under BUILD_ROOT itself. See
# `patched_sources` for the whole argument.
PATCHED_ROOT = BUILD_ROOT / "patched"
SHTNS_PREFIX = PROJECT_ROOT / "vendor" / "shtns-install"

# The markers a verification arm writes into the model source before it corrupts
# it. `scripts/smoke_test.py:check_no_control_patch` refuses to let either be
# committed; this file refuses to let either reach the registry's build tag.
CONTROL_MARKERS = ("CONTROL PATCH IN PROGRESS", "! CONTROL:")

# Resolution name -> latitudes, from `lib/gridding.py`, which is the one place
# the ladder is written down and checks each rung's name against its own grid.
# A value outside it must be an ERROR: `-r 170` used to build T21 and say
# nothing. SPAT-2.
RESOLUTIONS = dict(rungs.RUNGS)


def declared() -> dict:
    return yaml.safe_load(CONFIG.read_text(encoding="utf-8"))["model"]


def flag_line(profile: str) -> tuple[list[str], int]:
    """The declared Fortran flags for a profile, and the declared precision."""
    model = declared()
    build = model["compile_flags"]
    profiles = build["profiles"]
    if profile not in profiles:
        raise SystemExit(
            f"--profile {profile!r} is not one of {sorted(profiles)}; "
            f"they are declared in config/planet.yaml under model.compile_flags.profiles")
    flags = [str(f) for f in build["f90_opts"]] + [str(f) for f in profiles[profile]]
    for f in flags:
        if not f.startswith("-") or f.split() != [f]:
            raise SystemExit(
                f"config/planet.yaml declares {f!r}, which is not a single token "
                f"starting with '-'. One flag per entry.")
    return flags, int(model["precision_bytes"])


def resolve(res: str) -> tuple[str, int]:
    """Accept a name or a latitude count; refuse anything else."""
    if res in RESOLUTIONS:
        return res, RESOLUTIONS[res]
    by_lat = {v: k for k, v in RESOLUTIONS.items()}
    if res.isdigit() and int(res) in by_lat:
        return by_lat[int(res)], int(res)
    raise SystemExit(
        f"--res {res!r} is not a resolution. Names: {' '.join(RESOLUTIONS)}. "
        f"Latitude counts: {' '.join(str(v) for v in RESOLUTIONS.values())}.")


def patched_sources() -> list[Path]:
    """The model source files a verification arm has deliberately corrupted.

    WHY THE BUILD HAS TO ASK. The tag below names every input that changes a
    byte of the executable EXCEPT the state of the model source, so an arm that
    patches `plasim/src` in place builds into the same directory the registry
    entry claims. That happened: `verify_shtns_model.sh` strips the spectral
    filter out of `shtnsmod.f90` for its control arm, and
    `build/t21_l10_p16_production/plasim.x` was left holding a binary with
    zero filter sites while the live source had eleven. `--no-publish` keeps an
    arm out of the model run directory; it does nothing about the build
    directory, and `--print-path` hands the caller exactly that path. The
    `drop+extra` hash exists to keep flag arms apart, and a source patch is an
    arm the tag could not see. world-70k.

    THE MARKER IS THE SIGNAL, and it is the one that already exists: a control
    patch writes `CONTROL PATCH IN PROGRESS` or `! CONTROL:` into the file
    BEFORE it corrupts it, precisely so a tree in that state is recognisable,
    and `smoke_test.check_no_control_patch` refuses to let one be committed.
    Reading the same marker here means there is one convention rather than two.
    It cannot see a patch that only deletes code, which is why the rule is that
    a control adds its marker first.
    """
    src = PLASIM / "src"
    if not src.is_dir():
        return []
    out = []
    for path in sorted(src.glob("*.f90")):
        text = path.read_text(errors="replace")
        if any(m in text for m in CONTROL_MARKERS):
            out.append(path)
    return out


def source_state(patched: list[Path]) -> str:
    """A short hash of the patched files, so two control arms do not collide.

    The same shape as the extra-flag hash: an arm's identity is in its
    directory's name, so two arms cannot land in one directory and quietly
    reuse each other's objects.
    """
    h = hashlib.sha256()
    for path in patched:
        h.update(path.name.encode())
        h.update(path.read_bytes())
    return h.hexdigest()[:8]


def tag(res: str, levels: int, ranks: int, profile: str,
        frame_pointers: bool, extra: list[str], source: str = "") -> str:
    """The build directory's name, and it names every input that changes a byte.

    Extra flags are hashed rather than spelled out, because they contain
    characters a directory name cannot carry -- but they ARE in the name, so two
    flag arms cannot land in one directory and quietly reuse each other's
    objects. `source` is the same thing for a patched model source, and it goes
    under `PATCHED_ROOT` rather than beside the registry's builds."""
    parts = [res.lower(), f"l{levels}", f"p{ranks}", profile]
    if frame_pointers:
        parts.append("fp")
    if extra:
        parts.append("x" + hashlib.sha256(" ".join(extra).encode()).hexdigest()[:8])
    if source:
        parts.append("s" + source)
    return "_".join(parts)


def executable_name(res: str, levels: int, ranks: int,
                    frame_pointers: bool) -> str:
    """The registry's name for this binary.

    `most_plasim_<res>_l<levels>_p<ranks>.x`, and the `p` is the parallel width
    -- threads here -- not the precision. It carries no parallel-mode suffix
    because there is one parallel mode; `vendor/exoplasim/exoplasim/__init__.py`
    composes the same name when it goes looking for the binary to run.
    """
    suffix = "_fp" if frame_pointers else ""
    return f"most_plasim_{res.lower()}_l{levels}_p{ranks}{suffix}.x"


def build(res_arg: str, levels: int, ranks: int, profile: str,
          frame_pointers: bool, jobs: int | None, verbose: bool,
          extra: list[str] | None = None, drop: list[str] | None = None,
          publish: bool = True) -> Path:
    extra = list(extra or [])
    drop = list(drop or [])
    res, nlat = resolve(res_arg)
    if ranks < 1:
        raise SystemExit(f"--ranks must be at least 1, got {ranks}")
    if nlat % ranks:
        raise SystemExit(
            f"--ranks {ranks} does not divide the {nlat} latitudes of {res}. "
            f"NLPP is NLAT/NPRO and is a parameter, so the bands would not tile "
            f"the globe.")

    flags, precision = flag_line(profile)
    # DROP FIRST, so an arm that removes a declared flag and adds a replacement
    # gets both. A flag named here that is not in the declared line is an error
    # rather than a no-op: an arm that silently dropped nothing would be
    # measuring the control and reporting it as the treatment.
    for d in drop:
        if d not in flags:
            raise SystemExit(
                f"--drop-flag {d!r} is not in the declared line for profile "
                f"{profile!r}: {' '.join(flags)}")
        flags = [f for f in flags if f != d]
    if frame_pointers:
        flags = flags + ["-fno-omit-frame-pointer"]
    flags = flags + ["-fopenmp", "-DOMPSHARED"]
    # LAST, so an arm can override a declared flag by restating it. This is the
    # only route by which a flag reaches the compiler without being declared in
    # config/planet.yaml, and it exists for the verification arms that vary a
    # flag ON PURPOSE -- verify_fold_exactness.sh's -ffp-contract=off, and the
    # flag sweep. An arm's flags are part of its build directory's name.
    flags = flags + extra

    # WHICH FFT, from the ladder registry and not from a set named here: the
    # rule is which radix chain covers this rung's longitude count, and it is
    # derived in `lib/rungs.py:fft_module` rather than listed. A hand-kept set
    # was how the postprocessor came to route T63 and T106 into the module
    # that cannot transform them. world-i38.
    fft = rungs.fft_module(res)
    compiler = "gfortran"
    if shutil.which(compiler) is None:
        raise SystemExit(f"{compiler} is not on PATH, and the model is built with it.")

    # A PATCHED SOURCE IS AN ARM, and an arm is not a registry entry. It gets
    # its own root and its own hash, and it never publishes -- the registry's
    # tag has to keep meaning "built from the committed model source", and a
    # caller that reads --print-path gets the arm it asked for either way.
    patched = patched_sources()
    root, source = BUILD_ROOT, ""
    if patched:
        names = ", ".join(f.name for f in patched)
        if publish:
            # LOUD RATHER THAN DIVERTED, because a caller that asked to publish
            # is a registry caller: `rebuild_binaries.py` builds the whole
            # matrix this way, and a control patch left behind by an
            # interrupted verification would otherwise be what the registry
            # records the sha of.
            raise SystemExit(
                f"the model source carries a control patch marker in {names}, "
                f"so this build would not be the model. An arm built from a "
                f"patched source cannot publish under the registry's naming: "
                f"pass --no-publish, or restore the source with "
                f"`git checkout -- vendor/exoplasim/exoplasim/plasim/src`.")
        root = PATCHED_ROOT
        source = source_state(patched)
        print(f"model source carries a control patch marker in {names}: "
              f"building under {PATCHED_ROOT.name}/ so the registry's tag keeps "
              f"meaning the committed source.", file=sys.stderr)

    bdir = root / tag(res, levels, ranks, profile, frame_pointers,
                      drop + extra, source)
    bdir.mkdir(parents=True, exist_ok=True)

    configure = [
        "cmake", "-G", "Ninja", "-S", str(PLASIM), "-B", str(bdir),
        f"-DCMAKE_Fortran_COMPILER={compiler}",
        f"-DPLASIM_NLAT={nlat}",
        f"-DPLASIM_NLEV={levels}",
        f"-DPLASIM_NPRO={ranks}",
        f"-DPLASIM_PRECISION={precision}",
        f"-DPLASIM_FFT={fft}",
        f"-DPLASIM_FFLAGS={' '.join(flags)}",
        f"-DPLASIM_SHTNS_PREFIX={SHTNS_PREFIX}",
    ]
    run(configure, bdir, verbose)

    compile_cmd = ["cmake", "--build", str(bdir)]
    if jobs:
        compile_cmd += ["--parallel", str(jobs)]
    run(compile_cmd, bdir, verbose)

    built = bdir / "plasim.x"
    if not built.is_file():
        raise SystemExit(f"the build reported success but {built} is not there")

    # AN ARM IS NOT A REGISTRY ENTRY. A caller building a comparison arm --
    # verify_shtns_model's two, the flag sweep's, a thread count outside MATRIX
    # -- wants the executable, not the registry's copy of it. Publishing anyway
    # overwrites the shipped binary with an arm and leaves thread counts nobody
    # asked for beside it, so `check_consistency` reports the executable a run
    # would pick up as having unknown provenance. That is rule 4's failure mode
    # reached from inside a verification. world-v3d.
    if not publish:
        return built

    # THE NAME IS WRITTEN LAST AND ONLY ON SUCCESS, and the old file is removed
    # first. The build this replaced removed only the name it was about to write,
    # so a build that silently became another resolution left the caller's
    # expected name pointing at a binary from a previous session -- which is how
    # three arms of a flag comparison turned out to be one file.
    MODEL_RUN.mkdir(parents=True, exist_ok=True)
    out = MODEL_RUN / executable_name(res, levels, ranks, frame_pointers)
    out.unlink(missing_ok=True)
    shutil.copy2(built, out)
    return out


def run(cmd: list[str], cwd: Path, verbose: bool) -> None:
    r = subprocess.run(cmd, cwd=cwd, capture_output=not verbose, text=True)
    if r.returncode != 0:
        if not verbose:
            sys.stderr.write(r.stdout or "")
            sys.stderr.write(r.stderr or "")
        raise SystemExit(f"failed: {' '.join(cmd)}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Build one PlaSim executable. Every argument is checked; "
                    "an unrecognised value is an error, never a quiet fallback.")
    ap.add_argument("--res", required=True,
                    help="resolution name (T21..T170) or latitude count (32..256)")
    ap.add_argument("--levels", type=int, default=10, help="vertical levels (default 10)")
    ap.add_argument("--ranks", type=int, required=True,
                    help="OpenMP threads; NLAT must divide by it, and it is the "
                         "`p` in the executable's name")
    ap.add_argument("--profile", default=None,
                    help="build profile from config/planet.yaml "
                         "(default: model.compile_flags.profile)")
    ap.add_argument("--frame-pointers", action="store_true",
                    help="a PROFILING build: -fno-omit-frame-pointer, named _fp. "
                         "Costs -1.17%% at T170 with an identical restart sha, so a "
                         "profile taken this way measures the same model.")
    ap.add_argument("--extra-flag", action="append", default=[], metavar="FLAG",
                    help="an extra compiler flag, appended after the declared ones. "
                         "For verification arms that vary a flag deliberately; it is "
                         "part of the build directory's identity. Repeatable.")
    ap.add_argument("--drop-flag", action="append", default=[], metavar="FLAG",
                    help="remove a flag from the declared line. Errors if it is not "
                         "there, so an arm cannot quietly drop nothing. Repeatable.")
    ap.add_argument("--jobs", type=int, default=None,
                    help="parallel compile jobs (default: ninja's own choice)")
    ap.add_argument("--print-path", action="store_true",
                    help="print only the path of the executable produced")
    ap.add_argument("--no-publish", action="store_true",
                    help="leave the executable in its build directory instead "
                         "of copying it into the model run directory under the "
                         "registry's naming. For a comparison ARM: an arm that "
                         "publishes overwrites the shipped binary and gives "
                         "the next run unknown provenance")
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()

    profile = a.profile or str(declared()["compile_flags"]["profile"])
    for f in a.extra_flag + a.drop_flag:
        if not f.startswith("-") or f.split() != [f]:
            raise SystemExit(f"{f!r} must be one token starting with '-'")
    out = build(a.res, a.levels, a.ranks, profile,
                a.frame_pointers, a.jobs, a.verbose, a.extra_flag, a.drop_flag,
                publish=not a.no_publish)
    if a.print_path:
        print(out)
    else:
        print(f"built {out}")


if __name__ == "__main__":
    main()
