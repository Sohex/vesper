#!/usr/bin/env python3
"""Build one PlaSim executable, and refuse anything it does not understand.

    python exoplasim/scripts/build_model.py --res T170 --ranks 16 --parmode omp
    python exoplasim/scripts/build_model.py --res T21 --ranks 8 --parmode mpi --profile checked

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
`most_compiler_mpi`, `most_compiler_omp` -- of which only the MPI one was
written from the declaration, so the threaded build, which is what this project
now measures, took whatever was in a file nothing derived. That was CONS-14.

WHERE IT BUILDS. One directory per configuration under `vendor/exoplasim/build/`,
so two configurations cannot delete each other's objects. The old build had one
`plasim/bld` that it emptied on entry, which is why `rebuild_binaries.py` had to
be serial across its twelve executables and why the tree needed marker files to
notice it held the wrong configuration. Both go away.
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
SHTNS_PREFIX = PROJECT_ROOT / "vendor" / "shtns-install"

# Resolution name -> latitudes, from `lib/gridding.py`, which is the one place
# the ladder is written down and checks each rung's name against its own grid.
# A value outside it must be an ERROR: `-r 170` used to build T21 and say
# nothing. SPAT-2.
RESOLUTIONS = dict(rungs.RUNGS)
# T63 and T106 are not powers of two in longitude and need the other FFT.
NEEDS_FFT991 = {"T63", "T106"}

PARMODES = ("serial", "mpi", "omp")


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


def tag(res: str, levels: int, ranks: int, parmode: str, profile: str,
        frame_pointers: bool, extra: list[str]) -> str:
    """The build directory's name, and it names every input that changes a byte.

    Extra flags are hashed rather than spelled out, because they contain
    characters a directory name cannot carry -- but they ARE in the name, so two
    flag arms cannot land in one directory and quietly reuse each other's
    objects."""
    parts = [res.lower(), f"l{levels}", f"p{ranks}", parmode, profile]
    if frame_pointers:
        parts.append("fp")
    if extra:
        parts.append("x" + hashlib.sha256(" ".join(extra).encode()).hexdigest()[:8])
    return "_".join(parts)


def executable_name(res: str, levels: int, ranks: int, parmode: str,
                    frame_pointers: bool) -> str:
    """The registry's name for this binary. Unchanged from the old build, because
    `binary_manifest.json` and every bed manifest already carry these names."""
    suffix = ""
    if parmode == "omp":
        suffix += "_omp"
    if frame_pointers:
        suffix += "_fp"
    return f"most_plasim_{res.lower()}_l{levels}_p{ranks}{suffix}.x"


def build(res_arg: str, levels: int, ranks: int, parmode: str, profile: str,
          frame_pointers: bool, jobs: int | None, verbose: bool,
          extra: list[str] | None = None, drop: list[str] | None = None,
          publish: bool = True) -> Path:
    extra = list(extra or [])
    drop = list(drop or [])
    res, nlat = resolve(res_arg)
    if parmode not in PARMODES:
        raise SystemExit(f"--parmode {parmode!r} is not one of {' '.join(PARMODES)}")
    if ranks < 1:
        raise SystemExit(f"--ranks must be at least 1, got {ranks}")
    if nlat % ranks:
        raise SystemExit(
            f"--ranks {ranks} does not divide the {nlat} latitudes of {res}. "
            f"NLPP is NLAT/NPRO and is a parameter, so the bands would not tile "
            f"the globe.")
    if parmode == "serial" and ranks != 1:
        raise SystemExit("--parmode serial takes --ranks 1")

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
    if parmode == "omp":
        flags = flags + ["-fopenmp", "-DOMPSHARED"]
    # LAST, so an arm can override a declared flag by restating it. This is the
    # only route by which a flag reaches the compiler without being declared in
    # config/planet.yaml, and it exists for the verification arms that vary a
    # flag ON PURPOSE -- verify_fold_exactness.sh's -ffp-contract=off, and the
    # flag sweep. An arm's flags are part of its build directory's name.
    flags = flags + extra

    fft = "fft991mod" if res in NEEDS_FFT991 else "fftmod"
    compiler = "mpif90" if parmode == "mpi" else "gfortran"
    if shutil.which(compiler) is None:
        raise SystemExit(f"{compiler} is not on PATH, and --parmode {parmode} needs it.")

    bdir = BUILD_ROOT / tag(res, levels, ranks, parmode, profile, frame_pointers,
                            drop + extra)
    bdir.mkdir(parents=True, exist_ok=True)

    configure = [
        "cmake", "-G", "Ninja", "-S", str(PLASIM), "-B", str(bdir),
        f"-DCMAKE_Fortran_COMPILER={compiler}",
        f"-DPLASIM_NLAT={nlat}",
        f"-DPLASIM_NLEV={levels}",
        f"-DPLASIM_NPRO={ranks}",
        f"-DPLASIM_PRECISION={precision}",
        f"-DPLASIM_PARMODE={parmode}",
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
    out = MODEL_RUN / executable_name(res, levels, ranks, parmode, frame_pointers)
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
                    help="MPI ranks, or OpenMP threads for --parmode omp")
    ap.add_argument("--parmode", required=True, choices=PARMODES)
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
    out = build(a.res, a.levels, a.ranks, a.parmode, profile,
                a.frame_pointers, a.jobs, a.verbose, a.extra_flag, a.drop_flag,
                publish=not a.no_publish)
    if a.print_path:
        print(out)
    else:
        print(f"built {out}")


if __name__ == "__main__":
    main()
