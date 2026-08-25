#!/usr/bin/env python3
"""Every model source the executable is built from passes the compiler's front end.

    python exoplasim/scripts/verify_model_compiles.py

Worldbuilding frame: this compiles the Vesper climate model's Fortran. Nothing
here concerns the simulated planet; it is a property of the source text.

WHY THIS EXISTS. Every other check of the model source in this project is a grep
or a text parse, so the one property that makes the source usable at all -- that
a compiler accepts it -- was checked by nothing until somebody ran
`rebuild_binaries.py` by hand. A dropped continuation marker in `landmod.f90`
therefore survived a merge and five parallel batches with `smoke_test` reporting
every check green, and no binary could be built from the tip at all. world-adts.

WHAT IT IS NOT. It does not build a binary and it must not: whether an
executable is CURRENT is rule 4's question and `check_consistency.py` owns it.
This asks the separate question of whether the source in the tree compiles, and
those are different failures with different fixes. `-fsyntax-only` runs the
front end and stops: no optimisation, no code generation, no link.

THE FLAGS ARE THE DECLARED ONES, and that is load-bearing rather than tidy.
Without `-fdefault-real-8` the kind selectors in `plasimmod.f90` report a
REAL(4)-to-REAL(8) mismatch that the real build does not have, so a syntax gate
run on a flag line of its own INVENTS defects. `config/planet.yaml` declares the
line, `build_model.py:flag_line` reads it, CMake appends `-fdefault-real-8` for
`precision_bytes: 8`, and this reproduces exactly that composition rather than
restating any part of it.

THE SOURCE SET AND ITS ORDER. `restart_schema.compiled_modules()` already parses
the `_sources` list out of `plasim/CMakeLists.txt`, so the set of translation
units is read from the build file and not listed here. That list is the LINK
order, which is not a compile order -- `plasim.f90` precedes the `plasimmod.f90`
it uses, and under Ninja the edges come from scanning `use` statements. This
scans them the same way and topologically sorts, seeded by the link order so the
result stays as close to it as the dependencies allow.

TWO THINGS THE BUILD FILE SUPPLIES THAT THE SOURCE TREE DOES NOT. `resmod.f90`
is generated from the requested geometry and lives only in a build directory, so
it is written here from the CONFIGURED geometry -- this gate does not iterate the
resolution ladder and has no opinion about which rungs exist. And `shtnsmod.f90`
includes `shtns.f03` from the SHTns prefix, so the prefix must be present; a
missing one is reported rather than skipped, for the same reason CMake refuses
it.

THE GATE OWNS ITS WORKING DIRECTORY, and that is a correctness property rather
than tidiness. gfortran searches the CURRENT DIRECTORY for module files ahead of
both `-I` and `-J`, so a single `.mod` left lying where the caller happened to be
standing replaces one of the modules this chain just wrote -- and `*.mod` is
gitignored, so nothing shows it is there. The report is not "stale module": it is
a shape mismatch in some third file at a line that is entirely correct, which
reads as a defect in the model. `run_front_end` has the case that cost a day.

`PLASIM_FFT` selects between two sources that both declare `module fftmod`, so
they cannot share a module directory. The configured rung's module goes in the
chain and the other variant is checked standalone afterwards, because a defect
in the variant this rung does not select is still a defect in the tree.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import CONFIG, MODEL_SRC, PROJECT_ROOT  # noqa: E402
import rungs  # noqa: E402  -- the ladder registry; _paths put lib on the path
import build_model  # noqa: E402  -- the declared flag line, read once
import restart_schema  # noqa: E402  -- the _sources parser, read once

PLASIM = MODEL_SRC / "plasim"
SRC = PLASIM / "src"
SHTNS_PREFIX = build_model.SHTNS_PREFIX

# Modules the compiler supplies. `use omp_lib` needs -fopenmp, which the flag
# line already carries; the rest need nothing and have no defining source here.
INTRINSIC = {"iso_c_binding", "iso_fortran_env", "omp_lib", "omp_lib_kinds",
             "ieee_arithmetic", "ieee_exceptions", "ieee_features"}

_MODULE = re.compile(r"^\s*module\s+([a-z_]\w*)\s*(?:!.*)?$", re.IGNORECASE)
_USE = re.compile(r"^\s*use\s*(?:::)?\s*([a-z_]\w*)", re.IGNORECASE)


def declared_flags(profile: str | None) -> tuple[list[str], str]:
    """The flag line a build of this profile hands the compiler.

    Composed exactly as `build_model.build` composes it: the declared
    `f90_opts`, then the profile's additions, then the threading flags, then
    `-fdefault-real-8` where CMake appends it for `precision_bytes: 8`.
    """
    model = build_model.declared()
    profile = profile or str(model["compile_flags"]["profile"])
    flags, precision = build_model.flag_line(profile)
    flags = flags + ["-fopenmp", "-DOMPSHARED"]
    if precision == 8:
        flags = flags + ["-fdefault-real-8"]
    return flags, profile


def geometry() -> tuple[str, int, int, int]:
    """The configured rung and its geometry. This gate iterates nothing."""
    import yaml
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    rung, nlat, _nlon = rungs.model_grid(config)
    return rung, nlat, int(config["model"]["layers"]), int(config["model"]["ncpus"])


def write_resmod(into: Path, nlat: int, nlev: int, npro: int) -> Path:
    """The generated geometry module, as `plasim/CMakeLists.txt` writes it."""
    path = into / "resmod.f90"
    path.write_text(
        "      module resmod ! generated by verify_model_compiles.py\n"
        f"      parameter(NLAT_ATM = {nlat})\n"
        f"      parameter(NLEV_ATM = {nlev})\n"
        f"      parameter(NPRO_ATM = {npro})\n"
        "      end module resmod\n\n", encoding="utf-8")
    return path


def scan(path: Path) -> tuple[set[str], set[str]]:
    """The modules a source DEFINES and the ones it USES, lowercased."""
    defines, uses = set(), set()
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = _MODULE.match(line)
        if m and m[1].lower() != "procedure":
            defines.add(m[1].lower())
            continue
        m = _USE.match(line)
        if m:
            uses.add(m[1].lower())
    return defines, uses


def compile_order(paths: list[Path]) -> list[Path]:
    """The link order the build file gives, resorted so every `use` is defined first.

    Seeded by the given order, so the result is the link order wherever the
    dependencies permit it and differs only where they do not.
    """
    defines, uses, owner = {}, {}, {}
    for p in paths:
        d, u = scan(p)
        defines[p], uses[p] = d, u
        for name in d:
            if name in owner:
                # TWO FILES, ONE MODULE NAME. The second `.mod` written
                # overwrites the first, so which one a later unit binds depends
                # on compile order rather than on anything the build declares --
                # and the diagnostic that comes out is a shape or type mismatch
                # in a THIRD file, which reads as a defect in the model. The
                # `fftmod` pair is the known case and `sources()` splits it out
                # before the chain is built; anything else reaching here is a
                # question about the build file, not something to pick a winner
                # for.
                raise SystemExit(
                    f"{owner[name].name} and {p.name} both declare `module "
                    f"{name}`, and both are in the compiled set. One of them "
                    f"would overwrite the other's module file, so there is no "
                    f"chain to check. Split them the way PLASIM_FFT's pair is "
                    f"split in `sources()`.")
            owner[name] = p
    pending, done, out = list(paths), set(), []
    while pending:
        for p in pending:
            needed = {owner[n] for n in uses[p]
                      if n in owner and n not in INTRINSIC}
            if needed <= done | {p}:
                out.append(p)
                done.add(p)
                pending.remove(p)
                break
        else:
            names = ", ".join(p.name for p in pending)
            raise SystemExit(
                f"the `use` graph over {names} has a cycle, so there is no "
                f"compile order. Ninja would fail on the same graph.")
    return out


def sources(rung: str) -> tuple[list[str], list[str]]:
    """The compiled `.f90` set, split into the configured chain and the variants.

    `compiled_modules` expands every value a configurable slot admits, because
    a restart record must be the same across configurations. A COMPILE cannot
    take both `fftmod.f90` and `fft991mod.f90`: they declare the same module.
    """
    every = list(restart_schema.compiled_modules(PLASIM))
    fft = rungs.fft_module(rung) + ".f90"
    variants = [f for f in every if f.startswith("fft") and f != fft]
    return [f for f in every if f not in variants], variants


def run_front_end(path: Path, flags: list[str], moddir: Path,
                  search: list[Path], cwd: Path, verbose: bool) -> str:
    """One translation unit. `search` is the module path IN ORDER, `cwd` is owned.

    `flags` is the flag LINE, a flat list of strings -- not the `(flags,
    profile)` pair `declared_flags` returns. Checked rather than concatenated,
    because the pair is what a caller reaches for and the failure it produces
    otherwise is a TypeError from inside `subprocess`.

    THE WORKING DIRECTORY IS THE GATE'S OWN, and that is the whole point of the
    argument. gfortran searches the CURRENT DIRECTORY for module files BEFORE
    both `-I` and `-J`, so one stale `.mod` lying in whatever directory the
    caller happened to be in silently replaces one of the modules this chain
    just wrote. What comes out is not "stale module": it is a shape mismatch in
    a file that uses it, at a line that is perfectly correct, and it reads
    exactly like a defect in the model. A T42 `glaciermod.mod` in the caller's
    directory made this gate report `outmod.f90:979` as a 512-against-128
    assignment on a tree that compiles clean. A gate that inherits the cwd
    invents defects, so this one does not inherit it.
    """
    bad = [f for f in flags if not isinstance(f, str)]
    if bad:
        raise TypeError(
            f"run_front_end wants the flag line, a list of strings, and got "
            f"{bad[0]!r} in it. `declared_flags` returns (flags, profile): "
            f"pass the first element.")
    cmd = (["gfortran", "-fsyntax-only", "-J", str(moddir)]
           + [f"-I{d}" for d in search] + list(flags) + [str(path)])
    if verbose:
        print("  " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    return "" if r.returncode == 0 else (r.stderr.strip() or "(no diagnostic)")


def verify(profile: str | None, verbose: bool) -> list[str]:
    """Every translation unit, in compile order. Returns the failures."""
    if shutil.which("gfortran") is None:
        return ["gfortran is not on PATH, and the model is built with it"]
    include = SHTNS_PREFIX / "include" / "shtns.f03"
    if not include.is_file():
        return [f"no shtns.f03 under {SHTNS_PREFIX}/include, and shtnsmod.f90 "
                f"includes it. Run exoplasim/scripts/build_shtns.sh."]

    flags, profile = declared_flags(profile)
    rung, nlat, nlev, npro = geometry()
    chain, variants = sources(rung)
    problems = []
    with tempfile.TemporaryDirectory(prefix="plasim-syntax-") as tmp:
        tmp = Path(tmp)
        moddir = tmp / "modules"
        moddir.mkdir()
        # Empty, and it stays empty: it is what every invocation runs in, so
        # that no `.mod` outside this directory tree can be reached at all.
        work = tmp / "cwd"
        work.mkdir()
        resmod = write_resmod(tmp, nlat, nlev, npro)
        paths = compile_order([resmod] + [SRC / f for f in chain])
        # The chain's own modules FIRST. The SHTNs prefix is here for
        # `shtnsmod.f90`'s `include 'shtns.f03'` and holds no module files.
        includes = [moddir, SHTNS_PREFIX / "include", tmp]
        print(f"{len(paths) + len(variants)} translation units, "
              f"{rung} l{nlev} p{npro} profile {profile}")
        print(f"gfortran -fsyntax-only {' '.join(flags)}")
        for path in paths:
            err = run_front_end(path, flags, moddir, includes, work, verbose)
            if err:
                problems.append(f"{path.name} does not compile:\n{err}")
                # The rest of the chain would fail for want of this one's
                # module file, and those are not findings.
                break
        # A variant redeclares a module the chain already wrote, so it gets a
        # module directory of its own, named FIRST so its own copy wins over the
        # chain's. That ordering is not cosmetic: fftmod.f90's `fftmod` has no
        # `ifax`, so with the chain's directory ahead of it fft991mod.f90's own
        # subroutines were checked against the other variant's module and
        # reported three rank mismatches that no build has.
        for name in variants:
            vdir = tmp / f"modules-{name}"
            vdir.mkdir()
            err = run_front_end(SRC / name, flags, vdir,
                                [vdir] + includes, work, verbose)
            if err:
                problems.append(f"{name} does not compile:\n{err}")
    return problems


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--profile", default=None,
                    help="build profile from config/planet.yaml "
                         "(default: model.compile_flags.profile)")
    ap.add_argument("--verbose", action="store_true",
                    help="print each compiler invocation")
    a = ap.parse_args()
    problems = verify(a.profile, a.verbose)
    for p in problems:
        print(p, file=sys.stderr)
    print("the model source does not compile" if problems
          else "the model source compiles")
    raise SystemExit(1 if problems else 0)


if __name__ == "__main__":
    main()
