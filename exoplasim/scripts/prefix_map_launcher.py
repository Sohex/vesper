#!/usr/bin/env python3
"""Apply the declared -ffile-prefix-map to the cpp line markers GCC will not.

    <compiler> ... -c foo-pp.f90 -o foo.o      # invoked as a CMake launcher

Worldbuilding frame: this is part of the build of the Vesper project's climate
model, a hard fork of ExoPlaSim. Nothing here concerns the simulated planet.

WHY IT EXISTS. `gcc(1)` on -ffile-prefix-map: "Directories referenced by
directives are not affected by these options." A cpp line marker is such a
directive, and CMake's Ninja generator compiles Fortran in two stages -- it
preprocesses the source into `<build>/CMakeFiles/.../<name>-pp.f90` so it can
scan the module graph, then compiles that with -fpreprocessed. The absolute
source path CMake hands the first stage is written into the line markers, and
gfortran takes the line marker verbatim for the string it stores beside every
I/O statement so the runtime can say which file an error came from. Those
strings sit in .rodata and are not debug information: they survive without -g
and no prefix map reaches them.

WHAT THAT COST. The executable then differs byte for byte between two checkouts
of the same commit, so a build made in a git worktree is absent from
`exoplasim/binary_manifest.json` under any name and every consumer check refuses
it. A worktree could compile the model -- the one way to verify a Fortran change
-- and could not run what it built through the registry. world-ynpx.

WHAT THIS DOES, and it is the whole of it: before the second stage runs, rewrite
the line markers of the file it is about to compile by applying the
-ffile-prefix-map options that are ALREADY ON THAT COMMAND LINE. The mapping
stays declared in `plasim/CMakeLists.txt` and is not restated here; this only
extends it to the one place the compiler documents as out of scope. Nothing else
about the command line is touched, and the compiler is then exec'd unchanged.

The rewrite is idempotent: a marker that has already been mapped no longer
carries the prefix, so an incremental build that recompiles without
re-preprocessing maps nothing a second time.

Registered by `build_model.py` as CMAKE_Fortran_COMPILER_LAUNCHER. CMake places
the launcher on the compile rules only, never on the preprocess rule -- which is
why the fix is here and not a flag handed to the preprocessor.
"""
from __future__ import annotations

import os
import re
import sys

# `# <line> "<file>"` and `#line <line> "<file>"`, which is the whole of what a
# preprocessed Fortran file carries a path in.
MARKER = re.compile(rb'^([ \t]*#(?:line)?[ \t]+[0-9]+[ \t]+")([^"]*)(")', re.M)

PREPROCESSED_SUFFIXES = (".f90", ".f", ".f95", ".f03", ".f08")


def prefix_maps(args: list[str]) -> list[tuple[bytes, bytes]]:
    """The -ffile-prefix-map pairs on this command line, longest prefix first.

    Split on the LAST '=' because that is what GCC's own option parser does, so
    a map this reads is the map the compiler received and not a second reading
    of it.
    """
    out = []
    for a in args:
        if a.startswith("-ffile-prefix-map="):
            old, sep, new = a[len("-ffile-prefix-map="):].rpartition("=")
            if sep and old:
                out.append((old.encode(), new.encode()))
    out.sort(key=lambda pair: -len(pair[0]))
    return out


def compile_input(args: list[str]) -> str | None:
    """The already-preprocessed file this invocation compiles, if it is one.

    -fpreprocessed is the signal, and it is the one CMake's scanned-compile rule
    carries: without it the compiler is reading a source whose path it was
    handed on the command line, and there are no line markers to map.
    """
    if "-fpreprocessed" not in args:
        return None
    for i, a in enumerate(args):
        if a == "-c" and i + 1 < len(args):
            path = args[i + 1]
            if path.lower().endswith(PREPROCESSED_SUFFIXES) and os.path.isfile(path):
                return path
    return None


def remap(path: str, maps: list[tuple[bytes, bytes]]) -> None:
    data = open(path, "rb").read()

    def one(m: re.Match) -> bytes:
        name = m.group(2)
        for old, new in maps:
            if name.startswith(old):
                return m.group(1) + new + name[len(old):] + m.group(3)
        return m.group(0)

    mapped = MARKER.sub(one, data)
    if mapped != data:
        tmp = path + ".prefixmap"
        with open(tmp, "wb") as fh:
            fh.write(mapped)
        os.replace(tmp, path)


def main(argv: list[str]) -> None:
    if len(argv) < 2:
        raise SystemExit("usage: prefix_map_launcher.py <compiler> [args...]")
    # A launcher takes the compiler as its first argument, so no compiler can be
    # named `--help` and answering it here cannot shadow a real invocation.
    if argv[1] in ("-h", "--help"):
        print(__doc__)
        return
    compiler, args = argv[1], argv[2:]
    maps = prefix_maps(args)
    src = compile_input(args) if maps else None
    if src is not None:
        remap(src, maps)
    os.execvp(compiler, [compiler] + args)


if __name__ == "__main__":
    main(sys.argv)
