#!/bin/bash
# Build SHTns for this model, into a project-local prefix.
#
#   exoplasim/scripts/build_shtns.sh [--force]
#
# Worldbuilding frame: a build dependency of the Vesper climate model. Nothing
# here is about the simulated planet.
#
# WHY NOT A SUBTREE, and why not the system. `vendored-upstreams.md` describes
# two git subtrees, both of which this project EDITS -- the Orogen fork carries
# our generator changes and the ExoPlaSim fork is compiled in place. SHTns is
# neither: it is used unmodified, and vendoring an unmodified upstream buys
# nothing but a large tree to keep in step. Nor is it installed system-wide,
# because that puts a build dependency outside the repo where no manifest
# records it and no other machine reproduces it.
#
# So it is built from a PINNED revision into `vendor/shtns-install/`, untracked
# the way `.venv` is untracked, by a script that is registered the way
# `requirements.txt` is. The revision is here, in one place, and moving it is a
# deliberate edit.
#
# -march=znver4 matches most_compiler_omp so the library and the model agree on
# the instruction set; it is what gives SHTns its AVX-512, which is most of why
# it beats a hand-written transform. --enable-openmp is required even though the
# model calls it SINGLE-threaded, because that is the build that ships the
# thread-safe internals the parallel-over-field-and-level design relies on.
set -euo pipefail

# PINNED. Not master: a build dependency that moves under you is a build you
# cannot reproduce, and every convention this model reconciles against -- the
# Condon-Shortley phase, the Robert form, the sign of the toroidal potential --
# is a property of a REVISION. Moving this is a deliberate edit, and
# verify_shtns_equivalence.sh is what says the new one still agrees.
REV="${SHTNS_REV:-4e69ceb}"          # SHTns 3.7.5
URL="https://bitbucket.org/nschaeff/shtns.git"

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
PREFIX="$REPO/vendor/shtns-install"
SRC="$REPO/vendor/shtns-src"

if [ "${1:-}" != "--force" ] && [ -f "$PREFIX/include/shtns.f03" ]; then
    echo "already built at $PREFIX"
    echo "  $(ls "$PREFIX"/lib/libshtns*.a 2>/dev/null | head -1)"
    echo "  pass --force to rebuild"
    exit 0
fi

if [ ! -d "$SRC/.git" ]; then
    echo "cloning SHTns"
    rm -rf "$SRC"
    git clone --quiet "$URL" "$SRC"
fi
( cd "$SRC" && git checkout --quiet "$REV" )
ver=$(grep -oE 'AC_INIT\(\[SHTns\],\[[0-9.]+\]' "$SRC/configure.ac" | grep -oE '[0-9][0-9.]*')
echo "SHTns $ver at $(cd "$SRC" && git rev-parse --short HEAD)"

rm -rf "$PREFIX"; mkdir -p "$PREFIX"
cd "$SRC"
make distclean >/dev/null 2>&1 || true
CFLAGS="-O3 -march=znver4" ./configure --prefix="$PREFIX" --enable-openmp \
    > "$SRC/configure.log" 2>&1 \
    || { echo "configure failed, see $SRC/configure.log" >&2; tail -20 "$SRC/configure.log"; exit 1; }
make -j"$(nproc)" > "$SRC/build.log" 2>&1 \
    || { echo "build failed, see $SRC/build.log" >&2; tail -20 "$SRC/build.log"; exit 1; }
make install >> "$SRC/build.log" 2>&1

lib="$(ls "$PREFIX"/lib/libshtns*.a | head -1)"
n512=$(objdump -d "$lib" 2>/dev/null | grep -cE "%zmm[0-9]+" || true)
echo "installed $lib"
echo "  AVX-512 instructions: $n512"
[ "${n512:-0}" -gt 0 ] || { echo "  NO AVX-512 -- the march flag did not take" >&2; exit 1; }
echo "  fortran interface: $PREFIX/include/shtns.f03"
