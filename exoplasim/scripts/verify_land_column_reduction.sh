#!/bin/bash
# Build and run the land column reduction check.
#
#   exoplasim/scripts/verify_land_column_reduction.sh
#
# LSHY-3 puts the scalar soil water bucket behind a SELECTABLE reduction, and
# that word is the requirement: the bucket must stay reachable and must be the
# EXACT reduction of the column that replaces it. This compiles the model's own
# `landcolumn.f90` -- not a lifted copy -- against the checking program and
# requires bitwise equality of store and runoff over twenty thousand steps.
#
# It compiles ONE dependency-free module and one program, so it costs a
# fraction of a model build and needs no grid, no namelist and no run
# directory. That is why the kernels live in their own file.
#
# The program carries three negative controls that must break the equality. A
# control that does not break it fails the check just as a failed reduction
# does: a comparison that cannot fail proves only that both sides were run.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
SRC="$REPO/vendor/exoplasim/exoplasim/plasim/src"
W="${TMPDIR:-/tmp}/verify_land_column.$$"

FC="${FC:-gfortran}"
# No -ffast-math and no -Ofast anywhere near this. Both licence the compiler to
# reassociate the very expressions whose ORDER is what makes the reduction
# exact, and a check that passes only at -O0 is not a check on what ships.
#
# BOTH PRECISIONS, and the second is the one that matters. The model builds
# with -fdefault-real-8, so a reduction verified only at default real-4 has
# been verified at a precision nothing runs at; the single-precision arm stays
# because an exactness that depends on having spare mantissa is not exactness.
FFLAGS="${FFLAGS:--O2 -std=legacy}"

mkdir -p "$W"
trap 'rm -rf "$W"' EXIT

cp "$SRC/landcolumn.f90" "$W/"
cp "$HERE/verify_land_column_reduction.f90" "$W/"

cd "$W"
for precision in "" "-fdefault-real-8"; do
  label="${precision:-default real}"
  echo "=== $label ==="
  rm -f ./*.o ./*.mod ./verify
  $FC $FFLAGS $precision -c landcolumn.f90
  $FC $FFLAGS $precision -o verify verify_land_column_reduction.f90 landcolumn.o
  ./verify
  echo ""
done
