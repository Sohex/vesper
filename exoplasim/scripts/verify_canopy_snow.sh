#!/bin/bash
# Does the canopy snow store conserve what falls into it, and reduce to no
# store where the model has no canopy?
#
#   exoplasim/scripts/verify_canopy_snow.sh
#
# Worldbuilding frame: a correctness check on the snow held in the canopy of a
# simulated planet's vegetated land.
#
# cansnowmod reads no model state, so this compiles it alone against the
# checking program and needs no grid, no namelist and no build tree. The flag
# line matches the model's real default, -fdefault-real-8, because one of the
# checks is a bitwise reduction and a different working precision is a
# different check.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC="$ROOT/vendor/exoplasim/exoplasim/plasim/src"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
W="$(mktemp -d)"
trap 'rm -rf "$W"' EXIT

gfortran -O2 -fdefault-real-8 -J"$W" -o "$W/verify_canopy_snow" \
    "$SRC/cansnowmod.f90" "$HERE/verify_canopy_snow.f90"

"$W/verify_canopy_snow"
