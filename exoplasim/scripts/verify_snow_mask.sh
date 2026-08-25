#!/bin/bash
# Does the canopy/snow mask reduce to the scalar it replaced?
#
#   exoplasim/scripts/verify_snow_mask.sh
#
# Worldbuilding frame: a correctness check on the surface albedo of a simulated
# planet's snow-covered vegetated land.
#
# snowmaskmod reads no model state, so this compiles it alone against the
# checking program and needs no grid, no namelist and no build tree. The flag
# line matches the model's real default, -fdefault-real-8, because the
# reduction it checks is a bitwise one and a different working precision is a
# different check.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC="$ROOT/vendor/exoplasim/exoplasim/plasim/src"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
W="$(mktemp -d)"
trap 'rm -rf "$W"' EXIT

gfortran -O2 -fdefault-real-8 -J"$W" -o "$W/verify_snow_mask" \
    "$SRC/snowmaskmod.f90" "$HERE/verify_snow_mask.f90"

"$W/verify_snow_mask"
