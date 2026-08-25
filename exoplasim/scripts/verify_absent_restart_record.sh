#!/bin/bash
# Is an absent restart record reported to the caller, or scattered into it?
#
#   exoplasim/scripts/verify_absent_restart_record.sh
#
# Worldbuilding frame: a correctness check on the restart reader of the Vesper
# climate model. Nothing here is about the simulated planet.
#
# restartmod reads no model state, so this compiles it alone against the
# checking program and needs no grid, no namelist and no build tree.
#
# -finit-real=snan IS THE INSTRUMENT and it is what makes the control arm.
# The defect is a buffer that is never written and is copied out anyway; at the
# shipped flag line that buffer holds ordinary stack, which looks like a number
# and cannot be told from one. Poisoning every local real turns the same read
# into a NaN the check can see. -Og and not -O2 for the same reason
# verify_uninitialised_reads.sh gives for its trapping arm: at -O2 gfortran can
# fold a signalling NaN away before it is read. No -ffpe-trap here, because the
# check WANTS the NaN to propagate into the caller's array and be tested there.
#
# -fdefault-real-8 matches the model's real default: the restart records are
# written and read at the model's working precision and a different one is a
# different check.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC="$ROOT/vendor/exoplasim/exoplasim/plasim/src"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
W="$(mktemp -d)"
trap 'rm -rf "$W"' EXIT

gfortran -Og -fdefault-real-8 -finit-real=snan -J"$W" \
    -o "$W/verify_absent_restart_record" \
    "$SRC/restartmod.f90" "$HERE/verify_absent_restart_record.f90"

cd "$W"
"$W/verify_absent_restart_record"
