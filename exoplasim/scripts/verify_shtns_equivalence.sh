#!/bin/bash
# Build and run the SHTns equivalence gate.
#
#   exoplasim/scripts/verify_shtns_equivalence.sh [res] [shtns-prefix]
#
# The recipe and the reasoning are in the driver's header. This compiles it
# against the model's OWN modules at one process, so the globe is local and the
# comparison is of transforms rather than of decompositions, and runs both arms
# from one binary: the real one, and the control with the Condon-Shortley phase
# dropped, which must fail.
#
# SHTns is not a system package here. Point this at a prefix built with
# --enable-openmp; the second argument, or SHTNS_PREFIX, or the default below.
set -euo pipefail

res="${1:-T42}"
PREFIX="${2:-${SHTNS_PREFIX:-/usr/local}}"

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
SRC="$REPO/vendor/exoplasim/exoplasim/plasim/src"
W="${TMPDIR:-/tmp}/verify_shtns_equiv.$$"
trap 'rm -rf "$W"' EXIT
mkdir -p "$W"; cd "$W"

case "$res" in
  T21) nlat=32 ;; T31) nlat=48 ;; T42) nlat=64 ;;
  T85) nlat=128 ;; T127) nlat=192 ;; T170) nlat=256 ;;
  *) echo "unknown resolution $res" >&2; exit 2 ;;
esac

if [ ! -f "$PREFIX/include/shtns.f03" ]; then
    echo "no shtns.f03 under $PREFIX -- build SHTns with --enable-openmp and" >&2
    echo "pass its prefix, or set SHTNS_PREFIX" >&2
    exit 2
fi

cp "$SRC"/plasimmod.f90 "$SRC"/legmod.f90 "$SRC"/fftmod.f90 \
   "$SRC"/gaussmod.f90 "$SRC"/specblock.f90 .
cp "$HERE"/verify_shtns_equivalence.f90 drive.f90

cat > resmod.f90 <<EOF
      module resmod
      parameter(NLAT_ATM = $nlat)
      parameter(NLEV_ATM = 10)
      parameter(NPRO_ATM = 1)
      end module resmod
EOF

F="-c -O2 -cpp -ffixed-line-length-132 -ffpe-summary=none -finit-real=zero -fdefault-real-8"
for f in resmod plasimmod gaussmod specblock fftmod legmod; do
    gfortran $F -J . "$f.f90" -o "$f.o" 2>"$f.err" \
        || { echo "compile failed: $f"; head -12 "$f.err"; exit 1; }
done
gfortran -O2 -fdefault-real-8 -ffixed-line-length-132 -J . -I"$PREFIX/include" \
    -c drive.f90 -o drive.o 2>drive.err \
    || { echo "compile failed: drive"; head -25 drive.err; exit 1; }
gfortran -o equiv.x drive.o legmod.o fftmod.o gaussmod.o specblock.o \
    plasimmod.o resmod.o "$PREFIX"/lib/libshtns*.a -lfftw3_omp -lfftw3 -lm -fopenmp \
    2>link.err || { echo "link failed"; head -12 link.err; exit 1; }

echo "==== SHTns against legmod, $res ===="
./equiv.x | grep -vE "^ \*|^$"

echo
echo "==== the control, Condon-Shortley dropped: must FAIL ===="
SHTNS_DROP_CS=1 ./equiv.x | grep -vE "^ \*|^$" || true
if SHTNS_DROP_CS=1 ./equiv.x >/dev/null 2>&1; then
    echo "  [ FAIL ] the control passed"
    exit 1
fi
