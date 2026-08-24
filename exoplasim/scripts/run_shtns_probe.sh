#!/bin/bash
# Build and run one of the SHTns convention probes.
#
#   exoplasim/scripts/run_shtns_probe.sh <probe.f90> [res] [shtns-prefix]
#
# Worldbuilding frame: a measurement on the Vesper climate model's spectral
# transform. Nothing here is about the simulated planet.
#
# WHY THIS EXISTS. The four probe drivers are registered source files that were
# each built by hand when they were written, which makes them findable but not
# runnable: the recipe for standing legmod up outside the model -- the resmod
# it needs generated, the abort stub that keeps the restart I/O chain out, the
# compile order -- lived in a shell history. A probe that cannot be re-run is a
# claim rather than a measurement.
#
# A PROBE IS NOT A CHECK. These print ratios and have no verdict; they exist to
# DISCOVER a convention. What enforces it afterwards is
# verify_shtns_equivalence.sh, which compares against legmod with controls, and
# verify_shtns_model.sh, which runs the model. Do not read a probe's output as
# a pass.
#
# The build matches the model's: double precision and OMPSHARED
# defined. A probe compiled at real*4 disagrees with legmod by 1e-5 and worse
# at high total wavenumber -- that is the weight factorisation losing precision
# in legmod, not a convention -- and a probe built without OMPSHARED is refused
# by shtns_setup outright.
set -euo pipefail

drv="${1:?usage: run_shtns_probe.sh <probe.f90> [res] [shtns-prefix]}"
res="${2:-T42}"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
PREFIX="${3:-${SHTNS_PREFIX:-$REPO/vendor/shtns-install}}"
SRC="$REPO/vendor/exoplasim/exoplasim/plasim/src"
WORK="$REPO/exoplasim/bench/_shtnsprobe"

[ -f "$drv" ] || drv="$HERE/$drv"
[ -f "$drv" ] || { echo "no such probe: $1" >&2; exit 2; }

# The ladder is `lib/rungs.py` and nowhere else; this asks it rather than
# carrying a `case` that has already drifted twice -- the two probes disagreed
# about whether T31 and T63 exist. SPAT-2.
PY="$REPO/.venv/bin/python"; [ -x "$PY" ] || PY=python3
nlat=$("$PY" -c "import sys; sys.path.insert(0, '$REPO/lib'); import rungs; print(rungs.geometry('$res')[0])") || {
  echo "unknown resolution $res: it is not a rung in lib/rungs.py" >&2; exit 2; }

if [ ! -f "$PREFIX/include/shtns.f03" ]; then
    echo "no shtns.f03 under $PREFIX -- exoplasim/scripts/build_shtns.sh" >&2
    echo "installs it, or set SHTNS_PREFIX" >&2
    exit 2
fi

rm -rf "$WORK"; mkdir -p "$WORK"
cd "$WORK"
cp "$SRC"/plasimmod.f90 "$SRC"/legmod.f90 "$SRC"/fftmod.f90 \
   "$SRC"/gaussmod.f90 "$SRC"/specblock.f90 "$SRC"/shtnsmod.f90 .
cp "$drv" drive.f90

# shtnsmod calls mpabort on a configuration it refuses. That is the whole of its
# dependency on the MPI layer, and mpimod_stub drags the restart I/O chain in
# behind it, so a two-line stub stands in for the model's file handling.
cat > abortstub.f90 <<'EOF'
      subroutine mpabort(ytext)
      character(len=*) :: ytext
      write(*,*) 'mpabort: ', ytext
      stop 1
      end subroutine mpabort
EOF

cat > resmod.f90 <<EOF
      module resmod
      parameter(NLAT_ATM = $nlat)
      parameter(NLEV_ATM = 10)
      parameter(NPRO_ATM = 1)
      end module resmod
EOF

F="-c -O2 -cpp -DOMPSHARED -ffixed-line-length-132 -ffpe-summary=none -finit-real=zero -fdefault-real-8"
for f in resmod plasimmod abortstub gaussmod specblock fftmod legmod; do
    gfortran $F -J . "$f.f90" -o "$f.o" 2>"$f.err" \
        || { echo "compile failed: $f"; head -15 "$f.err"; exit 1; }
done
gfortran $F -J . -I"$PREFIX/include" shtnsmod.f90 -o shtnsmod.o 2>shtnsmod.err \
    || { echo "compile failed: shtnsmod"; head -20 shtnsmod.err; exit 1; }
gfortran -O2 -cpp -DOMPSHARED -fdefault-real-8 \
    -ffixed-line-length-132 -J . -I"$PREFIX/include" -c drive.f90 -o drive.o \
    2>drive.err || { echo "compile failed: $(basename "$drv")"; head -30 drive.err; exit 1; }
gfortran -o probe.x drive.o shtnsmod.o legmod.o fftmod.o gaussmod.o \
    specblock.o plasimmod.o abortstub.o resmod.o \
    "$PREFIX"/lib/libshtns*.a -lfftw3_omp -lfftw3 -lm -fopenmp \
    2>link.err || { echo "link failed"; head -15 link.err; exit 1; }

./probe.x
