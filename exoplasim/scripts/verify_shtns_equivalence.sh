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
# SHTns lives in vendor/shtns-install, built from a pinned revision by
# exoplasim/scripts/build_shtns.sh, which is where this looks by default. A
# second argument or SHTNS_PREFIX overrides it, which is how a candidate upgrade
# is checked before it is pinned.
set -euo pipefail

res="${1:-T42}"
PREFIX="${2:-${SHTNS_PREFIX:-$(cd "$(dirname "$0")/../.." && pwd)/vendor/shtns-install}}"

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
   "$SRC"/gaussmod.f90 "$SRC"/specblock.f90 "$SRC"/shtnsmod.f90 \
   .
cp "$HERE"/verify_shtns_equivalence.f90 drive.f90

# shtnsmod calls mpabort on a configuration it refuses. That is the whole of
# its dependency on the MPI layer, and mpimod_stub drags the restart I/O chain
# in behind it, so the transform test gets a two-line abort instead of the
# model's file handling.
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

F="-c -O2 -cpp -DOMPSHARED -DNOPAIRLAT -ffixed-line-length-132 -ffpe-summary=none -finit-real=zero -fdefault-real-8"

build () {
    local tag="$1"
    rm -f ./*.o ./*.mod "equiv_$tag.x"
    for f in resmod plasimmod abortstub gaussmod specblock fftmod legmod; do
        gfortran $F -J . "$f.f90" -o "$f.o" 2>"$f.err" \
            || { echo "compile failed: $f ($tag)"; head -12 "$f.err"; exit 1; }
    done
    gfortran $F -J . -I"$PREFIX/include" shtnsmod.f90 -o shtnsmod.o 2>shtnsmod.err \
        || { echo "compile failed: shtnsmod ($tag)"; head -20 shtnsmod.err; exit 1; }
    gfortran -O2 -cpp -DOMPSHARED -DNOPAIRLAT -fdefault-real-8 -ffixed-line-length-132 -J . -I"$PREFIX/include" \
        -c drive.f90 -o drive.o 2>drive.err \
        || { echo "compile failed: drive ($tag)"; head -25 drive.err; exit 1; }
    gfortran -o "equiv_$tag.x" drive.o shtnsmod.o legmod.o fftmod.o gaussmod.o \
        specblock.o plasimmod.o abortstub.o resmod.o \
        "$PREFIX"/lib/libshtns*.a -lfftw3_omp -lfftw3 -lm -fopenmp \
        2>link.err || { echo "link failed ($tag)"; head -12 link.err; exit 1; }
}

echo "==== the shipped wrappers against legmod, $res ===="
build good
./equiv_good.x | grep -vE "^ \*|^$"

echo
echo "==== the control, Condon-Shortley dropped: must FAIL ===="
# Patched into shtnsmod ITSELF, because that is the code that ships. Dropping
# the phase negates every odd zonal wavenumber, which a comparison of
# magnitudes cannot see -- which is how it survived once already.
sed -i 's/knorm = SHT_ORTHONORMAL  /knorm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE  /' shtnsmod.f90
grep -q "SHT_NO_CS_PHASE" shtnsmod.f90 || { echo "control patch missed"; exit 1; }
build bad
if ./equiv_bad.x >/dev/null 2>&1; then
    echo "  [ FAIL ] the control passed, so this check cannot see a dropped phase"
    exit 1
else
    echo "  [  ok  ] control rejected, so the check has teeth"
fi

echo
echo "==== the control, spectral filter dropped: must FAIL ===="
# The failure this check was blind to for a whole session. legini folds
# skspgp(n+1) into fsp, so every conversion legmod makes is filtered; a wrapper
# that omits it is a different operator, and the error GROWS with total
# wavenumber rather than announcing itself at n=1. It is invisible at the
# default nfilter=0, which is why the driver above now runs the beds' nfilter=2.
cp -f "$SRC"/shtnsmod.f90 shtnsmod.f90
before=$(grep -cE 'real\(fsp\(|real\(fgp\(|real\(pfil\(' shtnsmod.f90)
[ "$before" -ge 6 ] || { echo "control patch: expected the filter in every wrapper, found $before"; exit 1; }
sed -i 's/ \* real(fsp(jm),8)//g; s/ \* real(fsp(2),8)//g; s/real(fgp(jm),8)/1.0_8/g; s/real(pfil(jm),8)/1.0_8/g' shtnsmod.f90
after=$(grep -cE 'real\(fsp\(|real\(fgp\(|real\(pfil\(' shtnsmod.f90 || true)   # grep -c exits 1 on no matches, which here is success
[ "$after" = 0 ] || { echo "control patch missed: $after sites still filtered"; exit 1; }
echo "  dropped the filter from $before sites"
build nofilt
if ./equiv_nofilt.x >/dev/null 2>&1; then
    echo "  [ FAIL ] the control passed, so this check cannot see an unfiltered"
    echo "           wrapper -- which is the defect it exists to catch."
    exit 1
else
    echo "  [  ok  ] control rejected, so the filter is actually under test"
fi
