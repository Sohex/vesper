#!/bin/bash
# Build and run the transform round-trip identity check.
#
#   exoplasim/scripts/verify_transform_roundtrip.sh [res] [levels]
#
# The identity and why it outlives every other check on this transform are in
# the header of verify_transform_roundtrip.f90. This script compiles it against
# the model's OWN modules -- not a lifted copy -- at one process, so the whole
# globe is local and analysis is a complete quadrature rather than a partial
# sum needing a reduction.
#
# THE NEGATIVE CONTROL drops the Gaussian weight from the analysis direction.
# That is the mistake the pair invites: synthesis and analysis are inverse only
# BECAUSE the quadrature weight is there, it is easy to fold into the wrong
# matrix, and a model with it missing still runs and still looks like weather.
# The control must fail, or the check is not measuring the inverse property.
set -euo pipefail

res="${1:-T42}"
lev="${2:-10}"

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
SRC="$REPO/vendor/exoplasim/exoplasim/plasim/src"
W="${TMPDIR:-/tmp}/verify_roundtrip.$$"
trap 'rm -rf "$W"' EXIT
mkdir -p "$W"; cd "$W"

case "$res" in
  T21)  nlat=32  ;;  T31)  nlat=48  ;;  T42) nlat=64 ;;
  T85)  nlat=128 ;;  T127) nlat=192 ;;  T170) nlat=256 ;;
  *) echo "unknown resolution $res" >&2; exit 2 ;;
esac

cp "$SRC"/plasimmod.f90 "$SRC"/legmod.f90 "$SRC"/fftmod.f90 \
   "$SRC"/gaussmod.f90 "$SRC"/specblock.f90 . 2>/dev/null
cp "$HERE"/verify_transform_roundtrip.f90 drive.f90

# One process: the whole globe is local, so analysis is the complete quadrature
# and the identity holds without a reduction across threads.
cat > resmod.f90 <<EOF
      module resmod
      parameter(NLAT_ATM = $nlat)
      parameter(NLEV_ATM = $lev)
      parameter(NPRO_ATM = 1)
      end module resmod
EOF

F="-c -O2 -cpp -ffixed-line-length-132 -ffpe-summary=none -finit-real=zero"
F="$F -fdefault-real-8"

build () {
    local tag="$1"
    rm -f ./*.o ./*.mod
    for f in resmod plasimmod gaussmod specblock fftmod legmod; do
        gfortran $F -J . "$f.f90" -o "$f.o" 2>"$f.err" \
            || { echo "compile failed: $f ($tag)"; head -12 "$f.err"; exit 1; }
    done
    gfortran -O2 -fdefault-real-8 -ffixed-line-length-132 -J . -c drive.f90 -o drive.o \
        2>drive.err || { echo "compile failed: drive ($tag)"; head -20 drive.err; exit 1; }
    gfortran -o "run_$tag.x" drive.o legmod.o fftmod.o gaussmod.o specblock.o \
        plasimmod.o resmod.o 2>link.err \
        || { echo "link failed ($tag)"; head -12 link.err; exit 1; }
}

echo "==== the transform pair, $res ===="
build good
./run_good.x | grep -vE "^ \*|^$" || { echo "the identity does not hold"; exit 1; }

echo
echo "==== the control, analysis without its quadrature weight: must FAIL ===="
# The weight is applied inline in fc2sp since the factorisation removed the
# stored matrices, so it is dropped there. Synthesis is untouched, which is the
# point: only the inverse property breaks.
before=$(grep -c 'gwd(l)' legmod.f90)
sed -i 's/pmat(w,l)\*fgp(w)\*gwd(l)/pmat(w,l)*fgp(w)/g' legmod.f90
after=$(grep -c 'gwd(l)' legmod.f90 || true)   # grep -c exits 1 on zero matches, which here means success
[ "$after" -lt "$before" ] || { echo "control patch missed: gwd(l) count unchanged at $before"; exit 1; }
echo "  dropped the quadrature weight from $((before-after)) accumulations"

build bad
if ./run_bad.x >/dev/null 2>&1; then
    echo "  [ FAIL ] the round trip still held without the quadrature weight,"
    echo "           so it is not testing the inverse property."
    exit 1
else
    echo "  [  ok  ] control rejected, so the identity has teeth"
fi
