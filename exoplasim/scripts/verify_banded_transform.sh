#!/bin/bash
# Build and run the banded spectral analysis against an independent reference.
#
#   exoplasim/scripts/verify_banded_transform.sh [res] [threads] [levels]
#
# Worldbuilding frame: a correctness check on the Vesper climate model's
# spectral transform. Nothing here is about the simulated planet.
#
# WHAT THIS IS. world-38b left this component one build, and
# verify_threaded_numerics.sh was written as that build against an independent
# one. world-d5l is the decision that replaced the missing side: a standalone
# driver that computes the second side rather than a second registered runtime.
# verify_banded_transform.f90 is the driver and its header carries the argument;
# banded_transform_reference.py is the reference and its header states exactly
# what is independent of the model and what is not. This script is the wiring.
#
# WHERE IT WORKS. exoplasim/bench/_bandedtransform, or wherever BANDED_WORK
# points. A git worktree reaches exoplasim/bench through a symlink into the main
# checkout, which is shared, so a worktree sets BANDED_WORK somewhere of its own.
#
# WHAT IT BUILDS. The model's own modules, out of the model tree, at the
# declared flag line with the thread count under test compiled in -- so the
# arithmetic here is the arithmetic the model does. The sources are COPIED into
# the work directory and the controls patch the copies, so nothing here can
# leave the repository's model source modified, however it exits.
#
# THE THREE CASES the reference is built for, and each is a different way for a
# band error to fail to cancel:
#   dense   every mode driven, which is what crosses every band
#   corner  the (NTRU,NTRU) mode alone, the last the recurrence reaches and the
#           one the weights are smallest for; a transform can be right on a
#           dense spectrum and wrong on one coefficient
#   zonal   every m = 0 mode, the column the model's own global sums are taken
#           through
#
# THE TWO CONTROLS, and both must FAIL or the arms are not measuring what they
# say. They are chosen to break DIFFERENT halves:
#
#   noweight  drops the Gaussian weight from fc2sp's accumulation. That is the
#             weight pre-scaling defect's shape -- the quadrature weight in the
#             wrong place -- and it is invisible to a round trip, because
#             synthesis and analysis would drop it together. ARM B must reject
#             it; ARM A must not, since the weights themselves are untouched.
#
#   slideband slides every thread's band down one latitude row, so neighbouring
#             bands share a row. This is the same patch verify_threaded_numerics
#             .sh uses, deliberately: one control for one property, in both
#             gates. ARM B must reject it; ARM A must not, since a thread's
#             weight matrices come from the scatter and not from the band
#             pointer, and that split is the point of having two arms.
#
#   samelats  makes the scatter hand every thread the master's latitudes. ARM A
#             must reject it: this is the case ARM A exists for, a thread whose
#             arithmetic is perfect on the wrong rows of the globe. Without this
#             one, ARM A would be an arm with no demonstration that it can fail.
set -euo pipefail

res="${1:-T21}"
threads="${2:-4}"
lev="${3:-10}"

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
SRC="$REPO/vendor/exoplasim/exoplasim/plasim/src"
WORK="${BANDED_WORK:-$REPO/exoplasim/bench/_bandedtransform}"
PY="$REPO/.venv/bin/python"; [ -x "$PY" ] || PY=python3

CASES="dense corner zonal"

# The ladder is lib/rungs.py and nowhere else.
nlat=$("$PY" -c "import sys; sys.path.insert(0, '$REPO/lib'); import rungs; print(rungs.geometry('$res')[0])") || {
  echo "unknown resolution $res: it is not a rung in lib/rungs.py" >&2; exit 2; }

if [ $((nlat % threads)) -ne 0 ]; then
    echo "refusing: $threads threads does not divide the $nlat latitudes of $res." >&2
    echo "  NLPP is NLAT/NPRO and is a parameter, so the bands would not tile" >&2
    echo "  the globe and part of it would be owned by no thread." >&2
    exit 2
fi

# The model source has to be what is committed, or the arms measure a working
# copy and the result cannot be attributed to anything.
dirty="$(cd "$REPO" && git status --porcelain -- vendor/exoplasim/exoplasim/plasim/src)"
if [ -n "$dirty" ]; then
    echo "refusing: the model source is not what is committed --" >&2
    echo "$dirty" >&2
    exit 1
fi

# The declared flag line, asked of build_model.py rather than restated here, so
# this measures the compiler the model is built with. `checked` is the profile
# for verification arms: config/planet.yaml says so and gives the reason.
FFLAGS="$("$PY" -c "
import sys; sys.path.insert(0, '$REPO/exoplasim/scripts')
import build_model
flags, _prec = build_model.flag_line('checked')
print(' '.join(flags))")"

rm -rf "$WORK"; mkdir -p "$WORK"
cd "$WORK"
cp "$SRC"/plasimmod.f90 "$SRC"/legmod.f90 "$SRC"/gaussmod.f90 \
   "$SRC"/specblock.f90 "$SRC"/mpimod_omp.f90 .
cp "$HERE"/verify_banded_transform.f90 drive.f90

# mpimod_omp's restart and surface transfers call three routines that live in
# the model's file handling. The driver reads no restart and writes no output,
# so they stand in as refusals: reached, they mean the driver did something it
# was not written to do.
cat > iostub.f90 <<'EOF'
      subroutine get_restart_array(yn,pa,k1,k2,klev)
      character (len=*) :: yn
      real :: pa(k1,klev)
      integer :: k1, k2, klev
      write(*,*) 'get_restart_array reached in a driver that reads no restart: ', yn
      stop 3
      end subroutine get_restart_array

      subroutine put_restart_array(yn,pa,k1,k2,klev)
      character (len=*) :: yn
      real :: pa(k1,klev)
      integer :: k1, k2, klev
      write(*,*) 'put_restart_array reached in a driver that writes no restart: ', yn
      stop 3
      end subroutine put_restart_array

      subroutine get_surf_array(yn,pa,k1,klev,kread)
      character (len=*) :: yn
      real :: pa(k1,klev)
      integer :: k1, klev, kread
      write(*,*) 'get_surf_array reached in a driver that reads no surface file: ', yn
      kread = 0
      stop 3
      end subroutine get_surf_array
EOF

cat > resmod.f90 <<EOF
      module resmod ! generated by verify_banded_transform.sh
      parameter(NLAT_ATM = $nlat)
      parameter(NLEV_ATM = $lev)
      parameter(NPRO_ATM = $threads)
      end module resmod
EOF

# -fopenmp and -DOMPSHARED are what build_model.py adds to the declared line for
# every model compile, and they are what makes the bands bands: without them the
# thread pointers are ordinary arrays and there is nothing to check.
F="-c -cpp -DOMPSHARED -fopenmp -ffixed-line-length-132 $FFLAGS"

build () {
    local tag="$1"
    rm -f ./*.o ./*.mod
    for f in resmod plasimmod iostub gaussmod specblock mpimod_omp legmod; do
        # shellcheck disable=SC2086
        gfortran $F -fdefault-real-8 -J . "$f.f90" -o "$f.o" 2>"$f.err" \
            || { echo "compile failed: $f ($tag)"; head -20 "$f.err"; exit 1; }
    done
    # shellcheck disable=SC2086
    gfortran $F -fdefault-real-8 -J . drive.f90 -o drive.o 2>drive.err \
        || { echo "compile failed: drive ($tag)"; head -30 drive.err; exit 1; }
    gfortran -o "run_$tag.x" drive.o legmod.o mpimod_omp.o specblock.o \
        gaussmod.o iostub.o plasimmod.o resmod.o -fopenmp -g 2>link.err \
        || { echo "link failed ($tag)"; head -20 link.err; exit 1; }
}

echo "==== the reference ===="
for c in $CASES; do
    "$PY" "$HERE/banded_transform_reference.py" "$res" --npro "$threads" \
        --case "$c" --out "$WORK/ref_$c.bin"
done

export OMP_NUM_THREADS="$threads" OMP_PROC_BIND=close OMP_PLACES=cores
export OMP_STACKSIZE=512M
ulimit -s unlimited || true

rc=0

echo
echo "==== the model's banded analysis, $res on $threads threads ===="
build good
for c in $CASES; do
    echo
    echo "-- case $c --"
    if ! ./run_good.x "$WORK/ref_$c.bin"; then
        echo "[ FAIL ] case $c did not return the reference answer."
        rc=1
    fi
done

# ---------------------------------------------------------------------------
# The controls. Each patches the COPIED source, rebuilds, and must be rejected.
# A control that is not applied is a control that passes for the wrong reason,
# so each patch is counted rather than assumed.
# ---------------------------------------------------------------------------

control () {
    local tag="$1" want="$2"
    echo
    echo "==== control $tag: must NOT agree ===="
    build "$tag"
    if ./run_"$tag".x "$WORK/ref_dense.bin" >"ctl_$tag.log" 2>&1; then
        echo "[ FAIL ] the control passed, so the arms cannot see $tag."
        sed -n '1,40p' "ctl_$tag.log"
        rc=1
    elif grep -q "\[ FAIL \] $want" "ctl_$tag.log"; then
        echo "[  ok  ] control rejected by $want, so that arm has teeth."
    else
        echo "[ FAIL ] the control was rejected, but not by $want. A control has"
        echo "         to fail the arm it was built to fail, or the rejection"
        echo "         says nothing about that arm."
        sed -n '1,40p' "ctl_$tag.log"
        rc=1
    fi
}

# The Gaussian weight, dropped from the analysis accumulation only. Synthesis is
# untouched, which is the point: a round trip would not notice.
cp plasimmod.f90 plasimmod.orig.f90
cp legmod.f90 legmod.orig.f90
cp mpimod_omp.f90 mpimod_omp.orig.f90
# Six lines carry that product: two in fc2sp, which is the routine under test,
# and two each in qtend and mktend, which this driver never calls. Patching all
# six is the robust spelling and the four spare are inert here.
before=$(grep -c 'pmat(w,l)\*fgp(w)\*gwd(l)' legmod.f90 || true)
[ "$before" -eq 6 ] || { echo "control patch site moved: $before hits, expected 6" >&2; exit 1; }
sed -i 's/pmat(w,l)\*fgp(w)\*gwd(l)/pmat(w,l)*fgp(w)/g' legmod.f90
control noweight "ARM B"
cp legmod.orig.f90 legmod.f90

# The band, slid down one latitude row so neighbouring bands share one. The same
# patch verify_threaded_numerics.sh applies, on the same line, deliberately.
before=$(grep -c '^      lo = mypid \* NHOR + 1$' plasimmod.f90 || true)
[ "$before" -eq 1 ] || { echo "band patch site moved: $before hits, expected 1" >&2; exit 1; }
sed -i 's/^      lo = mypid \* NHOR + 1$/      lo = max(1, mypid * NHOR + 1 - NLON)   ! CONTROL: bands overlap by a row/' \
    plasimmod.f90
control slideband "ARM B"
cp plasimmod.orig.f90 plasimmod.f90

# The scatter, made to hand every thread the master's latitudes. This is ARM A's
# own control: without it ARM A would report a number and never a failure.
before=$(grep -c 'if (mypid /= NROOT) p(1:n) = zbufd(mypid\*n+1:mypid\*n+n)' mpimod_omp.f90 || true)
[ "$before" -eq 1 ] || { echo "scatter patch site moved: $before hits, expected 1" >&2; exit 1; }
sed -i 's|if (mypid /= NROOT) p(1:n) = zbufd(mypid\*n+1:mypid\*n+n)|if (mypid /= NROOT) p(1:n) = zbufd(1:n)   ! CONTROL: every thread takes the same latitudes|' \
    mpimod_omp.f90
control samelats "ARM A"
cp mpimod_omp.orig.f90 mpimod_omp.f90

echo
echo "work kept at $WORK"
exit $rc
