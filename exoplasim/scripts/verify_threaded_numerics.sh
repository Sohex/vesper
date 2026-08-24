#!/bin/bash
# Does the threaded build still compute the reference build's answer?
#
#   exoplasim/scripts/verify_threaded_numerics.sh <bed> <res> <n> [reference]
#
# Worldbuilding frame: a correctness check on the Vesper climate model's
# threaded build. Nothing here is about the simulated planet.
#
# WHY THIS EXISTS AS A REGISTERED CHECK, and why it did not need to before.
# Until the grid fields became bands of a shared globe, the threaded build was
# BIT IDENTICAL to the MPI build at T21 on two, and that was the sharpest tool
# in this component -- it caught the dv2uv planetary vorticity race and the
# weight pre-scaling bug, both of which a tolerance would have let through.
#
# That standard is gone and cannot come back while the layout serves SHTns. A
# thread's band is contiguous within a level and strided between them, so the
# compiler cannot assume contiguity, vectorises differently, and moves last
# bits. What replaces it is this: agreement at the scale of a regrouped sum, AT
# SHORT RANGE, where the seed has not yet amplified. That is a weaker check, so
# it is written down, given a control, and run rather than remembered.
#
# WHY SHORT RANGE IS THE WHOLE POINT. A last-bit difference in this model grows
# by roughly three decades every twenty steps: measured on this change, 3.2e-13
# at one step, 1.4e-13 at twenty, 3.0e-10 at sixty. A comparison at sixty steps
# therefore FAILS a 1e-10 tolerance on a change that is perfectly sound, and
# a comparison at three hundred fails on anything at all. Length is not a
# detail here; it is what separates a defect from Lyapunov growth.
#
# THE CONTROL must fail, and choosing it took a wrong turn worth recording.
# The obvious one -- every thread takes its NEIGHBOUR's band -- PASSES, because
# it is not a mistake. The physics is per-point and does not consult mypid, and
# a band's global position is set by the transfer routines through ilatperm
# rather than by the pointer offset, so if thread 0 works band 1 and thread 1
# works band 0 both bands are still computed correctly. The offset is a
# relabelling, which is exactly why the conversion was inert.
#
# The property that MATTERS is that the bands are DISJOINT. So the control
# breaks disjointness and nothing else: every thread's band slides down by one
# latitude row, so thread 0 and thread 1 share a row, thread 1 and thread 2
# share a row, and so on. An offset that is one row out is the mistake this
# design actually invites.
#
# IT SLIDES BY A ROW RATHER THAN COLLAPSING TO THREAD 0'S BAND, and the
# difference is the whole reliability of the control. Pointing every thread at
# lo = 1 makes four threads apply gp = exp(gp) to the same row, which overflows
# and takes the trap -- SOMETIMES, because it is a race. It ran and disagreed
# when it was written and it crashed the next time it was asked, and a crash is
# not a demonstration that the comparison works. A control whose outcome is up
# to the scheduler is not a control. The one-row slide corrupts deterministically
# and stays in range: the top thread's band ends one row short of the globe, so
# no thread reads past NUGP.
#
# WHERE THESE BOUNDS COME FROM, AND WHICH OF THEM IS DERIVED.
#
# THE REASSOCIATION FLOOR is arithmetic and is derived per rung. Every spectral
# restart record is a `fc2sp` analysis, a sum over the whole grid, so the
# rigorous upper bound on how much a regrouped float64 version of that sum can
# differ, in ANY association order, is NLAT*NLON*eps. It is 4.6e-13 at T21 and
# 2.9e-11 at T170 -- a factor of 64 across the ladder, because the ladder is
# what sets the length of the sum. `floor_at_rung` computes it, and the gate
# REFUSES to run where the declared birth bound is below it: a bound under the
# floor fails every sound change, and 1e-11 is under it from T106 up. That is
# the rung-dependence made loud rather than silent.
#
# THE BIRTH SCALE the floor implies is larger than the floor, because the birth
# norm is taken over every restart record and one of them applies a gain of its
# own. `rainmod.f90:92` sets rcrit = MAX(0.85, MAX(sigma, 1-sigma)) and the
# cloud fraction carries 1/(1-rcrit)^2, which on the linear sigma grid is
# (2*NLEV)^2 = 400 at the top and bottom levels. So a difference born at the
# floor can present as 400 times the floor and still be nothing but
# reassociation. Two measurements exist: 3.2e-13 for the grid-band change here,
# and 3.35e-12 for the SHTns change on the baseline bed, ten times larger for
# exactly this reason (exoplasim/notes/shtns-viability.md).
#
# TOL AND JUMP ARE NOT DERIVED, and the rung they were measured at is not on
# record. Both encode the rate at which the model amplifies a last-bit
# difference, which is a property of the flow and not of the arithmetic, and is
# rung-dependent by construction: a finer truncation resolves faster-growing
# modes and reaches a given norm in fewer steps for entirely sound reasons.
# `world-2ic` is that question, and it names the experiment that settles it:
# run this gate's growth curve unchanged at three rungs on the same terrain and
# fit the per-step factor at each. Until then the two are what they have been,
# and this header is what they rest on.
set -euo pipefail

bed="$(cd "${1:?usage: verify_threaded_numerics.sh <bed> <res> <n> [reference]}" && pwd)"
res="${2:?}"
n="${3:?}"
# THE REFERENCE IS THE SERIAL BUILD. This project builds and runs the threaded
# OpenMP parmode and nothing else; the MPI path is being removed under
# world-38b, so defaulting the reference arm to it would default this gate to a
# configuration that will not exist. `mpi` is still accepted as an explicit
# argument for as long as that build does, because a gate that can compare
# against two independent references is worth more than one that can compare
# against one -- but nothing should reach for it by accident.
reference="${4:-serial}"  # serial or mpi; see the note on the hard fork

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
# shellcheck source=_bed_guard.sh
. "$HERE/_bed_guard.sh"
PKG="$REPO/vendor/exoplasim/exoplasim"
BUILD="$REPO/.venv/bin/python $REPO/exoplasim/scripts/build_model.py"
SRC="$PKG/plasim/src"
low="$(echo "$res" | tr 'A-Z' 'a-z')"
WORK="$REPO/exoplasim/bench/_tnumerics"
TOL=1e-10
STEPS="1 2 5 10 20 40"
BIRTH=1e-11      # the norm at ONE step; above this the change is wrong at birth
JUMP=1e4         # the largest ratio allowed between adjacent samples

require_settled_bed "$bed"
require_bed_grid "$bed" "$res"

# The reassociation floor at THIS rung, and the refusal that goes with it. See
# the header. `floor` is NLAT*NLON*eps, and a birth bound below it cannot be met
# by correct arithmetic.
PYX="$REPO/.venv/bin/python"; [ -x "$PYX" ] || PYX=python3
floor="$("$PYX" -c "import sys; sys.path.insert(0, '$REPO/lib'); import rungs, numpy as np; nlat, nlon, _ = rungs.geometry('$res'); print(repr(nlat * nlon * float(np.finfo(np.float64).eps)))")" || {
    echo "unknown resolution $res: it is not a rung in lib/rungs.py" >&2; exit 2; }
gained="$(awk -v f="$floor" 'BEGIN{printf "%.3g", f*400}')"
if awk -v b="$BIRTH" -v f="$floor" 'BEGIN{exit !(b < f)}'; then
    echo "refusing: the birth bound $BIRTH is below the reassociation floor at" >&2
    echo "  $res, which is $floor = NLAT*NLON * float64 eps. A regrouped sum over" >&2
    echo "  this grid can differ by that much in any association order, so the" >&2
    echo "  bound fails every correct change at this rung. It has to be derived" >&2
    echo "  per rung before this gate runs here -- world-2ic." >&2
    exit 2
fi


dirty="$(cd "$REPO" && git status --porcelain -- vendor/exoplasim/exoplasim/plasim/src)"
if [ -n "$dirty" ]; then
    echo "refusing: the model source is not what is committed --" >&2
    echo "$dirty" >&2
    exit 1
fi

rm -rf "$WORK"; mkdir -p "$WORK/ref"
restore() { ( cd "$REPO" && git checkout -- vendor/exoplasim/exoplasim/plasim/src ); }
trap restore EXIT

build_arm() {
    local arm="$1" stamp="$WORK/.stamp" name flags
    restore
    if [ "$arm" = "wrongband" ]; then
        sed -i 's/^      lo = mypid \* NHOR + 1$/      lo = max(1, mypid * NHOR + 1 - NLON)   ! CONTROL: bands overlap by a row/' \
            "$SRC/plasimmod.f90"
        grep -q "CONTROL: bands overlap by a row" "$SRC/plasimmod.f90" || {
            echo "control patch missed" >&2; exit 1; }
    fi
    case "$arm" in
      reference)
        if [ "$reference" = "serial" ]; then
            flags="--ranks 1 --parmode serial"; name="most_plasim_${low}_l10_p1.x"
        else
            flags="--ranks $n --parmode mpi";  name="most_plasim_${low}_l10_p${n}.x"
        fi ;;
      *) flags="--ranks $n --parmode omp"; name="most_plasim_${low}_l10_p${n}_omp.x" ;;
    esac
    : > "$stamp"
    # shellcheck disable=SC2086
    ( $BUILD --res "$res" $flags ) >"$WORK/build_$arm.log" 2>&1 || true
    [ -f "$PKG/plasim/run/$name" ] && [ "$PKG/plasim/run/$name" -nt "$stamp" ] || {
        echo "build failed or stale: $arm (see $WORK/build_$arm.log)" >&2; exit 1; }
    cp -f "$PKG/plasim/run/$name" "$WORK/ref/$arm.x"
    echo "  built $arm  $(sha256sum "$WORK/ref/$arm.x" | cut -c1-16)"
}

run_arm() {
    local arm="$1" steps="$2" tag="$3"
    local d="$WORK/run_$tag"
    rm -rf "$d"; mkdir -p "$d"
    cp -a "$bed"/. "$d"/
    ( cd "$d"
      rm -f ./*.x plasim_status Abort_Message
      sed -i "s/^ *N_RUN_STEPS *=.*/ N_RUN_STEPS = $steps /" plasim_namelist
      # NSHTNS=0 ON BOTH ARMS, deliberately. This check is about the shared grid
      # bands, not the transform: the MPI build cannot run SHTns at all, so
      # leaving the threaded arm on its new default would compare two different
      # transforms and read the difference as a band problem. It would also
      # blunt the control, because SHTns rewrites the whole globe every step and
      # masks an overlapping band -- which is exactly how a committed control
      # patch once hid in plain sight.
      sed -i "/^ *NSHTNS *=/d" plasim_namelist
      sed -i "2i\\ NSHTNS      =     0" plasim_namelist
      cp -f "$WORK/ref/$arm.x" ./probe.x
      if [ "$arm" = "reference" ] && [ "$reference" != "serial" ]; then
          mpiexec -np "$n" ./probe.x >run.log 2>&1
      elif [ "$arm" = "reference" ]; then
          ./probe.x >run.log 2>&1
      else
          export OMP_NUM_THREADS="$n" OMP_PROC_BIND=close OMP_PLACES=cores
          export OMP_STACKSIZE=512M
          ulimit -s unlimited
          ./probe.x >run.log 2>&1
      fi ) >/dev/null 2>&1 || true
    [ -f "$d/plasim_status" ]
}

compare() {
    "$REPO"/.venv/bin/python "$REPO"/exoplasim/scripts/compare_restarts.py \
        "$WORK/run_$1/plasim_status" "$WORK/run_$2/plasim_status" \
        --tol "$TOL" --exact dls --exact doro --exact darea --quiet
}

norm() {
    "$REPO"/.venv/bin/python "$REPO"/exoplasim/scripts/compare_restarts.py \
        "$WORK/run_$1/plasim_status" "$WORK/run_$2/plasim_status" \
        --tol "$TOL" --norm 2>/dev/null
}

echo "$res, $n threads against the $reference build"
echo "declared before the arms ran: tolerance $TOL, lengths [$STEPS],"
echo "derived at $res: reassociation floor $floor = NLAT*NLON * float64 eps,"
echo "  and $gained once rainmod.f90:92's (2*NLEV)^2 record gain is allowed for"
echo "  birth bound $BIRTH at one step, jump bound ${JUMP}x between samples"
echo
build_arm reference
build_arm threaded
build_arm wrongband

rc=0

# HOW THE DIFFERENCE GROWS, which is the part two isolated tolerance checks
# cannot tell you. A sound change starts at rounding scale and grows smoothly as
# the model's own sensitivity amplifies it; a defect that only fires under some
# condition -- a guard that opens on a particular step, a branch reached once a
# field crosses a threshold -- puts a STEP in the curve instead.
#
# NOT gated on monotonicity, and that is measured rather than assumed. The
# statistic is the worst record's relative difference, and which record is worst
# changes with length, so the curve legitimately goes DOWN: on the grid-band
# change it ran 3.2e-13 at one step, 1.4e-13 at twenty, 3.0e-10 at sixty. A
# monotonicity gate would have failed a sound change.
#
# So the curve is the diagnostic and the gates are deliberately coarse: wrong at
# birth, or a jump too large to be amplification.
echo
echo "==== how the difference grows ===="
printf '  %8s  %14s  %s\n' steps norm "worst record"
prev=""
for s in $STEPS; do
    if run_arm reference "$s" "r$s" && run_arm threaded "$s" "t$s"; then
        read -r v rec <<<"$(norm "r$s" "t$s")"
        printf '  %8s  %14s  %s\n' "$s" "$v" "$rec"
        if [ "$s" = 1 ]; then
            over=$(awk -v a="$v" -v b="$BIRTH" 'BEGIN{print (a>b)?1:0}')
            if [ "$over" = 1 ]; then
                echo "  FAIL: $v at one step is above the birth bound $BIRTH."
                echo "        A change that is wrong at step one is not amplification."
                rc=1
            fi
        fi
        if [ -n "$prev" ]; then
            big=$(awk -v a="$v" -v b="$prev" -v j="$JUMP" \
                  'BEGIN{print (b>0 && a/b>j)?1:0}')
            if [ "$big" = 1 ]; then
                echo "  FAIL: the norm jumped by more than ${JUMP}x into $s steps."
                echo "        Amplification is smooth; a step in the curve is a"
                echo "        condition being met, not a seed growing."
                rc=1
            fi
        fi
        prev="$v"
    else
        echo "  $s steps: an arm produced no restart"; rc=1
    fi
done

for s in 1 20; do
    echo
    echo "==== threaded against $reference, $s step(s): must agree ===="
    compare "r$s" "t$s" || rc=1
done

echo
# WHAT COUNTS AS THE CONTROL BEING REJECTED, and the distinction is not the one
# this block used to draw. It refused to count a crash, on the grounds that a
# crash could come from anything -- which is right about the SUBJECT of a check
# and wrong about its CONTROL. A subject arm that does not run leaves the
# question open. A control arm that does not run has been DETECTED, loudly, and
# a detection is what a control is for. The failure it guards against is the
# control passing quietly, not the control failing in a way that is hard to
# read.
#
# The reason a crash is the expected outcome here rather than a surprise: the
# band offset slides the whole band, so every value in it is misaligned against
# the latitude arrays, not just the overlapping row. That state does not survive
# the shortwave radiation with -ffpe-trap on, and it should not.
#
# So the crash counts, but only once the alternatives are excluded: the patch
# must have applied and the build must have succeeded (build_arm exits
# otherwise), and the SAME binary configuration must have run to completion
# unpatched, which the threaded arm above did. Without those, a crash means the
# harness broke rather than the control being caught.
echo "==== the shared-band control, 1 step: must NOT agree ===="
if run_arm wrongband 1 w1; then
    if compare r1 w1; then
        echo "FAIL: the bands overlapped by a row and the comparison did not"
        echo "      notice, so it cannot see bands that are not disjoint."
        rc=1
    else
        echo "control rejected on the numbers, so the comparison has teeth."
    fi
else
    echo "control rejected by refusing to run: the misaligned band did not"
    echo "survive the physics. The build succeeded and the patch applied, and"
    echo "the same binary ran unpatched above, so this is the control being"
    echo "caught rather than the harness breaking."
fi
echo
echo "work kept at $WORK"
exit $rc
