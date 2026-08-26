#!/bin/bash
# Does the threaded build compute the right answer over its grid bands?
#
#   exoplasim/scripts/verify_threaded_numerics.sh <bed> <res> <n>
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
# WHAT REPLACED THAT COMPARISON, TWICE. The layout change took bit identity
# away: a thread's band is contiguous within a level and strided between them,
# so the compiler cannot assume contiguity, vectorises differently and moves
# last bits. What stood in its place was agreement between the threaded and the
# MPI build at short range, and world-38b then removed the MPI build, leaving
# that arm with nothing on the other side. world-d5l took the decision this file
# now implements: build a standalone INDEPENDENT driver rather than keep only
# the control arm or retire the gate.
#
# SO THE SECOND SIDE IS A DRIVER AND NOT A BUILD, and that is stronger than what
# it replaces rather than weaker. `verify_banded_transform.sh` computes the band
# decomposition's answer from a table scipy produces and a quadrature pinned by
# two identities, and requires the model's own scatter, weight matrices, banded
# partial sums and cross-thread reduction to return it, at rounding scale, with
# three controls that must fail. Two builds of one source can agree and both be
# wrong; a right answer cannot be agreed with wrongly. Read that script's header
# and the driver's for what is independent, what is not, and what the choice
# makes invisible.
#
# WHAT THE DRIVER CANNOT DO is integrate the model. It says nothing about a
# defect that needs a timestep to appear, which is why the arm below still runs
# the model itself.
#
# THE GROWTH CURVE THAT USED TO BE HERE NEEDS TWO MODEL ARMS THAT DIFFER AT
# LAST-BIT SCALE, and world-38b left one build, so it is not run here. Nothing
# is lost: `verify_shtns_model.sh` carries the same curve, the same
# reassociation floor, the same birth and jump bounds and the same reference to
# `world-2ic` on a pair that is still live -- NSHTNS=0 against NSHTNS=1 on one
# binary. That is where the rung-dependence experiment world-2ic names is run.
#
# A SECOND MODEL ARM IS AVAILABLE AND IS DELIBERATELY NOT TAKEN. WORLD-1YQR
# found that `build_model.py` still accepts `--ranks 1`: it refuses only ranks
# below one and counts that do not divide NLAT, so `most_plasim_<res>_l10_p1.x`
# is a buildable, differently named binary with NPRO=1, one band covering the
# globe and mpimod_omp's reductions collapsed to a single slot. That is a
# genuinely different execution path, and it is NOT the trap world-d5l's refusal
# was written against -- that trap was the mpi arm building the threaded binary
# under the threaded binary's own name, and p1 has a name of its own.
#
# The decision, taken 2026-08-26, is to keep ONE model arm here. Two builds of
# one source can agree and both be wrong: a defect in code both configurations
# run is invisible to their comparison, which is the whole reason the
# independent driver was built instead of a second build. What a p1 arm would
# restore is the growth curve at last-bit scale on the band question, and
# `verify_shtns_model.sh` already carries that curve with its birth bound, its
# jump bound and its reassociation floor. Buying a second arm that shares the
# source to re-run a curve that is already run elsewhere adds a build
# configuration to maintain and no coverage. Reopen this only if the shtns pair
# stops being live.
#
# THE CONTROL must fail, and choosing it took a wrong turn worth recording.
# The obvious one -- every thread takes its NEIGHBOUR's band -- PASSES, because
# it is not a mistake. The physics is per-point and does not consult mypid, and
# a band's global position is set by the transfer routines' own
# `jg = mypid*NLPP + jlat` arithmetic -- mpimod_omp's mpscgp and mpgagp --
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
# IT IS THE SAME PATCH `verify_banded_transform.sh` APPLIES, on the same line,
# deliberately. One property, one control, in both gates: there the driver
# catches it in the transform alone, here the model catches it after a timestep
# of physics, and a control that differed between them would leave open which
# of the two had changed.
set -euo pipefail

bed="$(cd "${1:?usage: verify_threaded_numerics.sh <bed> <res> <n>}" && pwd)"
res="${2:?}"
n="${3:?}"
# THE FOURTH ARGUMENT IS GONE. It named a reference BUILD, and `serial` and
# `mpi` are not build modes any more. Dropping it silently would have been the
# trap: the name the mpi arm looked for is now the threaded binary's own name,
# so an unflagged reference arm would have built the thing under test, compared
# it against itself and passed. A caller who still names one is asking for a
# comparison this gate does not perform, so it is refused rather than ignored.
if [ $# -gt 3 ]; then
    echo "refusing: there is no reference BUILD to name any more." >&2
    echo "  world-38b left one build. The independent side is now" >&2
    echo "  verify_banded_transform.sh, which this gate runs itself." >&2
    exit 2
fi

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
# shellcheck source=_bed_guard.sh
. "$HERE/_bed_guard.sh"
PKG="$REPO/vendor/exoplasim/exoplasim"
BUILD="$REPO/.venv/bin/python $REPO/exoplasim/scripts/build_model.py"
SRC="$PKG/plasim/src"
# exoplasim/bench/_tnumerics, or wherever TNUMERICS_WORK points. A git worktree
# reaches exoplasim/bench through a symlink into the MAIN checkout, which
# every worktree shares, so a worktree that runs this without setting the
# variable deletes and rewrites work another one is using -- this script
# opens with `rm -rf "$WORK"`.
WORK="${TNUMERICS_WORK:-$REPO/exoplasim/bench/_tnumerics}"

require_settled_bed "$bed"
require_bed_grid "$bed" "$res"

dirty="$(cd "$REPO" && git status --porcelain -- vendor/exoplasim/exoplasim/plasim/src)"
if [ -n "$dirty" ]; then
    echo "refusing: the model source is not what is committed --" >&2
    echo "$dirty" >&2
    exit 1
fi

rm -rf "$WORK"; mkdir -p "$WORK/ref"
restore() { ( cd "$REPO" && git checkout -- vendor/exoplasim/exoplasim/plasim/src ); }
trap restore EXIT

rc=0

# ---------------------------------------------------------------------------
# The independent side. It builds and runs on its own and needs no bed: its
# subject is the transform chain rather than an integration, which is exactly
# why it can have a right answer at all.
# ---------------------------------------------------------------------------
echo "==== the independent side: the banded analysis against its reference ===="
if "$HERE/verify_banded_transform.sh" "$res" "$n"; then
    echo "the model's banded analysis returns the independently computed answer."
else
    echo "FAIL: the banded analysis does not return the reference answer."
    rc=1
fi

build_arm() {
    local arm="$1" stamp="$WORK/.stamp" built
    restore
    if [ "$arm" = "wrongband" ]; then
        sed -i 's/^      lo = mypid \* NHOR + 1$/      lo = max(1, mypid * NHOR + 1 - NLON)   ! CONTROL: bands overlap by a row/' \
            "$SRC/plasimmod.f90"
        grep -q "CONTROL: bands overlap by a row" "$SRC/plasimmod.f90" || {
            echo "control patch missed" >&2; exit 1; }
    fi
    : > "$stamp"
    # --no-publish AND --print-path, for two reasons that both bite here.
    # An ARM must not land in the model run directory under the registry's
    # naming: publishing overwrites the shipped executable and leaves
    # check_consistency reporting the binary a run would pick up as having
    # unknown provenance, which is rule 4 reached from inside a check
    # (world-v3d). And the control arm CANNOT publish at all -- it writes a
    # `! CONTROL:` marker into plasimmod.f90 and build_model.py refuses to put
    # a patched source under the registry's tag (world-70k), so this gate's
    # control build failed outright until it stopped asking to. The patched
    # build lands under build/patched/ with a hash of what it patched, so the
    # two arms cannot share objects either.
    built=$( $BUILD --res "$res" --ranks "$n" --no-publish --print-path \
        2>"$WORK/build_$arm.log" ) || {
        echo "build failed: $arm (see $WORK/build_$arm.log)" >&2; exit 1; }
    [ -f "$built" ] && [ "$built" -nt "$stamp" ] || {
        echo "build produced nothing newer than the stamp: $arm "\
             "(see $WORK/build_$arm.log)" >&2; exit 1; }
    cp -f "$built" "$WORK/ref/$arm.x"
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
      # NSHTNS=0 ON EVERY ARM, deliberately. This check is about the shared
      # grid bands, not the transform, so the arms have to differ in the bands
      # and in nothing else. It also keeps the control sharp: SHTns rewrites
      # the whole globe every step and masks an overlapping band, which is
      # exactly how a committed control patch once hid in plain sight.
      sed -i "/^ *NSHTNS *=/d" plasim_namelist
      sed -i "2i\\ NSHTNS      =     0" plasim_namelist
      cp -f "$WORK/ref/$arm.x" ./probe.x
      export OMP_NUM_THREADS="$n" OMP_PROC_BIND=close OMP_PLACES=cores
      export OMP_STACKSIZE=512M
      ulimit -s unlimited
      ./probe.x >run.log 2>&1 ) >/dev/null 2>&1 || true
    [ -f "$d/plasim_status" ]
}

compare() {
    "$REPO"/.venv/bin/python "$REPO"/exoplasim/scripts/compare_restarts.py \
        "$WORK/run_$1/plasim_status" "$WORK/run_$2/plasim_status" \
        --tol 1e-10 --exact dls --exact doro --exact darea --quiet
}

echo
echo "==== the model itself, $res on $n threads ===="
build_arm threaded
build_arm wrongband

if ! run_arm threaded 1 t1; then
    echo "FAIL: the threaded build produced no restart at one step, so the"
    echo "      control below has nothing to be compared against."
    rc=1
fi

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
    if compare t1 w1; then
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
