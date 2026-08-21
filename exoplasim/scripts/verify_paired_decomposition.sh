#!/bin/bash
# Does the paired latitude decomposition compute the same model as the stock one?
#
#   exoplasim/scripts/verify_paired_decomposition.sh <bed> <res> <ranks> [steps]
#
# Worldbuilding frame: a correctness check on the Vesper climate model's MPI
# decomposition. Nothing here is about the simulated planet.
#
# WHAT IS BEING CHECKED. LPAIRLAT permutes the latitude scatter so that a
# latitude and its mirror land on the same process, which is what lets legmod
# fold a mirror pair together before it multiplies by the Legendre weight. The
# permutation is undone on every gather, so the model above mpimod is supposed
# to be unable to tell which layout it is running under.
#
# WHY THIS CANNOT BE A CHECKSUM. Each process sums the transform over the
# latitudes it holds and mpsumsc completes the partial sums, so changing which
# latitudes a process holds regroups a floating point sum. The last bits move
# legitimately. What separates a regrouped sum from a different computation is
# SIZE -- rounding scale against order unity -- and that is what
# compare_restarts.py reports.
#
# Two things in it can still fail hard, and they are the point:
#
#   1. the land-sea mask and the orography are scattered in and gathered back
#      out with no arithmetic between, so they must come back BIT IDENTICAL. A
#      permutation applied on the way in and not on the way out shows up here
#      and nowhere else. darea is checked the same way and is the sharpest of
#      the three: it is the Gaussian weight of a process's OWN latitude, so it
#      only survives the round trip if the sid/gwd scatter is permuted the same
#      way the grid fields are.
#   2. the NEGATIVE CONTROL arm must FAIL the same comparison. Its permutation
#      hands the southern blocks round the processes by one, so it is still a
#      bijection and every process still holds NLPP latitudes -- they are just
#      not each other's mirrors. That is the mistake this change actually
#      invites, and a control that merely corrupts the map would not stand in
#      for it. Without a control the tolerance could be passing everything.
set -euo pipefail

bed="$(cd "${1:?usage: verify_paired_decomposition.sh <bed> <res> <ranks> [steps]}" && pwd)"
res="${2:?}"
ranks="${3:?}"
steps="${4:-20}"

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck source=_bed_guard.sh
. "$(dirname "${BASH_SOURCE[0]}")/_bed_guard.sh"
require_settled_bed "$bed"
PKG="$ROOT/vendor/exoplasim/exoplasim"
SRC="$PKG/plasim/src"
low=$(echo "$res" | tr 'A-Z' 'a-z')
name="most_plasim_${low}_l10_p${ranks}.x"
WORK="$ROOT/exoplasim/bench/_paircheck"

# The tolerance is declared HERE, before any arm runs. A sum of NLAT terms in
# double precision regroups at about NLAT * 2.2e-16, which is 6e-14 at the
# top of the ladder; a few tens of timesteps amplify that by a couple of orders.
# 1e-10 is that bound rounded up. It is not a fitted number and the test does
# not hinge on it: a wrong permutation differs by order unity, thirteen decades
# above the bar, which is why the control arm is here to demonstrate the gap.
TOL=1e-10

dirty="$(git -C "$ROOT" status --porcelain "$SRC" || true)"
if [ -n "$dirty" ]; then
    echo "refusing: the model source is not what is committed --" >&2
    echo "$dirty" >&2
    echo "  The arms patch this tree in place and restore it afterwards, so an" >&2
    echo "  uncommitted edit here is both untested and destroyable." >&2
    exit 1
fi

rm -rf "$WORK"; mkdir -p "$WORK/ref"

# Both arms are patches to plasimmod.f90, so that is the only file to save. It
# is patched IN PLACE, in the tree the model is compiled from, which is why the
# check above refuses to start on a dirty source: an edit landing here while
# the arms are building would be silently compiled into one of them, and the
# restore would then write the wrong thing back.
restore() { cp -f "$WORK/orig_plasimmod.f90" "$SRC/plasimmod.f90"; }
cp -f "$SRC/plasimmod.f90" "$WORK/orig_plasimmod.f90"
trap restore EXIT

build_arm() {
    local arm="$1"
    restore
    case "$arm" in
      paired) ;;
      contig)
        sed -i 's/^      logical,parameter :: LPAIRLAT = .*/      logical,parameter :: LPAIRLAT = .false./' \
            "$SRC/plasimmod.f90"
        grep -q "LPAIRLAT = .false." "$SRC/plasimmod.f90" || { echo "contig patch missed" >&2; exit 1; }
        ;;
      broken)
        sed -i 's/ilatperm = NLAT - ir \* NLHP - NLPP + il/ilatperm = NLAT - mod(ir+1,NPRO) * NLHP - NLPP + il/' \
            "$SRC/plasimmod.f90"
        grep -q "mod(ir+1,NPRO)" "$SRC/plasimmod.f90" || { echo "control patch missed" >&2; exit 1; }
        ;;
    esac
    ( cd "$PKG" && ./compile.sh -n "$ranks" -p 8 -r "$res" -v 10 -O march=znver4 ) >"$WORK/build_$arm.log" 2>&1 || true
    [ -f "$PKG/plasim/run/$name" ] || { echo "build failed: $arm (see $WORK/build_$arm.log)" >&2; exit 1; }
    cp -f "$PKG/plasim/run/$name" "$WORK/ref/$arm.x"
    echo "built $arm  $(sha256sum "$WORK/ref/$arm.x" | cut -c1-16)"
}

run_arm() {
    local arm="$1" nsteps="$2" tag="$3"
    rm -rf "$WORK/run_$tag"; mkdir -p "$WORK/run_$tag"
    cp -a "$bed"/. "$WORK/run_$tag"/
    ( cd "$WORK/run_$tag"
      rm -f MOST_REST.* MOST_DIAG.* plasim_status Abort_Message
      sed -i "s/^ *N_RUN_STEPS *=.*/ N_RUN_STEPS = $nsteps /" plasim_namelist
      cp -f "$WORK/ref/$arm.x" ./probe.x
      mpiexec -np "$ranks" ./probe.x >run.log 2>&1 ) || return 1
    if [ -f "$WORK/run_$tag/Abort_Message" ]; then return 1; fi
    if [ ! -f "$WORK/run_$tag/plasim_status" ]; then return 1; fi
    return 0
}

require_arm() {
    local arm="$1" nsteps="$2"
    if ! run_arm "$arm" "$nsteps" "$arm"; then
        echo "$arm: the MODEL failed, not the comparison. Last lines:" >&2
        tail -8 "$WORK/run_$arm/run.log" >&2
        exit 1
    fi
    echo "ran $arm over $nsteps steps"
}

for arm in paired contig broken; do build_arm "$arm"; done
restore

require_arm paired "$steps"
require_arm contig "$steps"

echo
echo "==== paired against contiguous: must agree ===="
ok=0
"$ROOT/.venv/bin/python" "$ROOT/exoplasim/scripts/compare_restarts.py" \
    "$WORK/run_contig/plasim_status" "$WORK/run_paired/plasim_status" \
    --tol "$TOL" --exact dls --exact doro --exact darea --quiet || ok=$?

# The control, over ONE step rather than the full length. A permutation that
# pairs the wrong latitudes puts each process's Gaussian weights against
# somebody else's fields, and over any distance that goes to a floating point
# exception rather than to a comparable restart. One step keeps it inside the
# window where it still produces one, so the comparison itself is exercised
# rather than merely bypassed by a crash. If it cannot manage even that, the
# rejection stands on the crash and is stated as such.
echo
echo "==== the control against contiguous: must NOT agree ===="
control_rejected=0
control_how=""
if run_arm broken 1 broken1; then
    if run_arm contig 1 contig1; then
        if ! "$ROOT/.venv/bin/python" "$ROOT/exoplasim/scripts/compare_restarts.py" \
             "$WORK/run_contig1/plasim_status" "$WORK/run_broken1/plasim_status" \
             --tol "$TOL" --exact dls --exact doro --exact darea --quiet; then
            control_rejected=1
            control_how="the comparison rejected it after one step"
        fi
    else
        echo "the contiguous arm failed at one step, which is not about the control" >&2
        exit 1
    fi
else
    control_rejected=1
    control_how="the model would not integrate one step with it"
    tail -4 "$WORK/run_broken1/run.log" 2>/dev/null | sed 's/^/    /'
fi

if [ "$control_rejected" -eq 0 ]; then
    echo
    echo "FAIL: the control arm PASSES the same comparison, so the comparison"
    echo "      proves nothing as it stands. Its permutation hands the southern"
    echo "      blocks round the processes by one and cannot be right."
    exit 1
fi
echo "control rejected: $control_how"

echo
if [ "$ok" -eq 0 ]; then
    echo "PASS: paired and contiguous agree at the scale of a regrouped sum,"
    echo "      the round-tripped fields are bit identical, and a wrong"
    echo "      permutation is rejected."
else
    exit 1
fi

echo
echo "NOTE: $PKG/plasim/run/$name is now whichever arm was built last."
echo "      Run exoplasim/scripts/rebuild_binaries.py before anything else uses it."
