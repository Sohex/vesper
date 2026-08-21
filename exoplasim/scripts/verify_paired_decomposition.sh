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
    local arm="$1"
    rm -rf "$WORK/run_$arm"; mkdir -p "$WORK/run_$arm"
    cp -a "$bed"/. "$WORK/run_$arm"/
    ( cd "$WORK/run_$arm"
      rm -f MOST_REST.* MOST_DIAG.* plasim_status Abort_Message
      sed -i "s/^ *N_RUN_STEPS *=.*/ N_RUN_STEPS = $steps /" plasim_namelist
      cp -f "$WORK/ref/$arm.x" ./probe.x
      if ! mpiexec -np "$ranks" ./probe.x >run.log 2>&1; then
          echo "$arm: the MODEL failed, not the comparison. Last lines:" >&2
          tail -8 run.log >&2
          exit 1
      fi
      if [ -f Abort_Message ]; then echo "$arm ABORTED" >&2; exit 1; fi
      [ -f plasim_status ] || { echo "$arm produced no plasim_status" >&2; exit 1; } )
    echo "ran $arm"
}

for arm in paired contig broken; do build_arm "$arm"; done
restore
for arm in paired contig broken; do run_arm "$arm"; done

echo
echo "==== paired against contiguous: must agree ===="
ok=0
"$ROOT/.venv/bin/python" "$ROOT/exoplasim/scripts/compare_restarts.py" \
    "$WORK/run_contig/plasim_status" "$WORK/run_paired/plasim_status" \
    --tol "$TOL" --exact dls --exact doro --exact darea --quiet || ok=$?

echo
echo "==== the control against contiguous: must NOT agree ===="
if "$ROOT/.venv/bin/python" "$ROOT/exoplasim/scripts/compare_restarts.py" \
    "$WORK/run_contig/plasim_status" "$WORK/run_broken/plasim_status" \
    --tol "$TOL" --exact dls --exact doro --exact darea --quiet; then
    echo
    echo "FAIL: the control arm, whose permutation strides wrong, PASSES the"
    echo "      same comparison. The comparison proves nothing as it stands."
    exit 1
fi

echo
if [ "$ok" -eq 0 ]; then
    echo "PASS: paired and contiguous agree at the scale of a regrouped sum,"
    echo "      the round-tripped fields are bit identical, and a wrong"
    echo "      permutation is rejected by the same test."
else
    exit 1
fi

echo
echo "NOTE: $PKG/plasim/run/$name is now whichever arm was built last."
echo "      Run exoplasim/scripts/rebuild_binaries.py before anything else uses it."
