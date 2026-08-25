#!/bin/bash
# Does storing P and Q instead of eight products compute the same model?
#
#   exoplasim/scripts/verify_weight_factorisation_model.sh <bed> <res> <ranks>
#
# Worldbuilding frame: a correctness check on the Vesper climate model's
# Legendre transform. Nothing here is about the simulated planet.
#
# WHAT IS BEING CHECKED. legini used to build eight NCSP x NLPP matrices, each P
# or Q times a per-mode factor and a per-latitude factor. It now builds P and Q
# and applies the factors where they are cheap. `verify_weight_factorisation.py`
# is the other half: it compiles legini out of legmod.f90 and checks the eight
# weights against the same products formed from scipy, so the tables and the six
# per-mode factors are pinned to a reference rather than to each other. This one
# runs the MODEL, which is the part no table comparison reaches: that the
# separated form, applied inside the transforms, still integrates the same
# climate.
#
# WHY THIS CANNOT BE BIT IDENTICAL. a*(b*c) is not (a*b)*c in floating point, and
# the whole change is moving where the multiplication happens. So the bar is a
# regrouped sum -- rounding scale against order unity -- which is what
# compare_restarts.py separates from a different computation. The tolerance is
# declared here, before any arm runs, and is the same 1e-10 the paired
# decomposition and the symmetry work used.
#
# THE CONTROL must fail. It corrupts ONE per-mode factor, fmm, from m*skgpsp to
# n*skgpsp -- a factorisation that is still separable, still the right shape,
# still gives a model that runs, and is wrong. A control that merely zeroed a
# matrix would not stand in for the mistake this change actually invites.
set -euo pipefail

bed="$(cd "${1:?usage: verify_weight_factorisation_model.sh <bed> <res> <ranks>}" && pwd)"
res="${2:?}"
ranks="${3:?}"

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
# shellcheck source=_bed_guard.sh
. "$HERE/_bed_guard.sh"
PKG="$REPO/vendor/exoplasim/exoplasim"
BUILD="$REPO/.venv/bin/python $REPO/exoplasim/scripts/build_model.py"
SRC="$PKG/plasim/src"
low="$(echo "$res" | tr 'A-Z' 'a-z')"
name="most_plasim_${low}_l10_p${ranks}.x"
# exoplasim/bench/_wfcheck, or wherever WFCHECK_WORK points. A git worktree
# reaches exoplasim/bench through a symlink into the MAIN checkout, which
# every worktree shares, so a worktree that runs this without setting the
# variable deletes and rewrites work another one is using -- this script
# opens with `rm -rf "$WORK"`.
WORK="${WFCHECK_WORK:-$REPO/exoplasim/bench/_wfcheck}"
TOL=1e-10

require_settled_bed "$bed"
require_bed_grid "$bed" "$res"

dirty="$(cd "$REPO" && git status --porcelain -- vendor/exoplasim/exoplasim/plasim/src)"
if [ -n "$dirty" ]; then
    echo "refusing: the model source is not what is committed --" >&2
    echo "$dirty" >&2
    echo "  The arms rewrite legmod.f90 in place and restore it afterwards." >&2
    exit 1
fi

rm -rf "$WORK"; mkdir -p "$WORK/ref"
restore() { ( cd "$REPO" && git checkout -- vendor/exoplasim/exoplasim/plasim/src ); }
trap restore EXIT

build_arm() {
    local arm="$1"
    restore
    case "$arm" in
      factored) ;;
      eight)
        # The newest revision that still declares the eight matrices, found
        # rather than counted back to. A fixed HEAD~n is wrong the moment
        # anything else is committed, which it was within the hour.
        local rev
        rev="$(cd "$REPO" && git rev-list HEAD -- vendor/exoplasim/exoplasim/plasim/src/legmod.f90 \
               | while read -r r; do
                   if git show "$r:vendor/exoplasim/exoplasim/plasim/src/legmod.f90" \
                      2>/dev/null | grep -q "qq(NCSP,NLPP)"; then echo "$r"; break; fi
                 done)"
        [ -n "$rev" ] || { echo "no revision of legmod.f90 has the eight-matrix form" >&2; exit 1; }
        echo "  reference is ${rev:0:8}"
        ( cd "$REPO" && git show "$rev:vendor/exoplasim/exoplasim/plasim/src/legmod.f90" \
              > "$SRC/legmod.f90" )
        ;;
      broken)
        sed -i 's/^      fmm(lm) = m \* skgpsp(n+1)$/      fmm(lm) = n * skgpsp(n+1)/' "$SRC/legmod.f90"
        grep -q "fmm(lm) = n \* skgpsp" "$SRC/legmod.f90" || {
            echo "control patch missed" >&2; exit 1; }
        ;;
    esac
    # FRESHNESS: the name lost its parallel mode with world-38b, so it is now
    # the name the MPI builds published. An existence test alone would compare
    # arms against a binary this script did not build.
    local stamp="$WORK/.stamp"; : > "$stamp"
    ( $BUILD --res "$res" --ranks "$ranks" ) \
        >"$WORK/build_$arm.log" 2>&1 || true
    [ -f "$PKG/plasim/run/$name" ] && [ "$PKG/plasim/run/$name" -nt "$stamp" ] || {
        echo "build failed or left an older binary of that name: $arm "\
             "(see $WORK/build_$arm.log)" >&2; exit 1; }
    cp -f "$PKG/plasim/run/$name" "$WORK/ref/$arm.x"
    echo "built $arm  $(sha256sum "$WORK/ref/$arm.x" | cut -c1-16)"
}

run_arm() {
    local arm="$1" steps="$2" tag="$3"
    local d="$WORK/run_$tag"
    rm -rf "$d"; mkdir -p "$d"
    cp -a "$bed"/. "$d"/
    ( cd "$d"
      rm -f ./*.x plasim_status Abort_Message
      sed -i "s/^ *N_RUN_STEPS *=.*/ N_RUN_STEPS = $steps /" plasim_namelist
      cp -f "$WORK/ref/$arm.x" ./probe.x
      export OMP_NUM_THREADS="$ranks" OMP_PROC_BIND=close OMP_PLACES=cores
      export OMP_STACKSIZE=512M
      ulimit -s unlimited
      ./probe.x >run.log 2>&1 ) >/dev/null 2>&1 || true
    [ -f "$d/plasim_status" ]
}

compare() {
    "$REPO"/.venv/bin/python "$REPO"/exoplasim/scripts/compare_restarts.py \
        "$WORK/run_$1/plasim_status" "$WORK/run_$2/plasim_status" \
        --tol "$TOL" --exact dls --exact doro --exact darea --quiet
}

echo "res $res  threads $ranks  tolerance $TOL, declared before the arms ran"
echo
build_arm eight
build_arm factored
build_arm broken

rc=0
for n in 1 20; do
    echo
    echo "==== eight matrices against two, $n step(s): must agree ===="
    if run_arm eight "$n" "e$n" && run_arm factored "$n" "f$n"; then
        compare "e$n" "f$n" || rc=1
    else
        echo "an arm produced no restart"; rc=1
    fi
done

echo
echo "==== the control against the reference, 1 step: must NOT agree ===="
if run_arm broken 1 b1; then
    if compare e1 b1; then
        echo "FAIL: the control agreed, so this check cannot see a wrong factor."
        rc=1
    else
        echo "control rejected, so the comparison has teeth."
    fi
else
    echo "the control did not run; a crash is not a demonstration that the"
    echo "comparison works, so this is not counted as a pass."
    rc=1
fi
echo
echo "work kept at $WORK"
exit $rc
