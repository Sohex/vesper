#!/bin/bash
# Does the SHTns path compute the same model as legmod?
#
#   exoplasim/scripts/verify_shtns_model.sh <bed> <res> <n>
#
# Worldbuilding frame: a correctness check on the Vesper climate model's
# spectral transform. Nothing here is about the simulated planet.
#
# WHY THIS EXISTS ALONGSIDE verify_shtns_equivalence.sh. That one drives fields
# through both transforms and compares the arrays, which is the sharper tool and
# the right place for the convention controls -- normalisation, phase, the sign
# of the toroidal potential. It is also the one that passed at 5e-14 while the
# model disagreed by 100% at the first step, because it set plavor to zero and
# left the spectral filter at its default of off. Both are now switched on
# there. But the general lesson does not go away by fixing two instances of it:
# an array comparison certifies the transform, not the MODEL, and the model
# reaches the transform through a namelist, a build configuration and a call
# site that the array comparison never touches.
#
# So this check runs the model. One binary, two settings of NSHTNS, same bed,
# and the restarts have to agree.
#
# ONE BINARY, NOT TWO. NSHTNS is a namelist switch precisely so that the two
# transforms can be compared without a second build in the comparison, which
# would put compiler differences in the same column as transform differences.
#
# WHY IT NEEDS THE UNPAIRED BUILD. SHTns wants the grid in latitude order and
# LPAIRLAT permutes it, so the binary is built with -u. That is not a
# concession: LPAIRLAT exists to let legmod fold a mirror pair together, and
# SHTns replaces legmod. The nshtns=0 arm therefore runs legmod's contiguous
# branch, which is the honest reference for a build that has no paired layout.
#
# WHY THE LENGTHS ARE SHORT, and the bounds: see the long argument in
# verify_threaded_numerics.sh. A last-bit difference in this model grows by
# roughly three decades every twenty steps, so a comparison at three hundred
# steps fails on anything at all. The tolerance and both bounds are taken
# UNCHANGED from that check, which fixed them before this one existed.
#
# THE CONTROL is the defect that got through. legini folds skspgp(n+1) into
# fsp, so every conversion legmod performs is filtered; a wrapper that omits it
# is a different operator whose error grows with total wavenumber instead of
# announcing itself at n=1. Dropping fsp from the wrappers must break this.
set -euo pipefail

bed="$(cd "${1:?usage: verify_shtns_model.sh <bed> <res> <n>}" && pwd)"
res="${2:?}"
n="${3:?}"

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
# shellcheck source=_bed_guard.sh
. "$HERE/_bed_guard.sh"
PKG="$REPO/vendor/exoplasim/exoplasim"
SRC="$PKG/plasim/src"
low="$(echo "$res" | tr 'A-Z' 'a-z')"
WORK="$REPO/exoplasim/bench/_shtnsmodel"
TOL=1e-10
STEPS="1 2 5 10 20 40"
BIRTH=1e-11      # the norm at ONE step; above this the change is wrong at birth
JUMP=1e4         # the largest ratio allowed between adjacent samples

require_settled_bed "$bed"

if [ ! -f "$REPO/vendor/shtns-install/include/shtns.f03" ]; then
    echo "refusing: no SHTns at vendor/shtns-install." >&2
    echo "  exoplasim/scripts/build_shtns.sh installs it; nshtns=1 needs it." >&2
    exit 1
fi

dirty="$(cd "$REPO" && git status --porcelain -- vendor/exoplasim/exoplasim/plasim/src)"
if [ -n "$dirty" ]; then
    echo "refusing: the model source is not what is committed --" >&2
    echo "$dirty" >&2
    exit 1
fi

rm -rf "$WORK"; mkdir -p "$WORK/bin"
restore() { ( cd "$REPO" && git checkout -- vendor/exoplasim/exoplasim/plasim/src ); }
trap restore EXIT

build_arm() {
    local arm="$1"
    local stamp="$WORK/.stamp"
    local name="most_plasim_${low}_l10_p${n}_omp_np.x"
    restore
    if [ "$arm" = "nofilter" ]; then
        local before after
        before=$(grep -cE 'real\(fsp\(|real\(fgp\(|real\(pfil\(' "$SRC/shtnsmod.f90")
        [ "$before" -ge 6 ] || {
            echo "control patch: expected the filter in every wrapper, found $before" >&2
            exit 1; }
        sed -i 's/ \* real(fsp(jm),8)//g; s/ \* real(fsp(2),8)//g; s/real(fgp(jm),8)/1.0_8/g; s/real(pfil(jm),8)/1.0_8/g' "$SRC/shtnsmod.f90"
        # grep -c exits 1 on no matches, which here is exactly success
        after=$(grep -cE 'real\(fsp\(|real\(fgp\(|real\(pfil\(' "$SRC/shtnsmod.f90" || true)
        [ "$after" = 0 ] || {
            echo "control patch missed: $after sites still filtered" >&2; exit 1; }
    fi
    : > "$stamp"
    ( cd "$PKG" && ./compile.sh -j -u -p 8 -n "$n" -r "$res" -v 10 ) \
        >"$WORK/build_$arm.log" 2>&1 || true
    [ -f "$PKG/plasim/run/$name" ] && [ "$PKG/plasim/run/$name" -nt "$stamp" ] || {
        echo "build failed or stale: $arm (see $WORK/build_$arm.log)" >&2; exit 1; }
    cp -f "$PKG/plasim/run/$name" "$WORK/bin/$arm.x"
    echo "  built $arm  $(sha256sum "$WORK/bin/$arm.x" | cut -c1-16)"
}

run_arm() {
    local arm="$1" sw="$2" steps="$3" tag="$4"
    local d="$WORK/run_$tag"
    rm -rf "$d"; mkdir -p "$d"
    cp -a "$bed"/. "$d"/
    ( cd "$d"
      rm -f ./*.x plasim_status Abort_Message
      sed -i "s/^ *N_RUN_STEPS *=.*/ N_RUN_STEPS = $steps /" plasim_namelist
      sed -i "/^ *NSHTNS *=/d" plasim_namelist
      sed -i "2i\\ NSHTNS      =     $sw" plasim_namelist
      cp -f "$WORK/bin/$arm.x" ./probe.x
      export OMP_NUM_THREADS="$n" OMP_PROC_BIND=close OMP_PLACES=cores
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

norm() {
    "$REPO"/.venv/bin/python "$REPO"/exoplasim/scripts/compare_restarts.py \
        "$WORK/run_$1/plasim_status" "$WORK/run_$2/plasim_status" \
        --tol "$TOL" --norm 2>/dev/null
}

echo "$res, $n thread(s), unpaired build: NSHTNS=1 against NSHTNS=0"
echo "declared before the arms ran: tolerance $TOL, lengths [$STEPS],"
echo "  birth bound $BIRTH at one step, jump bound ${JUMP}x between samples"
echo "  (all three taken unchanged from verify_threaded_numerics.sh)"
echo
build_arm shipped

rc=0

echo
echo "==== how the difference grows ===="
printf '  %8s  %14s  %s\n' steps norm "worst record"
prev=""
for s in $STEPS; do
    if run_arm shipped 0 "$s" "l$s" && run_arm shipped 1 "$s" "s$s"; then
        read -r v rec <<<"$(norm "l$s" "s$s")"
        printf '  %8s  %14s  %s\n' "$s" "$v" "$rec"
        if [ "$s" = 1 ]; then
            over=$(awk -v a="$v" -v b="$BIRTH" 'BEGIN{print (a>b)?1:0}')
            if [ "$over" = 1 ]; then
                echo "  FAIL: $v at one step is above the birth bound $BIRTH."
                echo "        A transform that is wrong at step one is not amplification."
                rc=1
            fi
        fi
        if [ -n "$prev" ]; then
            big=$(awk -v a="$v" -v b="$prev" -v j="$JUMP" \
                  'BEGIN{print (b>0 && a/b>j)?1:0}')
            if [ "$big" = 1 ]; then
                echo "  FAIL: the norm jumped by more than ${JUMP}x into $s steps."
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
    echo "==== NSHTNS=1 against NSHTNS=0, $s step(s): must agree ===="
    compare "l$s" "s$s" || rc=1
done

# RUN TO RUN, which agreement with legmod does not imply. Every wrapper is
# called by the whole team at once -- the worksharing is over levels -- and
# concurrent calls on one configuration were measured safe in isolation
# (probe_shtns_concurrency.f90). This is that measurement made again with the
# model around it, where the reduction order and the barriers are also in play.
# Archive CLIM-44 is this model failing to be reproducible and it has been paid
# for once.
# FOUR RUNS, NOT TWO, and the number was bought. SHTns's default grid setup
# BENCHMARKS its algorithm variants and keeps the winner, so which one it uses
# depends on machine timing at startup and different variants round differently.
# The runs then cluster: three identical, then a different one, then three of
# those. Two runs agreeing said nothing, and it said nothing for long enough to
# send a whole afternoon after a data race that did not exist -- the model was
# nondeterministic on ONE thread, which no race explains. shtns_setup asks for
# SHT_QUICK_INIT now, which skips the benchmark.
REPEATS=4
echo
echo "==== NSHTNS=1 $REPEATS times, 20 steps: must be BIT identical ===="
first=""
for r in $(seq 2 "$REPEATS"); do
    if run_arm shipped 1 20 "s20r$r"; then
        h="$(sha256sum "$WORK/run_s20r$r/plasim_status" | cut -c1-16)"
    else
        echo "  [ FAIL ] repeat $r produced no restart"; rc=1; continue
    fi
    [ -z "$first" ] && first="$(sha256sum "$WORK/run_s20/plasim_status" | cut -c1-16)"
    if [ "$h" != "$first" ]; then
        echo "  [ FAIL ] $first then $h at repeat $r: not reproducible"
        rc=1
    fi
done
[ -n "$first" ] && echo "  [  ok  ] $first, $REPEATS runs" 

echo
echo "==== the control, spectral filter dropped: must NOT agree ===="
build_arm nofilter
if run_arm nofilter 1 1 c1; then
    if compare l1 c1; then
        echo "  [ FAIL ] the control agreed, so this check cannot see an"
        echo "           unfiltered wrapper -- the defect it exists to catch."
        rc=1
    else
        echo "  [  ok  ] control rejected, so the filter is under test"
    fi
else
    # A control that cannot produce a restart has not been shown to differ for
    # the right reason, so it is not evidence either way.
    echo "  [ FAIL ] the control produced no restart, so it tests nothing"
    rc=1
fi

echo
if [ "$rc" = 0 ]; then
    echo "PASS: the SHTns path computes this model."
else
    echo "FAIL: see above."
fi
exit "$rc"
