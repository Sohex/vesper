#!/bin/bash
# Does the threaded build compute the same thing twice?
#
#   exoplasim/scripts/verify_shared_determinism.sh <bed> <res> <threads> [steps] [repeats]
#
# Worldbuilding frame: a correctness check on the Vesper climate model's
# threaded build. Nothing here is about the simulated planet.
#
# WHAT IS BEING CHECKED. The threaded build shares the spectral state between
# threads instead of giving each one a private copy, which is what removes the
# staging traffic that made threads slower than ranks. Sharing is safe only if
# every write is either the owning thread writing its own slice or a root-only
# write with a barrier either side. A write that is neither does not usually
# crash and does not usually give an obviously wrong answer: it gives a
# DIFFERENT answer each time, because it depends on which thread got there
# first. So the check is repetition, not comparison against a reference.
#
# WHY BIT IDENTICAL IS THE RIGHT BAR HERE, when compare_restarts.py exists
# precisely because decomposition changes move the last bits. Nothing is being
# decomposed differently between these runs. Same binary, same bed, same thread
# count, same fixed reduction order in mpsumsc. Every floating point operation
# is the same operation in the same order, so the only thing that can move a
# bit is a race. A tolerance here would be a way of not noticing one.
#
# THE CONTROL, and why it is not the obvious one. The obvious control is a
# deliberately racy build that must fail the comparison. That control cannot be
# built usefully here: the race this check exists for corrupts planetary
# vorticity on the first step, so the pre-fix build does not survive five steps,
# and a control that crashes has bypassed the comparison rather than failed it.
# The arrival-order reduction the plan proposed has the opposite problem -- a
# sha detects a single flipped bit with certainty, so a comparison of hashes has
# no sensitivity to demonstrate.
#
# What CAN go wrong, and what the third arm therefore tests, is that the hash
# stops tracking the computation: three runs sharing one output directory, or a
# restart written before the model touched it, would all agree for reasons that
# have nothing to do with the model being reproducible. So the third arm runs
# the same binary on a bed that differs ONLY in the random seed and requires a
# DIFFERENT hash. If it comes back identical, the first arm's agreement means
# nothing and the script says so.
set -euo pipefail

bed="$(cd "${1:?usage: verify_shared_determinism.sh <bed> <res> <threads> [steps] [repeats]}" && pwd)"
res="${2:?}"
threads="${3:?}"
steps="${4:-200}"
repeats="${5:-3}"

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
# shellcheck source=_bed_guard.sh
. "$HERE/_bed_guard.sh"

PKG="$REPO/vendor/exoplasim/exoplasim"
BUILD="$REPO/.venv/bin/python $REPO/exoplasim/scripts/build_model.py"
WORK="${TMPDIR:-/tmp}/verify_shared_determinism.$$"
low="$(echo "$res" | tr 'A-Z' 'a-z')"
name="most_plasim_${low}_l10_p${threads}_omp.x"

dirty="$(cd "$REPO" && git status --porcelain -- "vendor/exoplasim/exoplasim/plasim/src")"
if [ -n "$dirty" ]; then
    echo "refusing: the model source is not what is committed --" >&2
    echo "$dirty" >&2
    exit 1
fi

require_settled_bed "$bed"

rm -rf "$WORK"; mkdir -p "$WORK/ref"

build() {
    # -j is a FLAG and takes no argument: it selects the threaded build. The
    # thread count is -n, as it is for ranks. Passing the count to -j leaves
    # the thread count is explicit: build_model.py has no default for it, so
    # whatever p2 was lying in plasim/run untouched. The freshness check below
    # is here because that is not hypothetical: it happened, and a build from
    # before the fix was run under the name of the one after it.
    local stamp="$WORK/stamp"; : > "$stamp"
    ( $BUILD --res "$res" --ranks "$threads" --parmode omp ) \
        >"$WORK/build.log" 2>&1 || true
    [ -f "$PKG/plasim/run/$name" ] || {
        echo "build failed: no $name (see $WORK/build.log)" >&2; exit 1; }
    [ "$PKG/plasim/run/$name" -nt "$stamp" ] || {
        echo "build failed: $name is older than this build started, so" >&2
        echo "  build_model.py did not produce it. See $WORK/build.log" >&2; exit 1; }
    cp -f "$PKG/plasim/run/$name" "$WORK/ref/probe.x"
    echo "built  $(sha256sum "$WORK/ref/probe.x" | cut -c1-16)"
}

# tag [seed-line] -> prints the restart sha, or nothing if the run died
run_once() {
    local tag="$1" seedline="${2:-}"
    local d="$WORK/run_$tag"
    rm -rf "$d"; mkdir -p "$d"
    cp -a "$bed"/. "$d"/
    ( cd "$d"
      rm -f ./*.x plasim_status Abort_Message MOST_REST.* MOST_DIAG.*
      sed -i "s/^ *N_RUN_STEPS *=.*/ N_RUN_STEPS = $steps /" plasim_namelist
      [ -n "$seedline" ] && sed -i "s/^ *SEED *=.*/ $seedline/" plasim_namelist
      cp -f "$WORK/ref/probe.x" ./probe.x
      # Bound the threads. Unbound, the runtime places them differently run to
      # run, which changes nothing about the arithmetic but changes the timing,
      # and timing is the only thing a race is sensitive to. Binding makes a
      # passing run mean MORE, not less: it is the harder case for detection.
      export OMP_NUM_THREADS="$threads" OMP_PROC_BIND=close OMP_PLACES=cores
      export OMP_STACKSIZE=512M
      ulimit -s unlimited
      ./probe.x >run.log 2>&1 ) 2>/dev/null || true
    [ -f "$d/plasim_status" ] && sha256sum "$d/plasim_status" | cut -d' ' -f1
}

echo "res $res  threads $threads  steps $steps  repeats $repeats"
echo
build
echo

echo "repetition (must all be identical):"
first=""; n_ok=0; n_diff=0; n_fail=0
for i in $(seq 1 "$repeats"); do
    s="$(run_once "rep$i")"
    if [ -z "$s" ]; then
        n_fail=$((n_fail+1)); echo "  run $i: DID NOT RUN"; continue
    fi
    n_ok=$((n_ok+1))
    [ -z "$first" ] && first="$s"
    [ "$s" = "$first" ] || n_diff=$((n_diff+1))
    echo "  run $i: ${s:0:16}"
done
echo "  $n_ok of $repeats ran; $n_diff of those differed from the first"

echo
echo "seed sensitivity (must DIFFER from the runs above):"
alt="$(run_once alt " SEED = 5,11,13,19,23,31,37,43")"
if [ -z "$alt" ]; then
    echo "  DID NOT RUN"
else
    echo "  altered seed: ${alt:0:16}"
fi

echo
rc=0
if [ "$n_fail" -gt 0 ]; then
    echo "FAIL: a run did not complete, so the build has been shown neither"
    echo "      reproducible nor unreproducible. Nothing was measured."
    echo "      Read run.log under the work directory."
    rc=1
elif [ "$n_diff" -gt 0 ]; then
    echo "FAIL: the threaded build is not reproducible run to run. That is a"
    echo "      race on shared state; find it before measuring anything."
    rc=1
else
    echo "pass: $n_ok runs, one answer."
fi
if [ -z "$alt" ]; then
    echo "FAIL: the seed arm did not run, so nothing shows the hash still"
    echo "      tracks the computation."
    rc=1
elif [ "$alt" = "$first" ]; then
    echo "FAIL: a different random seed gave the SAME restart hash. The hash is"
    echo "      not tracking the computation, so the agreement above is an"
    echo "      artifact of the harness and means nothing."
    rc=1
else
    echo "pass: a different seed gives a different hash, so the hash tracks the"
    echo "      computation and the agreement above is about the model."
fi
echo "work kept at $WORK"
exit $rc
