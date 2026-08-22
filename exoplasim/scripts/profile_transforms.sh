#!/bin/bash
# Profile a bed under perf, every rank, and leave the samples for scoring.
#
#   exoplasim/scripts/profile_transforms.sh <bed_dir> <binary> <ranks> [repeat]
#   PERF_LAUNCH=omp ... <ranks is the thread count>   for the threaded build
#
# Worldbuilding frame: a COMPUTE measurement of the Vesper climate model on this
# desktop. Nothing it produces is about the simulated planet.
#
# Method, and the two halves of it that are not optional. The machine must be
# quiet -- `exoplasim/notes/rank-layout-benchmark.md` -- AND a warm-up run must
# precede the first recorded one, because whichever arm runs first after an idle
# stretch gets the boost clock and the rest run warm
# (`notes/audits/aocl-and-model-build-flags.md`). The warm-up here is unrecorded
# and its samples are discarded.
#
# No `--map-by` and no `--bind-to`: production launches with Open MPI's defaults
# and the profile has to be of production. Adding pinning would measure a
# configuration nobody runs.
set -euo pipefail

bed="$(cd "${1:?usage: profile_transforms.sh <bed_dir> <binary> <ranks> [repeat]}" && pwd)"
binary="${2:?binary name inside the bed}"
ranks="${3:?rank count}"
repeat="${4:-1}"

[ -x "$bed/$binary" ] || { echo "no executable $bed/$binary" >&2; exit 1; }
[ -f "$bed/bed_manifest.json" ] || { echo "$bed is not a bed (no bed_manifest.json)" >&2; exit 1; }

here="$(cd "$(dirname "$0")" && pwd)"
export PERF_EVENT="${PERF_EVENT:-cpu-clock}"
export PERF_FREQ="${PERF_FREQ:-997}"

# PERF_LAUNCH=omp profiles the THREADED build, which is one process rather than
# one per rank -- so there is no filename collision to solve and no mpiexec to
# lose the ring buffer inside. The sample file is still named perf.rank00.data,
# because that is what score_transform_profile.py reads and a threaded profile
# is the same measurement of the same model.
#
# The two exports and the ulimit are not tuning. libgomp does not bind by
# default, and unbound threads on this processor land on SMT siblings and across
# both dies, which changes what any profile of them means. The master thread
# runs on the process stack rather than on OMP_STACKSIZE, and at T127 it
# overruns the default.
export PERF_LAUNCH="${PERF_LAUNCH:-mpi}"
run_model () {
    if [ "$PERF_LAUNCH" = "omp" ]; then
        OMP_PROC_BIND=close OMP_PLACES=cores OMP_STACKSIZE=512M \
            bash -c "ulimit -s unlimited; exec $*"
    else
        mpiexec -np "$ranks" $*
    fi
}

cd "$bed"

echo "== warm-up, unrecorded =="
run_model "./$binary" >/dev/null 2>&1 || {
    echo "warm-up failed; check $bed for Abort_Message" >&2; exit 1; }
[ -e Abort_Message ] && { echo "model aborted during warm-up" >&2; exit 1; }

for r in $(seq 1 "$repeat"); do
    # The EVENT is part of the directory name. It was not, and a cache-miss
    # pass then overwrote a cpu-clock pass of the same index in the same bed,
    # destroying it: two profiles of the same bed are different measurements and
    # must not collide on a filename.
    outdir="$bed/perf_${PERF_EVENT}_$(printf '%02d' "$r")"
    rm -rf "$outdir"; mkdir -p "$outdir"
    export PERF_OUTDIR="$outdir"
    echo "== recorded pass $r of $repeat -> $(basename "$outdir") =="
    start=$(date +%s.%N)
    if [ "$PERF_LAUNCH" = "omp" ]; then
        run_model "perf record --quiet --freq $PERF_FREQ --event $PERF_EVENT --call-graph fp --mmap-pages 32 --output $outdir/perf.rank00.data ./$binary"
    else
        mpiexec -np "$ranks" "$here/perf_rank.sh" "./$binary"
    fi
    end=$(date +%s.%N)
    [ -e Abort_Message ] && { echo "model aborted in pass $r" >&2; exit 1; }
    # The restart sha is the numerics tripwire: every pass integrates the same
    # steps from the same state, so a pass that disagrees with its siblings is
    # not a slow measurement, it is a different model.
    sha=$(sha256sum plasim_status 2>/dev/null | cut -c1-16 || echo "(no plasim_status)")
    printf '{"pass":%d,"wall_s":%.3f,"ranks":%d,"launch":"%s","binary":"%s","event":"%s","freq":%s,"status_sha256_16":"%s"}\n' \
        "$r" "$(echo "$end - $start" | bc)" "$ranks" "$PERF_LAUNCH" "$binary" "$PERF_EVENT" "$PERF_FREQ" "$sha" \
        > "$outdir/pass.json"
    cat "$outdir/pass.json"
done

echo
echo "score with:"
echo "  .venv/bin/python exoplasim/scripts/score_transform_profile.py $bed/perf_${PERF_EVENT}_*"
