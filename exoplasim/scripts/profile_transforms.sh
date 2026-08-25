#!/bin/bash
# Profile a bed under perf and leave the samples for scoring.
#
#   exoplasim/scripts/profile_transforms.sh <bed_dir> <binary> <threads> [repeat]
#
# Worldbuilding frame: a COMPUTE measurement of the Vesper climate model on this
# desktop. Nothing it produces is about the simulated planet.
#
# ONE PROCESS, ONE SAMPLE FILE. The model is a threaded build: <threads> is the
# thread count the binary was compiled for, the whole team lives in one address
# space, and one `perf record` around it sees every thread. The file is still
# named perf.rank00.data because that is what score_transform_profile.py reads.
#
# Method, and the two halves of it that are not optional. The machine must be
# quiet -- `exoplasim/notes/rank-layout-benchmark.md` -- AND a warm-up run must
# precede the first recorded one, because whichever arm runs first after an idle
# stretch gets the boost clock and the rest run warm
# (`notes/audits/aocl-and-model-build-flags.md`). The warm-up here is unrecorded
# and its samples are discarded.
#
# THE LAUNCH LINE IS PRODUCTION'S. `exoplasim/__init__.py` starts the model as
# `OMP_NUM_THREADS=n OMP_PLACES=cores OMP_PROC_BIND=close ./<binary>`, and every
# verification shell in this directory exports the same three, so the profile is
# of the placement the model actually runs under. They are not tuning and not
# optional: libgomp does not bind by default, and unbound threads on this
# processor land on SMT siblings and across both dies, which changes what any
# profile of them means. OMP_STACKSIZE sizes the NON-MASTER threads only; the
# master runs on the process stack, which at T127 overruns the default, hence
# the ulimit.
set -euo pipefail

bed="$(cd "${1:?usage: profile_transforms.sh <bed_dir> <binary> <threads> [repeat]}" && pwd)"
binary="${2:?binary name inside the bed}"
threads="${3:?thread count}"
repeat="${4:-1}"

[ -x "$bed/$binary" ] || { echo "no executable $bed/$binary" >&2; exit 1; }
[ -f "$bed/bed_manifest.json" ] || { echo "$bed is not a bed (no bed_manifest.json)" >&2; exit 1; }

export PERF_EVENT="${PERF_EVENT:-cpu-clock}"
export PERF_FREQ="${PERF_FREQ:-997}"
# PERF_MMAP_PAGES IS AN INPUT, because score_transform_profile.py's advice on a
# lost-sample count is to re-record with a larger one and that has to be
# followable. 32 pages is 128 KB and was verified at sixteen ranks with zero
# lost samples at 997 Hz. The constraint that forced the number is host-specific
# and is recorded because it recurs wherever perf runs inside something that has
# already locked pages: perf's DEFAULT ring buffer was refused outright when the
# model ran under mpiexec here -- `perf_event_mlock_kb` is 516 and RLIMIT_MEMLOCK
# is 8 MB, and Open MPI's transport had locked pages against the latter before
# perf asked. It reported "Permission error mapping pages", exited 255 and left a
# zero-byte perf.data, which is not a permissions problem despite the wording,
# and it did not reproduce outside the job. One process locks nothing first, so
# the default would serve here; 32 is kept because it is the size this workload
# has been sampled at.
export PERF_MMAP_PAGES="${PERF_MMAP_PAGES:-32}"

run_model () {
    OMP_NUM_THREADS="$threads" OMP_PROC_BIND=close OMP_PLACES=cores \
        OMP_STACKSIZE=512M bash -c "ulimit -s unlimited; exec $*"
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
    run_model "perf record --quiet --freq $PERF_FREQ --event $PERF_EVENT --call-graph fp --mmap-pages $PERF_MMAP_PAGES --output $outdir/perf.rank00.data ./$binary"
    end=$(date +%s.%N)
    [ -e Abort_Message ] && { echo "model aborted in pass $r" >&2; exit 1; }
    # The restart sha is the numerics tripwire: every pass integrates the same
    # steps from the same state, so a pass that disagrees with its siblings is
    # not a slow measurement, it is a different model.
    sha=$(sha256sum plasim_status 2>/dev/null | cut -c1-16 || echo "(no plasim_status)")
    printf '{"pass":%d,"wall_s":%.3f,"threads":%d,"binary":"%s","event":"%s","freq":%s,"mmap_pages":%s,"status_sha256_16":"%s"}\n' \
        "$r" "$(echo "$end - $start" | bc)" "$threads" "$binary" "$PERF_EVENT" "$PERF_FREQ" "$PERF_MMAP_PAGES" "$sha" \
        > "$outdir/pass.json"
    cat "$outdir/pass.json"
done

echo
echo "score with:"
echo "  .venv/bin/python exoplasim/scripts/score_transform_profile.py $bed/perf_${PERF_EVENT}_*"
