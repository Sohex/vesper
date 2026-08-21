#!/bin/sh
# Wrap one MPI rank in `perf record`, one output file per rank.
#
# `mpiexec -np 16 perf record -o perf.data ...` gives sixteen ranks one filename
# and fifteen of them lose. perf has no rank variable to interpolate, so the
# name has to be built inside the rank, which is what this exists for. Open MPI
# sets OMPI_COMM_WORLD_RANK; MPICH sets PMI_RANK, accepted here so the bed is
# not tied to one launcher.
#
# Driven by profile_transforms.sh, which sets PERF_OUTDIR, PERF_EVENT and
# PERF_FREQ. Not useful on its own.
set -eu

rank="${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-0}}"
outdir="${PERF_OUTDIR:?PERF_OUTDIR not set}"
event="${PERF_EVENT:-cpu-clock}"
freq="${PERF_FREQ:-997}"
# perf's default ring buffer does not fit inside an mpiexec'd rank on this host.
# `perf_event_mlock_kb` is 516 and RLIMIT_MEMLOCK is 8 MB, and Open MPI's
# transport has already locked pages against the latter by the time perf asks,
# so the mapping is refused -- "Permission error mapping pages", exit 255, and
# a zero-byte perf.data per rank. It is not a permissions problem despite the
# wording, and it does NOT reproduce outside mpiexec, which is what makes it
# worth a comment: perf works when you test it by hand and fails in the job.
# 32 pages is 128 KB per rank, verified at 16 ranks with zero lost samples.
mmap="${PERF_MMAP_PAGES:-32}"

# No --call-graph: flat symbol attribution is the whole question, callers are
# not, and this perf rejects `--call-graph=none` rather than treating it as the
# default it already is.
exec perf record --quiet \
    --event "$event" --freq "$freq" --mmap-pages "$mmap" \
    --output "$outdir/perf.rank$(printf '%02d' "$rank").data" \
    -- "$@"
