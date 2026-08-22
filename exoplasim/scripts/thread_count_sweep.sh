#!/bin/bash
# CLIM-65: is sixteen threads the right count at every resolution?
# The arms and the verdict are exoplasim/notes/thread-count-by-resolution.md.
#
# The thread count is COMPILED IN -- `!$omp parallel num_threads(NPRO)` -- so a
# sweep needs one binary per count and OMP_NUM_THREADS does nothing. NPRO must
# divide NLAT, which at T42 (NLAT 64) allows 2, 4, 8, 16 and more.
#
# Each count is benched PAIRWISE against sixteen with bench_ab, so every
# comparison is interleaved against a run seconds away and inherits the
# discipline that harness already has. Wall time is the decision variable: the
# question is which configuration finishes first, not which scales best.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PKG="$REPO/vendor/exoplasim/exoplasim"
SCRATCH="${SCRATCH:-${TMPDIR:-/tmp}/thread_count_sweep}"
mkdir -p "$SCRATCH"
RES="${1:-T42}"
BED="${2:-bed_t42}"
COUNTS="${3:-2 4 8 16}"
low=$(echo "$RES" | tr 'A-Z' 'a-z')

for n in $COUNTS; do
    ( cd "$PKG" && ./compile.sh -j -p 8 -r "$RES" -v 10 -n "$n" ) \
        > "$SCRATCH/ts_build_${RES}_$n.log" 2>&1
    b="$PKG/plasim/run/most_plasim_${low}_l10_p${n}_omp.x"
    [ -f "$b" ] || { echo "build failed at $n threads"; exit 1; }
    cp "$b" "$SCRATCH/ts_${RES}_$n.x"
    echo "built $RES $n threads: $(sha256sum "$SCRATCH/ts_${RES}_$n.x" | cut -c1-12)"
done

d="$SCRATCH/ts_bed_$RES"
rm -rf "$d"; mkdir -p "$d"
cp -a "$REPO/exoplasim/bench/$BED/." "$d/"
rm -f "$d"/*.x "$d/plasim_status" "$d/Abort_Message"; rm -rf "$d"/perf_*
for n in $COUNTS; do cp "$SCRATCH/ts_${RES}_$n.x" "$d/t$n.x"; done

cd "$REPO" || exit 1
for n in $COUNTS; do
    [ "$n" = 16 ] && continue
    echo
    echo "===== $RES: $n threads against 16 ====="
    .venv/bin/python exoplasim/scripts/bench_ab.py \
        --bed "$d" --a "$d/t16.x" --b "$d/t$n.x" \
        --label-a "16thr" --label-b "${n}thr" \
        --a-launch omp --b-launch omp --ranks "$n" \
        --rounds 4 --out "$SCRATCH/ts_${RES}_$n.json" 2>&1 | tail -5
done
