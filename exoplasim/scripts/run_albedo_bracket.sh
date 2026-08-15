#!/usr/bin/env bash
# Albedo-endmember bracket: two land-surface assumptions at two stellar fluxes.
#
# The question is whether this world has one stable climate or more than one.
# Bare rock and fully vegetated are the two physical endmembers, 0.315 and 0.197
# land-mean albedo, about 15 W/m2 apart in absorbed flux, which is comparable to
# the entire 0.85-to-0.95 stellar sweep. If both endmembers converge on a similar
# climate the vegetation feedback is a correction; if they diverge, the world is
# bistable and the state has to be chosen deliberately.
#
# Two fluxes because the old flux calibration is void: lithology albedo and the
# K2 spectrum together shift absorbed flux by about -9 W/m2 at fixed stellar
# input, so the transition has moved and we do not know where to.
#
# Each run is prepared fresh. The run directory name carries the surface-input
# digest, so the four cannot collide.

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

ORBITS=${ORBITS:-49}
LOG=${LOG:-/tmp/albedo_bracket.log}

run_case() {
    local mode=$1 flux=$2
    echo "=== ${mode} @ ${flux} S-Earth :: $(date +%H:%M:%S) ===" | tee -a "$LOG"
    python exoplasim/scripts/build_surface_albedo.py --mode "$mode" >>"$LOG" 2>&1 || return 1
    python exoplasim/scripts/run_exoplasim.py --flux-ratio "$flux" --run-years 1 >>"$LOG" 2>&1 || return 1
    python exoplasim/scripts/continue_exoplasim.py --flux-ratio "$flux" --orbits "$ORBITS" >>"$LOG" 2>&1 || return 1
    echo "    done $(date +%H:%M:%S)" | tee -a "$LOG"
}

for flux in ${FLUXES:-0.90 0.95}; do
    for mode in lithology vegetated; do
        if ! run_case "$mode" "$flux"; then
            echo "FAILED: ${mode} @ ${flux}; see $LOG" | tee -a "$LOG"
            exit 1
        fi
    done
done

echo "all four cases complete $(date +%H:%M:%S)" | tee -a "$LOG"
