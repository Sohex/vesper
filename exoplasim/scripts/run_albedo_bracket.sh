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
# Each case is prepared fresh, so each gets its own run. Run ids are UUIDs, so
# the four cannot collide whatever they share; what tells them apart afterwards
# is the physical fingerprint in `exoplasim/runs/INDEX.json`, whose `geography`
# digest covers the albedo file and so separates the two modes at one flux.
#
# HOW THIS SCRIPT KNOWS WHICH RUN TO EXTEND. It is told, by the process that
# created it. `run_exoplasim.py` announces `RUN_ID=<id>` on its own line before
# it does anything expensive, and that line is threaded straight into
# `continue_exoplasim.py --run`. The id is never reconstructed and never
# discovered: a UUID cannot be recomputed from the physics (CLAUDE.md rule 6),
# and searching INDEX.json for a run matching this configuration would be that
# derived name again with extra steps -- it would match every earlier bracket at
# the same flux on the same surface, which is precisely the collision the UUID
# exists to make impossible. Exactly one announced id continues; zero or more
# than one is a refusal, because continuing the wrong run is the failure all of
# this exists to prevent.
#
# `--self-test` exercises that threading against a stub interpreter, with no
# model and nothing written outside a temporary directory.

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

# Overridable so --self-test can substitute a stub. Everything else about the
# script is the real path.
PYTHON=${PYTHON:-python}
ORBITS=${ORBITS:-49}
LOG=${LOG:-/tmp/albedo_bracket.log}

run_case() {
    local mode=$1 flux=$2
    local out status run_id
    local ids=()
    out=$(mktemp "${TMPDIR:-/tmp}/albedo_bracket_case.XXXXXX") || return 1
    echo "=== ${mode} @ ${flux} S-Earth :: $(date +%H:%M:%S) ===" | tee -a "$LOG"

    "$PYTHON" exoplasim/scripts/build_surface_albedo.py --mode "$mode" >"$out" 2>&1
    status=$?
    cat "$out" >>"$LOG"
    if [ "$status" -ne 0 ]; then rm -f "$out"; return 1; fi

    # Captured rather than appended, so the id below is read from THIS
    # invocation's output and not from a shared log that four cases write into.
    : >"$out"
    "$PYTHON" exoplasim/scripts/run_exoplasim.py --flux-ratio "$flux" --run-years 1 >"$out" 2>&1
    status=$?
    cat "$out" >>"$LOG"
    # The id is announced before the prepare guards run, so a failed prepare can
    # leave a RUN_ID line behind with no run under it. Status first, always.
    if [ "$status" -ne 0 ]; then rm -f "$out"; return 1; fi

    mapfile -t ids < <(sed -n 's/^RUN_ID=//p' "$out")
    rm -f "$out"
    if [ "${#ids[@]}" -ne 1 ] || [ -z "${ids[0]}" ]; then
        echo "REFUSING ${mode} @ ${flux}: run_exoplasim.py announced ${#ids[@]}" \
             "RUN_ID lines, expected exactly one. Nothing was continued." | tee -a "$LOG"
        return 1
    fi
    run_id=${ids[0]}
    echo "    run ${run_id}" | tee -a "$LOG"

    # --flux-ratio is redundant for the model and is not redundant here: the
    # continuation refuses when it disagrees with the manifest of the run it was
    # handed, so it is a free assertion that the captured id is the run this
    # case prepared.
    "$PYTHON" exoplasim/scripts/continue_exoplasim.py --run "$run_id" \
        --flux-ratio "$flux" --orbits "$ORBITS" --purpose spinup >>"$LOG" 2>&1 || return 1
    echo "    done $(date +%H:%M:%S)" | tee -a "$LOG"
}

# Class 17: every case has an answer that is right rather than merely different,
# and each negative is paired with the positive that proves the setup was real.
# What is under test is the run-id threading, which is the whole of the
# interface, because a wrong id here continues someone else's climate.
self_test() {
    local tmp stub calls failures=0

    # name, RUN_ID lines the stub emits, its exit status, expected run_case
    # status, expected argument to continue_exoplasim.py --run.
    _case() {
        local name=$1 want_rc=$4 want_run=$5 rc got
        export STUB_RUN_IDS=$2 STUB_RUN_EXIT=$3
        : >"$calls"
        run_case lithology 0.90 >/dev/null 2>&1
        rc=$?
        got=$(sed -n 's/.*continue_exoplasim\.py --run \([^ ]*\).*/\1/p' "$calls")
        if [ "$rc" -ne "$want_rc" ]; then
            echo "FAIL ${name}: exit ${rc}, expected ${want_rc}"
            failures=$((failures + 1))
        elif [ "$got" != "$want_run" ]; then
            echo "FAIL ${name}: continued '${got}', expected '${want_run}'"
            failures=$((failures + 1))
        else
            echo "ok   ${name}"
        fi
    }

    tmp=$(mktemp -d) || return 1
    stub="$tmp/interpreter"
    calls="$tmp/calls"
    cat >"$stub" <<'STUB'
#!/usr/bin/env bash
# Stands in for the interpreter: records every invocation, and emits as many
# RUN_ID lines as the case under test asks for.
printf '%s\n' "$*" >>"$STUB_CALLS"
case "$1" in
  *run_exoplasim.py)
    i=0
    while [ "$i" -lt "$STUB_RUN_IDS" ]; do echo "RUN_ID=run_stub$i"; i=$((i + 1)); done
    exit "$STUB_RUN_EXIT" ;;
esac
exit 0
STUB
    chmod +x "$stub"
    PYTHON="$stub"
    ORBITS=1
    LOG="$tmp/log"
    export STUB_CALLS="$calls"

    # Positive: one announced id, and it is the one continued.
    _case "one id is threaded through"  1 0 0 "run_stub0"
    # Nothing announced, so there is no id to continue and nothing may be.
    _case "no id refuses"               0 0 1 ""
    # Two announced: an ambiguous match is a refusal and not a choice.
    _case "two ids refuse"              2 0 1 ""
    # Announced, but the prepare failed. The id exists and the run does not.
    _case "failed prepare refuses"      1 1 1 ""

    rm -rf "$tmp"
    unset STUB_CALLS STUB_RUN_IDS STUB_RUN_EXIT
    if [ "$failures" -ne 0 ]; then
        echo "${failures} self-test failure(s)"
        return 1
    fi
    echo "self-test passed"
}

if [ "${1:-}" = "--self-test" ]; then
    self_test
    exit $?
fi

for flux in ${FLUXES:-0.90 0.95}; do
    for mode in lithology vegetated; do
        if ! run_case "$mode" "$flux"; then
            echo "FAILED: ${mode} @ ${flux}; see $LOG" | tee -a "$LOG"
            exit 1
        fi
    done
done

echo "all four cases complete $(date +%H:%M:%S)" | tee -a "$LOG"
