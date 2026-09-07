#!/usr/bin/env bash
# The flux-to-kelvin bracket: two runs differing in NOTHING but the stellar flux.
#
# `lib/sensitivity.py` declares one sensitivity for the whole project, dT/df near
# the baseline flux, and it is a secant between two converged points. This script
# buys that pair. It exists because the pair has to be RE-bought whenever the
# staged surface moves: `_verify_bracket` refuses two arms whose `geography`
# digests differ, and that digest covers every staged surface field, so a change
# to any surface builder makes the standing pair unusable rather than merely old.
#
# WHY BOTH ARMS ARE BOUGHT TOGETHER AND NEVER ONE OF THEM. The secant's arms have
# to be one INSTRUMENT and one WINDOW. The instrument is the compiled model, the
# staged surface, the I/O regime and the thread count; the window is the orbit
# range the asymptote is fitted over. A slope taken across a change in any of
# them measures that change as well as the flux. Reusing a standing run as one
# arm looks cheaper and silently spends the measurement: the model source, the
# surface fields and the assessor's default window have all moved under this
# project inside a single week, and none of those movements is visible in a
# convergence report.
#
# THE STRUCTURE IS IDENTICAL ON BOTH ARMS BY CONSTRUCTION.
#
#   orbits 0..SPINUP-1        low I/O, purpose spinup
#   orbits SPINUP..TOTAL-1    clean I/O, purpose post_equilibrium_climatology
#   assessed window           the last WINDOW orbits, inside the clean block
#
# so the two reports cover the same absolute orbit indices in the same regime.
# The clean block is longer than the window on purpose: a window that reached
# back to the join would span two instruments and `segments.py` refuses it.
#
# HOW THIS SCRIPT KNOWS WHICH RUN TO EXTEND. It is told. `run_exoplasim.py`
# announces `RUN_ID=<id>` on its own line before it does anything expensive, and
# that line is threaded into every later command for the same arm. A UUID cannot
# be recomputed from the physics (CLAUDE.md rule 6) and searching INDEX.json for
# a run matching this configuration would match every earlier bracket at the same
# flux on the same surface. Exactly one announced id continues; zero or more than
# one is a refusal.
#
# THE FIRST ORBIT IS BOUGHT SEPARATELY so the run directory can be relocated
# before anything records a path inside it. In a worktree `exoplasim/runs/` is a
# real directory of per-file links, so a run created there lives only in the
# worktree and dies with it, while the declaration that cites it is tracked and
# does not. `RUNS_HOME`, when set, is the directory the payload is moved into,
# with a symlink left behind; the run then resolves to the same path from every
# tree. Leave it unset in the main checkout, where the run is already home.
#
# It runs nothing itself: wrap the whole invocation in `qrun`.
#
#     qrun -- \
#         exoplasim/scripts/run_flux_bracket.sh
#
# `--self-test` exercises the id threading and the relocation against a stub
# interpreter, with no model and nothing written outside a temporary directory.

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

PYTHON=${PYTHON:-.venv/bin/python}
# The two fluxes. COLD is the design flux `config/planet.yaml` declares; WARM is
# one unit-round step above it, wide enough that the response is far larger than
# either fit's uncertainty and narrow enough to stay on one side of the ice
# transition. Both are arguments so a re-bracket can move them together.
COLD=${COLD:-0.945}
WARM=${WARM:-1.000}
# Fixed BEFORE any result is seen. SPINUP is the cold-start convergence length
# this project has recorded for T21; WINDOW is assess_convergence.py's derived
# default, the shortest window at which the offset criterion discriminates
# rather than flips; CLEAN exceeds WINDOW so the window clears the I/O join.
SPINUP=${SPINUP:-70}
CLEAN=${CLEAN:-60}
WINDOW=${WINDOW:-55}
LOG=${LOG:-/tmp/flux_bracket.log}
RUNS_HOME=${RUNS_HOME:-}

relocate() {
    # Move a freshly prepared run out of this tree and link it back, so the
    # payload outlives the worktree the declaration was written in.
    local run_id=$1 here there
    [ -n "$RUNS_HOME" ] || return 0
    here="exoplasim/runs/${run_id}"
    there="${RUNS_HOME%/}/${run_id}"
    if [ -L "$here" ]; then return 0; fi
    if [ -e "$there" ]; then
        echo "REFUSING to relocate ${run_id}: ${there} already exists" | tee -a "$LOG"
        return 1
    fi
    mv "$here" "$there" || return 1
    ln -s "$there" "$here" || return 1
    # The prepare wrote its own absolute paths into the manifest before the move.
    # They name the same files and would otherwise name a directory that no
    # longer exists, which reads as a lost artifact rather than a moved one.
    if [ -f "${there}/run_manifest.json" ]; then
        "$PYTHON" - "$there" "$run_id" <<'PY' || return 1
import json, sys
from pathlib import Path
there, run_id = sys.argv[1], sys.argv[2]
path = Path(there) / "run_manifest.json"
text = path.read_text(encoding="utf-8")
old = str(Path.cwd() / "exoplasim" / "runs" / run_id)
if old != there and old in text:
    path.write_text(text.replace(old, there), encoding="utf-8")
    print(f"    manifest paths repointed from {old}")
PY
    fi
    echo "    payload at ${there}, linked from ${here}" | tee -a "$LOG"
}

run_arm() {
    local flux=$1
    local out status run_id
    local ids=()
    out=$(mktemp "${TMPDIR:-/tmp}/flux_bracket_arm.XXXXXX") || return 1
    echo "=== arm @ ${flux} S-Earth :: $(date +%H:%M:%S) ===" | tee -a "$LOG"

    "$PYTHON" exoplasim/scripts/run_exoplasim.py --flux-ratio "$flux" \
        --purpose spinup --run-years 1 >"$out" 2>&1
    status=$?
    cat "$out" >>"$LOG"
    # The id is announced before the prepare guards run, so a failed prepare can
    # leave a RUN_ID line behind with no run under it. Status first, always.
    if [ "$status" -ne 0 ]; then rm -f "$out"; return 1; fi

    mapfile -t ids < <(sed -n 's/^RUN_ID=//p' "$out")
    rm -f "$out"
    if [ "${#ids[@]}" -ne 1 ] || [ -z "${ids[0]}" ]; then
        echo "REFUSING arm @ ${flux}: run_exoplasim.py announced ${#ids[@]}" \
             "RUN_ID lines, expected exactly one. Nothing was continued." | tee -a "$LOG"
        return 1
    fi
    run_id=${ids[0]}
    echo "ARM_RUN_ID=${run_id}" | tee -a "$LOG"
    relocate "$run_id" || return 1

    # --flux-ratio is redundant for the model and is not redundant here: the
    # continuation refuses when it disagrees with the manifest of the run it was
    # handed, so it is a free assertion that the captured id is this arm's.
    "$PYTHON" exoplasim/scripts/continue_exoplasim.py --run "$run_id" \
        --flux-ratio "$flux" --orbits "$((SPINUP - 1))" --purpose spinup \
        --low-io >>"$LOG" 2>&1 || return 1
    # The cutoff a post-equilibrium segment is refused without. It is a claim
    # about where this arm has got to, not yet the bracket's reading.
    "$PYTHON" exoplasim/scripts/assess_convergence.py "exoplasim/runs/${run_id}" \
        >>"$LOG" 2>&1 || return 1
    "$PYTHON" exoplasim/scripts/continue_exoplasim.py --run "$run_id" \
        --flux-ratio "$flux" --orbits "$CLEAN" \
        --purpose post_equilibrium_climatology --clean-io >>"$LOG" 2>&1 || return 1
    # The bracket's reading: the last WINDOW orbits, wholly inside the clean
    # block, the same indices on both arms.
    "$PYTHON" exoplasim/scripts/assess_convergence.py "exoplasim/runs/${run_id}" \
        --window "$WINDOW" >>"$LOG" 2>&1 || return 1
    echo "    done $(date +%H:%M:%S)" | tee -a "$LOG"
}

# Class 17: what is under test is the id threading and the relocation, because a
# wrong id continues someone else's climate and a run left in a worktree is
# evidence that evaporates. Each negative is paired with the positive that shows
# the setup was real.
self_test() {
    local tmp stub calls failures=0 home

    _case() {
        local name=$1 want_rc=$4 want_run=$5 rc got
        export STUB_RUN_IDS=$2 STUB_RUN_EXIT=$3
        : >"$calls"
        run_arm 0.945 >/dev/null 2>&1
        rc=$?
        got=$(sed -n 's/.*continue_exoplasim\.py --run \([^ ]*\).*/\1/p' "$calls" | head -1)
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
printf '%s\n' "$*" >>"$STUB_CALLS"
case "$1" in
  *run_exoplasim.py)
    i=0
    while [ "$i" -lt "$STUB_RUN_IDS" ]; do
        mkdir -p "exoplasim/runs/run_stub$i"
        echo "RUN_ID=run_stub$i"
        i=$((i + 1))
    done
    exit "$STUB_RUN_EXIT" ;;
esac
exit 0
STUB
    chmod +x "$stub"
    PYTHON="$stub"
    SPINUP=2
    CLEAN=1
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

    # The relocation moves the payload and leaves a link that resolves to it.
    home="$tmp/home"
    mkdir -p "$home"
    RUNS_HOME="$home"
    rm -rf exoplasim/runs/run_stub0
    _case "relocated arm is still threaded" 1 0 0 "run_stub0"
    if [ ! -d "$home/run_stub0" ]; then
        echo "FAIL relocation: payload is not in RUNS_HOME"
        failures=$((failures + 1))
    elif [ ! -L "exoplasim/runs/run_stub0" ]; then
        echo "FAIL relocation: no link left behind"
        failures=$((failures + 1))
    else
        echo "ok   relocation moves the payload and links it back"
    fi
    # An occupied destination is a collision, and a collision is a refusal
    # rather than an overwrite: the payload already there is somebody's run.
    rm -f exoplasim/runs/run_stub0
    _case "occupied destination refuses" 1 0 1 ""

    rm -rf exoplasim/runs/run_stub0 exoplasim/runs/run_stub1
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

for flux in "$COLD" "$WARM"; do
    if ! run_arm "$flux"; then
        echo "FAILED: arm @ ${flux}; see $LOG" | tee -a "$LOG"
        exit 1
    fi
done

"$PYTHON" exoplasim/scripts/index_runs.py >>"$LOG" 2>&1 || exit 1
echo "both arms complete $(date +%H:%M:%S)" | tee -a "$LOG"
