#!/bin/bash
# Shared precondition for any check that compares two model runs.
#
# Worldbuilding frame: a guard on the Vesper climate model's verification beds.
# Nothing here is about the simulated planet.
#
# WHAT THIS EXISTS FOR. A bed that starts cold with KICK non-zero has white
# noise added to log surface pressure before the first step, and plasim.f90
# seeds that generator from the CLOCK unless the namelist sets SEED. Two runs of
# ONE binary on such a bed differ from step one, so every comparison built on it
# reports a broken model no matter what is being compared -- and reports it in
# the loud, unambiguous way that invites you to believe it. It cost a correct
# transform rewrite an afternoon of being disbelieved: the check said 94 of 199
# records were wrong at one step, and what was wrong was the bed.
#
# Two ways to be settled, and either is enough:
#   - the bed carries a restart, so there is no cold start and no kick, or
#   - the bed fixes SEED, so the kick is the same kick every time.
#
# This is a precondition and not a test. It cannot tell you the bed is
# reproducible; it can only tell you the one way it is known to fail. The
# comparison scripts that source this still run a binary against ITSELF first,
# because that is the check with an answer.
require_settled_bed() {
    local bed="$1"
    local nl="$bed/plasim_namelist"
    [ -f "$nl" ] || { echo "no plasim_namelist in $bed" >&2; return 1; }

    # A restart means the run does not start from rest, so noise() is not called.
    if [ -f "$bed/plasim_restart" ] || [ -f "$bed/plasim_status" ]; then
        return 0
    fi

    local kick seed
    kick="$(grep -iE "^ *KICK *=" "$nl" | head -1 | tr -dc '0-9')"
    if [ -z "$kick" ] || [ "$kick" = "0" ]; then
        return 0
    fi

    seed="$(grep -iE "^ *SEED *=" "$nl" | head -1 | tr -dc '0-9')"
    if [ -n "$seed" ] && [ "$seed" != "0" ]; then
        return 0
    fi

    echo "refusing: $bed is a cold start with KICK=$kick and no fixed SEED." >&2
    echo "  plasim.f90 seeds the kick from the clock when seed(1) is zero, so two" >&2
    echo "  runs of one binary on this bed do not agree and nothing compared on it" >&2
    echo "  means anything. Add a SEED line to plasim_namelist, or give the bed a" >&2
    echo "  restart:" >&2
    echo "    SEED = 17,29,41,53,67,79,91,103" >&2
    return 1
}


# WHY THE GRID IS A PRECONDITION TOO. Every gate that sources this takes the
# rung as an argument and builds its arms at that rung, but nothing looked at
# what the BED is. A T42 bed handed "T85" builds T85 binaries, hands them N064
# surface files, and dies inside the model on an SRA header -- which every one
# of these gates reports as "an arm produced no restart", the same line it
# prints for a genuine crash of the thing under test. The bed says its own grid
# in the names of the files the model reads, so this reads it there.
#
# The executable the bed was made with says it too, and is checked when it is
# present: `most_plasim_t42_l10_p16.x` and an N085 surface set cannot both be
# right, and a bed assembled from two runs is worth catching before a build.
require_bed_grid() {
    local bed="$1" res="$2"
    local here repo want got binary

    here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    repo="$(cd "$here/../.." && pwd)"
    local py="$repo/.venv/bin/python"; [ -x "$py" ] || py=python3
    want="$("$py" -c "import sys; sys.path.insert(0, '$repo/lib'); import rungs; print(rungs.geometry('$res')[0])")" || {
        echo "unknown resolution $res: it is not a rung in lib/rungs.py" >&2
        return 1; }

    got=""
    for f in "$bed"/N???_surf_*.sra; do
        [ -e "$f" ] || continue
        got="$(basename "$f")"; got="${got#N}"; got="${got%%_*}"
        got="$((10#$got))"
        break
    done
    if [ -z "$got" ]; then
        echo "refusing: no N???_surf_*.sra in $bed, so its grid cannot be read" >&2
        echo "  and nothing here can tell whether it is a $res bed." >&2
        return 1
    fi
    if [ "$got" != "$want" ]; then
        echo "refusing: $bed carries N$(printf '%03d' "$got") surface files and $res is $want latitudes." >&2
        echo "  The arms would be built at $res, handed this bed's inputs, and die" >&2
        echo "  inside the model on an SRA header -- which reads here as 'an arm" >&2
        echo "  produced no restart' and says nothing about the thing under test." >&2
        return 1
    fi

    for binary in "$bed"/most_plasim_*.x; do
        [ -e "$binary" ] || continue
        case "$(basename "$binary")" in
          most_plasim_"$(echo "$res" | tr 'A-Z' 'a-z')"_*) ;;
          *) echo "refusing: $bed was made with $(basename "$binary") and $res was asked for." >&2
             echo "  The surface files say $want latitudes, so the bed is assembled" >&2
             echo "  from two runs and one of the two is wrong." >&2
             return 1 ;;
        esac
        break
    done
    return 0
}
