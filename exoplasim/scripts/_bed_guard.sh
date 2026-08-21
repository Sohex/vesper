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
