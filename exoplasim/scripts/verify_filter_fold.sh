#!/bin/bash
# Prove the filter fold is exact, against something that can fail.
#
#   exoplasim/scripts/verify_filter_fold.sh <bed> <unpatched.x> <patched.x> <ranks>
#
# Worldbuilding frame: a correctness check on the Vesper climate model's
# spectral transform. Nothing here is about the simulated planet.
#
# THE TEST. Folding skgpsp/skspgp into the weight matrices replaces
# `q * field * filter` with `(q * filter) * field`. When the filters are
# DISABLED they are exactly 1.0, and multiplying an IEEE double by 1.0 is the
# identity for every finite value. So with NGPTFILTER=0 and NSPVFILTER=0 the
# patched model must reproduce the unpatched one BIT FOR BIT.
#
# That is a test with a right answer rather than a comparison that can only
# differ: if the restart shas disagree, the fold is indexing the wrong mode and
# no amount of "the climate looks similar" rescues it. The obvious alternative
# check -- run both with the filter ON and see whether the fields look close --
# cannot fail informatively, because reassociation makes them differ anyway.
#
# The filters-ON case is a separate question and is NOT settled here: there the
# shas legitimately differ in the last bits and the claim is only that the
# change is a reassociation, which this script does not attempt to prove.
set -euo pipefail

bed="$(cd "${1:?usage: verify_filter_fold.sh <bed> <unpatched.x> <patched.x> [ranks] [steps]}" && pwd)"
unpatched="$(readlink -f "${2:?}")"
patched="$(readlink -f "${3:?}")"
ranks="${4:-16}"
steps="${5:-20}"

work="$bed/../_foldcheck"
rm -rf "$work"; mkdir -p "$work"
cp -a "$bed"/. "$work"/
cd "$work"
rm -f MOST_REST.* MOST_DIAG.* plasim_status Abort_Message

# Disable both filters. This is the whole basis of the test: with them on, the
# fold is a reassociation and the shas are expected to differ.
sed -i 's/^ *NGPTFILTER *=.*/ NGPTFILTER = 0 /; s/^ *NSPVFILTER *=.*/ NSPVFILTER = 0 /' plasim_namelist

# Short, and that is not laziness. The filter is load-bearing for stability:
# with it off, an unpatched T127 cold start reaches a NaN in lwr_ within a few
# hundred steps and -ffpe-trap kills it. One timestep already exercises every
# transform routine in the model, so a few tens of steps prove the indexing
# while staying inside the window where filters-off still integrates.
sed -i "s/^ *N_RUN_STEPS *=.*/ N_RUN_STEPS = $steps /" plasim_namelist
grep -E "NGPTFILTER|NSPVFILTER|N_RUN_STEPS" plasim_namelist

run_one() {
    local exe="$1" tag="$2"
    rm -f plasim_status Abort_Message
    cp -f "$exe" ./probe.x
    if ! mpiexec -np "$ranks" ./probe.x >run.log 2>&1; then
        echo "$tag: the MODEL failed, not the comparison. Last lines:" >&2
        grep -v "Fortran runtime warning\|^At line" run.log | tail -6 >&2
        echo "  (if this is a floating point exception with the filters off, the" >&2
        echo "   run is too long for an unfiltered integration -- lower <steps>.)" >&2
        exit 1
    fi
    if [ -f Abort_Message ]; then echo "$tag ABORTED" >&2; exit 1; fi
    [ -f plasim_status ] || { echo "$tag produced no plasim_status" >&2; exit 1; }
    sha256sum plasim_status | cut -d' ' -f1
    rm -f probe.x
}

echo
echo "running unpatched ..."; a=$(run_one "$unpatched" unpatched); echo "  $a"
echo "running patched   ..."; b=$(run_one "$patched"   patched);   echo "  $b"

echo
if [ "$a" = "$b" ]; then
    echo "PASS: filters off, patched and unpatched restarts are BIT IDENTICAL."
    echo "      The fold indexes the same spectral mode the transforms do."
else
    echo "FAIL: restarts differ with the filters DISABLED, where the fold"
    echo "      multiplies by exactly 1.0 and must be the identity."
    echo "      The mode indexing is wrong -- legini counts m,n from 0 and the"
    echo "      transforms from 1, so check the skgpsp(n+1) offset."
    exit 1
fi
