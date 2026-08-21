#!/bin/bash
# Prove the filter fold applies each mode's OWN filter value, with a test that
# can fail -- including a deliberately wrong arm that must fail.
#
#   exoplasim/bench/verify_fold_indexing.sh
#
# WHY THIS EXISTS. The earlier test ran with the filters DISABLED, where every
# coefficient is exactly 1.0. Multiplying by 1.0 is the identity wherever you
# put it, so that test passes whether the fold indexes skspgp(n+1), skspgp(n),
# or anything else. It proved the fold did not corrupt the transform. It did NOT
# prove the mode indexing, and the PR claimed it did.
#
# THE TEST. Make the filter a function of the wavenumber whose values are exact
# POWERS OF TWO. Scaling an IEEE double by 2^k is exact, and rounding commutes
# with it: round(x * 2^k) == round(x) * 2^k. So
#     (qi * sk) * sp   ==   (qi * sp) * sk     bit for bit, when sk = 2^k
# and the folded and unfolded builds must agree EXACTLY -- while a fold that
# picked the wrong n would apply a different power of two to most modes and
# diverge immediately.
#
# Three arms, so the test has a negative control:
#   unfolded   stock
#   folded     skspgp(n+1) / skgpsp(n+1), the shipped fold
#   wrong      skspgp(n)   / skgpsp(n),   deliberately off by one
# PASS requires unfolded == folded AND wrong != folded. An arm that cannot fail
# is not a test.
set -euo pipefail
WT=/home/cfutro/docs/world/.claude/worktrees/spectral-optimisation
PKG="$WT/vendor/exoplasim/exoplasim"
SRC="$PKG/plasim/src/legmod.f90"
SP=/tmp/claude-1000/-home-cfutro-docs-world/ba0fe715-2dde-47df-beb2-5bee30f6d5c0/scratchpad
REF="$WT/exoplasim/bench/ref"
RES=T21; LOW=t21; NAME="most_plasim_t21_l10_p16.x"

cp -f "$SRC" "$SP/legmod_folded.f90"
git -C "$WT" show HEAD:vendor/exoplasim/exoplasim/plasim/src/legmod.f90 > "$SP/legmod_unfolded.f90"
cp -f "$PKG/most_compiler_mpi" "$PKG/most_compiler_mpi.orig"
restore() {
    cp -f "$SP/legmod_folded.f90" "$SRC"
    cp -f "$PKG/most_compiler_mpi.orig" "$PKG/most_compiler_mpi"
}
trap restore EXIT

# The wrong arm: same fold, off by one in the mode index.
sed 's/skspgp(n+1)/skspgp(n)/g; s/skgpsp(n+1)/skgpsp(n)/g' \
    "$SP/legmod_folded.f90" > "$SP/legmod_wrong.f90"
if cmp -s "$SP/legmod_folded.f90" "$SP/legmod_wrong.f90"; then
    echo "the wrong arm is identical to the folded arm; the sed did not bite" >&2
    exit 1
fi

# Inject the power-of-two test filter into every arm, as nfilter == 9.
inject () {
python - "$1" <<'PY'
import sys, pathlib
p = pathlib.Path(sys.argv[1]); s = p.read_text()
anchor = "   else\n     skgpsp(n) = 1.0\n     skspgp(n) = 1.0\n   endif"
assert anchor in s, "filter else-branch not found"
s = s.replace(anchor,
"""   else if (nfilter .eq. 9) then !TEST ONLY: exact powers of two, see
     !verify_fold_indexing.sh. Alternating so an off-by-one in the mode index
     !changes the factor on EVERY mode rather than only at a boundary.
     skgpsp(n) = 2.0**(-mod(n,2))
     skspgp(n) = 2.0**(-mod(n,2))

   else
     skgpsp(n) = 1.0
     skspgp(n) = 1.0
   endif""")
p.write_text(s)
PY
}

# Contraction off: with it on, deleting the multiply changes which FMA the
# compiler forms and the arms differ for reasons unrelated to indexing.
python - "$PKG/most_compiler_mpi" <<'PY'
import sys, pathlib
p = pathlib.Path(sys.argv[1])
p.write_text("\n".join(
    l + " -ffp-contract=off" if l.startswith("MOST_F90_OPTS=") else l
    for l in p.read_text().splitlines()) + "\n")
PY

for arm in unfolded folded wrong; do
    cp -f "$SP/legmod_${arm}.f90" "$SRC"
    inject "$SRC"
    ( cd "$PKG" && ./compile.sh -n 16 -p 8 -r $RES -v 10 -O march=znver4 >/dev/null 2>&1 || true )
    [ -f "$PKG/plasim/run/$NAME" ] || { echo "build failed: $arm" >&2; exit 1; }
    cp -f "$PKG/plasim/run/$NAME" "$REF/idx_${arm}_${LOW}.x"
    echo "built $arm"
done

work="$WT/exoplasim/bench/_idxcheck"
rm -rf "$work"; mkdir -p "$work"
cp -a "$WT/exoplasim/bench/bed_t21/." "$work/"
cd "$work"
rm -f MOST_REST.* plasim_status Abort_Message
sed -i 's/^ *NFILTER *=.*/ NFILTER = 9 /; s/^ *NGPTFILTER *=.*/ NGPTFILTER = 1 /; s/^ *NSPVFILTER *=.*/ NSPVFILTER = 1 /; s/^ *N_RUN_STEPS *=.*/ N_RUN_STEPS = 20 /' plasim_namelist
grep -E "NFILTER|NGPTFILTER|NSPVFILTER|N_RUN_STEPS" plasim_namelist
echo

declare -A SHA
for arm in unfolded folded wrong; do
    rm -f plasim_status Abort_Message
    cp -f "$REF/idx_${arm}_${LOW}.x" ./probe.x
    mpiexec -np 16 ./probe.x >run.log 2>&1 || true
    if [ ! -f plasim_status ]; then
        # The `wrong` arm failing is a PASS: an off-by-one makes legini read
        # skspgp(0), which -fcheck=all catches as an out-of-bounds and aborts.
        # It never reaches an integration to diverge in. Any other arm failing
        # is a real failure.
        if [ "$arm" = "wrong" ]; then
            SHA[$arm]="(aborted in legini: $(grep -c 'legini_' run.log) ranks)"
            echo "$arm  ABORTED -- out-of-bounds skspgp(0), caught by -fcheck=all"
            rm -f probe.x
            continue
        fi
        echo "$arm: no plasim_status (model failed under the test filter)" >&2
        grep -v "Fortran runtime warning\|^At line" run.log | tail -4 >&2
        exit 1
    fi
    SHA[$arm]=$(sha256sum plasim_status | cut -c1-32)
    echo "$arm  ${SHA[$arm]}"
    rm -f probe.x
done

echo
ok=0
if [ "${SHA[unfolded]}" = "${SHA[folded]}" ]; then
    echo "PASS  folded == unfolded, bit for bit, with a filter that VARIES by mode."
    echo "      Each mode receives its own filter value; the n+1 offset is right."
else
    echo "FAIL  folded != unfolded. The fold applies the wrong filter per mode."; ok=1
fi
if [ "${SHA[wrong]}" != "${SHA[folded]}" ]; then
    echo "PASS  the deliberately off-by-one arm DIVERGES, so the test can fail."
else
    echo "FAIL  off-by-one is indistinguishable; the test proves nothing."; ok=1
fi
exit $ok
