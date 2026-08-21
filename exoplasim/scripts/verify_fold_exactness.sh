#!/bin/bash
# Exactness test with FMA contraction DISABLED.
#
#   exoplasim/scripts/verify_fold_exactness.sh T21
#
# Why this exists. The first exactness test compared patched against unpatched
# with the filters off, expecting bit-identical results because the fold then
# multiplies by exactly 1.0. It failed everywhere. The premise was wrong:
# multiplying by 1.0 IS exact, but DELETING the multiply changes what the
# compiler can contract. gfortran defaults to -ffp-contract=fast and znver4 has
# FMA, so
#     sp + qc*fc*sk   ->  t = round(qc*fc); t2 = t*sk; round(sp + t2)
#     sp + qc*fc      ->  fma(qc, fc, sp)     -- ONE rounding, not two
# and the two differ in the last bits for reasons that have nothing to do with
# whether the fold indexes the right spectral mode.
#
# With -ffp-contract=off neither form can contract, so the filters-off
# comparison becomes a real test of the ALGEBRA again: pass means the fold picks
# the same mode the transform does, fail means the n+1 offset is wrong.
set -euo pipefail
res="${1:?usage: verify_nofma.sh <T21|T42|T85>}"
low=$(echo "$res" | tr 'A-Z' 'a-z')
WT=/home/cfutro/docs/world/.claude/worktrees/spectral-optimisation
PKG="$WT/vendor/exoplasim/exoplasim"
SRC="$PKG/plasim/src/legmod.f90"
KEEP="$WT/exoplasim/bench/legmod_patched.f90"
REF="$WT/exoplasim/bench/ref"
name="most_plasim_${low}_l10_p16.x"

cp -f "$SRC" "$KEEP"                      # the patched source, to restore
# cp, not $(cat) + printf: command substitution strips trailing newlines and
# printf '%s' does not put one back, so restoring that way left
# most_compiler_mpi without its final newline. compile.sh cats these files
# together, so the next build produced a makefile target literally named
# `pumax_stubMOST_PREC=-fdefault-real-8.c` and failed with no hint of why.
cp -f "$PKG/most_compiler_mpi" "$PKG/most_compiler_mpi.orig"
restore() {
    cp -f "$KEEP" "$SRC"
    cp -f "$PKG/most_compiler_mpi.orig" "$PKG/most_compiler_mpi"
}
trap restore EXIT

# -ffp-contract=off for both arms
python - "$PKG/most_compiler_mpi" <<'PY'
import sys, pathlib
p = pathlib.Path(sys.argv[1]); out = []
for line in p.read_text().splitlines():
    if line.startswith("MOST_F90_OPTS=") and "-ffp-contract" not in line:
        line += " -ffp-contract=off"
    out.append(line)
p.write_text("\n".join(out) + "\n")
PY
grep MOST_F90_OPTS "$PKG/most_compiler_mpi"

cd "$PKG"
echo "building UNPATCHED $res with contraction off ..."
git -C "$WT" checkout -- vendor/exoplasim/exoplasim/plasim/src/legmod.f90
./compile.sh -n 16 -p 8 -r "$res" -v 10 -O march=znver4 >/dev/null 2>&1 || true
[ -f "plasim/run/$name" ] || { echo "unpatched build failed" >&2; exit 1; }
cp -f "plasim/run/$name" "$REF/unpatched_nofma_${low}.x"

echo "building PATCHED   $res with contraction off ..."
cp -f "$KEEP" "$SRC"
./compile.sh -n 16 -p 8 -r "$res" -v 10 -O march=znver4 >/dev/null 2>&1 || true
[ -f "plasim/run/$name" ] || { echo "patched build failed" >&2; exit 1; }
cp -f "plasim/run/$name" "$REF/patched_nofma_${low}.x"

cd "$WT"
bed="exoplasim/bench/bed_$low"
exoplasim/scripts/verify_filter_fold.sh "$bed" \
    "$REF/unpatched_nofma_${low}.x" "$REF/patched_nofma_${low}.x" 16 20
