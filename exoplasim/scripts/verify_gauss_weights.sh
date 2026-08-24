#!/bin/bash
# Are the model's Gaussian nodes and weights right?
#
#   exoplasim/scripts/verify_gauss_weights.sh [res ...]
#
# Worldbuilding frame: a correctness check on the Vesper climate model's
# spectral transform. Nothing here is about the simulated planet.
#
# WHY IT EXISTS. The weights enter the FORWARD transform and nothing else --
# fc2sp carries gwd, uv2dv and mktend carry gwd/cos^2, and the synthesis
# direction never touches them. So an error here is invisible to every check
# that compares synthesis, which is how the previous inigau carried weights
# wrong by 1e-10 at the pole-most latitude for as long as this model has run.
# It surfaced only as an unexplained floor on the ANALYSIS arms of
# verify_shtns_equivalence.sh. CLIM-73.
#
# THE REFERENCE IS INDEPENDENT AND HIGHER PRECISION, which is what makes this a
# test rather than a comparison: gauss_weight_reference.py runs the same
# textbook algorithm in x86 longdouble, eps about 1.1e-19, eight orders of
# margin on the disagreement being judged. Three double-precision
# implementations disagreed here -- inigau, SHTns and numpy's leggauss -- and
# double precision cannot settle which is right.
#
# THE BARS ARE DERIVED FROM THE ALGORITHM, not chosen to be met. Evaluating P_n
# by the three-term recurrence accumulates about n*eps, and the weight squares
# the derivative, so the floor is about 2*n*eps -- 1.1e-13 at NLAT 256. The bar
# is 1e-12, an order above that. Beating it needs a different algorithm rather
# than a rearrangement of this one; converging the angle instead of the node,
# which avoids forming 1-z^2 near the pole, was tried and changed nothing
# outside noise.
#
# NODES ARE JUDGED IN ABSOLUTE TERMS and weights in relative. A relative bar on
# the nodes is not a well-posed criterion: they are cosines and pass through
# zero at the equator, so an absolute error of 1e-17 there reads as a relative
# error of 3e-15 and the number says more about which latitude happens to sit
# nearest the equator than about the algorithm. They live in [-1,1], so absolute
# error is the meaningful measure and 1e-15 is the bar.
#
# The sum of the weights must be 2 to 1e-15. That one is an identity, and it is
# the only check here that does not depend on the reference being right.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
SRC="$REPO/vendor/exoplasim/exoplasim/plasim/src"
WORK="$REPO/exoplasim/bench/_gaussweights"
TOL=1e-12       # weights, relative
NODETOL=1e-15   # nodes, ABSOLUTE: they pass through zero
SUMTOL=1e-15
# The commit before the rewrite: the control is the code that actually shipped.
CONTROL_REV="${CONTROL_REV:-447fdac9}"

resolutions=("$@")
# The default sweep is the ladder, taken from `lib/rungs.py` rather than
# written out: the five counts that stood here were the ladder with T31,
# T63 and T106 missing. SPAT-2.
if [ ${#resolutions[@]} -eq 0 ]; then
    PY="$REPO/.venv/bin/python"; [ -x "$PY" ] || PY=python3
    # shellcheck disable=SC2207
    resolutions=($("$PY" -c "import sys; sys.path.insert(0, '$REPO/lib'); import rungs; print(' '.join(str(rungs.geometry(r)[0]) for r in rungs.RUNGS))"))
fi

rm -rf "$WORK"; mkdir -p "$WORK"
cd "$WORK"
cp "$SRC/gaussmod.f90" .

cat > drive.f90 <<'EOF'
      program gaussdump
      implicit none
      integer :: klat, j
      real (kind=8), allocatable :: z0(:), zw(:)
      character(len=32) :: yarg
      call get_command_argument(1, yarg)
      read(yarg,*) klat
      allocate(z0(klat), zw(klat))
      call inigau(klat, z0, zw)
      do j = 1 , klat
         write(*,'(i6,a,e26.17,a,e26.17)') j, ',', z0(j), ',', zw(j)
      enddo
      end program gaussdump
EOF

gfortran -O2 -o dump.x drive.f90 gaussmod.f90 2>build.err \
    || { echo "build failed"; head -20 build.err; exit 1; }

echo "the model's inigau against extended precision"
echo "declared before the run: weights to $TOL relative, nodes to $NODETOL absolute,"
echo "  sum of weights to $SUMTOL of 2"
echo
rc=0
for n in "${resolutions[@]}"; do
    ./dump.x "$n" > "model_$n.csv"
    "$REPO"/.venv/bin/python - "$n" "$WORK/model_$n.csv" "$TOL" "$SUMTOL" "$NODETOL" <<'PY' || rc=1
import sys, importlib.util
import numpy as np

n, path = int(sys.argv[1]), sys.argv[2]
tol, sumtol, nodetol = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
here = __import__("pathlib").Path(path).parents[2] / "scripts" / "gauss_weight_reference.py"
spec = importlib.util.spec_from_file_location("gref", here)
gref = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gref)

nodes, wts = gref.gauss_legendre(n)
mx, mw = [], []
for line in open(path):
    _, x, w = line.split(",")
    mx.append(float(x)); mw.append(float(w))
mx = np.array(mx); mw = np.array(mw)
if mx.size != n:
    print(f"  [ FAIL ] NLAT {n}: dumped {mx.size} latitudes"); raise SystemExit(1)

dn = np.max(np.abs(mx - nodes.astype(float)))
dw = np.max(np.abs((mw - wts.astype(float)) / wts.astype(float)))
ds = abs(float(np.sum(mw)) - 2.0)
bad = (dn > nodetol) or (dw > tol) or (ds > sumtol)
tag = "[ FAIL ]" if bad else "[  ok  ]"
print(f"  {tag} NLAT {n:>4}   nodes {dn:.3e}   weights {dw:.3e}   sum-2 {ds:.3e}")
raise SystemExit(1 if bad else 0)
PY
done

# THE CONTROL IS THE IMPLEMENTATION THIS REPLACED, which is the sharpest one
# available: not a synthetic corruption but the code that actually shipped, and
# it must fail. If it passes, the check cannot see the defect it exists for.
echo
echo "==== the control, the trigonometric-series inigau: must FAIL ===="
if git -C "$REPO" show "$CONTROL_REV:vendor/exoplasim/exoplasim/plasim/src/gaussmod.f90" > old_gaussmod.f90 2>/dev/null && [ -s old_gaussmod.f90 ]; then
    if gfortran -O2 -o dump_old.x drive.f90 old_gaussmod.f90 2>build_old.err; then
        ./dump_old.x 192 > model_old_192.csv
        worst=$("$REPO"/.venv/bin/python "$HERE/gauss_weight_reference.py" 192 --compare "$WORK/model_old_192.csv")
        over=$(awk -v a="$worst" -v b="$TOL" 'BEGIN{print (a>b)?1:0}')
        if [ "$over" = 1 ]; then
            echo "  [  ok  ] control rejected at $worst, so the quadrature is under test"
        else
            echo "  [ FAIL ] the shipped implementation passed at $worst"
            rc=1
        fi
    else
        echo "  [ FAIL ] the control did not build"; head -10 build_old.err; rc=1
    fi
else
    echo "  [ FAIL ] could not retrieve the control from $CONTROL_REV"
    rc=1
fi

echo
if [ "$rc" = 0 ]; then
    echo "PASS: the quadrature is right to $TOL."
else
    echo "FAIL: see above."
fi
exit "$rc"
