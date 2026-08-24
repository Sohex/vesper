#!/usr/bin/env python3
"""The independent second side for the banded spectral analysis, as a file.

    python exoplasim/scripts/banded_transform_reference.py T21 --case dense --out ref.bin

Worldbuilding frame: a numerical reference for the Vesper climate model's
spectral transform. Nothing here is about the simulated planet.

WHAT THIS IS FOR. `verify_banded_transform.f90` drives the model's OWN analysis
chain -- the scatter that gives a thread its latitudes, `legini`, `fc2sp` over
that thread's band, and the reduction that adds the bands up -- and needs
something to compare the answer against. Everything it links is the subject of
the check, so the reference cannot come from any of it. It comes from here, as
a table of numbers the driver reads and never recomputes.

WHAT IS INDEPENDENT AND WHAT IS NOT, stated rather than implied.

INDEPENDENT: the associated Legendre functions. `legini` evaluates them by its
own three-term recurrence in mu; this takes them from `scipy.special.
sph_legendre_p_all`, a different implementation of a different formulation,
which shares no line of code and no recurrence with the model. The mapping
between the two is one constant and one factor, and both are DERIVED rather
than fitted:

    pmat(w,l) = sqrt(2*pi) * Pbar_n^m(theta_l)
    qmat(w,l) = (1 - mu^2) * d/dmu [ sqrt(2*pi) * Pbar_n^m ]
              = -sqrt(1 - mu^2) * d/dtheta [ sqrt(2*pi) * Pbar_n^m ]

sqrt(2*pi) is the ratio between the model's normalisation, which puts
P(0,0) = sqrt(1/2), and the fully normalised spherical harmonic, which puts it
at 1/sqrt(4*pi). The (1 - mu^2) is the factor the model carries because its
wind variables are u*cos(phi) rather than u.

NOT INDEPENDENT, and this is the judgement call: the Gaussian nodes and
weights. `gauss_weight_reference.gauss_legendre` is the SAME textbook algorithm
`inigau` runs -- Newton on P_n by the three-term recurrence, then
w = 2/((1-x^2) P'^2) -- carried out in x86 longdouble instead of float64. That
adjudicates ROUNDING, by eight orders, and it does not adjudicate METHOD: if
the textbook construction were itself the wrong thing to compute, both sides
would be wrong together and this file would not notice.

Reproducing the quadrature by a genuinely different route -- Golub-Welsch on
the Jacobi matrix -- would close that, and needs an extended-precision
eigensolver this project does not carry. What closes it instead are two
identities that involve no method at all and are checked below before anything
is written: the weights sum to exactly 2, and the rule integrates every
Legendre polynomial of degree 1 to 2*NLAT-1 to zero. A construction that passes
both is a Gauss-Legendre rule on NLAT points, whatever produced it, because
that rule is unique. So the method hole is closed by an identity even though
the implementations are related.

WHAT THE CHOICE MAKES INVISIBLE, in one line: nothing, given those two
identities hold; without them it would hide a shared misreading of what
Gauss-Legendre quadrature is.

THE RIGHT ANSWER, fixed before any arm runs. The driver is handed Fourier
coefficients

    fc(c,m,l) = sum over n of pmat_ref(w(m,n), l) * spdrv(c, w(m,n))

and the model's analysis of them must return `spdrv`, because

    sum over l of gwd(l) * pmat(w,l) * pmat(w',l) = delta(w,w')

for two modes of the same m. That is orthonormality under the quadrature, it is
a property of the functions and not of either implementation, and it is what
makes `spdrv` a right answer rather than a second opinion. This file measures
its own residual on that identity and prints it, so the reference's own error
can be compared against the bound the gate runs at rather than assumed smaller.

WHY THE FOURIER HALF IS ABSENT. The band decomposition splits LATITUDES. The
Fourier transform runs along a latitude row and is local to whichever thread
owns the row, so it carries no band structure, no partial sum and no reduction,
and putting it in would only add a convention -- the sign of the exponent --
that has to be agreed between the two sides for no gain. The driver starts from
Fourier coefficients for exactly that reason.

THE FILE. Little-endian float64, no record markers, in this order:

    header  8 values: nlat, nlon, ntru, ncsp, npro, nlpp, scale, version
    sid     nlat            sine of latitude, descending from +1
    gwd     nlat            Gaussian weights, same order
    pmat    ncsp*nlat       w fastest, then l
    qmat    ncsp*nlat       w fastest, then l
    spdrv   2*ncsp          component fastest, then mode
    fc      nlon*nlat       the model's fc(2,NLON/2,NLAT) layout exactly

`scale` is what the driving vector was divided by so that the largest Fourier
coefficient is 1. Fixing that magnitude is what lets the driver's tolerance be
an absolute number derived from the length of the sum.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import sph_legendre_p_all

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gauss_weight_reference import gauss_legendre  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import rungs  # noqa: E402

LD = np.longdouble
VERSION = 1.0

# The identities the quadrature is pinned by. Both are exact statements, so the
# bars are rounding bars and nothing else, and both are checked at the rung the
# reference is being built for rather than argued in general.
SUM_BAR = 1e-17        # |sum(w) - 2|, in longdouble; eps there is 1.1e-19
MOMENT_BAR = 1e-16     # |sum(w P_k(x))| for 1 <= k <= 2*NLAT-1


def triangular_modes(ntru: int) -> list[tuple[int, int]]:
    """The (m,n) modes in the model's packed order: m outer, n from m to NTRU.

    `legini` fills `pmat(lm,jlat)` walking this order, so the mode index here
    and the model's `w` are the same integer. Anything else here would show up
    as a gross mismatch rather than a subtle one.
    """
    return [(m, n) for m in range(ntru + 1) for n in range(m, ntru + 1)]


def legendre_moments(nodes: LD, weights: LD, degree: int) -> float:
    """The worst |sum(w P_k(x))| over 1 <= k <= degree, by the recurrence.

    The exact value is zero for every k >= 1, because P_0 = 1 is orthogonal to
    all of them. Evaluated in longdouble, so a residual above MOMENT_BAR is the
    rule being wrong and not the arithmetic.
    """
    p0 = np.ones_like(nodes)
    p1 = nodes.copy()
    worst = float(abs(np.sum(weights * p1)))
    for k in range(2, degree + 1):
        p0, p1 = p1, ((2 * k - 1) * nodes * p1 - (k - 1) * p0) / k
        worst = max(worst, float(abs(np.sum(weights * p1))))
    return worst


def legendre_table(ntru: int, nodes: LD):
    """pmat and qmat for every mode at every node, from scipy.

    Returned as (ncsp, nlat) arrays in longdouble. scipy computes in float64,
    so the table itself is good to about eps and no better; that is three
    hundred times inside the bound the driver runs at, and the orthonormality
    residual printed by `main` is the measurement of it.
    """
    modes = triangular_modes(ntru)
    ms = np.array([m for m, _n in modes])
    ns = np.array([n for _m, n in modes])
    nlat = nodes.size
    pmat = np.zeros((len(modes), nlat), dtype=LD)
    qmat = np.zeros((len(modes), nlat), dtype=LD)
    s = np.sqrt(LD(2) * LD(np.pi))
    for l in range(nlat):
        mu = float(nodes[l])
        theta = math.acos(mu)
        sin_theta = math.sqrt(1.0 - mu * mu)
        pt, dt = sph_legendre_p_all(ntru, ntru, theta, diff_n=1)
        pmat[:, l] = np.asarray(pt)[ns, ms].astype(LD) * s
        # d/dmu = -(1/sin theta) d/dtheta, and the model carries (1-mu^2)
        # times that.
        qmat[:, l] = -LD(sin_theta) * np.asarray(dt)[ns, ms].astype(LD) * s
    return pmat, qmat


def orthonormality_residual(pmat: LD, weights: LD, modes) -> float:
    """The worst |sum_l w_l P_w P_w' - delta| over modes sharing an m.

    This is the identity the driver's right answer rests on, measured on the
    reference alone. It is the reference's own error bar, and it is printed so
    that "the bound is larger than the instrument" is a comparison rather than
    an assumption.
    """
    worst = 0.0
    by_m: dict[int, list[int]] = {}
    for k, (m, _n) in enumerate(modes):
        by_m.setdefault(m, []).append(k)
    # float64 for the Gram matrices: the residual being measured is the
    # reference table's own, which is a float64 table, so carrying the product
    # in longdouble would measure nothing extra and costs an unvectorised
    # matmul that runs to minutes at the top of the ladder.
    w64 = weights.astype(np.float64)
    p64 = pmat.astype(np.float64)
    for _m, ks in by_m.items():
        block = p64[ks, :]
        gram = (block * w64) @ block.T
        eye = np.eye(len(ks), dtype=np.float64)
        worst = max(worst, float(np.max(np.abs(gram - eye))))
    return worst


def driving_vector(case: str, ntru: int, modes) -> LD:
    """The spectral vector the analysis has to give back, as (2, ncsp).

    The imaginary component of every m = 0 mode is zero. The zonal mean has no
    imaginary part, the model discards it, and driving one would make the right
    answer depend on how that discard is spelled rather than on the transform.
    """
    ncsp = len(modes)
    sp = np.zeros((2, ncsp), dtype=LD)
    if case == "dense":
        # A fixed seed, so the vector is a property of the case and not of the
        # run. Every mode is driven: the bands are what is under test and a
        # dense spectrum is what crosses all of them.
        rng = np.random.default_rng(20260824)
        draw = rng.standard_normal((2, ncsp))
        sp[:, :] = draw.astype(LD)
    elif case == "corner":
        # One mode, the (NTRU,NTRU) corner: the last the recurrence reaches and
        # the one whose weights are smallest, so a transform that is right on a
        # dense spectrum and wrong on one coefficient shows up here. The
        # planetary vorticity mode was that shape of defect.
        sp[0, ncsp - 1] = LD(1)
        sp[1, ncsp - 1] = LD(1)
    elif case == "zonal":
        # Every m = 0 mode and nothing else. The zonal mean is the one the
        # model's own global sums are taken through, and it is the column of
        # the transform a band error cannot cancel within a latitude.
        for k, (m, _n) in enumerate(modes):
            if m == 0:
                sp[0, k] = LD(1)
    else:
        raise SystemExit(f"unknown case {case!r}: dense, corner or zonal")
    for k, (m, _n) in enumerate(modes):
        if m == 0:
            sp[1, k] = LD(0)
    return sp


def synthesise_fc(sp: LD, pmat: LD, modes, nlon: int, nlat: int) -> LD:
    """The Fourier coefficients the driver is handed, in the model's layout.

    fc(c,m,l) = sum over the modes of this m of pmat(w,l) * sp(c,w). This is
    the synthesis half done on the REFERENCE table, which is the whole point:
    if the model's own synthesis were used instead, an error in `pmat` would
    cancel between the two halves and the round trip would hold on a wrong
    table. That is the failure this arm exists to be immune to.
    """
    fc = np.zeros((2, nlon // 2, nlat), dtype=LD)
    for k, (m, _n) in enumerate(modes):
        fc[0, m, :] += pmat[k, :] * sp[0, k]
        fc[1, m, :] += pmat[k, :] * sp[1, k]
    return fc


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("res", help="a rung in lib/rungs.py, e.g. T21")
    ap.add_argument("--npro", type=int, default=1,
                    help="the thread count the driver will run at; written to "
                         "the header so the driver refuses a mismatched file")
    ap.add_argument("--case", default="dense",
                    choices=("dense", "corner", "zonal"))
    ap.add_argument("--out", required=True, help="the binary the driver reads")
    args = ap.parse_args()

    nlat, nlon, ntru = rungs.geometry(args.res)
    if nlat % args.npro:
        raise SystemExit(f"--npro {args.npro} does not divide the {nlat} "
                         f"latitudes of {args.res}")
    modes = triangular_modes(ntru)
    ncsp = len(modes)

    nodes, weights = gauss_legendre(nlat)

    # The two identities that pin the quadrature, before anything else is
    # computed from it. Either one failing means the reference is not a
    # Gauss-Legendre rule and nothing below it means anything.
    dsum = float(abs(np.sum(weights) - LD(2)))
    dmom = legendre_moments(nodes, weights, 2 * nlat - 1)
    if dsum > SUM_BAR or dmom > MOMENT_BAR:
        print(f"refusing: the reference quadrature is not a Gauss-Legendre rule "
              f"on {nlat} nodes.", file=sys.stderr)
        print(f"  |sum(w) - 2| = {dsum:.3e}, bar {SUM_BAR:.0e}", file=sys.stderr)
        print(f"  worst |sum(w P_k)| for 1 <= k <= {2 * nlat - 1} = {dmom:.3e}, "
              f"bar {MOMENT_BAR:.0e}", file=sys.stderr)
        return 2

    pmat, qmat = legendre_table(ntru, nodes)
    resid = orthonormality_residual(pmat, weights, modes)

    sp = driving_vector(args.case, ntru, modes)
    fc = synthesise_fc(sp, pmat, modes, nlon, nlat)
    scale = float(np.max(np.abs(fc)))
    if scale <= 0.0:
        raise SystemExit("the driving vector synthesised to nothing")
    sp = sp / LD(scale)
    fc = fc / LD(scale)

    header = np.array([nlat, nlon, ntru, ncsp, args.npro, nlat // args.npro,
                       scale, VERSION], dtype=np.float64)
    with open(args.out, "wb") as fh:
        header.tofile(fh)
        nodes.astype(np.float64).tofile(fh)
        weights.astype(np.float64).tofile(fh)
        # Fortran order: w fastest, matching pmat(NCSP,NLPP).
        pmat.astype(np.float64).T.reshape(-1).tofile(fh)
        qmat.astype(np.float64).T.reshape(-1).tofile(fh)
        # component fastest, matching sp(2,NCSP)
        sp.astype(np.float64).T.reshape(-1).tofile(fh)
        # component fastest, then m, then l: fc(2,NLON/2,NLAT)
        fc.astype(np.float64).transpose(2, 1, 0).reshape(-1).tofile(fh)

    print(f"{args.res}  NLAT {nlat}  NLON {nlon}  NTRU {ntru}  modes {ncsp}  "
          f"NPRO {args.npro}  case {args.case}")
    print(f"  quadrature identities: |sum(w) - 2| = {dsum:.3e}, "
          f"worst Legendre moment {dmom:.3e}")
    print(f"  orthonormality residual of the reference table: {resid:.3e}")
    print(f"  the driving vector was scaled by {scale:.6e}, so the largest "
          f"Fourier coefficient is 1")
    print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
