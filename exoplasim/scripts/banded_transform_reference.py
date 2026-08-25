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

THE SYNTHESIS HALF, AND WHY IT IS HERE TOO. `dv2uv` turns spectral divergence
and ABSOLUTE vorticity into the two wind components on the latitudes a thread
holds. Its only band dependence is which `pmat` and `qmat` rows the thread has
and which rows of the shared wind globe it writes, and both of those are
questions only a banded driver can ask -- at one thread the local latitude
index and the global one are the same integer, so a term applied at the wrong
latitude is applied at the right one. `verify_inverse_transform.py` drives the
same loop at one thread against the same arithmetic and cannot see that class
at all.

PLANETARY VORTICITY, and why the reference needs no formula for it. `dv2uv`
takes ABSOLUTE vorticity and removes the planetary part from its RESULT rather
than from the input, because the input aliases the shared spectral state. The
two are the same thing by an identity that owes nothing to how the removal is
coded:

    dv2uv(pd, pz) with plavor = P  ==  dv2uv(pd, pz - P at mode w=2) with P = 0

so `pu` and `pv` below are computed from the RELATIVE vorticity with no
planetary term at all, and the vorticity handed to the driver is that plus P at
the first component of mode w=2. The mode, the component, the sign and the
factor are then all pinned by where P enters, and none of them is read off the
model. P is 1.7, the same declared value `verify_inverse_transform.py` uses, so
one number means one thing in both gates.

NOT INDEPENDENT, the second one: the per-mode factors `fmu` and `fmv`. The
driver refuses a filtered build, so `legini` leaves the physics filter at one
and the factors are m/(n(n+1)) and 1/(n(n+1)); those are restated here rather
than derived from anything else, exactly as `verify_inverse_transform.py`
restates them. A shared misreading of what the factor IS would be invisible to
both. What this arm is for is where the factor lands, not what it is.

`fmu(2)` is zero, because mode w=2 is m = 0 and `fmu` carries a factor of m.
So the model's second planetary-vorticity line, the one that writes `pv`, adds
exactly nothing whatever P is, and no arm anywhere can give it teeth. It is
stated here rather than left for a reader to rediscover from a control that
passes.

THE FILE. Little-endian float64, no record markers, in this order:

    header  12 values: nlat, nlon, ntru, ncsp, npro, nlpp, scale, version,
                       nlev, plavor, tscale, dvfmax
    sid     nlat            sine of latitude, descending from +1
    gwd     nlat            Gaussian weights, same order
    pmat    ncsp*nlat       w fastest, then l
    qmat    ncsp*nlat       w fastest, then l
    spdrv   2*ncsp          component fastest, then mode
    fc      nlon*nlat       the model's fc(2,NLON/2,NLAT) layout exactly
    fmu     ncsp
    fmv     ncsp
    spd     2*ncsp*nlev     component fastest, then mode, then level
    spz     2*ncsp*nlev     ABSOLUTE vorticity, same layout
    pu      2*(nlon/2)*nlat*nlev   the model's pu(2,NLON/2,NLAT,NLEV) layout
    pv      2*(nlon/2)*nlat*nlev   the same

`scale` is what the analysis driving vector was divided by so that the largest
Fourier coefficient is 1. Fixing that magnitude is what lets the driver's
tolerance be an absolute number derived from the length of the sum.

`tscale` and `dvfmax` do the same job for the synthesis arm and are measured
rather than declared. `tscale` is the largest sum of TERM MAGNITUDES that
`dv2uv` accumulates into any one output element, the planetary correction
included, which is what bounds the rounding of that accumulation. `dvfmax` is
the largest |factor| times the largest |coefficient| over the driving fields,
which is what the weight-matrix error of one mode gets multiplied by. The
synthesis fields are NOT normalised: the bound is stated in terms of these two
numbers and both travel in the file.
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
VERSION = 2.0

# The planetary vorticity the synthesis arm drives, the same declared value
# verify_inverse_transform.py uses. It is a number this reference chooses, not
# one read off the planet: what is under test is where the term lands, and a
# value shared between the two gates keeps one number meaning one thing.
PLAVOR = 1.7

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


def per_mode_factors(modes) -> tuple[LD, LD]:
    """fmu and fmv, at the unfiltered build the driver refuses to run without.

    `legini` sets fmu(w) = m/(n(n+1)) * skspgp(n+1) and fmv(w) = 1/(n(n+1)) *
    skspgp(n+1), and skspgp is one at nfilter = 0. These are RESTATED, not
    derived from anything independent; the reference header says so and says
    what that leaves invisible.
    """
    ncsp = len(modes)
    fmu = np.zeros(ncsp, dtype=LD)
    fmv = np.zeros(ncsp, dtype=LD)
    for k, (m, n) in enumerate(modes):
        if n > 0:
            znn1 = LD(1) / LD(n * (n + 1))
            fmu[k] = znn1 * LD(m)
            fmv[k] = znn1
    return fmu, fmv


def driving_field(case: str, modes, nlev: int) -> LD:
    """A (2, ncsp, nlev) spectral field shaped by the same three cases.

    The levels differ from one another rather than repeating, so a level index
    pinned to one inside the routine is a wrong answer here and not a passing
    one. The seed is fixed, so the field is a property of the case.
    """
    ncsp = len(modes)
    out = np.zeros((2, ncsp, nlev), dtype=LD)
    if case == "dense":
        rng = np.random.default_rng(20260825)
        out[:, :, :] = rng.standard_normal((2, ncsp, nlev)).astype(LD)
    elif case == "corner":
        for v in range(nlev):
            out[0, ncsp - 1, v] = LD(1 + v)
            out[1, ncsp - 1, v] = LD(1 + v)
    elif case == "zonal":
        for k, (m, _n) in enumerate(modes):
            if m == 0:
                for v in range(nlev):
                    out[0, k, v] = LD(1 + v)
    else:
        raise SystemExit(f"unknown case {case!r}: dense, corner or zonal")
    for k, (m, _n) in enumerate(modes):
        if m == 0:
            out[1, k, :] = LD(0)
    return out


def synthesise_uv(spd: LD, spzrel: LD, pmat: LD, qmat: LD, fmu: LD, fmv: LD,
                  modes, nlon: int, nlat: int, nlev: int):
    """dv2uv's four outputs, and the term-magnitude bound the driver needs.

    The four sign patterns are the transform's definition and are written out
    once here:

        pu1 = sum  qmat*fmv*pz1 + pmat*fmu*pd2
        pu2 = sum  qmat*fmv*pz2 - pmat*fmu*pd1
        pv1 = sum  pmat*fmu*pz2 - qmat*fmv*pd1
        pv2 = sum -pmat*fmu*pz1 - qmat*fmv*pd2

    with pz the RELATIVE vorticity, which is what makes this a reference for
    the planetary term rather than a restatement of the model's formula for it.

    `tscale` is taken over the ABSOLUTE vorticity, the field the model actually
    accumulates, and carries the planetary correction as a term of its own,
    because that accumulation is what rounds.
    """
    nm = nlon // 2
    pu = np.zeros((2, nm, nlat, nlev), dtype=LD)
    pv = np.zeros((2, nm, nlat, nlev), dtype=LD)
    au = np.zeros((2, nm, nlat, nlev), dtype=LD)
    av = np.zeros((2, nm, nlat, nlev), dtype=LD)
    spz = spzrel.copy()
    spz[0, 1, :] += LD(PLAVOR)          # mode w=2, first component
    for k, (m, _n) in enumerate(modes):
        qv = qmat[k, :, None] * fmv[k]          # (nlat, 1)
        pu_ = pmat[k, :, None] * fmu[k]
        zz1, zz2 = spzrel[0, k, :], spzrel[1, k, :]
        zd1, zd2 = spd[0, k, :], spd[1, k, :]
        pu[0, m] += qv * zz1 + pu_ * zd2
        pu[1, m] += qv * zz2 - pu_ * zd1
        pv[0, m] += pu_ * zz2 - qv * zd1
        pv[1, m] += -pu_ * zz1 - qv * zd2
        # The magnitudes the model's own accumulation carries: absolute
        # vorticity, not relative.
        az1, az2 = np.abs(spz[0, k, :]), np.abs(spz[1, k, :])
        ad1, ad2 = np.abs(spd[0, k, :]), np.abs(spd[1, k, :])
        aq, ap = np.abs(qv), np.abs(pu_)
        au[0, m] += aq * az1 + ap * ad2
        au[1, m] += aq * az2 + ap * ad1
        av[0, m] += ap * az2 + aq * ad1
        av[1, m] += ap * az1 + aq * ad2
    # The planetary correction is a term of the accumulation too.
    au[0, 0] += np.abs(qmat[1, :, None] * fmv[1] * LD(PLAVOR))
    av[1, 0] += np.abs(pmat[1, :, None] * fmu[1] * LD(PLAVOR))
    tscale = float(max(np.max(au), np.max(av)))
    return pu, pv, spz, tscale


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("res", help="a rung in lib/rungs.py, e.g. T21")
    ap.add_argument("--npro", type=int, default=1,
                    help="the thread count the driver will run at; written to "
                         "the header so the driver refuses a mismatched file")
    ap.add_argument("--case", default="dense",
                    choices=("dense", "corner", "zonal"))
    ap.add_argument("--nlev", type=int, default=10,
                    help="the level count the driver was compiled with; "
                         "written to the header so the driver refuses a "
                         "mismatched file")
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

    if args.nlev < 1:
        raise SystemExit("--nlev must be at least 1")
    fmu, fmv = per_mode_factors(modes)
    spd = driving_field(args.case, modes, args.nlev)
    spzrel = driving_field(args.case, modes, args.nlev)
    if args.case == "dense":
        # A second draw, so divergence and vorticity are not the same field:
        # they enter dv2uv through different matrices with different signs, and
        # driving them equal lets a swap between the two halves cancel.
        rng = np.random.default_rng(20260826)
        spzrel = rng.standard_normal((2, ncsp, args.nlev)).astype(LD)
        for k, (m, _n) in enumerate(modes):
            if m == 0:
                spzrel[1, k, :] = LD(0)
    pu, pv, spz, tscale = synthesise_uv(spd, spzrel, pmat, qmat, fmu, fmv,
                                        modes, nlon, nlat, args.nlev)
    dvfmax = float(max(np.max(np.abs(fmu)), np.max(np.abs(fmv)))
                   * max(np.max(np.abs(spd)), np.max(np.abs(spz))))

    header = np.array([nlat, nlon, ntru, ncsp, args.npro, nlat // args.npro,
                       scale, VERSION, args.nlev, PLAVOR, tscale, dvfmax],
                      dtype=np.float64)
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
        fmu.astype(np.float64).tofile(fh)
        fmv.astype(np.float64).tofile(fh)
        # component fastest, then mode, then level: pd(2,NCSP,NLEV)
        spd.astype(np.float64).transpose(2, 1, 0).reshape(-1).tofile(fh)
        spz.astype(np.float64).transpose(2, 1, 0).reshape(-1).tofile(fh)
        # component fastest, then m, then l, then v: pu(2,NLON/2,NLAT,NLEV)
        pu.astype(np.float64).transpose(3, 2, 1, 0).reshape(-1).tofile(fh)
        pv.astype(np.float64).transpose(3, 2, 1, 0).reshape(-1).tofile(fh)

    print(f"{args.res}  NLAT {nlat}  NLON {nlon}  NTRU {ntru}  modes {ncsp}  "
          f"NPRO {args.npro}  case {args.case}")
    print(f"  quadrature identities: |sum(w) - 2| = {dsum:.3e}, "
          f"worst Legendre moment {dmom:.3e}")
    print(f"  orthonormality residual of the reference table: {resid:.3e}")
    print(f"  the driving vector was scaled by {scale:.6e}, so the largest "
          f"Fourier coefficient is 1")
    print(f"  synthesis arm: {args.nlev} levels, plavor {PLAVOR}, "
          f"largest |u| {float(np.max(np.abs(pu))):.3e}, "
          f"largest |v| {float(np.max(np.abs(pv))):.3e}")
    print(f"  worst term-magnitude sum into one output {tscale:.3e}, "
          f"largest |factor|*|coefficient| {dvfmax:.3e}")
    print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
