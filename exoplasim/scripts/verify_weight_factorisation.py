#!/usr/bin/env python
"""Are legini's eight weight matrices really P and Q times separable factors?

CLIM-48 rests on this. `legini` builds eight `NCSP x NLPP` matrices and the
model touches every one each timestep, which is 48.4 MB a die at T127 and
114.9 MB at T170 against CCD1's 32 MB of L3 -- the measured cause of the load
imbalance in `notes/rank-imbalance-and-weight-traffic.md`. If all eight are one
of two arrays scaled by a per-MODE factor and a per-LATITUDE factor, then two
matrices plus a handful of vectors hold the same information and the resident
footprint falls fourfold.

This reproduces legini's recurrence, builds the eight matrices the way legini
does, and rebuilds each from the unfiltered P and Q plus its factors. Exact
agreement for the two stored raw; about one ulp for the six that are the same
product taken in a different order.

UNFILTERED on purpose. The filter fold put skspgp/skgpsp INSIDE the matrices,
and a filter can be exactly zero -- Cesaro at n = NTP1 -- so P is not
recoverable from qi by division. The factorisation has to keep P and Q raw, and
this runs with a filter that has a genuine zero in it, so a division-based
shortcut fails here rather than in production.

A NEGATIVE CONTROL drops the zonal wavenumber from qm's per-mode factor, which
is the shape of mistake this invites, and must be caught.
"""
import argparse
import sys

import numpy as np

NAMES = ("qi", "qj", "qc", "qu", "qv", "qe", "qq", "qm")


def legini_recurrence(zsin: float, zgwd: float, ntru: int, ncsp: int):
    """legmod.f90:legini's own recurrence, for one latitude."""
    zpli = np.zeros(ncsp)
    zpld = np.zeros(ncsp)
    zcsq = 1.0 - zsin * zsin
    f1m = np.sqrt(1.5)
    zpli[0] = np.sqrt(0.5)
    zpli[1] = f1m * zsin
    lm = 2
    for m in range(0, ntru + 1):
        if m > 0:
            lm += 1
            f2m = -f1m * np.sqrt(zcsq / (m + m))
            f1m = f2m * np.sqrt(m + m + 3.0)
            zpli[lm - 1] = f2m
            if lm < ncsp:
                lm += 1
                zpli[lm - 1] = f1m * zsin
                zpld[lm - 2] = -m * f2m * zsin
        amsq = m * m
        for n in range(m + 2, ntru + 1):
            lm += 1
            z1 = np.sqrt(((n - 1) ** 2 - amsq) / (4 * (n - 1) ** 2 - 1))
            z2 = zsin * zpli[lm - 2] - z1 * zpli[lm - 3]
            zpli[lm - 1] = z2 * np.sqrt((4 * n * n - 1) / (n * n - amsq))
            zpld[lm - 2] = (1 - n) * z2 + n * z1 * zpli[lm - 3]
        if lm < ncsp:
            z3 = np.sqrt((ntru * ntru - amsq) / (4 * ntru * ntru - 1))
            zpld[lm - 1] = (-ntru * zsin * zpli[lm - 1]
                            + (ntru + ntru + 1) * zpli[lm - 2] * z3)
        else:
            zpld[lm - 1] = -ntru * zsin * zpli[lm - 1]
    return zpli, zpld, zgwd / zcsq


def compare(ntru: int, nlat: int, break_qm: bool):
    ntp1 = ntru + 1
    ncsp = ntp1 * (ntp1 + 1) // 2
    mu, gw = np.polynomial.legendre.leggauss(nlat)
    mu, gw = mu[::-1], gw[::-1]
    modes = [(m, n) for m in range(ntru + 1) for n in range(m, ntru + 1)]

    # skspgp with an exact zero at the top mode, and skgpsp exponential
    a = np.array([1.0 - (n + 1) / ntp1 for (_, n) in modes])
    b = np.array([np.exp(-((n + 1) / ntru) ** 4) for (_, n) in modes])
    znn1 = np.array([0.0 if n == 0 else 1.0 / (n * (n + 1)) for (_, n) in modes])
    mm = np.array([float(m) for (m, _) in modes])
    nn = np.array([n * (n + 1) * 0.5 for (_, n) in modes])

    per_mode = {"qi": a, "qj": a, "qc": b, "qe": b,
                "qu": a * znn1 * mm, "qv": a * znn1,
                "qq": b * nn, "qm": b * (np.ones_like(mm) if break_qm else mm)}
    base_is_p = {"qi": True, "qj": False, "qc": True, "qu": True,
                 "qv": False, "qe": False, "qq": True, "qm": True}

    worst = {k: 0.0 for k in NAMES}
    for l in range(nlat):
        zpli, zpld, zgwdcsq = legini_recurrence(mu[l], gw[l], ntru, ncsp)
        zgwd = gw[l]
        legini = {
            "qi": zpli * a, "qj": zpld * a,
            "qc": (zpli * zgwd) * b,
            "qu": (zpli * znn1 * mm) * a,
            "qv": (zpld * znn1) * a,
            "qe": (zpld * zgwdcsq) * b,
            "qq": (zpli * zgwdcsq * nn) * b,
            "qm": (zpli * zgwdcsq * mm) * b,
        }
        per_lat = {"qi": 1.0, "qj": 1.0, "qu": 1.0, "qv": 1.0,
                   "qc": zgwd, "qe": zgwdcsq, "qq": zgwdcsq, "qm": zgwdcsq}
        for k in NAMES:
            got = (zpli if base_is_p[k] else zpld) * per_mode[k] * per_lat[k]
            scale = max(np.abs(legini[k]).max(), 1.0)
            worst[k] = max(worst[k], float(np.abs(legini[k] - got).max() / scale))
    return worst, ncsp


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ntru", type=int, default=12)
    ap.add_argument("--nlat", type=int, default=16)
    ap.add_argument("--tol", type=float, default=1e-14,
                    help="the six rebuilt matrices are the same products in a "
                         "different order, so a few ulp is the expected answer "
                         "and anything larger is a different quantity")
    args = ap.parse_args()

    worst, ncsp = compare(args.ntru, args.nlat, break_qm=False)
    print(f"NTRU {args.ntru}, {ncsp} modes, {args.nlat} latitudes, "
          f"a filter carrying an exact zero")
    failures = 0
    for k in NAMES:
        d = worst[k]
        state = "exact" if d == 0.0 else ("reassociated" if d <= args.tol else "DIFFERENT")
        if d > args.tol:
            failures += 1
        print(f"  {k}  {d:.3e}  {state}")

    print()
    print("the negative control, which must fail")
    broken, _ = compare(args.ntru, args.nlat, break_qm=True)
    if broken["qm"] > args.tol:
        print(f"  [  ok  ] dropping m from qm's per-mode factor is caught "
              f"({broken['qm']:.3e})")
    else:
        failures += 1
        print("  [ FAIL ] the control PASSES, so this check proves nothing")

    print()
    print("stored: 2 matrices of NCSP x NLPP, 6 vectors of NCSP, 3 scalars a latitude")
    print("against 8 matrices of NCSP x NLPP -- a quarter of the resident traffic")
    print(f"{failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
