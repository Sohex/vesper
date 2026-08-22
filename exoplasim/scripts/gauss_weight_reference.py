#!/usr/bin/env python3
"""Gauss-Legendre nodes and weights in extended precision, as an adjudicator.

Worldbuilding frame: a numerical reference for the Vesper climate model's
spectral transform. Nothing here is about the simulated planet.

WHY THIS EXISTS. `inigau` and SHTns disagree about the Gaussian quadrature
weights by up to 1e-10 near the poles, and numpy's `leggauss` disagrees with
both -- it sits between them. Three double-precision implementations that
disagree cannot settle which is right, and the weights enter every forward
transform the model performs, so the question is not academic. CLIM-61.

x86 `longdouble` carries a 64-bit mantissa, eps about 1.1e-19, which is eight
orders of margin on a disagreement at 1e-11. That is enough to adjudicate and
it needs no package this project does not already have; arbitrary precision
would be tidier and is not necessary.

THE METHOD is the textbook one, and the point is that it is textbook: Newton on
P_n evaluated by the three-term recurrence, then w = 2/((1-x^2) P'^2). What
`inigau` does instead is evaluate P_n as a TRIGONOMETRIC SERIES and take the
weight as z4(1-z^2)/z5^2, which squares the series' error into the weight.
"""
from __future__ import annotations

import argparse

import numpy as np

LD = np.longdouble


def legendre_pn(n: int, x):
    """P_n(x) and P'_n(x) by the three-term recurrence, in extended precision."""
    p0 = np.ones_like(x)
    p1 = x.copy()
    if n == 0:
        return p0, np.zeros_like(x)
    for k in range(2, n + 1):
        p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
    # derivative from the standard identity, exact where |x| < 1
    dp = n * (x * p1 - p0) / (x * x - LD(1))
    return p1, dp


def gauss_legendre(n: int, iters: int = 100):
    """Nodes descending from +1, and weights, matching the model's ordering."""
    j = np.arange(1, n // 2 + 1, dtype=np.int64)
    # Tricomi's asymptotic start, good to about 1e-5 and plenty for Newton
    x = np.cos(LD(np.pi) * (LD(4) * j - LD(1)) / (LD(4) * n + LD(2)))
    for _ in range(iters):
        p, dp = legendre_pn(n, x)
        step = p / dp
        x = x - step
        if np.max(np.abs(step)) < LD(1e-20):
            break
    _, dp = legendre_pn(n, x)
    w = LD(2) / ((LD(1) - x * x) * dp * dp)
    nodes = np.concatenate([x, -x[::-1]])
    wts = np.concatenate([w, w[::-1]])
    return nodes, wts


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("nlat", type=int, nargs="?", default=0,
                    help="number of Gaussian latitudes")
    ap.add_argument("--csv", action="store_true",
                    help="print every node and weight rather than a summary")
    ap.add_argument("--compare", metavar="CSV",
                    help="a dumped 'index,node,weight' file; print the worst "
                         "relative weight difference and nothing else, so a "
                         "shell can act on it")
    args = ap.parse_args()

    if args.compare:
        _, wts = gauss_legendre(args.nlat)
        mw = np.array([float(line.split(",")[2]) for line in open(args.compare)])
        ref = wts.astype(float)
        print(f"{float(np.max(np.abs((mw - ref) / ref))):.6e}")
        return 0

    if args.nlat <= 0:
        raise SystemExit("give a latitude count")
    if args.nlat % 2:
        raise SystemExit("this model's grids have an even latitude count")
    nodes, wts = gauss_legendre(args.nlat)

    if args.csv:
        for k, (x, w) in enumerate(zip(nodes, wts), start=1):
            print(f"{k},{x:.20e},{w:.20e}")
        return 0

    # A quadrature that integrates 1 exactly is a necessary check on the
    # weights, and one that is independent of how they were computed.
    print(f"NLAT {args.nlat}")
    print(f"  sum of weights - 2 = {float(np.sum(wts) - LD(2)):.3e}   (exact is 0)")
    ref = np.polynomial.legendre.leggauss(args.nlat)
    npw = ref[1][::-1]
    rel = np.abs((npw - wts.astype(float)) / wts.astype(float))
    print(f"  numpy leggauss, worst relative difference {float(np.max(rel)):.3e} "
          f"at latitude {int(np.argmax(rel)) + 1}")
    print(f"  first weight  {float(wts[0]):.17e}")
    print(f"  last weight   {float(wts[args.nlat // 2 - 1]):.17e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
