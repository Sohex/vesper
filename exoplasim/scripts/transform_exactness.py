#!/usr/bin/env python3
"""Is the model's Gaussian quadrature exact for the terms `calcgp` actually forms?

    python exoplasim/scripts/transform_exactness.py <run_dir> --time 20 --level 8

Worldbuilding frame: a numerical check on the Vesper project's climate model.
The fields are model arrays.

THE QUESTION. `world-bxr` measures an energy sink in the dynamical core. One
candidate is that the transform back to spectral space mishandles the nonlinear
terms. A spectral model on a Gaussian grid with NLON = 3*NTRU + 1 -- which is
every rung on this ladder -- transforms a product of TWO band-limited fields
exactly. It does NOT dealias a product of three, and `calcgp` forms products of
three: `zsdotp` is itself a sum of products and multiplies a vertical difference
to make `gtd`, `gud` and `gvd`. It also forms terms that are not polynomials at
all, `rcsq = 1/cos^2(phi)` and `gp = exp(gp)`, which no dealiasing rule covers.

WHAT IS TESTED, AND WHY IT IS THE GLOBAL MEAN. The sink is a global-mean
enthalpy loss, and the global-mean temperature tendency comes from ONE gridpoint
field: in `mktend` the (0,0) mode picks up neither the `fmm` term, which carries
a factor m, nor the `qmat` term, whose derivative vanishes there. So the sink can
only come through the transform if the quadrature that forms that mean is
inexact. That is an identity with a right answer.

HOW. Every field is analysed to T42 from the model's own gridpoint output, then
synthesised EXACTLY onto a grid three times finer in each direction, where the
quadrature is exact to far higher degree. The same term is formed on both grids
from the same coefficients and its global mean compared. The difference is the
quadrature error, measured rather than bounded.

THE CONTROLS COME FIRST AND CAN FAIL. A pure quadratic must agree to round-off,
because the grid is built to dealias it; a cubic must also agree in the GLOBAL
MEAN, because the mean needs only degree 126 against a 64-point rule exact to
127. If either misses, the apparatus is wrong and no other row here means
anything. The round trip is checked too: `ua` is NOT band-limited -- the model
carries U = u*cos(phi) and the postprocessor divides -- so the transform is
verified on U, where the error must fall to the float32 the output is stored in.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.special import sph_harm_y

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _paths  # noqa: F401

ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "exoplasim" / "analysis"

ROUNDOFF = 1e-10      # a control above this means the apparatus is wrong
STORAGE = 1e-6        # the output is written as float32; a round trip cannot beat it


class SHT:
    """Triangular spherical-harmonic transform on a Gaussian grid, in its own basis.

    Its own basis and not the model's, deliberately: it acts on gridpoint output,
    so no coefficient convention has to be matched and none can be got wrong.
    """

    def __init__(self, ntru: int, nlat: int, nlon: int):
        self.ntru, self.nlat, self.nlon = ntru, nlat, nlon
        mu, w = np.polynomial.legendre.leggauss(nlat)
        self.mu, self.w = mu[::-1].copy(), w[::-1].copy()   # north to south
        theta = np.arccos(self.mu)
        self.mmax = min(ntru, nlon // 2)
        self.L = [np.array([np.real(sph_harm_y(n, m, theta, 0.0))
                            for n in range(m, ntru + 1)])
                  for m in range(self.mmax + 1)]

    def analyse(self, f):
        F = np.fft.rfft(f, axis=-1) / self.nlon
        return [2.0 * np.pi * np.einsum("...y,ny->...n", F[..., m] * self.w, self.L[m])
                for m in range(self.mmax + 1)]

    def synthesise(self, c, nlon=None):
        nlon = nlon or self.nlon
        F = np.zeros(c[0].shape[:-1] + (self.nlat, nlon // 2 + 1), dtype=complex)
        for m in range(len(c)):
            F[..., m] = np.einsum("...n,ny->...y", c[m], self.L[m])
        return np.fft.irfft(F, n=nlon, axis=-1) * nlon


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--time", type=int, default=0, help="output sample index")
    ap.add_argument("--level", type=int, default=8, help="model level index")
    ap.add_argument("--refine", type=int, default=3,
                    help="how much finer the reference grid is, per direction")
    ap.add_argument("--out", type=Path, default=ANALYSIS / "transform_exactness.json")
    args = ap.parse_args()

    import netCDF4 as nc
    run = args.run_dir if args.run_dir.is_absolute() else ROOT / args.run_dir
    manifest = json.loads((run / "run_manifest.json").read_text())
    ntru = int(str(manifest["physical"]["resolution"]).lstrip("Tt"))
    files = sorted(run.glob("MOST.*.nc"))
    if not files:
        raise SystemExit(f"no MOST.*.nc in {run}")
    d = nc.Dataset(files[0])
    nlat = d.dimensions["lat"].size
    nlon = d.dimensions["lon"].size
    if nlon != 3 * ntru + 1 + (3 * ntru + 1) % 2:
        print(f"note: NLON={nlon} against 3*NTRU+1={3*ntru+1}; the dealiasing "
              f"statement in this script assumes the model's own rule")

    coarse = SHT(ntru, nlat, nlon)
    fine = SHT(ntru, nlat * args.refine, nlon * args.refine)
    wc, wf = coarse.w, fine.w
    nlon_f = nlon * args.refine
    muc, muf = coarse.mu, fine.mu
    cosc = np.cos(np.arcsin(muc))[:, None]

    def both(field):
        c = coarse.analyse(field)
        return coarse.synthesise(c), fine.synthesise(c, nlon=nlon_f)

    def gmean(f, w):
        return float((f.mean(axis=-1) * w).sum() / w.sum())

    k, lev = args.time, args.level
    ua = np.asarray(d.variables["ua"][:], float)[k, lev]
    va = np.asarray(d.variables["va"][:], float)[k, lev]
    ta = np.asarray(d.variables["ta"][:], float)[k, lev]
    ps = np.asarray(d.variables["ps"][:], float)[k] * 100.0
    d.close()

    U, V, lnps = ua * cosc, va * cosc, np.log(ps)
    # THE ROUND TRIP, on the field the model actually carries.
    trip = float(np.abs(coarse.synthesise(coarse.analyse(U)) - U).max()
                 / np.abs(U).max())
    trip_u = float(np.abs(coarse.synthesise(coarse.analyse(ua)) - ua).max()
                   / np.abs(ua).max())

    Uc, Uf = both(U); Vc, Vf = both(V); Tc, Tf = both(ta); Lc, Lf = both(lnps)
    mv = np.fft.rfftfreq(nlon, 1 / nlon)
    cl = coarse.analyse(lnps)
    cd = [1j * m * cl[m] for m in range(len(cl))]
    GLc = coarse.synthesise(cd, nlon=nlon)
    GLf = fine.synthesise(cd, nlon=nlon_f)
    rc, rf = 1.0 / (1.0 - muc ** 2), 1.0 / (1.0 - muf ** 2)

    rows = [
        ("CONTROL quadratic  T*U", Tc * Uc, Tf * Uf, True),
        ("CONTROL cubic      U*V*T", Uc * Vc * Tc, Uf * Vf * Tf, True),
        ("gke      U*U + V*V", Uc * Uc + Vc * Vc, Uf * Uf + Vf * Vf, False),
        ("zvgpg    rcsq*U*dlnps/dlam",
         rc[:, None] * Uc * GLc, rf[:, None] * Uf * GLf, False),
        ("gp       exp(lnps)", np.exp(Lc), np.exp(Lf), False),
        ("rcsq*(U*U+V*V)", rc[:, None] * (Uc * Uc + Vc * Vc),
         rf[:, None] * (Uf * Uf + Vf * Vf), False),
    ]

    print(f"{run.name}, T{ntru} on {nlat}x{nlon}, reference {fine.nlat}x{nlon_f}")
    print(f"  round trip on U = ua*cos(phi) : {trip:.2e}  "
          f"(must reach the float32 the output is stored in, {STORAGE:g})")
    print(f"  round trip on ua itself       : {trip_u:.2e}  "
          f"(ua is NOT band-limited; this one is expected to be large)")
    print()
    out, failed = [], []
    for name, fc, ff, control in rows:
        a, b = gmean(fc, wc), gmean(ff, wf)
        rms = float(np.sqrt((ff ** 2).mean()))
        rel = abs(a - b) / rms if rms else float("nan")
        flag = ""
        if control and rel > ROUNDOFF:
            flag = "  <-- CONTROL FAILED"
            failed.append(name)
        print(f"  {name:28s} coarse {a:>13.6e}  fine {b:>13.6e}  "
              f"error/rms {rel:>9.2e}{flag}")
        out.append({"term": name, "coarse": a, "fine": b,
                    "relative_error": rel, "control": control})

    verdict = ("APPARATUS INVALID: a control missed round-off, so no row here "
               "means anything" if failed else
               "the quadrature is exact for every term formed here")
    print(f"\n  {verdict}")
    if trip > STORAGE:
        print("  NOTE: the round trip on U did not reach storage precision either")

    payload = {"note": "quadrature exactness of the model's nonlinear terms, for "
                       "world-bxr. exoplasim/scripts/transform_exactness.py",
               "generated": datetime.now(timezone.utc).isoformat(),
               "run": run.name, "truncation": ntru, "grid": [nlat, nlon],
               "reference_grid": [fine.nlat, nlon_f],
               "time_index": k, "level_index": lev,
               "round_trip_U": trip, "round_trip_ua": trip_u,
               "roundoff_threshold": ROUNDOFF, "terms": out,
               "controls_failed": failed, "verdict": verdict}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"  wrote {args.out}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
