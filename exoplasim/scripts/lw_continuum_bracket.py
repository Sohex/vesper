#!/usr/bin/env python
"""How large `th2oc`, the model's entire longwave water vapour continuum, can be.

WORLDBUILDING FRAME. Vesper is a simulated super-Earth. Every number below is a
property of a toy GCM's longwave scheme as configured for that simulated planet.
Nothing here is a measurement of anything outside the simulation.

THE QUESTION. `radmod`'s `lwr` builds its clear-sky absorptance from Sasamori
(1968), whose broadband fits cover water vapour, CO2 and ozone bands and carry
NO window absorption at all. The single line

    if (th2oc > 0.) zah2o = min(zah2o + (1 - exp(-th2oc*zsumwv)), 1)

is therefore the whole of the modelled longwave continuum. `th2oc = 0.024` is
declared with the comment "absorption coefficient h2o continuum (lwr)" and no
unit, no source and no derivation anywhere in the tree; it is the fourth member
of the per-truncation `jtune` table upstream carried and the only one of the
four this project has never re-weighted or bracketed. world-2esd.

WHAT THIS SCRIPT DOES, and what it does NOT do. It does not source the value --
that needs a correlated-k or line-by-line calculation with the MT_CKD continuum
on this path, and `references/INDEX.md` records that Mlawer et al. (2012) does
not supply an evaluable continuum because the coefficients ship as data with
LBLRTM. What it DOES is derive a hard upper bound and measure where the
inherited value sits inside it, which is what turns an opaque constant into a
declared one with a bracket that can be swept.

THE BOUND, and it is physics rather than a fit. The term adds broadband
absorptance that Sasamori's fits leave out, and what they leave out is the
window. So the addition cannot exceed the share of the emitted Planck flux that
the window carries: above that, the continuum absorbs more than the spectral
region it stands in for contains. At the largest pressure-weighted water path
the model reaches,

    1 - exp(-th2oc * w_max) <= f_window(T)   =>   th2oc <= -ln(1 - f_window)/w_max

THE PATH IS THE MODEL'S OWN, not an assumption. `zsumwv` is not the geometric
column: `zq = zfh2o * sigma*dsigma/ga/1e5 * ps**2 * q` carries `ps**2`, so what
Sasamori's fits receive is an amount times a pressure. This script rebuilds that
sum from a climatology's own `hus` and `ps` on the model's own sigma grid.

Run it as `python exoplasim/scripts/lw_continuum_bracket.py`.

IT WRITES NOTHING. What it produces is a bracket, and a bracket belongs in the
declaration it bounds and in the arms registered against it -- `radmod.f90`'s
`th2oc` declaration and `exoplasim/notes/forcing-bundle-predictions.md`. An
artifact no step in `config/pipeline.yaml` generates does not exist.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset
from scipy.integrate import quad

import _paths  # noqa: F401  anchors every path on this file and adds lib/

RADMOD = _paths.MODEL_SRC / "plasim/src/radmod.f90"
DEFAULT_CLIM = _paths.ANALYSIS / "climatology" / "bootstrap_regular_climatology.nc"

H = 6.62607015e-34
C = 2.99792458e8
KB = 1.380649e-23
SIGMA_SB = 5.670374419e-8

# Sasamori (1968) fits the 6.3 um band and the rotational band and stops. The
# window between them is where the continuum lives, and 8 to 12 um is the
# interval the continuum literature uses for it; 8 to 13 um is the wider reading
# and is carried as the looser bound.
WINDOWS_UM = ((8.0, 12.0), (8.0, 13.0))


def read_constant(name: str) -> float:
    """Take a literal out of the model source rather than restating it here."""
    text = RADMOD.read_text(encoding="utf-8", errors="replace")
    for pattern in (rf"^\s*real\s*::\s*{name}\s*=\s*([0-9.eEdD+-]+)",
                    rf"^\s*parameter\s*\(\s*{name}\s*=\s*([0-9.eEdD+-]+)\s*\)"):
        m = re.search(pattern, text, re.MULTILINE)
        if m:
            return float(m.group(1).replace("d", "e").replace("D", "e"))
    raise SystemExit(f"{name} not found in {RADMOD}")


def planck_fraction(lo_um: float, hi_um: float, t_k: float) -> float:
    def b(lam: float) -> float:
        return 2 * H * C ** 2 / lam ** 5 / (math.exp(H * C / (lam * KB * t_k)) - 1)
    v, _ = quad(lambda lam: math.pi * b(lam), lo_um * 1e-6, hi_um * 1e-6, limit=200)
    return v / (SIGMA_SB * t_k ** 4)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--climatology", type=Path, default=DEFAULT_CLIM,
                    help="the field the model's own water path comes from")
    args = ap.parse_args()

    th2oc = read_constant("th2oc")
    zfh2o = read_constant("zfh2o")
    cfg = yaml.safe_load(_paths.CONFIG.read_text(encoding="utf-8"))
    ga = None
    for block in cfg.values():
        if isinstance(block, dict) and "gravity_m_s2" in block:
            ga = float(block["gravity_m_s2"])
    if ga is None:
        raise SystemExit("gravity_m_s2 not found in config/planet.yaml")

    ds = Dataset(args.climatology)
    hus = np.asarray(ds.variables["hus"][:])
    ps_hpa = np.asarray(ds.variables["ps"][:])
    lev = np.asarray(ds.variables["lev"][:], dtype=float)

    sig = lev / lev.max() if lev.max() > 1.5 else lev
    half = np.empty(len(sig) + 1)
    half[0], half[-1] = 0.0, 1.0
    half[1:-1] = 0.5 * (sig[:-1] + sig[1:])
    dsig = np.diff(half)

    # THE CHECK THAT CAN FAIL: the layer thicknesses of a sigma coordinate sum
    # to exactly 1. A grid rebuilt wrongly, or a `lev` axis that is not sigma,
    # would not, and every path below would then be wrong by that factor.
    closure = abs(dsig.sum() - 1.0)
    if closure > 1e-6:
        raise SystemExit(f"sigma thicknesses sum to {dsig.sum()}, not 1")

    ps_pa = ps_hpa * 100.0
    w = np.zeros(hus.shape[0:1] + hus.shape[2:])
    for j in range(len(sig)):
        w += zfh2o * sig[j] * dsig[j] / ga / 1.0e5 * ps_pa ** 2 * hus[:, j]

    w_max = float(w.max())
    w_mean = float(w.mean())

    print(f"model source   {RADMOD.relative_to(_paths.PROJECT_ROOT)}")
    print(f"climatology    {args.climatology.relative_to(_paths.PROJECT_ROOT)}")
    print(f"th2oc          {th2oc}   zfh2o {zfh2o}   ga {ga} m/s2")
    print(f"sigma closure  |sum(dsigma) - 1| = {closure:.2e}   PASS")
    print()
    print("THE MODEL'S OWN PRESSURE-WEIGHTED WATER PATH, g/cm2 at a bar")
    print(f"  min {w.min():.4f}   mean {w_mean:.4f}   "
          f"p99 {np.percentile(w, 99):.4f}   max {w_max:.4f}")
    print()
    print("WHAT THE INHERITED COEFFICIENT CLAIMS, as broadband absorptance")
    print(f"  at the mean path  {1 - math.exp(-th2oc * w_mean):.4f}")
    print(f"  at the max path   {1 - math.exp(-th2oc * w_max):.4f}")
    print()
    print("THE CEILING, from the window's share of the emitted Planck flux")
    for lo, hi in WINDOWS_UM:
        print(f"  window {lo:g}-{hi:g} um")
        for t in (250.0, 273.0, 288.0, 300.0):
            f = planck_fraction(lo, hi, t)
            ceiling = -math.log(1.0 - f) / w_max
            print(f"    T = {t:5.0f} K   f_window {f:.4f}   "
                  f"th2oc <= {ceiling:.4f}   inherited is {th2oc / ceiling:.2f} of it")


if __name__ == "__main__":
    main()
