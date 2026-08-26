#!/usr/bin/env python3
"""Reconcile the two grounded values for `playa_clastic` albedo.

    python analysis/playa_albedo.py

## The disagreement

`playa_clastic`, "Playa mud / alluvial fan fill", is about a quarter of this
planet's land and the largest single lever on its energy balance. Two sources in
`references/INDEX.md`, both read, give different answers:

  Post et al. (2000)              52 pyranometer measurements over 26 US soils,
                                  0.3-2.8 um: mean 0.189, range 0.048-0.402.
                                  `vendor/orogen/js/lithology.js` cites this for
                                  moving the class from 0.30 to 0.19.

  Henderson-Sellers, Wilson (1983) salt playas and light sand deserts 0.28-0.44.

They are not measuring the same thing, and neither can be used raw.

## What each source is good for, and what it is not

**Post is the right KIND of measurement and the wrong POPULATION.** A field
pyranometer over a real surface is exactly what a model albedo wants; 26 US
agricultural soils, many of them dark and organic-rich, are not a sample of
playa mud and desert-varnished fan gravel. Their mean is a population statistic
being applied to a specific material.

**Henderson-Sellers is a compilation, and "salt playa" is our OTHER class.** A
salt-encrusted playa surface is `evaporite` here, which carries 0.50 already.
The clastic apron around and under the crust is what `playa_clastic` names.

**ECOSTRESS is the right population and the wrong preparation.** It carries
saline desert soils measured as directional hemispherical reflectance, which is
the geometry a model albedo wants, but as prepared laboratory samples rather than
crusted field surfaces -- the same trap `analysis/rock_albedo.py` documents for
powders against slabs.

## The reconciliation

Use each for what it measures well, and let the third correct the second.

1. Within ECOSTRESS, measured consistently, saline desert soils are BRIGHTER
   than the soil population mean. That ratio is a property of the material, not
   of the laboratory, because numerator and denominator share a preparation.
2. Post gives the field level for an average soil.
3. Multiply: field level times material ratio is the field level for this
   material. Then re-weight from the Sun to this star, which the spectra allow
   directly.

The preparation offset cancels in step 1 and never enters the answer. What does
not cancel is that Post's population and ECOSTRESS's population are not the same
soils, so the ratio carries whatever difference there is between "average US
agricultural soil" and "average ECOSTRESS soil". That is the residual
uncertainty and it is reported rather than buried.

## A third preparation axis that does NOT cancel: water content

Both populations are dry. ECOSTRESS's saline desert soils are prepared
laboratory samples, and Post's pyranometer runs were over field surfaces whose
water content is not reported here, so the answer this script returns is the DRY
endmember for `playa_clastic` and is staged as one. Playa mud is the class where
that matters most on this world: it is the largest single surface class on the
simulated land and it is by definition ground that floods and dries.

The wet endmember is not measured in either population. It has to come from a
wetting RATIO applied to the dry answer, and the one held ratio -- Penndorf
(1956) Table 1, clay soil 7.5/15 -- is a luminous-reflectance quantity that
brackets band 1 and leaves band 2 open. The sign is not safe to assume on the
neighbouring class either: the twenty-year MODIS record over the largest halite
pan in `references/INDEX.md` reads BRIGHTER in wet years than in dry ones.
`exoplasim/notes/soil-albedo-moisture.md` carries the argument, the magnitude,
and the four preconditions that would arm a moisture-dependent soil albedo.
"""

from __future__ import annotations

import argparse
import glob
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ECOSTRESS = ROOT / "references" / "ecospeclib-all"
SPECTRUM = ROOT / "exoplasim" / "inputs" / "stellarspectra" / "k25v_hr.dat"

H, C, KB = 6.626e-34, 2.998e8, 1.381e-23
SUN_TEFF = 5772.0

# Post et al. (2000), as recorded in references/INDEX.md.
POST_MEAN = 0.189
POST_RANGE = (0.048, 0.402)
POST_BAND = (0.3, 2.8)          # the band their pyranometer integrated over

# Aridisol suborders that are the material in question. Salorthids are the
# saline ones -- salt-affected soils of lake plains, which is a playa in the
# taxonomy of soil science.
PLAYA_SUBORDERS = ("salorthid", "gypsiorthid", "calciorthid", "camborthid",
                   "haplargid")


def planck(microns: np.ndarray, teff: float) -> np.ndarray:
    lam = microns * 1e-6
    return (2 * H * C ** 2 / lam ** 5) / (np.exp(H * C / (lam * KB * teff)) - 1)


def read_ecostress(path: Path):
    """Header fields and an ascending (wavelength_um, reflectance_percent) table."""
    head, rows = {}, []
    for line in path.read_text(errors="ignore").splitlines():
        if re.match(r"^\s*[\d.]+\s+[-\d.]+\s*$", line):
            a, b = line.split()[:2]
            rows.append((float(a), float(b)))
        elif ":" in line:
            key, _, value = line.partition(":")
            head.setdefault(key.strip(), value.strip())
    if not rows:
        return head, None
    arr = np.array(rows)
    return head, arr[arr[:, 0].argsort()]


def weighted(arr, star_wl, flux, band) -> float | None:
    """Reflectance integrated against a spectrum over a band. Percent to fraction."""
    wl, refl = arr[:, 0], arr[:, 1] / 100.0
    keep = (wl >= band[0]) & (wl <= band[1])
    if keep.sum() < 10:
        return None
    wl, refl = wl[keep], refl[keep]
    f = np.interp(wl, star_wl, flux)
    return float(np.trapezoid(refl * f, wl) / np.trapezoid(f, wl))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--spectrum", type=Path, default=SPECTRUM)
    ap.add_argument("--library", type=Path, default=ECOSTRESS)
    args = ap.parse_args()
    if not args.library.is_dir():
        raise SystemExit(
            f"{args.library} is absent. It is bulk reference data, excluded from "
            "git; unzip references/ecospeclib_*.zip beside this path.")

    data = np.loadtxt(args.spectrum, skiprows=1)
    star_wl, star_k = data[:, 0], data[:, 1]
    star_sun = planck(star_wl, SUN_TEFF)

    rows = []
    for path in sorted(glob.glob(str(args.library / "soil.*spectrum.txt"))):
        head, arr = read_ecostress(Path(path))
        if arr is None:
            continue
        # The Sun over Post's own band, so the population comparison is
        # like for like, and this star over the full shortwave for the answer.
        sun = weighted(arr, star_wl, star_sun, POST_BAND)
        k25v = weighted(arr, star_wl, star_k, (0.3, 4.0))
        if sun is None or k25v is None:
            continue
        suborder = Path(path).name.split(".")[2]
        rows.append((suborder, sun, k25v, head.get("Name", "")))

    if not rows:
        raise SystemExit("no soil spectra parsed")
    sun_all = np.array([r[1] for r in rows])
    playa = [r for r in rows if r[0] in PLAYA_SUBORDERS]
    sun_playa = np.array([r[1] for r in playa])
    k_playa = np.array([r[2] for r in playa])

    ratio = float(sun_playa.mean() / sun_all.mean())
    star_gain = float(k_playa.mean() / sun_playa.mean())
    reconciled = POST_MEAN * ratio * star_gain

    print(f"ECOSTRESS soils, sun-weighted over Post's {POST_BAND[0]}-{POST_BAND[1]} um band")
    print(f"  all soils        n={len(rows):3d}  mean {sun_all.mean():.3f}  "
          f"range {sun_all.min():.3f}-{sun_all.max():.3f}")
    print(f"  Post, field      n= 26  mean {POST_MEAN:.3f}  "
          f"range {POST_RANGE[0]:.3f}-{POST_RANGE[1]:.3f}")
    print(f"  lab over field   {sun_all.mean()/POST_MEAN:.2f}x, which is the "
          f"preparation offset that must NOT enter the answer")
    print()
    print(f"  desert/saline    n={len(playa):3d}  mean {sun_playa.mean():.3f}  "
          f"range {sun_playa.min():.3f}-{sun_playa.max():.3f}")
    print(f"  material ratio   {ratio:.3f}x the soil population, same preparation")
    print(f"  this star vs Sun {star_gain:.3f}x, from the same spectra")
    print()
    print(f"  reconciled       {POST_MEAN:.3f} x {ratio:.3f} x {star_gain:.3f} "
          f"= {reconciled:.3f}")
    print(f"  generator has    0.190")
    print()
    for suborder, sun, k25v, name in sorted(playa, key=lambda r: -r[1]):
        print(f"    {sun:.3f} sun  {k25v:.3f} K2.5V  {suborder:12s} {name[:44]}")


if __name__ == "__main__":
    main()
