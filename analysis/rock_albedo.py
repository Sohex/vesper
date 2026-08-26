#!/usr/bin/env python3
"""Solar-weighted rock albedo under this world's star, from measured spectra.

    python analysis/rock_albedo.py

This is the derivation behind the `albedo` column of `ROCK_CLASSES` in
`vendor/orogen/js/lithology.js`, which `build_surface_albedo.py` turns into an
ExoPlaSim boundary condition. The numbers were computed once and pasted into that
table; without this script the table is a list of values whose provenance is a
sentence in an audit. Anything that changes -- a new spectrum, a different star,
a rock class added -- wants to re-run it rather than interpolate by eye.

## Why a K dwarf changes the answer

An albedo is a reflectance integrated against the light falling on it. Published
rock albedos are integrated against the Sun. This star is a K2.5V at 4965 K, so
more of its output sits at longer wavelengths, and a mineral with absorption
features in the near infrared is darker here than its published value. The
integration is against `exoplasim/inputs/stellarspectra/k25v_hr.dat`.

## Powders are not surfaces, and the correction is lithology-dependent

The trap this script exists to avoid. Logan (1973) and Hunt (1982) measure
0-74 um powders in a fairy-castle packing, not solid rock, and a powder is far
brighter than the slab it came from. Paragas (2025) measures the SAME rock as
slab, crushed and powder, so the conversion is measurable rather than assumed --
and it is not a constant. It runs about 2.2x for basalt against 4.1x for granite,
so powder data compresses the felsic/mafic contrast and must not be carried
across it. See `notes/audits/orogen-lithology.md`.

The data lives in `references/poseidon_surface_albedo/`, from the POSEIDON
retrieval code's surface reflectivity database, which reports
directional-hemispherical reflectance -- the geometry a model albedo actually
wants.

## Every number here is a DRY endmember, and that is a scope, not a caveat

Each spectrum integrated below was measured on a prepared laboratory sample --
slab, crushed or powder -- so each albedo is the reflectance of that material at
or near zero water content. That is exactly the state the modelled land column's
surface layer reaches when it empties, which is why the two are staged as one
pair of ends rather than as a level and an offset.

The SATURATED endmember of the same material is
`analysis/soil_albedo_wetting.py`. It takes these same spectra and applies the
closed-form wet-from-dry relation of Lekner and Dorf (1988) and of Twomey,
Bohren and Mergenthaler (1986) per WAVELENGTH, which is what gives a band-2 wet
endmember from band-2 dry data rather than from a broadband ratio. The two
mechanisms are the two arms of a bracket and not two estimates of one number.
`exoplasim/scripts/build_surface_albedo.py` stages both ends and
`exoplasim/notes/soil-albedo-moisture.md` carries the argument.
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
POSEIDON = ROOT / "references" / "poseidon_surface_albedo"
SPECTRUM = ROOT / "exoplasim" / "inputs" / "stellarspectra" / "k25v_hr.dat"

H, C, KB = 6.626e-34, 2.998e8, 1.381e-23
SUN_TEFF = 5772.0


def planck(microns: np.ndarray, teff: float) -> np.ndarray:
    lam = microns * 1e-6
    return (2 * H * C ** 2 / lam ** 5) / (np.exp(H * C / (lam * KB * teff)) - 1)


def weighted_albedo(path: Path, star_wl: np.ndarray, star_flux: np.ndarray
                    ) -> tuple[float, float, float]:
    """Reflectance integrated against a stellar spectrum, over their overlap."""
    table = np.loadtxt(path)
    wl, refl = table[:, 0], table[:, 1]
    order = np.argsort(wl)
    wl, refl = wl[order], refl[order]
    lo, hi = max(wl.min(), 0.3), min(wl.max(), 4.0)
    keep = (wl >= lo) & (wl <= hi)
    wl, refl = wl[keep], refl[keep]
    flux = np.interp(wl, star_wl, star_flux)
    return float(np.trapezoid(refl * flux, wl) / np.trapezoid(flux, wl)), lo, hi


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spectrum", type=Path, default=SPECTRUM)
    args = ap.parse_args()

    if not POSEIDON.is_dir():
        raise SystemExit(
            f"{POSEIDON} is absent. It is bulk reference data, excluded from git "
            f"by .gitignore; refetch it from the POSEIDON repository.")
    data = np.loadtxt(args.spectrum, skiprows=1)
    star_wl, star_flux = data[:, 0], data[:, 1]
    sun_flux = planck(star_wl, SUN_TEFF)

    base = POSEIDON / "Paragas2025-P25"
    print("Paragas 2025: the same rock as slab, crushed and powder, K2.5V weighted")
    print(f"{'rock':40}{'slab':>8}{'crushed':>9}{'powder':>8}{'powder/slab':>13}")
    ratios = []
    rocks = sorted({os.path.basename(p).replace("_slab_P25.txt", "")
                    for p in glob.glob(str(base / "*_slab_P25.txt"))})
    for rock in rocks:
        got = {}
        for form in ("slab", "crushed", "powder"):
            path = base / f"{rock}_{form}_P25.txt"
            if path.is_file():
                got[form] = weighted_albedo(path, star_wl, star_flux)[0]
        if "slab" in got and "powder" in got:
            ratio = got["powder"] / got["slab"]
            ratios.append(ratio)
            print(f"{rock[:39]:40}{got['slab']:8.3f}"
                  f"{got.get('crushed', float('nan')):9.3f}"
                  f"{got['powder']:8.3f}{ratio:13.2f}")
    if ratios:
        arr = np.array(ratios)
        print(f"\npowder/slab: median {np.median(arr):.2f}, "
              f"range {arr.min():.2f}-{arr.max():.2f}, n={len(arr)}")
        print("NOT a constant. It is lithology-dependent, which is why powder\n"
              "ratios must not be carried across the felsic/mafic divide.")

    print("\nSurface types, Sun against this star")
    print(f"{'surface':34}{'Sun':>8}{'K2.5V':>8}{'delta':>8}")
    for source, pattern in (("Hu2012-H12", "*_H12.txt"),
                            ("Hammond2025-H25", "*_H25.txt")):
        for path in sorted(glob.glob(str(POSEIDON / source / pattern))):
            name = os.path.basename(path).rsplit("_", 1)[0]
            a_sun, _, _ = weighted_albedo(Path(path), star_wl, sun_flux)
            a_star, _, _ = weighted_albedo(Path(path), star_wl, star_flux)
            print(f"{name[:33]:34}{a_sun:8.3f}{a_star:8.3f}{a_star - a_sun:+8.3f}")


if __name__ == "__main__":
    main()
