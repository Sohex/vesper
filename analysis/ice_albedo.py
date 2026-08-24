#!/usr/bin/env python3
"""Verify ExoPlaSim's own star-weighting of the snow, glacier and sea ice albedo.

    python analysis/ice_albedo.py

`radmod.f90:600-790` already does what PHYS-14 was opened to do: it carries the
ECOSTRESS snow, ice, water and frost blends as Fortran arrays over 965
wavelengths, weights each by the configured stellar spectrum, splits at 0.75 um
and assigns the two-band `dsnowalb`, `dsnowalbmn`, `dsnowalbmx`, `dicealbmn`,
`dicealbmx` and `dglacalbmn`. `radini` runs at `plasim.f90:513`, BEFORE
`surfini` reads any surface namelist, so the flat two-element declarations at
`plasimmod.f90:428-432` are defaults that never reach the radiation.

So there is nothing here to derive and no key to add. What there is, is a
normalisation that mixes two grids, and this script exists to measure it.

## The defect

The numerator is integrated on the LOW-resolution grid and the denominator on
the HIGH-resolution one:

- `a1`, `a2` sum `bb3(k)*blend(k)` over `wavelengths`, the 965-point grid the
  surface blends live on, which starts at 0.34 um and ends at 14.01 um.
- `z1`, `z2` sum `bb1`, `bb2` over `wv1`, `wv2`, the 2048-point grid, band 1
  running from the file's start with everything below `minwavel` zeroed, band 2
  running to 100 um, and the band-edge interval added to `z1`.
- `zdenom1 = 0.01/z1` then divides one by the other.

A ratio of integrals is a reflectance only when both cover the same interval
with the same spectrum. These do not, so the quotient is not one. It is the
same class of defect as CLAUDE.md rule 3 -- two labellings of one axis combined
without reconciliation -- on wavelength rather than longitude, and the same
family as the grid mismatch already found in this module's Rayleigh reference.

Band 2 nearly escapes: this star has little flux beyond 14.01 um, so the two
band-2 intervals almost coincide. Band 1 does not escape, and every surface is
one-signed DARKER than a self-consistent integration gives.

## The check that can fail

`RECORDED` holds what the model itself printed under "Finalized Albedos". This
script reproduces those numbers from source, and refuses if it cannot. That is
what makes the alternative column below a measurement of the defect rather than
a second opinion about it: the reproduction is exact, so the difference between
the two columns is the normalisation and nothing else.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib import stellar  # noqa: E402
from lib.paths import PROJECT_ROOT, rel  # noqa: E402

ROOT = PROJECT_ROOT
SURFACESPECS = ROOT / "vendor" / "exoplasim" / "exoplasim" / "surfacespecs.py"
OUTPUT = ROOT / "analysis" / "ice_albedo.json"

BAND1_LAST_INDEX = 41
"""`radmod.f90:663`, `do k=1,41`. On this grid that is 0.34 to 0.74 um."""

REPRODUCTION_TOLERANCE = 1.0e-9
"""How far the reproduction may sit from what the model printed.

Not a physical tolerance. The arithmetic is the same operations in the same
order on the same inputs, so it either reproduces to floating-point noise or
the source has changed under this script and the comparison means nothing.
Set at the print precision of `MOST_DIAG`, which is far looser than the 1e-16
actually achieved, so a real change cannot hide inside it.
"""

RECORDED = {
    # "Finalized Albedos" from exoplasim/runs/run_2b20e3324bb0/MOST_DIAG.00001,
    # a T42 run on the k25v spectrum. Band 1 then band 2.
    "dsnowalbmx": (0.96819269058784208, 0.58528588178241947),
    "dsnowalbmn": (0.50026688807985986, 0.27258907520749881),
    "dicealbmx": (0.87855514368045551, 0.47288429368749935),
    "dicealbmn": (0.62630199666337949, 0.33947282395248929),
    "dglacalbmn": (0.75242868443925126, 0.40617853940042364),
}

SURFACES = [
    ("iceblendmax", "dsnowalbmx", 0.8, "snow, maximum"),
    ("iceblendmin", "dsnowalbmn", 0.4, "snow, minimum"),
    ("seaicemax", "dicealbmx", 0.7, "sea ice, maximum"),
    ("seaicemin", "dicealbmn", 0.5, "sea ice, minimum"),
    ("glacalbmin", "dglacalbmn", 0.6, "glacier ice, minimum"),
]


def load_surfacespecs():
    spec = importlib.util.spec_from_file_location("surfacespecs", SURFACESPECS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _trapz_pairs(x, y, lo, hi):
    """`radmod`'s own loop: sum of trapezoids over the 1-based range [lo, hi)."""
    return float(np.trapezoid(y[lo - 1:hi], x[lo - 1:hi]))


def low_res_weighting(wl_um, bb3, blend_percent):
    """`a1`, `a2` and their self-consistent denominators, all on the 965 grid."""
    wl_m = wl_um * 1.0e-6
    k = BAND1_LAST_INDEX
    n = wl_um.size
    a1 = _trapz_pairs(wl_m, bb3 * blend_percent, 1, k + 1)
    a2 = _trapz_pairs(wl_m, bb3 * blend_percent, k + 1, n)
    z1 = _trapz_pairs(wl_m, bb3, 1, k + 1)
    z2 = _trapz_pairs(wl_m, bb3, k + 1, n)
    return a1, a2, z1, z2


def hi_res_denominators(path: Path):
    """`z1`, `z2` exactly as `radmod.f90:519-545` builds them."""
    wavelength, flux = stellar.read_hires(path)
    wv1, bb1, wv2, bb2 = stellar._bands(wavelength, flux)
    z1 = float(np.trapezoid(bb1, wv1))
    z2 = float(np.trapezoid(bb2, wv2))
    z1 += 0.5 * (bb1[-1] + bb2[0]) * (wv2[0] - wv1[-1])
    return z1, z2


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--spectrum", default=None)
    ap.add_argument("--output", type=Path, default=OUTPUT)
    args = ap.parse_args()

    specs = load_surfacespecs()
    wvl = specs.wvl
    lowres, hires = stellar.spectrum_paths(args.spectrum)
    if wvl.size != 965:
        raise SystemExit(f"surfacespecs ships {wvl.size} wavelengths, not 965")
    if abs(wvl[BAND1_LAST_INDEX - 1] - 0.74) > 1e-9 or abs(wvl[BAND1_LAST_INDEX] - 0.75) > 1e-9:
        raise SystemExit(
            f"the blend grid no longer changes band at index {BAND1_LAST_INDEX}: "
            f"{wvl[BAND1_LAST_INDEX - 1]} and {wvl[BAND1_LAST_INDEX]} um. "
            "radmod's hardcoded do-loop bounds would split somewhere else.")

    low = np.loadtxt(lowres, skiprows=1)
    if not np.allclose(low[:, 0], wvl):
        raise SystemExit(
            "the stellar spectrum and the surface blends are on different "
            "wavelength grids, so radmod is multiplying bb3(k) by blend(k) at "
            "mismatched wavelengths. That is a larger defect than this script "
            "measures and must be settled first.")
    bb3 = low[:, 1]
    zh1, zh2 = hi_res_denominators(hires)

    surfaces, failures = {}, []
    for blend, key, default, label in SURFACES:
        pct = getattr(specs, blend)
        a1, a2, zl1, zl2 = low_res_weighting(wvl, bb3, pct)
        model = (a1 * 0.01 / zh1, a2 * 0.01 / zh2)      # what radmod computes
        consistent = (a1 * 0.01 / zl1, a2 * 0.01 / zl2)  # one grid throughout
        want = RECORDED[key]
        for i, (got, expect) in enumerate(zip(model, want)):
            if abs(got - expect) > REPRODUCTION_TOLERANCE:
                failures.append(f"  {key} band {i + 1}: {got:.12f} against the "
                                f"recorded {expect:.12f}")
        independent = stellar.band_reflectances(
            wvl, np.clip(pct / 100.0, 0.0, 1.0), name=args.spectrum)
        surfaces[key] = {
            "blend": blend,
            "surface": label,
            "declared_default": default,
            "model_bands": [round(model[0], 6), round(model[1], 6)],
            "recorded_bands": [round(want[0], 6), round(want[1], 6)],
            "single_grid_bands": [round(consistent[0], 6), round(consistent[1], 6)],
            "band1_error": round(model[0] - consistent[0], 6),
            "band2_error": round(model[1] - consistent[1], 6),
            "band1_error_relative": round(model[0] / consistent[0] - 1.0, 6),
            "hi_res_reference_bands": [round(independent["band1"], 6),
                                       round(independent["band2"], 6)],
        }

    if failures:
        raise SystemExit(
            "this script no longer reproduces what the model printed, so the "
            "comparison below measures nothing. radmod's integration has "
            "changed, or the blends or the spectrum have:\n" + "\n".join(failures))

    band1_fraction, _ = stellar.band_fractions(name=args.spectrum)
    broadband_error = {
        k: round(band1_fraction * v["band1_error"]
                 + (1 - band1_fraction) * v["band2_error"], 6)
        for k, v in surfaces.items()}

    report = {
        "generated": datetime.date.today().isoformat(),
        "generator": "analysis/ice_albedo.py",
        "finding": "radmod.f90 normalises a numerator integrated on the "
                   "965-point blend grid by a denominator integrated on the "
                   "2048-point spectrum grid. The two do not span the same "
                   "band-1 interval, so the quotient is not a reflectance.",
        "spectrum_lowres": rel(lowres),
        "spectrum_hires": rel(hires),
        "spectrum_sha256": hashlib.sha256(hires.read_bytes()).hexdigest(),
        "library": rel(SURFACESPECS),
        "library_sha256": hashlib.sha256(SURFACESPECS.read_bytes()).hexdigest(),
        "reproduction_tolerance": REPRODUCTION_TOLERANCE,
        "band1_flux_fraction": round(band1_fraction, 6),
        "surfaces": surfaces,
        "broadband_error": broadband_error,
        "note": "band1_error is what the mixed normalisation costs; it is "
                "one-signed DARKER on every surface, so the simulation's "
                "ice-albedo feedback is slightly stronger than a consistent "
                "integration gives, not weaker. The flat declarations at "
                "plasimmod.f90:428-432 never reach the radiation and are not "
                "the defect.",
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")

    w = max(len(k) for k in surfaces)
    print(f"{'constant':<{w}}  {'model b1':>9}  {'1-grid b1':>9}  {'error':>8}  "
          f"{'rel':>7}  {'broadband':>9}")
    for k, s in surfaces.items():
        print(f"{k:<{w}}  {s['model_bands'][0]:9.6f}  "
              f"{s['single_grid_bands'][0]:9.6f}  {s['band1_error']:+8.5f}  "
              f"{s['band1_error_relative']:+7.2%}  {broadband_error[k]:+9.5f}")
    print(f"\nreproduced all {len(RECORDED)} recorded band pairs to within "
          f"{REPRODUCTION_TOLERANCE}")
    print(f"band 2 error is at most "
          f"{max(abs(s['band2_error']) for s in surfaces.values()):.2e}, "
          "because this star has little flux beyond the blend grid's 14.01 um")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
