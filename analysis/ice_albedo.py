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

## The defect this was written to measure, and its repair

`radmod` used to normalise a numerator integrated on the LOW-resolution grid by
a denominator integrated on the HIGH-resolution one:

- `a1`, `a2` sum `bb3(k)*blend(k)` over `wavelengths`, the 965-point grid the
  surface blends live on, which starts at 0.34 um and ends at 14.01 um.
- `z1`, `z2` sum `bb1`, `bb2` over `wv1`, `wv2`, the 2048-point grid, band 1
  running from the file's start with everything below `minwavel` zeroed, band 2
  running to 100 um, and the band-edge interval added to `z1`.
- `zdenom1 = 0.01/z1` then divided the one by the other.

A ratio of integrals is a reflectance only when both cover the same interval
with the same spectrum. These did not, so the quotient implicitly assigned ZERO
reflectance to the flux outside the blend grid: 1.48 percent of band 1 for this
star and 0.06 percent of band 2, one-signed DARK on every surface. It was
CLAUDE.md rule 3's class of defect -- two labellings of one axis combined
without reconciliation -- on wavelength rather than longitude.

The repair extends the numerator rather than shrinking the denominator, because
the albedo multiplies flux over the whole band and so must be averaged over the
whole band. Each blend is held at its endpoint value across the gap. The
alternative -- normalising by the blend grid's own flux and accepting an average
over a sub-interval -- is reported below as `single_grid_bands`; it differs by
at most 0.001, on ground, where an endpoint value is least representative.

## The check that can fail

`RECORDED` holds what the model itself printed under "Finalized Albedos". This
script reproduces those numbers from source and refuses if it cannot, which is
what keeps it a verification rather than a second opinion: the reproduction is
exact, so anything that moves the model's integration breaks it here first.

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
    # "Finalized Albedos" from a two-step T21 integration on the k25v spectrum
    # with the repaired radmod, band 1 then band 2. The pre-repair values were
    # 0.968193/0.585286, 0.500267/0.272589, 0.878555/0.472884, 0.626302/0.339473
    # and 0.752429/0.406179; the difference is what the mixed grids cost.
    "dsnowalb": (0.76382640121202905, 0.40611515565148087),
    "dsnowalbmx": (0.98278401164203333, 0.58529460781319775),
    "dsnowalbmn": (0.50782484289747332, 0.27261180703024718),
    "dicealbmx": (0.89180978408407729, 0.47289790510791729),
    "dicealbmn": (0.63575483791955212, 0.33949058411682098),
    "dglacalbmn": (0.76378241785467849, 0.40619422519279852),
    "dgroundalb": (0.19354492153968911, 0.23002255831102814),
    "doceanalb": (0.075569967429263696, 0.063668098143514609),
}

SURFACES = [
    ("iceblend", "dsnowalb", 0.75, "snow, base"),
    ("iceblendmax", "dsnowalbmx", 0.8, "snow, maximum"),
    ("iceblendmin", "dsnowalbmn", 0.4, "snow, minimum"),
    ("seaicemax", "dicealbmx", 0.7, "sea ice, maximum"),
    ("seaicemin", "dicealbmn", 0.5, "sea ice, minimum"),
    ("glacalbmin", "dglacalbmn", 0.6, "glacier ice, minimum"),
    ("groundblend", "dgroundalb", 0.2, "bare ground"),
    ("oceanblend", "doceanalb", 0.069, "open ocean"),
]

# All eight surfaces radmod weights in this block. The list is complete on
# purpose: the same zdenom pair normalises every one of them, so a check that
# covered only the ice surfaces would have left the same defect live on ground
# and ocean.


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


def hi_res_terms(path: Path, edge_lo_um: float, edge_hi_um: float):
    """`z1`, `z2`, `zout1`, `zout2` exactly as `radmod.f90:519-565` builds them.

    `zout1` and `zout2` are the band flux lying outside the blend grid, summed
    over whole intervals as the Fortran loop does. They are what the repair adds
    to each numerator so the quotient covers one interval.
    """
    wavelength, flux = stellar.read_hires(path)
    wv1, bb1, wv2, bb2 = stellar._bands(wavelength, flux)
    z1 = float(np.trapezoid(bb1, wv1))
    z2 = float(np.trapezoid(bb2, wv2))
    z1 += 0.5 * (bb1[-1] + bb2[0]) * (wv2[0] - wv1[-1])
    lo, hi = edge_lo_um * 1.0e-6, edge_hi_um * 1.0e-6
    zout1 = zout2 = 0.0
    for k in range(1, 1024):
        if wv1[k] <= lo:
            zout1 += 0.5 * (bb1[k - 1] + bb1[k]) * (wv1[k] - wv1[k - 1])
        if wv2[k - 1] >= hi:
            zout2 += 0.5 * (bb2[k - 1] + bb2[k]) * (wv2[k] - wv2[k - 1])
    return z1, z2, zout1, zout2


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
    zh1, zh2, zout1, zout2 = hi_res_terms(hires, wvl[0], wvl[-1])

    surfaces, failures = {}, []
    for blend, key, default, label in SURFACES:
        pct = getattr(specs, blend)
        a1, a2, zl1, zl2 = low_res_weighting(wvl, bb3, pct)
        # The repaired radmod: numerator extended across the flux the blend grid
        # does not reach, holding the blend at its endpoint value.
        model = ((a1 + pct[0] * zout1) * 0.01 / zh1,
                 (a2 + pct[-1] * zout2) * 0.01 / zh2)
        # The alternative repair, kept so the choice stays priced: normalise by
        # the blend grid's own flux and accept an average over a sub-interval.
        consistent = (a1 * 0.01 / zl1, a2 * 0.01 / zl2)
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
            "band1_vs_alternative": round(model[0] - consistent[0], 6),
            "band2_vs_alternative": round(model[1] - consistent[1], 6),
            "hi_res_reference_bands": [round(independent["band1"], 6),
                                       round(independent["band2"], 6)],
        }

    if failures:
        raise SystemExit(
            "this script no longer reproduces what the model printed, so the "
            "comparison below measures nothing. radmod's integration has "
            "changed, or the blends or the spectrum have:\n" + "\n".join(failures))

    band1_fraction, _ = stellar.band_fractions(name=args.spectrum)
    broadband_difference = {
        k: round(band1_fraction * v["band1_vs_alternative"]
                 + (1 - band1_fraction) * v["band2_vs_alternative"], 6)
        for k, v in surfaces.items()}

    report = {
        "generated": datetime.date.today().isoformat(),
        "generator": "analysis/ice_albedo.py",
        "finding": "radmod.f90's surface albedo quotient covers one interval "
                   "in numerator and denominator. Before the repair the "
                   "numerator stopped at the blend grid's 0.34 um while the "
                   "denominator ran from minwavel, which assigned zero "
                   "reflectance to 1.48 percent of band 1.",
        "spectrum_lowres": rel(lowres),
        "spectrum_hires": rel(hires),
        "spectrum_sha256": hashlib.sha256(hires.read_bytes()).hexdigest(),
        "library": rel(SURFACESPECS),
        "library_sha256": hashlib.sha256(SURFACESPECS.read_bytes()).hexdigest(),
        "reproduction_tolerance": REPRODUCTION_TOLERANCE,
        "band1_flux_fraction": round(band1_fraction, 6),
        "surfaces": surfaces,
        "band1_flux_outside_blend_grid": round(zout1 / zh1, 6),
        "band2_flux_outside_blend_grid": round(zout2 / zh2, 6),
        "broadband_vs_alternative": broadband_difference,
        "note": "band1_vs_alternative is the residual choice between the two "
                "self-consistent repairs, not an error: extending the "
                "numerator against normalising by the blend grid's own flux. "
                "It is largest on ground, where holding an endpoint value "
                "across the gap is least representative. The flat "
                "declarations at plasimmod.f90:428-432 never reach the "
                "radiation and were never the defect.",
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")

    w = max(len(k) for k in surfaces)
    print(f"{'constant':<{w}}  {'band 1':>9}  {'band 2':>9}  {'alt b1':>9}  "
          f"{'vs alt':>8}")
    for k, v in surfaces.items():
        print(f"{k:<{w}}  {v['model_bands'][0]:9.6f}  {v['model_bands'][1]:9.6f}  "
              f"{v['single_grid_bands'][0]:9.6f}  "
              f"{v['band1_vs_alternative']:+8.5f}")
    print(f"\nreproduced all {len(RECORDED)} recorded band pairs to within "
          f"{REPRODUCTION_TOLERANCE}")
    print(f"flux outside the blend grid: band 1 {zout1 / zh1:.4%}, "
          f"band 2 {zout2 / zh2:.4%}")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
