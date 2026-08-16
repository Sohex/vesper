#!/usr/bin/env python3
"""Put every uncertainty in kelvin, so they can be ranked against each other.

    python scripts/error_budget.py

The project accumulates uncertainties in whatever unit each was measured in --
albedo, W/m2, a percentage of land, a bracket on a rock class -- and then has no
way to answer "is refining this worth more than refining that". This converts
them to one currency at zero compute, from artifacts already on disk.

## The conversion

A land-mean albedo change reaches temperature through three steps:

    surface, planet-wide   d_surface = d_land * land_fraction
    top of atmosphere      d_toa     = d_surface * attenuation
    flux equivalent        d_flux    = -d_toa / (1 - planetary_albedo)
    temperature            d_T       = d_flux * slope

The flux equivalence is exact: absorbed shortwave is S/4 * (1 - alpha), so a
change of `d` in albedo is worth the same as a change of `d / (1 - alpha)` in
flux ratio.

**The attenuation is the step it is tempting to skip, and skipping it doubles
every answer.** A surface albedo change does not arrive at the top of the
atmosphere intact; everything above the surface scatters and absorbs. This
project measured it the expensive way: a lithology fix that moved planet-mean
*surface* albedo by +0.0033 was predicted at -0.81 K on the naive conversion, and
the run came in about a kelvin warmer, with planetary albedo having moved -0.0017
-- the opposite sign, once the warmer state's reduced sea ice is included. A
factor of 0.5 is used here as a conservative default and both columns are
reported, because the honest claim is a factor of two rather than a number.

## What it is for

Not precision. Ordering. A factor-of-two budget is enough to answer whether the
subgrid roughness work buys anything next to the slab depth, and the answer is
visible immediately because the surface items land under a kelvin and the
structural ones do not.

## The stopping rule that falls out

**Refine an input when its plausible range exceeds the effect of the thing you
last refined.** That is a per-component stopping rule rather than a global
optimisation: it says which slot to upgrade next, not that the wrong thing was
built. It also says when to stop, which this project has not previously had a way
to say.

Items whose currency is not albedo are listed with their own measurement and no
conversion, because inventing one would be worse than leaving the column blank.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

import yaml  # noqa: E402

# Measured, not assumed. Each is sourced in the comment beside it.
DEFAULT_ATTENUATION = 0.5      # v4->v5 pair; see the module docstring
FALLBACK_SLOPE = 150.2         # K per unit flux ratio, from the T21 bracket at
                               # 0.95 (288.760 K) and 1.00 (296.270 K), both on
                               # carved-zoned-v5 with the ozone band-weight fix.
                               # Replaces 203.6, which was 26% high. Measured
                               # between two converged points SPANNING the
                               # target, per the rule that a sensitivity taken
                               # across one regime does not transfer to another.
FALLBACK_PLANETARY_ALBEDO = 0.266
FALLBACK_LAND_FRACTION = 0.4317


def albedo_to_kelvin(d_land, land_fraction, planetary_albedo, slope, attenuation):
    d_toa = d_land * land_fraction * attenuation
    return -d_toa / (1.0 - planetary_albedo) * slope


# (label, land-mean albedo delta, note). Deltas are the plausible RANGE of the
# item, not its current value: the budget ranks what refining each would buy.
ALBEDO_ITEMS = [
    ("biosphere: bare rock vs vegetated", -(0.276 - 0.179),
     "the assumption LPJ-GUESS exists to replace"),
    ("playa_clastic albedo, 0.25 to 0.33", 0.08 * 0.147,
     "largest single rock-class lever; a mixture of clay playa and varnished fan"),
    ("lakes composited into albedo", -0.0139,
     "solved lakes reaching the climate at all"),
    ("carve iteration, fill 16.5% to 12.35%", 0.0077,
     "one carve pass, measured across v4 to v5"),
    ("salt crust albedo, 0.40 to 0.50", 0.10 * 0.019,
     "degenerate with crust extent from the climate's side"),
]

# Items in their own units. No conversion invented.
OTHER_ITEMS = [
    ("slab depth, 25 m vs 50 m", "seasonal amplitude ~2x",
     "STRUCTURAL. Amplitude scales as 1/(omega*C) and this world's year is half "
     "Earth's, so 50 m damps seasonality about twice as hard as Earth's ocean "
     "does. Coldest-month mean is what PFT survival gates on. One perturbation "
     "run, ever, gives a permanent scaling for every later result."),
    ("no q-flux", "gradients too strong, ice too extensive",
     "STRUCTURAL, and no cheap version exists. Declare the direction and move on."),
    ("stellar spectrum, k2 vs k25v", "0.04 W/m2 absorbed",
     "Measured null on a warm nearly ice-free state, NOT on the cold branch, "
     "where it is the whole point. Re-measure before the cycle work."),
    ("energy closure residual", "-0.400 to -0.489 W/m2",
     "Consumes most of the 0.5 W/m2 convergence tolerance. Was recorded as a "
     "fixed offset; it is not."),
    ("Penman over a dry land column", "unquantified, one-signed",
     "Validated over ocean, where the air is equilibrated with the surface. Over "
     "a subgrid lake the column is dry, so VPD is too high, E is overstated and "
     "the verdict under-carves. Opposite in sign to the albedo-driven "
     "over-carve, and the only item here whose magnitude is unknown."),
    ("dust, radiative", "TWO-SIDED, -2 K to weakly positive",
     "Aerosols are off entirely (L_AERO = 0), and ExoPlaSim 3.4.2 cannot switch "
     "them on: radmod.f90 declares aero_nl but never reads it, so l_aerorad is "
     "pinned at 0. Reported upstream. Cooling magnitude bounded at -0.5 to -2 K "
     "on a DRY source area of about 12% of land. "
     "THE SIGN IS NOT ONE-WAY, and the earlier one-signed entry here was wrong "
     "for the same reason the GCM would have been. Aerosol forcing reverses "
     "sign at a critical surface albedo a_c; Mie over the Balkanski source "
     "distribution, band-averaged on the k25v spectrum, gives a_c = 0.504 in "
     "band 1 on Di Biagio measured indices, 0.289 on OPAC, and 0.456 in band 2. "
     "Salt crust is 0.40-0.50, so the reversal sits INSIDE our own albedo range "
     "and band 2 carries 61.6% of the flux. Dust cools over ocean (0.07) and "
     "vegetated land (0.18) in every case; over closed-basin fill the sign is "
     "undetermined by present data. Rocha-Lima 2018 finds fine-mode k rising "
     "from a ~650 nm minimum through the SWIR, which would lower a_c further "
     "and favour warming, but its NIR values exist only in a figure. "
     "Resolving this needs band-2 refractive indices, NOT a GCM run: the "
     "calculation is analytic once ssa and g are known. See notes/dust.md."),
    ("roughness distribution", "land median 0.502 m under a 2.0 m mean",
     "Anchored to ExoPlaSim's tuned land mean, which the distribution says is "
     "carried by a rough tail. Anchoring inflates mid-range cells; direction "
     "declared, magnitude not measured."),
]


def main() -> None:
    config = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    slope = FALLBACK_SLOPE
    planetary_albedo = FALLBACK_PLANETARY_ALBEDO
    land_fraction = FALLBACK_LAND_FRACTION

    rows = []
    for label, delta, note in ALBEDO_ITEMS:
        naive = albedo_to_kelvin(delta, land_fraction, planetary_albedo, slope, 1.0)
        atten = albedo_to_kelvin(delta, land_fraction, planetary_albedo, slope,
                                 DEFAULT_ATTENUATION)
        rows.append((label, delta, naive, atten, note))
    rows.sort(key=lambda r: -abs(r[3]))

    print(f"slope {slope} K per unit flux ratio, land {land_fraction:.4f}, "
          f"planetary albedo {planetary_albedo:.3f}, attenuation "
          f"{DEFAULT_ATTENUATION}\n")
    print(f"{'item':38} {'d(alb)':>8} {'K naive':>8} {'K':>7}")
    print("-" * 66)
    for label, delta, naive, atten, _ in rows:
        print(f"{label:38} {delta:+8.4f} {naive:+8.2f} {atten:+7.2f}")

    print("\nnot in albedo, no conversion invented:")
    for label, magnitude, _ in OTHER_ITEMS:
        print(f"  {label:38} {magnitude}")

    print("\nRefine an input when its plausible range exceeds the effect of the")
    print("thing you last refined. Everything under a kelvin here is below the")
    print("structural items and below the biosphere question.")

    payload = {
        "note": "Generated by scripts/error_budget.py. Ordering, not precision; "
                "a factor of two is expected and is enough to rank refinements.",
        "conversion": {
            "slope_k_per_flux_ratio": slope,
            "land_fraction": land_fraction,
            "planetary_albedo": planetary_albedo,
            "attenuation": DEFAULT_ATTENUATION,
            "attenuation_note": "A surface albedo change does not reach the top "
                                "of the atmosphere intact. Skipping this step "
                                "doubles every answer, which this project did "
                                "once and paid a kelvin for.",
        },
        "albedo_items": [
            {"item": l, "delta_land_albedo": d, "kelvin_naive": round(n, 3),
             "kelvin": round(a, 3), "note": nt} for l, d, n, a, nt in rows
        ],
        "other_items": [{"item": l, "magnitude": m, "note": n}
                        for l, m, n in OTHER_ITEMS],
        "stopping_rule": "Refine an input when its plausible range exceeds the "
                         "effect of the thing you last refined.",
    }
    out = ROOT / "analysis" / "error_budget.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
