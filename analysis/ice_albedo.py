#!/usr/bin/env python3
"""Star-weighted snow, glacier ice and sea ice albedo, in the model's two bands.

    python analysis/ice_albedo.py

This is the derivation behind `model.snow_albedo_*`, `model.glacier_albedo_min`
and `model.sea_ice_albedo_*` in `config/planet.yaml`, and it writes
`analysis/ice_albedo.json`.

## Why these surfaces and not the others

`plasimmod.f90:428-432` declares `dsnowalbmx(2)`, `dsnowalbmn(2)` and
`dglacalbmn(2)` as two-element BAND arrays initialised from a scalar, so
element 1 equals element 2 and the simulation's snow reflects the same fraction
either side of 0.75 um. `dicealbmx` and `dicealbmn` in `icemod` are the same
and have no configuration key at all. That is the defect the vegetation work
already fixed for the canopy, where `config/planet.yaml` records "one identical
field used to assert that a canopy reflects equally either side of the split".

It also means the stellar spectrum reaches these surfaces not at all. The
two-band scheme carries a star's shape by moving flux BETWEEN the bands; moving
flux between two equal numbers changes nothing.

## Where the numbers come from

`vendor/exoplasim/exoplasim/surfacespecs.py` ships the ECOSTRESS snow, ice,
water and frost blends the model's own constants were made from, interpolated
to the 965 wavelengths ExoPlaSim uses, and only `pRT.py` reads them. Nothing
weights them against the configured star. This does.

## The check that can fail

Weighting each blend by a 5772 K Planck must reproduce the constant it was
built to yield. Six for six is what establishes that these spectra ARE the
provenance of the model's constants, and it validates the integral before
anything is concluded from the K dwarf column. `SOLAR_TOLERANCE` says what the
check is worth and why it cannot be tighter.
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
SOURCE = ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "src"
OUTPUT = ROOT / "analysis" / "ice_albedo.json"

SOLAR_TOLERANCE = 0.02
"""How far the 5772 K column may sit from the constant the blend was built for.

Set from what the comparison is worth rather than from the size of any error it
is meant to catch, and it cannot be tight: the model states these constants to
ONE DECIMAL PLACE, so a blend reproducing 0.7 exactly would still be read back
as anything in 0.65 to 0.75. Rounding therefore dominates, and 0.02 is that
floor with the integration's own spread on top -- two independent integrations
of the same quantity, one on `solarini`'s grids and one on the blend's own,
differ by about 0.004.

What it can still catch is everything it exists for: a substituted or
re-interpolated blend, a band convention that splits somewhere other than where
`solarini` splits, a spectrum file at the wrong resolution, and a key paired
with the wrong blend. It cannot certify the constants to better than their own
stated precision, and it does not claim to.
"""

# Blend, the model constant it was built to yield, and where that constant is
# declared. `dsnowalb` and `doceanalb` are excluded deliberately: the first is
# a base array the min/max pair supersedes in use, and the second is a
# declaration `radmod.f90:2892-2903` overwrites with a zenith formula over open
# water, leaving it in force only under sea ice. OCN-7 and OCN-22 own that.
SURFACES = [
    ("iceblendmax", "snow_albedo_max",     0.8, "plasimmod.f90 dsnowalbmx"),
    ("iceblendmin", "snow_albedo_min",     0.4, "plasimmod.f90 dsnowalbmn"),
    ("glacalbmin",  "glacier_albedo_min",  0.6, "plasimmod.f90 dglacalbmn"),
    ("seaicemax",   "sea_ice_albedo_max",  0.7, "icemod.f90 dicealbmx"),
    ("seaicemin",   "sea_ice_albedo_min",  0.5, "icemod.f90 dicealbmn"),
]

REFERENCE = ("oceanblend", 0.069, "seamod.f90 doceanalb")
"""Carried as a NEGATIVE result rather than a key.

Water is dark and nearly flat across the split, so re-weighting it for this
host moves the broadband value by about 0.002 and the band pair barely
separates. It does not earn a key, and recording that here is what stops the
question being reopened. `notes/external-model-survey.md` section 11c.
"""


def load_surfacespecs():
    """Import the spectra without importing the ExoPlaSim package around them."""
    spec = importlib.util.spec_from_file_location("surfacespecs", SURFACESPECS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def evaluate(wvl, percent, spectrum_name):
    """Sun and star band integrals for one blend, reflectance given in percent."""
    reflectance = np.clip(percent / 100.0, 0.0, 1.0)
    sun = stellar.band_reflectances(wvl, reflectance, temperature_k=5772.0)
    star = stellar.band_reflectances(wvl, reflectance, name=spectrum_name)
    return sun, star


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--spectrum", default=None,
                    help="stellar spectrum name; the configured one by default")
    ap.add_argument("--output", type=Path, default=OUTPUT)
    args = ap.parse_args()

    specs = load_surfacespecs()
    wvl = specs.wvl
    if wvl.size != 965:
        raise SystemExit(
            f"surfacespecs ships {wvl.size} wavelengths, not the 965 its own "
            "docstring names. The blends have been re-interpolated and the "
            "solar check below would be comparing something else.")

    stellar.solar_partition_identity()

    hires = stellar.spectrum_paths(args.spectrum)[1]
    surfaces, failures = {}, []
    for blend, key, constant, declared_at in SURFACES:
        sun, star = evaluate(wvl, getattr(specs, blend), args.spectrum)
        miss = sun["broadband"] - constant
        if abs(miss) > SOLAR_TOLERANCE:
            failures.append(f"  {blend}: 5772 K gives {sun['broadband']:.4f} "
                            f"against {declared_at} = {constant}, off by {miss:+.4f}")
        surfaces[key] = {
            "blend": blend,
            "declared_at": declared_at,
            "model_constant": constant,
            "solar_broadband": round(sun["broadband"], 4),
            "solar_check_residual": round(miss, 4),
            "bands": [round(star["band1"], 3), round(star["band2"], 3)],
            "broadband": round(star["broadband"], 3),
            "shift_from_model_constant": round(star["broadband"] - constant, 3),
            "recombination_residual": float(f"{star['recombination_residual']:.3e}"),
        }

    if failures:
        raise SystemExit(
            "the 5772 K column no longer reproduces the model's own constants, "
            "so these blends are not their provenance and nothing should be "
            "concluded from the star column:\n" + "\n".join(failures))

    ref_blend, ref_constant, ref_declared = REFERENCE
    ref_sun, ref_star = evaluate(wvl, getattr(specs, ref_blend), args.spectrum)

    band1_fraction, _ = stellar.band_fractions(name=args.spectrum)
    report = {
        "generated": datetime.date.today().isoformat(),
        "generator": "analysis/ice_albedo.py",
        "spectrum": rel(hires),
        "spectrum_sha256": hashlib.sha256(hires.read_bytes()).hexdigest(),
        "library": rel(SURFACESPECS),
        "library_sha256": hashlib.sha256(SURFACESPECS.read_bytes()).hexdigest(),
        "wavelengths": int(wvl.size),
        "band_split_um": stellar.BAND_SPLIT_UM,
        "band1_flux_fraction": round(band1_fraction, 4),
        "solar_tolerance": SOLAR_TOLERANCE,
        "flux_outside_measured_range": float(
            f"{ref_star['flux_outside_measured']:.3e}"),
        "surfaces": surfaces,
        "ocean_reference": {
            "blend": ref_blend,
            "declared_at": ref_declared,
            "model_constant": ref_constant,
            "solar_broadband": round(ref_sun["broadband"], 4),
            "bands": [round(ref_star["band1"], 3), round(ref_star["band2"], 3)],
            "broadband": round(ref_star["broadband"], 3),
            "verdict": "no key: water is flat across the split, so the "
                       "re-weighting is worth about 0.002 and the band pair "
                       "barely separates",
        },
        "note": "Bands are the star-flux-weighted integrals either side of the "
                "model's band edge, so their recombination at "
                "band1_flux_fraction reproduces the broadband value by "
                "construction; the residual is the band-edge interval "
                "solarini books to band 1. Every surface here is currently "
                "too BRIGHT in the model, one-signed against the ice-albedo "
                "feedback at the cold end.",
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")

    width = max(len(k) for k in surfaces)
    print(f"{'key':<{width}}  {'5772K':>7}  {'const':>6}  "
          f"{'band1':>6}  {'band2':>6}  {'broad':>6}  {'shift':>7}")
    for key, s in surfaces.items():
        print(f"{key:<{width}}  {s['solar_broadband']:7.4f}  "
              f"{s['model_constant']:6.3f}  {s['bands'][0]:6.3f}  "
              f"{s['bands'][1]:6.3f}  {s['broadband']:6.3f}  "
              f"{s['shift_from_model_constant']:+7.3f}")
    print(f"\nsolar check passed for all {len(surfaces)} surfaces within "
          f"{SOLAR_TOLERANCE}")
    print(f"ocean, for reference and NOT a key: {ref_star['broadband']:.3f} "
          f"against {ref_sun['broadband']:.3f} under the Sun")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
