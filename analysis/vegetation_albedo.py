#!/usr/bin/env python3
"""Star-weighted vegetated land albedo, from measured leaf spectra.

    python analysis/vegetation_albedo.py

This is the derivation behind `model.vegetation_albedo` in
`config/planet.yaml`, which `build_surface_albedo.py` paints over every land
region outside `barren_rock_classes` in `--mode vegetated` -- about four
fifths of this planet's land, and the largest single item in the albedo half
of `analysis/error_budget.json`. Until 2026-08-19 that value was a bare 0.15
on a CLI default: an Earth canopy figure integrated against the Sun, never
re-weighted to this star, while every rock class beside it went through
`analysis/rock_albedo.py` and `playa_clastic` through `analysis/playa_albedo.py`.
`notes/audits/inherited-earth-constants.md` finding 1 is the finding; this
script is the derivation made re-runnable, and it writes
`analysis/vegetation_albedo.json`.

## Why the sign is opposite to rock

An albedo is a reflectance integrated against the light falling on it. A leaf
is dark in the visible and bright in the near infrared -- the opposite slope
to a mineral with near-infrared absorption features -- so a redder star makes
a canopy BRIGHTER where it makes such a mineral darker. `k25v` puts 0.382 of
its shortwave below 0.75 um where the Sun puts 0.517.

## What transfers, and what does not

Leaf reflectance is not canopy albedo: a canopy traps light between elements
and shows ground through gaps, both of which pull the level down. What
transfers is the RATIO of the star-weighted to the Sun-weighted integral, by
the same argument `playa_albedo.py` uses to cancel its preparation offset --
numerator and denominator share whatever the measurement does not represent.
A canopy treatment would give a smaller ratio than a leaf one, so the leaf
ratio UPPER-bounds the correction, and the bracket below runs from no
correction at all to the leaf ratio with the spectral tails carried.

## The band pair

Codes 174, 175 and 176 historically carried one identical field, which was a
fair statement about a rock table and asserts of a canopy the one thing a
canopy certainly does not do: reflect equally either side of 0.75 um. The
same integrals split at the model's band edge give the vegetated band pair,
anchored so the star-flux-weighted combination reproduces the broadband
value by construction.

Data: ECOSTRESS spectral library, `references/ecospeclib-all/`, the green
vegetation VSWIR set (`vegetation.*` with `vswir` in the name, which excludes
the non-photosynthetic material). Directional-hemispherical reflectance, the
geometry a model albedo wants; laboratory leaves, hence the ratio-only use.
"""

from __future__ import annotations

import argparse
import datetime
import glob
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ECOSTRESS = ROOT / "references" / "ecospeclib-all"
SPECTRUM = ROOT / "exoplasim" / "inputs" / "stellarspectra" / "k25v_hr.dat"
OUTPUT = ROOT / "analysis" / "vegetation_albedo.json"

H, C, KB = 6.626e-34, 2.998e8, 1.381e-23
SUN_TEFF = 5772.0
BAND_SPLIT_UM = 0.75            # radmod's two shortwave bands meet here
TRUNCATED = (0.35, 2.5)         # the range every spectrum actually measures
TAILS = (0.2, 4.0)              # padded range; tail reflectance for a leaf
TAIL_REFLECTANCE = 0.05         # right order in both tails: dark in UV and SWIR
EARTH_SUN_ENDMEMBER = 0.15      # the value this derivation replaces

# The audit's own numbers, reproduced as a check that can fail: if this script
# stops matching notes/audits/inherited-earth-constants.md finding 1, one of
# the two is wrong and the disagreement must not be averaged over.
AUDIT = {"sun_truncated": 0.2492, "k25v_truncated": 0.2742,
         "ratio_truncated": 1.101, "ratio_tails": 1.141}


def planck(microns: np.ndarray, teff: float) -> np.ndarray:
    lam = microns * 1e-6
    return (2 * H * C ** 2 / lam ** 5) / (np.exp(H * C / (lam * KB * teff)) - 1)


def read_spectrum(path: Path):
    """Ascending (wavelength um, reflectance fraction) from an ECOSTRESS file."""
    import re
    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        if re.match(r"^\s*[\d.]+\s+[-\d.]+\s*$", line):
            a, b = line.split()[:2]
            rows.append((float(a), float(b)))
    if not rows:
        return None
    arr = np.array(rows)
    arr = arr[arr[:, 0].argsort()]
    return arr[:, 0], np.clip(arr[:, 1] / 100.0, 0.0, 1.0)


def band_mean(wl, refl, star_wl, star_flux, lo, hi, pad=None):
    """Reflectance integrated against a stellar spectrum over [lo, hi].

    With pad set, the reflectance is held at that value outside the measured
    range instead of the range being truncated to the overlap.
    """
    grid = star_wl[(star_wl >= lo) & (star_wl <= hi)]
    if pad is None:
        grid = grid[(grid >= wl.min()) & (grid <= wl.max())]
        if len(grid) < 10:
            return None
        r = np.interp(grid, wl, refl)
    else:
        r = np.where((grid >= wl.min()) & (grid <= wl.max()),
                     np.interp(grid, wl, refl), pad)
    f = np.interp(grid, star_wl, star_flux)
    return float(np.trapezoid(r * f, grid) / np.trapezoid(f, grid))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--spectrum", type=Path, default=SPECTRUM)
    ap.add_argument("--library", type=Path, default=ECOSTRESS)
    args = ap.parse_args()
    if not args.library.is_dir():
        raise SystemExit(f"{args.library} is absent. It is bulk reference data, "
                         "excluded from git; unzip references/ecospeclib_*.zip "
                         "beside this path.")

    data = np.loadtxt(args.spectrum, skiprows=1)
    star_wl, star_k = data[:, 0], data[:, 1]
    star_sun = planck(star_wl, SUN_TEFF)

    files = sorted(p for p in glob.glob(str(args.library / "vegetation.*spectrum.txt"))
                   if "vswir" in Path(p).name)
    per_class: dict[str, list] = {}
    sun_t, k_t, sun_p, k_p, k_b1, k_b2 = [], [], [], [], [], []
    for path in files:
        parsed = read_spectrum(Path(path))
        if parsed is None:
            continue
        wl, refl = parsed
        st = band_mean(wl, refl, star_wl, star_sun, *TRUNCATED)
        kt = band_mean(wl, refl, star_wl, star_k, *TRUNCATED)
        sp = band_mean(wl, refl, star_wl, star_sun, *TAILS, pad=TAIL_REFLECTANCE)
        kp = band_mean(wl, refl, star_wl, star_k, *TAILS, pad=TAIL_REFLECTANCE)
        b1 = band_mean(wl, refl, star_wl, star_k, TRUNCATED[0], BAND_SPLIT_UM)
        b2 = band_mean(wl, refl, star_wl, star_k, BAND_SPLIT_UM, TRUNCATED[1])
        if None in (st, kt, sp, kp, b1, b2):
            continue
        sun_t.append(st); k_t.append(kt); sun_p.append(sp); k_p.append(kp)
        k_b1.append(b1); k_b2.append(b2)
        per_class.setdefault(Path(path).name.split(".")[1], []).append((st, kt))

    n = len(sun_t)
    if n < 500:
        raise SystemExit(f"only {n} spectra parsed; the audit used 553, so the "
                         "library or the filter has changed and the comparison "
                         "below would not mean anything")
    means = {k: float(np.mean(v)) for k, v in
             (("sun_truncated", sun_t), ("k25v_truncated", k_t),
              ("sun_tails", sun_p), ("k25v_tails", k_p),
              ("k25v_band1", k_b1), ("k25v_band2", k_b2))}
    ratio_t = means["k25v_truncated"] / means["sun_truncated"]
    ratio_p = means["k25v_tails"] / means["sun_tails"]

    # The check against the audit. 0.005 of ratio is far below anything the
    # numbers are used for and far above interpolation noise.
    for key, got in (("ratio_truncated", ratio_t), ("ratio_tails", ratio_p)):
        if abs(got - AUDIT[key]) > 0.005:
            raise SystemExit(f"{key} = {got:.4f} against the audit's "
                             f"{AUDIT[key]}; reconcile before using either")

    value = round(EARTH_SUN_ENDMEMBER * ratio_t, 3)
    bracket = [EARTH_SUN_ENDMEMBER, round(EARTH_SUN_ENDMEMBER * ratio_p, 3)]
    shape1 = means["k25v_band1"] / means["k25v_truncated"]
    shape2 = means["k25v_band2"] / means["k25v_truncated"]
    bands = [round(value * shape1, 3), round(value * shape2, 3)]

    # Reconstruction with the model's own band fractions, computed from the
    # same spectrum file over its full range the way lib/stellar.py does. The
    # residual is the range mismatch (full file against 0.35-2.5 um here) and
    # is reported so nobody discovers it as a surprise.
    f1 = float(np.trapezoid(star_k[star_wl <= BAND_SPLIT_UM],
                            star_wl[star_wl <= BAND_SPLIT_UM])
               / np.trapezoid(star_k, star_wl))
    recombined = f1 * bands[0] + (1 - f1) * bands[1]

    report = {
        "generated": datetime.date.today().isoformat(),
        "generator": "analysis/vegetation_albedo.py",
        "spectrum": str(args.spectrum.relative_to(ROOT)),
        "spectrum_sha256": hashlib.sha256(args.spectrum.read_bytes()).hexdigest(),
        "library": str(args.library.relative_to(ROOT)),
        "spectra_used": n,
        "per_class": {k: {"n": len(v),
                          "sun": round(float(np.mean([a for a, _ in v])), 4),
                          "k25v": round(float(np.mean([b for _, b in v])), 4)}
                      for k, v in sorted(per_class.items())},
        "leaf_means": {k: round(v, 4) for k, v in means.items()},
        "ratio_truncated_035_25um": round(ratio_t, 4),
        "ratio_tails_02_40um": round(ratio_p, 4),
        "earth_sun_endmember": EARTH_SUN_ENDMEMBER,
        "vegetation_albedo": value,
        "vegetation_albedo_bracket": bracket,
        "vegetation_albedo_bands": bands,
        "band_shape_from_leaf": [round(shape1, 4), round(shape2, 4)],
        "band1_flux_fraction_full_spectrum": round(f1, 4),
        "band_recombination_residual": round(recombined - value, 4),
        "note": ("The ratio, not the level, transfers from leaf to canopy; the "
                 "leaf ratio upper-bounds the correction, so the bracket runs "
                 "from the uncorrected Earth-Sun endmember to the tails-carried "
                 "leaf ratio. The band pair is the same integrals split at "
                 "0.75 um, anchored to the broadband value."),
    }
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"{n} green-vegetation VSWIR spectra against {args.spectrum.name}")
    for cls, row in report["per_class"].items():
        print(f"  {cls:6s} n={row['n']:4d}  sun {row['sun']:.4f}  k25v {row['k25v']:.4f}")
    print(f"  ratio {ratio_t:.4f} truncated, {ratio_p:.4f} with tails "
          f"(audit: {AUDIT['ratio_truncated']}, {AUDIT['ratio_tails']})")
    print(f"  vegetation_albedo {value} bracket {bracket} bands {bands}")
    print(f"  band recombination residual {recombined - value:+.4f} "
          f"at f1 = {f1:.4f}")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
