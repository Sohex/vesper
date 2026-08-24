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

## The band pair, and what it is anchored on

Codes 174, 175 and 176 historically carried one identical field, which was a
fair statement about a rock table and asserts of a canopy the one thing a
canopy certainly does not do: reflect equally either side of 0.75 um. The
same integrals split at the model's band edge give the vegetated band pair.

What transfers from the spectra is again a RATIO and not a level: `band2/band1`,
how much brighter the canopy is beyond 0.75 um than below it. The level is set
by anchoring, and the anchor is the model's own band weights. `radmod.f90` forms
`zsolars(1)*dsalb(1) + zsolars(2)*dsalb(2)`, and `lib/stellar.band_fractions`
reproduces `solarini` exactly -- the `minwavel` cut and the band-edge interval
included -- so

    band1 = broadband / (z1 + z2 * rho),    band2 = rho * band1

makes the recombination an IDENTITY in the weights the radiation actually uses.

Anchoring instead on the flux share of the range the leaves were MEASURED over
does not, because 0.35-2.5 um is not the star's whole shortwave: that share is
0.399 against the model's 0.382, and the pair it produces recombines about
0.0025 too bright on vegetated ground, one-signed. That number is
`band_recombination_residual_naive_anchor` below, it is what this project
shipped until the anchoring was corrected, and it is the residual `bio-18` asks
to see reported.

## Whether tree and grass need separate band ratios: not resolved, and by how much

The broadband correction takes the POPULATION ratio, because tree and grass
differ there by 0.006 in ratio and 0.001 in albedo against a grass sample of
four spectra. The BAND ratio is a much larger difference -- tree 3.15 against
grass 2.56 -- and the question is whether four grass spectra resolve it.

The test was fixed before it was run: the classes are separated if the whole
grass sample lies below the trees' tenth percentile. It misses, by less than
the last digit either number is quoted at. So the endmembers take the
population ratio and carry the per-class ratio as the other end of the bracket,
which is the same answer the broadband correction reached and for the same
reason. `band_ratio_separates_tree_from_grass` in the output is that test, and
what it takes to settle it is more grass spectra, not more argument: the trees'
own ratios run from 1.78 to 4.92, so a class mean here is only ever as good as
its sample.

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
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from lib.paths import rel  # noqa: E402
from lib import stellar  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ECOSTRESS = ROOT / "references" / "ecospeclib-all"
SPECTRUM = ROOT / "exoplasim" / "inputs" / "stellarspectra" / "k25v_hr.dat"
OUTPUT = ROOT / "analysis" / "vegetation_albedo.json"

H, C, KB = 6.626e-34, 2.998e8, 1.381e-23
SUN_TEFF = 5772.0
BAND_SPLIT_UM = stellar.BAND_SPLIT_UM   # radmod's two shortwave bands meet here
TRUNCATED = (0.35, 2.5)         # the range every spectrum actually measures
TAILS = (0.2, 4.0)              # padded range; tail reflectance for a leaf
TAIL_REFLECTANCE = 0.05         # right order in both tails: dark in UV and SWIR
EARTH_SUN_ENDMEMBER = 0.15      # the value this derivation replaces

# The other two endmembers named in the same sentence of
# `notes/audits/inherited-earth-constants.md` finding 1, and left behind when
# 0.15 was re-weighted: `build_surface_albedo.py`'s full-tree-cover and
# full-grass-cover albedos, which are what `--mode modelled` blends. Same kind
# of quantity as 0.15, so they take the same correction.
#
# The POPULATION ratio and not the per-class one. The per-class truncated
# ratios are 1.105 for tree and 1.099 for grass against 1.100 for the whole
# leaf set, a spread of 0.006 in ratio and 0.001 in albedo, and the grass class
# has n = 4. A per-class ratio is not better resolved than the effect it would
# be claiming to separate.
COVER_ENDMEMBERS = {"tree": 0.13, "grass": 0.19}

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
        per_class.setdefault(Path(path).name.split(".")[1], []).append(
            (st, kt, b1, b2))

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

    # The model's own weights, not a second integration of the same file. This
    # is the one place the star's colour enters the surface energy balance and
    # `lib/stellar` is the only module allowed to compute it; a hand-rolled
    # split here misses the minwavel cut and the band-edge interval, which is
    # worth 0.0022 of z1 and was the whole of the residual reported before.
    z1, z2 = (float(v) for v in stellar.band_fractions(args.spectrum))

    def anchored(level: float, rho: float) -> list[float]:
        """The band pair for a broadband `level` at band ratio `rho`.

        Anchored so `z1*band1 + z2*band2 == level` exactly, which is what the
        radiation forms. Rounded to three decimals like every other albedo this
        project writes; the rounding is worth at most 5e-4 and the writer's
        recombination check is set well above it.
        """
        band1 = level / (z1 + z2 * rho)
        return [round(band1, 3), round(rho * band1, 3)]

    rho_all = means["k25v_band2"] / means["k25v_band1"]
    bands = anchored(value, rho_all)

    # Per-class band ratios, and the pre-fixed test of whether the class split
    # is resolved: the grass sample must lie entirely below the trees' tenth
    # percentile. Whichever way it comes out, the two ratios are the two ends
    # of the endmembers' bracket.
    class_rho = {}
    for cls, rows in per_class.items():
        b1 = float(np.mean([r[2] for r in rows]))
        b2 = float(np.mean([r[3] for r in rows]))
        class_rho[cls] = b2 / b1
    tree_p10 = float(np.percentile(
        [r[3] / r[2] for r in per_class["tree"]], 10))
    grass_max = max(r[3] / r[2] for r in per_class["grass"])
    rho_resolved = grass_max < tree_p10

    # The central value takes the ratio the test above licenses, and the other
    # ratio is the bracket. Both endmembers therefore get their own PAIR even
    # when they share a ratio, because their levels differ; that is what makes
    # codes 175 and 176 distinct fields under `--mode modelled`.
    cover_bands = {}
    cover_bands_bracket = {}
    cover_band_ratio = {}
    for cls, level in ((k, round(v * ratio_t, 3))
                       for k, v in COVER_ENDMEMBERS.items()):
        chosen = class_rho[cls] if rho_resolved else rho_all
        other = rho_all if rho_resolved else class_rho[cls]
        cover_band_ratio[cls] = chosen
        cover_bands[cls] = anchored(level, chosen)
        alt = anchored(level, other)
        cover_bands_bracket[cls] = [sorted((cover_bands[cls][0], alt[0])),
                                    sorted((cover_bands[cls][1], alt[1]))]

    # What the pair would have come to under the anchoring this script used
    # before: normalised on the flux share of the MEASURED range rather than on
    # the star's whole shortwave. Reported rather than deleted, because it is
    # the size of the error that anchoring carried and bio-18 asks for it.
    g1 = ((means["k25v_band2"] - means["k25v_truncated"])
          / (means["k25v_band2"] - means["k25v_band1"]))
    naive = [value * means["k25v_band1"] / means["k25v_truncated"],
             value * means["k25v_band2"] / means["k25v_truncated"]]
    recombined = z1 * bands[0] + z2 * bands[1]
    naive_residual = z1 * naive[0] + z2 * naive[1] - value

    report = {
        "generated": datetime.date.today().isoformat(),
        "generator": "analysis/vegetation_albedo.py",
        "spectrum": rel(args.spectrum),
        "spectrum_sha256": hashlib.sha256(args.spectrum.read_bytes()).hexdigest(),
        "library": rel(args.library),
        "spectra_used": n,
        "per_class": {k: {"n": len(v),
                          "sun": round(float(np.mean([r[0] for r in v])), 4),
                          "k25v": round(float(np.mean([r[1] for r in v])), 4),
                          "band2_over_band1": round(class_rho[k], 4)}
                      for k, v in sorted(per_class.items())},
        "leaf_means": {k: round(v, 4) for k, v in means.items()},
        "ratio_truncated_035_25um": round(ratio_t, 4),
        "ratio_tails_02_40um": round(ratio_p, 4),
        "earth_sun_endmember": EARTH_SUN_ENDMEMBER,
        "vegetation_albedo": value,
        "vegetation_albedo_bracket": bracket,
        "cover_albedo": {k: round(v * ratio_t, 3)
                         for k, v in COVER_ENDMEMBERS.items()},
        "cover_albedo_earth_sun": dict(COVER_ENDMEMBERS),
        "cover_albedo_bracket": {k: [v, round(v * ratio_p, 3)]
                                 for k, v in COVER_ENDMEMBERS.items()},
        "vegetation_albedo_bands": bands,
        "cover_albedo_bands": cover_bands,
        "cover_albedo_bands_bracket": cover_bands_bracket,
        "cover_band_ratio": {k: round(v, 4)
                             for k, v in sorted(cover_band_ratio.items())},
        "tree_band_ratio_p10": round(tree_p10, 4),
        "grass_band_ratio_max": round(grass_max, 4),
        "band2_over_band1": round(rho_all, 4),
        "band2_over_band1_per_class": {k: round(v, 4)
                                       for k, v in sorted(class_rho.items())},
        "band_ratio_separates_tree_from_grass": rho_resolved,
        "model_band_flux_fractions": [round(z1, 6), round(z2, 6)],
        "model_band_flux_source": ("lib/stellar.band_fractions, which reproduces "
                                   "radmod.f90 solarini including the minwavel "
                                   "cut and the band-edge interval"),
        "measured_range_band1_flux_share": round(g1, 4),
        "band_recombination_residual": round(recombined - value, 5),
        "band_recombination_residual_naive_anchor": round(naive_residual, 5),
        "note": ("The ratio, not the level, transfers from leaf to canopy; the "
                 "leaf ratio upper-bounds the correction, so the bracket runs "
                 "from the uncorrected Earth-Sun endmember to the tails-carried "
                 "leaf ratio. The band pair is the same integrals split at "
                 "0.75 um, anchored on the model's own band weights so that "
                 "z1*band1 + z2*band2 returns the broadband value as an "
                 "identity. cover_albedo_bands gives tree and grass their own "
                 "band ratios, which unlike the broadband ratio are resolved "
                 "against the within-class spread."),
    }
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"{n} green-vegetation VSWIR spectra against {args.spectrum.name}")
    for cls, row in report["per_class"].items():
        print(f"  {cls:6s} n={row['n']:4d}  sun {row['sun']:.4f}  k25v {row['k25v']:.4f}")
    print(f"  ratio {ratio_t:.4f} truncated, {ratio_p:.4f} with tails "
          f"(audit: {AUDIT['ratio_truncated']}, {AUDIT['ratio_tails']})")
    print(f"  vegetation_albedo {value} bracket {bracket} bands {bands}")
    for cls, v in COVER_ENDMEMBERS.items():
        print(f"  {cls}_albedo {report['cover_albedo'][cls]} "
              f"bracket {report['cover_albedo_bracket'][cls]} "
              f"bands {cover_bands[cls]} at rho {cover_band_ratio[cls]:.3f}")
    print(f"  band ratio {rho_all:.3f} population, "
          + ", ".join(f"{k} {v:.3f}" for k, v in sorted(class_rho.items()))
          + f"; grass max {grass_max:.3f} against tree p10 {tree_p10:.3f}, "
          f"separated: {rho_resolved}")
    print(f"  recombination residual {recombined - value:+.5f} at "
          f"z1 = {z1:.5f}; the naive anchor on the measured range's "
          f"{g1:.4f} would have left {naive_residual:+.5f}")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
