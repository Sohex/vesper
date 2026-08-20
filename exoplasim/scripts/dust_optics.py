#!/usr/bin/env python3
"""Two-band mineral dust optics for this world, and the sign of its forcing.

    python exoplasim/scripts/dust_optics.py            # from committed indices
    python exoplasim/scripts/dust_optics.py --digitize # re-extract from the PDF

Produces `analysis/dust_optics.json`: band-averaged single-scattering albedo,
asymmetry parameter, mass extinction efficiency and the critical surface albedo
at which dust forcing changes sign, for each refractive-index dataset.

## Why this exists as a script rather than a note

`notes/dust.md` quotes these numbers as measured. They were originally computed
in a session scratchpad, which meant a note in the repository cited values that
existed nowhere in the repository. That is the failure this file closes: the
inputs are committed, the computation is committed, and the product carries its
provenance.

## The question it answers

Aerosol top-of-atmosphere forcing reverses sign at a critical surface albedo:
below it the layer scatters and cools, above it the layer is darker than the
ground and warms. This world's closed-basin fill is bright -- salt crust at 0.40
to 0.50 -- so which side of that threshold it sits on decides whether dust warms
or cools over the surfaces that make this world unusual.

## Band-2 indices are the hard part

ExoPlaSim's shortwave split is at 0.75 um, and most of this star's flux is above
it -- `lib/stellar.py` has the share, computed the way `solarini` computes it
rather than from a blackbody. That is the band where measured dust indices
effectively run out:

  Di Biagio 2019   0.37-0.95 um   measured, 19 samples          -> band 1
  Rocha-Lima 2018  0.95-2.45 um   measured, but FIGURE ONLY     -> band 2
  OPAC             0.25-40 um     continuous, but 1-2x more
                                  absorbing than measured

Rocha-Lima's spectral k exists only as plotted points in its Fig. 10. It is
digitised here rather than read off by eye, and validated against the four values
the paper states in its text.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dust_indices import NAMES, OPAC_PATH, ROCHALIMA_PATH, indices  # noqa: E402
from mie_dust import lognormal_integrate  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
import sys as _sys
if str(ROOT / "lib") not in _sys.path:
    _sys.path.insert(0, str(ROOT / "lib"))
from paths import rel  # noqa: E402
from stellar import band1_fraction, spectrum_paths  # noqa: E402
DATA = ROOT / "exoplasim" / "data" / "dust"
OUT = ROOT / "analysis" / "dust_optics.json"
SPECTRUM = ROOT / "exoplasim" / "inputs" / "stellarspectra" / "k25v.dat"
PDF = ROOT / "references" / "rochalima2018-fennec-saharan-dust.pdf"

# The refractive-index datasets themselves live in `dust_indices.py`, which is
# also what `aeolian/config/dust.yaml` selects between and what every other
# consumer resolves a dataset name through. This file computes the table; it
# does not own the tables it computes from.

# Balkanski et al. (2007) source distribution: modal (number-median) diameter
# 0.59 um, sigma 2.0. Density 2.6 g/cm3 as OPAC uses for all mineral components.
R_MOD_UM, SIGMA_G, RHO_G_CM3 = 0.295, 2.0, 2.6

# Surfaces on this world, for the sign test. Salt crust is the decisive one.
# Vegetated and playa come from config/planet.yaml (rule 2); ocean is the
# model's open-water value and salt crust the export rock table's evaporite.
import yaml as _yaml
_planet = _yaml.safe_load(
    (Path(__file__).resolve().parents[2] / "config" / "planet.yaml")
    .read_text(encoding="utf-8"))
SURFACES = {
    "ocean": 0.07,
    "vegetated land": float(_planet["model"]["vegetation_albedo"]),
    "playa fill": float(_planet["model"]["lithology_albedo_overrides"]
                        ["playa_clastic"]["albedo"]),
    "salt crust (bright)": 0.50,
}

BAND_SPLIT_UM = 0.75
BAND_LO_UM, BAND_HI_UM = 0.34, 4.00   # 0.34 um is the spectrum file's floor;
                                      # 4.0 um is ExoPlaSim's shortwave ceiling


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stellar_weights():
    """Flux per wavelength bin, for shape WITHIN a band. Not the band split.

    The k25v grid is NOT uniform -- 0.01 um in the visible, 0.04 um in the
    infrared -- so summing flux without multiplying by the bin width overstates
    the infrared's share of a band.

    This weights the Mie integration across each band and nothing else. The
    band SPLIT comes from `lib/stellar.py`, which reproduces `solarini`'s own
    integration of the hi-res file; integrating this low-resolution file across
    the split instead gave 0.3777, which is the 0.34-14.01 um truncation
    speaking rather than the star.
    """
    d = np.loadtxt(SPECTRUM, skiprows=1)
    return d[:, 0], d[:, 1] * np.gradient(d[:, 0])


def critical_surface_albedo(ssa: float, beta: float) -> float:
    """Surface albedo at which TOA aerosol forcing changes sign.

    Two-stream condition (1-a)^2 = 2a(1-w)/(b w), taking the physical root.
    Above `a_c` the layer is darker than the ground and warms.
    """
    b = 2.0 + 2.0 * (1.0 - ssa) / (beta * ssa)
    return (b - np.sqrt(b * b - 4.0)) / 2.0


def band_average(lo, hi, k_of_lam, n_of_lam, lam_s, f_s, n_grid=16, n_r=240):
    """Flux-weighted band mean of extinction, single-scattering albedo and g."""
    grid = np.linspace(lo, hi, n_grid)
    w = np.interp(grid, lam_s, f_s)
    w = w / w.sum()
    ext = sca = gsc = 0.0
    for lam, ww in zip(grid, w):
        _, ssa, g, mee = lognormal_integrate(
            lam, n_of_lam(lam), k_of_lam(lam), R_MOD_UM, SIGMA_G,
            0.01, 10.0, RHO_G_CM3, n_r=n_r)
        ext += mee * ww
        sca += mee * ssa * ww
        gsc += g * mee * ssa * ww
    return ext, sca / ext, gsc / sca


def digitize_fig10(pdf: Path, out: Path) -> None:
    """Extract Rocha-Lima Fig. 10 by colour-thresholding the marker series.

    The page is rendered at 300 dpi, the plot frame located as the longest runs
    of dark pixels, and each series masked by colour. Error-bar stems are 2-3 px
    wide and markers and caps are 6+, so keeping only pixels in horizontal runs
    of 6 or more drops the stems; caps sit symmetrically about the marker, so the
    per-column median then lands on the marker itself.

    Validated against the four values the paper states in text: Algeria mixed
    mode reads 0.0032 at 450 nm against a stated 0.0030, and 0.0005 at 850 nm
    against 0.0005. Accuracy is about +/-0.0005 in k.
    """
    import subprocess
    from PIL import Image

    tmp = out.parent / "_fig10"
    subprocess.run(["pdftoppm", "-f", "13", "-l", "13", "-r", "300", "-png",
                    str(pdf), str(tmp)], check=True)
    im = np.array(Image.open(f"{tmp}-13.png").convert("RGB")).astype(int)

    rows = []
    for panel, (x0, y0) in (("algeria", (540, 280)), ("mauritania", (1280, 280))):
        sub = im[y0:y0 + 550, x0:x0 + 730]
        dark = sub.sum(axis=2) < 250
        cc = np.where(dark.sum(axis=0) > 0.6 * dark.shape[0])[0]
        rr = np.where(dark.sum(axis=1) > 0.6 * dark.shape[1])[0]
        xl, xr, yt, yb = cc[0], cc[-1], rr[0], rr[-1]
        R, G, B = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]
        series = {"fine_mie": (R > 140) & (G < 110) & (B < 110),
                  "fine_tmatrix": (R < 90) & (G < 90) & (B < 90),
                  "mixed_mie": (B > 140) & (R < 110) & (G < 110)}
        for name, mask in series.items():
            m = mask.copy()
            m[:, :xl + 2] = False
            m[:, xr - 1:] = False
            m[:yt + 2, :] = False
            m[yb - 1:, :] = False
            m[50:195, 175:600] = False          # legend box
            keep = np.zeros_like(m)
            for j in range(m.shape[0]):
                row = m[j]
                if not row.any():
                    continue
                d = np.diff(np.concatenate(([0], row.view(np.int8), [0])))
                for a, b in zip(np.where(d == 1)[0], np.where(d == -1)[0]):
                    if b - a >= 6:
                        keep[j, a:b] = True
            for i in range(xl + 2, xr - 1):
                ys = np.where(keep[:, i])[0]
                if len(ys) < 3:
                    continue
                lam_nm = 250.0 + (i - xl) / (xr - xl) * 2400.0
                k = 0.016 + (np.median(ys) - yt) / (yb - yt) * (-0.016)
                rows.append((panel, name, round(lam_nm, 1), round(float(k), 6)))

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="ascii") as fh:
        fh.write("# Rocha-Lima et al. (2018) ACP 18, 1023-1043, Fig. 10.\n")
        fh.write("# Imaginary refractive index of Saharan dust, digitised from\n")
        fh.write("# the published figure -- the paper tabulates no k anywhere,\n")
        fh.write("# and its data are available only on request from the author.\n")
        fh.write("# Real part is an ASSUMED constant 1.56 in that paper, not a\n")
        fh.write("# retrieval. Accuracy about +/-0.0005 in k; validated against\n")
        fh.write("# the four values the text states. Regenerate with --digitize.\n")
        fh.write("panel,series,wavelength_nm,k\n")
        for r in rows:
            fh.write(f"{r[0]},{r[1]},{r[2]},{r[3]}\n")
    print(f"wrote {out} ({len(rows)} points)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--digitize", action="store_true",
                    help="re-extract Rocha-Lima Fig. 10 from the PDF")
    ap.add_argument("--output", type=Path, default=OUT)
    args = ap.parse_args()

    if args.digitize:
        if not PDF.is_file():
            raise SystemExit(f"{PDF} not present; see references/INDEX.md")
        digitize_fig10(PDF, ROCHALIMA_PATH)

    for p in (SPECTRUM, OPAC_PATH, ROCHALIMA_PATH):
        if not p.is_file():
            raise SystemExit(f"missing input {p}")

    lam_s, f_s = stellar_weights()
    b1 = band1_fraction()

    # Every dataset in `dust_indices.NAMES`, over the band it is a candidate
    # for. The config picks one per band from exactly this table, so a dataset
    # added to `dust_indices` and not here would be selectable and unpriced.
    bands = {"band 1": (BAND_LO_UM, BAND_SPLIT_UM),
             "band 2": (BAND_SPLIT_UM, BAND_HI_UM)}
    candidates = [("Di Biagio 2019 measured", "band 1"),
                  ("Rocha-Lima Algeria fine", "band 2"),
                  ("Rocha-Lima Mauritania fine", "band 2"),
                  ("OPAC", "band 1"),
                  ("OPAC", "band 2")]
    unpriced = set(NAMES) - {name for name, _ in candidates}
    if unpriced:
        raise SystemExit(
            f"dust_indices knows {sorted(unpriced)} and this file does not band "
            f"average them, so aeolian/config/dust.yaml could select an entry "
            f"that analysis/dust_optics.json has no row for.")
    cases = []
    for name, band in candidates:
        n_of, k_of = indices(name)
        lo, hi = bands[band]
        cases.append((name, band, lo, hi, k_of, n_of))

    results = []
    print(f"band 1 carries {b1*100:.1f}% of stellar flux, band 2 {100-b1*100:.1f}%\n")
    print(f"{'indices':30}{'band':8}{'MEE':>8}{'ssa':>8}{'g':>8}{'beta':>8}{'a_crit':>9}")
    for label, band, lo, hi, kf, nf in cases:
        mee, ssa, g = band_average(lo, hi, kf, nf, lam_s, f_s)
        beta = (1.0 - g) / 2.0
        a_c = critical_surface_albedo(ssa, beta)
        print(f"{label:30}{band:8}{mee:>8.3f}{ssa:>8.4f}{g:>8.4f}"
              f"{beta:>8.4f}{a_c:>9.3f}")
        results.append({"indices": label, "band": band,
                        "wavelength_um": [lo, hi],
                        "mass_extinction_efficiency_m2_g": round(mee, 4),
                        "single_scattering_albedo": round(ssa, 4),
                        "asymmetry_parameter": round(g, 4),
                        "backscatter_fraction": round(beta, 4),
                        "critical_surface_albedo": round(a_c, 4),
                        "warms_over": sorted(
                            [s for s, a in SURFACES.items() if a > a_c])})

    payload = {
        "note": "Band-averaged dust optics and the surface albedo at which dust "
                "forcing changes sign. Generated by exoplasim/scripts/"
                "dust_optics.py; do not edit.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "stellar_flux_fraction_band1": round(float(b1), 4),
        "stellar_flux_fraction_band1_source":
            "lib/stellar.py, reproducing radmod.f90:solarini on "
            + spectrum_paths()[1].name,
        "size_distribution": {
            "source": "Balkanski et al. 2007",
            "number_median_radius_um": R_MOD_UM,
            "sigma_g": SIGMA_G, "density_g_cm3": RHO_G_CM3},
        "surfaces": SURFACES,
        "inputs": {p.name: sha256(p) for p in (SPECTRUM, OPAC_PATH, ROCHALIMA_PATH)},
        "results": results,
        "caveats": [
            "Rocha-Lima k is digitised from a figure, about +/-0.0005.",
            "n taken as 1.52 (Di Biagio, measured); Rocha-Lima assumed 1.56.",
            "k held flat above 2.45 um, which carries about 5% of the flux.",
            "a_crit depends on the assumed size distribution.",
            "Rocha-Lima measured Saharan silicate; this world's source is "
            "evaporite, which is less absorbing and pushes a_crit higher.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {rel(args.output)}")


if __name__ == "__main__":
    main()
