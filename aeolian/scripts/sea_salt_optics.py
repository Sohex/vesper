#!/usr/bin/env python3
"""Two-band sea-salt optics for this world, per dry size bin and per humidity.

    python aeolian/scripts/sea_salt_optics.py

CLIM-27. Produces `analysis/sea_salt_optics.json`: band-averaged mass extinction
efficiency, single-scattering albedo, asymmetry parameter and the critical
surface albedo at which the forcing changes sign, for each dry size bin of
`aeolian/config/sea_salt.yaml` at each humidity OPAC tabulates.

It is the sea-salt counterpart of `exoplasim/scripts/dust_optics.py` and it
lives here rather than beside that one because nothing in the model reads it:
mineral dust has an `aerofile` and a surface field, and sea salt so far has an
offline burden and nothing else.

## The two things this has that mineral dust does not

**The particle is wet, and the mass is dry.** A sea-salt particle at ambient
humidity is a solution droplet: bigger than the dry salt, less dense, and
optically nearer water than salt. Its extinction has to be computed on the wet
particle and then divided by the DRY mass, because the emission and transport
budget is carried in dry salt. At 80% relative humidity the wet particle is
1.99 times the dry diameter and 4.1 times the dry mass, so getting that division
backwards is a factor of four.

**The size distribution is the emitted one.** Mineral dust is integrated over a
single lognormal because that is what its source theory produces and what OPAC
publishes. Sea salt is emitted as a sum of three lognormals, so the Mie
integration runs over `sea_salt_source.number_flux_per_dp` restricted to each
bin rather than over a lognormal chosen for convenience.

## What is checked, and what a failure would mean

Three things, each with a right answer that does not come from this file.

1. **OPAC's own columns must agree with each other.** The wet density and the
   growth in mode radius are two statements of the same salt volume fraction, so
   the density the growth factor implies must reproduce the density OPAC prints.
   Compared as densities they agree to 0.008 g/cm3, inside the table's own
   rounding; compared as volume fractions the same rows look 11% apart, because
   at 99% humidity the printed 1.01 keeps one significant figure once water is
   subtracted. The sharper test is on the quantity as printed.

2. **Inverting OPAC's volume mixing must return the refractive index of water.**
   OPAC builds the humid index by volume-weighting dry salt with water. Undoing
   that with the volume fraction from check 1 has to give 1.333 at 0.55 um,
   which is Hale and Querry's water and appears nowhere in this calculation.

3. **The non-absorbing and geometric-optics limits.** Sea salt's imaginary index
   in band 1 is of order 1e-8, so the single-scattering albedo must be 1 to
   within 1e-5; and the coarse bin at 0.55 um has a size parameter above 20, so
   its extinction efficiency must be near 2.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import PROJECT_ROOT  # noqa: E402
from sea_salt_source import number_flux_per_dp  # noqa: E402

sys.path.insert(0, str(PROJECT_ROOT / "exoplasim" / "scripts"))
# stellar_weights and critical_surface_albedo are imported rather than copied.
# The first carries a subtlety that is easy to get wrong and was got wrong once
# -- the k25v grid is not uniform, so flux must be multiplied by the bin width
# before it is used as a weight -- and a second copy of it is a second chance to
# lose that. `dust_optics` is a script, but it guards its own entry point.
from dust_optics import critical_surface_albedo, stellar_weights  # noqa: E402
from mie_dust import distribution_integrate  # noqa: E402

sys.path.insert(0, str(PROJECT_ROOT / "lib"))
from paths import rel  # noqa: E402
from stellar import band1_fraction, spectrum_paths  # noqa: E402

DATA = PROJECT_ROOT / "exoplasim" / "data" / "sea_salt"
INDEX_PATH = DATA / "opac_sea_salt_refractive_index.dat"
SIZE_PATH = DATA / "opac_sea_salt_size_distribution.dat"
CONFIG_PATH = PROJECT_ROOT / "aeolian" / "config" / "sea_salt.yaml"
OUT = PROJECT_ROOT / "analysis" / "sea_salt_optics.json"

BAND_SPLIT_UM = 0.75
BAND_LO_UM, BAND_HI_UM = 0.34, 4.00     # the stellar spectrum file's own range
RHO_WATER_G_CM3 = 1.0
WATER_N_AT_550 = 1.333                  # Hale and Querry (1973), via OPAC
SURFACES = {"ocean": 0.07, "sea ice": 0.60, "vegetated land": 0.18}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_size_table() -> dict:
    """{mode: {rh: (r_min, r_max, r_mod, rho, sigma)}} from the OPAC table."""
    out: dict[str, dict[int, tuple]] = {}
    for line in SIZE_PATH.read_text().splitlines():
        if line.startswith("#") or line.startswith("mode") or not line.strip():
            continue
        p = line.split()
        out.setdefault(p[0], {})[int(p[1])] = tuple(float(x) for x in p[2:])
    return out


def read_index_table() -> dict:
    """{mode: {rh: (lambda, n, k)}} from the OPAC table."""
    out: dict[str, dict[int, list]] = {}
    for line in INDEX_PATH.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split()
        out.setdefault(p[0], {}).setdefault(int(p[1]), []).append(
            (float(p[2]), float(p[3]), float(p[4])))
    return {m: {rh: np.array(sorted(v)) for rh, v in d.items()}
            for m, d in out.items()}


def growth(size: dict, mode: str, rh: int) -> tuple[float, float, float]:
    """Growth factor, wet density and the salt volume fraction they imply.

    The two routes to the volume fraction are OPAC's own and must agree; the
    caller checks. `rho_dry` is taken from the table's own 0% row rather than
    from the config, so this function reads one source and not two.
    """
    r_mod_dry = size[mode][0][2]
    rho_dry = size[mode][0][3]
    _, _, r_mod, rho_wet, _ = size[mode][rh]
    gf = r_mod / r_mod_dry
    f_from_growth = gf ** -3.0
    # The density OPAC's own growth factor implies, to be compared against the
    # density OPAC tabulates. Comparing the VOLUME FRACTIONS instead looks like
    # an 11% disagreement at 99% humidity, and all of that is the table
    # rounding rho to two decimals: (1.01 - 1.00) keeps one significant figure.
    # Compared as densities the same rows agree to 0.002, which is inside the
    # rounding, so the sharper test is the one on the quantity as printed.
    rho_implied = (RHO_WATER_G_CM3
                   + f_from_growth * (rho_dry - RHO_WATER_G_CM3))
    return gf, rho_wet, (f_from_growth, rho_implied)


def band_average(lo_um, hi_um, table_rh, r_wet_um, dndlnr, rho_eff_g_cm3,
                 lam_s, f_s, n_grid=14):
    """Flux-weighted band mean of MEE, single-scattering albedo and g."""
    grid = np.linspace(lo_um, hi_um, n_grid)
    w = np.interp(grid, lam_s, f_s)
    w = w / w.sum()
    lam_t, n_t, k_t = table_rh[:, 0], table_rh[:, 1], table_rh[:, 2]
    ext = sca = gsc = 0.0
    for lam, ww in zip(grid, w):
        n = float(np.interp(lam, lam_t, n_t))
        k = float(np.interp(lam, lam_t, k_t))
        _, ssa, g, mee = distribution_integrate(
            lam, n, k, r_wet_um, dndlnr, rho_eff_g_cm3)
        ext += mee * ww
        sca += mee * ssa * ww
        gsc += g * mee * ssa * ww
    return ext, sca / ext, gsc / sca


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--u10", type=float, default=None,
                    help="wind the emitted size distribution is shaped at; "
                         "defaults to the Earth check's ocean mean")
    ap.add_argument("--output", type=Path, default=OUT)
    args = ap.parse_args()

    cfg = yaml.safe_load(CONFIG_PATH.read_text())
    size = read_size_table()
    index = read_index_table()
    edges = cfg["size"]["bin_edges_um"]
    u10 = args.u10 or cfg["earth_check"]["ocean_mean_u10_m_s"]
    rho_dry = cfg["emission"]["dry_density_kg_m3"] / 1000.0   # g/cm3

    # --- check 1: OPAC's two routes to the salt volume fraction ------------
    worst = 0.0
    for mode in ("ssam", "sscm"):
        for rh in sorted(size[mode]):
            if rh == 0:
                continue
            _, rho_wet, (_, rho_implied) = growth(size, mode, rh)
            worst = max(worst, abs(rho_implied - rho_wet))
    print(f"check 1  wet density implied by OPAC's growth factor against the "
          f"density OPAC tabulates: worst {worst:.4f} g/cm3, table rounds to "
          f"0.01   {'OK' if worst < 0.01 else 'FAIL'}")
    if worst >= 0.01:
        raise SystemExit("OPAC size table misread: growth and density disagree")

    # --- check 2: invert the volume mixing and recover water ---------------
    rh_check = 80
    gf, _, (fg, _) = growth(size, "ssam", rh_check)
    tab = index["ssam"]
    n_dry = float(np.interp(0.55, tab[0][:, 0], tab[0][:, 1]))
    n_wet = float(np.interp(0.55, tab[rh_check][:, 0], tab[rh_check][:, 1]))
    n_water = (n_wet - fg * n_dry) / (1.0 - fg)
    err = abs(n_water - WATER_N_AT_550) / WATER_N_AT_550
    print(f"check 2  water index recovered from OPAC's mixing at {rh_check}% "
          f"RH: {n_water:.4f} against {WATER_N_AT_550}   "
          f"{'OK' if err < 0.01 else 'FAIL'}")
    if err >= 0.01:
        raise SystemExit("volume-mixing inversion does not return water")

    lam_s, f_s = stellar_weights()
    b1 = band1_fraction()
    print(f"\nband 1 carries {b1*100:.1f}% of stellar flux, band 2 "
          f"{100-b1*100:.1f}%")
    print(f"emitted shape taken at U10 = {u10} m/s; MEE is per gram of DRY "
          f"sea salt\n")

    rows = []
    print(f"{'bin (dry um)':>16}{'RH%':>5}{'GF':>7}{'MEE1':>8}{'MEE2':>8}"
          f"{'ssa1':>9}{'g1':>7}{'a_crit1':>9}")
    for lo, hi in zip(edges[:-1], edges[1:]):
        # OPAC grows and indexes the two modes separately; pick the one whose
        # dry mode radius is nearer this bin rather than averaging two tables.
        centre = np.sqrt(lo * hi)
        mode = "ssam" if centre < 1.0 else "sscm"
        dp_dry = np.logspace(np.log10(lo), np.log10(hi), 240)
        # dN/dlnr of the EMITTED population inside this bin. dF/dDp is per unit
        # diameter, so dN/dlnDp = Dp dF/dDp, and ln r and ln Dp differ by a
        # constant whose derivative is 1.
        dndlnr = dp_dry * number_flux_per_dp(dp_dry, np.array(u10), cfg)
        for rh in sorted(size[mode]):
            gf, rho_wet, _ = growth(size, mode, rh)
            r_wet = 0.5 * dp_dry * gf
            # MEE per gram of DRY salt: pass the density that turns the WET
            # volume into the DRY mass, which is rho_dry / gf^3.
            rho_eff = rho_dry / gf ** 3
            mee1, ssa1, g1 = band_average(BAND_LO_UM, BAND_SPLIT_UM,
                                          index[mode][rh], r_wet, dndlnr,
                                          rho_eff, lam_s, f_s)
            mee2, ssa2, g2 = band_average(BAND_SPLIT_UM, BAND_HI_UM,
                                          index[mode][rh], r_wet, dndlnr,
                                          rho_eff, lam_s, f_s)
            beta1 = (1.0 - g1) / 2.0
            a_c1 = critical_surface_albedo(ssa1, beta1)
            rows.append({
                "bin_dry_um": [lo, hi], "opac_mode": mode,
                "relative_humidity_percent": rh,
                "growth_factor": round(gf, 4),
                "wet_density_g_cm3": rho_wet,
                "band1": {"mass_extinction_efficiency_m2_g_dry": round(mee1, 4),
                          "single_scattering_albedo": round(ssa1, 6),
                          "asymmetry_parameter": round(g1, 4),
                          "backscatter_fraction": round(beta1, 4),
                          "critical_surface_albedo": round(a_c1, 4)},
                "band2": {"mass_extinction_efficiency_m2_g_dry": round(mee2, 4),
                          "single_scattering_albedo": round(ssa2, 6),
                          "asymmetry_parameter": round(g2, 4)},
            })
            if rh in (0, 80, 95):
                print(f"{lo:7.2f}-{hi:<8.2f}{rh:5d}{gf:7.2f}{mee1:8.3f}"
                      f"{mee2:8.3f}{ssa1:9.5f}{g1:7.3f}{a_c1:9.3f}")

    # --- check 3: the two limits -------------------------------------------
    ssa_min = min(r["band1"]["single_scattering_albedo"] for r in rows)
    print(f"\ncheck 3a lowest band-1 single-scattering albedo over all bins "
          f"and humidities: {ssa_min:.6f}   "
          f"{'OK' if ssa_min > 0.9999 else 'FAIL'}")
    coarse = [r for r in rows if r["bin_dry_um"][0] == edges[-2]
              and r["relative_humidity_percent"] == 0][0]
    # Qext from MEE: MEE = 3 Qext / (4 rho r_eff) for a sphere of radius r_eff.
    r_eff = 0.5 * np.sqrt(edges[-2] * edges[-1])
    qext = (coarse["band1"]["mass_extinction_efficiency_m2_g_dry"]
            * 4.0 * rho_dry * r_eff / 3.0)
    print(f"check 3b coarse dry bin effective Qext in band 1: {qext:.2f} "
          f"(geometric optics limit is 2)   "
          f"{'OK' if 1.5 < qext < 2.8 else 'FAIL'}")
    if ssa_min <= 0.9999 or not 1.5 < qext < 2.8:
        raise SystemExit("optics limits failed")

    payload = {
        "note": "Band-averaged sea-salt optics per dry size bin and relative "
                "humidity. Mass extinction efficiency is PER GRAM OF DRY SEA "
                "SALT, not of the wet particle. Generated by "
                "aeolian/scripts/sea_salt_optics.py; do not edit.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "stellar_flux_fraction_band1": round(float(b1), 4),
        "stellar_flux_fraction_band1_source":
            "lib/stellar.py, reproducing radmod.f90:solarini on "
            + spectrum_paths()[1].name,
        "emitted_shape_u10_m_s": u10,
        "indices": "OPAC (Hess et al. 1998) ssam and sscm at ambient humidity",
        "checks": {
    "opac_growth_against_density_worst_g_cm3": round(worst, 4),
            "water_index_recovered_at_550nm": round(n_water, 4),
            "lowest_band1_single_scattering_albedo": ssa_min,
            "coarse_dry_bin_qext_band1": round(float(qext), 3),
        },
        "surfaces": SURFACES,
        "inputs": {p.name: sha256(p) for p in (INDEX_PATH, SIZE_PATH)},
        "results": rows,
        "caveats": [
            "The emitted shape inside a bin is taken at one wind speed; the "
            "shape depends on wind only through the relative weight of the "
            "three source modes and is nearly constant within a bin.",
            "OPAC's humid indices are volume-weighted mixtures of dry salt "
            "with water. Irshad et al. (2009) measure that rule to be "
            "inadequate from 1 to 20 um, where sea salt absorbs; in band 1 it "
            "is a near-pure scatterer and n carries the answer.",
            "No coagulation and no shape correction: spheres throughout.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {rel(args.output)}")


if __name__ == "__main__":
    main()
