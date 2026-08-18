#!/usr/bin/env python3
"""What the dust burden is worth radiatively, shortwave AND longwave.

    python exoplasim/scripts/dust_forcing.py

DUST-2. `dust_optics.py` answers whether dust warms or cools over a given
surface; this answers by how much, at the burden DUST-1 actually produces, and
prices the thermal-infrared term that the shortwave-only picture leaves out.

## Why the longwave is the point

Mineral dust scatters sunlight, which cools, and absorbs and re-emits terrestrial
infrared, which warms. Only the first of those exists in ExoPlaSim: `radmod.f90`
puts the aerosol in its two shortwave bands and the longwave solver has no
aerosol term at all. So before DUST-3 puts an emission scheme in the model, this
has to say how large the missing half is -- because a model that can cool with
dust and cannot warm with it may be further from the truth than one with no dust
at all, and this project already has a narrow flux window to protect.

## The two things that decide the answer

**The size distribution.** Longwave interaction needs particles comparable to a
10 um wavelength, so it is carried almost entirely by the coarse mode, while the
shortwave is carried by the fine one. This is where our own components disagree
with each other: `analysis/dust_optics.json` is computed for a Balkanski
number-median radius of 0.295 um, and `aeolian/config/dust.yaml` emits a Kok
distribution with a volume-median DIAMETER of 3.4 um and a geometric sigma of 3.
The atmospheric burden sits between them -- coarse settles out fastest -- so both
are computed here and reported as a bracket rather than one being chosen.

**The surface albedo.** Shortwave forcing changes sign at a critical albedo, and
this world's closed-basin fill is bright. Over salt crust the shortwave term is
already small or positive, so the longwave decides the sign outright.

## What it does not do

No spectral overlap with water vapour or CO2: the longwave term here is computed
against a clear-sky window approximation, which OVERSTATES it, because part of
the band the dust absorbs in is already opaque. Treat the longwave number as an
upper bound and the net as a bracket. A line-by-line treatment is not available
here and would be precision theatre against a burden uncertain by its own factor
of several.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mie_dust import lognormal_integrate  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "lib") not in sys.path:
    sys.path.insert(0, str(ROOT / "lib"))
from paths import rel  # noqa: E402

OPAC = ROOT / "exoplasim" / "data" / "dust" / "opac_mineral_refractive_index.dat"
OUT = ROOT / "analysis" / "dust_forcing.json"

# Declared geometry, all of it stated rather than fitted.
#
# DIFFUSIVITY is the standard two-stream factor for a diffuse longwave field.
# LAPSE_K_PER_KM and the dust scale height together set the temperature the layer
# radiates at, which is what the longwave term is proportional to; the scale
# height is `aeolian/config/dust.yaml`'s own, and for an exponential profile the
# mass-weighted mean height IS the scale height.
# TRANSMISSION is the clear-sky shortwave transmission of the atmosphere above
# the layer, which the two-stream forcing expression squares.
DIFFUSIVITY = 1.66
LAPSE_K_PER_KM = 6.5
TRANSMISSION = 0.79
LW_LO_UM, LW_HI_UM = 4.0, 40.0
# The atmospheric window, and how much of it survives between a dust layer at the
# scale height and the ground. Outside 8-12 um the air below the dust is already
# opaque, so downward emission from the dust adds nothing there; inside it, water
# vapour still absorbs some. 0.7 is declared, bracketed 0.5 to 0.9, and it scales
# the surface longwave term linearly.
WINDOW_LO_UM, WINDOW_HI_UM = 8.0, 12.0
WINDOW_TRANSMITTANCE = 0.7
SW_BANDS_UM = ((0.34, 0.75), (0.75, 4.0))

# The two ends of the size bracket. `fine` is what the optical depth in
# `dust_baseline.json` was computed with; `emitted` is what the emission scheme
# actually produces, converted from a volume-median diameter to the number-median
# radius a lognormal integrator wants by ln D_n = ln D_v - 3 ln^2 sigma.
SIZE_ENDS = {
    "fine (Balkanski, what the reported AOD assumes)":
        {"r_mod_um": 0.295, "sigma_g": 2.0},
    "emitted (Kok, volume-median diameter 3.4 um)":
        {"r_mod_um": 0.5 * 3.4 * np.exp(-3.0 * np.log(3.0) ** 2), "sigma_g": 3.0},
}
RHO_G_CM3 = 2.6
SURFACES = {"ocean": 0.07, "vegetated land": 0.18,
            "playa fill": 0.40, "salt crust (bright)": 0.50}


def opac_indices():
    """OPAC mineral n and k, 0.25 to 40 um. k ships negative; sign is flipped."""
    d = np.loadtxt(OPAC)
    lam, n, k = d[:, 0], d[:, 1], np.abs(d[:, 2])
    order = np.argsort(lam)
    return lam[order], n[order], k[order]


def planck(lam_um, temperature_k):
    """Spectral radiant EXITANCE, W/m2 per metre of wavelength.

    `c1` is 2*pi*h*c^2, so the pi is already in it and this is a flux out of a
    surface rather than a radiance. Multiplying by pi again overstates every
    longwave number by 3.14, which is how the first version of this read 39.7
    W/m2 against a hard ceiling of sigma(Ts^4 - Td^4) times the emissivity.
    `_check_planck` asserts the integral against Stefan-Boltzmann.
    """
    lam = np.asarray(lam_um, dtype=float) * 1e-6
    c1, c2 = 3.7418e-16, 1.4388e-2
    return c1 / (lam ** 5 * (np.exp(c2 / (lam * temperature_k)) - 1.0))


def _check_planck(temperature_k=289.8, tol=0.02) -> float:
    """Integral of `planck` against sigma T^4. Raises if the normalisation slips.

    Over 0.1 to 500 um essentially all of a 290 K body's emission is captured, so
    this is a real check on the constant and on the units of the integration
    variable, not a formality.
    """
    lam = np.exp(np.linspace(np.log(0.1), np.log(500.0), 4000))
    total = float(np.trapezoid(planck(lam, temperature_k), lam * 1e-6))
    sigma_t4 = 5.670374e-8 * temperature_k ** 4
    if abs(total / sigma_t4 - 1.0) > tol:
        raise SystemExit(
            f"planck() integrates to {total:.2f} W/m2 against a Stefan-Boltzmann "
            f"{sigma_t4:.2f} at {temperature_k} K. The normalisation is wrong.")
    return total / sigma_t4


def band_optics(lo_um, hi_um, weights_of_lam, dist, n_grid=24, n_r=200):
    """Weighted band-mean extinction, single-scattering albedo and asymmetry."""
    lam_t, n_t, k_t = opac_indices()
    grid = np.linspace(lo_um, hi_um, n_grid)
    w = weights_of_lam(grid)
    w = w / w.sum()
    ext = sca = gsc = 0.0
    for lam, ww in zip(grid, w):
        n = float(np.interp(lam, lam_t, n_t))
        k = float(np.interp(lam, lam_t, k_t))
        _, ssa, g, mee = lognormal_integrate(
            lam, n, k, dist["r_mod_um"], dist["sigma_g"],
            0.01, 25.0, RHO_G_CM3, n_r=n_r)
        ext += mee * ww
        sca += mee * ssa * ww
        gsc += g * mee * ssa * ww
    return ext, sca / ext, gsc / max(sca, 1e-30)


def backscatter_fraction(g):
    """Hemispheric upscatter fraction from the asymmetry parameter, Wiscombe-Grams."""
    return 0.5 * (1.0 - g)


def shortwave_forcing(tau_ext, ssa, beta, albedo, insolation):
    """Chylek and Coakley (1974) two-stream TOA forcing for a thin layer.

    Negative is cooling. The scattering term cools in proportion to how dark the
    ground is, and the absorbing term warms in proportion to how bright it is,
    which is the whole reason a critical albedo exists.
    """
    tau_sca, tau_abs = tau_ext * ssa, tau_ext * (1.0 - ssa)
    return -insolation * TRANSMISSION ** 2 * (
        (1.0 - albedo) ** 2 * beta * tau_sca - 2.0 * albedo * tau_abs)


def longwave_forcing(column_kg_m2, dist, t_surface, t_layer, n_grid=40):
    """Greybody TOA longwave warming from an absorbing layer aloft.

    Spectrally resolved rather than band-averaged, because the dust absorption
    and the Planck function both vary strongly across the thermal infrared and
    their product is what matters:

        dF = integral [ B(T_s) - B(T_d) ] * (1 - exp(-D * tau_abs)) dlam

    CLEAR-SKY WINDOW APPROXIMATION, so this is an UPPER bound: part of the band
    is already opaque to water vapour and CO2, and dust absorbing where the
    atmosphere already absorbs adds nothing.
    """
    lam_t, n_t, k_t = opac_indices()
    grid = np.linspace(LW_LO_UM, LW_HI_UM, n_grid)
    total = 0.0
    tau_abs_mean = 0.0
    b_s, b_d = planck(grid, t_surface), planck(grid, t_layer)
    for i, lam in enumerate(grid):
        n = float(np.interp(lam, lam_t, n_t))
        k = float(np.interp(lam, lam_t, k_t))
        _, ssa, _, mee = lognormal_integrate(
            lam, n, k, dist["r_mod_um"], dist["sigma_g"],
            0.01, 25.0, RHO_G_CM3, n_r=200)
        tau_abs = mee * 1e3 * column_kg_m2 * (1.0 - ssa)   # mee is m2/g
        emissivity = 1.0 - np.exp(-DIFFUSIVITY * tau_abs)
        total += (b_s[i] - b_d[i]) * emissivity
        tau_abs_mean += tau_abs
    dlam = (LW_HI_UM - LW_LO_UM) / (n_grid - 1) * 1e-6
    return float(total * dlam), float(tau_abs_mean / n_grid)


def surface_forcing(column_kg_m2, dist, albedo, t_layer, insolation,
                    sw_bands, n_grid=24):
    """Shortwave and longwave forcing AT THE SURFACE, which is what Penman reads.

    Not the top-of-atmosphere numbers, and the difference decides the sign. An
    absorbing layer removes energy from the beam that never reaches the ground,
    so the surface shortwave term is negative even where the TOA term is positive
    over bright ground -- the TOA sign there is energy retained in the ATMOSPHERE.

    Shortwave, per Chylek and Coakley's surface form: the ground loses what the
    layer scatters back and what it absorbs, scaled by how much of the beam it
    would have kept.

        dF_sw,surf = -S T^2 (1 - a) (beta tau_sca + tau_abs)

    Longwave: the ground gains downward emission from the layer, but only in the
    window, because outside it the air between is already opaque and adds nothing.
    """
    tau_sca = tau_abs = beta_eff = 0.0
    for (lo, hi), share in sw_bands:
        ext, ssa, g = band_optics(lo, hi, lambda x: np.ones_like(x), dist)
        tau = ext * 1e3 * column_kg_m2
        tau_sca += tau * ssa * share
        tau_abs += tau * (1.0 - ssa) * share
        beta_eff += backscatter_fraction(g) * share
    sw = -insolation * TRANSMISSION ** 2 * (1.0 - albedo) * (
        beta_eff * tau_sca + tau_abs)

    lam_t, n_t, k_t = opac_indices()
    grid = np.linspace(WINDOW_LO_UM, WINDOW_HI_UM, n_grid)
    b = planck(grid, t_layer)
    lw = 0.0
    for i, lam in enumerate(grid):
        n = float(np.interp(lam, lam_t, n_t))
        k = float(np.interp(lam, lam_t, k_t))
        _, ssa, _, mee = lognormal_integrate(
            lam, n, k, dist["r_mod_um"], dist["sigma_g"], 0.01, 25.0,
            RHO_G_CM3, n_r=200)
        emissivity = 1.0 - np.exp(-DIFFUSIVITY * mee * 1e3 * column_kg_m2 * (1.0 - ssa))
        lw += b[i] * emissivity
    dlam = (WINDOW_HI_UM - WINDOW_LO_UM) / (n_grid - 1) * 1e-6
    return float(sw), float(lw * dlam * WINDOW_TRANSMITTANCE)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dust", type=Path,
                    default=ROOT / "aeolian" / "analysis" / "dust_baseline.json")
    ap.add_argument("--config", type=Path, default=ROOT / "config" / "planet.yaml")
    ap.add_argument("--output", type=Path, default=OUT)
    args = ap.parse_args()

    ratio = _check_planck()
    print(f"planck normalisation checked: {ratio:.4f} of sigma T^4\n")
    import yaml
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    dust = json.loads(args.dust.read_text(encoding="utf-8"))
    dust_cfg = yaml.safe_load(
        (ROOT / "aeolian" / "config" / "dust.yaml").read_text(encoding="utf-8"))

    flux = float(config["orbit"]["baseline_flux_earth"])
    insolation = 1361.0 * flux / 4.0
    t_surface = 289.80
    scale_h_m = float(dust_cfg["transport"]["dust_scale_height_m"])
    t_layer = t_surface - LAPSE_K_PER_KM * scale_h_m / 1000.0

    mee_chain = float(dust["mass_extinction_efficiency_m2_kg"])
    central = dust["shelter_bracket"]["central"]
    aod_land = float(central["land_mean_aod"])
    column = aod_land / mee_chain          # kg/m2, the size-independent quantity

    lam_s = None
    d = np.loadtxt(ROOT / "exoplasim" / "inputs" / "stellarspectra" / "k25v.dat",
                   skiprows=1)
    lam_s, f_s = d[:, 0], d[:, 1] * np.gradient(d[:, 0])

    print(f"burden: land-mean AOD {aod_land:.4f} at {mee_chain:.1f} m2/kg")
    print(f"     -> column mass {column * 1e3:.3f} g/m2")
    print(f"layer: {t_surface:.1f} K surface, {t_layer:.1f} K at "
          f"{scale_h_m / 1000:.1f} km, insolation {insolation:.1f} W/m2\n")

    results = []
    for name, dist in SIZE_ENDS.items():
        sw = {}
        tau_sw_total = 0.0
        for (lo, hi), weight in zip(SW_BANDS_UM, (None, None)):
            ext, ssa, g = band_optics(
                lo, hi, lambda x: np.interp(x, lam_s, f_s), dist)
            share = float(np.trapezoid(np.interp(np.linspace(lo, hi, 64), lam_s, f_s),
                                       np.linspace(lo, hi, 64)))
            sw[f"{lo}-{hi} um"] = {
                "mass_extinction_efficiency_m2_g": round(ext, 4),
                "single_scattering_albedo": round(float(ssa), 4),
                "asymmetry_parameter": round(float(g), 4),
                "flux_share_raw": share}
            tau_sw_total += ext * 1e3 * column * share
        norm = sum(v["flux_share_raw"] for v in sw.values())
        tau_sw = 0.0
        ssa_eff = beta_eff = 0.0
        for v in sw.values():
            wgt = v["flux_share_raw"] / norm
            t = v["mass_extinction_efficiency_m2_g"] * 1e3 * column
            tau_sw += t * wgt
            ssa_eff += v["single_scattering_albedo"] * wgt
            beta_eff += backscatter_fraction(v["asymmetry_parameter"]) * wgt
            v.pop("flux_share_raw")
            v["flux_share"] = round(wgt, 4)

        lw, tau_lw = longwave_forcing(column, dist, t_surface, t_layer)
        per_surface = {}
        for surface, albedo in SURFACES.items():
            sw_f = shortwave_forcing(tau_sw, ssa_eff, beta_eff, albedo, insolation)
            per_surface[surface] = {
                "shortwave_w_m2": round(float(sw_f), 3),
                "longwave_w_m2": round(lw, 3),
                "net_w_m2": round(float(sw_f) + lw, 3)}
        results.append({
            "size_distribution": name,
            "r_mod_um": round(dist["r_mod_um"], 4),
            "sigma_g": dist["sigma_g"],
            "shortwave_optical_depth": round(float(tau_sw), 4),
            "shortwave_bands": sw,
            "longwave_absorption_optical_depth": round(tau_lw, 5),
            "longwave_forcing_w_m2": round(lw, 3),
            "by_surface": per_surface})
        print(f"{name}")
        print(f"   SW optical depth {tau_sw:.4f}   LW absorption depth {tau_lw:.5f}")
        print(f"   LONGWAVE  {lw:+7.3f} W/m2   (clear-sky window, an upper bound)")
        for surface, v in per_surface.items():
            print(f"   {surface:<24} SW {v['shortwave_w_m2']:+8.3f}  "
                  f"net {v['net_w_m2']:+8.3f} W/m2")
        print()

    # The global mean is what the error budget wants, and it is not the land
    # mean: the burden thins over ocean, and ocean is where the shortwave cools
    # hardest because the surface underneath is dark. Area-weighted with the
    # planet's own land fraction, and scaled by the ratio of global to land
    # optical depth the transport actually produced.
    aod_global = float(central["global_mean_aod"])
    land_fraction = float(dust["source_map"]["land_fraction"])
    scale = aod_global / aod_land
    global_mean = []
    for r in results:
        land_net = 0.5 * (r["by_surface"]["vegetated land"]["net_w_m2"]
                          + r["by_surface"]["playa fill"]["net_w_m2"])
        ocean_net = r["by_surface"]["ocean"]["net_w_m2"]
        gm = scale * (land_fraction * land_net + (1.0 - land_fraction) * ocean_net)
        r["global_mean_net_w_m2"] = round(float(gm), 3)
        global_mean.append(gm)
        # And what ExoPlaSim would actually produce, which is the shortwave
        # alone. This is the number DUST-3 turns on.
        land_sw = 0.5 * (r["by_surface"]["vegetated land"]["shortwave_w_m2"]
                         + r["by_surface"]["playa fill"]["shortwave_w_m2"])
        ocean_sw = r["by_surface"]["ocean"]["shortwave_w_m2"]
        gm_sw = scale * (land_fraction * land_sw + (1.0 - land_fraction) * ocean_sw)
        r["global_mean_shortwave_only_w_m2"] = round(float(gm_sw), 3)
        r["shortwave_only_error_w_m2"] = round(float(gm_sw - gm), 3)
        print(f"{r['size_distribution'][:34]:<36} global mean net "
              f"{gm:+7.3f}   shortwave only {gm_sw:+7.3f}   "
              f"error {gm_sw - gm:+7.3f} W/m2")
    print(f"\nreopening threshold is 1.5 W/m2 in the global mean -> "
          f"{'CROSSES' if max(abs(g) for g in global_mean) > 1.5 else 'below'}")

    payload = {
        "note": "DUST-2. Shortwave and longwave forcing of the DUST-1 burden. "
                "The longwave is a clear-sky window estimate and is an UPPER "
                "bound; the size distribution is a bracket because this world's "
                "own components disagree about it. Generated by "
                "exoplasim/scripts/dust_forcing.py; do not edit.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "burden": {"land_mean_aod": aod_land,
                   "mass_extinction_efficiency_m2_kg_used_by_chain": mee_chain,
                   "column_mass_kg_m2": column,
                   "source": rel(args.dust)},
        "geometry": {"surface_temperature_k": t_surface,
                     "layer_temperature_k": round(t_layer, 2),
                     "dust_scale_height_m": scale_h_m,
                     "lapse_k_per_km": LAPSE_K_PER_KM,
                     "diffusivity": DIFFUSIVITY,
                     "shortwave_transmission": TRANSMISSION,
                     "mean_insolation_w_m2": round(insolation, 2)},
        "reopening_threshold_w_m2": 1.5,
        "global_mean": {
            "global_mean_aod": aod_global,
            "land_fraction": land_fraction,
            "note": "Area-weighted over ocean and a half-and-half land mix of "
                    "vegetated and playa, scaled by the global-to-land optical "
                    "depth ratio the transport produced. Crude, and the sign is "
                    "what matters rather than the third digit.",
            "net_w_m2_by_size_end": [round(float(g), 3) for g in global_mean]},
        "results": results,
        "caveats": [
            "The longwave term is a clear-sky window estimate. Part of the band "
            "is already opaque to water vapour and CO2, so this OVERSTATES it.",
            "The size distribution is a bracket, not a measurement: the reported "
            "optical depth assumes a fine Balkanski distribution and the emission "
            "scheme produces a much coarser one. The atmospheric burden is "
            "between them and nothing here says where.",
            "OPAC indices run 1 to 2x more absorbing than the measured datasets "
            "in the shortwave, so the shortwave absorption here is an upper bound "
            "too, which cuts the cooling.",
            "Two-stream, optically-thin forcing expressions. At an optical depth "
            "near 1 they are being used past where they are strictly valid.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()


def write_per_cell(dust_nc, albedo, out_path, dist, t_layer, insolation, sw_bands):
    """Per-cell surface forcing, as perturbations to `rss` and `rls`.

    Written so the carve verdict can read them straight into its own Penman
    without re-deriving any optics. The optical properties depend only on the
    size distribution, so they are computed once and applied across the grid;
    only the longwave emissivity is nonlinear in column mass.
    """
    from netCDF4 import Dataset as DS
    with DS(dust_nc) as ds:
        aod = np.asarray(ds["aod_central"][:], dtype=float)
        lat = np.asarray(ds["lat"][:]); lon = np.asarray(ds["lon"][:])
    mee = None
    tau_sca_c = tau_abs_c = beta_eff = 0.0
    for (lo, hi), share in sw_bands:
        ext, ssa, g = band_optics(lo, hi, lambda x: np.ones_like(x), dist)
        if mee is None:
            mee = ext * 1e3                      # m2/kg, for column from AOD
        tau_sca_c += ext * 1e3 * ssa * share
        tau_abs_c += ext * 1e3 * (1.0 - ssa) * share
        beta_eff += backscatter_fraction(g) * share
    column = aod / mee                            # kg/m2 per cell

    drss = -insolation * TRANSMISSION ** 2 * (1.0 - albedo) * (
        beta_eff * tau_sca_c + tau_abs_c) * column

    lam_t, n_t, k_t = opac_indices()
    grid = np.linspace(WINDOW_LO_UM, WINDOW_HI_UM, 24)
    b = planck(grid, t_layer)
    drls = np.zeros_like(column)
    for i, lam in enumerate(grid):
        n = float(np.interp(lam, lam_t, n_t)); k = float(np.interp(lam, lam_t, k_t))
        _, ssa, _, m = lognormal_integrate(lam, n, k, dist["r_mod_um"],
                                           dist["sigma_g"], 0.01, 25.0,
                                           RHO_G_CM3, n_r=200)
        drls += b[i] * (1.0 - np.exp(-DIFFUSIVITY * m * 1e3 * column * (1.0 - ssa)))
    drls *= (WINDOW_HI_UM - WINDOW_LO_UM) / 23 * 1e-6 * WINDOW_TRANSMITTANCE

    with DS(out_path, "w") as ds:
        ds.createDimension("lat", lat.size); ds.createDimension("lon", lon.size)
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        for name, arr, note in (("drss", drss, "surface net shortwave"),
                                ("drls", drls, "surface net longwave")):
            v = ds.createVariable(name, "f8", ("lat", "lon"))
            v[:] = arr; v.units = "W m-2"; v.long_name = f"dust perturbation to {note}"
        ds.note = ("DUST-2/DUST-10. Add to rss and rls before Penman. SURFACE "
                   "forcing, not top-of-atmosphere: the layer absorbs, so the "
                   "ground loses even where TOA gains. Longwave is restricted to "
                   "the 8-12 um window at a declared transmittance of "
                   f"{WINDOW_TRANSMITTANCE}, which scales it linearly.")
        ds.window_transmittance = WINDOW_TRANSMITTANCE
    return drss, drls, column
