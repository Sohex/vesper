#!/usr/bin/env python3
"""What the dust burden is worth radiatively, shortwave AND longwave.

    python exoplasim/scripts/dust_forcing.py
    python exoplasim/scripts/dust_forcing.py --indices bracket   # the absorbing end

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

**The refractive indices**, which are DECLARED rather than chosen here.
`aeolian/config/dust.yaml` says which dataset this world's dust is, per spectral
region, and this file and the aerofile ExoPlaSim reads are built from that one
declaration. Until DUST-12 they were not: this computed everything from OPAC
while the aerofile used the measured datasets, so the two priced the same burden
at absorption optical depths a factor of 2.5 apart. `--indices bracket` prices
the absorbing end the config also declares.

## What it writes

`analysis/dust_forcing.json`, the per-surface pricing, and
`analysis/dust_surface_forcing.nc`, the per-cell perturbation to `rss` and `rls`
that `carve_verdict.py --dust-forcing` reads. Both are declared against the
`dust_forcing` step in `config/pipeline.yaml`.

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
from dust_indices import check_coverage, indices, load_config, selection  # noqa: E402
from mie_dust import lognormal_integrate  # noqa: E402
from sra import read_sra  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "lib") not in sys.path:
    sys.path.insert(0, str(ROOT / "lib"))
from lapse import environmental_lapse_k_per_km  # noqa: E402
import climatology  # noqa: E402  from lib/, via _paths
from paths import climatology_path, rel  # noqa: E402
from provenance import staged_surface_field  # noqa: E402
from stellar import band1_fraction  # noqa: E402

OUT = ROOT / "analysis" / "dust_forcing.json"
OUT_NC = ROOT / "analysis" / "dust_surface_forcing.nc"

# Declared geometry, all of it stated rather than fitted.
#
# DIFFUSIVITY is the standard two-stream factor for a diffuse longwave field.
# The lapse rate and the dust scale height together set the temperature the
# layer radiates at, which is what the longwave term is proportional to; the
# scale height is `aeolian/config/dust.yaml`'s own, and for an exponential
# profile the mass-weighted mean height IS the scale height. The rate itself is
# NOT declared here: it was Earth's 6.5 K/km until PHYS-12, and it is now
# measured from the baseline climatology at run time by `lib/lapse.py`
# (season="annual", because t_surface below is an annual mean).
# TRANSMISSION is the clear-sky shortwave transmission of the atmosphere above
# the layer, which the two-stream forcing expression squares.
DIFFUSIVITY = 1.66
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
# The per-cell field can only be written at ONE end, and this is which. The
# emitted end is the smaller perturbation of the two and is the end DUST-10's
# lake result was computed through, so it is kept rather than quietly swapped.
# The other end is worth about 1.14x on the surface shortwave. What is NOT a
# bracket is the column mass: the chain reports a mass and an optical depth
# formed from it, so `aod / mee_chain` recovers the mass whichever end is then
# applied to it.
SIZE_END_KEYS = {"fine": "fine (Balkanski, what the reported AOD assumes)",
                 "emitted": "emitted (Kok, volume-median diameter 3.4 um)"}
RHO_G_CM3 = 2.6
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


def band_optics(lo_um, hi_um, weights_of_lam, dist, n_of, k_of,
                n_grid=24, n_r=200):
    """Weighted band-mean extinction, single-scattering albedo and asymmetry.

    `n_of` and `k_of` come from `dust_indices` for the dataset the config
    declares, so this function has no opinion about which dust it is describing.
    """
    grid = np.linspace(lo_um, hi_um, n_grid)
    w = weights_of_lam(grid)
    w = w / w.sum()
    ext = sca = gsc = 0.0
    for lam, ww in zip(grid, w):
        _, ssa, g, mee = lognormal_integrate(
            lam, n_of(lam), k_of(lam), dist["r_mod_um"], dist["sigma_g"],
            0.01, 25.0, RHO_G_CM3, n_r=n_r)
        ext += mee * ww
        sca += mee * ssa * ww
        gsc += g * mee * ssa * ww
    return ext, sca / ext, gsc / max(sca, 1e-30)


def shortwave_bands(sel, band1_share):
    """[(lo, hi), flux share, dataset name, n_of, k_of] for the two SW bands.

    One place the config's declared choice becomes something the integrations
    below can use, so a band and the indices it is integrated with cannot drift
    apart inside this file.
    """
    out = []
    for (lo, hi), share, name in zip(SW_BANDS_UM,
                                     (band1_share, 1.0 - band1_share),
                                     (sel["band1"], sel["band2"])):
        check_coverage(name, lo, hi)
        n_of, k_of = indices(name)
        out.append(((lo, hi), share, name, n_of, k_of))
    return out


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


def longwave_forcing(column_kg_m2, dist, t_surface, t_layer, lw_indices,
                     n_grid=40):
    """Greybody TOA longwave warming from an absorbing layer aloft.

    Spectrally resolved rather than band-averaged, because the dust absorption
    and the Planck function both vary strongly across the thermal infrared and
    their product is what matters:

        dF = integral [ B(T_s) - B(T_d) ] * (1 - exp(-D * tau_abs)) dlam

    CLEAR-SKY WINDOW APPROXIMATION, so this is an UPPER bound: part of the band
    is already opaque to water vapour and CO2, and dust absorbing where the
    atmosphere already absorbs adds nothing.
    """
    n_of, k_of = lw_indices
    grid = np.linspace(LW_LO_UM, LW_HI_UM, n_grid)
    total = 0.0
    tau_abs_mean = 0.0
    b_s, b_d = planck(grid, t_surface), planck(grid, t_layer)
    for i, lam in enumerate(grid):
        _, ssa, _, mee = lognormal_integrate(
            lam, n_of(lam), k_of(lam), dist["r_mod_um"], dist["sigma_g"],
            0.01, 25.0, RHO_G_CM3, n_r=200)
        tau_abs = mee * 1e3 * column_kg_m2 * (1.0 - ssa)   # mee is m2/g
        emissivity = 1.0 - np.exp(-DIFFUSIVITY * tau_abs)
        total += (b_s[i] - b_d[i]) * emissivity
        tau_abs_mean += tau_abs
    dlam = (LW_HI_UM - LW_LO_UM) / (n_grid - 1) * 1e-6
    return float(total * dlam), float(tau_abs_mean / n_grid)


def surface_forcing(column_kg_m2, dist, albedo, t_layer, insolation,
                    sw_bands, lw_indices, n_grid=24):
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
    for (lo, hi), share, _name, n_of, k_of in sw_bands:
        ext, ssa, g = band_optics(lo, hi, lambda x: np.ones_like(x), dist,
                                  n_of, k_of)
        tau = ext * 1e3 * column_kg_m2
        tau_sca += tau * ssa * share
        tau_abs += tau * (1.0 - ssa) * share
        beta_eff += backscatter_fraction(g) * share
    sw = -insolation * TRANSMISSION ** 2 * (1.0 - albedo) * (
        beta_eff * tau_sca + tau_abs)

    n_of, k_of = lw_indices
    grid = np.linspace(WINDOW_LO_UM, WINDOW_HI_UM, n_grid)
    b = planck(grid, t_layer)
    lw = 0.0
    for i, lam in enumerate(grid):
        _, ssa, _, mee = lognormal_integrate(
            lam, n_of(lam), k_of(lam), dist["r_mod_um"], dist["sigma_g"],
            0.01, 25.0, RHO_G_CM3, n_r=200)
        emissivity = 1.0 - np.exp(-DIFFUSIVITY * mee * 1e3 * column_kg_m2 * (1.0 - ssa))
        lw += b[i] * emissivity
    dlam = (WINDOW_HI_UM - WINDOW_LO_UM) / (n_grid - 1) * 1e-6
    return float(sw), float(lw * dlam * WINDOW_TRANSMITTANCE)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dust", type=Path,
                    default=ROOT / "aeolian" / "analysis" / "dust_baseline.json")
    ap.add_argument("--dust-field", type=Path,
                    default=ROOT / "aeolian" / "analysis" / "dust_baseline.nc",
                    help="per-cell optical depth, for the surface forcing field")
    ap.add_argument("--config", type=Path, default=ROOT / "config" / "planet.yaml")
    ap.add_argument("--indices", default="central", choices=("central", "bracket"),
                    help="which refractive-index set from aeolian/config/dust.yaml: "
                         "`central` is the declared choice, `bracket` the absorbing "
                         "OPAC end. DUST-12")
    ap.add_argument("--output", type=Path, default=OUT)
    ap.add_argument("--size-end", default="emitted", choices=("fine", "emitted"),
                    help="which end of the size bracket the per-cell field is "
                         "written at; the JSON always reports both")
    ap.add_argument("--output-nc", type=Path, default=OUT_NC,
                    help="per-cell surface forcing, which carve_verdict.py reads")
    args = ap.parse_args()

    ratio = _check_planck()
    print(f"planck normalisation checked: {ratio:.4f} of sigma T^4\n")
    import yaml
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    dust = json.loads(args.dust.read_text(encoding="utf-8"))
    dust_cfg = load_config()

    # WHICH DUST THIS IS. Declared once in aeolian/config/dust.yaml and read
    # here, so this file and the aerofile the model reads describe the same
    # particle. They did not until DUST-12: this computed everything from OPAC
    # while the aerofile was built from the measured datasets, and the absorption
    # optical depths differed by a factor of 2.5 for the same burden.
    sel = selection(dust_cfg, end=args.indices)
    check_coverage(sel["longwave"], LW_LO_UM, LW_HI_UM)
    lw_indices = indices(sel["longwave"])
    print(f"indices ({args.indices}): band 1 {sel['band1']}, "
          f"band 2 {sel['band2']}, thermal {sel['longwave']}")

    flux = float(config["orbit"]["baseline_flux_earth"])
    insolation = float(config["orbit"]["earth_solar_constant_w_m2"]) * flux / 4.0
    # Global-mean surface temperature from the climatology the config names,
    # bin-weighted. Was a hardcoded 289.80, a current value living in code --
    # the same defect family as the 6.5 K/km this file just shed, flagged in
    # the PHYS-12 pass. Planck emission goes as T^4, so drift here reaches the
    # longwave term directly.
    from netCDF4 import Dataset as _Dataset
    from numpy.polynomial.legendre import leggauss as _leggauss
    with _Dataset(climatology_path()) as _ds:
        _ts = climatology.annual_mean_of(_ds, "ts")
        _gw = _leggauss(_ts.shape[0])[1][::-1]
        t_surface = float((_ts * _gw[:, None]).sum() / (_gw.sum() * _ts.shape[1]))
    scale_h_m = float(dust_cfg["transport"]["dust_scale_height_m"])
    lapse_k_per_km = environmental_lapse_k_per_km(config)
    print(f"lapse rate {lapse_k_per_km:.2f} K/km, measured (lib/lapse.py, annual)")
    t_layer = t_surface - lapse_k_per_km * scale_h_m / 1000.0

    mee_chain = float(dust["mass_extinction_efficiency_m2_kg"])
    central = dust["shelter_bracket"]["central"]
    aod_land = float(central["land_mean_aod"])
    column = aod_land / mee_chain          # kg/m2, the size-independent quantity

    # Shape WITHIN each band, for the Mie integration. The band SPLIT is not
    # taken from this file: integrating it across 0.75 um gives 0.3777, which is
    # its own 0.34-14.01 um truncation rather than the star. `lib/stellar.py`
    # reproduces what `radmod.f90:solarini` does with the hi-res spectrum, which
    # is the weight the model itself applies.
    d = np.loadtxt(ROOT / "exoplasim" / "inputs" / "stellarspectra" / "k25v.dat",
                   skiprows=1)
    lam_s, f_s = d[:, 0], d[:, 1] * np.gradient(d[:, 0])
    band1_share = band1_fraction()
    sw_bands = shortwave_bands(sel, band1_share)

    print(f"burden: land-mean AOD {aod_land:.4f} at {mee_chain:.1f} m2/kg")
    print(f"     -> column mass {column * 1e3:.3f} g/m2")
    print(f"layer: {t_surface:.1f} K surface, {t_layer:.1f} K at "
          f"{scale_h_m / 1000:.1f} km, insolation {insolation:.1f} W/m2\n")

    results = []
    for name, dist in SIZE_ENDS.items():
        sw = {}
        tau_sw_total = 0.0
        for (lo, hi), share, iname, n_of, k_of in sw_bands:
            ext, ssa, g = band_optics(
                lo, hi, lambda x: np.interp(x, lam_s, f_s), dist, n_of, k_of)
            sw[f"{lo}-{hi} um"] = {
                "indices": iname,
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

        lw, tau_lw = longwave_forcing(column, dist, t_surface, t_layer, lw_indices)
        per_surface = {}
        for surface, albedo in SURFACES.items():
            sw_f = shortwave_forcing(tau_sw, ssa_eff, beta_eff, albedo, insolation)
            # AT THE SURFACE, which is a different quantity from the two above
            # and is the one Penman reads. Reported here rather than left to an
            # ad-hoc call, because notes/dust.md quotes it and the carve verdict
            # is decided through it.
            sw_s, lw_s = surface_forcing(column, dist, albedo, t_layer,
                                         insolation, sw_bands, lw_indices)
            per_surface[surface] = {
                "shortwave_w_m2": round(float(sw_f), 3),
                "longwave_w_m2": round(lw, 3),
                "net_w_m2": round(float(sw_f) + lw, 3),
                "surface_shortwave_w_m2": round(float(sw_s), 3),
                "surface_longwave_w_m2": round(float(lw_s), 3),
                "surface_net_w_m2": round(float(sw_s + lw_s), 3)}
        results.append({
            "size_distribution": name,
            "r_mod_um": round(dist["r_mod_um"], 4),
            "sigma_g": dist["sigma_g"],
            "shortwave_optical_depth": round(float(tau_sw), 4),
            "shortwave_bands": sw,
            "shortwave_absorption_optical_depth": round(
                float(tau_sw * (1.0 - ssa_eff)), 5),
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

    # -- the per-cell surface field, which is what the carve verdict reads ----
    #
    # Written here rather than by an ad-hoc call, because config/pipeline.yaml
    # declares this step as its generator and an artifact no step writes is one
    # nothing can invalidate. See SIZE_END_KEYS for which end of the size
    # bracket it is written at and why that is a choice rather than a default.
    surface_field = None
    if args.dust_field.is_file():
        model = config["model"]
        nlat, nlon = int(model["latitudes"]), int(model["longitudes"])
        # Through the one door. That path is keyed by the RUNG alone while
        # `surface_albedo` rewrites it per BUILD, so it could not say which
        # build's field it held, and the per-cell forcing this writes is what
        # carve_verdict.py divides by (1 - albedo). No cross-build read: the
        # forcing is formed on the build the config names. world-z7bu, rule 5.
        staged_albedo = staged_surface_field(174, config)
        albedo_path = ROOT / staged_albedo["path"]
        albedo = read_sra(albedo_path, nlat, nlon)
        end_name = SIZE_END_KEYS[args.size_end]
        drss, drls, col = write_per_cell(
            args.dust_field, albedo, args.output_nc, SIZE_ENDS[end_name],
            t_layer, insolation, sw_bands, lw_indices, mee_chain,
            meta={"indices_band1": sel["band1"], "indices_band2": sel["band2"],
                  "indices_longwave": sel["longwave"],
                  "indices_end": args.indices,
                  "size_distribution": end_name,
                  "dust_field": rel(args.dust_field),
                  "background_albedo": staged_albedo["path"],
                  "background_albedo_build": staged_albedo["build"],
                  "background_albedo_sha256": staged_albedo["sha256"],
                  "generated": datetime.now(timezone.utc).isoformat()})
        from netCDF4 import Dataset as _DS
        with _DS(args.dust_field) as _ds:
            _lat = np.asarray(_ds["lat"][:], dtype=float)
        wt = np.cos(np.deg2rad(_lat))[:, None] * np.ones((1, drss.shape[1]))
        surface_field = {
            "output": rel(args.output_nc),
            "size_distribution": end_name,
            "background_albedo_field": staged_albedo,
            "area_weighted_drss_w_m2": round(float(np.average(drss, weights=wt)), 4),
            "area_weighted_drls_w_m2": round(float(np.average(drls, weights=wt)), 4),
            "min_drss_w_m2": round(float(drss.min()), 3),
            "max_drls_w_m2": round(float(drls.max()), 3),
            "note": "Per-cell perturbations to rss and rls, formed with the "
                    "per-cell background albedo because carve_verdict.py "
                    "recovers the downward shortwave as rss / (1 - albedo).",
        }
        print(f"\nsurface field: area-weighted drss "
              f"{surface_field['area_weighted_drss_w_m2']:+.4f}, drls "
              f"{surface_field['area_weighted_drls_w_m2']:+.4f} W/m2 "
              f"-> {rel(args.output_nc)}")
    else:
        print(f"\n{rel(args.dust_field)} is absent; the per-cell surface "
              f"forcing field was NOT written.")

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
                     "lapse_k_per_km": round(lapse_k_per_km, 3),
                     "lapse_source": "lib/lapse.py, measured from the baseline "
                                     "climatology, season=annual",
                     "diffusivity": DIFFUSIVITY,
                     "shortwave_transmission": TRANSMISSION,
                     "mean_insolation_w_m2": round(insolation, 2)},
        "reopening_threshold_w_m2": 1.5,
        "indices": {"end": args.indices, **sel,
                    "declared_in": "aeolian/config/dust.yaml, optics.indices",
                    "note": "The same declaration the aerofile ExoPlaSim reads "
                            "is built from, so the offline forcing and the "
                            "in-model forcing describe one particle. DUST-12."},
        "global_mean": {
            "global_mean_aod": aod_global,
            "land_fraction": land_fraction,
            "note": "Area-weighted over ocean and a half-and-half land mix of "
                    "vegetated and playa, scaled by the global-to-land optical "
                    "depth ratio the transport produced. Crude, and the sign is "
                    "what matters rather than the third digit.",
            "net_w_m2_by_size_end": [round(float(g), 3) for g in global_mean]},
        "surface_field": surface_field,
        "results": results,
        "caveats": [
            "The longwave term is a clear-sky window estimate. Part of the band "
            "is already opaque to water vapour and CO2, so this OVERSTATES it.",
            "The size distribution is a bracket, not a measurement: the reported "
            "optical depth assumes a fine Balkanski distribution and the emission "
            "scheme produces a much coarser one. The atmospheric burden is "
            "between them and nothing here says where.",
            "The shortwave and the thermal infrared are not on the same "
            "dataset and cannot be: the measured indices stop at 2.45 um, so "
            "the thermal term is OPAC's, which runs 1 to 2x more absorbing than "
            "the measurements where the two overlap. The mixture biases the net "
            "warm and the surface shortwave reduction small. aeolian/config/dust.yaml "
            "says "
            "why, and --indices bracket prices the all-OPAC end.",
            "Two-stream, optically-thin forcing expressions. At an optical depth "
            "near 1 they are being used past where they are strictly valid.",
        ],
    }
    args.output.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {rel(args.output)}")


def write_per_cell(dust_nc, albedo, out_path, dist, t_layer, insolation,
                   sw_bands, lw_indices, mee_chain_m2_kg, variant="central",
                   meta=None):
    """Per-cell surface forcing, as perturbations to `rss` and `rls`.

    Written so the carve verdict can read them straight into its own Penman
    without re-deriving any optics. The optical properties depend only on the
    size distribution, so they are computed once and applied across the grid;
    only the longwave emissivity is nonlinear in column mass.

    `albedo` is the BACKGROUND SURFACE ALBEDO PER CELL, and it has to be, because
    of what the consumer does with the result. `carve_verdict.py` recovers the
    downward shortwave as `rss / (1 - land_albedo)` before re-absorbing it at
    water's albedo, so a perturbation added to `rss` is only recovered correctly
    if it was formed with that same per-cell albedo. A scalar here is a per-cell
    error of `(1 - a_scalar) / (1 - a_cell)`, which is 1.55x over playa fill if
    the scalar is water's.

    The column mass comes from the chain's own mass extinction efficiency, the
    same quantity the reported optical depth was formed with, so `column` is the
    burden the aeolian component actually produced rather than a re-derivation
    of it at one band's efficiency.
    """
    from netCDF4 import Dataset as DS
    with DS(dust_nc) as ds:
        aod = np.asarray(ds[f"aod_{variant}"][:], dtype=float)
        lat = np.asarray(ds["lat"][:]); lon = np.asarray(ds["lon"][:])
    albedo = np.asarray(albedo, dtype=float)
    if albedo.shape != aod.shape:
        raise SystemExit(
            f"albedo is {albedo.shape} and the dust field is {aod.shape}. They "
            f"are the same grid index for index; do not reconcile them.")
    tau_sca_c = tau_abs_c = beta_eff = 0.0
    for (lo, hi), share, _name, n_of, k_of in sw_bands:
        ext, ssa, g = band_optics(lo, hi, lambda x: np.ones_like(x), dist,
                                  n_of, k_of)
        tau_sca_c += ext * 1e3 * ssa * share
        tau_abs_c += ext * 1e3 * (1.0 - ssa) * share
        beta_eff += backscatter_fraction(g) * share
    column = aod / mee_chain_m2_kg                # kg/m2 per cell

    drss = -insolation * TRANSMISSION ** 2 * (1.0 - albedo) * (
        beta_eff * tau_sca_c + tau_abs_c) * column

    n_of, k_of = lw_indices
    grid = np.linspace(WINDOW_LO_UM, WINDOW_HI_UM, 24)
    b = planck(grid, t_layer)
    drls = np.zeros_like(column)
    for i, lam in enumerate(grid):
        _, ssa, _, m = lognormal_integrate(lam, n_of(lam), k_of(lam),
                                           dist["r_mod_um"], dist["sigma_g"],
                                           0.01, 25.0, RHO_G_CM3, n_r=200)
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
                   f"{WINDOW_TRANSMITTANCE}, which scales it linearly. Formed "
                   "with the per-cell background albedo, which is the convention "
                   "carve_verdict.py inverts.")
        ds.window_transmittance = WINDOW_TRANSMITTANCE
        for key, value in (meta or {}).items():
            setattr(ds, key, value)
    return drss, drls, column


if __name__ == "__main__":
    main()
