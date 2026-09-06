#!/usr/bin/env python3
"""Emit, transport and deposit sea salt, offline, from a climatology.

    python aeolian/scripts/build_sea_salt.py
    python aeolian/scripts/build_sea_salt.py --bracket    # every declared end

CLIM-27. `notes/audits/unpriced-terms.md` finding 2 says why this exists: this
project priced mineral dust carefully and had never asked what any other aerosol
was worth, on a planet whose ocean is more than half the surface. The answer has
to be a number, because the argument that a scattering aerosol over a dark ocean
matters is available for free and settles nothing.

It is `build_dust.py`'s counterpart and deliberately its sibling rather than its
generalisation: the two share `settling_velocity` and `advect_to_steady_state`
and nothing else, because everything upstream of transport differs.

## What differs from dust, and why each difference is not cosmetic

**The source is the wind and only the wind.** No lithology, no soil moisture, no
threshold, no drag partition, no vegetation. Whitecaps break where the wind
blows, so the erodibility machinery that is most of `build_dust.py` has no
counterpart here. What survives is the subgrid wind distribution, and it
survives in a better form: production goes as U10^3.5 with no threshold, so the
subgrid moment is analytic rather than quadrature.

**The particle is wet.** It settles as a solution droplet, twice the dry
diameter at 80% relative humidity and half the density, and it extinguishes as
one. Only the mass budget is dry. `sea_salt_optics.py` carries the optics and
the check that OPAC's growth table has been read correctly.

**Gravity enters twice, in opposite directions.** Settling is 1.31 times faster
at 12.81 m/s2, which shortens the lifetime and lowers the burden. And the
Charnock roughness of the sea surface is 0.766 times Earth's for the same
friction velocity, which steepens the wind profile and lowers U10 slightly for
the same level wind, which lowers emission. Both push the same way.

**Sea ice suppresses it.** There is no whitecap under ice, and this is the term
that makes this world's cold branch differ from its warm one in sea salt.

## What it does not do

No feedback: the climatology is an input. No coagulation, no in-cloud
processing, no source dependence on salinity -- `config/planet.yaml` declares a
salinity but the source function used here carries no salinity weight, and
Grythe et al. find the one function that does is not improved by it.

Every constant is declared in `aeolian/config/sea_salt.yaml` with its source,
and anything that could not be taken from a source carries a bracket this script
reports rather than collapsing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset

from _paths import ANALYSIS, CONFIG, PROJECT_ROOT  # noqa: E402

import climatology  # noqa: E402  from lib/, via _paths
from build_dust import advect_to_steady_state, flag_anomalous_bins, settling_velocity
from build_dust import steering_wind, weibull_shape_from_samples
from gridding import gaussian_area_weights
from paths import best_available_climatology, rel, snapshot_beside

sys.path.insert(0, str(PROJECT_ROOT / "exoplasim" / "scripts"))
# The two-stream expression and the upscatter fraction are imported, not
# copied, so mineral dust and sea salt are priced through one formula.
from dust_forcing import backscatter_fraction, shortwave_forcing  # noqa: E402
from sea_salt_source import mass_flux, moment_ratio
from aerosol_deposition import write_deposition

SEA_SALT_CONFIG = PROJECT_ROOT / "aeolian" / "config" / "sea_salt.yaml"
OPTICS = PROJECT_ROOT / "analysis" / "sea_salt_optics.json"
OUT_JSON = ANALYSIS / "sea_salt_baseline.json"
OUT_NC = ANALYSIS / "sea_salt_baseline.nc"
# The nutrient carrier, and it is a SEPARATE FILE on purpose. See
# aerosol_deposition.py: the abiotic nutrient ledger's screen refuses a carrier
# whose path names an optics product, and `sea_salt_baseline.nc` holds optical
# depth. A file that is half optics cannot be the mass carrier.
OUT_DEP_JSON = ANALYSIS / "sea_salt_deposition.json"
OUT_DEP_NC = ANALYSIS / "sea_salt_deposition.nc"

R_DRY = 287.05
EARTH_YEAR_S = 365.25 * 86400.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_commit() -> str | None:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT,
                              capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:
        return None


def u10_from_level(u_level, z_level_m, cfg, gravity):
    """10 m wind from the lowest model level, over a Charnock sea surface.

    Iterated because the roughness depends on the friction velocity it sets:
    z0 = 0.11 nu / u* + a u*^2 / g (Smith 1988). Gravity is this world's, which
    is the only place it reaches the emission other than through settling.

    Returns U10 and the roughness, so the caller can report what it used.
    """
    s = cfg["surface_layer"]
    kappa, a, nu = s["von_karman"], s["charnock"], s["kinematic_viscosity_m2_s"]
    z0 = np.full_like(u_level, 1.5e-4)
    for _ in range(int(s["iterations"])):
        ustar = kappa * u_level / np.log(np.maximum(z_level_m / z0, 1.001))
        z0 = np.maximum(s["smooth_flow_coefficient"] * nu
                        / np.maximum(ustar, 1e-6) + a * ustar ** 2 / gravity,
                        1e-7)
    ustar = kappa * u_level / np.log(np.maximum(z_level_m / z0, 1.001))
    return (ustar / kappa) * np.log(s["reference_height_m"] / z0), z0


def optics_table(path: Path) -> dict:
    """{(bin_lo, bin_hi): (rh[], mee1[], mee2[], ssa1[], g1[], gf[])}."""
    payload = json.loads(path.read_text())
    out: dict[tuple[float, float], dict[str, list]] = {}
    for row in payload["results"]:
        key = tuple(row["bin_dry_um"])
        d = out.setdefault(key, {"rh": [], "mee1": [], "mee2": [], "ssa1": [],
                                 "g1": [], "gf": []})
        d["rh"].append(row["relative_humidity_percent"])
        d["mee1"].append(row["band1"]["mass_extinction_efficiency_m2_g_dry"])
        d["mee2"].append(row["band2"]["mass_extinction_efficiency_m2_g_dry"])
        d["ssa1"].append(row["band1"]["single_scattering_albedo"])
        d["g1"].append(row["band1"]["asymmetry_parameter"])
        d["gf"].append(row["growth_factor"])
    for d in out.values():
        order = np.argsort(d["rh"])
        for k in d:
            d[k] = np.asarray(d[k], dtype=float)[order]
    return out, payload


def annual_mean(field, weights):
    return np.tensordot(weights, field, axes=(0, 0))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--climatology", type=Path, default=None,
                    help="defaults to the configured bootstrap_climatology, "
                         "which is what the sea_salt step declares it needs. "
                         "It moves the wind tail with it")
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--sea-salt-config", type=Path, default=SEA_SALT_CONFIG)
    ap.add_argument("--optics", type=Path, default=OPTICS)
    ap.add_argument("--gust-samples", type=Path, default=None,
                    help="high-cadence wind extract; the Weibull shape over "
                         "OCEAN cells is fitted from it")
    ap.add_argument("--bracket", action="store_true",
                    help="run every declared bracket end and report the spread")
    ap.add_argument("--output", type=Path, default=OUT_JSON)
    ap.add_argument("--output-nc", type=Path, default=OUT_NC)
    ap.add_argument("--output-deposition", type=Path, default=OUT_DEP_JSON)
    ap.add_argument("--output-deposition-nc", type=Path, default=OUT_DEP_NC)
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text())
    cfg = yaml.safe_load(args.sea_salt_config.read_text())
    # `advect_to_steady_state` is shared with build_dust.py and reads the
    # planet radius off the config dict it is handed, because the CFL step and
    # the polar cos(lat) floor are both lengths on this world's sphere. It is
    # injected the same way build_dust.py injects it. Without it the transport
    # raises KeyError before producing anything, which is what it did.
    cfg["_planet_radius_earth"] = config["planet"]["radius_earth"]
    gravity = float(config["planet"]["gravity_m_s2"])
    # THE BEST AVAILABLE, which is the baseline once one is named and the
    # bootstrap before that. The source function is the near-surface wind to a
    # high power over open water, gated by sea ice, and the burden is grown at
    # the modelled humidity and washed out by the modelled precipitation: every
    # one of those is the climate STATE. The `needs: bootstrap_climatology`
    # edge in config/pipeline.yaml is unchanged and states what must EXIST for
    # the first pass to run, not what a later pass should read.
    clim_path, clim_stage = best_available_climatology(args.climatology)
    if not args.optics.is_file():
        raise SystemExit(f"{args.optics} missing; run "
                         f"aeolian/scripts/sea_salt_optics.py first")
    optics, optics_payload = optics_table(args.optics)

    with Dataset(clim_path) as ds:
        # The bin centres, not a bin index. `lib/climatology.py` recovers the
        # record count each bin holds from the SPACING of these, and evenly
        # spaced integers say "the count divides the bin count" whatever the
        # file was written from.
        bin_centres = np.asarray(ds["time"][:], dtype=float)
        lat = np.asarray(ds["lat"][:], dtype=float)
        lon = np.asarray(ds["lon"][:], dtype=float)
        lev = np.asarray(ds["lev"][:], dtype=float)
        spd_all = np.asarray(ds["spd"][:], dtype=float)
        ua = np.asarray(ds["ua"][:], dtype=float)
        va = np.asarray(ds["va"][:], dtype=float)
        ta = np.asarray(ds["ta"][:], dtype=float)[:, -1]
        hur = np.asarray(ds["hur"][:], dtype=float)[:, -1]
        ps = np.asarray(ds["ps"][:], dtype=float)
        ts = np.asarray(ds["ts"][:], dtype=float)
        sic = np.asarray(ds["sic"][:], dtype=float)
        pr = np.asarray(ds["pr"][:], dtype=float)
        rst = np.asarray(ds["rst"][:], dtype=float)
        rsut = np.asarray(ds["rsut"][:], dtype=float)
        alb1 = np.asarray(ds["alb1"][:], dtype=float)
        alb2 = np.asarray(ds["alb2"][:], dtype=float)
        lsm = np.asarray(ds["lsm"][:], dtype=float)[0]

    # The same corrupted first output bin `build_dust.py` found, and for the
    # same reason: this reads `spd` from the same regular climatology.
    bad = flag_anomalous_bins(spd_all)
    spd = spd_all[:, -1]
    weights = np.asarray(climatology.bin_weights(bin_centres), dtype=float)
    if bad:
        weights[bad] = 0.0
    weights = weights / weights.sum()

    ocean = lsm < 0.5
    # THE AREA WEIGHT and THE SOLVER'S METRIC are two different quantities and
    # this is the one place they part. `area` is each cell's share of the
    # sphere, the Gauss-Legendre partition the model's own budget is taken over,
    # and every mean and every emission total below is weighted by it.
    # `coslat_solver` is the cosine `advect_to_steady_state` divides its
    # east-west width by, floored near the poles so the explicit step survives
    # the converging longitude grid; the mass balance is checked against that
    # because it is the metric the solver actually conserved over.
    area = gaussian_area_weights(lat, lon.size, what=str(clim_path))
    coslat_solver = np.maximum(
        np.cos(np.deg2rad(lat))[:, None] * np.ones((1, lon.size)),
        cfg["transport"].get("polar_coslat_floor", 1e-3))
    area_w = area * ocean

    spd_a = annual_mean(spd, weights)
    ta_a = annual_mean(ta, weights)
    ps_a = annual_mean(ps, weights) * 100.0                 # hPa -> Pa
    ts_a = annual_mean(ts, weights)
    sic_a = np.clip(annual_mean(sic, weights), 0.0, 1.0)
    pr_a = annual_mean(pr, weights)
    # Incident TOA shortwave, from the model's own net and outgoing: the
    # climatology carries no downward flux, and rsdt = rst + rsut by definition.
    rsdt_a = annual_mean(rst, weights) + annual_mean(rsut, weights)
    alb1_a = np.clip(annual_mean(alb1, weights), 0.0, 1.0)
    alb2_a = np.clip(annual_mean(alb2, weights), 0.0, 1.0)
    # `hur` IS IN PERCENT, and code 157 now declares it. Upstream ExoPlaSim's
    # pyburn table labelled it "1" while PlaSim writes a percentage, so a reader
    # that trusted the attribute and multiplied by 100 got 100% everywhere,
    # which through the growth curve is a factor of five on the optical depth;
    # WORLD-HF12 fixed the attribute rather than the values. A climatology
    # postprocessed before that fix still carries the "1" attribute over
    # percentage values, so the range is checked below rather than assumed.
    hur_a = np.clip(annual_mean(hur, weights), 0.0, 100.0)
    if not 1.0 < float(np.nanmax(hur_a)) <= 105.0:
        raise SystemExit(
            f"relative humidity peaks at {float(np.nanmax(hur_a)):.3f}, which "
            f"is neither a percentage nor this reader's assumption; check the "
            f"climatology's `hur` before trusting anything downstream")
    rho_a = ps_a * lev[-1] / (R_DRY * ta_a)

    # Reference height of the lowest level, from the model's own sigma through
    # the hypsometric relation. Not an assumed 10 m: it is near 109 m here.
    z_level = (R_DRY * ta_a / gravity) * np.log(1.0 / lev[-1])
    sst_c = ts_a - 273.15

    def wmean(f, w):
        return float((f * w).sum() / w.sum())

    def solve(cfg):
        """Everything downstream of the climatology, for one configuration.

        Taken as a function of the config so that `--bracket` can re-run it at
        each declared bracket end rather than scaling the central answer, which
        would assume a linearity none of these terms has: the wet lifetime acts
        through a steady state, the scale height through the settling loss, and
        the Charnock coefficient through a logarithm.
        """
        u10, z0_sea = u10_from_level(spd_a, z_level, cfg, gravity)

        # Weibull shape over OCEAN, fitted the way dust fits it over erodible land.
        # Beside the regular product this run resolved, so the wind tail and
        # the binned fields cannot come from two different climatologies and
        # `--climatology` moves both.
        samples = args.gust_samples or snapshot_beside(clim_path)
        fit = weibull_shape_from_samples(samples, ocean)
        k_fit, n_samples = fit if fit else (None, 0)
        k = k_fit if k_fit else cfg["subgrid_wind"]["weibull_shape"]

        edges = cfg["size"]["bin_edges_um"]
        _, _, clamped = mass_flux(u10, sst_c, cfg, weibull_k=k)
        gate = ocean * (1.0 - sic_a if cfg["emission"]["sea_ice_suppression"]
                        else 1.0)

        rh_cell = np.clip(hur_a, cfg["growth"]["rh_min_percent"],
                          cfg["growth"]["rh_max_percent"])

        # Steering wind: LAYER-MASS weighted over the layers at or below the
        # steering sigma, which for a boundary-layer aerosol is most of what
        # carries it. The weight is `dsigma` renormalised over those layers,
        # from `build_dust.steering_weights`. The two layers it selects here
        # are 0.087 and 0.034 thick, so a plain mean over them is a different
        # quantity and puts the weight 0.22 wrong on each.
        u_lev, v_lev = steering_wind(ua, va, lev, cfg, config)
        u_steer = annual_mean(u_lev, weights)
        v_steer = annual_mean(v_lev, weights)

        # Wet removal, anchored to Jaegle et al. (2011)'s measured accumulation-mode
        # lifetime at Earth's ocean-mean precipitation and carried to this world's
        # own precipitation with Sportisse's exponent. Size-independent, which is
        # the point: in-cloud scavenging removes whatever activated, and sea salt
        # activates at every size.
        scav = cfg["removal"]
        p_mm_h = np.maximum(pr_a, 0.0) * 1000.0 * 3600.0
        lam_ref = 1.0 / (scav["wet_lifetime_hours"] * 3600.0)
        lam_wet = lam_ref * (p_mm_h
                             / scav["reference_precipitation_mm_per_hour"]
                             ) ** scav["scavenging_b"]

        planet_area = 4.0 * np.pi * (
            float(config["planet"]["radius_earth"]) * 6.371e6) ** 2
        b1 = optics_payload["stellar_flux_fraction_band1"]
        # Optical properties for the forcing, taken at the ocean-mean humidity and
        # weighted across bins by their optical depth rather than their mass, which
        # is what a single-layer two-stream expression is asking for.
        rh_mean = float((rh_cell * area_w).sum() / area_w.sum())
        ssa_band1 = float(np.mean([np.interp(rh_mean, optics[kk]["rh"],
                                             optics[kk]["ssa1"])
                                   for kk in optics]))
        g_band1 = float(np.mean([np.interp(rh_mean, optics[kk]["rh"],
                                           optics[kk]["g1"]) for kk in optics]))
        ssa_band2, g_band2 = ssa_band1, g_band1
        beta_band1 = backscatter_fraction(g_band1)
        beta_band2 = backscatter_fraction(g_band2)

        height = cfg["transport"]["scale_height_m"]

        def run(per_bin_emission):
            """Transport every size bin to steady state and accumulate the column.

            Returns the burden, the two band optical depths, a per-bin report and
            a mass residual, which is a DIAGNOSTIC and not an identity.

            At steady state the area integral of emission would equal the area
            integral of loss times burden if advection only moved mass. The shared
            `advect_to_steady_state` is written in ADVECTIVE form -- it steps
            `u dm/dx` rather than `d(um)/dx` -- so it conserves mass only where the
            steering wind is non-divergent, and a horizontal wind on a sigma
            surface is not. The residual therefore measures the divergence of the
            steering field, and a few percent is expected. What it still catches,
            loudly, is a bin that has not relaxed or a sign error in the loss term,
            which show up as tens of percent.
            """
            burden = np.zeros_like(u10)
            aod1 = np.zeros_like(u10)
            aod2 = np.zeros_like(u10)
            emitted = removed = 0.0
            report = []
            # The two removal terms, kept apart and per bin, because that is the
            # carrier the abiotic nutrient ledger needs and it is a PARTITION of
            # a rate the solver has already balanced against emission rather
            # than a second calculation beside it. ANUT-7.
            dry_dep, wet_dep = [], []
            for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
                tab = optics[(lo, hi)]
                gf = np.interp(rh_cell, tab["rh"], tab["gf"])
                # Wet density from the same table the growth came from: water plus
                # the salt volume fraction the growth implies.
                rho_p = (1000.0 + (cfg["emission"]["dry_density_kg_m3"] - 1000.0)
                         / gf ** 3)
                d_dry = np.sqrt(lo * hi)
                v_s = settling_velocity(d_dry * gf, rho_p, rho_a, ta_a, gravity)
                loss = v_s / height + lam_wet
                m, steps, converged = advect_to_steady_state(
                    per_bin_emission[i], u_steer, v_steer, loss, lat, lon, cfg)
                if not converged:
                    print(f"  WARNING bin {lo}-{hi} um did not converge in "
                          f"{steps} steps")
                dry_dep.append(v_s / height * m)
                wet_dep.append(lam_wet * m)
                burden = burden + m
                # Weighted by the metric the SOLVER uses, not by the true one.
                # `advect_to_steady_state` floors cos(lat) near the poles to keep
                # the explicit step from collapsing on the converging longitude
                # grid, so the cell areas it conserves mass over are the floored
                # ones. Checking against true areas measures that declared
                # distortion instead of the solver, and reads as a 2% leak.
                emitted += float((per_bin_emission[i] * coslat_solver).sum())
                removed += float((loss * m * coslat_solver).sum())
                d_aod1 = m * 1000.0 * np.interp(rh_cell, tab["rh"], tab["mee1"])
                d_aod2 = m * 1000.0 * np.interp(rh_cell, tab["rh"], tab["mee2"])
                aod1 = aod1 + d_aod1
                aod2 = aod2 + d_aod2
                report.append({
                    "bin_dry_um": [lo, hi],
                    "emission_tg_per_earth_year": round(float(
                        (per_bin_emission[i] * area).sum() / area.sum()
                        * planet_area * EARTH_YEAR_S / 1e9), 1),
                    "burden_mg_m2_ocean_mean": round(float(
                        (m * area_w).sum() / area_w.sum() * 1e6), 4),
                    "optical_depth_ocean_mean": round(float(
                        ((b1 * d_aod1 + (1 - b1) * d_aod2) * area_w).sum()
                        / area_w.sum()), 5),
                    "mean_wet_diameter_um": round(float(
                        (d_dry * gf * area_w).sum() / area_w.sum()), 3),
                    "mean_settling_velocity_mm_s": round(float(
                        (v_s * area_w).sum() / area_w.sum() * 1000.0), 3),
                    "mean_lifetime_days": round(float(
                        (1.0 / loss * area_w).sum() / area_w.sum() / 86400.0), 3),
                    "transport_steps": steps, "converged": bool(converged),
                })
            residual = abs(emitted - removed) / max(emitted, 1e-30)
            # IDENTITY, and it can fail. The deposition fields are a PARTITION
            # of the same `loss * m` the mass balance above already summed, so
            # summing them back over the same metric has to reproduce `removed`
            # to rounding. A sign error, a dropped bin or a scale height applied
            # to the wrong half all break it, and none of them would show in the
            # mass residual, which is a diagnostic of the steering field rather
            # than of this split.
            dry_a, wet_a = np.array(dry_dep), np.array(wet_dep)
            split_sum = float(((dry_a + wet_a).sum(axis=0)
                               * coslat_solver).sum())
            if abs(split_sum / max(removed, 1e-30) - 1.0) > 1e-10:
                raise SystemExit(
                    f"the wet/dry deposition split sums to {split_sum:.6e} "
                    f"against a removal of {removed:.6e}. It is meant to be the "
                    f"same quantity partitioned, so this is an implementation "
                    f"error and not a physical residual.")
            return (burden, aod1, aod2, report, residual, dry_a, wet_a)


        # TWO VARIANTS, AND THE SECOND IS THE ONE ANY EARTH NUMBER IS COMPARABLE
        # WITH. Grythe's third mode is centred at Dp = 30 um and its lower tail
        # dominates the mass below the 10 um cut; those drops live about five hours
        # here and are not what a PM10 network or an optical depth retrieval sees.
        # `all_modes` is equation 7 as written. `no_spume` drops the third mode and
        # is what reproduces Grythe's own stated production. Neither is a
        # correction to the other and both are reported.
        variants = {}
        deposition = {}
        for label, modes in (("all_modes", None), ("no_spume", (0, 1))):
            _, per_bin_v, _ = mass_flux(u10, sst_c, cfg, weibull_k=k,
                                        bin_edges_um=edges, modes=modes)
            per_bin_v = per_bin_v * gate
            (burden_v, aod1_v, aod2_v, report_v, residual_v,
             dry_v, wet_v) = run(per_bin_v)
            if residual_v >= 0.10:
                raise SystemExit(
                    f"{label}: mass residual {residual_v*100:.2f}% is far beyond "
                    f"what steering-wind divergence explains; check the per-bin "
                    f"convergence flags and the loss term")
            aod_v = b1 * aod1_v + (1.0 - b1) * aod2_v
            # Chylek and Coakley two-stream, the same expression dust_forcing.py
            # prices mineral dust with, so the two aerosols are compared through
            # one formula rather than two. Sea salt's single-scattering albedo is 1
            # to within 1e-5 in both bands, so the absorbing term is identically
            # zero and the sign cannot come out positive over any surface: this is
            # a scatterer, and over an ocean albedo of 0.07 it can only cool.
            f1 = shortwave_forcing(aod1_v, ssa_band1, beta_band1, alb1_a,
                                   b1 * rsdt_a)
            f2 = shortwave_forcing(aod2_v, ssa_band2, beta_band2, alb2_a,
                                   (1.0 - b1) * rsdt_a)
            forcing = f1 + f2
            variants[label] = {
                "toa_shortwave_forcing_w_m2_global": round(
                    wmean(forcing, area), 4),
                "toa_shortwave_forcing_w_m2_ocean": round(
                    wmean(forcing, area_w), 4),
                "emission_tg_per_earth_year": round(float(
                    (per_bin_v.sum(axis=0) * area).sum() / area.sum()
                    * planet_area * EARTH_YEAR_S / 1e9), 1),
                "burden_mg_m2_ocean_mean": round(wmean(burden_v, area_w) * 1e6, 4),
                "burden_mg_m2_global_mean": round(wmean(burden_v, area) * 1e6, 4),
                "optical_depth_global_mean": round(wmean(aod_v, area), 5),
                "optical_depth_ocean_mean": round(wmean(aod_v, area_w), 5),
                "optical_depth_band1_ocean_mean": round(wmean(aod1_v, area_w), 5),
                "optical_depth_band2_ocean_mean": round(wmean(aod2_v, area_w), 5),
                "conservation_residual": round(residual_v, 6),
                "bins": report_v,
            }
            deposition[label] = (dry_v, wet_v)
            if label == "all_modes":
                burden, aod = burden_v, aod_v
                emission_field = per_bin_v.sum(axis=0)

        return {"variants": variants, "u10": u10, "z0": z0_sea,
                "rh": rh_cell, "clamped": clamped, "burden": burden,
                "aod": aod, "emission_field": emission_field,
                "deposition": deposition,
                "k": k, "k_fit": k_fit, "n_samples": n_samples}

    r = solve(cfg)
    variants = r["variants"]
    u10, z0_sea, rh_cell, clamped = r["u10"], r["z0"], r["rh"], r["clamped"]
    burden, aod, emission_field = r["burden"], r["aod"], r["emission_field"]
    deposition = r["deposition"]
    k, k_fit, n_samples = r["k"], r["k_fit"], r["n_samples"]

    # THE BRACKET, re-solved rather than scaled. Each end is a full solve
    # because none of these terms acts linearly on the answer: the wet lifetime
    # acts through a steady state, the scale height through the settling loss,
    # and the Charnock coefficient through a logarithm. The spread is reported
    # and never collapsed, as `dust.yaml` requires of its own constants.
    bracket = {}
    if args.bracket:
        import copy
        ends = [("removal", "wet_lifetime_hours", "wet_lifetime_hours_bracket"),
                ("transport", "scale_height_m", "scale_height_bracket_m"),
                ("surface_layer", "charnock", "charnock_bracket")]
        if not k_fit:
            ends.append(("subgrid_wind", "weibull_shape",
                         "weibull_shape_bracket"))
        for block, key, bkey in ends:
            for value in cfg[block][bkey]:
                if value == cfg[block][key]:
                    continue
                alt = copy.deepcopy(cfg)
                alt[block][key] = value
                v = solve(alt)["variants"]
                bracket[f"{block}.{key}={value}"] = {
                    label: {
                        "optical_depth_global_mean":
                            d["optical_depth_global_mean"],
                        "toa_shortwave_forcing_w_m2_global":
                            d["toa_shortwave_forcing_w_m2_global"],
                    } for label, d in v.items()}
        if k_fit:
            bracket["note"] = (
                "the Weibull shape is FITTED from instantaneous samples over "
                "ocean, so its declared bracket is not run; refit it rather "
                "than bracketing it")

    print(f"climatology       {rel(clim_path)}")
    if bad:
        print(f"  excluded time bins {bad} as anomalous in `spd`")
    print(f"Weibull shape over ocean  k = {k:.3f}"
          + (f" fitted from {n_samples} samples" if k_fit
             else " (config fallback)")
          + f", subgrid moment x{moment_ratio(3.5, k):.3f} at p = 3.5")
    print(f"ocean-mean U10 {wmean(u10, area_w):.3f} m/s from a level wind of "
          f"{wmean(spd_a, area_w):.3f} m/s at {wmean(z_level, area_w):.1f} m, "
          f"Charnock z0 {wmean(z0_sea, area_w):.2e} m")
    print(f"ocean-mean SST {wmean(sst_c, area_w):.2f} C, "
          f"{clamped*100:.1f}% of cells outside the temperature fit range; "
          f"relative humidity {wmean(rh_cell, area_w):.1f}%\n")
    print(f"{'variant':>10}{'emission Tg/yr':>16}{'burden mg/m2':>14}"
          f"{'AOD ocean':>11}{'AOD global':>12}{'mass resid':>12}"
          f"{'W/m2 TOA':>10}")
    for label, v in variants.items():
        print(f"{label:>10}{v['emission_tg_per_earth_year']:>16,.0f}"
              f"{v['burden_mg_m2_ocean_mean']:>14.2f}"
              f"{v['optical_depth_ocean_mean']:>11.4f}"
              f"{v['optical_depth_global_mean']:>12.4f}"
              f"{v['conservation_residual']*100:>11.3f}%"
              f"{v['toa_shortwave_forcing_w_m2_global']:>10.3f}")
    for name, v in bracket.items():
        if name == "note":
            continue
        print(f"  bracket {name:38s} "
              f"AOD {v['no_spume']['optical_depth_global_mean']:.4f} / "
              f"{v['all_modes']['optical_depth_global_mean']:.4f}   "
              f"W/m2 {v['no_spume']['toa_shortwave_forcing_w_m2_global']:+.3f} / "
              f"{v['all_modes']['toa_shortwave_forcing_w_m2_global']:+.3f}")

    payload = {
        "note": "Offline sea-salt emission, transport and optical depth. "
                "Generated by aeolian/scripts/build_sea_salt.py; do not edit.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(),
        "climatology": str(rel(clim_path)),
        # WHICH STAGE this burden was emitted from. lib/paths.py.
        "climatology_stage": clim_stage,
        "excluded_time_bins": bad,
        "gravity_m_s2": gravity,
        "subgrid_wind": {
            "weibull_shape": round(float(k), 4),
            "fitted": bool(k_fit), "samples": n_samples,
            "moment_ratio_p3.5": round(moment_ratio(3.5, k), 4),
        },
        "surface_layer": {
            "ocean_mean_level_wind_m_s": round(wmean(spd_a, area_w), 4),
            "level_height_m": round(wmean(z_level, area_w), 2),
            "ocean_mean_u10_m_s": round(wmean(u10, area_w), 4),
            "ocean_mean_charnock_z0_m": float(f"{wmean(z0_sea, area_w):.4e}"),
        },
        "forcing_note": "Emission, burden and optical depth are reported for "
                        "two variants. `all_modes` is Grythe equation 7 as "
                        "written; `no_spume` drops the 30 um mode, whose lower "
                        "tail dominates the mass below the 10 um cut and is "
                        "what Grythe's own reported production leaves out. Any "
                        "comparison with an Earth number is against no_spume.",
        "conditions": {
            "ocean_mean_sst_c": round(wmean(sst_c, area_w), 3),
            "fraction_outside_temperature_fit": round(clamped, 4),
            "ocean_mean_sea_ice_fraction": round(wmean(sic_a, area_w), 4),
            "ocean_mean_relative_humidity_percent": round(
                wmean(rh_cell, area_w), 2),
        },
        "variants": variants,
        "bracket": bracket,
        "inputs": {p.name: sha256(p) for p in
                   (args.sea_salt_config, args.optics, args.config)},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    with Dataset(args.output_nc, "w") as out:
        out.createDimension("lat", lat.size)
        out.createDimension("lon", lon.size)
        for name, data, units in (
                ("lat", lat, "deg"), ("lon", lon, "deg")):
            v = out.createVariable(name, "f8", (name,))
            v.units = units
            v[:] = data
        for name, data, units, long_name in (
                ("burden", burden, "kg m-2", "sea salt column burden, dry mass"),
                ("aod", aod, "1", "sea salt optical depth, band combined"),
                ("emission", emission_field, "kg m-2 s-1",
                 "sea salt emission, dry mass, all modes"),
                ("u10", u10, "m s-1", "10 m wind over the Charnock sea surface")):
            v = out.createVariable(name, "f8", ("lat", "lon"), zlib=True)
            v.units = units
            v.long_name = long_name
            v[:] = data
        out.note = payload["note"]
        out.generated = payload["generated"]
        out.climatology = payload["climatology"]
        out.climatology_stage = payload["climatology_stage"]
    print(f"\nwrote {rel(args.output)} and {rel(args.output_nc)}")

    # -- the nutrient carrier, a SEPARATE file holding mass and only mass -----
    #
    # ANUT-7 retains marine aerosol as a base-cation and sulfur source to the
    # abiotic nutrient ledger and had no carrier, because what this component
    # produced for it was optics. `biosphere/config/abiotic_nutrients.yaml`
    # refuses an optical quantity as a carrier by name and the refusal is
    # enforced, so the mass cannot ride in `sea_salt_baseline.nc` beside the
    # optical depth: it needs a file whose whole content is mass.
    #
    # ALL MODES, not `no_spume`. `no_spume` exists so that a total can be
    # compared with Grythe's own reported production, which cuts the same mode;
    # what lands on the ground is what is emitted, and the spume mode's fate is
    # settling rather than truncation. The JSON reports the land mean under both
    # so the difference is visible rather than assumed, and over land it is
    # small because a 30 um drop does not reach the coast.
    dry_all, wet_all = deposition["all_modes"]
    dry_ns, wet_ns = deposition["no_spume"]
    land_w = area * (~ocean)
    comp = cfg["composition"]
    dep = write_deposition(
        args.output_deposition_nc, args.output_deposition,
        lat=lat, lon=lon,
        bins_um=list(zip(cfg["size"]["bin_edges_um"][:-1],
                         cfg["size"]["bin_edges_um"][1:])),
        dry_per_bin=dry_all, wet_per_bin=wet_all, comp=comp,
        land_weight=land_w, ocean_weight=area_w,
        payload={
            "note": "Sea-salt wet and dry DEPOSITION MASS, the abiotic nutrient "
                    "ledger's carrier for the marine_aerosol candidate. "
                    "Generated by aeolian/scripts/build_sea_salt.py; do not edit.",
            "generated": payload["generated"],
            "git_commit": payload["git_commit"],
            "climatology": payload["climatology"],
            "climatology_stage": payload["climatology_stage"],
            "source_function": "Grythe et al. (2014) equation 7, all three "
                               "modes, transported and deposited by the same "
                               "steady state that produces the burden",
            "variant": "all_modes",
            "no_spume_land_mean_total_mg_m2_earth_year": round(float(
                ((dry_ns.sum(axis=0) + wet_ns.sum(axis=0)) * land_w).sum()
                / land_w.sum() * EARTH_YEAR_S * 1e6), 4),
            "inputs": payload["inputs"],
        })
    anchor = cfg["composition"].get("magnitude_anchor", {})
    ca = dep["land_mean_deposition_mg_m2_earth_year"].get("Ca", {}).get("total")
    print(f"land-mean deposition, mg per m2 per Earth year:")
    for element, v in dep["land_mean_deposition_mg_m2_earth_year"].items():
        print(f"   {element:3s} {v['total']:10.3f}   (dry {v['dry']:.3f} "
              f"wet {v['wet']:.3f})")
    if anchor and ca is not None:
        lo, hi = anchor["calcium_mg_m2_earth_year"]
        print(f"   Ca against {anchor['source']}: {lo}-{hi} mg/m2/yr   "
              f"{'inside' if lo <= ca <= hi else 'OUTSIDE'}")
    print(f"wrote {rel(args.output_deposition)} and "
          f"{rel(args.output_deposition_nc)}")


if __name__ == "__main__":
    main()
