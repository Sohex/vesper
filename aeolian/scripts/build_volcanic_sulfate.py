#!/usr/bin/env python3
"""Volcanic sulfate: source from this world's outgassing, then the same chain.

    python aeolian/scripts/build_volcanic_sulfate.py
    python aeolian/scripts/build_volcanic_sulfate.py --bracket

CLIM-28, the third aerosol. `notes/audits/unpriced-terms.md` finding 2 asked
what the aerosols that are not mineral dust are worth; sea salt turned out to be
dust's size with the opposite sign, which retired the prior that a second
species is negligible. This asks the same question of the third.

## What makes this one different, and why it is the last computable one

Dust and sea salt are both wind acting on a vast surface, so both are large.
Volcanic sulfur is set by how much magma degasses, and this project already has
that number: `weathering_fluxes.py` computes the CO2 outgassing this world needs
at steady state to balance its own silicate weathering, and a magmatic S/C ratio
converts it. So the source is a genuine coupling to the pedology thermostat
rather than a scaled guess, which is the whole reason this can be computed.

Everything downstream is the machinery the other two already use: OPAC optics at
ambient humidity through `sea_salt_optics.band_average`, the shared
`advect_to_steady_state`, and `dust_forcing.shortwave_forcing`, so three
aerosols are priced through one two-stream expression rather than three.

## What this deliberately does not include, and cannot

**Explosive eruptions.** Stratospheric sulfate from a large eruption is where
volcanic aerosol does its climatic work on Earth, and it is episodic. An
episodic source needs a frequency-magnitude distribution, which needs an
eruption history; `docs/src/reference/no-time-axis.md` says why this project has none. What
the arc classes carry is a place, not a rate in time. The passive flux here is
the continuously degassing part and the config's upper bracket carries roughly
what Earth's sporadic contribution adds to it.

**Everything sulphur does that is not volcanic.** On Earth the largest natural
sulphur source is marine biogenic, and this project has no marine biosphere at
all: LPJ-GUESS is terrestrial. That is a bigger gap than this component fills
and it is recorded in the audit rather than papered over here.
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
from builds import grid_export, mesh_export
from gridding import land_fraction_of_class
from orogen import Export
from paths import best_available_climatology, rel
from aerosol_deposition import write_deposition
from sea_salt_optics import band_average, growth, read_index_table, read_size_table

sys.path.insert(0, str(PROJECT_ROOT / "exoplasim" / "scripts"))
from dust_forcing import backscatter_fraction, shortwave_forcing  # noqa: E402
from dust_optics import critical_surface_albedo, stellar_weights  # noqa: E402

sys.path.insert(0, str(PROJECT_ROOT / "lib"))
from stellar import band1_fraction  # noqa: E402

SULFATE_CONFIG = PROJECT_ROOT / "aeolian" / "config" / "volcanic_sulfate.yaml"
WEATHERING = PROJECT_ROOT / "pedology" / "analysis" / "weathering_fluxes.json"
OUT_JSON = ANALYSIS / "volcanic_sulfate.json"
OUT_NC = ANALYSIS / "volcanic_sulfate.nc"
# The nutrient carrier, and a SEPARATE file on purpose: the abiotic nutrient
# ledger's screen names `aeolian/analysis/volcanic_sulfate.nc` in
# `forbidden_carriers` because that file is an optics product, so the mass has
# to live somewhere whose whole content is mass. See aerosol_deposition.py.
OUT_DEP_JSON = ANALYSIS / "volcanic_sulfate_deposition.json"
OUT_DEP_NC = ANALYSIS / "volcanic_sulfate_deposition.nc"

R_DRY = 287.05
EARTH_YEAR_S = 365.25 * 86400.0
BAND_SPLIT_UM = 0.75
BAND_LO_UM, BAND_HI_UM = 0.34, 4.00


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_commit() -> str | None:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT,
                              capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:
        return None


def sulfur_flux_kg_per_year(cfg: dict, outgassing_over_earth: float,
                            so2_tg: float) -> float:
    """Global sulphur mass emitted per Earth year, from the carbon anchor.

    Earth's measured passive SO2 flux, converted to sulphur, scaled by the
    outgassing this world needs relative to Earth's. The scaling assumes an
    Earth-like magmatic S/C ratio and an Earth-like subaerial share, both of
    which the config declares and neither of which this project constrains.
    """
    s = cfg["source"]
    s_per_so2 = s["sulfur_molar_mass"] / s["so2_molar_mass"]
    return so2_tg * 1e9 * s_per_so2 * outgassing_over_earth


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--climatology", type=Path, default=None,
                    help="defaults to the configured bootstrap_climatology, "
                         "which is what the volcanic_sulfate step declares it "
                         "needs")
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--sulfate-config", type=Path, default=SULFATE_CONFIG)
    ap.add_argument("--weathering", type=Path, default=WEATHERING)
    ap.add_argument("--bracket", action="store_true",
                    help="run every declared bracket end and report the spread")
    ap.add_argument("--output", type=Path, default=OUT_JSON)
    ap.add_argument("--output-nc", type=Path, default=OUT_NC)
    ap.add_argument("--output-deposition", type=Path, default=OUT_DEP_JSON)
    ap.add_argument("--output-deposition-nc", type=Path, default=OUT_DEP_NC)
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text())
    cfg = yaml.safe_load(args.sulfate_config.read_text())
    # `advect_to_steady_state` is shared with build_dust.py and reads the
    # planet radius off the config dict it is handed, because the CFL step and
    # the polar cos(lat) floor are both lengths on this world's sphere. It is
    # injected the same way build_dust.py injects it. Without it the transport
    # raises KeyError before producing anything, which is what it did.
    cfg["_planet_radius_earth"] = config["planet"]["radius_earth"]
    gravity = float(config["planet"]["gravity_m_s2"])
    # THE BEST AVAILABLE, on the same terms as its two sibling aerosol steps:
    # the burden is transported by the modelled winds, grown at the modelled
    # humidity and washed out by the modelled precipitation, so it is a
    # function of the climate STATE and the bootstrap is the best answer only
    # on the pass where it is the only one. The `needs: bootstrap_climatology`
    # edge in config/pipeline.yaml is unchanged and states what must EXIST.
    clim_path, clim_stage = best_available_climatology(args.climatology)

    if not args.weathering.is_file():
        raise SystemExit(
            f"{args.weathering} missing. The source term IS the outgassing "
            f"requirement, so this cannot fall back to a literal: run "
            f"pedology/scripts/weathering_fluxes.py first.")
    weathering = json.loads(args.weathering.read_text())
    outgassing = float(
        weathering["carbon_balance"]["implied_outgassing_over_earth"])

    size = read_size_table(PROJECT_ROOT / cfg["growth"]["size_table"])
    index = read_index_table(PROJECT_ROOT / cfg["growth"]["index_table"])
    mode = cfg["growth"]["mode"]

    with Dataset(clim_path) as ds:
        lat = np.asarray(ds["lat"][:], dtype=float)
        lon = np.asarray(ds["lon"][:], dtype=float)
        lev = np.asarray(ds["lev"][:], dtype=float)
        spd_all = np.asarray(ds["spd"][:], dtype=float)
        ua = np.asarray(ds["ua"][:], dtype=float)
        va = np.asarray(ds["va"][:], dtype=float)
        ta = np.asarray(ds["ta"][:], dtype=float)[:, -1]
        hur = np.asarray(ds["hur"][:], dtype=float)[:, -1]
        ps = np.asarray(ds["ps"][:], dtype=float)
        pr = np.asarray(ds["pr"][:], dtype=float)
        rst = np.asarray(ds["rst"][:], dtype=float)
        rsut = np.asarray(ds["rsut"][:], dtype=float)
        alb1 = np.asarray(ds["alb1"][:], dtype=float)
        alb2 = np.asarray(ds["alb2"][:], dtype=float)
        lsm = np.asarray(ds["lsm"][:], dtype=float)[0]

    bad = flag_anomalous_bins(spd_all)
    weights = np.asarray(climatology.bin_weights(np.arange(spd_all.shape[0])),
                         dtype=float)
    if bad:
        weights[bad] = 0.0
    weights = weights / weights.sum()

    def annual(field):
        return np.tensordot(weights, field, axes=(0, 0))

    coslat = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, lon.size))
    coslat_solver = np.maximum(
        coslat, cfg["transport"].get("polar_coslat_floor", 1e-3))
    ta_a, ps_a, pr_a = annual(ta), annual(ps) * 100.0, annual(pr)
    rho_a = ps_a * lev[-1] / (R_DRY * ta_a)
    # `hur` is in percent; see build_sea_salt for what the attribute used to say.
    rh_cell = annual(hur)
    if not 1.0 < float(np.nanmax(rh_cell)) <= 105.0:
        raise SystemExit("relative humidity is not a percentage; check `hur`")
    rh_cell = np.clip(rh_cell, cfg["growth"]["rh_min_percent"],
                      cfg["growth"]["rh_max_percent"])
    rsdt_a = annual(rst) + annual(rsut)
    alb1_a = np.clip(annual(alb1), 0.0, 1.0)
    alb2_a = np.clip(annual(alb2), 0.0, 1.0)

    below = lev >= cfg["transport"]["steering_sigma"]
    u_steer = annual(ua[:, below].mean(axis=1))
    v_steer = annual(va[:, below].mean(axis=1))

    # --- where the vents are ------------------------------------------------
    export = Export(mesh_export(config))
    grid_dir = grid_export(config)
    lit = json.loads((export.root / "manifest.json").read_text(
        encoding="utf-8"))["lithology"]
    ids = {r["code"]: int(r["id"]) for r in lit["rockClasses"]}
    substrate = export.field("substrate_class")
    missing = [c for c in cfg["emission"]["arc_classes"] if c not in ids]
    if missing:
        raise SystemExit(
            f"lithology has no class {missing}; the arc classes named in "
            f"aeolian/config/volcanic_sulfate.yaml are not in this build's "
            f"rock table, and guessing a substitute would invent a geography")
    arc = np.zeros((lat.size, lon.size))
    for name in cfg["emission"]["arc_classes"]:
        arc = arc + land_fraction_of_class(export, grid_dir,
                                           substrate == ids[name])
    if arc.sum() <= 0:
        raise SystemExit(
            f"no arc terrain found for {cfg['emission']['arc_classes']}; the "
            f"source has nowhere to sit and a uniform fallback would be an "
            f"invented geography")
    arc_area_fraction = float((arc * coslat).sum() / coslat.sum())

    planet_area = 4.0 * np.pi * (1.2 * 6.371e6) ** 2
    b1 = band1_fraction()
    lam_s, f_s = stellar_weights()

    # Optics: one lognormal per humidity, per gram of the 75% acid solution the
    # size table is denominated in.
    r_mod_dry = size[mode][0][2]
    rho_dry = size[mode][0][3]
    sigma_g = size[mode][0][4]
    r_min, r_max = size[mode][0][0], size[mode][0][1]
    optics_rh = sorted(size[mode])
    mee1_t, mee2_t, ssa1_t, g1_t, gf_t = [], [], [], [], []
    for rh in optics_rh:
        gf, _, _ = growth(size, mode, rh)
        r_grid = np.logspace(np.log10(r_min * gf), np.log10(r_max * gf), 400)
        dndlnr = np.exp(-0.5 * (np.log(r_grid / (r_mod_dry * gf))
                                / np.log(sigma_g)) ** 2)
        rho_eff = rho_dry / gf ** 3      # wet volume -> reference dry mass
        m1, s1, g1 = band_average(BAND_LO_UM, BAND_SPLIT_UM, index[mode][rh],
                                  r_grid, dndlnr, rho_eff, lam_s, f_s)
        m2, s2, g2 = band_average(BAND_SPLIT_UM, BAND_HI_UM, index[mode][rh],
                                  r_grid, dndlnr, rho_eff, lam_s, f_s)
        mee1_t.append(m1)
        mee2_t.append(m2)
        ssa1_t.append(s1)
        g1_t.append(g1)
        gf_t.append(gf)
    optics_rh = np.asarray(optics_rh, dtype=float)

    def solve(cfg):
        s = cfg["source"]
        s_kg = sulfur_flux_kg_per_year(cfg, outgassing, s["earth_so2_tg_per_year"])
        per_s = (s["h2so4_molar_mass"] / s["sulfur_molar_mass"]
                 / s["acid_mass_fraction"])
        aerosol_kg_per_year = s_kg * s["aerosol_yield"] * per_s
        # Distributed over the arc, per unit area, as a steady flux.
        cell_area = coslat / coslat.sum() * planet_area
        share = arc * cell_area
        emission = (aerosol_kg_per_year / EARTH_YEAR_S
                    * share / share.sum() / cell_area)

        gf = np.interp(rh_cell, optics_rh, gf_t)
        rho_p = 1000.0 + (rho_dry * 1000.0 - 1000.0) / gf ** 3
        v_s = settling_velocity(2.0 * r_mod_dry * gf, rho_p, rho_a, ta_a,
                                gravity)
        scav = cfg["removal"]
        p_mm_h = np.maximum(pr_a, 0.0) * 1000.0 * 3600.0
        lam_ref = 1.0 / (scav["wet_lifetime_hours"] * 3600.0)
        lam_wet = lam_ref * (p_mm_h
                             / scav["reference_precipitation_mm_per_hour"]
                             ) ** scav["scavenging_b"]
        loss = v_s / cfg["transport"]["scale_height_m"] + lam_wet
        m, steps, converged = advect_to_steady_state(
            emission, u_steer, v_steer, loss, lat, lon, cfg)
        residual = abs(float((emission * coslat_solver).sum())
                       - float((loss * m * coslat_solver).sum())
                       ) / max(float((emission * coslat_solver).sum()), 1e-30)
        aod1 = m * 1000.0 * np.interp(rh_cell, optics_rh, mee1_t)
        aod2 = m * 1000.0 * np.interp(rh_cell, optics_rh, mee2_t)
        ssa = float(np.interp(np.mean(rh_cell), optics_rh, ssa1_t))
        beta = backscatter_fraction(float(np.interp(np.mean(rh_cell),
                                                    optics_rh, g1_t)))
        forcing = (shortwave_forcing(aod1, ssa, beta, alb1_a, b1 * rsdt_a)
                   + shortwave_forcing(aod2, ssa, beta, alb2_a,
                                       (1.0 - b1) * rsdt_a))
        return {
            "sulfur_tg_per_earth_year": round(s_kg / 1e9, 4),
            "aerosol_tg_per_earth_year": round(aerosol_kg_per_year / 1e9, 4),
            "burden_mg_m2_global_mean": round(
                float((m * coslat).sum() / coslat.sum()) * 1e6, 5),
            "lifetime_days_area_mean": round(
                float((1.0 / loss * coslat).sum() / coslat.sum()) / 86400.0, 3),
            "optical_depth_global_mean": round(
                float((b1 * aod1 + (1 - b1) * aod2).sum() * 0 +
                      ((b1 * aod1 + (1 - b1) * aod2) * coslat).sum()
                      / coslat.sum()), 6),
            "toa_shortwave_forcing_w_m2_global": round(
                float((forcing * coslat).sum() / coslat.sum()), 5),
            "single_scattering_albedo": round(ssa, 6),
            "critical_surface_albedo": round(
                critical_surface_albedo(ssa, beta), 4),
            "transport_steps": steps, "converged": bool(converged),
            "mass_residual": round(residual, 6),
        }, m, (b1 * aod1 + (1 - b1) * aod2), emission, (
            v_s / cfg["transport"]["scale_height_m"] * m, lam_wet * m)

    central, burden, aod, emission, (dry_dep, wet_dep) = solve(cfg)
    if not central["converged"] or central["mass_residual"] >= 0.10:
        raise SystemExit(f"transport failed: {central}")

    bracket = {}
    if args.bracket:
        import copy
        for block, key, bkey in (
                ("source", "earth_so2_tg_per_year",
                 "earth_so2_tg_per_year_bracket"),
                ("source", "aerosol_yield", "aerosol_yield_bracket"),
                ("removal", "wet_lifetime_hours", "wet_lifetime_hours_bracket"),
                ("transport", "scale_height_m", "scale_height_bracket_m")):
            for value in cfg[block][bkey]:
                if value == cfg[block][key]:
                    continue
                alt = copy.deepcopy(cfg)
                alt[block][key] = value
                r = solve(alt)[0]
                bracket[f"{block}.{key}={value}"] = {
                    "optical_depth_global_mean": r["optical_depth_global_mean"],
                    "toa_shortwave_forcing_w_m2_global":
                        r["toa_shortwave_forcing_w_m2_global"]}

    print(f"climatology            {rel(clim_path)}")
    print(f"outgassing requirement {outgassing:.4f} x Earth, from "
          f"{rel(args.weathering)}")
    print(f"arc terrain            {arc_area_fraction*100:.3f}% of the planet's "
          f"surface carries an arc class")
    print(f"sulphur                {central['sulfur_tg_per_earth_year']:.2f} "
          f"Tg S per Earth year, {central['aerosol_tg_per_earth_year']:.2f} Tg "
          f"of 75% acid aerosol")
    print(f"burden                 "
          f"{central['burden_mg_m2_global_mean']:.4f} mg/m2 global mean, "
          f"lifetime {central['lifetime_days_area_mean']:.2f} days")
    print(f"optical depth          "
          f"{central['optical_depth_global_mean']:.5f} global mean")
    print(f"TOA shortwave forcing  "
          f"{central['toa_shortwave_forcing_w_m2_global']:+.4f} W/m2 global "
          f"mean, single-scattering albedo "
          f"{central['single_scattering_albedo']:.5f}")
    for name, v in bracket.items():
        print(f"  bracket {name:42s} AOD "
              f"{v['optical_depth_global_mean']:.5f}  W/m2 "
              f"{v['toa_shortwave_forcing_w_m2_global']:+.4f}")

    payload = {
        "note": "Offline volcanic-sulfate source, burden, optical depth and "
                "shortwave forcing. PASSIVE degassing only: the eruptive part "
                "is episodic and needs an eruption history this project does "
                "not have, per docs/src/reference/no-time-axis.md. Generated by "
                "aeolian/scripts/build_volcanic_sulfate.py; do not edit.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(),
        "climatology": str(rel(clim_path)),
        # WHICH STAGE this burden was transported by. lib/paths.py.
        "climatology_stage": clim_stage,
        "excluded_time_bins": bad,
        "outgassing_over_earth": outgassing,
        "outgassing_source": str(rel(args.weathering)),
        "arc_area_fraction_of_planet": round(arc_area_fraction, 6),
        "central": central,
        "bracket": bracket,
        "inputs": {p.name: sha256(p) for p in
                   (args.sulfate_config, args.weathering, args.config)},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    with Dataset(args.output_nc, "w") as out:
        out.createDimension("lat", lat.size)
        out.createDimension("lon", lon.size)
        for name, data, units in (("lat", lat, "deg"), ("lon", lon, "deg")):
            v = out.createVariable(name, "f8", (name,))
            v.units = units
            v[:] = data
        for name, data, units, long_name in (
                ("burden", burden, "kg m-2", "sulfate column burden"),
                ("aod", aod, "1", "sulfate optical depth, band combined"),
                ("emission", emission, "kg m-2 s-1", "sulfate source"),
                ("arc_fraction", arc, "1", "arc-class land fraction")):
            v = out.createVariable(name, "f8", ("lat", "lon"), zlib=True)
            v.units = units
            v.long_name = long_name
            v[:] = np.asarray(data)
        out.note = payload["note"]
        out.generated = payload["generated"]
    print(f"\nwrote {rel(args.output)} and {rel(args.output_nc)}")

    # -- the nutrient carrier, mass and only mass ----------------------------
    #
    # ANUT-7 registers volcanic sulfate deposition as a sulfur source to
    # soil_solution and it had no carrier. One size, because the source is a
    # single lognormal mode, so the bin axis has length one rather than being
    # absent: the ledger reads both carriers through one convention.
    land_w = coslat * (lsm >= 0.5)
    ocean_w = coslat * (lsm < 0.5)
    dep = write_deposition(
        args.output_deposition_nc, args.output_deposition,
        lat=lat, lon=lon,
        bins_um=[[2.0 * r_mod_dry, 2.0 * r_mod_dry]],
        dry_per_bin=dry_dep[None, :, :], wet_per_bin=wet_dep[None, :, :],
        comp=cfg["composition"], land_weight=land_w, ocean_weight=ocean_w,
        payload={
            "note": "Volcanic sulfate wet and dry DEPOSITION MASS, the abiotic "
                    "nutrient ledger's carrier for the "
                    "volcanic_sulfate_deposition candidate. PASSIVE degassing "
                    "only. Generated by "
                    "aeolian/scripts/build_volcanic_sulfate.py; do not edit.",
            "generated": payload["generated"],
            "git_commit": payload["git_commit"],
            "climatology": payload["climatology"],
            "climatology_stage": payload["climatology_stage"],
            "outgassing_over_earth": outgassing,
            "outgassing_source": payload["outgassing_source"],
            "inputs": payload["inputs"],
        })
    print("land-mean deposition, mg per m2 per Earth year:")
    for element, v in dep["land_mean_deposition_mg_m2_earth_year"].items():
        print(f"   {element:3s} {v['total']:10.4f}   (dry {v['dry']:.4f} "
              f"wet {v['wet']:.4f})")
    print(f"wrote {rel(args.output_deposition)} and "
          f"{rel(args.output_deposition_nc)}")


if __name__ == "__main__":
    main()
