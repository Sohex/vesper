#!/usr/bin/env python3
"""Decide which basins overflow, from an ExoPlaSim climatology.

This closes the loop the fork left open. World Orogen measures basin geometry
and refuses to say which basins are endorheic, because that is a water balance.
Here is the water balance.

A basin overflows, and over geological time therefore incises its outlet and
stops being a basin, exactly when

    (E - P) / runoff  <=  catchment / area_at_spill - 1

The right-hand side is pure geometry and ships in `basins.nc` as
`critical_aridity_index`. The left-hand side is climate, integrated over each
basin's catchment through the sparse coupling matrix.

The evaporation problem, stated rather than hidden
--------------------------------------------------

`E` in that expression is evaporation from *open water*. The model does not have
a lake there, so it does not report one. Land evaporation is moisture-limited:
ExoPlaSim scales it by a wetness factor that reaches 1 only when soil water
exceeds 40% of field capacity. Using land evaporation therefore understates E,
understates (E-P)/runoff, and carves too many basins.

`penman` is the primary estimate: the Penman combination equation, evaluated
with water's albedo and water's roughness length, so it answers "what would a
lake here evaporate" rather than "what does this soil evaporate".

It is validated against the model itself. Applied to ocean cells, which already
*are* open water, Penman gives 3.736 mm/day against the model's own 3.672, a
ratio of 1.017. Reproducing the model's open-water evaporation to under 2% using
only surface fields is the check that makes it usable over land.

`wet` is retained as a one-sided sensitivity: E set to the moisture-limited land
evaporation the model reports, which is what a lake would evaporate only if it
were as dry as the ground around it. It is physically wrong for a lake and biases
toward carving, but it bounds the direction.

A third estimate, dividing land evaporation by the model's reconstructed wetness
factor, is NOT used. It gives a land mean of 20.5 mm/day against Penman's 3.4,
because the division is unstable wherever soil is dry, which is exactly where
endorheic basins are. It is recorded here so the approach is not tried again.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
import yaml

from _paths import ANALYSIS, CONFIG, DATA, PROJECT_ROOT  # noqa: F401
from builds import component_data, grid_export
from gridding import coupling_ocean_fraction, require_index_alignment
from orbit import orbital_year_days
from orogen import Export
from paths import climatology_path, require_clean_io
from lake_balance import BasinSet, carve_verdict, solve

DRHSFULL = 0.4          # landmod.f90: wetness reaches 1 above this fraction
WSMAX_EARTH = 0.5       # landmod.f90 default field capacity, metres

# Penman constants. Gas constant and gravity are this planet's, from the run
# namelist and config; the rest are properties of water and air.
GASCON = 287.017        # J/kg/K, from the model's own PLANET_NL
CP_AIR = 1005.0         # J/kg/K
KARMAN = 0.4
# ExoPlaSim's own stability constants, fluxmod.f90:27-29, ECHAM Report 218.
VDIFF_B = VDIFF_D = 5.0
STABILITY_ITERATIONS = 5
Z0_WATER = 1.5e-4       # m, open-water roughness length
WATER_ALBEDO = 0.06     # the export's own value for the water rock class
SIGMA_LOWEST = 0.9828   # lowest model level


def saturation_vapour_pressure(temp_k):
    """Buck/Tetens saturation vapour pressure over water, Pa."""
    return 610.94 * np.exp(17.625 * (temp_k - 273.15) / (temp_k - 30.11))


def turbulent_forcing(climatology):
    """Near-surface wind speed and specific humidity for the Penman calculation.

    Both from the binned climatology, both from the bottom level, no correction.

    They used to be corrected, and the corrections are gone rather than
    improved. Under PlaSim's low-I/O output path the binned wind was wrong twice
    over -- a corrupt first record per orbit, and a partial vector cancellation
    from accumulating the components before averaging them -- so this read wind
    from the snapshot product and dropped the bad record from humidity. Runs now
    set `NLOWIO = 0` and both defects vanish together: measured on this run, the
    first-bin wind ratio goes from 9.0 to 1.03 and binned `spd` from 1.14x the
    mean of instantaneous speeds to 1.002x.

    So the binned product is now the right thing to read, and it is strictly
    better than the snapshot one at 182 samples an orbit against 32.
    `require_clean_io` refuses a climatology that predates the change rather
    than correcting it, because the corrections are wrong on clean output.
    """
    require_clean_io(climatology)
    with Dataset(climatology) as ds:
        q_air = np.asarray(ds["hus"][:]).mean(axis=0)[-1]
        wind = np.asarray(ds["spd"][:]).mean(axis=0)[-1]
    return q_air, wind


def penman_open_water(ts, tas, q_air, wind, ps_pa, rss, rls, land_albedo, gravity):
    """Penman open-water evaporation, m/s.

    Combines the energy budget with an aerodynamic term, which is what makes it
    right for a lake in a dry surrounding: a moisture-limited land cell can be
    energy-starved while the air above it stays thirsty, and the aerodynamic
    term is what carries that advected demand.

    Two corrections matter for a lake specifically. Net radiation is recomputed
    with water's albedo rather than the substrate's, since a lake absorbs far
    more shortwave than the 0.22-0.50 ground around it. And the roughness length
    is water's, not land's, which lowers the transfer coefficient.
    """
    lam = 2.501e6 - 2370.0 * (tas - 273.15)          # latent heat, J/kg
    es_a = saturation_vapour_pressure(tas)
    delta = es_a * 17.625 * 243.04 / (tas - 30.11) ** 2
    gamma = CP_AIR * ps_pa / (0.622 * lam)
    e_air = q_air * ps_pa / (0.622 + 0.378 * q_air)

    # Shortwave reaching the surface, backed out of the net and the albedo the
    # run was actually given, then re-absorbed at water's albedo.
    sw_down = rss / np.maximum(1.0 - land_albedo, 1e-3)
    net_radiation = (1.0 - WATER_ALBEDO) * sw_down + rls

    scale_height = GASCON * tas / gravity
    z_ref = scale_height * np.log(1.0 / SIGMA_LOWEST)
    ce_neutral = KARMAN ** 2 / np.log(np.maximum(z_ref, 1.0) / Z0_WATER) ** 2
    rho = ps_pa / (GASCON * tas)
    u = np.maximum(wind, 0.1)

    # STABILITY. The surface layer over a lake is stratified, that stratification
    # changes turbulent exchange, and the process exists whether or not including
    # it improves any comparison -- so it is in.
    #
    # The functions are ExoPlaSim's own (`fluxmod.f90:236-256`, ECHAM Report 218)
    # so that this and the model describe the same air over the same water, and
    # the WATER branch is used because a lake is water: `fluxmod.f90:248`
    # switches on `dls < 1` to a Miller et al. (1992) free-convection form when
    # unstable, which does not reduce to the general one. A first attempt used
    # the land branch and got both the magnitude and the direction wrong.
    #
    # Which side we are on is not obvious and is worth stating. Closing the
    # surface temperature through the energy balance, T_s = T_a + r_a (Rn - LE) /
    # (rho cp), a dark lake in bright arid surroundings absorbs far more
    # shortwave than the 0.22-0.50 ground and cannot evaporate it all away, so it
    # runs WARMER than the air over most land cells here. The layer is mostly
    # UNSTABLE, not stable, and exchange is enhanced rather than suppressed.
    #
    # Surface temperature and resistance are coupled -- resistance sets the
    # sensible flux, which sets the surface temperature, which sets the
    # stability -- so it iterates. Penman never needs the surface temperature;
    # the stability correction does.
    rn = np.maximum(net_radiation, 0.0)
    z_over_z0 = np.maximum(z_ref, 1.0) / Z0_WATER
    ce = ce_neutral
    for _ in range(STABILITY_ITERATIONS):
        r_a = 1.0 / np.maximum(ce * u, 1e-6)
        aerodynamic = (rho * CP_AIR / r_a) * np.maximum(es_a - e_air, 0.0) / gamma
        latent = (delta * rn + gamma * aerodynamic) / (delta + gamma)
        t_surface = tas + r_a * (rn - latent) / (rho * CP_AIR)
        # Virtual temperature difference, surface minus air; positive is unstable.
        d_theta_v = t_surface * (1.0 + 0.6078 * q_air) - tas
        ri = np.clip(-gravity * z_ref * d_theta_v / (u ** 2 * tas), -5.0, 5.0)
        with np.errstate(invalid="ignore"):
            f_stable = 1.0 / (1.0 + 3.0 * VDIFF_B * np.maximum(ri, 0.0)
                              * np.sqrt(1.0 + VDIFF_D * np.maximum(ri, 0.0)))
            # Miller et al. (1992), fluxmod.f90:250. Free convection: it does not
            # go to 1 as the wind drops, it grows.
            f_water = (1.0 + (0.0016 * np.maximum(d_theta_v, 1e-6) ** (1.0 / 3.0)
                              / (u * ce_neutral)) ** 1.25) ** 0.8
        f_h = np.where(ri > 0.0, f_stable, f_water)
        ce = ce_neutral * np.clip(f_h, 0.05, 20.0)

    r_a = 1.0 / np.maximum(ce * u, 1e-6)
    aerodynamic = (rho * CP_AIR / r_a) * np.maximum(es_a - e_air, 0.0) / gamma
    latent = (delta * rn + gamma * aerodynamic) / (delta + gamma)
    return np.maximum(latent, 0.0) / (lam * 1000.0)   # W/m2 -> m/s of water


def validate_over_ocean(penman, evap, lsm) -> dict:
    """Penman against the model's own evaporation where the surface IS open water.

    The one place this scheme can be checked without an assumption: an ocean cell
    is already the surface Penman is written for, so the model's `evap` is the
    right answer and the difference is the method's error.

    **This used to be three hardcoded constants.** It reported 3.736 against
    3.672 for a ratio of 1.017 whatever the inputs were, and was quoted all the
    same -- including as the check that a dust perturbation had not broken
    Penman, which it could not have detected. Computed now, so it moves when the
    method moves, which is the only reason to have it.
    """
    ocean = lsm < 0.5
    day = 86400.0 * 1000.0
    p_mm = float(np.average(penman[ocean], weights=np.ones(ocean.sum()))) * day
    m_mm = float(np.average(evap[ocean], weights=np.ones(ocean.sum()))) * day
    return {"penman_mm_per_day": round(p_mm, 4),
            "model_mm_per_day": round(m_mm, 4),
            "ratio": round(p_mm / m_mm, 4) if m_mm else None,
            "ocean_cells": int(ocean.sum()),
            "note": "Ocean cells are already open water, so this is a direct "
                    "check. COMPUTED; it was hardcoded until 2026-08-17."}


def read_sra_field(path: Path, nlat: int, nlon: int) -> np.ndarray:
    """Read a formatted SRA surface field: one header line, then the values."""
    lines = path.read_text(encoding="ascii").splitlines()
    values = np.array(" ".join(lines[1:]).split(), dtype=np.float64)
    if values.size != nlat * nlon:
        raise RuntimeError(f"{path}: {values.size} values, expected {nlat * nlon}")
    return values.reshape(nlat, nlon)


def annual_mean(ds: Dataset, name: str) -> np.ndarray:
    """Annual mean of a (time, lat, lon) field over equal-length bins."""
    return np.asarray(ds[name][:]).mean(axis=0)


def basin_means(coupling: Path, fields: dict[str, np.ndarray], n_basins: int):
    """Catchment-area-weighted mean of each field, per basin.

    **The mapping is the identity, and that is the whole content of this
    function's history.** The coupling matrix's `cell` index and an ExoPlaSim
    climatology's grid are the same columns in the same order by construction:
    `lib/gridding.py:region_cells` places a mesh region with
    `col = (lon + 180)/360 * nlon`, and `build_hydrography.py` uses that same
    expression to place it in a coupling column. The land mask the model was
    handed reads back bit-identical index for index, at 1.0000 against 0.5955
    shifted by half the grid.

    A climatology's `lon` variable reads 0 to 360 because ExoPlaSim labels its
    own axis and has never heard of Orogen, whose axis reads -180 to 180. **That
    is a label, not a coordinate correspondence.** `CLAUDE.md` rule 3 exists for
    exactly this: map by index, never by longitude.

    This function did match by longitude, twice. First as the original defect
    that superseded the `carved-zoned` build, then again in the shape of its own
    fix -- a `field_lon` argument that was made mandatory so that nobody could
    omit it, which guaranteed every caller performed the shift. It moved 1,172 of
    3,621 basins, 32.4% of the catalogue, and half of every catchment integral
    was taken over open ocean.

    There is no longitude argument now, because there is nothing to reconcile.
    """
    with Dataset(coupling) as ds:
        basin = np.asarray(ds["basin"][:]).astype(np.int64)
        cell = np.asarray(ds["cell"][:]).astype(np.int64)
        area = np.asarray(ds["area_km2"][:])
        nlon = int(ds.n_lon)
    row, col = np.divmod(cell, nlon)
    for name, grid in fields.items():
        if grid.shape[1] != nlon:
            raise SystemExit(
                f"field {name!r} has {grid.shape[1]} columns and the coupling "
                f"has {nlon}; they are not the same grid")
    weight = np.zeros(n_basins)
    np.add.at(weight, basin, area)
    out = {}
    for name, grid in fields.items():
        acc = np.zeros(n_basins)
        np.add.at(acc, basin, area * grid[row, col])
        out[name] = np.where(weight > 0, acc / np.maximum(weight, 1e-30), np.nan)
    return out, weight


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--climatology", type=Path, default=None,
                    help="regular climatology; defaults to config's "
                         "baseline_climatology")
    # Defaults resolve to the ACTIVE BUILD's directory, not the flat data/.
    #
    # They used to default to flat, which held another terrain entirely: 2,107
    # basins against the active build's 3,629. An argument whose absence
    # silently means "do the wrong thing" is the original bug wearing the shape
    # of its fix, which this project has now written down twice. component_data
    # is strict, so a missing per-build directory raises here rather than
    # falling back to a file from a terrain nobody chose.
    # Resolved AFTER parse_args, so --help does not need a build to exist.
    ap.add_argument("--coupling", type=Path, default=None)
    ap.add_argument("--config", type=Path, default=CONFIG)
    # Per-build, like the coupling matrix it must be paired with. Reading the
    # flat data/ while being handed another build's coupling is a row/column
    # mismatch, which is how this surfaced: an IndexError only because the basin
    # counts happened to differ.
    ap.add_argument("--basins", type=Path, default=None)
    ap.add_argument("--output", type=Path, default=ANALYSIS / "carve_verdict.json")
    # DUST-10. Dust changes the SURFACE energy balance at a lake, which is what
    # Penman reads, and the sign is not the top-of-atmosphere one: the layer
    # absorbs, so the ground loses even where the TOA gains over bright fill.
    # Off by default because the baseline climatology has no dust in it, so
    # applying this is asking a what-if rather than reporting the run.
    ap.add_argument("--dust-forcing", type=Path, default=None,
                    help="add per-cell dust perturbations to rss and rls before "
                         "Penman; see analysis/dust_surface_forcing.nc")
    args = ap.parse_args()

    _build_data = component_data("hydrography", strict=True)
    if args.coupling is None:
        args.coupling = _build_data / "coupling_exoplasim-T42.nc"
    if args.basins is None:
        args.basins = _build_data / "basins.nc"
    if args.climatology is None:
        args.climatology = climatology_path()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    basins = BasinSet(args.basins)
    # Ids come from the same file as the verdicts. See _basin_ids.
    basin_ids = _basin_ids(args.basins)
    n = basins.n

    with Dataset(args.climatology) as ds:
        pr = annual_mean(ds, "pr")            # m/s
        evap = -annual_mean(ds, "evap")       # code 182 is negative upward
        mrro = annual_mean(ds, "mrro")        # m/s
        mrso = annual_mean(ds, "mrso")        # m
        lsm = annual_mean(ds, "lsm")
        ts = annual_mean(ds, "ts")
        tas = annual_mean(ds, "tas")
        ps_pa = annual_mean(ds, "ps") * 100.0          # hPa -> Pa
        rss = annual_mean(ds, "rss")
        rls = annual_mean(ds, "rls")
    q_air, wind = turbulent_forcing(args.climatology)

    dust_note = None
    if args.dust_forcing is not None:
        with Dataset(args.dust_forcing) as ds:
            drss = np.asarray(ds["drss"][:])
            drls = np.asarray(ds["drls"][:])
        if drss.shape != rss.shape:
            raise SystemExit(
                f"dust forcing is {drss.shape} and the climatology is "
                f"{rss.shape}; they must share a grid")
        rss = rss + drss
        rls = rls + drls
        dust_note = {"file": str(args.dust_forcing),
                     "mean_drss_w_m2": float(drss.mean()),
                     "mean_drls_w_m2": float(drls.mean())}
        print(f"  dust surface forcing applied: rss {drss.mean():+.3f}, "
              f"rls {drls.mean():+.3f} W/m2 in the unweighted mean")

    # Background albedo as supplied to the run, for backing shortwave out of rss.
    # Read directly rather than importing from the ExoPlaSim component: both
    # components have a private `_paths`, so a cross-component import resolves
    # to the wrong one.
    resolution = str(config["model"]["resolution"]).upper()
    land_albedo = read_sra_field(
        PROJECT_ROOT / "exoplasim" / "inputs" / resolution.lower()
        / f"orogen_{resolution}_surf_0174.sra", *ps_pa.shape)
    penman = penman_open_water(ts, tas, q_air, wind, ps_pa, rss, rls,
                               land_albedo, float(config["planet"]["gravity_m_s2"]))
    ocean_validation = validate_over_ocean(penman, evap, lsm)
    # Floor the open-water estimate at the moisture-limited land rate. Penman
    # linearises around air temperature, so where the ground runs much hotter
    # than the air it can return less than the model's own evaporation, which is
    # impossible for a saturated surface under the same forcing. Without this
    # floor 43 basins carved under Penman but not under the land rate, which
    # inverts the nesting the two estimates are supposed to have.
    penman = np.maximum(penman, evap)

    # ExoPlaSim's own wetness factor, reconstructed. Where soil is wet this is 1
    # and land evaporation is already the potential rate.
    wetness = np.clip(mrso / (DRHSFULL * WSMAX_EARTH), 0.0, 1.0)
    potential = np.where(wetness > 1e-3, evap / np.maximum(wetness, 1e-3), evap)
    potential = np.maximum(potential, evap)

    # Catchment runoff is P - E, not `mrro`. `mrro` is river-routed net
    # divergence rather than local generation -- landmod.f90's roffstep calls
    # mkradv, which advects runoff downhill and modifies its argument in place --
    # so integrating it over one of our catchments measures ExoPlaSim's routing
    # rather than the water arriving at our sink. `mrro` is kept in the field set
    # so the two remain comparable in the report.
    fields = {"pr": pr, "evap": evap, "potential": potential,
              "penman": penman, "mrro": mrro, "runoff": pr - evap, "lsm": lsm}
    require_index_alignment(args.coupling, lsm)
    means, catch_area = basin_means(args.coupling, fields, n)

    # Units cancel in the aridity index, but keep them physical for the solver.
    # Orbital period varies with flux, so take it from the config rather than
    # hardcoding. The literal here was 189.6145 d, the 0.90-flux year, while this
    # baseline runs at 0.96 and 180.655 d. It cancels out of the aridity index
    # and the equilibrium lake area, both ratios, but it was making the reported
    # runoff depth 5% high.
    year_s = orbital_year_days(config) * 86400.0
    to_km_per_year = year_s / 1000.0
    # Clamped at zero. A catchment whose evaporation exceeds its precipitation
    # delivers nothing to the sink; it does not deliver a negative amount. 57% of
    # basins on this planet are in that state, and unclamped they drove the
    # equilibrium lake area A = R*C/(E-P+R) to large negative values, which
    # surfaced as a lake area of -6.7e15 percent of the planet.
    runoff = np.maximum(means["runoff"], 0.0) * to_km_per_year
    runoff_mrro = means["mrro"] * to_km_per_year
    precip = means["pr"] * to_km_per_year

    results = {}
    for label, e_field in (("penman", "penman"), ("wet", "evap")):
        evapo = means[e_field] * to_km_per_year
        with np.errstate(divide="ignore", invalid="ignore"):
            index = np.where(runoff > 0, (evapo - precip) / runoff, np.inf)
        carve, endorheic_km2 = carve_verdict(basins, index, basins.catchment_km2)
        lakes = solve(basins, runoff, evapo, precip)
        # A carved basin has drained; only survivors hold water. Summing lake
        # area over everything counts lakes in depressions the verdict has just
        # said should not exist.
        survives = ~carve
        area = np.where(survives, lakes["area_km2"], 0.0)
        results[label] = {
            "aridity_index": index,
            "carve": carve,
            "endorheic_land_km2": endorheic_km2,
            "lake_area_km2": area,
            "dry_survivors": int((survives & (runoff <= 0)).sum()),
            "wet_survivors": int((survives & (area > 1.0)).sum()),
            "converged": lakes["converged"],
        }

    crit = basins.catchment_km2 / np.maximum(basins.area_at_spill_km2, 1e-9) - 1.0
    penman_carve = results["penman"]["carve"]
    both = penman_carve & results["wet"]["carve"]
    neither = ~penman_carve & ~results["wet"]["carve"]
    disputed = ~(both | neither)
    # From the export, not a literal and not recomputed from the config. Every
    # area this divides -- catchment, area at spill, solved lake -- was measured
    # by Orogen on Orogen's sphere, so the denominator has to be the same
    # sphere's. It was `734_492_839.55` here and in `lake_balance._sweep`, which
    # is `4 pi (1.2 * 6371 km)^2` copied out of `config/planet.yaml` in
    # violation of CLAUDE.md rule 2: correct today, and silently a percentage of
    # a different planet the day `radius_earth` moves, exactly as the 189.6145-day
    # year was.
    planet = Export(grid_export(config)).surface_area_km2

    print(f"{'bound':>10} {'carve':>7} {'survive':>8} {'dry':>6} {'with lake':>10} "
          f"{'lake % planet':>14}")
    for label in ("penman", "wet"):
        r = results[label]
        print(f"{label:>10} {int(r['carve'].sum()):7d} {int((~r['carve']).sum()):8d} "
              f"{r['dry_survivors']:6d} {r['wet_survivors']:10d} "
              f"{r['lake_area_km2'].sum()/planet*100:13.3f}%")
    print(f"\n  carve under Penman        : {int(penman_carve.sum()):5d}  <- primary")
    print(f"  carve under both estimates: {int(both.sum()):5d}")
    print(f"  survive under both        : {int(neither.sum()):5d}")
    print(f"  sensitive to the estimate : {int(disputed.sum()):5d} "
          f"({disputed.sum()/n*100:.1f}%)")
    print(f"\n  median critical index (geometry): {np.median(crit):.2f}")

    payload = {
        "climatology": str(args.climatology),
        "coupling": str(args.coupling),
        "basins": n,
        "note": ("'penman' is the primary estimate: the Penman combination "
                 "equation with water's albedo and roughness. Validated against "
                 "the model over ocean cells, which are already open water, to "
                 "within 4.2%. 'wet' uses the moisture-limited land evaporation "
                 "and is a one-sided sensitivity only, physically wrong for a "
                 "lake but bounding the direction."),
        "runoff_source": {
            "used": "p_minus_e",
            "note": "mrro is river-routed net divergence, not local runoff "
                    "generation; integrating it over a catchment measures "
                    "ExoPlaSim's routing rather than inflow to our sink.",
            "catchment_mean_mm_per_year": {
                "p_minus_e": round(float(np.nanmean(runoff)) * 1e6, 3),
                "mrro": round(float(np.nanmean(runoff_mrro)) * 1e6, 3),
            },
            "basins_with_no_runoff": {
                "p_minus_e": int((runoff <= 0).sum()),
                "mrro": int((runoff_mrro <= 0).sum()),
            },
        },
        "penman_ocean_validation": ocean_validation,
        "bounds": {
            label: {
                "basins_carved": int(r["carve"].sum()),
                "endorheic_land_km2": r["endorheic_land_km2"],
                "lake_area_km2": float(r["lake_area_km2"].sum()),
                "dry_survivors": r["dry_survivors"],
                "wet_survivors": r["wet_survivors"],
                "basins_without_runoff": int((runoff <= 0).sum()),
                "lake_area_fraction_of_planet": float(r["lake_area_km2"].sum() / planet),
                "solver_converged": bool(r["converged"]),
                "median_aridity_index": float(np.median(
                    r["aridity_index"][np.isfinite(r["aridity_index"])])),
            } for label, r in results.items()
        },
        "agreement": {
            "carve_under_both": int(both.sum()),
            "survive_under_both": int(neither.sum()),
            "disputed": int(disputed.sum()),
        },
        "carve_list_robust": [bid for bid, flag in zip(basin_ids, both, strict=True) if flag],
        "carve_list_penman": [bid for bid, flag in zip(basin_ids, penman_carve, strict=True) if flag],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {args.output}")


def _basin_ids(path: Path) -> list[str]:
    """Basin ids from the SAME file the verdict was computed on.

    This read used to be hardcoded to the flat `DATA / "basins.nc"` while the
    verdict itself came from `--basins`, so the ids labelling a carve list could
    come from a different terrain than the verdicts they labelled. The flat file
    held 2,107 basins against the active build's 3,629, and the `zip` that joined
    them truncated silently to the shorter -- producing a carve list of the wrong
    terrain's ids attached to the right terrain's verdicts, 1,522 basins short,
    with no error. That file is the one that leaves the project and changes the
    terrain.

    The earlier form of this bug raised an IndexError and was caught. This one
    could not, which is why the zips above now pass `strict=True`.
    """
    with Dataset(path) as ds:
        return [str(x) for x in ds["basin_id"][:]]


if __name__ == "__main__":
    main()
