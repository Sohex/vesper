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

It is validated against the model itself, on the one surface where the answer is
already known. Applied to ocean cells, which already *are* open water and which
the model gives water's own roughness, Penman reproduces the model's own
evaporation to within a few percent -- the run's number is computed by
`validate_over_ocean` and written into the report, never quoted from here.
Reproducing open-water evaporation over open water, from surface fields alone,
is what makes the estimate usable over a lake the model does not have.

`wet` is a SENSITIVITY, not a bound, and which side it falls on depends on the
regime. E is set to the moisture-limited land evaporation the model reports.
Where the ground is dry that is far below what a lake would evaporate, so it
biases toward carving. Where the ground is wet it is ABOVE it, because this
world's land carries a roughness field with a median z0 of 0.521 m against open
water's 1.5e-4 and therefore a transfer coefficient 6.4 times larger: a smooth
lake in a rough wet landscape evaporates less than the land around it. The two
estimates are consequently not nested, and a floor that asserted they were used
to stand in `main` -- it bound on 60.2% of land cells and decided 73% of the
overflowing basins by clamp rather than by climate.

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
from climatology import annual_mean as weighted_annual_mean, bin_weights
from gridding import (coupling_cells, coupling_ocean_fraction, coupling_path,
                      require_index_alignment)
from orbit import orbital_year_days
from orogen import Export
from paths import climatology_path, rel, require_clean_io
from provenance import staged_surface_field
from lake_balance import BasinSet, carve_verdict, solve
from lapse import reference_height_m

DRHSFULL = 0.4          # landmod.f90: wetness reaches 1 above this fraction
WSMAX_EARTH = 0.5       # landmod.f90 default field capacity, metres

# Penman constants. Gas constant and gravity are this planet's, from the run
# namelist and config; the rest are properties of water and air.
GASCON = 287.017        # J/kg/K, from the model's own PLANET_NL
CP_AIR = 1005.0         # J/kg/K
KARMAN = 0.4
# Quadrature points around one diurnal cycle. 24 is far more than the smoothness
# of a sinusoid needs; the cost is nothing and it removes the point count as a
# thing to wonder about.
DIURNAL_POINTS = 24
# ExoPlaSim's own stability constants, fluxmod.f90:27-29, ECHAM Report 218.
VDIFF_B = VDIFF_D = 5.0
STABILITY_ITERATIONS = 5
Z0_WATER = 1.5e-4       # m, open-water roughness length
WATER_ALBEDO = 0.06     # the export's own value for the water rock class
SIGMA_LOWEST = 0.9828   # lowest model level


def saturation_vapour_pressure(temp_k):
    """Buck/Tetens saturation vapour pressure over water, Pa."""
    return 610.94 * np.exp(17.625 * (temp_k - 273.15) / (temp_k - 30.11))


def reference_level_air(climatology, bin_index=None):
    """Temperature, humidity, wind and pressure of the air Penman reads.

    `bin_index` takes ONE of the climatology's time bins instead of the annual
    mean, which is what a seasonally resolved open-water evaporation is
    evaluated on. It changes nothing about the annual path: with `bin_index`
    absent this returns exactly the weighted annual mean it always did, and the
    bin weights are the ones `lib/climatology.py` recovers, so summing the bins
    under those weights returns the annual answer rather than approximating it.

    Returns `(t_air, q_air, wind, p_air)`, ALL FOUR AT THE LOWEST MODEL LEVEL,
    which `SIGMA_LOWEST` puts of order 300 m up.

    **One level, because the equation has one reference height.** This used to
    return only humidity and wind, and every caller then handed Penman the 2 m
    `tas` alongside them. So the saturation term was evaluated low and warm
    while the actual vapour pressure was taken high and dry, which inflates the
    deficit that drives the aerodynamic term. Measured on this climatology, the
    2 m air is 1.85 K warmer than the lowest level over land, worth 13.7% on
    `e_s` -- one-signed, and in the same direction as the sub-grid dry-column
    error, so the two compounded rather than offsetting.

    The lowest model level is the right height rather than merely a consistent
    one. The transfer coefficient this scheme uses is derived over `z_ref`, so
    that is the height the bulk formula's humidity difference belongs at, which
    is also how the model computes its own surface fluxes. And a sub-grid lake
    does not sit under the 2 m air the model reports, which is a diagnostic of
    the dry ground beside it; 300 m up is the nearest thing to the regional air
    a lake would actually see.

    The measurement that says this is a fix and not a preference is the ocean
    validation, which is the one place the answer is already known: mixing the
    levels put Penman 8.45% above the model's own open-water evaporation, and
    reading one level puts it 3.28% below.

    Everything comes from the binned climatology with no correction. The
    corrections that used to be here are gone rather than improved: under
    PlaSim's low-I/O output path the binned wind was wrong twice over, a corrupt
    first record per orbit and a partial vector cancellation from accumulating
    components before averaging them, so this read wind from the snapshot
    product and dropped the bad record from humidity. Runs now set `NLOWIO = 0`
    and both defects vanish together -- the first-bin wind ratio goes from 9.0
    to 1.03 and binned `spd` from 1.14x the mean of instantaneous speeds to
    1.002x. The binned product is now strictly better than the snapshot one, at
    182 samples an orbit against 32, and `require_clean_io` refuses a
    climatology that predates the change rather than correcting it.
    """
    require_clean_io(climatology)
    with Dataset(climatology) as ds:
        centres = np.asarray(ds["time"][:])
        if bin_index is None:
            def take(name):
                return weighted_annual_mean(np.asarray(ds[name][:]), centres)
        else:
            k = int(bin_index)
            if not 0 <= k < centres.size:
                raise IndexError(f"bin {k} outside the climatology's {centres.size}")

            def take(name):
                return np.asarray(ds[name][k])
        t_air = take("ta")[-1]
        q_air = take("hus")[-1]
        wind = take("spd")[-1]
        p_air = take("ps") * 100.0 * SIGMA_LOWEST
    return t_air, q_air, wind, p_air


def penman_open_water(t_air, q_air, wind, p_air, rss, rls, land_albedo,
                      gravity, diurnal_range=None, column_relative_humidity=None,
                      cfg=None):
    """Penman open-water evaporation, m/s.

    Combines the energy budget with an aerodynamic term, which is what makes it
    right for a lake in a dry surrounding: a moisture-limited land cell can be
    energy-starved while the air above it stays thirsty, and the aerodynamic
    term is what carries that advected demand.

    `t_air`, `q_air`, `wind` and `p_air` must all be at ONE height, and
    `reference_level_air` is what supplies them. The radiation terms are surface
    quantities and stay so, which is what Penman is: a surface energy budget
    closed by an aerodynamic term measured to a reference level.

    Two corrections matter for a lake specifically. Net radiation is recomputed
    with water's albedo rather than the substrate's, since a lake absorbs far
    more shortwave than the 0.22-0.50 ground around it. And the roughness length
    is water's, not land's, which lowers the transfer coefficient.

    `column_relative_humidity` floors the ambient vapour pressure at that
    fraction of saturation. It is OFF in every product and exists to bracket the
    sub-grid moistening a lake does to the air over it, which is real and which
    nothing here can compute: the internal boundary layer over a lake tens of
    kilometres across is deeper than this reference level, so the air a lake
    sees is wetter than the cell mean the model reports. Correcting it needs a
    fetch-dependent boundary-layer model, so it is bracketed rather than
    applied -- a floor chosen to move the answer would be a knob.

    `cfg` is the parsed `config/planet.yaml`, for the reference height's gas
    constant and gravity. It is optional only so the existing callers that pass
    `gravity` alone keep working; omitting it re-reads the same file.
    """
    lam = 2.501e6 - 2370.0 * (t_air - 273.15)          # latent heat, J/kg
    # DIURNAL INTEGRATION. Saturation vapour pressure is convex in temperature at
    # about 6.7% per kelvin, so the daily MEAN of e_s exceeds e_s of the daily
    # mean, and evaluating Penman at a 12-bin mean understates the vapour
    # pressure deficit that drives it. Jensen's inequality, and it is one-signed.
    #
    # The model carries `maxt` and `mint`, so the range is measured rather than
    # assumed: 8.08 K land-mean and 12.32 K at the ninetieth percentile. A
    # sinusoid between them is integrated by Gauss-Legendre over one cycle.
    # `delta` is the derivative of the same curve and takes the same treatment,
    # because using a diurnal-mean e_s with a point-evaluated slope would be
    # inconsistent in the direction that flatters the answer.
    if diurnal_range is None:
        es_a = saturation_vapour_pressure(t_air)
        delta = es_a * 17.625 * 243.04 / (t_air - 30.11) ** 2
    else:
        phase = 2.0 * np.pi * (np.arange(DIURNAL_POINTS) + 0.5) / DIURNAL_POINTS
        amp = np.maximum(diurnal_range, 0.0) / 2.0
        es_a = np.zeros_like(t_air)
        delta = np.zeros_like(t_air)
        for ph in phase:
            t = t_air + amp * np.sin(ph)
            e = saturation_vapour_pressure(t)
            es_a += e
            delta += e * 17.625 * 243.04 / (t - 30.11) ** 2
        es_a /= DIURNAL_POINTS
        delta /= DIURNAL_POINTS
    gamma = CP_AIR * p_air / (0.622 * lam)
    e_air = q_air * p_air / (0.622 + 0.378 * q_air)
    if column_relative_humidity is not None:
        # The bracket, not a correction. See the docstring.
        e_air = np.maximum(e_air, column_relative_humidity * es_a)

    # Shortwave reaching the surface, backed out of the net and the albedo the
    # run was actually given, then re-absorbed at water's albedo.
    sw_down = rss / np.maximum(1.0 - land_albedo, 1e-3)
    net_radiation = (1.0 - WATER_ALBEDO) * sw_down + rls

    # ONE derivation of the height the transfer coefficient is taken over.
    # `exoplasim/scripts/build_surface_roughness.py` needs the same z_ref and
    # carried the EARTH-gravity answer as a literal instead; the expression is
    # now `lib/lapse.py:reference_height_m`, whose R comes from the configured
    # composition rather than being retyped from the run namelist.
    z_ref = reference_height_m(t_air, cfg, SIGMA_LOWEST)
    ce_neutral = KARMAN ** 2 / np.log(np.maximum(z_ref, 1.0) / Z0_WATER) ** 2
    rho = p_air / (GASCON * t_air)
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
        t_surface = t_air + r_a * (rn - latent) / (rho * CP_AIR)
        # Virtual temperature difference, surface minus air; positive is unstable.
        d_theta_v = t_surface * (1.0 + 0.6078 * q_air) - t_air
        ri = np.clip(-gravity * z_ref * d_theta_v / (u ** 2 * t_air), -5.0, 5.0)
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
    p_mm = float(penman[ocean].mean()) * day
    m_mm = float(evap[ocean].mean()) * day
    return {"penman_mm_per_day": round(p_mm, 4),
            "model_mm_per_day": round(m_mm, 4),
            "ratio": round(p_mm / m_mm, 4) if m_mm else None,
            "ocean_cells": int(ocean.sum()),
            "note": "Ocean cells are already open water, so this is a direct "
                    "check."}


def read_sra_field(path: Path, nlat: int, nlon: int) -> np.ndarray:
    """Read a formatted SRA surface field: one header line, then the values."""
    lines = path.read_text(encoding="ascii").splitlines()
    values = np.array(" ".join(lines[1:]).split(), dtype=np.float64)
    if values.size != nlat * nlon:
        raise RuntimeError(f"{path}: {values.size} values, expected {nlat * nlon}")
    return values.reshape(nlat, nlon)


def annual_mean(ds: Dataset, name: str) -> np.ndarray:
    """Annual mean of a (time, lat, lon) field, weighted by records per bin.

    The bins are NOT equal length, which this used to assume. pyburn splits the
    raw stream with `np.linspace(0, ntimes, nbin+1).astype(int)`, so at the 182
    records per orbit of a clean-I/O run two bins in twelve hold sixteen and the
    other ten hold fifteen. `lib/climatology.py` recovers the counts from the
    bin centres. Worth +0.16% on land P - E, which is under a basin against the
    error budget's own response -- small, and no reason for the criterion to be
    computed on a mean nobody can defend. CLIM-13.
    """
    return weighted_annual_mean(np.asarray(ds[name][:]), np.asarray(ds["time"][:]))


def seasonal_rectification(climatology, year_s: float) -> dict:
    """What clamping `P - E` per bin would add, and where that water came from.

    The catchment integral uses the ANNUAL MEAN of `P - E`, clamped at zero
    after aggregation, and that is the conserved answer rather than a choice.
    Over one annual cycle at steady state a land cell's storage returns to where
    it started, so `annual(P - E)` equals the runoff the cell generated, exactly.
    Clamping the negative bins away instead counts the wet season's supply twice:
    once as runoff, and once as the water that refilled the soil the dry season
    emptied.

    This measures that directly. `deficit` is what per-bin clamping would add;
    `storage_range` is the seasonal swing of the model's own soil water plus
    snow. Cell by cell the median ratio of the two is 0.99 -- the store IS the
    deficit, and the water per-bin clamping would add is water the model already
    spent. The identity is what makes this a check rather than a comparison: a
    bucket cannot generate negative runoff, and it cannot supply water it never
    stored.
    """
    with Dataset(climatology) as ds:
        pme = np.asarray(ds["pr"][:]) - (-np.asarray(ds["evap"][:]))
        store = np.asarray(ds["mrso"][:]) + np.asarray(ds["snd"][:])
        centres = np.asarray(ds["time"][:])
        lsm = weighted_annual_mean(np.asarray(ds["lsm"][:]), centres)
        lat = np.asarray(ds["lat"][:])
    # Per-bin duration, from the same weights: the bins are not equal length, so
    # a single `year_s / nbin` mis-times the dry-bin deficit below.
    dt = year_s * bin_weights(centres)
    land = lsm > 0.5
    w = np.where(land, np.cos(np.deg2rad(lat))[:, None] * np.ones_like(lsm), 0.0)

    def land_mean(field):
        return float((field * w).sum() / w.sum()) * 1000.0     # mm

    deficit = (np.maximum(-pme, 0.0) * dt[:, None, None]).sum(axis=0)
    swing = store.max(axis=0) - store.min(axis=0)
    ratio = np.where(deficit > 1e-6, swing / np.maximum(deficit, 1e-9), np.nan)
    usable = np.isfinite(ratio) & land
    return {
        "annual_mean_p_minus_e_mm": round(land_mean(
            weighted_annual_mean(pme, centres) * year_s), 2),
        "sum_of_positive_bins_mm": round(land_mean(
            (np.maximum(pme, 0.0) * dt[:, None, None]).sum(axis=0)), 2),
        "dry_bin_deficit_mm": round(land_mean(deficit), 2),
        "seasonal_soil_and_snow_range_mm": round(land_mean(swing), 2),
        "median_store_over_deficit": round(float(np.median(ratio[usable])), 3),
        "note": "land means per Vesper year. The soil-and-snow store's seasonal "
                "swing IS the dry-bin deficit, cell by cell, so clamping P - E "
                "per bin adds water the model drew from storage and did not run "
                "off. The annual mean is the conserved quantity and is what the "
                "criterion uses.",
    }


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
    The decode itself lives in `lib/gridding.py:coupling_cells`, which is the
    module that owns the convention; this function is an area-weighted sum and
    nothing else.
    """
    basin, row, col, area, _, nlon = coupling_cells(coupling)
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
    ap.add_argument("--for-build", default=None,
                    help="declare a DELIBERATE cross-build read of the staged "
                         "background albedo by naming the build it was staged "
                         "from. The overshoot measurement needs this once the "
                         "carved build is staged and the verdict is being "
                         "re-taken against the pre-carve one. A name from "
                         "lib/orogen.py's registry; any other build is refused.")
    # DUST-10. Dust changes the SURFACE energy balance at a lake, which is what
    # Penman reads, and the sign is not the top-of-atmosphere one: the layer
    # absorbs, so the ground loses even where the TOA gains over bright fill.
    # Off by default because the baseline climatology has no dust in it, so
    # applying this is asking a what-if rather than reporting the run.
    ap.add_argument("--dust-forcing", type=Path, default=None,
                    help="add per-cell dust perturbations to rss and rls before "
                         "Penman; see analysis/dust_surface_forcing.nc")
    # DUST-11 and HYD-14 and HYD-11 all move catchment runoff, which is the
    # DENOMINATOR of the aridity index, and all of them are priced in percent.
    # This converts a percent into basins so the three can be compared and
    # ranked against each other in the currency the carve is decided in.
    #
    # It is a SENSITIVITY, not a climate: precipitation and evaporation are left
    # where they are, so nothing here should be read as a verdict. The verdict on
    # a dust climate is the one taken on that climate's own climatology.
    ap.add_argument("--runoff-scale", type=float, default=1.0,
                    help="scale catchment runoff by this factor before the "
                         "aridity index; a sensitivity, not a climate")
    # GW-14. The criterion's supply term was SURFACE runoff, because that was
    # the only supply a surface-only model had. It is not the only one now.
    #
    # 605 basins here have `max(P - E, 0)` of exactly zero over their whole
    # catchment, which sends the aridity index to infinity and makes them
    # incapable of carving whatever their geometry -- and that was never a
    # statement about geometry, it was a division by a supply term missing a
    # term. Handed a water table, the supply becomes
    #
    #     r_eff = (surface_runoff * C + Qg) / C
    #
    # and the infinite branch disappears for any basin that receives water by
    # some route. It is KEPT for `r_eff <= 0`, because a basin that genuinely
    # receives nothing -- including one exporting its entire recharge
    # underground -- still never overflows and still never incises, which is
    # correct physics rather than an artifact.
    #
    # OFF unless a field is given, and the reduction identity is the test:
    # without one, every verdict is bit-identical to the surface-only answer.
    ap.add_argument("--groundwater", type=Path, default=None,
                    help="water_table.nc, to add its per-basin net groundwater "
                         "exchange to the criterion's supply term. Omitted, the "
                         "criterion is exactly the surface-only one")
    args = ap.parse_args()

    _build_data = component_data("hydrography", strict=True)
    # LOADED BEFORE THE DEFAULTS, because one of them needs it: the coupling
    # matrix is per rung and `coupling_path` reads `model.resolution`.
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if args.coupling is None:
        args.coupling = coupling_path(_build_data, config)
    if args.basins is None:
        args.basins = _build_data / "basins.nc"
    if args.climatology is None:
        args.climatology = climatology_path()
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
        rss = annual_mean(ds, "rss")
        rls = annual_mean(ds, "rls")
        diurnal = annual_mean(ds, "maxt") - annual_mean(ds, "mint")
    t_air, q_air, wind, p_air = reference_level_air(args.climatology)

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
    # Resolved through the one door rather than by rebuilding the path from
    # `model.resolution`: that path is keyed by the RUNG alone while
    # `surface_albedo` rewrites it per BUILD, so it cannot say which build's
    # field is in it. `--for-build` is how the OVERSHOOT measurement declares
    # the cross-build read it makes on purpose -- the staged albedo is the right
    # one while the carved build is staged and the wrong one for anything re-run
    # on the pre-carve build afterwards -- and it declares it by NAMING the
    # build, which admits that one and refuses every other. world-z7bu, world-xgtj,
    # CLAUDE.md rule 5.
    staged_albedo = staged_surface_field(174, config, for_build=args.for_build)
    land_albedo = read_sra_field(PROJECT_ROOT / staged_albedo["path"],
                                 *p_air.shape)
    gravity = float(config["planet"]["gravity_m_s2"])
    penman = penman_open_water(t_air, q_air, wind, p_air, rss, rls, land_albedo,
                               gravity, diurnal_range=diurnal, cfg=config)
    ocean_validation = validate_over_ocean(penman, evap, lsm)
    # The sub-grid dry column, bracketed rather than applied. A lake moistens
    # the air over itself and the cell mean does not know it; how much is a
    # fetch problem this project cannot solve, so what is reported is what a
    # floor on ambient humidity would be worth. HYD-11.
    land = lsm > 0.5
    dry_column_bracket = {
        f"rh_floor_{int(rh * 100)}": round(float(np.mean(
            penman_open_water(t_air, q_air, wind, p_air, rss, rls, land_albedo,
                              gravity, diurnal_range=diurnal,
                              column_relative_humidity=rh, cfg=config)[land])) * 86400.0 * 1000.0, 4)
        for rh in (0.7, 0.8, 0.9)}
    dry_column_bracket["as_computed"] = round(
        float(np.mean(penman[land])) * 86400.0 * 1000.0, 4)
    # THERE IS NO FLOOR AT THE LAND RATE, AND THERE MUST NOT BE ONE.
    #
    # `penman = max(penman, evap)` stood here. Its reason was that a saturated
    # surface cannot evaporate less than the moisture-limited ground beside it,
    # so a Penman below the model's own `evap` had to be Penman failing. That is
    # not so, and the reason is roughness: this run's land carries a roughness
    # field with a median z0 of 0.521 m against open water's 1.5e-4, which is a
    # transfer coefficient 6.4 times larger. A smooth lake in a rough, wet,
    # vegetated landscape genuinely evaporates LESS than the land around it, and
    # the model's own numbers say so -- the floor bound on 60.2% of land cells,
    # and where it bound the mean soil wetness was 0.597 against 0.121 where it
    # did not. It fired on 1,324 of the 1,347 cells at full wetness. It was not
    # catching a failure; it was catching the wet regime.
    #
    # What it cost is that the lake estimate simply became the land estimate
    # across whole catchments: on 1,252 of the 2,543 basins with catchment
    # runoff the aridity index came out at exactly -1.0, so the criterion
    # reduced to `crit > -1` and the basin carved regardless of its climate.
    # That is 73% of the overflowing set decided by a clamp.
    #
    # The nesting it was defending is itself the error. `wet` bounds toward
    # carving only where the ground is moisture-limited; where the ground is wet
    # the roughness contrast reverses the ordering, and the two estimates are
    # not nested at all. The note on `wet` above says so now.
    #
    # The check that settles it is the ocean, where the model uses water's own
    # roughness and the answer is already known: `validate_over_ocean` puts this
    # scheme within a few percent there, in exactly the configuration a lake is
    # in. A scheme that reproduces open-water evaporation over open water does
    # not need rescuing over land.

    # Catchment runoff is P - E, not `mrro`. `mrro` is river-routed net
    # divergence rather than local generation -- landmod.f90's roffstep calls
    # mkradv, which advects runoff downhill and modifies its argument in place --
    # so integrating it over one of our catchments measures ExoPlaSim's routing
    # rather than the water arriving at our sink. `mrro` is kept in the field set
    # so the two remain comparable in the report.
    fields = {"pr": pr, "evap": evap,
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
    if args.runoff_scale <= 0.0:
        raise SystemExit("--runoff-scale must be positive")
    runoff = np.maximum(means["runoff"], 0.0) * to_km_per_year * args.runoff_scale
    if args.runoff_scale != 1.0:
        print(f"  catchment runoff scaled by {args.runoff_scale}: a sensitivity "
              "on the criterion's denominator, not a climate")

    # GW-14: the groundwater half of the supply, as a depth over the catchment
    # so it adds to a runoff depth. This is the SUPPLY the basin receives, so it
    # feeds the lake solver as well as the index -- they take the same quantity
    # and separating them would be two answers to one question.
    groundwater_depth = np.zeros_like(runoff)
    groundwater_info = None
    if args.groundwater is not None:
        with Dataset(args.groundwater) as ds:
            if getattr(ds, "terrain_hash", None) != basins.terrain_hash:
                raise SystemExit(
                    f"{args.groundwater} was built from a different terrain "
                    "than basins.nc; the basin indices do not correspond")
            qg = np.asarray(ds["basin_groundwater_m3_s"][:]).astype(np.float64)
        if qg.size != n:
            raise SystemExit(
                f"{args.groundwater} has {qg.size} basins and basins.nc has {n}")
        with np.errstate(divide="ignore", invalid="ignore"):
            groundwater_depth = np.where(
                catch_area > 0, qg / (catch_area * 1e6), 0.0) * to_km_per_year
        import hashlib
        _h = hashlib.sha256()
        with open(args.groundwater, "rb") as _fh:
            while _chunk := _fh.read(1 << 20):
                _h.update(_chunk)
        groundwater_info = {
            "file": rel(args.groundwater),
            "sha256": _h.hexdigest(),
            "basins_gaining": int((qg > 0).sum()),
            "basins_losing": int((qg < 0).sum()),
            "note": ("added to the criterion's supply term as a catchment "
                     "depth; the lake solver takes the same supply"),
        }
        print(f"  groundwater supply from {rel(args.groundwater)}: "
              f"{int((qg > 0).sum())} basins gain, {int((qg < 0).sum())} lose")
    # Clamped for the same reason the surface half is: a basin cannot deliver a
    # negative amount of water to its own sink.
    runoff = np.maximum(runoff + groundwater_depth, 0.0)
    runoff_mrro = means["mrro"] * to_km_per_year
    precip = means["pr"] * to_km_per_year

    # THE INTEGRAL, which is not the mean reported below, and the distinction
    # has already cost one wrong number. `catchment_mean_mm_per_year` is an
    # UNWEIGHTED mean over basins, so multiplying it by the total catchment area
    # weights a 400 km2 basin the same as a 4,000,000 km2 one. Doing that gave
    # 97% of all land runoff arriving off 76% of the land, which is not credible
    # on the arid share of a planet; `notes/audits/unpriced-terms.md` records it.
    #
    # The area-weighted integral is the water the endorheic system receives, and
    # it is the quantity the absent lake moisture source has to be checked
    # against: at steady state a basin's inflow leaves as evaporation or as
    # spill, and the model provides neither over land.
    earth_years_per_orbit = year_s / (365.25 * 86400.0)
    signed = means["runoff"] * to_km_per_year
    inflow_unclamped = float(np.nansum(catch_area * signed)) / earth_years_per_orbit
    inflow_per_basin = (catch_area * np.maximum(np.nan_to_num(signed), 0.0)
                        / earth_years_per_orbit)
    inflow_clamped = float(inflow_per_basin.sum())
    catchment_total_km2 = float(np.nansum(np.where(np.isnan(signed), 0.0, catch_area)))

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
            # PRE-CARVE, and deliberately not the `area` above. `area` zeroes
            # the basins the verdict carves, which describes the world AFTER
            # the carve; the moisture source that is missing from the climate
            # that produced this verdict is the one over the lakes standing
            # today, and an overflowing basin stands at its spill area.
            "lake_net_evaporation_km3": float(np.nansum(
                lakes["area_km2"] * (evapo - precip))),
            "spill_to_ocean_km3": float(np.nansum(
                lakes["overflow_km3_per_year"][basins.spill_target < 0])),
            # The moisture-source difference, which is NOT the line above.
            # What the atmosphere gains or loses by the lake existing is the
            # lake's evaporation against what the model already evaporates on
            # the same ground, not the lake's net water demand. On this world
            # the sign is not obvious and is measured rather than argued: land
            # here is several times rougher than open water, so the model's
            # own evaporation exceeds Penman on most land cells.
            "lake_minus_model_evaporation_km3": float(np.nansum(
                lakes["area_km2"] * (means["penman"] - means["evap"])
                * to_km_per_year)),
            "runoff_over_lake_footprint_km3": float(np.nansum(
                lakes["area_km2"] * runoff)),
            "dry_survivors": int((survives & (runoff <= 0)).sum()),
            "wet_survivors": int((survives & (area > 1.0)).sum()),
            "converged": lakes["converged"],
        }

    # The split that turns the integral above into a CHECK. A basin that does
    # not overflow retains everything its catchment delivers, and at steady
    # state retention leaves as evaporation over and above the precipitation
    # falling on the lake itself -- which the model already has, because the
    # bucket on the lake fraction evaporates local rainfall at the potential
    # rate. So retention, not gross lake evaporation, is the moisture source
    # ExoPlaSim is missing over land, and an independent estimate of it from
    # lake area times a Penman rate has to land inside this.
    retained = float(inflow_per_basin[~results["penman"]["carve"]].sum())
    spilling = float(inflow_per_basin[results["penman"]["carve"]].sum())

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
        "staged_background_albedo": staged_albedo,
        "groundwater": groundwater_info,
        "basins": n,
        "note": ("'penman' is the primary estimate: the Penman combination "
                 "equation with water's albedo and roughness, evaluated wholly "
                 "at the lowest model level. Its error against the model over "
                 "ocean cells, which are already open water, is in "
                 "penman_ocean_validation and is COMPUTED, never quoted. 'wet' "
                 "uses the moisture-limited land evaporation; it is a "
                 "sensitivity and not a bound, because this world's land is "
                 "several times rougher than open water and so evaporates more "
                 "than a lake wherever the ground is wet."),
        "runoff_source": {
            "used": "p_minus_e",
            "note": "mrro is river-routed net divergence, not local runoff "
                    "generation; integrating it over a catchment measures "
                    "ExoPlaSim's routing rather than inflow to our sink.",
            "catchment_mean_mm_per_year": {
                "per": "VESPER year, the orbital period, and an UNWEIGHTED mean "
                       "over basins. It is a per-basin typical depth and NOT an "
                       "integral: see endorheic_inflow below for the water the "
                       "endorheic system actually receives.",
                "p_minus_e": round(float(np.nanmean(runoff)) * 1e6, 3),
                "mrro": round(float(np.nanmean(runoff_mrro)) * 1e6, 3),
            },
            "endorheic_inflow": {
                "what": "area-weighted integral of catchment P - E over every "
                        "basin's catchment, in EARTH years. This is the water "
                        "the endorheic system receives, and at steady state it "
                        "leaves as lake evaporation or as spill.",
                "catchment_km2": round(catchment_total_km2, 1),
                "km3_per_earth_year": {
                    "unclamped": round(inflow_unclamped, 1),
                    "clamped": round(inflow_clamped, 1),
                    "retained_by_surviving_basins": round(retained, 1),
                    "delivered_to_overflowing_basins": round(spilling, 1),
                },
                "lake_net_evaporation_km3_per_earth_year": round(
                    results["penman"]["lake_net_evaporation_km3"]
                    / earth_years_per_orbit, 1),
                "spill_to_ocean_km3_per_earth_year": round(
                    results["penman"]["spill_to_ocean_km3"]
                    / earth_years_per_orbit, 1),
                "runoff_over_lake_footprint_km3_per_earth_year": round(
                    results["penman"]["runoff_over_lake_footprint_km3"]
                    / earth_years_per_orbit, 1),
                "identity": "inflow(clamped) = lake net evaporation + runoff "
                            "generated over the lake footprint + spill to "
                            "ocean, at steady state. The middle term is there "
                            "because a basin's catchment includes its own bed, "
                            "so the solver's demand carries it; the three "
                            "close the budget and that is what makes the "
                            "integral usable.",
                "moisture_source": "THE TERM THE CLIMATE MODEL LACKS is "
                                   "lake_minus_model_evaporation, not the net "
                                   "evaporation above: what the atmosphere "
                                   "gains from a lake is its evaporation "
                                   "against what the model already evaporates "
                                   "on that ground. It is signed, and on this "
                                   "world the sign is not obvious, because "
                                   "land roughness puts the model's own "
                                   "evaporation above Penman on most land "
                                   "cells. notes/audits/unpriced-terms.md.",
                "lake_minus_model_evaporation_km3_per_earth_year": round(
                    results["penman"]["lake_minus_model_evaporation_km3"]
                    / earth_years_per_orbit, 1),
                "retention": "the split is under the PENMAN bound and moves "
                             "with the verdict. Retention is the moisture "
                             "source the model has no way to return to the "
                             "atmosphere over land: a surviving basin's inflow "
                             "leaves as evaporation, and in ExoPlaSim it leaves "
                             "as runoff to the sea instead. "
                             "notes/audits/unpriced-terms.md finding 1.",
                "mm_per_earth_year_over_catchment": {
                    "unclamped": round(
                        inflow_unclamped / max(catchment_total_km2, 1e-9) * 1e6, 3),
                    "clamped": round(
                        inflow_clamped / max(catchment_total_km2, 1e-9) * 1e6, 3),
                },
                "clamp": "clamped is the delivery the criterion uses: a "
                         "catchment whose evaporation exceeds its precipitation "
                         "delivers nothing to its sink, not a negative amount. "
                         "The gap between the two is what the drying catchments "
                         "would have subtracted, and it is not a correction to "
                         "the delivery.",
            },
            "basins_with_no_runoff": {
                "p_minus_e": int((runoff <= 0).sum()),
                "mrro": int((runoff_mrro <= 0).sum()),
            },
            "clamped": "after aggregation over the catchment, which is where "
                       "the noise on a non-negative quantity has averaged down "
                       "furthest. A catchment cannot deliver a negative amount "
                       "and, at steady state, cannot deliver more than its "
                       "annual P - E either.",
            "seasonal_rectification": seasonal_rectification(
                args.climatology, year_s),
        },
        "penman_ocean_validation": ocean_validation,
        # What was done to the climatology before the verdict was taken. Null
        # under both keys means the run's own numbers, unperturbed.
        "perturbations": {
            "dust_surface_forcing": dust_note,
            "runoff_scale": args.runoff_scale,
            "note": "A perturbed verdict is a SENSITIVITY. The verdict on a "
                    "climate that contains dust is the one taken on that "
                    "climate's own climatology, and the two must not be added.",
        },
        "penman_dry_column_bracket_mm_per_day": dry_column_bracket,
        "penman_land_evaporation_ordering": {
            "land_cells_where_model_evap_exceeds_penman": int(
                (lsm > 0.5).sum() and ((lsm > 0.5) & (evap > penman)).sum()),
            "land_cells": int((lsm > 0.5).sum()),
            "note": "NOT an error and NOT floored. This world's land roughness "
                    "gives a transfer coefficient several times water's, so a "
                    "smooth lake in a rough wet landscape evaporates less than "
                    "the ground around it. A floor at the land rate stood here "
                    "and decided 73% of the overflowing basins by clamp; see "
                    "the note in main().",
        },
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
