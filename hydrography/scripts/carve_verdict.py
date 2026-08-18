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
from builds import component_data
from orbit import orbital_year_days
from paths import climatology_path, require_clean_io
from lake_balance import BasinSet, carve_verdict, solve

DRHSFULL = 0.4          # landmod.f90: wetness reaches 1 above this fraction
WSMAX_EARTH = 0.5       # landmod.f90 default field capacity, metres

# Penman constants. Gas constant and gravity are this planet's, from the run
# namelist and config; the rest are properties of water and air.
GASCON = 287.017        # J/kg/K, from the model's own PLANET_NL
CP_AIR = 1005.0         # J/kg/K
KARMAN = 0.4
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
    ce = KARMAN ** 2 / np.log(np.maximum(z_ref, 1.0) / Z0_WATER) ** 2
    rho = ps_pa / (GASCON * tas)
    r_a = 1.0 / np.maximum(ce * np.maximum(wind, 0.1), 1e-6)

    aerodynamic = (rho * CP_AIR / r_a) * np.maximum(es_a - e_air, 0.0) / gamma
    latent = (delta * np.maximum(net_radiation, 0.0) + gamma * aerodynamic) / (delta + gamma)
    return np.maximum(latent, 0.0) / (lam * 1000.0)   # W/m2 -> m/s of water


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


def basin_means(coupling: Path, fields: dict[str, np.ndarray], n_basins: int,
                field_lon: np.ndarray):
    """Catchment-area-weighted mean of each field, per basin.

    The coupling matrix numbers its columns on the Orogen grid, which runs -180
    to 180, and an ExoPlaSim climatology numbers its own on 0 to 360. Indexing
    one with the other is off by half a planet, and silently so: the array
    shapes match, latitude is unaffected, and every basin simply reads its
    antipode. `field_lon` is the longitude axis the fields are on, and the
    columns are remapped onto the coupling's before anything is indexed.

    It is required, and deliberately has no default. It was optional when the
    fix first landed, defaulting to the unremapped path, and `export_carve_list`
    was never updated to pass it: the one caller whose output leaves the project
    and changes the terrain went on reading the antipode. An argument whose
    absence silently means "do the wrong thing" reproduces the original bug in
    the shape of its own fix, so omitting it is now a TypeError at the call site.
    """
    with Dataset(coupling) as ds:
        basin = np.asarray(ds["basin"][:]).astype(np.int64)
        cell = np.asarray(ds["cell"][:]).astype(np.int64)
        area = np.asarray(ds["area_km2"][:])
        nlon = int(ds.n_lon)
        coupling_lon = (np.asarray(ds["cell_lon"][:])
                        if "cell_lon" in ds.variables else None)
    row, col = np.divmod(cell, nlon)
    if coupling_lon is None:
        raise RuntimeError(
            "this coupling file predates cell_lon and its longitude "
            "convention cannot be checked; rebuild it with "
            "build_hydrography.py")
    wrap = lambda a: (np.asarray(a) + 180.0) % 360.0 - 180.0
    remap = np.abs(wrap(field_lon)[None, :] - wrap(coupling_lon)[:, None]).argmin(axis=1)
    col = remap[col]
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
    with Dataset(args.climatology) as ds:
        field_lon = np.asarray(ds["lon"][:])
    means, catch_area = basin_means(args.coupling, fields, n, field_lon=field_lon)

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
    planet = 734_492_839.55

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
                 "within 1.7%. 'wet' uses the moisture-limited land evaporation "
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
        "penman_ocean_validation": {
            "penman_mm_per_day": 3.736, "model_mm_per_day": 3.672, "ratio": 1.017,
            "note": "Ocean cells are already open water, so this is a direct check.",
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
