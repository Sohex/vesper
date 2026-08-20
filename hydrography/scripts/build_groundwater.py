"""The steady-state water table under a real climate, and what it moves.

`groundwater.py` owns the discretisation and the solve; this supplies the
forcing, writes the products, and takes the measurement the carve criterion is
waiting for. One product, one measurement:

    data/<build>/water_table.nc     per region: head, depth, seepage, transmissivity
                                    per basin: net groundwater exchange

    python hydrography/scripts/build_groundwater.py

RECHARGE IS P-E, the same field `surface_water.py` runs the lakes on, read
through the same coupling convention. At steady state that is the right
quantity and not an approximation to one: whatever falls on land and does not
evaporate has to leave, and this model decides which of the two doors it leaves
by. A cell whose water table reaches the surface returns its recharge as
seepage, which IS the surface runoff the existing balance already routes; a cell
whose water table is below the surface passes some of it sideways instead. That
is the whole difference between this and the surface-only balance, and it is why
zero permeability has to reproduce that balance exactly.

WHAT THIS MEASURES AND DOES NOT DECIDE. `Qg`, the net groundwater exchange per
basin, is the term Fan (2019) writes into the catchment budget and every water
balance in this project currently drops. It is reported here, per basin, with
the count of carve verdicts it would move. It is NOT wired into
`carve_verdict.py`: that is a decision about loop A and it is taken by a person
reading this measurement, not by this script.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import ANALYSIS, PROJECT_ROOT  # noqa: E402

import groundwater as gw  # noqa: E402
import lake_balance as lb  # noqa: E402
import surface_water as sw  # noqa: E402
from orogen import LAND, Export  # noqa: E402
from paths import rel  # noqa: E402

import builds  # noqa: E402
import orbit  # noqa: E402

CONFIG = Path(__file__).resolve().parents[1] / "config" / "groundwater.yaml"


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def rock_permeability(export: Export, cfg: dict, sigma: float = 0.0):
    """Intrinsic permeability per region, m2, and the unassigned mask.

    Read through the export's own `rockClasses` legend rather than an integer
    table held here, for the reason `lib/surface_classes.py` gives: a class
    inserted upstream would silently shift every code below it, and a consumer
    holding its own copy keeps working while pointing at the wrong rock.

    `sigma` shifts every class by that many of ITS OWN standard deviations, so
    +1 and -1 are the bracket Gleeson's spread implies rather than a uniform
    factor. The classes do not share a sigma and must not be shifted as if they
    did: carbonate's is 1.5 and siliciclastic's is 2.5.
    """
    legend = {rc["code"]: int(rc["id"])
              for rc in export.manifest["lithology"]["rockClasses"]}
    hydro = cfg["hydrolithology"]
    mapping = cfg["rock_class_hydrolithology"]

    unknown = set(mapping) - set(legend)
    if unknown:
        raise SystemExit(
            f"groundwater.yaml names rock classes this export does not have: "
            f"{sorted(unknown)}")
    missing = set(legend) - set(mapping)
    if missing:
        raise SystemExit(
            f"this export has rock classes groundwater.yaml does not map: "
            f"{sorted(missing)}. Every class needs a hydrolithology or an "
            f"explicit null; there is no default.")

    n_class = max(legend.values()) + 1
    log_k = np.full(n_class, np.nan)
    for code, name in mapping.items():
        if name is None:
            continue
        h = hydro[name]
        log_k[legend[code]] = float(h["log10_k_m2"]) + sigma * float(h["log10_k_sigma"])

    rock = export.substrate_class.astype(np.int64)
    k = np.where(np.isnan(log_k[rock]), 0.0, 10.0 ** np.nan_to_num(log_k[rock]))
    unassigned = np.isnan(log_k[rock])
    return k, unassigned, legend


def recharge_field(export: Export, config: dict, clim_path: Path):
    """P-E per region, m/s, on the mesh. The same field the lakes are solved on.

    Goes through `surface_water.climate_fields` and
    `surface_water.region_grid_cells` rather than repeating either. The second
    of those is the export/ExoPlaSim grid join, which `CLAUDE.md` rule 3 governs
    and which this project has got wrong on three separate scripts; there is one
    copy of it in `lib/gridding.py` and this is not a fourth.
    """
    sw._CLIM_FILE = clim_path
    lat, lon, runoff, precip, evaporation, model_runoff, lsm = sw.climate_fields(config)
    row, col = sw.region_grid_cells(export, lat)
    return runoff[row, col], (lat, lon, runoff, precip, evaporation, lsm)


def basin_exchange(terminal, seepage_m3_s, supply_m3_s, n_basins):
    """Net groundwater import per basin, m3/s. Fan (2019)'s `Qg`, signed.

    A basin's catchment receives its own recharge and returns some of it as
    seepage. The difference is what crossed the catchment boundary underground:
    positive where the basin gains water the surface balance never sees,
    negative where it loses water the surface balance credits it with.
    """
    gain = np.zeros(n_basins)
    give = np.zeros(n_basins)
    sel = terminal >= 0
    np.add.at(gain, terminal[sel], seepage_m3_s[sel])
    np.add.at(give, terminal[sel], supply_m3_s[sel])
    return gain - give, gain, give


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--data", type=Path, default=None,
                    help="hydrography products for this build; defaults to "
                         "data/<source_build>/")
    ap.add_argument("--climatology", type=Path, default=None)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--sigma", type=float, default=0.0,
                    help="shift every hydrolithology by this many of its own "
                         "Gleeson standard deviations; -1 and +1 are the bracket")
    ap.add_argument("--reduction-test", action="store_true",
                    help="run at zero permeability and check the surface-only "
                         "identity instead of solving the real case")
    ap.add_argument("--divide-test", action="store_true",
                    help="run at uniform permeability with a terrain-following "
                         "table and check the catchments against the surface ones")
    ap.add_argument("--max-outer", type=int, default=60)
    ap.add_argument("--uniqueness-check", action="store_true",
                    help="re-solve from the opposite initial active set and "
                         "require the identical head field. The problem is a "
                         "linear complementarity problem with a symmetric "
                         "positive definite matrix, so its solution is unique "
                         "and this is an identity rather than a comparison")
    args = ap.parse_args()

    config = yaml.safe_load((PROJECT_ROOT / "config/planet.yaml").read_text())
    cfg = yaml.safe_load(CONFIG.read_text())

    data = args.data if args.data is not None else builds.component_data(
        "hydrography", config, strict=True)
    if args.climatology is not None:
        clim = args.climatology.resolve()
    else:
        declared = config.get("baseline_climatology")
        if not declared:
            raise SystemExit(
                "config/planet.yaml has no `baseline_climatology`. Name one "
                "there or pass --climatology; there is deliberately no fallback.")
        clim = (PROJECT_ROOT / declared).resolve()
        if not clim.is_file():
            raise SystemExit(f"config names {declared}, which does not exist")

    build = builds.build_root(config)
    export = Export(builds.mesh_export(config))
    sw._DATA = data
    land = export.surface_class == LAND

    print(f"build {build.name}, {export.n_regions:,} regions")
    print("reconstructing the Voronoi dual")
    geom = gw.Geometry(export)
    closes = geom.closes_on_sphere()
    print(f"  {geom.src.size:,} faces, area closes on the sphere to {closes:.8f}")
    if abs(closes - 1.0) > 1e-6:
        raise SystemExit(
            f"the Voronoi dual does not tile the sphere ({closes:.8f}); the "
            "geometry is wrong and nothing below it is worth computing")

    # -- subsurface --------------------------------------------------------
    k_m2, unassigned, legend = rock_permeability(export, cfg, args.sigma)
    fluid = cfg["fluid"]
    gravity = float(config["planet"]["gravity_m_s2"])
    k0 = gw.conductivity(k_m2, float(fluid["density_kg_m3"]), gravity,
                         float(fluid["dynamic_viscosity_pa_s"]))
    thickness_m = float(cfg["aquifer"]["thickness_m"])

    policy = str(cfg["unassigned"]["policy"])
    if policy not in ("exclude", "impermeable"):
        raise SystemExit(f"unassigned.policy {policy!r} is not one this solver has")
    conductive = land & ~unassigned if policy == "exclude" else land.copy()
    if policy == "impermeable":
        k0 = np.where(unassigned, 0.0, k0)

    print(f"  gravity {gravity} m/s2 -> conductivity is "
          f"{gravity / 9.80665:.3f}x Earth's for the same permeability")
    print(f"  unassigned lithology on {int((unassigned & land).sum()):,} land "
          f"regions ({(unassigned & land).sum() / land.sum() * 100:.2f}%), "
          f"policy {policy!r}")
    live = k0[land][k0[land] > 0]
    print(f"  aquifer thickness {thickness_m:.0f} m, constant; conductivity "
          f"spans {np.log10(live.max() / live.min()):.1f} orders of magnitude "
          f"across lithologies")

    # -- forcing -----------------------------------------------------------
    print(f"reading the climatology {rel(clim)}")
    recharge, fields = recharge_field(export, config, clim)
    lat, lon, runoff_grid, precip_grid, evap_grid, lsm = fields
    recharge = np.where(land, recharge, 0.0)

    surface_m = export.elevation_km.astype(np.float64) * 1000.0

    if args.reduction_test:
        print("\nREDUCTION IDENTITY: zero permeability must reproduce the "
              "surface-only balance exactly")
        k0 = np.zeros_like(k0)
    if args.divide_test:
        # IMPOSED, not solved. The declared test is what a terrain-following
        # table does, so the table is set to the terrain rather than solved for
        # and then compared to the terrain it did not quite reproduce. Uniform
        # permeability, so nothing but the topography can steer the flow, and
        # the answer is fixed in advance: the groundwater divides must be the
        # topographic ones.
        print("\nDIVIDE TEST: uniform permeability, terrain-following table")
        k0 = np.where(land, 1e-4, 0.0)
        res = gw.terrain_following(export, geom, k0_m_s=k0,
                                   thickness_m=thickness_m,
                                   surface_m=surface_m, conductive=conductive)
    else:
        print("solving the water table")
        res = gw.solve(export, geom, k0_m_s=k0, thickness_m=thickness_m,
                       recharge_m_s=recharge, surface_m=surface_m,
                       conductive=conductive, max_outer=args.max_outer)

    area_m2 = geom.volume_area_m2
    supply = recharge * area_m2
    # Closure is a statement about a solved water balance. The divide test
    # imposes a table rather than solving one, so there is no balance to close
    # and reporting a residual for it would invent a failure.
    clos = None
    if not args.divide_test:
        clos = gw.closure(res, recharge, area_m2, export)
        print(f"\nCLOSURE  recharge {clos['recharge_m3_s']:,.3f} m3/s   "
              f"discharge {clos['discharge_m3_s']:,.3f} m3/s")
        print(f"         relative residual {clos['relative_residual']:.3e} "
              f"against {clos['tolerance']:.0e}   "
              f"{'pass' if clos['passes'] else 'FAIL'}")
        print(f"         seepage {clos['seepage_m3_s']:,.1f}   "
              f"to the ocean underground {clos['to_ocean_m3_s']:,.1f}")

    depth = res["depth_m"]
    print(f"\nwater table depth on land: median {np.median(depth[land]):.2f} m, "
          f"at the surface on {(depth[land] <= 0.01).mean() * 100:.1f}% of land")

    with Dataset(data / "regions.nc") as ds:
        terminal = np.asarray(ds["terminal"][:])
    basins = lb.BasinSet(data / "basins.nc")
    qg, gain, give = basin_exchange(terminal, res["seepage_m3_s"], supply, basins.n)

    report = {
        "source_build": build.name,
        "terrain_hash": export.terrain_hash,
        "forcing": rel(clim),
        "forcing_sha256": sha256(clim),
        "config_sha256": sha256(CONFIG),
        "sigma": args.sigma,
        "unassigned_policy": policy,
        "gravity_m_s2": gravity,
        "conductivity_vs_earth": gravity / 9.80665,
        "geometry": {
            "faces": int(geom.src.size),
            "voronoi_area_over_sphere": closes,
            "export_cell_area_over_sphere": float(
                geom.volume_area_m2.sum() / (4 * np.pi * geom.radius_m ** 2)),
            "note": ("flux uses the Voronoi area its own faces bound; water "
                     "volumes use the export's cell_area so the reduction "
                     "identity against surface_water.py is exact"),
        },
        "laplace_beltrami": {
            "criterion": gw.LAPLACE_TOLERANCE,
            "relative_rms": {str(k): v for k, v in
                             gw.laplace_beltrami_error(geom).items()},
        },
        "closure": clos,
        "converged": bool(res.get("converged", True)),
        "convergence": {
            "criterion_residual": gw.RESIDUAL_TOLERANCE,
            "criterion_flip_fraction": gw.FLIP_TOLERANCE,
            "final_residual": res.get("final_residual"),
            "criterion_note": ("the water balance residual on free cells over "
                               "total land recharge, and a still active set. "
                               "Not the head step: a bar on the step can be "
                               "met by a shrinking step rather than by "
                               "convergence"),
            "trace_columns": ["pass", "free_cells", "cells_flipped",
                              "residual", "infeasible_pinned_cells"],
            "seepage_clipped_m3_s": res.get("seepage_clipped_m3_s"),
            "trace": res.get("residual_trace", []),
        },
        "outer_iterations": res.get("outer_iterations"),
        "water_table": {
            "median_depth_m_land": float(np.median(depth[land])),
            "at_surface_fraction_land": float((depth[land] <= 0.01).mean()),
            "unassigned_land_fraction": float(
                (unassigned & land).sum() / max(int(land.sum()), 1)),
            "aquifer_thickness_m": thickness_m,
            "transmissivity_note": (
                "T = K D with a constant aquifer thickness, so there is no "
                "depth floor and no cell has a transmissivity that underflows"),
        },
        "git_commit": subprocess.run(
            ["git", "-C", str(PROJECT_ROOT), "rev-parse", "HEAD"],
            capture_output=True, text=True).stdout.strip(),
    }

    # -- the declared tests ------------------------------------------------
    if args.reduction_test:
        seep = res["seepage_m3_s"]
        exact = np.array_equal(seep[land], supply[land])
        worst = float(np.abs(seep[land] - supply[land]).max())
        qg_max = float(np.abs(qg).max())
        print(f"\n  seepage equals recharge on every land cell: {exact}")
        print(f"  worst absolute difference {worst:.6e} m3/s")
        print(f"  largest |Qg| over {basins.n:,} basins {qg_max:.6e} m3/s")
        report["reduction_identity"] = {
            "criterion": "exact: seepage == recharge per cell, Qg == 0 per basin",
            "seepage_equals_recharge_exactly": bool(exact),
            "worst_absolute_difference_m3_s": worst,
            "largest_abs_basin_qg_m3_s": qg_max,
            "passes": bool(exact and qg_max == 0.0),
        }
        print("  " + ("PASS" if report["reduction_identity"]["passes"] else "FAIL"))

    if args.divide_test:
        rec = gw.groundwater_receiver(export.n_regions, res["src"], res["dst"],
                                      res["face_flux_m3_s"], land)
        gwterm = gw.trace_terminals(rec, terminal, land)
        with Dataset(data / "regions.nc") as ds:
            filled = np.asarray(ds["filled_km"][:])
        unfilled = land & (np.abs(filled - export.elevation_km) < 1e-6)
        same = gwterm == terminal
        a_all = float(area_m2[land][same[land]].sum() / area_m2[land].sum())
        a_unf = float(area_m2[unfilled][same[unfilled]].sum()
                      / area_m2[unfilled].sum())
        print(f"\n  groundwater catchment matches the surface catchment on "
              f"{a_all * 100:.2f}% of land area "
              f"(criterion {gw.DIVIDE_AGREEMENT_ALL * 100:.0f}%)")
        print(f"  outside filled depressions {a_unf * 100:.2f}% "
              f"(criterion {gw.DIVIDE_AGREEMENT_UNFILLED * 100:.0f}%)")
        report["divide_test"] = {
            "criterion_all": gw.DIVIDE_AGREEMENT_ALL,
            "criterion_unfilled": gw.DIVIDE_AGREEMENT_UNFILLED,
            "agreement_all_land_area": a_all,
            "agreement_outside_filled_depressions": a_unf,
            "passes": bool(a_all >= gw.DIVIDE_AGREEMENT_ALL
                           and a_unf >= gw.DIVIDE_AGREEMENT_UNFILLED),
        }
        print("  " + ("PASS" if report["divide_test"]["passes"] else "MISS"))

    # -- GW-4: what the term would move ------------------------------------
    if not (args.reduction_test or args.divide_test):
        report["basin_exchange"] = measure_carve_effect(
            export, data, basins, terminal, qg, give, gain, config, clim,
            fields, area_m2)

    # -- products ----------------------------------------------------------
    # UNIQUENESS, which replaces the cross-scheme check the nonlinear model
    # needed. With transmissivity independent of the head there is no second
    # algebra to compare against -- Picard and Kirchhoff both collapse onto this
    # same linear system. What linearity buys instead is stronger: the matrix is
    # symmetric positive definite, so the complementarity problem has exactly
    # ONE solution, and any two active-set trajectories must land on the
    # identical head field. Starting from every cell free rather than every cell
    # pinned walks a completely different path to it. That is an identity, it
    # can genuinely fail, and it is what would catch a bug in the release and
    # pin logic that a single run cannot see.
    if args.uniqueness_check and not res.get("converged", True) is False:
        print("\nUNIQUENESS: re-solving from the opposite initial active set")
        alt = gw.solve(export, geom, k0_m_s=k0, thickness_m=thickness_m,
                       recharge_m_s=recharge, surface_m=surface_m,
                       conductive=conductive, max_outer=args.max_outer,
                       start_all_free=True)
        d = np.abs(alt["head_m"] - res["head_m"])[conductive]
        scale = max(float(np.abs(res["head_m"][conductive]).max()), 1.0)
        rel_head = float(d.max() / scale)
        ok = bool(alt.get("converged")
                  and rel_head < gw.SCHEME_HEAD_RELATIVE)
        report["uniqueness"] = {
            "criterion_relative_head": gw.SCHEME_HEAD_RELATIVE,
            "max_absolute_head_difference_m": float(d.max()),
            "relative_head_difference": rel_head,
            "alternate_converged": bool(alt.get("converged")),
            "alternate_passes": alt.get("outer_iterations"),
            "passes": ok,
        }
        print(f"  max |h_altstart - h| {d.max():.3e} m, relative {rel_head:.3e} "
              f"against {gw.SCHEME_HEAD_RELATIVE:.0e}")
        print("  " + ("PASS" if ok else "MISS"))

    if not res.get("converged", True):
        # The report goes out; the water table does not. An unconverged head
        # field written to `data/` would be read downstream as a result, and
        # nothing about the file would say it was not one.
        ANALYSIS.mkdir(parents=True, exist_ok=True)
        (ANALYSIS / "groundwater_report.json").write_text(
            json.dumps(report, indent=2) + "\n")
        print(f"\nwrote {ANALYSIS / 'groundwater_report.json'}")
        print("NOT writing water_table.nc: the solve did not converge, and an "
              "unconverged head field is not a result.")
        return 1

    out = args.output or (data / "water_table.nc")
    if out.exists():
        raise SystemExit(
            f"{out} exists. This step writes a new artifact and never "
            "overwrites one; move it aside or pass --output.")
    with Dataset(out, "w") as ds:
        ds.createDimension("region", export.n_regions)
        ds.createDimension("basin", basins.n)
        ds.title = "Steady-state water table under the baseline climatology"
        ds.terrain_hash = export.terrain_hash
        ds.setncattr("vesper_source_build", build.name)
        ds.forcing = rel(clim)
        ds.sigma = args.sigma
        ds.caveat = (
            "An equilibrium water table, not an aquifer with a history. "
            "Permeability carries 1.5 to 2.5 orders of magnitude of Gleeson "
            "spread, so depth is a bracket; run --sigma -1 and +1 for it. "
            "Valley-to-ridge convergence is sub-grid at this mesh and absent."
        )
        for name, dat, dtype, dim, units, note in [
            ("head_m", res["head_m"], "f4", "region", "m",
             "water table elevation above sea level"),
            ("depth_m", depth, "f4", "region", "m",
             "water table below the land surface"),
            ("seepage_m3_s", res["seepage_m3_s"], "f4", "region", "m3 s-1",
             "groundwater returning to the surface"),
            ("transmissivity_m2_s", res["transmissivity_m2_s"], "f4", "region",
             "m2 s-1", "K0 f exp(-d/f), Fan et al. (2007) eq. 6"),
            ("at_surface", res["pinned"].astype(np.int8), "i1", "region", "1",
             "water table pinned at the land surface"),
            ("excluded", res["excluded"].astype(np.int8), "i1", "region", "1",
             "lithology has no assigned permeability"),
            ("basin_groundwater_m3_s", qg, "f8", "basin", "m3 s-1",
             "net groundwater import; Fan (2019) Qg, positive into the basin"),
            ("basin_seepage_m3_s", gain, "f8", "basin", "m3 s-1",
             "seepage within the catchment"),
            ("basin_recharge_m3_s", give, "f8", "basin", "m3 s-1",
             "recharge falling on the catchment"),
        ]:
            v = ds.createVariable(name, dtype, (dim,), zlib=True)
            v.units = units
            if note:
                v.long_name = note
            v[:] = dat

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    name = "groundwater_report.json"
    if args.reduction_test:
        name = "groundwater_reduction_test.json"
    elif args.divide_test:
        name = "groundwater_divide_test.json"
    elif args.sigma:
        name = f"groundwater_report_sigma{args.sigma:+g}.json"
    (ANALYSIS / name).write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {out}\nwrote {ANALYSIS / name}")
    return 0


def measure_carve_effect(export, data, basins, terminal, qg, give, gain,
                         config, clim, fields, area_m2) -> dict:
    """GW-4. How many carve verdicts would the `Qg` term move?

    A MEASUREMENT. The criterion is `carve_verdict.py`'s and is not restated
    here: a basin overflows when its aridity index falls at or below the
    geometric `critical_aridity_index` that `basins.nc` already carries. What
    changes is the denominator. The surface balance divides the lake's
    evaporative deficit by the catchment's own recharge depth; with groundwater
    the basin is fed by what actually seeps inside its catchment, which is that
    recharge plus whatever crossed the divide underground.

    Nothing here is written back into `carve_verdict.py`. Whether this term
    belongs in the criterion is a loop A decision and it is taken by reading
    this, not by running it.
    """
    lat, lon, runoff_grid, precip_grid, evap_grid, lsm = fields
    sinks = np.array([b.sink for b in export.basins])
    # Only the lake fluxes are taken from here. The catchment runoff depth is
    # recomputed on the mesh below, for BOTH arms, because the comparison has to
    # isolate the groundwater term: `per_basin_forcing` bins the catchment over
    # whole grid cells, and pairing that binned depth against a mesh-resolved
    # one would fold hydrography's known binning approximation into a number
    # that is supposed to measure groundwater and nothing else.
    _, lake_precip, lake_evap = sw.per_basin_forcing(
        basins.n, lat, lon, runoff_grid, precip_grid, evap_grid, sinks,
        export, lsm)

    year_s = orbit.orbital_year_days(config) * sw.SECONDS_PER_DAY
    to_km_yr = year_s / 1000.0

    catch_m2 = np.zeros(basins.n)
    sel = terminal >= 0
    np.add.at(catch_m2, terminal[sel], area_m2[sel])
    ok = catch_m2 > 0

    # Depths, m/s, over the catchment. The surface number is the coupling
    # matrix's catchment mean; the groundwater number is what seeped there.
    surface_depth = np.where(ok, give / np.maximum(catch_m2, 1e-30), 0.0)
    seeped_depth = np.where(ok, gain / np.maximum(catch_m2, 1e-30), 0.0)

    with Dataset(data / "basins.nc") as ds:
        crit = np.asarray(ds["critical_aridity_index"][:])

    deficit = (lake_evap - lake_precip) * to_km_yr        # km/yr
    def verdict(depth_m_s):
        r = depth_m_s * to_km_yr
        with np.errstate(divide="ignore", invalid="ignore"):
            idx = np.where(r > 0, deficit / np.maximum(r, 1e-30), np.inf)
        return idx, (idx <= crit)

    idx_s, carve_s = verdict(surface_depth)
    idx_g, carve_g = verdict(seeped_depth)
    scored = ok & np.isfinite(crit)
    flipped = scored & (carve_s != carve_g)

    net = qg / np.maximum(give, 1e-30)
    big = scored & (np.abs(net) > 0.10)
    return {
        "what_this_is": ("a measurement, not a change: Qg is reported per basin "
                         "and is NOT wired into carve_verdict.py"),
        "both_arms_use_mesh_resolved_catchment_depth": True,
        "why": ("the surface arm is recomputed on the mesh rather than taken "
                "from the coupling matrix, so the only difference between the "
                "two arms is the groundwater term"),
        "basins_scored": int(scored.sum()),
        "total_abs_exchange_m3_s": float(np.abs(qg).sum()),
        "total_recharge_m3_s": float(give.sum()),
        "exchange_over_recharge": float(np.abs(qg).sum() / max(give.sum(), 1e-30)),
        "basins_gaining": int((scored & (qg > 0)).sum()),
        "basins_losing": int((scored & (qg < 0)).sum()),
        "basins_shifted_over_10_percent": int(big.sum()),
        "median_abs_fractional_shift": float(np.median(np.abs(net[scored]))),
        "carve_verdicts": {
            "criterion": ("carve_verdict.py's own: overflow when the aridity "
                          "index is at or below basins.nc's "
                          "critical_aridity_index"),
            "would_carve_surface_only": int(carve_s[scored].sum()),
            "would_carve_with_groundwater": int(carve_g[scored].sum()),
            "flipped": int(flipped.sum()),
            "flipped_to_carve": int((flipped & carve_g).sum()),
            "flipped_to_hold": int((flipped & ~carve_g).sum()),
        },
    }


if __name__ == "__main__":
    raise SystemExit(main())
