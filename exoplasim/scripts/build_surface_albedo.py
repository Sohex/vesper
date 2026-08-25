#!/usr/bin/env python3
"""Write background land albedo from World Orogen lithology.

ExoPlaSim presets background albedo to a uniform `albland` of 0.22 and reads
codes 174 (broadband), 175 (<0.75 um) and 176 (>0.75 um) if the files exist, and
at `NSIMPLEALBEDO = 0` -- this project's setting -- the radiation reads the PAIR
and never 174. All three carry a real field here; see "The two bands" below. On
this planet the substrate is not uniform at all: area-weighted bare-rock albedo
over land is 0.315, and the spread runs 0.10 for basalt to 0.50 for evaporite.
Using 0.22 everywhere is a planetary albedo error of about 0.04, roughly
-12 W/m2 of absorbed flux before snow, cloud and vegetation feedbacks. That is a
large fraction of the entire 0.85-to-0.95 flux sweep, so it is not a detail.

`rock_albedo` in the export exists for this. Its own description calls it "a
surface boundary condition for the climate stage; vegetation and snow override
it where they exist. Orogen does not apply it to anything."

What this is and is not
-----------------------

This is *substrate* albedo, not land-surface albedo. A vegetated Earth-like
surface sits nearer 0.12 to 0.20, well below bare rock. Writing bare rock gives a
genuinely bare-rock planet, which is the honest starting point when no vegetation
model has run yet, but it will run cold relative to a vegetated world, and that
cold bias propagates into whatever vegetation model consumes the result.

`--mode` picks which risk to take:

  lithology  bare-rock albedo as exported. Physically what the surface is before
             anything grows on it. Correct for a first pass whose purpose is to
             feed a vegetation model, and biased cold.
  vegetated  the opposite physical endmember: everything that could carry
             vegetation set to `--vegetation-albedo`, with evaporite left bare
             because nothing grows on a salt pan. Land mean about 0.223. Use it
             with `lithology` to bracket the answer.
  scaled     the lithology *pattern*, rescaled so the land mean matches
             `--target-mean`. Kept for experiments, but not recommended as a
             production setting: it fixes the mean to a guess at the answer, and
             it puts every cell somewhere physically wrong, making basalt darker
             than any real rock and evaporite far too dark for a salt pan.
  modelled   the answer rather than an endmember: per-cell cover from an
             LPJ-GUESS run, blended against the lithology substrate. This is the
             mode that closes the loop, and the only one that is not an
             assumption. Needs --vegetation pointing at a run's fpc.out. Its two
             cover endmembers come from `model.tree_albedo` and
             `model.grass_albedo`, star-weighted by the same derivation as
             `model.vegetation_albedo`.
  uniform    write nothing and let ExoPlaSim default to 0.22.

The two endmembers are about 15 to 19 W/m2 apart in absorbed flux, against 21
W/m2 for the entire 0.85-to-0.95 stellar sweep that produced a 33 K range. The
land albedo is therefore not a bootstrap detail; it is comparable to the largest
forcing this project has varied on purpose, and with a positive
vegetation-albedo feedback on top it may select between distinct equilibria
rather than shifting one. Run both endmembers before trusting either.

Ocean cells are written with the water class value; ExoPlaSim computes ocean
albedo separately, so the value there is inert, but the array has to be full.

The two bands
-------------

Every surface written here carries a band pair, not just the vegetated paint.
The shape comes from a measured reflectance spectrum for the material -- rock
class by rock class from `analysis/rock_albedo_bands.json`, canopy from
`analysis/vegetation_albedo.json` -- and it is applied as a dimensionless
RATIO, so it composes with whatever level the class carries after overrides,
the derived evaporite split and lakes have moved it.

The pair is anchored on the model's own band weights, `lib/stellar`'s
reproduction of `solarini`, so that

    z1 * band1 + z2 * band2 == broadband

holds cell by cell. That identity is asserted on the three fields immediately
before they are written; it is the check that can fail, and it fails if a
repaint moves a level without moving the material with it.

`--flat-bands` writes one identical field to all three codes, which is what this
script did before the rock table had a split and is the only way to reproduce a
surface built then.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
import yaml

from _paths import CONFIG, INPUTS, PROJECT_ROOT

import climatology  # noqa: E402  from lib/, put on sys.path by _paths
import stellar
from sra import write_sra
from builds import resolution_of, grid_export, mesh_export
from gridding import land_fraction_of_class, land_weighted, region_cells
from provenance import config_stamp
from orogen import Export, LAND

# 174 broadband, 175 below 0.75 um, 176 above. With NSIMPLEALBEDO=0 the
# radiation uses the two-band pair; 174 is written too so the broadband
# diagnostic agrees rather than silently keeping 0.22.
ALBEDO_CODES = (174, 175, 176)

# The two derivations that supply the band SHAPE. Neither supplies a level:
# levels come from the export's rock table, the config overrides and the config
# cover endmembers, so a shape file and a level file can never disagree about a
# number they do not both hold.
ROCK_BANDS = PROJECT_ROOT / "analysis" / "rock_albedo_bands.json"
VEGETATION_BANDS = PROJECT_ROOT / "analysis" / "vegetation_albedo.json"

# The recombination identity's tolerance, in albedo, fixed before any field was
# written. Set from the INSTRUMENT: `.sra` is written at `%12.5f`, so a value
# reaches the model with a quantum of 5e-6 and a two-term recombination of such
# values cannot be held tighter. Twice that leaves the check able to fail on a
# construction error while ignoring the format's own rounding.
RECOMBINATION_TOLERANCE = 1.0e-5

# 212 is forest fraction. It is written here rather than left at its default
# because it is not independent of the albedo assumption: landmod blends snow
# albedo between forested and unforested endpoints by `dforest`, so leaving it at
# ExoPlaSim's uniform 0.5 asserts a half-forested planet while the background
# albedo asserts bare rock. That inconsistency darkens snow on the cold branch,
# which is exactly where the bistability question is decided.
FOREST_CODE = 212

# Forest fraction implied by each albedo mode. `modelled` is absent because it
# does not imply one: it reads tree cover per cell from LPJ-GUESS.
MODE_FOREST_FRACTION = {
    "lithology": 0.0,   # bare rock carries no canopy
    "vegetated": 0.5,   # Earth-like mixed cover where anything grows
    "scaled": 0.25,     # midpoint, matching that mode's compromise character
}

# LPJ-GUESS grass PFTs. Everything else in fpc.out except Lon, Lat, Year and
# Total is a tree, so this is the list that has to be right rather than a list of
# a dozen tree codes that would silently miss a newly added one.
GRASS_PFTS = ("C3G", "C4G")


def read_foliar_cover(path: Path) -> dict[tuple[float, float], tuple[float, float]]:
    """Tree and grass foliar projective cover per cell, from a run's fpc.out.

    Takes the last simulated year per cell, which is the equilibrium the run
    reached. Coordinates are rounded to two decimals because that is the
    precision LPJ-GUESS prints them at, coarser than the driver's four.
    """
    lines = path.read_text().splitlines()
    header = lines[0].split()
    pfts = [name for name in header[3:] if name != "Total"]
    grass = [i for i, name in enumerate(pfts) if name in GRASS_PFTS]
    trees = [i for i, name in enumerate(pfts) if name not in GRASS_PFTS]
    latest: dict[tuple[float, float], tuple[int, list[float]]] = {}
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 3 + len(pfts):
            continue
        key = (round(float(parts[0]), 2), round(float(parts[1]), 2))
        year = int(float(parts[2]))
        values = [float(v) for v in parts[3:3 + len(pfts)]]
        if key not in latest or year > latest[key][0]:
            latest[key] = (year, values)
    return {key: (sum(values[i] for i in trees), sum(values[i] for i in grass))
            for key, (_, values) in latest.items()}


def _rock_id(mesh_dir: Path, code: str) -> int:
    lit = json.loads((mesh_dir / "manifest.json").read_text(encoding="utf-8"))["lithology"]
    for r in lit["rockClasses"]:
        if r["code"] == code:
            return int(r["id"])
    raise KeyError(f"no rock class {code!r} in {mesh_dir}")


def band_shapes(rho: np.ndarray, z1: float, z2: float):
    """Dimensionless (band1, band2) multipliers for a band ratio `rho`.

    `z1*s1 + z2*s2 == 1` by construction, so multiplying any broadband albedo
    by the pair gives a pair that recombines to it in the model's own weights.
    """
    s1 = 1.0 / (z1 + z2 * rho)
    return s1, rho * s1


def rock_band_ratios(path: Path, mesh_root: Path) -> dict[int, float]:
    """Rock class id -> band2/band1, for the classes this export carries.

    Keyed by CODE in the file and resolved to this export's ids here, because
    ids move when Orogen adds a class and codes do not. A class the file does
    not cover raises rather than defaulting to 1.0: a silent 1.0 is exactly the
    spectrally flat substrate this table exists to remove.
    """
    table = json.loads(path.read_text(encoding="utf-8"))["classes"]
    lit = json.loads((mesh_root / "manifest.json").read_text(
        encoding="utf-8"))["lithology"]
    ratios = {}
    for entry in lit["rockClasses"]:
        code = entry["code"]
        if code not in table:
            raise SystemExit(
                f"{path.name} has no band ratio for rock class {code!r}, which "
                f"{mesh_root} carries. Add a proxy in "
                "analysis/rock_albedo_bands.py and re-run it; a class with no "
                "ratio would be written spectrally flat.")
        ratios[int(entry["id"])] = float(table[code]["band2_over_band1"])
    return ratios


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--mesh", type=Path, default=None,
                    help="export carrying raw/; defaults to the configured build")
    ap.add_argument("--grid", type=Path, default=None,
                    help="export whose grid to target; defaults to the config resolution")
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--mode",
                    choices=("lithology", "vegetated", "scaled", "modelled",
                             "uniform"),
                    default=None)
    ap.add_argument("--vegetation", type=Path, default=None,
                    help="fpc.out from an LPJ-GUESS run, for --mode modelled")
    ap.add_argument("--lakes", type=Path, default=None,
                    help="surface_water.nc from hydrography; paints solved lake "
                         "regions with open water's albedo before gridding, so "
                         "each cell gets an area-weighted composite")
    # Resolved from config.baseline_climatology, not hardcoded; the default
    # here named the superseded `climatology_s096` until 2026-08-17. Left
    # optional because the BOOTSTRAP run's albedo is built before any
    # climatology exists. See lib/paths.py:climatology_path.
    ap.add_argument("--climatology", type=Path, default=None,
                    help="coordinate source for --mode modelled, and the "
                         "evaporation field the derived evaporite split needs. "
                         "Must be the climatology the LPJ-GUESS driver was "
                         "built from when used for the former.")
    ap.add_argument("--tree-albedo", type=float, default=None,
                    help="albedo of full tree cover, for --mode modelled. "
                         "Defaults to model.tree_albedo in the config, derived "
                         "for THIS star by analysis/vegetation_albedo.py; 0.13 "
                         "is the Earth-Sun endmember it replaced, pass it "
                         "explicitly to reproduce the old surface")
    ap.add_argument("--grass-albedo", type=float, default=None,
                    help="albedo of full grass cover, for --mode modelled. "
                         "Defaults to model.grass_albedo in the config, same "
                         "derivation; 0.19 is the Earth-Sun endmember")
    ap.add_argument("--target-mean", type=float, default=0.20,
                    help="land-mean albedo for --mode scaled")
    ap.add_argument("--vegetation-albedo", type=float, default=None,
                    help="albedo of vegetated ground for --mode vegetated. "
                         "Defaults to model.vegetation_albedo in the config, "
                         "derived for THIS star by analysis/vegetation_albedo.py; "
                         "0.15 is the Earth-Sun endmember it replaced (SPEC-5), "
                         "pass it explicitly to reproduce the old surface")
    ap.add_argument("--forest-fraction", type=float, default=None,
                    help="override the forest fraction implied by --mode")
    ap.add_argument("--flat-bands", action="store_true",
                    help="write one identical field to 174, 175 and 176. What "
                         "this script did before the rock table had a band "
                         "split, and the only way to reproduce a surface built "
                         "then; the radiation reads 175 and 176 and not 174 at "
                         "NSIMPLEALBEDO=0, so this asserts every material on "
                         "the planet reflects equally either side of 0.75 um")
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model = config["model"]
    mode = args.mode or model.get("land_albedo_source", "lithology")
    veg_from_config = False
    if mode == "vegetated" and args.vegetation_albedo is None:
        # No silent 0.15: that default was an Earth-Sun value nothing had
        # re-weighted for this star, and reviving it by omission would be the
        # defect SPEC-5 closed. The config value is the derived one.
        try:
            args.vegetation_albedo = float(model["vegetation_albedo"])
            veg_from_config = True    # recorded in the report; see below
        except KeyError:
            raise SystemExit(
                "--mode vegetated needs model.vegetation_albedo in the config "
                "or an explicit --vegetation-albedo; the derivation is "
                "analysis/vegetation_albedo.py, and 0.15 is the un-reweighted "
                "Earth-Sun endmember if that is really what is meant")
    # Same rule as --vegetation-albedo above and for the same reason: no silent
    # Earth-Sun default. 0.13 and 0.19 were the two endmembers left behind when
    # 0.15 was re-weighted and moved to the config, and `modelled` is the mode
    # that blends them. world-9m5.
    cover_from_config = []
    if mode == "modelled":
        for name, earth_sun in (("tree_albedo", 0.13), ("grass_albedo", 0.19)):
            if getattr(args, name) is not None:
                continue
            try:
                setattr(args, name, float(model[name]))
                cover_from_config.append(name)
            except KeyError:
                flag = "--" + name.replace("_", "-")
                raise SystemExit(
                    f"--mode modelled needs model.{name} in the config or an "
                    f"explicit {flag}; the derivation is "
                    f"analysis/vegetation_albedo.py, and {earth_sun} is the "
                    "un-reweighted Earth-Sun endmember if that is really what "
                    "is meant")
    nlat, nlon = int(model["latitudes"]), int(model["longitudes"])

    # Deliberately from the grid, not from config: see builds.resolution_of.
    grid_dir = args.grid or grid_export(config)
    # From the grid, not from config: the two differ exactly when someone
    # builds for another resolution, which is when the filename matters.
    resolution = resolution_of(grid_dir)
    output = args.output or (INPUTS / resolution.lower())

    if mode == "uniform":
        for code in ALBEDO_CODES + (FOREST_CODE,):
            p = output / f"orogen_{resolution}_surf_{code:04d}.sra"
            if p.exists():
                p.unlink()
        print("mode=uniform: removed any albedo SRA; ExoPlaSim will use albland=0.22")
        return

    mesh = Export(args.mesh or mesh_export(config))
    rock = mesh.substrate_class
    # Classes that cannot carry a canopy. This was a single hardcoded reference
    # to evaporite, which silently became wrong when Orogen split that class:
    # playa_clastic is 18.8% of land and nothing roots in playa mud either, but
    # it would have been handed vegetation albedo because it is "not evaporite".
    # Driven from config so the next class addition is a config change.
    barren_codes = model.get("barren_rock_classes") or ["evaporite"]
    barren = np.zeros(rock.shape, dtype=bool)
    barren_applied = []
    for code in barren_codes:
        try:
            barren |= rock == _rock_id(mesh.root, code)
            barren_applied.append(code)
        except KeyError:
            pass    # class absent from this export; older builds lack playa_clastic
    is_land = mesh.surface_class == LAND

    # Per-region substrate albedo, then integrated over land only. Taking this
    # from the gridded export instead would average open water's 0.06 into every
    # coastal cell.
    region_albedo = mesh.rock_albedo.astype(np.float64).copy()

    # The band ratio travels beside the albedo and is repainted wherever the
    # albedo is, because a repaint changes the MATERIAL and the ratio is a
    # property of the material. The two arrays are only ever assigned together;
    # the recombination check at the end is what catches it if they are not.
    z1, z2 = (float(v) for v in stellar.band_fractions())
    rock_ratio = rock_band_ratios(ROCK_BANDS, mesh.root)
    region_rho = np.empty(region_albedo.shape, dtype=np.float64)
    for rid, value in rock_ratio.items():
        region_rho[rock == rid] = value
    veg_json = json.loads(VEGETATION_BANDS.read_text(encoding="utf-8"))
    canopy_rho = float(veg_json["band2_over_band1"])
    cover_rho = {k: float(v) for k, v in veg_json["cover_band_ratio"].items()}

    # Rock-class overrides, applied before anything else reads the field. These
    # exist because the export's albedo table was written for plausibility rather
    # than radiative accuracy, and one of its entries dominates this planet's
    # energy balance. Recorded in the report so no result is quoted without it.
    overrides = model.get("lithology_albedo_overrides") or {}
    applied = {}
    for code, spec in overrides.items():
        value = float(spec["albedo"] if isinstance(spec, dict) else spec)
        expected = float(spec["replaces"]) if isinstance(spec, dict) and "replaces" in spec else None
        rid = _rock_id(mesh.root, code)
        sel = rock == rid
        if not sel.any():
            raise RuntimeError(f"override names rock class {code!r}, absent from this export")
        exported = float(np.median(region_albedo[sel]))
        if expected is not None and abs(exported - expected) > 1e-6:
            raise RuntimeError(
                f"override for {code!r} was written against an exported albedo of "
                f"{expected}, but this export reports {exported}. The class has "
                "been redefined; re-derive the override or delete it rather than "
                "applying it to a class that no longer means the same thing."
            )
        applied[code] = {"rock_id": rid, "albedo": value,
                         "exported_albedo": exported, "regions": int(sel.sum())}
        region_albedo[sel] = value

    # The two endmembers the error budget's biosphere item is the difference of,
    # captured here because this is the only place both are unambiguous.
    #
    # They are taken on the SAME rock table -- overrides applied -- and both
    # PRE-LAKE. Neither is optional. `land_mean_bare_rock` in the report below
    # is a DIFFERENT quantity: the raw export mean that `scaled` mode divides
    # by, still carrying the exported playa albedo the override exists to
    # replace, so pairing it with a vegetated mean compares two rock tables and
    # understates the bare end by the override times playa's share of land. And
    # pairing a lake-composited vegetated mean with an uncomposited bare one
    # double-counts the lakes, which are their own line in the same budget.
    area_land = mesh.cell_area.astype(np.float64)[is_land]

    def _land_mean(field: np.ndarray) -> float:
        return float(np.average(field[is_land], weights=area_land))

    endmembers = {
        "note": "land-mean albedo of the two surface endmembers, on the "
                "overridden rock table and before lakes. Their difference is "
                "what the biosphere is worth; pair only these two.",
        "bare_rock": _land_mean(region_albedo),
    }
    veg_painted = None
    if mode == "vegetated":
        veg_painted = is_land & ~barren
        region_albedo[veg_painted] = args.vegetation_albedo
        region_rho[veg_painted] = canopy_rho
        endmembers["vegetated"] = _land_mean(region_albedo)

    # Lakes, last, because a lake covers whatever lithology is under it and no
    # vegetation grows on open water.
    #
    # Applied per region and then integrated, which is what makes this a
    # composite rather than a mask: a lake occupying part of a T42 cell leaves
    # that cell with an area-weighted albedo, and almost every lake here is far
    # below the grid. No mask is flipped, per notes/lake-representation.md --
    # these are surface-property perturbations on cells that stay land.
    #
    # This closes a gap rather than refining one. The lake solution existed and
    # was rendered on maps, but never reached the climate: the cells were being
    # given their dry substrate albedo, which on this planet is bright playa and
    # salt crust exactly where the water is.
    lake_report = None
    lake_mask = None
    # Bound here rather than inside the lake block, which is the only place it is
    # filled in: the BOOTSTRAP run's albedo is built with no lakes, because
    # they need a climatology this terrain does not have yet, and the report at
    # the end reads this whether or not that block ran.
    evap_report = None
    if args.lakes is not None:
        area_r = mesh.cell_area.astype(np.float64)
        water_albedo_value = float(next(
            r["albedo"] for r in mesh.manifest["lithology"]["rockClasses"]
            if r["code"] == "water"))
        with Dataset(args.lakes) as lds:
            lake = np.asarray(lds["lake"][:]).astype(bool)
            lake_terrain = getattr(lds, "terrain_hash", None)
        if lake.shape != region_albedo.shape:
            raise RuntimeError(
                f"{args.lakes} has {lake.shape[0]} regions, this mesh has "
                f"{region_albedo.shape[0]}; they are different builds")
        if lake_terrain and lake_terrain != mesh.terrain_hash:
            raise RuntimeError(
                f"lake solution was computed on terrain {lake_terrain[:16]} and "
                f"this mesh is {mesh.terrain_hash[:16]}. Lake extent is a "
                "property of the basin floor, so a solution from another build "
                "does not describe these lakes; re-run surface_water.py.")
        lake &= is_land

        # --- derived evaporite: salt crust is the EPHEMERAL zone -------------
        #
        # Orogen assigns `evaporite` geometrically -- the lowest quarter of each
        # basin's depth range -- because it has no water balance to consult. That
        # is the only rule it can apply, and it overstates salt: 2.85% of land
        # against Earth's modern salt flats at 0.13-0.34%.
        #
        # Salt crust is where brine repeatedly evaporates to saturation, which is
        # the shoreline zone that wets and dries, not the deep floor. A lake
        # shallower than the local annual evaporation depth dries within an orbit
        # and leaves crust; a deeper one persists and is painted as water below.
        # Deriving it here gives 0.47% of land, and two independent routes agree
        # the geometric rule is 6-7x too generous.
        #
        # This belongs in the derived-surface classifier once that exists; it is
        # here because albedo is what reaches the model and the classifier is not
        # built yet. See pedology/notes/derived-surface-classes.md.
        if args.climatology is not None and Path(args.climatology).is_file():
            import run_exoplasim as _rx
            year_s = float(_rx.derive(
                config, float(config["orbit"]["baseline_flux_earth"])
            )["orbital_year_seconds"])
            with Dataset(args.climatology) as cds:
                # Weighted by records per bin, which are unequal. CLIM-13.
                ev_grid = ((-climatology.annual_mean(
                    np.asarray(cds["evap"][:]),
                    np.asarray(cds["time"][:])) * year_s).ravel()
                    if "evap" in cds.variables else None)
            if ev_grid is not None:
                cellidx, _nla, _nlo = region_cells(mesh, grid_dir)
                evaporite_id = next(r["id"] for r in
                                    mesh.manifest["lithology"]["rockClasses"]
                                    if r["code"] == "evaporite")
                with Dataset(args.lakes) as lds:
                    depth_m = np.asarray(lds["lake_depth_km"][:]) * 1000.0
                e_local = ev_grid[cellidx]
                ephemeral = is_land & (depth_m > 0) & (depth_m <= e_local)
                salt_a = float(next(r["albedo"] for r in
                                    mesh.manifest["lithology"]["rockClasses"]
                                    if r["code"] == "evaporite"))
                playa_a = float(next(r["albedo"] for r in
                                     mesh.manifest["lithology"]["rockClasses"]
                                     if r["code"] == "playa_clastic"))
                geo_salt = is_land & (mesh.substrate_class == evaporite_id)
                before_e = float(np.average(region_albedo[is_land],
                                            weights=area_r[is_land]))
                # Geometric salt that is not ephemeral becomes playa; ephemeral
                # ground becomes salt crust regardless of what Orogen called it.
                region_albedo[geo_salt & ~ephemeral] = playa_a
                region_rho[geo_salt & ~ephemeral] = rock_ratio[
                    _rock_id(mesh.root, "playa_clastic")]
                region_albedo[ephemeral] = salt_a
                region_rho[ephemeral] = rock_ratio[evaporite_id]
                after_e = float(np.average(region_albedo[is_land],
                                           weights=area_r[is_land]))
                evap_report = {
                    "method": "salt crust = regions whose solved lake depth is at "
                              "or below the local annual evaporation depth, so "
                              "they dry within one orbit",
                    # The albedo now depends on a CLIMATOLOGY, which it did not
                    # before. Record which one: an evaporite split derived from
                    # another climate describes another world's salt flats.
                    "climatology": str(args.climatology),
                    "orbital_year_seconds": year_s,
                    "geometric_salt_fraction_of_land":
                        float(area_r[geo_salt].sum() / area_r[is_land].sum()),
                    "derived_salt_fraction_of_land":
                        float(area_r[ephemeral].sum() / area_r[is_land].sum()),
                    "salt_albedo": salt_a, "playa_albedo": playa_a,
                    "land_mean_albedo_before": round(before_e, 6),
                    "land_mean_albedo_after": round(after_e, 6),
                    "delta": round(after_e - before_e, 6),
                }

        # Ephemeral regions dry within an orbit and keep their salt crust;
        # only persistent lakes are painted as water. lake_mask stays the full
        # lake set because crust is as unvegetable as water downstream.
        lake_mask = lake
        try:
            paint_water = lake & ~ephemeral
        except NameError:
            paint_water = lake
        before = float(np.average(region_albedo[is_land], weights=area_r[is_land]))
        region_albedo[paint_water] = water_albedo_value
        region_rho[paint_water] = rock_ratio[_rock_id(mesh.root, "water")]
        after = float(np.average(region_albedo[is_land], weights=area_r[is_land]))
        lake_report = {
            "source": str(args.lakes),
            "terrain_hash": lake_terrain,
            "lake_regions": int(lake.sum()),
            "lake_fraction_of_land": float(area_r[lake].sum() / area_r[is_land].sum()),
            "water_albedo": water_albedo_value,
            "land_mean_albedo_before": round(before, 6),
            "land_mean_albedo_after": round(after, 6),
            "delta": round(after - before, 6),
        }

    fraction, alb_grid, _empty = land_weighted(mesh, grid_dir, region_albedo)
    land_cells = fraction >= float(model["geography_land_threshold"])

    # Codes 175 and 176, from the same regions and the same gridding as 174.
    # One identical field in all three was a fair statement about a rock table
    # with no spectral split and is false about every material on this planet:
    # a canopy is dark below 0.75 um and bright above it, a soil more so, and
    # halite runs the other way. The model reads the pair and not 174 at
    # NSIMPLEALBEDO = 0.
    #
    # `land_weighted` is an area-weighted mean, so it is linear in its argument
    # and the recombination identity survives gridding; the assertion before
    # the write is what proves that rather than this comment.
    band_grids = None
    if not args.flat_bands:
        shape1, shape2 = band_shapes(region_rho, z1, z2)
        band_grids = [land_weighted(mesh, grid_dir, region_albedo * shape1)[1],
                      land_weighted(mesh, grid_dir, region_albedo * shape2)[1]]

    raw_fraction, raw_alb, _ = land_weighted(mesh, grid_dir,
                                             mesh.rock_albedo.astype(np.float64))
    area = mesh.cell_area.astype(np.float64)
    raw_mean = float(np.average(mesh.rock_albedo[is_land], weights=area[is_land]))

    if mode == "scaled":
        scaled_grid = np.clip(alb_grid * (args.target_mean / raw_mean), 0.05, 0.80)
        # The bands take the SAME per-cell factor, clip included, so the
        # identity survives a mode whose whole point is to move the level.
        if band_grids is not None:
            with np.errstate(invalid="ignore", divide="ignore"):
                factor = np.where(alb_grid > 0, scaled_grid / alb_grid, 1.0)
            band_grids = [g * factor for g in band_grids]
        alb_grid = scaled_grid

    # --- the loop closure -----------------------------------------------------
    # Blend the lithology substrate against what actually grew, cell by cell,
    # instead of asserting a single land-surface albedo. This is the only mode
    # whose answer came from a model rather than from a choice.
    vegetation_summary = None
    tree_cover = None
    if mode == "modelled":
        if args.vegetation is None:
            raise SystemExit("--mode modelled needs --vegetation <run>/fpc.out")
        cover = read_foliar_cover(args.vegetation)
        nlat_g, nlon_g = alb_grid.shape

        # Coordinates come from the climatology, not from the export's
        # planet.nc, and this is not a preference.
        #
        # The two label the same grid differently. planet.nc runs its longitudes
        # from -178.5938, ExoPlaSim's output from 0, offset by half a cell and an
        # origin. The *grids* are identical: taking the land mask from each and
        # comparing gives 4106 cells both ways and 100% cellwise agreement. Only
        # the labels disagree.
        #
        # LPJ-GUESS's printed coordinates descend from the driver file, which was
        # built from the climatology, so matching against planet.nc's labels
        # returns zero cells out of 4106. That failure has now happened three
        # times in this project on three different scripts. Share the coordinate
        # source; never reconstruct one.
        if args.climatology is None or not args.climatology.is_file():
            raise SystemExit(
                "--mode modelled needs --climatology, the same one the LPJ-GUESS "
                "driver was built from, for its coordinates")
        with Dataset(args.climatology) as data:
            lat_axis = np.asarray(data["lat"][:], dtype=float)
            lon_axis = np.asarray(data["lon"][:], dtype=float)
        if (len(lat_axis), len(lon_axis)) != (nlat_g, nlon_g):
            raise SystemExit(
                f"climatology grid is {len(lat_axis)}x{len(lon_axis)} but the "
                f"surface grid is {nlat_g}x{nlon_g}")
        lon_signed = np.where(lon_axis > 180.0, lon_axis - 360.0, lon_axis)

        tree_cover = np.zeros_like(alb_grid)
        grass_cover = np.zeros_like(alb_grid)
        matched = 0
        for j in range(nlat_g):
            for i in range(nlon_g):
                if not land_cells[j, i]:
                    continue
                entry = cover.get((round(float(lon_signed[i]), 2),
                                   round(float(lat_axis[j]), 2)))
                if entry is None:
                    continue
                tree_cover[j, i], grass_cover[j, i] = entry
                matched += 1

        # LPJ-GUESS is not told which ground is salt crust or playa mud. It is
        # given a soil texture there and will happily grow on it if the climate
        # allows, so the barren classes are masked out here rather than trusted.
        # They are 12.65% of carved-zoned land and the brightest surfaces on the
        # planet, so letting vegetation cover them would be a real error.
        barren_fraction = land_fraction_of_class(mesh, grid_dir, barren)
        rootable = np.clip(1.0 - barren_fraction, 0.0, 1.0)
        tree_cover *= rootable
        grass_cover *= rootable

        total_cover = np.clip(tree_cover + grass_cover, 0.0, 1.0)

        # The same blend on the bands, with each endmember carrying its own
        # band pair. The substrate term keeps the pair the lithology gave it,
        # so a modelled cell is a real mixture of a canopy spectrum and a rock
        # spectrum rather than one broadband number copied three times. This is
        # the mode the docstring calls the only one that is not an assumption,
        # and until now it was the one mode with no split at all.
        #
        # The weights are whatever the cover says; they are not renormalised
        # here, so the identity holds through the clip above without this block
        # needing to agree with the broadband line about what the weights mean.
        cover_shapes = {}
        for name in ("tree", "grass"):
            cs1, cs2 = band_shapes(cover_rho[name], z1, z2)
            cover_shapes[name] = (cs1, cs2)
        if band_grids is not None:
            blended = []
            for index in (0, 1):
                blended.append(np.where(
                    land_cells,
                    tree_cover * args.tree_albedo * cover_shapes["tree"][index]
                    + grass_cover * args.grass_albedo * cover_shapes["grass"][index]
                    + (1.0 - total_cover) * band_grids[index],
                    band_grids[index]))
            band_grids = blended

        alb_grid = np.where(
            land_cells,
            tree_cover * args.tree_albedo
            + grass_cover * args.grass_albedo
            + (1.0 - total_cover) * alb_grid,
            alb_grid)

        missing = int(land_cells.sum()) - matched
        vegetation_summary = {
            "source": str(args.vegetation),
            "land_cells_matched": matched,
            "land_cells_without_vegetation_left_bare": missing,
            "tree_albedo": args.tree_albedo,
            "grass_albedo": args.grass_albedo,
            "cover_albedo_source": (
                "config model." + ", config model.".join(cover_from_config)
                if cover_from_config else "explicit on the command line"),
            "cover_albedo_derivation": "analysis/vegetation_albedo.json",
            "tree_albedo_bands": [round(args.tree_albedo * cover_shapes["tree"][0], 5),
                                  round(args.tree_albedo * cover_shapes["tree"][1], 5)],
            "grass_albedo_bands": [round(args.grass_albedo * cover_shapes["grass"][0], 5),
                                   round(args.grass_albedo * cover_shapes["grass"][1], 5)],
            "cover_band_ratio": {k: round(v, 4) for k, v in cover_rho.items()},
            "barren_masked": ("barren classes forced to zero cover; LPJ-GUESS is "
                              "not told which ground is salt crust or playa"),
            "coordinate_source": str(args.climatology),
        }

    # Ocean cells carry the water value; ExoPlaSim computes ocean albedo itself,
    # so it is inert, but the array has to be full.
    water_albedo = float(next(r["albedo"] for r in json.loads(
        (mesh.root / "manifest.json").read_text(encoding="utf-8"))["lithology"]["rockClasses"]
        if r["code"] == "water"))
    field = np.where(land_cells, alb_grid, water_albedo)

    if mode == "modelled" and args.forest_fraction is None:
        # Tree cover per cell, which is what 212 is for. No longer a constant.
        forest_value = None
        forest = np.where(land_cells, np.clip(tree_cover, 0.0, 1.0), 0.0)
    else:
        forest_value = (args.forest_fraction if args.forest_fraction is not None
                        else MODE_FOREST_FRACTION[mode])
        # Open water carries no canopy, so lake regions come out of the
        # vegetable fraction as well as out of the albedo. Leaving them in
        # would put forest on the lakes and, through dforest, walk the albedo
        # we just set back toward the forested endpoint.
        canopy = ~barren if lake_mask is None else (~barren & ~lake_mask)
        vegetable = land_fraction_of_class(mesh, grid_dir, canopy)
        forest = np.where(land_cells, vegetable * forest_value, 0.0)

    gw = np.fromfile(grid_dir / "grid" / "gauss_weights.bin", dtype="float64")
    def gmean(a):
        num = (np.where(land_cells, a, 0.0).mean(axis=1) * gw).sum()
        den = (land_cells.mean(axis=1) * gw).sum()
        return float(num / den) if den else 0.0
    final_mean = gmean(field)

    output.mkdir(parents=True, exist_ok=True)
    written = []
    band_fields = {code: field for code in ALBEDO_CODES}
    band_means = None
    band_check = {"split": False, "reason": "--flat-bands"}
    if band_grids is not None:
        # Ocean gets the water class's own pair rather than its broadband value
        # in both bands. It is inert -- the model computes ocean albedo itself
        # -- but a flat pair there would make the identity below pass over
        # ocean for the wrong reason, and a check that cannot fail on two
        # thirds of the grid is not much of a check.
        ws1, ws2 = band_shapes(rock_ratio[_rock_id(mesh.root, "water")], z1, z2)
        band_fields[175] = np.where(land_cells, band_grids[0], water_albedo * ws1)
        band_fields[176] = np.where(land_cells, band_grids[1], water_albedo * ws2)
        band_means = {"band1_land_mean": gmean(band_fields[175]),
                      "band2_land_mean": gmean(band_fields[176])}

        # The check that can fail. Every level in this script is reached by a
        # chain of repaints, and the band arrays follow that chain through a
        # separate variable; if any repaint moved one without the other, or if
        # a mode rescaled the broadband and not the pair, the recombination
        # stops returning the broadband field. Asserted on the arrays about to
        # be written, over the whole grid, and reported as a number rather than
        # a boolean so a later change can see it move.
        recombined = z1 * band_fields[175] + z2 * band_fields[176]
        residual = np.abs(recombined - band_fields[174])
        worst = float(residual.max())
        if worst > RECOMBINATION_TOLERANCE:
            j, i = np.unravel_index(int(residual.argmax()), residual.shape)
            raise SystemExit(
                f"the band pair does not recombine to the broadband field: "
                f"worst cell ({j}, {i}) is off by {worst:.3e} against a "
                f"tolerance of {RECOMBINATION_TOLERANCE:.1e}. "
                f"174 = {band_fields[174][j, i]:.6f}, "
                f"175 = {band_fields[175][j, i]:.6f}, "
                f"176 = {band_fields[176][j, i]:.6f}, "
                f"z = ({z1:.5f}, {z2:.5f}). A repaint has moved a level "
                "without moving its band ratio, or a mode has rescaled one and "
                "not the other.")
        band_check = {
            "split": True,
            "band_flux_fractions": [round(z1, 6), round(z2, 6)],
            "band_flux_source": "lib/stellar.band_fractions",
            "rock_band_ratios": "analysis/rock_albedo_bands.json",
            "canopy_band_ratios": "analysis/vegetation_albedo.json",
            "recombination_worst_residual": worst,
            "recombination_tolerance": RECOMBINATION_TOLERANCE,
            "note": "z1*175 + z2*176 == 174 over the whole grid, checked before "
                    "writing. The .sra format rounds at %12.5f, so the "
                    "tolerance is twice that quantum and not tighter.",
        }
    for code in ALBEDO_CODES:
        path = output / f"orogen_{resolution}_surf_{code:04d}.sra"
        write_sra(path, code, band_fields[code])
        written.append(str(path))
    forest_path = output / f"orogen_{resolution}_surf_{FOREST_CODE:04d}.sra"
    write_sra(forest_path, FOREST_CODE, forest)
    written.append(str(forest_path))

    report = {
        "mode": mode,
        "mesh": str(mesh.root),
        "grid": str(grid_dir),
        "resolution": resolution,
        "terrain_hash": mesh.terrain_hash,
        "codes": list(ALBEDO_CODES) + [FOREST_CODE],
        "forest_fraction_value": forest_value,
        "vegetation": vegetation_summary,
        "forest_fraction_land_mean": gmean(forest),
        "lithology_albedo_overrides": applied,
        "vegetation_albedo": (None if mode != "vegetated" else {
            "value": args.vegetation_albedo,
            "source": ("config model.vegetation_albedo" if veg_from_config
                       else "explicit --vegetation-albedo"),
            "derivation": "analysis/vegetation_albedo.json",
            "band_ratio": canopy_rho,
            "bands_175_176": (
                None if args.flat_bands
                else [round(args.vegetation_albedo * v, 5)
                      for v in band_shapes(canopy_rho, z1, z2)]),
        }),
        "bands": {**band_check, **(band_means or {})},
        "barren_rock_classes": barren_applied,
        "forest_note": ("dforest blends snow albedo between forested and "
                        "unforested endpoints, so it has to agree with the "
                        "background albedo assumption. ExoPlaSim's default is a "
                        "uniform 0.5."),
        "land_mean_bare_rock": raw_mean,
        "land_mean_written": final_mean,
        "endmembers": endmembers,
        "lakes": lake_report,
        "derived_evaporite": evap_report,
        "exoplasim_default_albland": 0.22,
        "land_min": float(field[land_cells].min()),
        "land_max": float(field[land_cells].max()),
        "files": written,
        "caveat": ("Substrate albedo, not land-surface albedo. Vegetation and "
                   "snow are applied by the model on top of this."),
    }
    # Provenance stamp; lib/provenance.py owns the shape and the inert set.
    # `inputs` is the half the config stamp cannot cover: every DERIVED FILE
    # this generator read whose content no configuration key names. A re-run of
    # either band-shape derivation moves 175 and 176 without moving a single
    # config key, and `check_consistency.py` compares these hashes against the
    # files as they now stand, on the same footing as the config stamp. The
    # optional per-run inputs are listed too, because a lake solution or a
    # climatology re-derived under the same name is the same failure.
    inputs = [ROCK_BANDS, VEGETATION_BANDS]
    inputs += [p for p in (args.vegetation, args.lakes, args.climatology)
               if p is not None]
    report.update(config_stamp(config, "exoplasim/scripts/build_surface_albedo.py",
                               inputs=inputs))
    (output / "albedo_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
