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

Wet ground and dry ground
-------------------------

Two ends, not one. Codes 174, 175 and 176 carry the DRY endmember, the
reflectance of each material at or near zero water content, which is where every
spectrum this project integrates was measured. Codes 1742, 1750 and 1760 carry
the SATURATED endmember of the same material in the same three surfaces, from
`analysis/soil_albedo_wetting.json`. `landmod`'s `getalb` under `nwetsoil = 1`
mixes between them on the degree of saturation of the land column's surface
layer, so the modelled soil albedo now moves through a wetting and drying cycle
instead of holding one value through it.

**The state comes from one place and this script does not define it.** It is the
top water layer of the layered land column, 0.02 m thick, declared in
`pedology/config/land_column_properties.yaml` under `surface_layer`: that block
is the ONE top-of-column profile, and DUST-17's emitting depth reads the same
profile at its own depth. A near-surface saturation defined here would be the
second central hydrology that row exists to prevent. Three things about that
layer matter to this field. It empties under evaporation in about a day where
the 0.5 m layer above it relaxed on the column's timescale, which is what makes
an albedo driven from it right in PHASE -- Craft and Horel measure a salt crust
recovering from 0.22 to 0.32 in eight days. Its capacity is cut from AIR DRY
rather than from the wilting point, so an empty layer and the dry field written
here are the same state and the mixing has no level shift at its dry end. And
the store it reads is LIQUID water, so a frozen surface reads as a dry one.

**The shape between the ends is not linear in the albedo.** Sadeghi, Jones and
Philpot (2015) derive the reflectance of a wetting soil from Kubelka-Munk
two-flux theory, and what interpolates linearly in the water content is the
transformed reflectance `r = (1-R)^2 / (2R)`, not `R`. The CLM form this row
started from, `alpha_dry*(1-S) + alpha_wet*S`, has the same two ends and the
wrong curve between them, and the difference is largest at low saturation, which
is where a drying skin spends its time.

**The saturated end closes band 2 without a band-2 measurement.** Lekner and Dorf
(1988) and Twomey, Bohren and Mergenthaler (1986) both give the wet reflectance
as a closed form in the DRY reflectance at the same wavelength, so applying
either per wavelength to the spectra `analysis/rock_albedo_bands.py` already
reads produces a band-2 wet endmember from band-2 dry data. The two bands come
out wetting by different amounts because the dry spectra differ between them.
The two papers describe two different surfaces -- a water film over a rough
solid, and interstitial water in a finely divided medium -- so they are the two
arms of a bracket rather than two estimates of one number; `--wetting-arm`
sweeps it, and the interstitial arm is staged because this world's land is
regolith and playa fill far more than it is bare outcrop.

**Two classes are REFUSED and staged wet equal to dry.** `water`, because open
water is not a surface that wets. And `evaporite`, because the sign is disputed
on it: Malek et al. (1990) and Craft and Horel (2019) both measure a wetted
halite crust much darker in situ, and a twenty-year MODIS series over the largest
halite pan reads wet or cool years brighter than dry or warm ones. One saturation
number cannot carry both, and the one physical route by which wetting brightens a
surface -- the Fresnel reflection of a continuous film, which Sadeghi's Eq (19)
adds -- is at most the normal-incidence reflectance of water times the porosity,
under 0.01, against an interannual difference near 0.075. So the brightening is a
crust-CONDITION effect, dissolution and reprecipitation and dust, which this
model carries no state for. The refusal is in the staged FIELD and not only in
the report: a refused class is written with its dry pair as its wet pair, and the
check before the write requires a cell made entirely of refused material to
carry identical dry and saturated fields.

`moisture_dependence` in the report is the magnitude, on two footings. The
ceiling is the substrate's own albedo minus open water's, because a wetted
modelled surface cannot be darker than the open water that would cover it if the
wetting went all the way. `band1_wetting` is tighter and covers less: the
band-1 loss at full wetting, from the wet/dry ratio bracket below, with band 2's
unbracketed flux share reported beside it. Both are land-mean albedo deltas,
which is what `scripts/error_budget.py` consumes. The ceiling's sign claim is
checked before anything is written: a land region at or below open water's
albedo would be one the ceiling does not bound, and the script refuses rather
than reporting a bound it has not got.
"""

from __future__ import annotations

import argparse
import json
import sys
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
from lpj_output import reduce_table, require_lpj_acceptance
from rootable import read_rootable_partition

# 174 broadband, 175 below 0.75 um, 176 above. With NSIMPLEALBEDO=0 the
# radiation uses the two-band pair; 174 is written too so the broadband
# diagnostic agrees rather than silently keeping 0.22.
ALBEDO_CODES = (174, 175, 176)

# The SATURATED pair, the far end of the moisture mixing, on the companion-code
# precedent 1730 and 2290 already set: 1742, 1750 and 1760 sit beside 174, 175
# and 176 and carry the same three surfaces at full saturation. `surfmod.f90`
# registers them and `landmod`'s `getalb` mixes toward them under
# `nwetsoil = 1`. ONE record each rather than fourteen: the wet reflectance of
# a material is a property of the material, and the season is carried by the
# saturation of the surface layer that the mixing reads.
WET_ALBEDO_CODES = (1742, 1750, 1760)

# Where the wetting ratios come from. Per class and per band, derived per
# WAVELENGTH from the same reflectance spectra the dry band ratios come from.
WETTING = PROJECT_ROOT / "analysis" / "soil_albedo_wetting.json"

# WHICH ARM IS STAGED, AND WHY IT IS NOT AN AVERAGE. `soil_albedo_wetting.json`
# carries two, from two mechanisms that describe two different surfaces: Lekner
# and Dorf's water FILM over a rough solid, and Twomey, Bohren and
# Mergenthaler's INTERSTITIAL water in a finely divided medium. This world's
# land is regolith, playa mud and fan fill far more than it is bare outcrop --
# `playa_clastic` alone is the largest single class -- and TBM's own paper
# names finely divided media as its case while Lekner and Dorf name rough solid
# surfaces as theirs. So the interstitial arm is staged and the film arm is the
# other end of the declared bracket, swept with `--wetting-arm`.
WETTING_ARMS = ("twomey", "lekner")

# THE CLASSES THE MOISTURE TERM REFUSES, and the refusal travels with the
# material the way every other per-class property here does.
#
# `evaporite` is refused because the SIGN is disputed on it and one saturation
# number cannot carry both signs. Malek et al. (1990) and Craft and Horel
# (2019) both measure a wetted halite crust much darker in situ; a twenty-year
# MODIS series over the largest halite pan on Earth reads wet or cool years
# BRIGHTER than dry or warm ones. They are not measuring the same thing -- a
# water film darkens a crust, and a wet year also dissolves and reprecipitates
# it and keeps dust off it -- and Sadeghi's Fresnel term, the one physical
# route by which wetting could brighten a surface, is an order of magnitude too
# small to carry the second: it adds at most the normal-incidence reflectance
# of water times the porosity, under 0.01, against an interannual difference
# near 0.075. So the brightening branch is a crust-CONDITION effect that this
# model has no state for, and staging a single-signed darkening on this class
# would assert against a held measurement.
#
# `water` is refused because open water is not a surface that wets.
#
# A refused class stages its DRY pair as its wet pair, so the mixing returns
# the dry albedo at every saturation and the refusal is visible in the staged
# field rather than only in a report.
WETTING_REFUSED = {
    "evaporite": "the sign is disputed: two in-situ pyranometer pairs make a "
                 "wetted halite crust much darker and a twenty-year MODIS "
                 "series makes wet years brighter, and the Fresnel term is too "
                 "small by an order of magnitude to be the second effect",
    "water": "open water is not a surface that wets",
}

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

# The only WETTING measurement this project holds, and it covers one of the two
# bands. Penndorf (1956) Table 1, the Sewing Handbook column, read by eye off
# the render because that file's text layer is machine-generated: clay soil
# 7.5/15, sand 18/31, bare rich soil 5.5/7.2. Ratios are taken WITHIN one column
# so the measurement's own preparation and geometry cancel, which is the same
# discipline analysis/playa_albedo.py applies to ECOSTRESS.
#
# It is LUMINOUS reflectance, eye-weighted over 0.38-0.77 um, so it bears on
# band 1 and says nothing whatever about band 2 -- and band 2 is where liquid
# water's absorption sits and where the larger share of this star's flux
# arrives. Penndorf's own text puts +/-25 percent on "the actual variations
# under natural conditions". The LEVELS do not transfer to this world and are
# not used; the ratio is what transfers, which is the reading
# references/INDEX.md already records for this source.
PENNDORF_WET_DRY_RATIO = (0.500, 0.764)

# The BROADBAND wetting of a salt crust, from the two in-situ pyranometer
# sources this project has read. Malek et al. (1990) on Pilot Valley playa,
# 0.24 wet against above 0.75 after three dry weeks; Craft and Horel (2019) on
# the Bonneville Salt Flats, 0.22 under 20-40 mm of flooding against 0.45 dry
# summer crust. Both cover BOTH bands, which Penndorf does not, and both are the
# darkening direction. A twenty-year MODIS series over Uyuni runs the other way
# on interannual composites -- see the note -- so this bracket describes the
# water film and not the state of the crust between floods.
SALT_CRUST_WET_DRY_RATIO = (0.320, 0.489)

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


def read_foliar_cover(path: Path, peers: list[Path] | None = None):
    """Tree and grass foliar projective cover per cell, from a run's fpc.out.

    BIO-12's shared reducer rejects a trending or incomplete end window and
    averages ten complete forcing cycles. Coordinates are rounded to two
    decimals by that reader because that is LPJ-GUESS's output precision.
    """
    require_lpj_acceptance(path)
    for peer in peers or ():
        require_lpj_acceptance(peer)
    reduced = reduce_table(path, peers or ())
    grass = [i for i, name in enumerate(reduced.names) if name in GRASS_PFTS]
    trees = [i for i, name in enumerate(reduced.names)
             if name not in GRASS_PFTS and name != "Total"]
    cover = {
        key: (sum(values[i] for i in trees), sum(values[i] for i in grass))
        for key, values in reduced.values.items()
    }
    return cover, reduced.report


def composite_rootable(total_surface: np.ndarray,
                       rootable_substrate_contribution: np.ndarray,
                       rootable_fraction: np.ndarray,
                       tree_fpc: np.ndarray, grass_fpc: np.ndarray,
                       tree_albedo: float, grass_albedo: float):
    """Replace canopy only inside rootable area, preserving all other surface.

    ``total_surface`` and ``rootable_substrate_contribution`` are contributions
    to the whole model-land cell, not conditional means. Thus subtracting the
    latter leaves solved water and dry barren ground exactly where the native
    mesh placed them. LPJ FPC is conditional on the rootable environment, so it
    is multiplied by BIO-11's rootable fraction exactly once.
    """
    arrays = [np.asarray(value, dtype=float) for value in (
        total_surface, rootable_substrate_contribution, rootable_fraction,
        tree_fpc, grass_fpc)]
    if any(value.shape != arrays[0].shape for value in arrays[1:]):
        raise ValueError("rootable compositing arrays have different shapes")
    total_surface, root_contribution, rootable, tree, grass = arrays
    if np.any((rootable < -1e-9) | (rootable > 1.0 + 1e-9)):
        raise ValueError("rootable fraction is outside [0,1]")
    if np.any(tree < -1e-9) or np.any(grass < -1e-9):
        raise ValueError("LPJ foliar cover is negative")
    cover = np.clip(tree + grass, 0.0, 1.0)
    tree_share = rootable * tree
    grass_share = rootable * grass
    # When overlapping FPCs sum above one, retain their relative shares while
    # enforcing the same total-cover ceiling the old code used.
    raw_share = tree_share + grass_share
    rescale = np.ones_like(raw_share)
    np.divide(rootable * cover, raw_share, out=rescale, where=raw_share > 0)
    tree_share *= rescale
    grass_share *= rescale
    nonrootable = total_surface - root_contribution
    result = (nonrootable + (1.0 - cover) * root_contribution
              + tree_share * tree_albedo + grass_share * grass_albedo)
    return result, tree_share, grass_share, cover


def _rock_id(mesh_dir: Path, code: str) -> int:
    lit = json.loads((mesh_dir / "manifest.json").read_text(encoding="utf-8"))["lithology"]
    for r in lit["rockClasses"]:
        if r["code"] == code:
            return int(r["id"])
    raise KeyError(f"no rock class {code!r} in {mesh_dir}")


LAND_COLUMN_CONTRACT = (PROJECT_ROOT / "pedology" / "config"
                        / "land_column_properties.yaml")


def _configured_column(config: dict) -> dict:
    """The land water column the run is configured with, beside the albedo depth.

    Three facts and one comparison, all READ. The mixing staged here needs the
    model's first water layer to BE the depth `land_column_properties.yaml`
    declares the albedo consumer reads at, because `landmod` has no access to
    that contract and cannot check it: `run_exoplasim.py` refuses a three-layer
    column cut anywhere else, and this is what makes the state visible in the
    report a reader has in front of them rather than only in a refusal they
    have to provoke.

    `agrees` is the comparison and not a gate. The gate is at the run; a
    surface build is legitimately made before the column is cut, and saying so
    is more use than refusing.
    """
    column = (config.get("surface", {}) or {}).get("land_water_column", {}) or {}
    moisture = (config.get("surface", {}) or {}).get("soil_albedo_moisture", {}) or {}
    thicknesses = column.get("layer_thickness_m")
    contract = yaml.safe_load(LAND_COLUMN_CONTRACT.read_text(encoding="utf-8"))
    albedo_depth = float(contract["surface_layer"]["consumers"]["albedo"]["depth_m"])
    first = float(thicknesses[0]) if thicknesses else None
    return {
        "scheme": column.get("scheme"),
        # The thickness list IS the layer count, so the count is not emitted
        # beside it. Whether the config's own declared count agrees with the
        # list is `run_exoplasim.py`'s refusal and not this report's claim.
        "layer_thickness_m": thicknesses,
        "albedo_depth_m": albedo_depth,
        "albedo_depth_source": "pedology/config/land_column_properties.yaml "
                               "surface_layer.consumers.albedo.depth_m",
        "first_layer_is_the_albedo_depth": (
            first is not None and abs(first - albedo_depth) <= 1.0e-9),
        "soil_albedo_moisture_enabled": bool(moisture.get("enabled", False)),
    }


def _full_layer_fall(config: dict, dry_fields: dict, wet_fields: dict,
                     gmean) -> dict:
    """What the mixing is worth at a PERMANENTLY FULL surface layer.

    The staged saturated pair is the endmember at a degree of saturation of
    one, and the modelled surface layer caps at field capacity, so the model
    never reaches it. The bound the term can actually be judged against is the
    mixing evaluated at `saturation_at_full_layer`, which is where a skin
    sitting at its capacity everywhere and always would put it. That is the
    number an arm's land-mean `alb` difference is compared with and the number
    `scripts/error_budget.py` prices, and neither should have to re-derive it
    from the two endmembers in a document.

    PER CELL AND THEN MEANED, never the other way round. The mixing is concave
    in the albedo pair, so mixing two land means and meaning the mixing of the
    pairs are different numbers, and only the second is what the model does.

    THE MIXING IS `analysis/soil_albedo_wetting.py`'s and is not restated here.
    That module holds Sadeghi Eq (13), reproduces the papers it comes from, and
    is checked bitwise against the compiled `wet_soil_albedo`; a second copy
    would be free to agree with neither.
    """
    sys.path.insert(0, str(PROJECT_ROOT / "analysis"))
    from soil_albedo_wetting import sadeghi_mix

    moisture = (config.get("surface", {}) or {}).get("soil_albedo_moisture", {}) or {}
    saturation = float(moisture["saturation_at_full_layer"])
    sigma = {
        1742: float(moisture.get("shape_sigma_broadband", 1.0)),
        1750: float(moisture.get("shape_sigma_band1", 1.0)),
        1760: float(moisture.get("shape_sigma_band2", 1.0)),
    }
    out = {
        "saturation": saturation,
        "saturation_source": "config/planet.yaml "
                             "surface.soil_albedo_moisture.saturation_at_full_layer",
        "note": "the land-mean albedo at a surface layer permanently at its "
                "own capacity, mixed per cell through Sadeghi Eq (13) and then "
                "meaned. A BOUND on the realised term and not an estimate of "
                "it: the realised fall is this times how wet the modelled skin "
                "is over the orbit, and it acts only where the surface is "
                "snow-free, ice-free land, so the realised area is below the "
                "land fraction and one-signed downward.",
    }
    for code, dry_code, name in ((1742, 174, "broadband"), (1750, 175, "band1"),
                                 (1760, 176, "band2")):
        mixed = gmean(sadeghi_mix(dry_fields[dry_code], wet_fields[code],
                                  saturation, sigma[code]))
        out[f"{name}_land_mean"] = mixed
        out[f"{name}_fall"] = gmean(dry_fields[dry_code]) - mixed
    return out


def band_shapes(rho: np.ndarray, z1: float, z2: float):
    """Dimensionless (band1, band2) multipliers for a band ratio `rho`.

    `z1*s1 + z2*s2 == 1` by construction, so multiplying any broadband albedo
    by the pair gives a pair that recombines to it in the model's own weights.
    """
    s1 = 1.0 / (z1 + z2 * rho)
    return s1, rho * s1


def rock_wetting_ratios(path: Path, mesh_root: Path, arm: str
                        ) -> tuple[dict[int, tuple[float, float]], dict]:
    """Rock class id -> (band1, band2) wet/dry ratio, for one arm.

    Keyed by CODE in the file and resolved to this export's ids here, on the
    same argument `rock_band_ratios` gives: ids move when Orogen adds a class
    and codes do not. A class the file does not cover raises rather than
    defaulting, because a silent 1.0 is a class asserted never to darken when
    wet, which is the omission this table exists to remove -- and the classes
    that ARE asserted not to darken are named in `WETTING_REFUSED` with their
    reason.
    """
    if not path.is_file():
        raise SystemExit(
            f"{path} is absent. Run analysis/soil_albedo_wetting.py, which "
            "derives the saturated endmember per class and per band; the wet "
            "pair cannot be staged without it.")
    table = json.loads(path.read_text(encoding="utf-8"))
    index = WETTING_ARMS.index(arm)
    lit = json.loads((mesh_root / "manifest.json").read_text(
        encoding="utf-8"))["lithology"]
    ratios: dict[int, tuple[float, float]] = {}
    refused: dict[str, str] = {}
    for entry in lit["rockClasses"]:
        code = entry["code"]
        if code in WETTING_REFUSED:
            ratios[int(entry["id"])] = (1.0, 1.0)
            refused[code] = WETTING_REFUSED[code]
            continue
        if code not in table["classes"]:
            raise SystemExit(
                f"{path.name} has no wetting ratio for rock class {code!r}, "
                f"which {mesh_root} carries. Add a proxy in "
                "analysis/soil_albedo_wetting.py and re-run it, or name the "
                "class in WETTING_REFUSED with its reason; a class with no "
                "ratio would be staged as never darkening when wet.")
        row = table["classes"][code]
        ratios[int(entry["id"])] = (float(row["band1_wet_over_dry"][index]),
                                    float(row["band2_wet_over_dry"][index]))
    return ratios, {
        "arm": arm,
        "arms_available": list(WETTING_ARMS),
        "source": str(path),
        "generated": table.get("generated"),
        "refused_classes": refused,
        "held_pair_bracket_passes": table["held_pair_bracket"]["passes"],
        "band2_liquid_absorption_direction":
            table["band2_liquid_absorption"]["direction"],
    }


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
    ap.add_argument("--vegetation-peer", type=Path, action="append", default=[],
                    help="comparable fpc.out with another root seed or patch "
                         "count; repeat to measure stochastic uncertainty")
    ap.add_argument("--rootable", type=Path, default=None,
                    help="BIO-11 rootable-fraction artifact; required by and "
                         "defaulted for --mode modelled")
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
    ap.add_argument("--wetting-arm", choices=WETTING_ARMS, default="twomey",
                    help="which arm of the wet-endmember bracket to stage. "
                         "`twomey` is the interstitial mechanism and is the "
                         "one this world's regolith and playa surfaces are; "
                         "`lekner` is the water-film mechanism and the other "
                         "end of the bracket. Not an average: they are two "
                         "surfaces, not two estimates of one")
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
    # The WETTING ratio travels with the material for the same reason the band
    # ratio does, and is repainted wherever that is: a repaint changes what the
    # ground is made of, and how much darker it gets when wet is a property of
    # the material and not of the cell. Two numbers per region, one per band,
    # because the two bands wet by different amounts -- which is the whole
    # result `soil_albedo_wetting.py` exists to produce.
    rock_wet, wetting_report = rock_wetting_ratios(
        WETTING, mesh.root, args.wetting_arm)
    region_wet1 = np.empty(region_albedo.shape, dtype=np.float64)
    region_wet2 = np.empty(region_albedo.shape, dtype=np.float64)
    for rid, (w1, w2) in rock_wet.items():
        region_wet1[rock == rid] = w1
        region_wet2[rock == rid] = w2
    # WHAT THE GROUND IS MADE OF AFTER THE REPAINTS, which is not what Orogen
    # called it. The derived evaporite split turns geometric salt that is not
    # ephemeral into playa and ephemeral ground into salt crust whatever it
    # started as, and a solved lake covers whatever lies under it. Every one of
    # those moves the wetting ratio, so the class the refusal is asserted about
    # has to move with them: `rock` is the geometry and this is the material.
    region_material = rock.copy()
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

    # The substrate as the moisture term would act on it: overrides applied,
    # nothing painted over it yet. Snapshot rather than recompute, because the
    # array below is mutated in place by the vegetation paint, the derived
    # evaporite split and the lakes, and the moisture term acts on the ground
    # under all three. See "Wet ground and dry ground" above. The band ratio is
    # snapshot beside it for the same reason the two are only ever assigned
    # together: a level without its material is not a spectrum.
    substrate_albedo = region_albedo.copy()
    substrate_rho = region_rho.copy()
    veg_painted = None
    if mode == "vegetated":
        veg_painted = is_land & ~barren
        region_albedo[veg_painted] = args.vegetation_albedo
        region_rho[veg_painted] = canopy_rho
        # A canopy is not the wetting surface. The moisture term acts on soil,
        # and what a wetted canopy does is a different question with different
        # sources; painting the substrate's ratio onto it would answer that
        # question by accident.
        region_wet1[veg_painted] = 1.0
        region_wet2[veg_painted] = 1.0
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
                region_material[geo_salt & ~ephemeral] = _rock_id(
                    mesh.root, "playa_clastic")
                pw1, pw2 = rock_wet[_rock_id(mesh.root, "playa_clastic")]
                region_wet1[geo_salt & ~ephemeral] = pw1
                region_wet2[geo_salt & ~ephemeral] = pw2
                region_albedo[ephemeral] = salt_a
                region_rho[ephemeral] = rock_ratio[evaporite_id]
                region_material[ephemeral] = evaporite_id
                # Salt crust, which the moisture term refuses: the sign is
                # disputed on this class and the refusal follows the material
                # through the derived split, not the geometric one.
                ew1, ew2 = rock_wet[evaporite_id]
                region_wet1[ephemeral] = ew1
                region_wet2[ephemeral] = ew2
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
        region_material[paint_water] = _rock_id(mesh.root, "water")
        ww1, ww2 = rock_wet[_rock_id(mesh.root, "water")]
        region_wet1[paint_water] = ww1
        region_wet2[paint_water] = ww2
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
    wet_grids = None
    if not args.flat_bands:
        shape1, shape2 = band_shapes(region_rho, z1, z2)
        band_grids = [land_weighted(mesh, grid_dir, region_albedo * shape1)[1],
                      land_weighted(mesh, grid_dir, region_albedo * shape2)[1]]
        # The saturated pair, through the SAME gridding and from the same
        # regions, so the two ends of the mixing describe one surface.
        wet_grids = [
            land_weighted(mesh, grid_dir,
                          region_albedo * shape1 * region_wet1)[1],
            land_weighted(mesh, grid_dir,
                          region_albedo * shape2 * region_wet2)[1]]

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
            wet_grids = [g * factor for g in wet_grids]
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
        cover, equilibrium_window = read_foliar_cover(
            args.vegetation, args.vegetation_peer)
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

        tree_fpc = np.zeros_like(alb_grid)
        grass_fpc = np.zeros_like(alb_grid)
        matched = 0
        for j in range(nlat_g):
            for i in range(nlon_g):
                if not land_cells[j, i]:
                    continue
                entry = cover.get((round(float(lon_signed[i]), 2),
                                   round(float(lat_axis[j]), 2)))
                if entry is None:
                    continue
                tree_fpc[j, i], grass_fpc[j, i] = entry
                matched += 1

        # BIO-11 is the population contract. It subtracts persistent solved
        # water and dry barren substrate once, including water over barren, and
        # LPJ's FPC is conditional on what remains. Re-derive the native masks
        # only to prove the surface inputs supplied to THIS builder are the
        # same partition; the artifact is the fraction used in the algebra.
        partition, rootable_provenance = read_rootable_partition(
            config, lat_axis, lon_axis, args.rootable, land=land_cells)
        rootable = partition["rootable"]
        # BIO-11 follows the solved-lake extent. In the albedo arm an
        # evaporation diagnostic may stage part of that extent as ephemeral
        # salt crust rather than open water, but both populations remain
        # non-rootable. Their optical levels already differ in region_albedo;
        # only their common exclusion from canopy belongs here.
        native_water = (np.zeros_like(is_land) if lake_mask is None
                        else lake_mask & is_land)
        native_barren = barren & is_land & ~native_water
        native_rootable = is_land & ~native_water & ~barren
        native_partition = {}
        for name, population in (("rootable", native_rootable),
                                 ("water", native_water),
                                 ("barren", native_barren)):
            native_partition[name] = land_weighted(
                mesh, grid_dir, population.astype(np.float64))[1]
        partition_residual = max(
            float(np.max(np.abs(native_partition[name][land_cells]
                                - partition[name][land_cells])))
            for name in partition)
        if partition_residual > 2.0e-6:
            raise SystemExit(
                "the lake/barren inputs supplied to surface albedo disagree "
                "with BIO-11's rootable partition by "
                f"{partition_residual:.3e}; rebuild both from the same derived "
                "surface rather than compositing two worlds")

        # Contribution of ROOTABLE substrate to the whole land cell. This is
        # deliberately not a conditional mean. Subtracting it from alb_grid
        # leaves persistent water and dry barren contributions untouched; the
        # old `(1-cover) * alb_grid` scaled those contributions down and then
        # painted canopy over the missing share.
        rootable_broadband = land_weighted(
            mesh, grid_dir,
            np.where(native_rootable, region_albedo, 0.0))[1]
        pre_modelled_broadband = alb_grid.copy()

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
            blended, wet_blended = [], []
            for index in (0, 1):
                shape = shape1 if index == 0 else shape2
                wet_ratio = region_wet1 if index == 0 else region_wet2
                rootable_band = land_weighted(
                    mesh, grid_dir, np.where(
                        native_rootable, region_albedo * shape, 0.0))[1]
                rootable_wet_band = land_weighted(
                    mesh, grid_dir, np.where(
                        native_rootable,
                        region_albedo * shape * wet_ratio, 0.0))[1]
                dry, _tree, _grass, _cover = composite_rootable(
                    band_grids[index], rootable_band, rootable,
                    tree_fpc, grass_fpc,
                    args.tree_albedo * cover_shapes["tree"][index],
                    args.grass_albedo * cover_shapes["grass"][index])
                wet, _tree, _grass, _cover = composite_rootable(
                    wet_grids[index], rootable_wet_band, rootable,
                    tree_fpc, grass_fpc,
                    args.tree_albedo * cover_shapes["tree"][index],
                    args.grass_albedo * cover_shapes["grass"][index])
                blended.append(np.where(land_cells, dry, band_grids[index]))
                wet_blended.append(np.where(land_cells, wet, wet_grids[index]))
            band_grids = blended
            wet_grids = wet_blended

        composite, tree_cover, grass_cover, conditional_cover = composite_rootable(
            alb_grid, rootable_broadband, rootable,
            tree_fpc, grass_fpc, args.tree_albedo, args.grass_albedo)
        alb_grid = np.where(land_cells, composite, alb_grid)
        total_cover = tree_cover + grass_cover
        pure_water = land_cells & (partition["water"] >= 1.0 - 2.0e-6)
        pure_barren = land_cells & (partition["barren"] >= 1.0 - 2.0e-6)
        nonrootable_error = 0.0
        pure_nonrootable = pure_water | pure_barren
        if pure_nonrootable.any():
            nonrootable_error = float(np.max(np.abs(
                alb_grid[pure_nonrootable]
                - pre_modelled_broadband[pure_nonrootable])))
            if nonrootable_error > RECOMBINATION_TOLERANCE:
                raise SystemExit(
                    "modelled cover moved a wholly non-rootable cell by "
                    f"{nonrootable_error:.3e}; water and barren ground cannot "
                    "receive canopy")
            if np.any(tree_cover[pure_nonrootable] != 0.0):
                raise SystemExit(
                    "modelled forest code 212 would be nonzero on wholly water "
                    "or barren ground")

        missing = int(land_cells.sum()) - matched
        vegetation_summary = {
            "source": str(args.vegetation),
            "equilibrium_window": equilibrium_window,
            "rootable_surface": rootable_provenance,
            "rootable_partition_max_absolute_residual": partition_residual,
            "pure_water_cells_unchanged": int(pure_water.sum()),
            "pure_barren_cells_unchanged": int(pure_barren.sum()),
            "pure_nonrootable_max_absolute_change": nonrootable_error,
            "cover_semantics": (
                "LPJ FPC is conditional on rootable ground; code 212 and the "
                "canopy albedo shares multiply it by f_rootable once. Solved "
                "water and dry barren contributions are retained exactly."),
            "conditional_rootable_cover_land_mean": float(
                np.mean(conditional_cover[land_cells])),
            "effective_cover_land_mean": float(np.mean(total_cover[land_cells])),
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
            "barren_and_lakes_masked": (
                "BIO-11 rootable fraction excludes dry barren substrate and "
                "persistent solved water, with overlap deducted once"),
            "coordinate_source": str(args.climatology),
        }

    # Ocean cells carry the water value; ExoPlaSim computes ocean albedo itself,
    # so it is inert, but the array has to be full.
    water_albedo = float(next(r["albedo"] for r in json.loads(
        (mesh.root / "manifest.json").read_text(encoding="utf-8"))["lithology"]["rockClasses"]
        if r["code"] == "water"))
    # --- the moisture term, declared absent -----------------------------------
    #
    # The ceiling on what a moisture-dependent soil albedo could be worth, and
    # the check that the ceiling means what it says. A wetted modelled surface
    # cannot be darker than the open water that would cover it if the wetting
    # went all the way, so the substrate's own albedo minus open water's bounds
    # the darkening branch per region. Taken on the substrate rather than on the
    # written field, because vegetation and snow override the soil albedo in the
    # model and the moisture term acts on what is underneath them.
    #
    # THE CHECK THAT CAN FAIL: a land region at or below open water's albedo is
    # one this bound does not bound, and on such a region the claim that wetting
    # darkens the ground is not something the export supports. Orogen's darkest
    # land class is basalt and open water is well below it, so this passes with
    # room today; it fails if a class is added or re-valued into that range,
    # which is exactly when the declaration would need rewriting.
    substrate_land = substrate_albedo[is_land]
    ceiling_land = substrate_land - water_albedo
    if float(ceiling_land.min()) <= 0.0:
        dark = int((ceiling_land <= 0.0).sum())
        raise SystemExit(
            f"{dark} land regions carry a substrate albedo at or below open "
            f"water's {water_albedo}, the darkest being "
            f"{float(substrate_land.min()):.4f}. The moisture-term ceiling in "
            "exoplasim/notes/soil-albedo-moisture.md bounds a DARKENING, and it "
            "does not bound ground already darker than the water that would "
            "cover it. Re-derive the bound for those classes before writing a "
            "report that claims it.")
    moisture_dependence = {
        "state": "staged",
        "what": "codes 174, 175 and 176 carry the DRY endmember and 1742, 1750 "
                "and 1760 the SATURATED one, and `landmod`'s `getalb` under "
                "`nwetsoil = 1` mixes between them on the degree of saturation "
                "of the land column's surface layer",
        "form": "Sadeghi, Jones and Philpot (2015) Eq (13): linear in the "
                "Kubelka-Munk transformed reflectance r = (1-R)^2/(2R) and "
                "therefore NOT in the albedo. The CLM form this row started "
                "from, alpha_dry*(1-S) + alpha_wet*S, has the same two ends and "
                "the wrong shape between them",
        "state_source": "pedology/config/land_column_properties.yaml "
                        "`surface_layer`, which is the ONE top-of-column "
                        "profile DUST-17's emitting depth reads at its own "
                        "depth. Not a second hydrology",
        "what_is_still_open": [
            "THE SIGN ON SALT CRUST, which is refused per class rather than "
            "resolved: see `wetting.refused_classes`. Two in-situ pyranometer "
            "pairs make a wetted halite crust much darker and a twenty-year "
            "MODIS series makes wet years brighter; the second is a "
            "crust-CONDITION effect this model has no state for, and the "
            "Fresnel term is too small by an order of magnitude to be it",
            "THE LIQUID'S OWN ABSORPTION IN BAND 2. Both mechanisms behind the "
            "saturated endmember treat the water as non-absorbing, so both are "
            "UPPER bounds there and the bracket is open at its dark end. "
            "Closing it needs an absolute scattering coefficient for the dry "
            "soil or a wet-and-dry spectrum pair on one sample",
            "THE SHAPE PARAMETER IN BAND 1. Sadeghi's sigma is one where the "
            "water's own scattering is negligible, which they verify in the "
            "short-wave infrared; their visible-band fits run 0.042 to 0.528, "
            "and a value below one darkens the modelled surface sooner in a "
            "wetting cycle without moving its ends",
        ],
        # WHAT THE RUN ACTUALLY CONFIGURES, read rather than asserted. This
        # used to be a sentence in the list above claiming the column was
        # still cut at 0.5 and 1.0 m, and it stayed there after the cut moved:
        # a report that states what a config SAYS goes stale the moment the
        # config does, and nothing re-runs to catch it. Emitting the column
        # instead puts the comparison in front of the reader, who can see the
        # surface cut against the albedo depth the contract declares.
        "configured_land_water_column": _configured_column(config),
        "argument": "exoplasim/notes/soil-albedo-moisture.md",
        "ceiling_note":
            "the most the darkening branch could be worth, as a land-mean albedo "
            "delta, which is the unit scripts/error_budget.py consumes. A "
            "CEILING and not an estimate: the realised term is this times a "
            "wet-area fraction and a wetting efficiency, neither of which is "
            "known here. It bounds a FULLY INUNDATED surface, which is a "
            "different quantity from the saturated endmember staged in 1742, "
            "1750 and 1760: the modelled surface layer caps at field capacity "
            "and never reaches it",
        "open_water_albedo": water_albedo,
        "ceiling_land_mean": round(
            float(np.average(ceiling_land, weights=area_land)), 6),
        "ceiling_min_region": round(float(ceiling_land.min()), 6),
        "ceiling_max_region": round(float(ceiling_land.max()), 6),
    }

    # The one band the magnitude can be bracketed on rather than merely bounded.
    # Penndorf's wet/dry ratio is a band-1 quantity, so this is the broadband
    # albedo delta a fully wetted modelled land surface would show from band 1
    # ALONE, at the two ends of that ratio. Band 2 carries no bracket at all and
    # carries the larger flux share, which is reported beside it so the number
    # is not mistaken for the whole term.
    sub_s1, _sub_s2 = band_shapes(substrate_rho, z1, z2)
    substrate_band1_mean = float(np.average(
        (substrate_albedo * sub_s1)[is_land], weights=area_land))
    moisture_dependence["band1_wetting"] = {
        "source": "Penndorf (1956) Table 1, Sewing Handbook column; ratios "
                  "taken within that column",
        "wet_over_dry_ratio_bracket": list(PENNDORF_WET_DRY_RATIO),
        "substrate_band1_land_mean": round(substrate_band1_mean, 6),
        "broadband_delta_if_all_land_wet": [
            round(z1 * substrate_band1_mean * (1.0 - r), 6)
            for r in sorted(PENNDORF_WET_DRY_RATIO, reverse=True)],
        "band2_flux_fraction_with_no_bracket": round(z2, 6),
        "note": "band 1 only, and a FULLY wetted simulated land surface, which "
                "no state here can say the area of. Luminous reflectance is "
                "eye-weighted over 0.38-0.77 um and band 1 is flat to 0.75 um, "
                "so the ratio is transferred across two weightings; Penndorf "
                "puts +/-25 percent on the natural variation. Band 2 has no "
                "measurement, and it is the band liquid water absorbs in. "
                "Taken on the substrate's own spectrum, so under --flat-bands "
                "it describes the material and not the flat pair staged.",
    }

    # Both bands, one class: salt crust is where this project holds in-situ
    # pyranometer pairs, so the darkening branch can be bracketed outright
    # there instead of merely bounded. The level shift is the number that
    # decides the verdict -- what arming the mixing on the current land column
    # would apply to every salt-crust cell in every season -- and it is stated
    # per unit of the reachable saturation floor rather than multiplied by it
    # here, because that floor belongs to
    # pedology/config/land_column_properties.yaml and is not restated in this
    # tree. See the note for the arithmetic.
    try:
        salt_sel = is_land & (rock == _rock_id(mesh.root, "evaporite"))
    except KeyError:
        salt_sel = np.zeros(rock.shape, dtype=bool)
    salt_share = float(area[salt_sel].sum() / area[is_land].sum())
    salt_level = float(np.median(substrate_albedo[salt_sel])) if salt_sel.any() else 0.0
    moisture_dependence["salt_crust_wetting"] = {
        "source": "Malek et al. (1990) and Craft, Horel (2019), both in-situ "
                  "and broadband; see references/INDEX.md",
        "wet_over_dry_ratio_bracket": list(SALT_CRUST_WET_DRY_RATIO),
        "class": "evaporite",
        "share_of_land": round(salt_share, 6),
        "dry_level": round(salt_level, 6),
        "per_cell_delta_if_wet": [
            round(salt_level * (1.0 - r), 6)
            for r in sorted(SALT_CRUST_WET_DRY_RATIO, reverse=True)],
        "land_mean_contribution_if_wet": [
            round(salt_share * salt_level * (1.0 - r), 6)
            for r in sorted(SALT_CRUST_WET_DRY_RATIO, reverse=True)],
        "note": "covers BOTH bands, unlike band1_wetting, and one class. The "
                "sign is disputed for this class on interannual timescales; "
                "this bracket is the water film, not the state of the crust "
                "between floods.",
    }

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

        # --- the saturated pair, codes 1742, 1750 and 1760 -----------------
        #
        # The broadband member is RECOMBINED from the pair rather than gridded
        # on its own, so the same identity holds at the wet end by construction
        # and there is no third array to keep in step.
        wet_fields = {
            1750: np.where(land_cells, wet_grids[0], water_albedo * ws1),
            1760: np.where(land_cells, wet_grids[1], water_albedo * ws2),
        }
        wet_fields[1742] = z1 * wet_fields[1750] + z2 * wet_fields[1760]

        # THREE CHECKS THAT CAN FAIL, on the arrays about to be written.
        #
        # The wet field must not be brighter than the dry one anywhere: the
        # term staged here is a DARKENING and a cell where it is not is a
        # repaint that moved one array and not the other. The tolerance is the
        # `.sra` write quantum, because the two fields reach the model through
        # that format and a difference below it is not a difference.
        #
        # It must not be negative, which the maps cannot produce but a
        # composite of them could if a cover weight went out of range.
        #
        # And the refused classes must be staged EQUAL to their dry pair, cell
        # by cell, because that is what a refusal is: the mixing returns the
        # dry albedo at every saturation. A refusal that is only a sentence in
        # a report is not a refusal.
        for code, dry_code in ((1750, 175), (1760, 176), (1742, 174)):
            excess = float((wet_fields[code] - band_fields[dry_code]).max())
            if excess > RECOMBINATION_TOLERANCE:
                j, i = np.unravel_index(
                    int((wet_fields[code] - band_fields[dry_code]).argmax()),
                    wet_fields[code].shape)
                raise SystemExit(
                    f"the saturated field {code} is brighter than the dry "
                    f"field {dry_code} at cell ({j}, {i}) by {excess:.3e}. "
                    "Wetting darkens; a cell where it does not is a repaint "
                    "that moved the level without moving the wetting ratio.")
            if float(wet_fields[code].min()) < 0.0:
                raise SystemExit(
                    f"the saturated field {code} is negative somewhere. The "
                    "maps cannot produce that, so a cover weight or a repaint "
                    "is out of range.")
        refused_ids = [_rock_id(mesh.root, code) for code in WETTING_REFUSED
                       if any(r["code"] == code for r in json.loads(
                           (mesh.root / "manifest.json").read_text(
                               encoding="utf-8"))["lithology"]["rockClasses"])]
        refused_regions = np.isin(region_material, refused_ids) & is_land
        # The refusal at the MATERIAL level, which is where it is declared and
        # the only place it is always testable: a T21 cell need not be made
        # entirely of one class, so the cell-level check below can have no
        # sample, and this one cannot.
        #
        # `region_material` AND NOT `rock`, which is the geometry. Reading the
        # geometric class made this refuse the first field ever built with
        # lakes: the derived evaporite split had turned non-ephemeral
        # geometric salt into PLAYA, which is a class that wets and is
        # correctly given a ratio, while the check still asked whether
        # evaporite carried one. The refusal is about what the ground is made
        # of after every repaint, which is the only thing the staged field
        # describes.
        if refused_regions.any():
            off = float(max(np.abs(region_wet1[refused_regions] - 1.0).max(),
                            np.abs(region_wet2[refused_regions] - 1.0).max()))
            if off != 0.0:
                raise SystemExit(
                    f"a refused rock class carries a wetting ratio {off:.3e} "
                    "away from 1. A refusal is staged as wet equal to dry, and "
                    "a repaint has given refused material a ratio.")
        refused_pure = None
        if refused_regions.any():
            # Cells made ENTIRELY of refused material, which are the only ones
            # where the staged fields must match exactly: a mixed cell carries
            # some material that does wet.
            share = land_weighted(mesh, grid_dir,
                                  refused_regions.astype(np.float64))[1]
            pure = land_cells & (share >= 1.0 - 1.0e-9)
            if pure.any():
                gap = float(np.abs(wet_fields[1750][pure]
                                   - band_fields[175][pure]).max())
                if gap > RECOMBINATION_TOLERANCE:
                    raise SystemExit(
                        f"a cell made entirely of refused material differs "
                        f"between its dry and saturated band-1 fields by "
                        f"{gap:.3e}. A refused class is staged wet EQUAL to "
                        "dry, so the mixing returns the dry albedo at every "
                        "saturation; this one would move.")
            refused_pure = int(pure.sum())
        wetting_report["refused_regions"] = int(refused_regions.sum())
        wetting_report["cells_entirely_refused_material"] = refused_pure
        wetting_report["band1_land_mean_saturated"] = gmean(wet_fields[1750])
        wetting_report["band2_land_mean_saturated"] = gmean(wet_fields[1760])
        wetting_report["broadband_land_mean_saturated"] = gmean(wet_fields[1742])
        wetting_report["broadband_land_mean_dry"] = final_mean
        wetting_report["codes"] = list(WET_ALBEDO_CODES)
        wetting_report["at_full_surface_layer"] = _full_layer_fall(
            config, band_fields, wet_fields, gmean)
        for code in WET_ALBEDO_CODES:
            path = output / f"orogen_{resolution}_surf_{code:04d}.sra"
            write_sra(path, code, wet_fields[code])
            written.append(str(path))
    else:
        wetting_report["staged"] = False
        wetting_report["why_not"] = ("--flat-bands writes one field to all "
                                     "three dry codes and there is no pair for "
                                     "the wet end to be a pair with")

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
        "codes": list(ALBEDO_CODES) + [FOREST_CODE]
                 + (list(WET_ALBEDO_CODES) if band_grids is not None else []),
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
        "moisture_dependence": moisture_dependence,
        "wetting": wetting_report,
        "lakes": lake_report,
        "derived_evaporite": evap_report,
        "exoplasim_default_albland": 0.22,
        "land_min": float(field[land_cells].min()),
        "land_max": float(field[land_cells].max()),
        "files": written,
        "caveat": ("Substrate albedo, not land-surface albedo. Vegetation and "
                   "snow are applied by the model on top of this. It is the DRY "
                   "substrate; `moisture_dependence` says what that omits and "
                   "what it is worth."),
    }
    # Provenance stamp; lib/provenance.py owns the shape and the inert set.
    # `inputs` is the half the config stamp cannot cover: every DERIVED FILE
    # this generator read whose content no configuration key names. A re-run of
    # either band-shape derivation moves 175 and 176 without moving a single
    # config key, and `check_consistency.py` compares these hashes against the
    # files as they now stand, on the same footing as the config stamp. The
    # optional per-run inputs are listed too, because a lake solution or a
    # climatology re-derived under the same name is the same failure.
    inputs = [ROCK_BANDS, VEGETATION_BANDS, WETTING]
    inputs += [p for p in (args.vegetation, args.lakes, args.climatology)
               if p is not None]
    report.update(config_stamp(config, "exoplasim/scripts/build_surface_albedo.py",
                               inputs=inputs))
    (output / "albedo_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
