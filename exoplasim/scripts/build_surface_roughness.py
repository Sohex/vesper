#!/usr/bin/env python3
"""Aerodynamic roughness length per cell, surface code 0173.

`configure()` clears every surface field when handed a landmap, so `dz0clim`
falls back to the uniform namelist default of `dz0land = 2.0 m` over all land.
That is the total roughness, and it is used directly: `NVEG = 0` in these runs,
so `landmod.f90:409` takes `dz0 = dz0clim` and SIMBA's separate vegetation and
orographic terms never run.

A single 2.0 m over all land is wrong in a way that matters here. It asserts
forest-scale roughness over the 16.5% of land that is salt crust and playa, where
real values are nearer 0.001 m, and those surfaces are not scattered: they are the
floors of closed basins, which are flat by construction and are exactly the cells
the carve verdict integrates evaporation over. Roughness enters the turbulent
exchange coefficient as

    ce = k^2 / ln(z_ref / z0)^2

with `k` von Karman's constant and `z_ref` the height of the lowest model level,
which `lib/lapse.py:reference_height_m` derives hypsometrically at THIS planet's
gravity. Over the liquid-water span the default gives about 8.1 to 8.5 times the
exchange a real playa surface would, and therefore that much too much evaporation
from precisely the basins whose water balance decides whether they survive. The
contrast is reported bracketed rather than at one temperature because `z_ref` is
linear in the air temperature and this step runs before any climatology exists to
measure one from.

## What this computes

Total roughness combines a surface term and a subgrid orographic term in
quadrature, which is the same combination `simba.f90:474` uses when it builds
`dz0` from vegetation and topography:

    z0 = sqrt(z0_surface^2 + z0_orographic^2)

`z0_surface` is per mesh region, from land cover: barren classes take a bare-ground
value, everything else takes a canopy value scaled by forest fraction. It is
integrated to the grid over land only, so a coastal cell is not dragged toward
open water.

**It is integrated in `ce`, not in the length.** The model applies one exchange
coefficient to a whole cell and the turbulent flux is linear in it, so what the
cell owes the atmosphere is the area mean of `ce` over its own surfaces, and the
roughness written out is the length that reproduces that mean. `ce` is
logarithmic in `z0`, so the two orders are different reductions, and this
world's land is the case where they part: bare ground and canopy are two orders
apart in `z0`, so a cell that is mostly playa with a canopy minority has its
exchange set by the minority once the lengths are mixed. SPAT-7 measured the
difference on the 10M mesh -- negligible in the land mean, 24 times inside the
`z_ref` bracket this file already declares, but past that bracket on a few
percent of land area which GROWS with refinement, and reaching 30% on closed
basin floors. Those cells are the reason this field exists. The inversion needs
`z_ref`, which is unknown before a climatology; it is taken at the midpoint of
the liquid-water bracket and the report carries what the whole bracket is worth,
which is 1.5e-4 of the land mean.

`z0_orographic` is the part worth retaining the native mesh for. The 10M
fine-support reference contributes about 1,221 regions per global T42 cell on
average (and about 76 even at T170), enough to measure the distribution of
elevation *inside* the cell rather than infer it from resolved slope. Counts
over land vary and are reported by the builder; SPAT-3 is the shared artifact
for comparisons across the full T21/T42/T85/T127/T170 ladder. That subgrid
relief is what a gridded mean elevation cannot supply.

## The one free parameter, what it is anchored on, and where that came from

The constant relating subgrid relief to an effective roughness is not something
this project can derive; formulations in the literature differ by more than an
order of magnitude, and Elvidge et al. (2019) records that the standard deviation
of subgrid orography enters a drag scheme "multiplied by a model-dependent tuning
constant", with the resulting zonal-mean orographic surface stress differing by a
factor of four between operational models of comparable resolution. So the
coefficient is not picked. It is SOLVED, so that the area-weighted land mean of
the result lands on a reference, and the whole question is what the reference is.

**THE REFERENCE IS DERIVED, NOT DECLARED.** The vendored model ships PlaSim's own
boundary dataset for the model's Earth configuration at two rungs, `N032` and
`N064`: code 172 the land mask, code 173 the total roughness `dz0clim`, code 1730
its topography-only part. Reduce code 173 over that dataset's land in `ce` -- the
same operator this file argues for below, at THIS planet's reference height --
and it returns the effective land roughness the model carries when it is given a
real roughness map instead of the uniform fallback. `earth_reference_land_z0()`
computes it here at build time from those files, so nothing about it is a
literal.

That is the derivation the namelist fallback `dz0land` never carried in the model
source, where its only documentation is the comment "roughness length land": the
`N032` answer reproduces the fallback to within a tenth, and `N032` is the
resolution PlaSim's own configuration defaults to. The residual has a direction
rather than being scatter, because the reduction falls as the reference height
rises and this planet's lowest model level sits lower than Earth's at the same
sigma, its gravity being the larger. Both reductions of the dataset are reported,
in `ce` and in the length, along with the fallback, so a reader sees the
agreement and its size rather than taking either on trust.

**WHAT REMAINS UNKNOWN IS THE TRANSFER, AND IT IS BRACKETED.** A reference derived
from the model's Earth dataset is a sound derivation for the wrong planet: this
world's land is not that land, and its relief and lithology are exactly what the
field exists to express. There is no way to close that from inside this project,
so the two ends are declared and both are run.

- `--orographic-arm reference`, the default: the land mean lands on the reference
  above. In a model with no orographic form drag scheme and no gravity wave drag
  scheme, `dz0` is the only channel subgrid orography has, and this arm carries
  as much of it as the model's own Earth dataset does.
- `--orographic-arm none`: the coefficient is zero and the field is the roughness
  its own lithology and land cover imply, with no orographic contribution. This
  is the arm in which subgrid orographic drag is not represented at all.

Both arms' land means are reported by every build whichever one is written, so
the anchored value and the field's own value are always visible together.
`--target-mean` and `--orographic-coefficient` still override, for a third arm or
for a fixed ladder.

What bounds that transfer is a comparison neither side was built to satisfy. The
reference dataset carries its surface and orographic parts separately, so its
land-cover roughness is recoverable; the low arm here IS this world's land-cover
roughness, derived from its own lithology. The two agree to well within a factor
of two at both shipped rungs, and the reference's cover part barely moves between
them while its total falls sharply. So the transfer is importing a level of
subgrid relief and not a land cover, and the report carries the ratio on every
build. It is reported and not gated: no threshold on it was fixed in advance.

**THE REFERENCE MOVES WITH THE RUNG, AND IT HAS TO.** A coarser cell folds more
relief into its own effective roughness, so the reference is a property of the
support as well as of the land: taken from the shipped dataset it falls sharply
between `N032` and `N064` while the surface part of the same dataset barely
moves, which is the signature of an orographic term and not of a land-cover one.
Anchoring every rung on ONE number therefore anchors the coarse rungs and the
fine rungs on different things.

**IT DOES NOT MAKE THE COEFFICIENT RUNG-INVARIANT, AND SPAT-8's CONSTRAINT
STANDS.** Both anchors were run on the same build at both rungs the model ships
a dataset for. Against a fixed target the solved coefficient RISES from T21 to
T42 by about a fifth; against the rung-matched reference it FALLS, by about two
fifths. The reference falls faster with the support than this world's own relief
does, so moving to it reverses the sign of the drift and leaves its size
comparable. A LADDER COMPARISON MUST STILL PASS ONE COEFFICIENT TO EVERY RUNG
with `--orographic-coefficient`, or it measures the calibration as well as the
terrain. `notes/audits/nonlinear-spatial-reductions.md` carries the ladder
numbers.

For a rung the model ships no dataset for the reference cannot be derived, and
the builder refuses rather than reaching for another rung's: pass
`--orographic-coefficient` or `--target-mean` there, which a ladder comparison
has to do anyway.

**This changes climate results.** It is off unless `model.roughness_source`
is set in `config/planet.yaml` (currently `lithology`); changing the key is a
configuration change `continue_exoplasim.py` compares for and refuses to
resume across, so moving it blocks any run in flight.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import yaml

from _paths import CONFIG, INPUTS, PROJECT_ROOT
from sra import read_sra, write_sra
from builds import resolution_of, grid_export, mesh_export
# ONE SOURCE for how much forest the land carries. `model.land_albedo_source`
# already decides it, and code 212 is derived from that decision, so importing
# the mapping is what keeps this field and the albedo field describing one land
# cover. Re-deriving it here would be two formulations of one quantity, which
# is failure class 17. The same idiom continue_exoplasim.py uses on
# run_exoplasim.py, and for the same reason.
from build_surface_albedo import MODE_FOREST_FRACTION
from gridding import (cell_expectation, cell_moments, gaussian_grid,
                      region_cells, transfer_ledger)
from provenance import config_stamp
from orogen import Export, LAND
# ONE derivation of the height a bulk transfer coefficient is taken over.
# hydrography/scripts/carve_verdict.py builds the same z_ref for the same
# reason; this file used to carry the EARTH-gravity answer as a literal.
from lapse import reference_height_m

ROUGHNESS_CODE = 173

KARMAN = 0.4                     # von Karman's constant; `ce` uses its square.

# The reference height is linear in the lowest-level air temperature and this
# step runs before any climatology exists to measure one from, so the reported
# exchange coefficients are BRACKETED over the range in which surface water is
# liquid -- which is the range the field's consumer, lake and playa
# evaporation, is defined over. Fixed here, before any field is built.
CE_BRACKET_K = (273.15, 313.15)

# landmod.f90:51, the uniform FALLBACK this replaces, and it is not the anchor.
# `configure()` falls back to it over all land when no code 173 is staged. It is
# reported beside the derived reference below, which the anchor is taken from.
EXOPLASIM_DZ0LAND_M = 2.0

# PlaSim's own boundary dataset for the model's Earth configuration, shipped with
# the model source, and the ONLY place in this tree where a roughness map the
# model itself uses can be read rather than assumed. Codes: 172 land mask, 173
# total roughness `dz0clim`, 1730 its topography-only part `dz0climo`. Rungs are
# whatever the vendored tree carries; `earth_reference_land_z0` refuses a rung
# that is absent rather than substituting another one, because the reference is
# support-dependent and substituting is what made a ladder comparison measure the
# calibration as well as the terrain.
EARTH_REFERENCE_DIR = PROJECT_ROOT / "vendor/exoplasim/exoplasim/plasim/run"
EARTH_REFERENCE_RUNGS = {"T21": ("N032", 32), "T42": ("N064", 64)}
EARTH_REFERENCE_CODES = {"land_mask": 172, "total": 173, "topography_only": 1730}

# Ocean roughness is set by the model from wind, so the value written on ocean
# cells is inert. It is written as open water's own scale rather than zero so the
# field is never mistaken for a land-only mask.
OCEAN_Z0_M = 1.5e-4

# Surface roughness by cover, metres. Bare ground spans 1e-4 (salt crust) to
# 1e-2 (rocky desert); 3e-3 sits in that range and well below anything vegetated.
DEFAULT_BARE_Z0_M = 0.003
DEFAULT_CANOPY_Z0_M = 0.25       # closed grass to open woodland
DEFAULT_FOREST_Z0_M = 1.0        # closed forest


def neutral_ce(z0, z_ref: float):
    """Neutral bulk exchange coefficient at the lowest model level.

    `fluxmod.f90` forms ONE `(karman / ln(z/z0 + 1))^2` from `dz0` and hands it
    to both `dtransm` and `dtransh`, so this single coefficient carries momentum,
    heat and moisture alike. That is why an orographic contribution folded into
    the roughness cannot be confined to the momentum budget here, and why the
    arm that carries it and the arm that does not are both run.
    """
    return KARMAN ** 2 / np.log(z_ref / np.maximum(z0, 1e-6)) ** 2


def effective_length(ce_value, z_ref: float):
    """The single roughness whose `ce` at `z_ref` is `ce_value`. Inverts above."""
    return z_ref * np.exp(-KARMAN / np.sqrt(np.maximum(ce_value, 1e-30)))


def earth_reference_land_z0(resolution: str, z_ref: float) -> dict:
    """The model's own Earth land roughness, reduced the way this file reduces.

    Reads PlaSim's shipped boundary dataset at the rung `resolution` names and
    returns the effective land roughness of its code 173, taken as the
    area-weighted mean of `ce` over its land and inverted at `z_ref`. This is
    what `dz0land` is a rounding of, and it is COMPUTED here rather than quoted
    so that a change under `vendor/exoplasim` moves the anchor instead of
    silently disagreeing with it.

    The reduction is in `ce` and not in the length for the same reason the cell
    reduction below is: the model applies the coefficient, not the length, and
    the two orders part company over a distribution that spans four decades.
    Both are returned, because their spread is what the choice of reduction is
    worth and it is not negligible.

    Refuses a rung the vendored tree ships no dataset for. The reference depends
    on the support -- a coarser cell folds more relief into its own effective
    roughness -- so answering from another rung's files would import a
    calibration from a grid the build is not on.
    """
    if resolution not in EARTH_REFERENCE_RUNGS:
        raise SystemExit(
            f"no Earth roughness reference for {resolution}: the vendored model "
            f"ships a boundary dataset at "
            f"{', '.join(sorted(EARTH_REFERENCE_RUNGS))} only, and the reference "
            f"is support-dependent, so another rung's files are not a stand-in. "
            f"Pass --orographic-coefficient or --target-mean for this rung, "
            f"which a ladder comparison has to do in any case.")
    prefix, nlat = EARTH_REFERENCE_RUNGS[resolution]
    nlon = 2 * nlat
    fields = {}
    for name, code in EARTH_REFERENCE_CODES.items():
        path = EARTH_REFERENCE_DIR / f"{prefix}_surf_{code:04d}.sra"
        if not path.exists():
            raise SystemExit(
                f"the Earth roughness reference needs {path}, which the vendored "
                f"model source does not carry. It is tracked in the subtree, so "
                f"a missing file means the tree is incomplete rather than that "
                f"the reference has moved.")
        fields[name] = read_sra(path, nlat, nlon)

    weight = gaussian_grid(nlat, nlon).cell_area_fraction()
    # The dataset's land mask is a fraction, not a class, so the population is
    # the cells the model treats as land in it.
    land = fields["land_mask"] > 0.5
    if not land.any():
        raise SystemExit(f"{prefix} land mask selects no cells; the reference "
                         f"cannot be taken over an empty population")
    share = weight[land] / weight[land].sum()

    def land_mean(field):
        return float((field[land] * share).sum())

    total = fields["total"]
    topography = fields["topography_only"]
    # The dataset carries the two parts in quadrature exactly as this file
    # builds them, so the surface part is recoverable and is the one quantity
    # the two planets can be compared on directly: it is land cover, not relief,
    # so it should NOT move with the support and does not.
    surface_part = np.sqrt(np.maximum(total ** 2 - topography ** 2, 0.0))
    ce_bar = land_mean(neutral_ce(total, z_ref))
    return {
        "rung": resolution,
        "dataset": prefix,
        "files": sorted(str((EARTH_REFERENCE_DIR
                             / f"{prefix}_surf_{code:04d}.sra").relative_to(PROJECT_ROOT))
                        for code in EARTH_REFERENCE_CODES.values()),
        "land_cells": int(land.sum()),
        "reference_height_m": round(z_ref, 2),
        "effective_land_z0_m": float(effective_length(ce_bar, z_ref)),
        "land_mean_ce": ce_bar,
        "land_mean_of_length_m": land_mean(total),
        "topography_only_land_mean_m": land_mean(topography),
        "surface_part_land_mean_m": land_mean(surface_part),
        "land_max_m": float(total[land].max()),
        "namelist_fallback_m": EXOPLASIM_DZ0LAND_M,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--mesh", type=Path, default=None)
    ap.add_argument("--grid", type=Path, default=None)
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--orographic-coefficient", type=float, default=None,
                    help="fix the coefficient relating subgrid relief to a "
                         "roughness instead of solving it on this grid. A "
                         "ladder comparison MUST pass one value to every rung: "
                         "the solve is per grid, so two rungs otherwise differ "
                         "by their terrain and by their calibration at once")
    ap.add_argument("--orographic-arm", choices=("reference", "none"),
                    default="reference",
                    help="which end of the declared bracket to WRITE. "
                         "'reference' solves the coefficient so the land mean "
                         "lands on the effective land roughness of the model's "
                         "own Earth boundary dataset at this rung; 'none' sets "
                         "the coefficient to zero, so the field is what this "
                         "world's lithology and land cover imply and subgrid "
                         "orographic drag is absent. Both arms are REPORTED "
                         "whichever is written")
    ap.add_argument("--target-mean", type=float, default=None,
                    help="area-weighted land mean to anchor to, as a third arm. "
                         "Overrides --orographic-arm; the default anchor is the "
                         "derived Earth reference and NOT the namelist fallback")
    ap.add_argument("--bare-z0", type=float, default=DEFAULT_BARE_Z0_M)
    ap.add_argument("--canopy-z0", type=float, default=DEFAULT_CANOPY_Z0_M)
    ap.add_argument("--forest-z0", type=float, default=DEFAULT_FOREST_Z0_M)
    ap.add_argument("--forest-fraction", type=float, default=None,
                    help="override the fraction implied by "
                         "model.land_albedo_source, the same override "
                         "build_surface_albedo.py takes, for a pair that "
                         "moves both fields together")
    ap.add_argument("--lakes", type=Path, default=None,
                    help="surface_water.nc; open water is smooth, so lake "
                         "regions take the ocean value before integration")
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    model = config["model"]
    # Deliberately from the grid, not from config: see builds.resolution_of.
    mesh = Export(args.mesh or mesh_export(config))
    grid_dir = args.grid or grid_export(config)
    # From the grid, not from config: the two differ exactly when someone
    # builds for another resolution, which is when the filename matters.
    resolution = resolution_of(grid_dir)

    # The reference height is linear in the lowest-level air temperature and this
    # step runs before any climatology exists to measure one from, so it is taken
    # at the MIDPOINT of the liquid-water bracket declared above and the
    # bracket's own effect on the answer is reported rather than hidden.
    z_ref_anchor = reference_height_m(0.5 * (CE_BRACKET_K[0] + CE_BRACKET_K[1]),
                                      config)

    # THE HIGH ARM OF THE DECLARED BRACKET, and the anchor. Derived, not
    # declared: the effective land roughness of the model's own Earth boundary
    # dataset at this rung, reduced in `ce` at this planet's reference height.
    # The namelist fallback travels with it so a reader sees what it is a
    # rounding of. Neither is consulted when the caller fixes the coefficient or
    # names its own target, and the reference cannot be derived at all for a rung
    # the model ships no dataset for -- so it is resolved HERE, before the mesh
    # is touched, and a rung it cannot serve refuses in a second rather than
    # after the reduction.
    reference = None
    if args.orographic_coefficient is None and args.target_mean is None:
        reference = earth_reference_land_z0(resolution, z_ref_anchor)

    rock = mesh.field("substrate_class").astype(int)
    is_land = mesh.surface_class == LAND
    area = mesh.cell_area.astype(np.float64)
    elev_km = mesh.field("elevation_km").astype(np.float64)

    codes = [r["code"] for r in mesh.manifest["lithology"]["rockClasses"]]
    barren = np.zeros(rock.shape, dtype=bool)
    for code in model.get("barren_rock_classes", []):
        if code in codes:
            barren |= rock == codes.index(code)

    # Surface term. The forest fraction comes from `model.land_albedo_source`,
    # the same key build_surface_albedo.py resolves its mode from, so roughness
    # and code 212 describe ONE land cover. It used to read
    # `model.forest_fraction_assumed`, a key the config does not carry, so it
    # silently used zero forest while 212 asserted half a canopy; that is
    # CLIM-36. An explicit --forest-fraction still overrides, for a sensitivity
    # pair. `modelled` is absent from the mapping because it does not imply a
    # fraction, it reads tree cover per cell, so it falls back to no blend here
    # until this reads that field too.
    mode = str(model.get("land_albedo_source", "lithology"))
    forest_fraction = (args.forest_fraction if args.forest_fraction is not None
                       else MODE_FOREST_FRACTION.get(mode))
    forest_fraction = float(forest_fraction) if forest_fraction else None
    canopy_z0 = args.canopy_z0
    if forest_fraction is not None:
        canopy_z0 = ((1.0 - forest_fraction) * args.canopy_z0
                     + forest_fraction * args.forest_z0)
    z0_surface = np.where(barren, args.bare_z0, canopy_z0)

    if args.lakes is not None:
        from netCDF4 import Dataset
        with Dataset(args.lakes) as lds:
            lake = np.asarray(lds["lake"][:]).astype(bool)
            lake_terrain = getattr(lds, "terrain_hash", None)
        if lake_terrain and lake_terrain != mesh.terrain_hash:
            raise SystemExit(
                f"lake solution is on terrain {lake_terrain[:16]}, mesh is "
                f"{mesh.terrain_hash[:16]}; re-run surface_water.py")
        z0_surface = np.where(lake & is_land, OCEAN_Z0_M, z0_surface)

    # Orographic term: the standard deviation of elevation among the mesh regions
    # inside each cell, which is subgrid relief by construction. `cell_moments`
    # in `lib/gridding.py` is the operator; the three bincounts used to be here.
    cells, nlat, nlon = region_cells(mesh, grid_dir)
    n = nlat * nlon
    sel = is_land
    _, var, counts, land_flat = cell_moments(cells, n, area, elev_km, sel)
    w = np.bincount(cells[sel], weights=area[sel], minlength=n)
    sigma_m = np.sqrt(var).reshape(nlat, nlon) * 1000.0
    counts = counts.reshape(nlat, nlon)

    land_cells = land_flat.reshape(nlat, nlon)
    weight = np.where(land_cells, w.reshape(nlat, nlon), 0.0)

    def land_mean(field):
        return float((field * weight).sum() / weight.sum())

    # THE SURFACE TERM IS REDUCED IN THE EXCHANGE COEFFICIENT, NOT IN THE LENGTH.
    #
    # The model applies one `ce` to the whole cell and the turbulent flux is
    # linear in it, so what a cell owes the atmosphere is the AREA MEAN OF `ce`
    # over its own surfaces. `ce` is not linear in `z0` -- it is
    # `k^2 / ln(z_ref/z0)^2` -- so averaging the lengths first and taking `ce`
    # of that is not the same reduction, and this world's land is exactly the
    # case where the two part company: the barren classes carry a roughness two
    # orders below the canopy, so a cell that is mostly playa with a canopy
    # minority has its exchange set by the minority once the lengths are mixed.
    # SPAT-7 measured it on the 10M mesh: negligible in the land mean, 24 times
    # inside the step's own `z_ref` bracket, but past that bracket on a few
    # percent of land area which GROWS with refinement, and reaching 30% on the
    # closed-basin floors the carve verdict integrates evaporation over. Those
    # cells are the reason this field exists. `analysis/spatial_reduction_gap.py`
    # is the measurement and `notes/audits/nonlinear-spatial-reductions.md` the
    # finding.
    #
    # So: `ce` per mesh region, area-mean over the cell's land, and the
    # roughness written out is the length that reproduces that mean, at the
    # `z_ref_anchor` set above.

    def region_z0(k_oro: float) -> np.ndarray:
        """Total roughness PER MESH REGION, the two terms in quadrature."""
        return np.sqrt(z0_surface ** 2
                       + (k_oro * sigma_m.reshape(-1)[cells]) ** 2)

    def reduce_ce(k_oro: float, z_ref: float):
        """`(cell z0, ce expectation, ce of the mean length, covered)`.

        The operator is `gridding.cell_expectation`, the NONLINEAR one, and it
        is called rather than reimplemented: it takes the LAW as a callable, so
        a caller cannot hand it a field that has already been reduced. The
        three bincounts that used to be here were the same reduction written a
        second time, which is what `lib/gridding.py` exists to stop.
        """
        ce_bar, z0_bar, covered = cell_expectation(
            cells, n, area, lambda v: neutral_ce(v, z_ref), region_z0(k_oro), sel)
        out = np.zeros(n)
        np.copyto(out, effective_length(ce_bar, z_ref), where=covered)
        return (out.reshape(nlat, nlon), ce_bar,
                neutral_ce(z0_bar, z_ref), covered)

    def effective_z0(k_oro: float, z_ref: float) -> np.ndarray:
        """Cell roughness whose `ce` is the area mean of the regions' own."""
        return reduce_ce(k_oro, z_ref)[0]

    # THE LOW ARM OF THE DECLARED BRACKET, computed whichever arm is written:
    # no orographic contribution, so the land mean is the one this world's own
    # lithology and land cover imply. It is the only land mean here that is
    # sourced end to end, and it is reported beside the anchored one so the two
    # are never seen apart.
    unrescaled_land_mean = land_mean(effective_z0(0.0, z_ref_anchor))

    if args.target_mean is not None:
        arm, target = "explicit", float(args.target_mean)
    elif args.orographic_coefficient is not None:
        arm, target = "fixed-coefficient", None
    elif args.orographic_arm == "none":
        arm, target = "none", None
    else:
        arm, target = "reference", reference["effective_land_z0_m"]

    # Solve the orographic coefficient so the land mean lands on the target.
    # 60 bisections resolve the coefficient to 1e-18 on [0, 1]; the 200 this
    # carried were free when the trial was arithmetic on the grid and are not
    # now that each one reduces the mesh.
    if args.orographic_coefficient is not None:
        k_oro = float(args.orographic_coefficient)
    elif target is None:
        k_oro = 0.0
    else:
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if land_mean(effective_z0(mid, z_ref_anchor)) < target:
                lo = mid
            else:
                hi = mid
        k_oro = 0.5 * (lo + hi)
        # A BISECTION THAT ENDS ON A BOUND HAS NOT SOLVED ANYTHING. The bracket
        # is [0, 1] and both ends are reachable: a target below the land mean
        # the cover alone gives drives it to zero, and a target the coefficient
        # cannot reach at all drives it to one. Either way the returned value is
        # the bound, the land mean is not the target, and the field is not what
        # the report says it is. Refuse instead of writing it.
        if not 1e-6 < k_oro < 1.0 - 1e-6:
            raise SystemExit(
                f"the orographic coefficient solved to {k_oro:.6g}, which is a "
                f"bound of the search bracket and not a solution: a land mean of "
                f"{target:.4f} m is not reachable on this grid, where the cover "
                f"alone gives {unrescaled_land_mean:.4f} m. Nothing was written.")

    z0, ce_expectation, ce_of_mean, covered = reduce_ce(k_oro, z_ref_anchor)
    field = np.where(land_cells, z0, OCEAN_Z0_M)

    # THE INVERSION IS AN IDENTITY AND IS CHECKED AS ONE. The cell roughness is
    # defined as the length whose `ce` is the area mean of the regions' own, so
    # putting it back through the law must return that mean. This has a right
    # answer rather than a plausible one, and it fails if the law and its
    # inverse ever stop being each other's -- which is the way a reduction
    # written twice goes wrong.
    round_trip = neutral_ce(z0.reshape(-1)[covered], z_ref_anchor)
    residual = float(np.max(np.abs(round_trip / ce_expectation[covered] - 1.0)))
    if residual > 1e-9:
        raise SystemExit(
            f"the ce inversion does not round-trip: worst relative residual "
            f"{residual:.3e} over {int(covered.sum())} covered cells. The "
            f"written roughness would not be the length that reproduces the "
            f"mean exchange coefficient. Nothing was written.")

    # THE JENSEN GAP IS THE FIELD'S OWN EVIDENCE THAT THE ORDER MATTERS.
    # `expectation - law(mean)` in `ce`, relative to the expectation, per cell.
    # It is what separates this reduction from the one this file used to do,
    # and reporting it means a later reader does not have to re-derive SPAT-7's
    # measurement to see whether the two orders part on THIS build and THIS
    # grid. A gap that reads zero everywhere would say the correction is inert
    # here, which is a finding and not a reason to drop the operator.
    gap = np.zeros(n)
    np.divide(ce_expectation - ce_of_mean, ce_expectation, out=gap,
              where=covered & (ce_expectation > 0.0))
    gap_land = gap[covered]

    # What the unknown reference height is worth in the field itself, reported
    # as a bracket because it cannot be verified before a climatology exists.
    z0_bracket = [land_mean(effective_z0(k_oro, reference_height_m(t, config)))
                  for t in CE_BRACKET_K]

    output = args.output or (INPUTS / resolution.lower()
                             / f"orogen_{resolution}_surf_{ROUGHNESS_CODE:04d}.sra")
    output.parent.mkdir(parents=True, exist_ok=True)
    write_sra(output, ROUGHNESS_CODE, field)

    z_ref = {t: reference_height_m(t, config) for t in CE_BRACKET_K}

    def ce(z0v, t_air):
        """Neutral bulk exchange coefficient at the lowest model level."""
        return neutral_ce(z0v, z_ref[t_air])

    def bracketed(z0v):
        lo, hi = (float(ce(z0v, t)) for t in CE_BRACKET_K)
        return [round(min(lo, hi), 6), round(max(lo, hi), 6)]

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "code": ROUGHNESS_CODE,
        "field": "dz0clim, total aerodynamic roughness length, metres",
        "terrain_hash": mesh.terrain_hash,
        "replaces_uniform": EXOPLASIM_DZ0LAND_M,
        "surface_z0_m": {"bare": args.bare_z0, "canopy": canopy_z0,
                         "forest": args.forest_z0,
                         "forest_fraction_used": forest_fraction},
        # THE ANCHOR, AND BOTH ENDS OF THE BRACKET IT SITS AT.
        #
        # `arm` is the end that was WRITTEN. `earth_reference` is the high end,
        # derived at build time from the model's own Earth boundary dataset at
        # this rung rather than declared, and it carries the namelist fallback
        # beside it so a reader can see what that number is a rounding of.
        # `unrescaled_land_mean_m` is the low end: what this world's lithology
        # and land cover give with no orographic contribution at all, which is
        # the only land mean here that is sourced end to end.
        #
        # The two are far apart and the difference is not cosmetic: `ce` is
        # logarithmic in the roughness, so the land turbulent exchange, and with
        # it land evaporation and land sensible heat, differs between the arms by
        # a factor the report states rather than implies. That is why both are
        # here on every build and why neither is presented as the answer.
        "anchor": {
            "arm": arm,
            "target_land_mean_m": target,
            "unrescaled_land_mean_m": round(unrescaled_land_mean, 5),
            "unrescaled_ce_bracket": bracketed(unrescaled_land_mean),
            "written_ce_bracket": bracketed(land_mean(z0)),
            "earth_reference": reference,
            # WHAT BOUNDS THE TRANSFER, and the one quantity the two planets can
            # be compared on directly. The reference dataset carries its surface
            # and orographic parts separately, so its land-cover roughness is
            # recoverable, and this world's low arm IS its land-cover roughness.
            # Neither is derived from the other. If they were decades apart the
            # transfer would be importing a land cover as well as a relief; that
            # they are close is what makes the high arm a bounded assumption
            # rather than an arbitrary one, and it is reported rather than
            # gated because no threshold on it was fixed in advance.
            "cover_roughness_ratio_to_earth": (
                None if reference is None
                else round(unrescaled_land_mean
                           / reference["surface_part_land_mean_m"], 4)),
            "note": "the high end is the effective land roughness of PlaSim's "
                    "own Earth boundary dataset, code 173 reduced over its land "
                    "in ce at this planet's reference height; the low end is "
                    "this world's own cover with no orographic term. The "
                    "namelist fallback dz0land is reported inside "
                    "earth_reference and is not the anchor.",
        },
        "orographic": {
            "coefficient_solved": round(k_oro, 6),
            "relation": "z0_oro = coefficient * stdev(elevation) within the cell",
            "subgrid_stdev_m": {
                "min": round(float(sigma_m[land_cells].min()), 2),
                "median": round(float(np.median(sigma_m[land_cells])), 2),
                "max": round(float(sigma_m[land_cells].max()), 2)},
            "mesh_regions_per_land_cell_median": int(np.median(counts[land_cells])),
            # SPAT-8 must hold this fixed across the ladder. It is SOLVED per
            # grid, so a difference between two rungs is part terrain and part
            # calibration, and the rung-matched reference does NOT close that:
            # measured on one build at both shipped rungs, the coefficient rises
            # from T21 to T42 against a fixed target and falls against the
            # reference, by comparable amounts in opposite directions.
            # `--orographic-coefficient` supplies one value for a whole ladder
            # comparison; the default still solves, so a single-rung build is
            # unchanged.
            "coefficient_was_solved": args.orographic_coefficient is None,
        },
        "surface_aggregation": {
            "operator": "gridding.cell_expectation, NONLINEAR, over the cell's "
                        "land population: ce evaluated per mesh region and then "
                        "area-averaged, inverted to the length that reproduces "
                        "that mean. lib/gridding.py owns the reduction "
                        "operators and this file calls one rather than "
                        "reimplementing it.",
            "law": "ce = karman^2 / ln(z_ref / z0)^2",
            "population": "surface_class == land",
            "reference_height_anchor_m": round(z_ref_anchor, 2),
            "land_mean_m_over_reference_air_bracket": [round(v, 5)
                                                       for v in z0_bracket],
            # What taking `ce` of the mean length instead would have cost,
            # per cell, at the anchor height. The sign is the law's curvature:
            # `d2ce/dz0^2` carries the factor `3 - ln(z_ref/z0)`, so `ce` is
            # CONCAVE in `z0` wherever the cell is more than about twenty
            # reference heights rougher than smooth, which is nearly all land
            # here, and the expectation therefore sits BELOW `ce` of the mean.
            # A cell that is mostly playa with a canopy minority is the extreme
            # of it, and that is the case this field exists for.
            "jensen_gap_in_ce_relative": {
                "min": float(gap_land.min()),
                "median": float(np.median(gap_land)),
                "max": float(gap_land.max()),
                "land_area_weighted_mean": land_mean(gap.reshape(nlat, nlon)),
            },
            "transfer": transfer_ledger(cells, n, area, sel),
        },
        "land_mean_m": round(land_mean(z0), 5),
        "land_min_m": round(float(z0[land_cells].min()), 6),
        "land_max_m": round(float(z0[land_cells].max()), 4),
        "exchange_coefficient": {
            "note": "ce = k^2 / ln(z_ref / z0)^2, the quantity roughness acts "
                    "through. z_ref is lib/lapse.py:reference_height_m at this "
                    "planet's gravity; each value is [low, high] over "
                    "reference_air_k, across which z_ref is linear.",
            "karman": KARMAN,
            "reference_air_k": list(CE_BRACKET_K),
            "reference_height_m": [round(z_ref[t], 2) for t in CE_BRACKET_K],
            "uniform_default": bracketed(EXOPLASIM_DZ0LAND_M),
            # ce rises with z0, so the roughest cell carries the largest ce.
            "this_field_land_max": bracketed(z0[land_cells].max()),
            "this_field_land_min": bracketed(z0[land_cells].min()),
        },
        "file": str(output),
    }
    report_path = output.parent / f"roughness_{resolution.lower()}_report.json"
    # Provenance stamp; lib/provenance.py owns the shape and the inert set.
    report.update(config_stamp(config, "exoplasim/scripts/build_surface_roughness.py"))
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"arm {arm}: land-mean z0 {land_mean(z0):.4f} m against the uniform "
          f"{EXOPLASIM_DZ0LAND_M} m fallback it replaces")
    print(f"  the field's own land mean with no orographic term is "
          f"{unrescaled_land_mean:.4f} m, the other end of the bracket")
    if reference is not None:
        print(f"  Earth reference from {reference['dataset']}: "
              f"{reference['effective_land_z0_m']:.4f} m in ce, "
              f"{reference['land_mean_of_length_m']:.4f} m in the length, "
              f"surface part {reference['surface_part_land_mean_m']:.4f} m")
    print(f"  range {z0[land_cells].min():.5f} to {z0[land_cells].max():.3f} m")
    print(f"  orographic coefficient solved to {k_oro:.5f}, subgrid stdev median "
          f"{np.median(sigma_m[land_cells]):.1f} m")
    lo_t, hi_t = CE_BRACKET_K
    print(f"  exchange coefficient {ce(z0[land_cells].min(), lo_t):.5f} "
          f"(flattest, barest) to {ce(z0[land_cells].max(), lo_t):.5f} "
          f"(roughest), against a uniform "
          f"{ce(EXOPLASIM_DZ0LAND_M, lo_t):.5f}, all at {lo_t} K")
    print(f"wrote {output}\n      {report_path}")


if __name__ == "__main__":
    main()
