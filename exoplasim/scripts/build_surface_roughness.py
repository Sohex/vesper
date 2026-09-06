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

**Under `model.land_albedo_source: modelled` the canopy is the one that grew.**
The forest fraction the other modes imply is an assumption about land cover,
and the albedo field stopped making it as soon as an LPJ-GUESS run existed to
read. While this field went on making it, one cell could be optically forest and
aerodynamically bare ground, and the surface energy balance was reporting on
two worlds at once. It now reads the same run's `fpc.out` through the same
reader and multiplies it by BIO-11's rootable fraction through the same
`rootable_cover_shares`, so both fields describe one cover. Solved water and dry
barren keep the roughness the native mesh gave them, and the orographic term is
taken on the cover roughness the mixture produces.

The cover mixture is taken in `ce` for the reason the region reduction is: tree,
grass and uncovered rootable ground are FRACTIONS of a cell's rootable area with
no position on the mesh, forest and bare ground are two orders apart in `z0`,
and `ce` is logarithmic in it, so mixing the lengths one scale further down
repeats the same error at the same place -- a canopy minority setting a
mostly-bare cell's exchange.

**It is integrated in `ce`, not in the length.** The model applies one exchange
coefficient to a whole cell and the turbulent flux is linear in it, so what the
cell owes the atmosphere is the area mean of `ce` over its own surfaces, and the
roughness written out is the length that reproduces that mean. `ce` is
logarithmic in `z0`, so the two orders are different reductions, and this
world's land is the case where they part: bare ground and canopy are two orders
apart in `z0`, so a cell that is mostly playa with a canopy minority has its
exchange set by the minority once the lengths are mixed. SPAT-7 measured the
difference on the 10M mesh: it moves the land-mean `ce` by about a tenth, past
this step's own `z_ref` bracket over most of land area, and reaches a quarter in
the top decile of barren share. Those cells are the closed-basin floors the
carve verdict integrates evaporation over, and they are the reason this field
exists. `notes/audits/nonlinear-spatial-reductions.md` section 3 carries the
measurement, and the height the average is taken at is Mason (1988)'s blending
height rather than the lowest model level -- see below.

`z0_orographic` is the part worth retaining the native mesh for. The 10M
fine-support reference contributes about 1,221 regions per global T42 cell on
average (and about 76 even at T170), enough to measure the distribution of
elevation *inside* the cell rather than infer it from resolved slope. Counts
over land vary and are reported by the builder; SPAT-3 is the shared artifact
for comparisons across the full T21/T42/T85/T127/T170 ladder. That subgrid
relief is what a gridded mean elevation cannot supply.
## The orographic term: DERIVED, with no free coefficient

The orographic contribution is not a coefficient this project picks and not one
it solves onto a reference. It is computed from this world's own subgrid slope
through the scheme large-scale models use for it, and every constant in that
scheme is read out of the papers that established it.

**The combination is additive in the DRAG COEFFICIENT, not in the length.**
Wood and Mason (1993) Eq. (33): the effective roughness of hilly ground is the
length whose drag coefficient at the pressure scale height equals the sum of
the pressure (form) drag coefficient of the hills and the skin-friction drag
coefficient of the ground between them,

    karman^2 / ln(z_m / z0_eff)^2 = Ca + karman^2 / ln(z_m / z0_cover)^2

so the two terms add where they act. `simba.f90:474` combines its vegetation and
topography terms in quadrature in the LENGTH instead; that is the model's own
convention and it is not the one the literature derives, so this file follows
the derivation.

**The form drag coefficient is quadratic in the subgrid SLOPE, not linear in the
relief.** Beljaars et al. (2004) Eq. (6), which is their Eq. (3) -- WM93's
surface pressure drag -- with the hills' base area set equal to the domain area
for mountainous terrain:

    Ca = 2 * alpha * beta * Cmd * theta_bar^2

`theta_bar^2` is the variance of the subgrid slope, `Cmd` the neutral drag
coefficient of the underlying cover at the pressure scale height, `beta` a shape
factor and `alpha` a shear parameter. This is the dependence that matters most
here: relief amplitude alone does not set the drag, because a given amount of
relief spread over a longer horizontal scale is a gentler and less resistive
surface. A scheme written as `coefficient * stdev(elevation)` has to carry that
missing horizontal scale inside its coefficient, which is why such a coefficient
cannot be sourced and drifts with the support.

**Nothing in it is fitted here.** `h_m`, the pressure scale height, and `l`, the
inner-layer depth, come from Beljaars Eq. (4); `alpha` from Eq. (5), which is
WM93's first-order-closure form; `Cmd` is the neutral drag coefficient at `h_m`
over the cell's own cover roughness; `beta = 1`, which WM93's Table 1 supports
for the parameter used here, since packed three-dimensional hills give effective
roughnesses of the same order as two-dimensional ones AT THE SAME `A/Sd`, and
`A/Sd` is the slope this scheme takes. `alpha` evaluated from Eq. (5) on this
world's land returns Beljaars' own assumed constant of 12 to two figures, which
is a check on the solvers rather than a coincidence.

**The slope is measured, not inferred.** `Export.local_slope_deg` fits a plane
through each mesh region and its neighbours, so the slope comes from the mesh at
its own spacing rather than from a gridded elevation. The variance per wind
direction is half the mean square gradient, the land being isotropic in slope to
the accuracy Beljaars section 2 finds for real terrain.

**The cell reduction is Mason (1988), and it fixes the height.** Mason defines
the area-average roughness as the one that reproduces the correct area-mean
surface stress, and shows it is obtained by averaging DRAG COEFFICIENTS at the
BLENDING HEIGHT `l_b` -- the height at which the flow is both in local
equilibrium with the surface beneath it and independent of horizontal position.
`l_b` follows from his Eq. (14), `l_b ln(l_b/z0)^2 = 2 karman^2 L_c` with
`L_c = L_D / 2pi` and `L_D` the horizontal scale over which the cover varies,
which here is the mesh spacing. That is an order of magnitude below the lowest
model level this file used to average at, and averaging too high understates the
weight a rough minority carries: Mason's whole result is that the average gives
extra weight to high roughnesses occupying small fractions of area. His
procedure is not sensitive to the exact `l_b`, only to its order, so one height
is solved for the grid and the per-cell spread is reported.

## What is declared and swept, and what is a limit rather than a bracket

**The declared bracket is on the scheme's inputs, and it is narrow.** Two ends,
both read out of the papers:

- The wavelength attributed to the mesh-resolved slope band. A plane fit through
  neighbours at spacing `d` reports the slope of features no shorter than about
  `2d`, so `2d` is the central attribution and `[d, 4d]` the bracket. It enters
  only through logarithms.
- The turbulence closure. WM93 measure the same pressure force at two closures:
  `5.9 theta^2` for the mixing-length (first-order) closure and `3.48 theta^2`
  for the second-order one, against Taylor et al.'s `3.5 theta^2`. Eq. (5) is
  the first-order form, so the second-order end scales `Ca` by `3.48/5.9`.

Both ends are computed and reported on every build. Together they move the land
mean by a few per cent, which is the point: the scheme's answer is insensitive
to everything about it that is not measured.

**What is NOT bracketed is the terrain the world does not have, and that is a
limit of the generator.** `Ca` is quadratic in slope, and slope variance is
dominated by the shortest wavelength present -- for a spectrum `F(k)` it is
`integral k^2 F(k) dk`, so the large scales contribute almost nothing. Orogen has
a measured terrain-information floor near 20 km, so this world carries no relief
below it; Beljaars' scheme, by contrast, is written for the band BELOW 5 km and
reaches it by extrapolating a power-law orographic spectrum fitted to United
States topography. Supplying that band here would be importing Earth's small-scale
topography under a derivation's name. So the derived roughness is what this
world's terrain supports, and the distance from Earth's is reported rather than
closed. `notes/audits/tuned-values.md` section 9 carries the argument.

**Lettau (1969) is carried as an independent cross-check and it does not agree at
this scale.** His `z0 = 0.5 h* s/S` gets a roughness from element geometry with
no fit, and his own worked example is orographic at exactly this scale --
Colorado's peaks at `h* = 1000 m`, `s = 5e6 m^2`, `S = 2e8 m^2`, giving 12.5 m.
Applied per mesh region here it runs one to two orders above the WM93 answer in
steep terrain, because it is purely geometric: at fixed slope it grows linearly
with the size of the elements, where the WM93 form saturates through the
logarithmic profile it is built on. Lettau flags his own continental example as
a borderline case for the surface Rossby numbers involved, and Beljaars assigns
scales above 5 km to gravity-wave and blocking schemes rather than to a
roughness at all. The builder computes both and reports the ratio; the WM93 form
is what is written, because it is the one with the boundary-layer physics in it
and the one the operational schemes descend from.

## The Earth reference is a COMPARISON, and no longer a target

The vendored model ships PlaSim's own boundary dataset for the model's Earth
configuration at two rungs, `N032` and `N064`: code 172 the land mask, code 173
the total roughness `dz0clim`, code 1730 its topography-only part. Reducing code
173 over that dataset's land in `ce` returns the effective land roughness the
model carries when it is given a real roughness map instead of the uniform
fallback, and at `N032` that reproduces the namelist fallback `dz0land` to within
a tenth. So `dz0land` is a rounding of the model's Earth land roughness, and the
derivation it never carried in the source lives in a boundary dataset.

`earth_reference_land_z0()` computes it at build time rather than quoting it, and
every build reports the derived land mean beside it. What it is NOT is something
to solve onto. Earth's land roughness is a sound derivation for the wrong planet:
this world's relief is scaled by its own gravity, its lithology and land cover
are its own, and its terrain has a coarser information floor. Anchoring on
Earth's number would report the distance as zero by construction and discard the
one land mean that is derived end to end.

The distance is reported instead. The Earth dataset's own cover term, recovered
from the quadrature, and this world's cover term agree closely, so the two
worlds' LAND COVER roughness is not what separates them; the whole of the gap is
in the relief term, which is where the terrain floor sits.

`--orographic-arm reference` still solves a coefficient on the old
`coefficient * stdev(elevation)` relation so the Earth-anchored field can be
rebuilt for a paired run that prices the change. It is a DIAGNOSTIC ARM: it is
not the default, and no climatology lineage may rest on it. `--orographic-arm
none` writes the cover term alone. `--orographic-coefficient` and `--target-mean`
override for a fixed ladder or a third arm.

**The derived arm needs no ladder discipline.** The old solve was per grid, so
two rungs differed by their terrain and by their calibration at once and SPAT-8
had to pass one coefficient to every rung. The derived orographic term is a
property of the land and the mesh, not of the grid, so it does not move with the
rung: the land means at the two rungs the model ships a reference for differ by a
couple of per cent, against a factor of nearly two for the Earth-anchored arm.

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
from build_surface_albedo import (MODE_FOREST_FRACTION, foliar_cover_on_grid,
                                  read_foliar_cover, rootable_cover_shares)
from rootable import read_rootable_partition
from gridding import (cell_expectation, cell_fraction, cell_mean, cell_moments,
                      gaussian_grid, region_cells, transfer_ledger)
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

# landmod.f90:51, the uniform FALLBACK this replaces, and it is not a target.
# `configure()` falls back to it over all land when no code 173 is staged. It is
# reported beside the Earth reference below, which it is a rounding of.
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

# -- the orographic form drag scheme, and every constant in it ---------------
#
# Wood and Mason (1993), QJRMS 119(514), 1233-1267, 10.1002/qj.49711951402, for
# the combination and the effective roughness; Beljaars, Brown and Wood (2004),
# QJRMS 130(599), 1327-1347, 10.1256/qj.03.73, for the operational form of the
# drag and for the constants of Eqs (4) to (6); Mason (1988), QJRMS 114(480),
# 399-420, 10.1002/qj.49711448007, for what an areally-averaged roughness IS and
# for the height its average is taken at; Lettau (1969), J. Appl. Meteorol. 8,
# 828-832, 10.1175/1520-0450(1969)008<0828:NOARPE>2.0.CO;2, for the independent
# geometric cross-check. `references/INDEX.md` records which sections were read.

# WM93's shape factor. 1 for two-dimensional hills, and WM93's Table 1 supports
# it for the parameter this scheme takes: packed three-dimensional hills give
# effective roughnesses "of the same order of magnitude" as two-dimensional ones
# at the SAME A/Sd, and A/Sd is the slope. It is the isolated-hill case that
# departs, and this world's land is packed rather than isolated everywhere the
# term is large.
WM93_SHAPE_FACTOR = 1.0

# THE CLOSURE BRACKET, and it is measured rather than assumed. WM93 give the
# same pressure force at two levels of turbulence closure for periodic
# two-dimensional hills at A/z0 = 1e4: 5.9 theta^2 for the mixing-length
# (first-order) closure and 3.48 theta^2 for the second-order one, against
# Taylor et al. (1989)'s 3.5 theta^2. Eq (5) below is WM93's first-order form,
# so the second-order end of the bracket scales the form drag by their ratio.
WM93_CLOSURE_FIRST_ORDER = 5.9
WM93_CLOSURE_SECOND_ORDER = 3.48

# THE BAND BRACKET. A plane fit through a mesh region and its neighbours at
# spacing `d` reports the slope of features no shorter than about 2d, which is
# the central attribution; the ends are one and four spacings. It enters the
# scheme only through logarithms, so the bracket is narrow by construction and
# is reported rather than argued.
RESOLVED_BAND_MULTIPLE = 2.0
RESOLVED_BAND_BRACKET = (1.0, 4.0)

# Beljaars' own simplification point, and the check that the solvers below
# reproduce it: with c_m = 0.1 and Cmd = 0.005 the paper fixes alpha = 12.
BELJAARS_CM = 0.1
BELJAARS_CMD = 0.005
BELJAARS_ALPHA = 12.0

# ECMWF limits effective roughness to 100 m to avoid numerical problems
# (Beljaars et al. 2004, section 1). It is the only cap in the literature this
# tree holds, so it is the one used, and the area it binds on is reported.
ECMWF_Z0_CAP_M = 100.0

# Lettau's numerical factor in `z0 = 0.5 h* s/S`, which he identifies as "the
# average drag coefficient of the characteristic individual obstacle of
# silhouette area s". Not a fit: his bushel-basket experiments put the relation
# within about 25% of wind-profile analysis.
LETTAU_DRAG_FACTOR = 0.5


def _fixed_point(step, guess, residual, name: str, tol: float = 1e-8):
    """Iterate `step` to convergence and REFUSE if it did not converge.

    Every length in Beljaars Eq (4) and in Mason Eq (14) is defined implicitly,
    by an equation with the unknown inside a logarithm. A fixed count of
    iterations is a guess that the recursion contracted; the residual is the
    equation itself and it has a right answer, so it is checked rather than
    assumed. A scheme whose lengths silently did not solve would still produce
    an ordinary-looking field.
    """
    x = np.asarray(guess, dtype=np.float64)
    for _ in range(200):
        nxt = step(x)
        if np.all(np.abs(nxt - x) <= tol * np.maximum(np.abs(nxt), 1e-12)):
            x = nxt
            break
        x = nxt
    worst = float(np.max(np.abs(residual(x))))
    if not worst < 1e-6:
        raise SystemExit(
            f"{name} did not solve: worst relative residual {worst:.3e} on its "
            f"own defining equation. The scheme's lengths would be wrong and the "
            f"field would still look ordinary. Nothing was written.")
    return x


def pressure_scale_height(wavelength_m, z0):
    """`h_m`, from Beljaars Eq (4): `h_m ln^(1/2)(h_m/z0) = lambda/4`.

    The height at which the square of the velocity scales the pressure field,
    and the level WM93 form the effective roughness at.
    """
    lam = np.asarray(wavelength_m, dtype=np.float64)
    z0 = np.asarray(z0, dtype=np.float64)

    def ln(x):
        return np.log(np.maximum(x / z0, 1.0 + 1e-9))

    return _fixed_point(lambda h: 0.25 * lam / np.sqrt(ln(h)),
                        BELJAARS_CM * lam,
                        lambda h: h * np.sqrt(ln(h)) / (0.25 * lam) - 1.0,
                        "the pressure scale height h_m")


def inner_layer_depth(wavelength_m, z0):
    """`l`, from Beljaars Eq (4): `l ln(l/z0) = karman^2 lambda / 2`."""
    lam = np.asarray(wavelength_m, dtype=np.float64)
    z0 = np.asarray(z0, dtype=np.float64)
    rhs = KARMAN ** 2 * lam / 2.0

    def ln(x):
        return np.log(np.maximum(x / z0, 1.0 + 1e-9))

    return _fixed_point(lambda l: rhs / ln(l), rhs / 10.0,
                        lambda l: l * ln(l) / rhs - 1.0,
                        "the inner-layer depth l")


def shear_parameter(h_m, l, z0):
    """`alpha`, Beljaars Eq (5), which is WM93's first-order-closure form.

    `alpha = 2 [ln(h_m/z0) / ln(l/z0)]^4 + ln(l/z0)`. Beljaars assumes it
    constant at 12; here it is evaluated, and it returns 12 on this world's land
    to two figures, which is what says the two lengths above solved correctly.
    """
    a = np.log(np.maximum(np.asarray(h_m) / z0, 1.0 + 1e-9))
    b = np.log(np.maximum(np.asarray(l) / z0, 1.0 + 1e-9))
    return 2.0 * (a / b) ** 4 + b


def blending_height(cover_scale_m, z0):
    """`l_b`, Mason (1988) Eq (14): `l_b ln(l_b/z0)^2 = 2 karman^2 L_c`.

    The height at which the flow is both in local equilibrium with the surface
    beneath it and independent of horizontal position, and therefore the height
    at which drag coefficients must be averaged to get the roughness that
    reproduces the correct area-mean surface stress. `L_c = L_D / 2pi` is
    Mason's own scaling of the periodicity length `L_D` to a differentiation
    scale; `L_D` here is the spacing at which the cover varies, the mesh's.
    """
    lc = np.asarray(cover_scale_m, dtype=np.float64) / (2.0 * np.pi)
    rhs = 2.0 * KARMAN ** 2 * lc

    def ln(x):
        return np.log(np.maximum(x / z0, 1.0 + 1e-9))

    return _fixed_point(lambda lb: rhs / ln(lb) ** 2,
                        np.full(np.shape(rhs), 20.0) if np.ndim(rhs) else 20.0,
                        lambda lb: lb * ln(lb) ** 2 / rhs - 1.0,
                        "the blending height l_b")


def orographic_form_drag(slope_variance, z_m, h_m, l, z0, closure: float = 1.0):
    """`Ca`, Beljaars Eq (6): `2 alpha beta Cmd theta_bar^2`.

    The pressure force per unit area on the hills, as a drag coefficient
    referenced to the velocity at `z_m`, which is the quantity WM93 Eq (33)
    adds to the skin friction there. `z_m` is the pressure scale height where
    the hills are lower than it and the hill top where they are higher, WM93's
    own definition; `alpha` is a property of the pressure scale height itself
    and so takes `h_m` in either case. `closure` scales the result to the
    second-order end of the bracket; 1 is WM93's first-order form.
    """
    cmd = neutral_ce(z0, z_m)
    alpha = shear_parameter(h_m, l, z0)
    return closure * 2.0 * alpha * WM93_SHAPE_FACTOR * cmd * slope_variance, cmd, alpha


def lettau_z0(obstacle_height_m, silhouette_area_m2, lot_area_m2):
    """Lettau (1969) Eq (1): `z0 = 0.5 h* s/S`.

    Roughness from element geometry, with the 0.5 identified as the obstacles'
    own drag coefficient rather than fitted. Carried as an INDEPENDENT check on
    the WM93 answer: it shares no constant with it and reaches a roughness by a
    different argument, so agreement is evidence and disagreement is a finding.
    """
    return (LETTAU_DRAG_FACTOR * np.asarray(obstacle_height_m, dtype=np.float64)
            * np.asarray(silhouette_area_m2, dtype=np.float64)
            / np.asarray(lot_area_m2, dtype=np.float64))


def check_against_published_values() -> dict:
    """The scheme against numbers its own papers print. Runs on every build.

    Three checks, each with a right answer rather than a plausible one:

    - Beljaars fixes `alpha = 12` and `Cmd = 0.005` with `h_m = 0.1 lambda`.
      Pick the wavelength at which `Cmd` is 0.005 over a given cover roughness,
      solve Eq (4) there, and both `h_m/lambda` and `alpha` from Eq (5) must
      come back at the paper's values. This fails if either implicit length is
      solved wrongly, if Eq (5) is transcribed wrongly, or if `Cmd` is formed
      with the wrong power of the logarithm.
    - Lettau's own worked orographic example: Colorado's peaks at `h* = 1000 m`,
      `s = 5e6 m^2`, `S = 2e8 m^2` give `z0 = 12.5 m`.
    - WM93 Eq (33) with no hills must return the cover roughness exactly, which
      is the identity the additive combination rests on.
    """
    z0 = 0.03
    # Cmd = karman^2 / ln^2(h_m/z0) = 0.005 fixes h_m/z0, and h_m = 0.1 lambda
    # then fixes the wavelength Beljaars' constants belong to.
    h_over_z0 = np.exp(KARMAN / np.sqrt(BELJAARS_CMD))
    lam = float(h_over_z0 * z0 / BELJAARS_CM)
    h_m = float(pressure_scale_height(np.array([lam]), np.array([z0]))[0])
    l = float(inner_layer_depth(np.array([lam]), np.array([z0]))[0])
    alpha = float(shear_parameter(np.array([h_m]), np.array([l]),
                                  np.array([z0]))[0])
    cmd = float(neutral_ce(np.array([z0]), h_m)[0])
    # 10%, and the tolerance is a statement about the paper rather than about
    # this code: `c_m = 0.1` and `Cmd = 0.005` are the round numbers Beljaars
    # SIMPLIFIES Eq (4) to, not its exact solution, so agreement closer than a
    # few per cent would say the wrong thing. A solver that has not converged,
    # or an equation transcribed with the wrong power or the wrong quarter,
    # misses by orders of magnitude and not by a tenth.
    if not abs(h_m / lam - BELJAARS_CM) < 0.10 * BELJAARS_CM:
        raise SystemExit(
            f"Beljaars Eq (4) returns h_m/lambda = {h_m / lam:.4f} where the "
            f"paper's simplification fixes {BELJAARS_CM}. The pressure scale "
            f"height is not being solved. Nothing was written.")
    if not abs(alpha / BELJAARS_ALPHA - 1.0) < 0.10:
        raise SystemExit(
            f"Beljaars Eq (5) returns alpha = {alpha:.3f} where the paper fixes "
            f"{BELJAARS_ALPHA}. Nothing was written.")
    colorado = float(lettau_z0(1000.0, 5.0e6, 2.0e8))
    if not abs(colorado - 12.5) < 0.05:
        raise SystemExit(
            f"Lettau Eq (1) returns {colorado:.3f} m on his own Colorado "
            f"example, which the paper works out at 12.5 m. Nothing was written.")
    # WM93 Eq (33) with Ca = 0 is an identity on the cover roughness.
    flat = float(effective_length(neutral_ce(np.array([z0]), h_m), h_m)[0])
    if not abs(flat / z0 - 1.0) < 1e-9:
        raise SystemExit(
            f"WM93 Eq (33) with no hills returns {flat:.6g} m over a cover "
            f"roughness of {z0} m; the combination is not an identity where it "
            f"has to be. Nothing was written.")
    return {"beljaars_hm_over_lambda": round(h_m / lam, 5),
            "beljaars_alpha": round(alpha, 4),
            "beljaars_cmd": round(cmd, 6),
            "lettau_colorado_z0_m": round(colorado, 4),
            "published": {"hm_over_lambda": BELJAARS_CM, "alpha": BELJAARS_ALPHA,
                          "cmd": BELJAARS_CMD, "lettau_colorado_z0_m": 12.5}}


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
    so that a change under `vendor/exoplasim` moves the comparison instead of
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
    ap.add_argument("--orographic-arm",
                    choices=("derived", "reference", "none"), default="derived",
                    help="which orographic term to WRITE. 'derived', the "
                         "default, computes it from this world's own subgrid "
                         "slope through WM93 Eq (33) and Beljaars Eq (6) with "
                         "no free coefficient; 'reference' solves the old "
                         "coefficient * stdev(elevation) relation so the land "
                         "mean lands on the model's Earth boundary dataset, "
                         "which is a DIAGNOSTIC arm for pricing the change and "
                         "not a lineage; 'none' writes the cover term alone. "
                         "All three land means are REPORTED whichever is written")
    ap.add_argument("--closure", choices=("first-order", "second-order"),
                    default="first-order",
                    help="which end of WM93's measured closure bracket the "
                         "derived arm writes. Both are reported")
    ap.add_argument("--band-multiple", type=float,
                    default=RESOLVED_BAND_MULTIPLE,
                    help="the wavelength attributed to the mesh-resolved slope, "
                         "in mesh spacings. Both ends of the declared bracket "
                         "are reported whichever is written")
    ap.add_argument("--target-mean", type=float, default=None,
                    help="area-weighted land mean to solve onto, as a further "
                         "arm. Overrides --orographic-arm, and like the "
                         "reference arm it is a DIAGNOSTIC: the default arm "
                         "derives the land mean rather than targeting one")
    ap.add_argument("--bare-z0", type=float, default=DEFAULT_BARE_Z0_M)
    ap.add_argument("--canopy-z0", type=float, default=DEFAULT_CANOPY_Z0_M)
    ap.add_argument("--forest-z0", type=float, default=DEFAULT_FOREST_Z0_M)
    ap.add_argument("--forest-fraction", type=float, default=None,
                    help="override the fraction implied by "
                         "model.land_albedo_source, the same override "
                         "build_surface_albedo.py takes, for a pair that "
                         "moves both fields together")
    ap.add_argument("--vegetation", type=Path, default=None,
                    help="fpc.out from an accepted LPJ-GUESS run. Required "
                         "under model.land_albedo_source: modelled, which is "
                         "the mode in which the canopy this field describes is "
                         "the canopy the albedo field describes")
    ap.add_argument("--vegetation-peer", type=Path, action="append", default=[],
                    help="comparable LPJ run for the equilibrium reducer's "
                         "stochastic spread; repeat")
    ap.add_argument("--rootable", type=Path, default=None,
                    help="BIO-11 rootable-fraction artifact; default the "
                         "active build and rung")
    ap.add_argument("--climatology", type=Path, default=None,
                    help="the climatology the LPJ-GUESS driver was built from, "
                         "for its coordinate labels. Required with --vegetation")
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

    # THE EARTH COMPARISON. Derived, not declared: the effective land roughness
    # of the model's own Earth boundary dataset at this rung, reduced in `ce` at
    # this planet's reference height. The namelist fallback travels with it so a
    # reader sees what it is a rounding of.
    #
    # It is a COMPARISON and not a target, so a rung the
    # model ships no dataset for is no longer a reason to refuse: the derived
    # arm needs nothing from it. It is still resolved HERE, before the mesh is
    # touched, so the `reference` arm -- which does solve onto it -- refuses in
    # a second rather than after the reduction.
    wants_reference = (args.orographic_arm == "reference"
                       and args.orographic_coefficient is None
                       and args.target_mean is None)
    reference = None
    if wants_reference or resolution in EARTH_REFERENCE_RUNGS:
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
    # pair.
    #
    # `modelled` is absent from the mapping because it does not imply a
    # fraction: it reads tree and grass cover PER CELL from the LPJ-GUESS run,
    # and this file reads the same field through the same reader below. That is
    # BIO-16 and it is what closes the aerodynamic half of the vegetation
    # feedback. Until it did, this field asserted one canopy over all land while
    # the albedo field carried another, so a cell could be optically forest and
    # aerodynamically bare ground, and the surface energy balance was reporting
    # on two different worlds.
    mode = str(model.get("land_albedo_source", "lithology"))
    modelled_cover = mode == "modelled"
    if args.vegetation is not None and not modelled_cover:
        raise SystemExit(
            f"--vegetation was given but model.land_albedo_source is {mode!r}, "
            "so the albedo field is not reading the LPJ run either. Deriving "
            "the roughness from a canopy the albedo does not have is the "
            "mismatch this arm exists to close, in the other direction.")
    if modelled_cover and args.vegetation is None:
        raise SystemExit(
            "model.land_albedo_source is 'modelled', so the roughness comes "
            "from the same LPJ-GUESS cover the albedo does. Pass --vegetation "
            "<run>/fpc.out and --climatology, the one its driver was built "
            "from. There is no per-mode fraction to fall back to and falling "
            "back is what put two land covers in one surface.")
    forest_fraction = (args.forest_fraction if args.forest_fraction is not None
                       else MODE_FOREST_FRACTION.get(mode))
    forest_fraction = float(forest_fraction) if forest_fraction else None
    canopy_z0 = args.canopy_z0
    if forest_fraction is not None:
        canopy_z0 = ((1.0 - forest_fraction) * args.canopy_z0
                     + forest_fraction * args.forest_z0)
    z0_surface = np.where(barren, args.bare_z0, canopy_z0)

    lake_mask = np.zeros(is_land.shape, dtype=bool)
    if args.lakes is not None:
        from netCDF4 import Dataset
        with Dataset(args.lakes) as lds:
            lake = np.asarray(lds["lake"][:]).astype(bool)
            lake_terrain = getattr(lds, "terrain_hash", None)
        if lake_terrain and lake_terrain != mesh.terrain_hash:
            raise SystemExit(
                f"lake solution is on terrain {lake_terrain[:16]}, mesh is "
                f"{mesh.terrain_hash[:16]}; re-run surface_water.py")
        lake_mask = lake & is_land
        z0_surface = np.where(lake_mask, OCEAN_Z0_M, z0_surface)

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

    # -- THE MODELLED CANOPY, read from the run the albedo field reads --------
    #
    # `cover_weights` is None in every mode that implies one cover for all land,
    # and the reduction below is then exactly what it always was. Under
    # `modelled` it is the four populations a land cell is made of -- tree,
    # grass, uncovered rootable ground, and the non-rootable water and barren
    # the mesh placed -- each with the share of the cell it occupies.
    #
    # The shares come from `build_surface_albedo.rootable_cover_shares`, the one
    # derivation of how much canopy a cell carries, not from a second reading of
    # the same run. LPJ FPC is conditional on the rootable environment and is
    # multiplied by BIO-11's rootable fraction exactly once, there.
    cover_weights = None
    # BIO-11's three populations, on the mesh. A lake over barren ground is a
    # lake, which is the order the albedo builder's native masks take too.
    sel_water = lake_mask & is_land
    sel_barren = barren & is_land & ~sel_water
    sel_root = is_land & ~barren & ~sel_water
    sel_non = sel_water | sel_barren
    vegetation_provenance = None
    if modelled_cover:
        if args.climatology is None or not args.climatology.is_file():
            raise SystemExit(
                "--vegetation needs --climatology, the same one the LPJ-GUESS "
                "driver was built from. Its coordinate labels are the run's "
                "own; the export's planet.nc labels the same grid differently "
                "and matching against it returns zero cells.")
        from netCDF4 import Dataset as _Dataset
        with _Dataset(args.climatology) as data:
            lat_axis = np.asarray(data["lat"][:], dtype=float)
            lon_axis = np.asarray(data["lon"][:], dtype=float)
        table, equilibrium_window = read_foliar_cover(
            args.vegetation, args.vegetation_peer)
        tree_fpc, grass_fpc, matched = foliar_cover_on_grid(
            table, lat_axis, lon_axis, land_cells)
        partition, rootable_provenance = read_rootable_partition(
            config, lat_axis, lon_axis, args.rootable, land=land_cells)
        rootable = partition["rootable"]

        # The same check the albedo builder makes, for the same reason: the
        # artifact is what the algebra uses, and the native masks supplied to
        # THIS builder have to be the population it was built from. Compositing
        # two worlds is how a cell ends up with more surface than it has.
        native = {
            "rootable": cell_fraction(cells, n, area, sel_root, sel)[0],
            "water": cell_fraction(cells, n, area, sel_water, sel)[0],
            "barren": cell_fraction(cells, n, area, sel_barren, sel)[0],
        }
        partition_residual = max(
            float(np.max(np.abs(native[name].reshape(nlat, nlon)[land_cells]
                                - partition[name][land_cells])))
            for name in partition)
        if partition_residual > 2.0e-6:
            raise SystemExit(
                "the lake/barren inputs supplied to surface roughness disagree "
                f"with BIO-11's rootable partition by {partition_residual:.3e}; "
                "rebuild both from the same derived surface rather than "
                "compositing two worlds")

        tree_share, grass_share, _ = rootable_cover_shares(
            rootable, tree_fpc, grass_fpc)
        bare_share = np.clip(rootable - tree_share - grass_share, 0.0, None)
        nonroot_share = np.clip(1.0 - rootable, 0.0, 1.0)
        cover_weights = [
            (tree_share.reshape(-1), args.forest_z0),
            (grass_share.reshape(-1), args.canopy_z0),
            (bare_share.reshape(-1), args.bare_z0),
            (nonroot_share.reshape(-1), None),
        ]
        # THE PARTITION IS AN IDENTITY AND IS CHECKED AS ONE. Tree, grass, what
        # neither covers, and the non-rootable remainder are the whole of a land
        # cell. If they do not sum to one the reduction below silently drops or
        # double-counts area, and the roughness written out is an average over a
        # cell that is not this one. There is a right answer here rather than a
        # plausible one.
        weight_sum = (tree_share + grass_share + bare_share + nonroot_share)
        worst = float(np.max(np.abs(weight_sum[land_cells] - 1.0)))
        if worst > 1e-9:
            raise SystemExit(
                f"the cover partition does not sum to one: worst residual "
                f"{worst:.3e} over {int(land_cells.sum())} land cells. Nothing "
                "was written.")

        vegetation_provenance = {
            "vegetation": str(args.vegetation),
            "climatology": str(args.climatology),
            "cells_matched": matched,
            "equilibrium_window": equilibrium_window,
            "rootable_surface": rootable_provenance,
            "rootable_partition_max_absolute_residual": partition_residual,
            "land_mean_tree_share": land_mean(tree_share),
            "land_mean_grass_share": land_mean(grass_share),
            "note": "tree, grass and uncovered rootable ground are mixed in ce "
                    "over the cell's rootable regions; solved water and dry "
                    "barren keep the roughness the native mesh gave them",
        }

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
    # roughness written out is the length that reproduces that mean.
    #
    # AT THE BLENDING HEIGHT, which is Mason (1988)'s result and not this
    # file's choice. The average of drag coefficients that reproduces the
    # correct area-mean surface stress is the one taken at the height where the
    # flow is both in equilibrium with the surface under it and independent of
    # horizontal position, and that height follows from the scale over which the
    # cover varies -- the mesh spacing -- not from where the model happens to
    # put its lowest level. Averaging an order of magnitude too high understates
    # the weight a rough minority carries, and giving that minority extra weight
    # is the whole of Mason's finding. The height depends on the roughness it is
    # averaging, so it is iterated, as he prescribes.

    def region_z0_from(surface_z0: np.ndarray, k_oro: float) -> np.ndarray:
        """Total roughness PER MESH REGION, the two terms in quadrature."""
        return np.sqrt(np.asarray(surface_z0, dtype=np.float64) ** 2
                       + (k_oro * sigma_m.reshape(-1)[cells]) ** 2)

    def region_z0(k_oro: float) -> np.ndarray:
        return region_z0_from(z0_surface, k_oro)

    def reduce_ce(k_oro: float, z_ref: float):
        """`(cell z0, ce expectation, ce of the mean length, covered)`.

        The operator is `gridding.cell_expectation`, the NONLINEAR one, and it
        is called rather than reimplemented: it takes the LAW as a callable, so
        a caller cannot hand it a field that has already been reduced. The
        three bincounts that used to be here were the same reduction written a
        second time, which is what `lib/gridding.py` exists to stop.

        WITH A MODELLED CANOPY THE MIXTURE IS ALSO IN `ce`. The cover shares
        inside a cell's rootable ground -- tree, grass and what neither covers
        -- are FRACTIONS with no position on the mesh, so mixing their lengths
        first repeats, one scale down, exactly the error mixing region lengths
        makes: forest and bare ground are two orders apart in `z0`, and `ce` is
        logarithmic in it, so a mostly-bare cell with a canopy minority has its
        exchange set by the minority once the lengths are mixed. Each cover
        takes the same reduction over the same rootable regions, and the four
        expectations are combined with the cover shares, which is linear and is
        where the mixture belongs.
        """
        law = lambda v: neutral_ce(v, z_ref)
        if cover_weights is None:
            ce_bar, z0_bar, covered = cell_expectation(
                cells, n, area, law, region_z0(k_oro), sel)
        else:
            ce_bar = np.zeros(n)
            z0_bar = np.zeros(n)
            for weight, surface in cover_weights:
                population = sel_non if surface is None else sel_root
                values = (region_z0(k_oro) if surface is None
                          else region_z0_from(np.full(z0_surface.shape, surface),
                                              k_oro))
                part, mean, part_covered = cell_expectation(
                    cells, n, area, law, values, population)
                ce_bar += weight * np.where(part_covered, part, 0.0)
                z0_bar += weight * np.where(part_covered, mean, 0.0)
            covered = covered_land
        out = np.zeros(n)
        np.copyto(out, effective_length(ce_bar, z_ref), where=covered)
        return (out.reshape(nlat, nlon), ce_bar,
                neutral_ce(z0_bar, z_ref), covered)

    def effective_z0(k_oro: float, z_ref: float) -> np.ndarray:
        """Cell roughness whose `ce` is the area mean of the regions' own."""
        return reduce_ce(k_oro, z_ref)[0]

    # THE MIXTURE IS CHECKED AGAINST THE REDUCTION IT GENERALISES. Give every
    # cover inside the rootable ground the same roughness and the four-part
    # expectation must equal a single-cover reduction over one array carrying
    # that roughness on rootable regions and the mesh's own value elsewhere.
    # That is an identity with a right answer -- the weights are a partition and
    # the law is the same law -- and it fails if the populations, the weights or
    # the combination ever stop being each other's complement. It is one extra
    # reduction on a step that costs minutes.
    if cover_weights is not None:
        control_z0 = 0.1
        control_surface = np.where(sel_root, control_z0, z0_surface)
        flat = [(weight, None if surface is None else control_z0)
                for weight, surface in cover_weights]
        mixed = np.zeros(n)
        for weight, surface in flat:
            population = sel_non if surface is None else sel_root
            values = (region_z0_from(control_surface, 0.0) if surface is None
                      else region_z0_from(np.full(z0_surface.shape, surface), 0.0))
            part, _, part_covered = cell_expectation(
                cells, n, area, lambda v: neutral_ce(v, 10.0), values, population)
            mixed += weight * np.where(part_covered, part, 0.0)
        direct, _, direct_covered = cell_expectation(
            cells, n, area, lambda v: neutral_ce(v, 10.0),
            region_z0_from(control_surface, 0.0), sel)
        worst = float(np.max(np.abs(
            mixed[direct_covered] / direct[direct_covered] - 1.0)))
        if worst > 1e-9:
            raise SystemExit(
                f"the cover mixture does not reproduce the single-cover "
                f"reduction on a uniform cover: worst relative residual "
                f"{worst:.3e}. The four populations are not a partition of the "
                "cell's land, or the weights are not their areas. Nothing was "
                "written.")

    # The mesh spacing, per cell: the scale over which the cover varies, which
    # is what sets the blending height, and the scale the resolved slope belongs
    # to, which is what sets the wavelength. Centre-to-centre distance of
    # hexagonal cells of the region's own area.
    spacing_region_m = 1000.0 * np.sqrt(2.0 / np.sqrt(3.0)) * np.sqrt(area)
    spacing_cell, cell_land_area, covered_land = cell_mean(
        cells, n, area, spacing_region_m, sel)
    cover_scale_m = float(np.average(spacing_cell[covered_land],
                                     weights=cell_land_area[covered_land]))

    # Mason (1988) Eq (14), iterated on the roughness it averages. ONE height
    # for the grid: his own result is that the procedure needs the ORDER of the
    # blending height and not its exact value, and this world is sampled on a
    # quasi-uniform mesh, so the per-cell spread is reported rather than carried
    # through a per-cell reduction the shared operator does not take.
    z_blend = 20.0
    for _ in range(40):
        trial = land_mean(effective_z0(0.0, z_blend))
        nxt = float(blending_height(cover_scale_m, trial))
        if abs(nxt - z_blend) <= 1e-9 * nxt:
            z_blend = nxt
            break
        z_blend = nxt
    else:
        raise SystemExit(
            f"the blending height did not settle: {z_blend:.4f} m still moving "
            f"after 40 passes. Mason's average is taken at it, so nothing was "
            f"written.")

    # THE COVER TERM, computed whichever arm is written: no orographic
    # contribution at all, so the land mean is the one this world's own
    # lithology and land cover imply.
    unrescaled_land_mean = land_mean(effective_z0(0.0, z_blend))

    # -- THE DERIVED OROGRAPHIC TERM -----------------------------------------
    #
    # This is the arm with no free coefficient in it. Every constant comes from
    # the papers; the only measured input is this world's own subgrid slope.
    published = check_against_published_values()

    # The slope from the mesh's own plane fit, and its variance PER WIND
    # DIRECTION, which is half the mean square gradient for terrain with no
    # preferred direction. Beljaars section 2 measures the directionality of
    # real orographic spectra and finds it second order, so the halving is the
    # isotropic case rather than an assumption of convenience.
    grad2_region = np.tan(np.radians(mesh.local_slope_deg.astype(np.float64))) ** 2
    theta2_cell, _, _ = cell_mean(cells, n, area, 0.5 * grad2_region, sel)

    # Safe stand-ins off land: these cells carry the ocean value in the end and
    # the logarithms below are undefined on a zero roughness or a zero spacing.
    z0_cover_flat = effective_z0(0.0, z_blend).reshape(-1)
    z0_cover_safe = np.where(covered_land, np.maximum(z0_cover_flat, 1e-6), 0.1)
    spacing_safe = np.where(covered_land, np.maximum(spacing_cell, 1.0), 1.0e4)

    def derived_z0(band_multiple: float, closure: float):
        """Cell roughness from WM93 Eq (33) with Beljaars Eq (6) for the hills.

        Returns `(z0, Ca, Cmd, alpha, z_m)`. `z_m` is WM93's pressure scale
        height where the hills are lower than it and the hill top where they are
        higher, their Eq (4); the peak-to-valley amplitude of the resolved band
        follows from the same slope variance and wavelength the drag does, so no
        second measurement of the relief enters.
        """
        lam = float(band_multiple) * spacing_safe
        h_m = pressure_scale_height(lam, z0_cover_safe)
        l_in = inner_layer_depth(lam, z0_cover_safe)
        h_peak_to_valley = lam * np.sqrt(2.0 * theta2_cell) / np.pi
        z_m = np.maximum(h_m, h_peak_to_valley)
        ca, cmd, alpha = orographic_form_drag(theta2_cell, z_m, h_m, l_in,
                                              z0_cover_safe, closure)
        z0 = effective_length(ca + cmd, z_m)
        # WM93 Eq (33) cannot return less than the cover roughness -- the form
        # drag is added, never subtracted -- and ECMWF's cap is the only limit
        # on the other end that this tree holds a source for.
        z0 = np.clip(z0, z0_cover_safe, ECMWF_Z0_CAP_M)
        return (np.where(covered_land, z0, 0.0).reshape(nlat, nlon),
                ca, cmd, alpha, z_m)

    # LETTAU (1969) AS AN INDEPENDENT CHECK, sharing no constant with the above.
    # One mesh region is one roughness element: its characteristic height is the
    # rise across it, `h* = rms|grad h| * d`, its silhouette the triangle
    # `s = h* d / 2` the wind sees, and its lot area `S = d^2`. That is his own
    # Colorado construction at his own scale, and `check_against_published_values`
    # holds the function to the answer he prints for it.
    lettau_height = np.sqrt(2.0 * theta2_cell) * spacing_safe
    lettau_cell = lettau_z0(lettau_height, 0.5 * lettau_height * spacing_safe,
                            spacing_safe ** 2)
    lettau_land_mean = float(
        (lettau_cell[covered_land] * cell_land_area[covered_land]).sum()
        / cell_land_area[covered_land].sum())

    closures = {"first-order": 1.0,
                "second-order": WM93_CLOSURE_SECOND_ORDER / WM93_CLOSURE_FIRST_ORDER}
    derived_field, ca_cell, cmd_cell, alpha_cell, zm_cell = derived_z0(
        args.band_multiple, closures[args.closure])
    derived_land_mean = land_mean(derived_field)
    # THE DECLARED BRACKET, both ends computed on every build whichever is
    # written: the wavelength the resolved slope is attributed to, and the
    # turbulence closure WM93 measure at two levels.
    derived_bracket = {
        "band_multiple": [land_mean(derived_z0(m, closures[args.closure])[0])
                          for m in RESOLVED_BAND_BRACKET],
        "closure": [land_mean(derived_z0(args.band_multiple, c)[0])
                    for c in (closures["second-order"], closures["first-order"])],
    }

    if args.target_mean is not None:
        arm, target = "explicit", float(args.target_mean)
    elif args.orographic_coefficient is not None:
        arm, target = "fixed-coefficient", None
    elif args.orographic_arm == "none":
        arm, target = "none", None
    elif args.orographic_arm == "derived":
        arm, target = "derived", None
    else:
        # The Earth reference is reduced at the lowest model level, which is the
        # height the model applies the field at, and Mason's blending height is
        # not transferable to it: its population is a coarse Gaussian grid's
        # cells and not mesh regions inside one cell, so the scale that sets the
        # height is hundreds of kilometres there and ten here.
        arm, target = "reference", (None if reference is None
                                    else reference["effective_land_z0_m"])

    # Solve the orographic coefficient so the land mean lands on the target.
    # 60 bisections resolve the coefficient to 1e-18 on [0, 1]; the 200 this
    # carried were free when the trial was arithmetic on the grid and are not
    # now that each one reduces the mesh.
    if args.orographic_coefficient is not None:
        k_oro = float(args.orographic_coefficient)
    elif target is None:
        k_oro = 0.0
    elif reference is None:
        raise SystemExit(
            f"the reference arm solves onto the model's own Earth boundary "
            f"dataset, which the vendored tree ships at "
            f"{', '.join(sorted(EARTH_REFERENCE_RUNGS))} only. Nothing was "
            f"written.")
    else:
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if land_mean(effective_z0(mid, z_blend)) < target:
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

    z0, ce_expectation, ce_of_mean, covered = reduce_ce(k_oro, z_blend)
    if arm == "derived":
        z0 = derived_field
    field = np.where(land_cells, z0, OCEAN_Z0_M)

    # THE INVERSION IS AN IDENTITY AND IS CHECKED AS ONE. The cell cover
    # roughness is defined as the length whose `ce` is the area mean of the
    # regions' own, so putting it back through the law must return that mean.
    # This has a right answer rather than a plausible one, and it fails if the
    # law and its inverse ever stop being each other's -- which is the way a
    # reduction written twice goes wrong.
    round_trip = neutral_ce(
        effective_z0(k_oro, z_blend).reshape(-1)[covered], z_blend)
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

    # What the unknown lowest-level air temperature is worth in the field, as a
    # bracket because it cannot be verified before a climatology exists. It
    # reaches the derived arm only through `z_ref` in the reported `ce`, since
    # the blending height and the pressure scale height are properties of the
    # land rather than of the model's vertical grid; the arms that solve a
    # coefficient carry it in the field itself.
    z0_bracket = ([derived_land_mean, derived_land_mean] if arm == "derived"
                  else [land_mean(effective_z0(k_oro, reference_height_m(t, config)))
                        for t in CE_BRACKET_K])

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
                         "cover_source": mode,
                         "forest_fraction_used": forest_fraction},
        # WHERE THE COVER CAME FROM, and it is the same place code 174 and code
        # 212 came from or it is a per-mode fraction. Null under every mode that
        # implies one cover for all land; under `modelled` it names the LPJ run,
        # the climatology whose labels matched it, BIO-11's rootable partition
        # and the equilibrium window the cover was reduced over. That edge is
        # the aerodynamic half of the vegetation feedback and it is a BACK edge
        # in loop C: config/pipeline.yaml carries it on this step's `needs`.
        "modelled_cover": vegetation_provenance,
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
            # WHERE THE DISTANCE IS, and the one quantity the two planets can
            # be compared on directly. The reference dataset carries its surface
            # and orographic parts separately, so its land-cover roughness is
            # recoverable, and this world's cover term IS its land-cover
            # roughness. Neither is derived from the other. That they are close
            # is what locates the whole of the distance between the two planets
            # in the RELIEF term, which is where the terrain-information floor
            # sits. Reported rather than gated: no threshold on it was fixed in
            # advance.
            "cover_roughness_ratio_to_earth": (
                None if reference is None
                else round(unrescaled_land_mean
                           / reference["surface_part_land_mean_m"], 4)),
            "note": "the derived land mean is what this world's own subgrid "
                    "slope gives through WM93 Eq (33) and Beljaars Eq (6); the "
                    "Earth reference beside it is the effective land roughness "
                    "of PlaSim's own boundary dataset, code 173 reduced over "
                    "its land in ce, and is a COMPARISON rather than a target. "
                    "The namelist fallback dz0land is reported inside "
                    "earth_reference and is a rounding of that reference.",
        },
        # THE DERIVED ARM, and the whole of what it rests on. No entry here is
        # fitted: the constants are the papers', the slope is measured, and the
        # bracket is the two ends of the two declared attributions.
        "derived_orographic": {
            "relation": "WM93 Eq (33): karman^2/ln(z_m/z0)^2 = Ca + "
                        "karman^2/ln(z_m/z0_cover)^2, with Beljaars Eq (6) "
                        "Ca = 2 alpha beta Cmd theta_bar^2",
            "sources": ["Wood and Mason (1993), 10.1002/qj.49711951402",
                        "Beljaars, Brown and Wood (2004), 10.1256/qj.03.73",
                        "Mason (1988), 10.1002/qj.49711448007",
                        "Lettau (1969), "
                        "10.1175/1520-0450(1969)008<0828:NOARPE>2.0.CO;2"],
            "land_mean_m": round(derived_land_mean, 5),
            "ce_bracket": bracketed(derived_land_mean),
            "written": arm == "derived",
            "shape_factor_beta": WM93_SHAPE_FACTOR,
            "closure_written": args.closure,
            "band_multiple_written": args.band_multiple,
            # BOTH ENDS OF THE DECLARED BRACKET, on every build. The wavelength
            # the mesh-resolved slope is attributed to enters only through
            # logarithms; the closure is WM93's own measurement of the same
            # pressure force at two levels of turbulence closure.
            "bracket_land_mean_m": {
                "band_multiple": {
                    "ends": list(RESOLVED_BAND_BRACKET),
                    "land_mean_m": [round(v, 5)
                                    for v in derived_bracket["band_multiple"]]},
                "closure": {
                    "ends": ["second-order", "first-order"],
                    "wm93_theta2_coefficients": [WM93_CLOSURE_SECOND_ORDER,
                                                 WM93_CLOSURE_FIRST_ORDER],
                    "land_mean_m": [round(v, 5)
                                    for v in derived_bracket["closure"]]},
            },
            "slope_variance_per_direction": {
                "median": float(np.median(theta2_cell[covered_land])),
                "p95": float(np.percentile(theta2_cell[covered_land], 95)),
                "max": float(theta2_cell[covered_land].max())},
            "form_drag_over_skin_friction": {
                "median": float(np.median((ca_cell / cmd_cell)[covered_land])),
                "p95": float(np.percentile((ca_cell / cmd_cell)[covered_land], 95)),
                "max": float((ca_cell / cmd_cell)[covered_land].max())},
            "shear_parameter_alpha_median": round(
                float(np.median(alpha_cell[covered_land])), 3),
            "pressure_scale_height_m_median": round(
                float(np.median(zm_cell[covered_land])), 1),
            "mesh_spacing_m": round(cover_scale_m, 1),
            "capped_land_area_fraction": float(
                (cell_land_area[covered_land]
                 * (derived_field.reshape(-1)[covered_land]
                    >= ECMWF_Z0_CAP_M - 1e-9)).sum()
                / cell_land_area[covered_land].sum()),
            # WHERE THE MODEL'S VERTICAL GRID CANNOT CARRY THE ANSWER. WM93's
            # effective roughness reaches hundreds of times the cover roughness
            # over steep ground -- their own Table 1 runs to 300 for periodic
            # two-dimensional ridges -- and ECMWF caps at 100 m for the same
            # reason. This model's lowest level sits at `z_ref`, so a cell whose
            # roughness is a large fraction of it has no surface layer left to
            # apply a logarithmic profile in. `fluxmod.f90` forms
            # `ln(z/z0 + 1)`, so nothing diverges; what happens instead is a
            # very large exchange coefficient over a very small area. Reported
            # with the share of land area and of land-mean `ce` those cells
            # carry, and NOT gated: no threshold on it was fixed in advance.
            "above_scale_separation": {
                "lowest_model_level_m": round(z_ref_anchor, 2),
                "land_area_fraction_z0_over_tenth_of_level": float(
                    (cell_land_area[covered_land]
                     * (derived_field.reshape(-1)[covered_land]
                        > 0.1 * z_ref_anchor)).sum()
                    / cell_land_area[covered_land].sum()),
                "land_area_fraction_z0_over_level": float(
                    (cell_land_area[covered_land]
                     * (derived_field.reshape(-1)[covered_land]
                        > z_ref_anchor)).sum()
                    / cell_land_area[covered_land].sum()),
                "share_of_land_mean_ce_from_those_cells": float(
                    (cell_land_area[covered_land]
                     * neutral_ce(derived_field.reshape(-1)[covered_land],
                                  z_ref_anchor)
                     * (derived_field.reshape(-1)[covered_land]
                        > 0.1 * z_ref_anchor)).sum()
                    / (cell_land_area[covered_land]
                       * neutral_ce(derived_field.reshape(-1)[covered_land],
                                    z_ref_anchor)).sum()),
            },
            # The scheme against numbers its own papers print. It refuses rather
            # than reporting when any of them misses.
            "published_value_checks": published,
            # THE INDEPENDENT CROSS-CHECK, and it does not agree at this scale.
            # Lettau's geometric route shares no constant with WM93's and grows
            # linearly with the size of the roughness elements at fixed slope,
            # where the WM93 form saturates through the logarithmic profile it
            # is built on. The two agree to within a factor of a few at the
            # kilometre scale the operational schemes were built for and part
            # company at the ten-kilometre scale this world's terrain floor
            # puts the resolved band at, which is Beljaars' own reason for
            # assigning scales above 5 km to gravity-wave and blocking schemes
            # rather than to a roughness. Reported and NOT gated: no threshold
            # on it was fixed in advance.
            "lettau_cross_check": {
                "relation": "Lettau (1969) Eq (1), z0 = 0.5 h* s/S, one mesh "
                            "region as one roughness element",
                "orographic_land_mean_m": round(lettau_land_mean, 5),
                "ratio_to_wm93_enhancement": (
                    None if derived_land_mean <= unrescaled_land_mean
                    else round(lettau_land_mean
                               / (derived_land_mean - unrescaled_land_mean), 3)),
            },
            # WHAT THE DERIVATION CANNOT REACH, stated rather than bracketed. Ca
            # is quadratic in slope and slope variance is dominated by the
            # shortest wavelength present, so the band below Orogen's terrain-
            # information floor is where an Earth-like world's orographic
            # roughness mostly comes from. Beljaars reaches it by extrapolating
            # a power-law orographic spectrum fitted to United States
            # topography; supplying that here would import Earth's small-scale
            # terrain under a derivation's name.
            "terrain_information_floor_note":
                "the derived value is what this world's terrain supports. "
                "Orogen's measured terrain-information floor is near 20 km and "
                "Beljaars' turbulent orographic form drag lives below 5 km, so "
                "the band that dominates slope variance on Earth is absent "
                "here rather than unmeasured.",
        },
        "orographic": {
            "coefficient_solved": round(k_oro, 6),
            "relation": "z0_oro = coefficient * stdev(elevation) within the "
                        "cell. SUPERSEDED by derived_orographic: the drag is "
                        "quadratic in slope, not linear in relief amplitude, so "
                        "this coefficient was carrying the missing horizontal "
                        "scale. Retained for the reference and fixed-coefficient "
                        "arms, which exist to price the change.",
            "subgrid_stdev_m": {
                "min": round(float(sigma_m[land_cells].min()), 2),
                "median": round(float(np.median(sigma_m[land_cells])), 2),
                "max": round(float(sigma_m[land_cells].max()), 2)},
            "mesh_regions_per_land_cell_median": int(np.median(counts[land_cells])),
            # SPAT-8's constraint applies to THESE arms and not to the derived
            # one. The solve is per grid, so two rungs differ by their terrain
            # and by their calibration at once and a ladder comparison must pass
            # one coefficient to every rung. The derived orographic term is a
            # property of the land and the mesh rather than of the grid, so it
            # needs no such discipline.
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
            # MASON (1988) FIXES THE HEIGHT THE AVERAGE IS TAKEN AT, and it is
            # not the model's lowest level. His Eq (14) from the scale over
            # which the cover varies, solved on the roughness it averages.
            "blending_height_m": round(z_blend, 3),
            "blending_height_source": "Mason (1988) Eq (14), l_b ln(l_b/z0)^2 = "
                                      "2 karman^2 L_c with L_c = mesh spacing / "
                                      "2pi. The average of drag coefficients "
                                      "that reproduces the correct area-mean "
                                      "surface stress is taken at it, not at "
                                      "the lowest model level.",
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
    print(f"  derived from this world's own subgrid slope: "
          f"{derived_land_mean:.4f} m, bracket "
          f"{min(derived_bracket['band_multiple'] + derived_bracket['closure']):.4f}"
          f" to "
          f"{max(derived_bracket['band_multiple'] + derived_bracket['closure']):.4f}"
          f" m over the band and closure ends")
    print(f"  the cover term alone is {unrescaled_land_mean:.4f} m; Mason "
          f"blending height {z_blend:.1f} m against a lowest model level at "
          f"{z_ref_anchor:.1f} m")
    if reference is not None:
        print(f"  Earth comparison from {reference['dataset']}: "
              f"{reference['effective_land_z0_m']:.4f} m in ce, "
              f"{reference['land_mean_of_length_m']:.4f} m in the length, "
              f"surface part {reference['surface_part_land_mean_m']:.4f} m "
              f"-- the derived land mean is "
              f"{derived_land_mean / reference['effective_land_z0_m']:.3f} times it")
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
