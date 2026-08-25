#!/usr/bin/env python3
"""Verify the prescribed ocean heat-transport channel end to end, on a field whose answer is known.

WORLDBUILDING CONTEXT, stated first because this file's vocabulary invites
misreading: Vesper is a FICTIONAL planet and this is engineering work on the
simulation of it. Ocean, sea ice and heat transport below name modelled
quantities of an invented world. Nothing here is Earth science.

WHAT THE CHANNEL IS. `oceanmod_nl` carries `nfluko`, and at `nfluko = 1`
`oceanmod.f90:oceanini` reads surface code 903 as a fourteen-month W/m2 field
into `yfsst`, `getflxco` interpolates it to the model timestep as `yfsst2`,
`addfc` applies it to the slab under its own header comment "add oceanic flux
correction (prescribed advection)", `oceanstep` accumulates it into `yfssta`,
and `oceanout` writes that accumulation back out as code 903 of the
`ocean_output` stream. `exoplasim/scripts/close_ocean_energy.py` reads the same
code off the other end. Both ends existed and both were null: the field was
never written and the diagnostic was recorded as identically zero.

`notes/audits/ocean-and-marine-biosphere.md` section 7 makes this channel the
one the adopted offline ocean architecture depends on, because the
atmosphere-ocean transport partition is iterated rather than solved and code
903 is the only path by which a correction returns to the next atmospheric run.
So the channel is verified BEFORE anything computes a real field for it, and it
is verified with a field whose transported answer can be stated in advance --
otherwise the first real field would arrive with no way to separate a wrong
field from a wrong channel.

THE FIELD, AND WHY ITS ANSWER IS KNOWN.

    f(phi, lambda) = A * L(mu) * G(lambda),  mu = sin(phi)
    L(mu) = P2(mu) + 0.3 * mu = 0.5 * (3 mu^2 - 1) + 0.3 mu
    G(lambda) = 1 + 0.4 * cos(lambda + pi/3)

Four properties, each of which makes one prediction that the channel could fail:

1. IT IS CONSTANT IN TIME. Only one record is written to the `.sra`, and
   `surfmod.f90:get_surf_array` expands a single record to all fourteen months
   by copy. `getflxco`'s month weights `1 - zgw2` and `zgw2` sum to one whatever
   the calendar says, so `yfsst2` equals the prescribed field at EVERY ocean
   step, and `yfssta`, which is their mean over `nout` steps, equals it exactly.
   The reported field must therefore reproduce the written field cell for cell,
   with no dependence on the calendar. That is deliberate: this verifies the
   channel, not the monthly interpolation, and separating the two is the point.

2. ITS AREA-WEIGHTED GLOBAL INTEGRAL IS ZERO ANALYTICALLY. `L` is a
   second-degree polynomial in `mu` orthogonal to the constant, and the grid's
   latitudes and weights are Gauss-Legendre, which integrates a degree-2
   polynomial exactly, so `sum_j w_j L(mu_j) = 0` to roundoff. `G` has
   longitudinal mean one. A prescribed advection that does not integrate to zero
   is a heat source, so this is the conservation law the channel must not break.

3. IT IS NEITHER SYMMETRIC IN LATITUDE NOR IN LONGITUDE. The `0.3 mu` term
   breaks the equatorial symmetry `P2` alone would have, and the phase in `G`
   breaks the longitudinal one. A north-south flip, a longitudinal roll and a
   transpose all produce a DIFFERENT field, so the cell-by-cell comparison names
   which of them happened rather than reporting a difference.

4. IT IS NONZERO OVER LAND. `getflxco` interpolates over the whole grid and only
   `addfc` restricts to `yls < 1`, so the accumulated diagnostic carries the
   prescribed field over land as well as over ocean. That is stated here as a
   prediction because the code could have masked it, and because the ocean-only
   integral of a globally zero-integral field is NOT zero -- which is what a
   real forcing artifact will have to reckon with, since the ocean is what the
   correction is applied to.

THE CRITERIA ARE FIXED HERE, IN THE SOURCE, BEFORE ANY RUN IS LOOKED AT. They
are the module constants below. The reproduction bar is four orders of magnitude
above the float32 precision of the stream and four orders below the smallest
mistake that could be called one.

WHAT IS DELIVERED, NOT JUST WHAT IS REPORTED. Reproducing the diagnostic proves
the field reached `yfssta`; it does not prove `addfc` put the energy into the
slab. On an ocean cell that carries no sea ice, `mksst` and `addfc` between them
give exactly

    CRHOS * CPS * mld * (SST_k - SST_{k-1}) = (yheata_k + yfssta_k) * dt_record

with CRHOS and CPS taken through `lib/sea_water.py` from `icemod.f90`'s
declaration and the run's own `icemod_namelist` rather than written down here -- they are namelist keys, and a
copy of them is how the recovered interval can miss by five per cent while every
other criterion passes.

so across every ice-free cell and every record the temperature change must be
ONE constant times the total flux, and that constant must be
`dt_record / (CRHOS * CPS * mld)`. The constant is fitted rather than assumed,
which makes it two tests: the scatter about a single proportionality can fail,
and the fitted value can disagree with the declared mixed-layer depth and record
interval. This is `close_ocean_energy.py`'s D4 with the flux-correction term
added, computed from the ocean stream alone so it needs no cross-stream
alignment and no postprocessed output.

FOUR ARMS, AND THE FIRST ONE ISOLATES THE OTHER THREE. The `channel` arm above
was built to hold the calendar and the sea ice OUT of the answer, and that was
the reason it left two paths untested rather than an argument that they should
stay that way: a constant field makes `getflxco`'s month weights irrelevant, so
a failure cannot be attributed to the calendar, and an ice-free ocean keeps
`addfc` in branch (a), so a failure cannot be attributed to the ice. With the
channel verified, each of those can now be turned on ALONE.

    monthly   a field that varies month to month, on an ice-free ocean. Exercises
              getflxco's interpolation and nothing else.
    ice       a cold declared cold-start profile so the model grows sea ice,
              carrying the `channel` arm's constant field. Exercises addfc's
              branch (c) -- model ice, no climatological ice -- and nothing else.
    ice_clim  the same with a climatological ice thickness at code 211 staged
              beside it, so `ycliced > 0` and addfc takes branch (b) instead.

THE MONTHLY FIELD IS SEPARABLE, and that is what makes its answer known. It is
the `channel` arm's spatial pattern times a declared per-month SCALAR:

    f_m(phi, lambda) = A * L(mu) * G(lambda) * s_m

Twelve records are written and `surfmod.f90:get_surf_array` expands them
cyclically to the fourteen `yfsst` holds, so index 0 is December, 1 to 12 are
the months and 13 is January. Four predictions follow, and each is a way the
interpolation could fail:

1. SEPARABILITY. Whatever weights `momint` supplies, `yfsst2` is
   `(1-w) f_{m1} + w f_{m2}`, and both terms are the SAME spatial pattern, so
   the reported record must be that pattern times ONE scalar at every cell. A
   per-cell mixing error, a month/space transposition and a partially updated
   field all break this and none of them breaks the envelope.
2. THE ENVELOPE. `w` lies in [0, 1] whatever the calendar says, so the
   interpolated value is a convex combination of two months and must lie inside
   the envelope of the fourteen at every step and every cell. This is the
   identity the calendar cannot influence, and it is the weaker of the two.
3. THE SERIES ACTUALLY VARIES. A model that ignored the month index entirely and
   used one record forever would pass both statements above. So the recovered
   scalar series must span at least a declared fraction of the declared one, and
   the arm fails if it does not: an interpolation that was never exercised is
   not an interpolation that was verified.
4. THE GLOBAL INTEGRAL IS STILL ZERO, at every record and not only in the mean.
   Every month is a scale of one pattern orthogonal to the constant, so a
   prescribed advection carrying it adds no net energy at any point in the
   cycle.

A stronger statement -- that the series is piecewise linear with breakpoints
where `momint` changes month -- needs the calendar, so the recovered series and
its extremes are REPORTED and no bar is put on them here. That is a check on the
calendar port and belongs to whoever owns it.

THE ICE ARMS CLOSE AN ENERGY IDENTITY, and it is the `channel` arm's delivery
identity with the term that arm required to be zero put back. `mksst` withholds
the atmospheric flux from an iced slab, `mkiflux` charges the freezing clamp,
`addfc` applies part of the correction and hands the rest to `yifluxr`, and
`mkiflx` adds `yifluxr` into `yiflux`. Every one of those is a flux either
applied to the slab or charged to `yiflux`, with `yiflux` carrying the NEGATIVE
of what was withheld, so across any ocean cell and any record

    CRHOS * CPS * mld * (SST_k - SST_{k-1}) = (yheata_k + yfssta_k + yifluxa_k) * dt_record

"what addfc withholds from the slab equals what mkiflx adds to yiflux, cell by
cell", written in the three quantities the ocean stream actually carries --
`yifluxr` is not one of them. The three signs are DECLARED here before any run;
if the identity fails, the residual under each of the four sign assignments of
the `yiflux` term is reported so the failure names itself, but the pass is on
the declared one and on nothing else.

Two guards keep an ice arm from passing vacuously. The modelled ice has to
actually appear on a declared fraction of the ocean, and for `ice_clim` it has
to overlap the staged climatological ice; and `yiflux` has to be materially
nonzero, because an identity evaluated where the new term is zero is the arm
that was already run.

USAGE, in two phases, because the model runs between them:

    python exoplasim/scripts/verify_ocean_flux_channel.py stage RUNDIR --rung T42
    ...run the model in RUNDIR with nfluko = 1...
    python exoplasim/scripts/verify_ocean_flux_channel.py check RUNDIR

    python exoplasim/scripts/verify_ocean_flux_channel.py stage RUNDIR --rung T42 --arm monthly
    python exoplasim/scripts/verify_ocean_flux_channel.py stage RUNDIR --rung T42 --arm ice
    python exoplasim/scripts/verify_ocean_flux_channel.py stage RUNDIR --rung T42 --arm ice_clim

`stage` writes `N<nlat>_surf_0903.sra`, a companion sea surface temperature
climatology at code 169 that `oceanini` refuses `nfluko` without, for `ice_clim`
a sea-ice thickness climatology at code 211, and
`ocean_flux_channel_prediction.json`, which carries the arm, the field's digest,
the grid, the predicted integrals, the namelist the run must set and the
criteria. `check` reads that file and the `ocean_output` stream and answers each
criterion of the arm it names. It exits non-zero if any criterion fails.

THE INSTRUMENT IS CHECKED BEFORE THE MODEL IS RUN, and that is what `selftest`
is for:

    python exoplasim/scripts/verify_ocean_flux_channel.py selftest RUNDIR

Every criterion here is a statement about a stream, so a stream satisfying all of
them exactly is synthesised for the staged arm and each criterion is then
perturbed in one named way. A criterion that passes the clean stream and fails
its own break has been shown to have a right answer. This verifies the
INSTRUMENT and says nothing about the model; it is what stands between staging an
arm and running it.

IT ALSO CAUGHT A BAR THAT COULD NOT DISCRIMINATE. `TOL_DELIVERY_RELATIVE` is a
fraction of the largest temperature change in the window, and `ysst` is stored
float32, so on a window whose changes are hundredths of a kelvin the bar sits
BELOW the storage resolution and the criterion reports the format rather than the
model. `MIN_DELIVERY_CHANGE_K` is the floor that makes it answerable, it applies
to every arm including `channel`, and the segment section 11b already ran clears
it by a factor of five.

EACH ARM IS ITS OWN RUN DIRECTORY. The arms differ in the namelist as well as in
the staged fields, and the prediction file is written under one name, so staging
a second arm over a first one replaces it.

The finding this produces belongs in
`notes/audits/ocean-and-marine-biosphere.md`; this script is the instrument.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import struct
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from _paths import PROJECT_ROOT  # noqa: F401  (also puts lib/ on sys.path)

import rungs  # noqa: E402  from lib/, via _paths
import sea_water  # noqa: E402  from lib/, via _paths
from sra import read_sra, write_sra  # noqa: E402  same directory

# THE SEA WATER PAIR COMES FROM THE MODEL, never from a literal here.
# `icemod.f90` declares density and specific heat as `icemod_nl` keys and passes
# them to `oceanini`, so a run can set them and two runs need not agree.
# `lib/sea_water.py` reads the declaration and then the run's own namelist. This verification's first pass carried the pair as
# literals copied from that file, which had them at fresh water's 4180 J/kg/K
# while the model had moved to sea water's 3990.34, and the fitted record
# interval missed the declared one by 4.75 per cent for that reason alone.

# ---------------------------------------------------------------------------
# THE FIELD AND THE CRITERIA. Declared before any run exists, and not to be
# adjusted after one does: a criterion chosen after the run it judges is not a
# criterion. CLAUDE.md, `docs/src/practice/conventions.md`.
# ---------------------------------------------------------------------------

# W/m2. The magnitude of an ocean heat-transport convergence, chosen so the
# signal is far above the float32 resolution of the output stream and far above
# the five-decimal write of the `.sra`, and so that a sign error or a factor of
# two is unmistakable.
AMPLITUDE_W_M2 = 20.0
MERIDIONAL_ASYMMETRY = 0.3   # weight of the P1 term that breaks equatorial symmetry
ZONAL_AMPLITUDE = 0.4        # weight of the zonal wave that breaks zonal symmetry
ZONAL_PHASE_RAD = np.pi / 3  # so a roll of +k and -k are different fields

# The sea surface temperature climatology staged beside it. `oceanini` aborts
# under `nfluko` without one, because the correction is a relaxation toward a
# climatology. Uniform and well above any freezing point, so the ocean starts
# ice-free everywhere and the delivery identity is testable on the whole ocean.
CLIMATOLOGICAL_SST_K = 290.0

# Criteria.
#
# TOL_REPRODUCTION: the stream stores float32, which resolves 20 W/m2 to about
# 1.2e-6, and the `.sra` is written to five decimals, which the comparison
# absorbs by comparing against the round-tripped field rather than the ideal
# one. 1e-4 W/m2 is two orders above that noise floor and five orders below the
# smallest structural error -- a sign, a roll, a mask, a unit -- that could
# occur.
TOL_REPRODUCTION_W_M2 = 1.0e-4
# TOL_STAGE_INTEGRAL_IDEAL: the field BEFORE it is written, where the integral
# is exact arithmetic on a Gauss-Legendre quadrature and lands at roundoff.
TOL_STAGE_INTEGRAL_IDEAL_W_M2 = 1.0e-9
# TOL_STAGE_INTEGRAL_WRITTEN: the field AS WRITTEN, which is a different object
# and needs a bar derived from the instrument rather than from the mathematics.
# `sra.py` writes five decimals, so each cell carries a quantisation error
# uniform on +-5e-6 W/m2, standard deviation 2.9e-6. Averaged over the nlat*nlon
# cells of a rung -- 8192 at T42, more above it -- the area-weighted mean of
# those errors has a standard deviation near 3e-8 W/m2. 1e-6 is thirty times
# that and is still five parts in 1e8 of the field's own amplitude, so it says
# "zero to the precision the file format can carry" and nothing weaker.
TOL_STAGE_INTEGRAL_WRITTEN_W_M2 = 1.0e-6
# TOL_OUTPUT_INTEGRAL: the same integral taken on each output record, so it
# carries the float32 storage of every cell.
TOL_OUTPUT_INTEGRAL_W_M2 = 1.0e-4
# TOL_DELIVERY: the residual about the single fitted proportionality, as a
# fraction of the largest temperature change in the window. The identity is
# exact in the model's own arithmetic; what limits it is that `ysst` is stored
# float32 and the differences are small, so the bar is loose enough to survive
# that and tight enough that a missing or doubled term cannot pass.
TOL_DELIVERY_RELATIVE = 1.0e-3
# TOL_RECORD_SECONDS: the fitted constant against the declared mixed-layer
# depth and record interval.
TOL_RECORD_SECONDS_RELATIVE = 1.0e-3
# MIN_DELIVERY_CHANGE_K: the smallest largest-temperature-change a window can
# carry and still be able to answer the delivery question at all.
#
# THE INSTRUMENT HAS TO BE CHECKED AGAINST THE SIZE OF THE EFFECT. `ysst` is
# stored float32, so near 290 K it resolves to about 1.7e-5 K and the difference
# of two such values carries about twice that. TOL_DELIVERY_RELATIVE is a
# fraction of the LARGEST change in the window, so the bar in kelvin is that
# fraction times the largest change: at 0.469 K, which is what a cold-start
# segment gives, the bar is 4.7e-4 K and sits thirty times above the
# quantisation; at 0.02 K it is 2e-5 K and sits BELOW it, and the criterion then
# reports the storage format rather than the model. 0.1 K puts the bar at 1e-4 K
# and three times the quantisation, which is the smallest window this criterion
# can discriminate on. Below it the answer is "not answerable", and that is a
# failure of the segment rather than a pass. CLAUDE.md and
# docs/src/practice/failure-modes.md class 34.
MIN_DELIVERY_CHANGE_K = 0.1

# ---------------------------------------------------------------------------
# THE ARMS, AND THE CRITERIA EACH ADDS. Declared here, before any run of any of
# them exists, for the reason the block above gives.
# ---------------------------------------------------------------------------

ARMS = ("channel", "monthly", "ice", "ice_clim")

# --- monthly ---------------------------------------------------------------
# The per-month scalar the spatial pattern is multiplied by, months 1 to 12.
# A cosine and a sine of the same period, so the sequence is asymmetric under a
# month reversal and a reversed calendar is a different series rather than the
# same one read backwards. The mean is one, so the pattern's own amplitude is
# what the orbit mean of the correction carries.
MONTHLY_SCALES = tuple(
    1.0 + 0.5 * np.cos(2.0 * np.pi * (m - 1) / 12.0)
        + 0.2 * np.sin(2.0 * np.pi * (m - 1) / 12.0)
    for m in range(1, 13))
# TOL_SEPARABILITY_RELATIVE: the reported record against one scalar times the
# pattern, as a fraction of the record's own largest value. The stream stores
# float32, which resolves the field to about one part in 1e7, and the `.sra`'s
# five-decimal write is absorbed by comparing against the round-tripped pattern.
# 1e-5 is two orders above that and four below any structural error.
TOL_SEPARABILITY_RELATIVE = 1.0e-5
# MONTHLY_MIN_OBSERVED_SPAN_FRACTION: how much of the declared scalar span the
# recovered series has to cover before the arm counts as having exercised the
# interpolation at all. Half, which a segment covering a full orbit reaches
# comfortably and a segment that ignored the month index reaches never.
MONTHLY_MIN_OBSERVED_SPAN_FRACTION = 0.5

# --- ice and ice_clim ------------------------------------------------------
# The declared cold sea surface temperature climatology, as the same
# hemispherically symmetric form icemod.f90's own cold start uses:
# T_pol + (T_eq - T_pol) * cos(lat)**2. The pole value is BELOW any freezing
# point the salinity bracket reaches, so the model grows ice at high latitude,
# and the equator value is well above it, so ice-free ocean and the branch (a)
# control survive in the same segment as the ice branches.
COLD_SST_EQUATOR_K = 300.0
COLD_SST_POLE_K = 265.0
# The climatological sea-ice thickness staged at code 211 for the `ice_clim`
# arm, on the cells whose declared climatological temperature is at or below the
# freezing point. Thick enough that ycliced > 0 is unambiguous in float32.
CLIMATOLOGICAL_ICE_THICKNESS_M = 2.0
# MIN_ICED_OCEAN_FRACTION: how much of the ocean has to carry modelled ice in at
# least one record before an ice arm counts as having reached its branch. Two
# per cent of the ocean at the cold end of a 35 K profile is a low bar to clear
# and an unmissable one to fail.
MIN_ICED_OCEAN_FRACTION = 0.02
# MIN_ICE_FREE_OCEAN_FRACTION: and how much has to stay ice-free in every
# record, so the branch (a) control is present in the same run and a delivery
# failure can be localised to the ice branches rather than to the slab.
MIN_ICE_FREE_OCEAN_FRACTION = 0.10
# MIN_ICE_FLUX_W_M2: the ice term has to be materially nonzero somewhere, or the
# identity below degenerates to the one the `channel` arm already closed. A
# tenth of the correction amplitude.
MIN_ICE_FLUX_W_M2 = 0.1 * AMPLITUDE_W_M2
# The signs of the three-term slab budget, declared before any run. See the
# module docstring: yiflux carries the NEGATIVE of what was withheld from the
# slab, so it enters with the same sign as the two fluxes offered.
SLAB_BUDGET_SIGNS = {"yheat": 1.0, "yfsst": 1.0, "yiflux": 1.0}
# TOL_FREEZING_CLAMP_K: addfc's final clamp and mkiflux both hold the modelled
# sea surface at or above the freezing point, so a reported value below it is a
# clamp that did not fire. The bar is float32 storage of a temperature near
# 271 K, which resolves to about 3e-5 K; 1e-3 is thirty times that and far below
# any missed clamp, which would be tenths of a kelvin at least.
TOL_FREEZING_CLAMP_K = 1.0e-3

OCEAN_CODES = {901: "yheat", 902: "yiflux", 903: "yfsst", 904: "ydsst",
               905: "yqhd", 906: "yfldo", 910: "yicec", 939: "ysst",
               972: "yls", 990: "yclsst"}
# Terms this verification requires to be off, so that the slab integration it
# checks has exactly two sources. Each is a namelist switch or a compiled
# dimension, and a nonzero value means the run is not the one described.
MUST_BE_ZERO = {"ydsst": "ocean vertical diffusion, NLEV_OCE = 1",
                "yqhd": "ocean horizontal diffusion, nhdiff = 0",
                "yfldo": "deep-ocean flux, nlsg = 0"}

PREDICTION_NAME = "ocean_flux_channel_prediction.json"


def grid(nlat: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gaussian latitudes north to south, their weights, and the longitudes.

    `lib/gridding.py` is the one grid convention and it runs rows north to
    south; `leggauss` returns nodes ascending in sin(latitude), which is south
    to north, so both are reversed together. Reversing one and not the other is
    the error `predict_ocean_terms.py` records having already made once.
    """
    nlon = 2 * nlat
    mu, w = leggauss(nlat)
    mu, w = mu[::-1], w[::-1]
    lon = 2.0 * np.pi * np.arange(nlon) / nlon
    return mu, w, lon


def prescribed_field(nlat: int) -> np.ndarray:
    """The field, as (nlat, nlon), north to south. See the module docstring."""
    mu, _, lon = grid(nlat)
    lat_part = 0.5 * (3.0 * mu ** 2 - 1.0) + MERIDIONAL_ASYMMETRY * mu
    lon_part = 1.0 + ZONAL_AMPLITUDE * np.cos(lon + ZONAL_PHASE_RAD)
    return AMPLITUDE_W_M2 * np.outer(lat_part, lon_part)


def write_sra_months(path: Path, code: int, fields: np.ndarray) -> None:
    """Write a twelve-month `.sra`: one `write_sra` record per month, in order.

    `sra.py` writes one record and this needs twelve, so the records are
    appended rather than a second writer being introduced -- the header, the
    fixed date stamp and the five-decimal field format all stay `sra.py`'s.
    `surfmod.f90:get_surf_array` expands twelve records cyclically into the
    fourteen `yfsst` holds, which is the model's own path and is exercised here
    rather than bypassed by writing fourteen.
    """
    if fields.shape[0] != 12:
        raise ValueError(f"{fields.shape[0]} months, expected 12")
    parts = []
    for month in range(12):
        write_sra(path, code, fields[month])
        parts.append(path.read_text(encoding="ascii"))
    path.write_text("".join(parts), encoding="ascii")


def read_sra_months(path: Path, nlat: int, nlon: int) -> np.ndarray:
    """Read a twelve-month `.sra` back, as (12, nlat, nlon)."""
    lines = path.read_text(encoding="ascii").splitlines()
    per_record = 1 + (nlat * nlon) // 8
    if len(lines) % per_record:
        raise ValueError(f"{path}: {len(lines)} lines is not a whole number of "
                         f"{per_record}-line records for {nlat}x{nlon}")
    months = len(lines) // per_record
    out = np.empty((months, nlat, nlon))
    for m in range(months):
        block = lines[m * per_record + 1:(m + 1) * per_record]
        out[m] = np.array(" ".join(block).split(),
                          dtype=np.float64).reshape(nlat, nlon)
    return out


def expand_to_fourteen(months: np.ndarray) -> np.ndarray:
    """The twelve written months as the fourteen `yfsst` holds.

    `get_surf_array`'s own cyclic expansion, reproduced so the envelope is taken
    over what the model has rather than over what was written: index 0 is
    December, 1 to 12 are the months and 13 is January.
    """
    return np.concatenate([months[11:12], months, months[0:1]], axis=0)


def latitude_field(nlat: int, nlon: int, pole: float, equator: float) -> np.ndarray:
    """`pole + (equator - pole) * cos(lat)**2`, as (nlat, nlon), north to south.

    icemod.f90's own declared cold-start form, reproduced so the staged sea
    surface temperature climatology and the profile the model would build from
    icemod_nl are the same shape. Hemispherically symmetric by construction, for
    the reason config/planet.yaml gives: an asymmetric initial condition puts a
    difference into the answer that nothing in the world put there.
    """
    mu, _, _ = grid(nlat)
    cos2 = 1.0 - mu ** 2
    return np.repeat((pole + (equator - pole) * cos2)[:, None], nlon, axis=1)


def global_integral(field: np.ndarray, w: np.ndarray) -> float:
    """Area-weighted global mean, W/m2. The Gaussian weights sum to two."""
    return float((field.mean(axis=1) * w).sum() / w.sum())


def masked_mean(field: np.ndarray, w: np.ndarray, mask: np.ndarray) -> float:
    """Area-weighted mean over a mask, per unit area of the mask."""
    weights = np.broadcast_to(w[:, None], field.shape)
    area = float((weights * mask).sum())
    if area <= 0.0:
        return float("nan")
    return float((field * weights * mask).sum()) / area


def digest(field: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(field, dtype=np.float64)
                          .tobytes()).hexdigest()


def namelist_value(path: Path, key: str) -> float | None:
    """One real from a Fortran namelist file, or None if it is not set there."""
    if not path.is_file():
        return None
    m = re.search(rf"^\s*{key}\s*=\s*([-+0-9.eEdD]+)", path.read_text(),
                  re.IGNORECASE | re.MULTILINE)
    return float(m.group(1).replace("d", "e").replace("D", "e")) if m else None


def read_stream(paths: list[Path], keep: dict[int, str]) -> tuple[dict, list]:
    """Records of a PlaSim service stream, plus each record's eight-word header.

    Same layout `close_ocean_energy.py:read_service` reads: one Fortran
    unformatted record of eight int32 header words, then one of NLON*NLAT
    float32 values. The headers are kept here because words three and four are
    the model's own date and time and are how the record interval is checked
    against the namelist rather than assumed from it.
    """
    out: dict[int, list[np.ndarray]] = {}
    heads: dict[int, list[np.ndarray]] = {}
    shape = None
    for path in paths:
        size = os.path.getsize(path)
        with open(path, "rb") as handle:
            while handle.tell() < size:
                length = struct.unpack("i", handle.read(4))[0]
                head = np.frombuffer(handle.read(length), dtype=np.int32).copy()
                handle.read(4)
                length = struct.unpack("i", handle.read(4))[0]
                raw = handle.read(length)
                handle.read(4)
                nlon, nlat = int(head[4]), int(head[5])
                if shape is None:
                    shape = (nlat, nlon)
                elif shape != (nlat, nlon):
                    raise SystemExit(f"{path} changes grid mid-file")
                code = int(head[0])
                if code in keep:
                    out.setdefault(code, []).append(
                        np.frombuffer(raw, dtype=np.float32).reshape(shape))
                    heads.setdefault(code, []).append(head)
    missing = [name for code, name in keep.items() if code not in out]
    if missing:
        raise SystemExit(f"{[str(p) for p in paths]} does not carry {missing}. "
                         "The ocean stream is written whenever nocean and "
                         "noutput are on; code 903 needs nfluko = 1 as well.")
    fields = {keep[code]: np.stack(v).astype(float) for code, v in out.items()}
    return fields, heads[903]


def stream_paths(run_dir: Path) -> list[Path]:
    """The ocean stream, per-call files first, then the in-place name."""
    per_call = sorted(run_dir.glob("MOST_OCEAN.[0-9][0-9][0-9][0-9][0-9]"))
    if per_call:
        return per_call
    single = run_dir / "ocean_output"
    if single.is_file():
        return [single]
    raise SystemExit(f"{run_dir} carries neither MOST_OCEAN.NNNNN nor "
                     "ocean_output, so no ocean stream was written")


# ---------------------------------------------------------------------------
# stage
# ---------------------------------------------------------------------------

def stage(run_dir: Path, rung: str, arm: str = "channel") -> dict:
    if arm != "channel":
        return stage_extra(run_dir, rung, arm)
    nlat = rungs.RUNGS[rung]
    nlon = 2 * nlat
    _, w, _ = grid(nlat)
    ideal = prescribed_field(nlat)

    run_dir.mkdir(parents=True, exist_ok=True)
    flux_path = run_dir / f"N{nlat:03d}_surf_0903.sra"
    sst_path = run_dir / f"N{nlat:03d}_surf_0169.sra"
    write_sra(flux_path, 903, ideal)
    write_sra(sst_path, 169, np.full((nlat, nlon), CLIMATOLOGICAL_SST_K))

    # THE PREDICTION IS THE FIELD AS WRITTEN, not as intended. The `.sra` holds
    # five decimals, so the model reads the rounded field and can only be asked
    # to reproduce that. Reading it back is also the only check that the writer
    # and the reader agree about the grid.
    written = read_sra(flux_path, nlat, nlon)
    if written.shape != (nlat, nlon):
        raise SystemExit(f"{flux_path} read back as {written.shape}, not "
                         f"{(nlat, nlon)}")
    ideal_integral = global_integral(ideal, w)
    if abs(ideal_integral) > TOL_STAGE_INTEGRAL_IDEAL_W_M2:
        raise SystemExit(
            f"the field's area-weighted global integral is {ideal_integral:.3e} "
            f"W/m2 before it is written, against a bar of "
            f"{TOL_STAGE_INTEGRAL_IDEAL_W_M2:.0e}. It is zero by construction "
            "on a Gauss-Legendre quadrature, so either the quadrature or the "
            "field is not what this file says it is.")
    integral = global_integral(written, w)
    if abs(integral) > TOL_STAGE_INTEGRAL_WRITTEN_W_M2:
        raise SystemExit(
            f"the staged field's area-weighted global integral is {integral:.3e} "
            f"W/m2 once written, against a bar of "
            f"{TOL_STAGE_INTEGRAL_WRITTEN_W_M2:.0e}, which is the five-decimal "
            "write of sra.py averaged over the grid. The file format cannot "
            "lose this much of the integral, so something else did.")

    prediction = {
        "what": "prescribed ocean heat-transport channel, surface code 903 "
                "under oceanmod_nl nfluko = 1",
        "rung": rung, "nlat": nlat, "nlon": nlon,
        "field": {
            "form": "A * (P2(mu) + b*mu) * (1 + c*cos(lambda + p)), mu = sin(latitude)",
            "amplitude_w_m2": AMPLITUDE_W_M2,
            "meridional_asymmetry": MERIDIONAL_ASYMMETRY,
            "zonal_amplitude": ZONAL_AMPLITUDE,
            "zonal_phase_rad": ZONAL_PHASE_RAD,
            "months_written": 1,
            "expanded_to_months": 14,
            "sha256_as_written": digest(written),
            "min_w_m2": float(written.min()), "max_w_m2": float(written.max()),
        },
        "companion_sst_climatology_k": CLIMATOLOGICAL_SST_K,
        "files": {"flux": flux_path.name, "sst": sst_path.name},
        "predicted": {
            "global_integral_w_m2": integral,
            "global_integral_before_writing_w_m2": ideal_integral,
            "reported_field_equals_written_field": True,
            "reported_over_land_equals_written_field": True,
        },
        "criteria": {
            "reproduction_w_m2": TOL_REPRODUCTION_W_M2,
            "stage_integral_ideal_w_m2": TOL_STAGE_INTEGRAL_IDEAL_W_M2,
            "stage_integral_written_w_m2": TOL_STAGE_INTEGRAL_WRITTEN_W_M2,
            "output_integral_w_m2": TOL_OUTPUT_INTEGRAL_W_M2,
            "delivery_relative": TOL_DELIVERY_RELATIVE,
            "record_seconds_relative": TOL_RECORD_SECONDS_RELATIVE,
        },
        "requires_zero": MUST_BE_ZERO,
        "namelist_required": {"oceanmod_nl": {"nfluko": 1, "nocean": 1,
                                              "nhdiff": 0, "nlsg": 0}},
        "arm": "channel",
    }
    (run_dir / PREDICTION_NAME).write_text(
        json.dumps(prediction, indent=2) + "\n", encoding="utf-8")
    return prediction


def stage_extra(run_dir: Path, rung: str, arm: str) -> dict:
    """The three arms the `channel` arm held the calendar and the ice out of.

    Each turns on exactly one of them and leaves the other off, so a failure
    names a path rather than a combination.
    """
    nlat = rungs.RUNGS[rung]
    nlon = 2 * nlat
    _, w, _ = grid(nlat)
    run_dir.mkdir(parents=True, exist_ok=True)
    pattern = prescribed_field(nlat)
    flux_path = run_dir / f"N{nlat:03d}_surf_0903.sra"
    sst_path = run_dir / f"N{nlat:03d}_surf_0169.sra"
    files = {"flux": flux_path.name, "sst": sst_path.name}
    prediction: dict = {"arm": arm, "rung": rung, "nlat": nlat, "nlon": nlon,
                        "requires_zero": MUST_BE_ZERO}

    if arm == "monthly":
        months = np.stack([pattern * s for s in MONTHLY_SCALES])
        write_sra_months(flux_path, 903, months)
        # Ice-free and well above any freezing point, so addfc stays in branch
        # (a) and this arm carries the calendar and nothing else.
        write_sra(sst_path, 169, np.full((nlat, nlon), CLIMATOLOGICAL_SST_K))
        written = read_sra_months(flux_path, nlat, nlon)
        if written.shape != (12, nlat, nlon):
            raise SystemExit(f"{flux_path} read back as {written.shape}, not "
                             f"{(12, nlat, nlon)}")
        envelope = expand_to_fourteen(written)
        integrals = [global_integral(envelope[m], w) for m in range(14)]
        worst = float(np.abs(integrals).max())
        if worst > TOL_STAGE_INTEGRAL_WRITTEN_W_M2 * max(MONTHLY_SCALES):
            raise SystemExit(
                f"the staged field's worst monthly area-weighted global integral "
                f"is {worst:.3e} W/m2 once written. Every month is a scale of one "
                "pattern orthogonal to the constant, so the file format cannot "
                "lose this much of the integral and something else did.")
        prediction.update({
            "what": "the MONTHLY interpolation of surface code 903: getflxco's "
                    "convex combination of two months, on an ice-free ocean",
            "field": {
                "form": "A * (P2(mu) + b*mu) * (1 + c*cos(lambda + p)) * s_m",
                "amplitude_w_m2": AMPLITUDE_W_M2,
                "months_written": 12,
                "expanded_to_months": 14,
                "expansion": "surfmod.f90:get_surf_array, cyclic: 0 is December, "
                             "1 to 12 the months, 13 January",
                "monthly_scales": [float(x) for x in MONTHLY_SCALES],
                "scale_span": float(max(MONTHLY_SCALES) - min(MONTHLY_SCALES)),
                "sha256_pattern_as_written": digest(written[0] / MONTHLY_SCALES[0]),
                "sha256_months_as_written": digest(written),
            },
            "companion_sst_climatology_k": CLIMATOLOGICAL_SST_K,
            "files": files,
            "predicted": {
                "reported_is_one_scalar_times_the_pattern": True,
                "reported_lies_inside_the_fourteen_month_envelope": True,
                "global_integral_every_record_w_m2": 0.0,
                "recovered_series_spans_the_declared_one": True,
            },
            "criteria": {
                "separability_relative": TOL_SEPARABILITY_RELATIVE,
                "envelope_w_m2": TOL_REPRODUCTION_W_M2,
                "output_integral_w_m2": TOL_OUTPUT_INTEGRAL_W_M2,
                "observed_span_fraction": MONTHLY_MIN_OBSERVED_SPAN_FRACTION,
                "delivery_relative": TOL_DELIVERY_RELATIVE,
                "record_seconds_relative": TOL_RECORD_SECONDS_RELATIVE,
            },
            "namelist_required": {"oceanmod_nl": {"nfluko": 1, "nocean": 1,
                                                  "nhdiff": 0, "nlsg": 0}},
            "run_must_cover": "a full orbit, or the recovered series cannot span "
                              "the declared one and the arm fails by design",
        })
    elif arm in ("ice", "ice_clim"):
        # The `channel` arm's field, unchanged and constant in time, so this arm
        # carries the ice branches and nothing else.
        write_sra(flux_path, 903, pattern)
        written = read_sra(flux_path, nlat, nlon)
        cold = latitude_field(nlat, nlon, COLD_SST_POLE_K, COLD_SST_EQUATOR_K)
        write_sra(sst_path, 169, cold)
        water = sea_water.constants(None)
        freezing = float(water["TFREEZE"])
        below = cold <= freezing
        if not below.any():
            raise SystemExit(
                f"the declared cold profile reaches {cold.min():.2f} K and the "
                f"model's freezing point is {freezing:.2f} K, so no cell would "
                "grow ice and neither ice arm can reach its branch")
        if arm == "ice_clim":
            ice_path = run_dir / f"N{nlat:03d}_surf_0211.sra"
            thickness = np.where(below, CLIMATOLOGICAL_ICE_THICKNESS_M, 0.0)
            write_sra(ice_path, 211, thickness)
            files["ice_climatology"] = ice_path.name
        prediction.update({
            "what": ("addfc branch (b), model ice and climatological ice"
                     if arm == "ice_clim" else
                     "addfc branch (c), model ice and no climatological ice"),
            "field": {
                "form": "A * (P2(mu) + b*mu) * (1 + c*cos(lambda + p)), constant in time",
                "amplitude_w_m2": AMPLITUDE_W_M2,
                "months_written": 1,
                "expanded_to_months": 14,
                "sha256_as_written": digest(written),
            },
            "cold_start": {
                "form": "T_pol + (T_eq - T_pol) * cos(lat)**2, icemod.f90's own",
                "sst_equator_k": COLD_SST_EQUATOR_K,
                "sst_pole_k": COLD_SST_POLE_K,
                "freezing_point_k": freezing,
                "sea_water": water,
                "cells_at_or_below_freezing": int(below.sum()),
                "climatological_ice_thickness_m":
                    CLIMATOLOGICAL_ICE_THICKNESS_M if arm == "ice_clim" else 0.0,
                "note": "oceanini clamps the climatology up to the freezing "
                        "point, so yclsst2 is TFREEZE wherever this profile is "
                        "below it, and that is the temperature addfc branch (b) "
                        "drives the modelled sea surface toward",
            },
            "files": files,
            "predicted": {
                "slab_budget_signs": SLAB_BUDGET_SIGNS,
                "sea_surface_stays_at_or_above_freezing": True,
                "ice_branch_is_reached": True,
                "an_ice_free_control_survives_in_the_same_run": True,
            },
            "criteria": {
                "reproduction_w_m2": TOL_REPRODUCTION_W_M2,
                "output_integral_w_m2": TOL_OUTPUT_INTEGRAL_W_M2,
                "delivery_relative": TOL_DELIVERY_RELATIVE,
                "record_seconds_relative": TOL_RECORD_SECONDS_RELATIVE,
                "freezing_clamp_k": TOL_FREEZING_CLAMP_K,
                "iced_ocean_fraction": MIN_ICED_OCEAN_FRACTION,
                "ice_free_ocean_fraction": MIN_ICE_FREE_OCEAN_FRACTION,
                "ice_flux_w_m2": MIN_ICE_FLUX_W_M2,
            },
            "namelist_required": {
                "oceanmod_nl": {"nfluko": 1, "nocean": 1, "nhdiff": 0, "nlsg": 0},
                # icemod's OWN nfluko stays off. It is a different key in a
                # different namelist, and at nfluko = 1 icemod would relax the
                # ice toward a climatology as well, which is a second correction
                # this arm is not testing.
                "icemod_nl": {"nice": 1, "nseaice": 1, "nfluko": 0, "ntskin": 1},
            },
            "run_must_start": "COLD. A restart carries an ocean that is already "
                              "warm and the declared profile never runs",
        })
    else:
        raise SystemExit(f"unknown arm {arm!r}")

    (run_dir / PREDICTION_NAME).write_text(
        json.dumps(prediction, indent=2) + "\n", encoding="utf-8")
    return prediction


# ---------------------------------------------------------------------------
# check
# ---------------------------------------------------------------------------

def _terms_that_must_be_off(fields: dict) -> dict:
    off = {}
    for name, why in MUST_BE_ZERO.items():
        peak = float(np.abs(fields[name]).max())
        off[name] = {"peak_w_m2": peak, "requires": why, "pass": peak == 0.0}
    return off


def _fit_delivery(dT: np.ndarray, flux: np.ndarray, mask: np.ndarray,
                  run_dir: Path, mixed_layer_m: float | None,
                  record_seconds: float | None) -> dict:
    """One fitted proportionality between a temperature change and a flux.

    Shared by every arm, because every arm closes the same slab budget and only
    the flux differs: the `channel` and `monthly` arms offer two terms on cells
    that carry no ice, the ice arms offer three on every ocean cell.
    """
    m = np.broadcast_to(mask[None, :, :], dT.shape)
    x, y = flux[m], dT[m]
    slope = float((x * y).sum() / (x * x).sum()) if x.size and (x * x).sum() else float("nan")
    resid = y - slope * x
    scale = float(np.abs(y).max()) if y.size else float("nan")
    rel = float(np.abs(resid).max() / scale) if scale else float("nan")
    fitted_seconds = slope * sea_water.slab_heat_capacity(
        mixed_layer_m if mixed_layer_m else float("nan"), run_dir)
    out = {
        "cells": int(mask.sum()), "record_pairs": int(dT.shape[0]),
        "samples": int(x.size),
        "fitted_k_per_w_m2": slope,
        "max_abs_residual_k": float(np.abs(resid).max()) if resid.size else None,
        "largest_change_k": scale,
        "relative_residual": rel,
        "criterion_relative": TOL_DELIVERY_RELATIVE,
        "criterion_largest_change_k": MIN_DELIVERY_CHANGE_K,
        "discriminates": bool(np.isfinite(scale) and scale >= MIN_DELIVERY_CHANGE_K),
        "discrimination_note":
            "the relative bar is a fraction of the largest change, so below "
            "MIN_DELIVERY_CHANGE_K it sits under the float32 resolution of ysst "
            "and reports the storage format rather than the model. A window that "
            "small has not answered the delivery question",
        "pass": bool(x.size and rel <= TOL_DELIVERY_RELATIVE
                     and np.isfinite(scale) and scale >= MIN_DELIVERY_CHANGE_K),
        "mixed_layer_depth_m": mixed_layer_m,
        "sea_water": sea_water.constants(run_dir),
        "fitted_record_seconds": fitted_seconds,
        "declared_record_seconds": record_seconds,
    }
    if record_seconds and np.isfinite(fitted_seconds):
        off_by = abs(fitted_seconds / record_seconds - 1.0)
        out["record_seconds_relative_error"] = float(off_by)
        out["record_seconds_pass"] = bool(off_by <= TOL_RECORD_SECONDS_RELATIVE)
    return out


def check(run_dir: Path, record_seconds: float | None,
          mixed_layer_m: float | None) -> dict:
    prediction = json.loads((run_dir / PREDICTION_NAME).read_text(encoding="utf-8"))
    arm = prediction.get("arm", "channel")
    if arm != "channel":
        return check_extra(run_dir, prediction, record_seconds, mixed_layer_m)
    nlat, nlon = int(prediction["nlat"]), int(prediction["nlon"])
    _, w, _ = grid(nlat)
    written = read_sra(run_dir / prediction["files"]["flux"], nlat, nlon)
    if digest(written) != prediction["field"]["sha256_as_written"]:
        raise SystemExit(
            f"{prediction['files']['flux']} is not the file this prediction was "
            "written for. The channel cannot be verified against a field that "
            "changed between staging and the run.")

    fields, heads = read_stream(stream_paths(run_dir), OCEAN_CODES)
    if fields["yfsst"].shape[1:] != (nlat, nlon):
        raise SystemExit(f"the stream is on {fields['yfsst'].shape[1:]} and the "
                         f"staged field on {(nlat, nlon)}")
    nrec = fields["yfsst"].shape[0]
    results: dict[str, dict] = {}

    # Terms that must be off, so the slab has exactly two sources.
    off = {}
    for name, why in MUST_BE_ZERO.items():
        peak = float(np.abs(fields[name]).max())
        off[name] = {"peak_w_m2": peak, "requires": why, "pass": peak == 0.0}
    results["terms_that_must_be_off"] = off

    # C1 REPRODUCTION, cell by cell and record by record.
    err = np.abs(fields["yfsst"] - written[None, :, :])
    worst = float(err.max())
    where = np.unravel_index(int(err.argmax()), err.shape)
    results["C1_reproduction"] = {
        "statement": "every record of code 903 equals the written field, cell "
                     "for cell, because the field is constant in time and the "
                     "month weights sum to one",
        "max_abs_error_w_m2": worst,
        "at_record_lat_lon": [int(x) for x in where],
        "criterion_w_m2": TOL_REPRODUCTION_W_M2,
        "pass": worst <= TOL_REPRODUCTION_W_M2,
    }

    # C2 ORIENTATION. Diagnostic rather than a second bar: it names the failure
    # C1 would report as a difference.
    reported = fields["yfsst"][0]
    alternatives = {
        "as_written": written,
        "north_south_flipped": written[::-1, :],
        "sign_reversed": -written,
        "transposed_if_square": (written.T if written.shape[0] == written.shape[1]
                                 else None),
    }
    for k in (1, nlon // 4, nlon // 2):
        alternatives[f"rolled_east_{k}"] = np.roll(written, k, axis=1)
    results["C2_orientation"] = {
        "statement": "the field is asymmetric in both latitude and longitude, "
                     "so a flip, a roll, a sign error and a transpose are "
                     "different fields and the nearest one names the failure",
        "max_abs_difference_w_m2": {
            name: (float(np.abs(reported - alt).max()) if alt is not None
                   else None)
            for name, alt in alternatives.items()},
    }

    # C3 GLOBAL INTEGRAL, on every record.
    integrals = [global_integral(fields["yfsst"][k], w) for k in range(nrec)]
    results["C3_global_integral"] = {
        "statement": "the field is orthogonal to the constant on a Gauss-"
                     "Legendre quadrature, so a prescribed advection carrying "
                     "it adds no net energy",
        "written_field_w_m2": global_integral(written, w),
        "max_abs_over_records_w_m2": float(np.abs(integrals).max()),
        "criterion_w_m2": TOL_OUTPUT_INTEGRAL_W_M2,
        "pass": float(np.abs(integrals).max()) <= TOL_OUTPUT_INTEGRAL_W_M2,
    }

    # C4 LAND. `getflxco` does not mask; `addfc` does.
    land = fields["yls"][0] >= 0.5
    ocean = ~land
    land_err = (float(np.abs(fields["yfsst"][:, land] - written[None, land]).max())
                if land.any() else float("nan"))
    results["C4_land_is_not_masked"] = {
        "statement": "getflxco interpolates over the whole grid and only addfc "
                     "restricts to yls < 1, so the reported accumulation "
                     "carries the field over land as well",
        "land_cells": int(land.sum()),
        "max_abs_error_over_land_w_m2": land_err,
        "criterion_w_m2": TOL_REPRODUCTION_W_M2,
        "pass": bool(land.any() and land_err <= TOL_REPRODUCTION_W_M2),
        "ocean_only_integral_of_written_field_w_m2": masked_mean(written, w, ocean),
        "note": "the ocean-only integral is NOT zero, and a real forcing "
                "artifact has to have zero integral over the OCEAN rather than "
                "over the globe, because the ocean is what addfc applies it to",
    }

    # C5 DELIVERY. One fitted proportionality across every ice-free ocean cell
    # and every record, and the fitted constant against the declared model.
    strict = ocean & (fields["yicec"].max(0) <= 0.0) \
        & (np.abs(fields["yiflux"]).max(0) == 0.0)
    dT = fields["ysst"][1:] - fields["ysst"][:-1]
    flux = fields["yheat"][1:] + fields["yfsst"][1:]
    delivery = _fit_delivery(dT, flux, strict, run_dir, mixed_layer_m,
                             record_seconds)
    delivery["statement"] = (
        "on an ice-free ocean cell mksst and addfc give "
        "CRHOS*CPS*mld*dSST = (yheat + yfsst)*dt_record, so the temperature "
        "change is ONE constant times the total flux and that constant is "
        "dt_record/(CRHOS*CPS*mld)")
    results["C5_delivery"] = delivery

    # The record clock the model itself stamped, so the declared interval is
    # checked against the stream rather than taken on trust. Words three and
    # four of the header are the date and the time of day.
    stamps = [(int(h[2]), int(h[3])) for h in heads]
    results["record_headers"] = {"first": stamps[0], "last": stamps[-1],
                                 "records": nrec}

    checks = [results["C1_reproduction"]["pass"],
              results["C3_global_integral"]["pass"],
              results["C4_land_is_not_masked"]["pass"],
              results["C5_delivery"]["pass"],
              results["C5_delivery"].get("record_seconds_pass", True),
              all(v["pass"] for v in off.values())]
    results["verdict"] = "PASS" if all(checks) else "FAIL"
    results["prediction"] = prediction
    return results


def check_extra(run_dir: Path, prediction: dict, record_seconds: float | None,
                mixed_layer_m: float | None) -> dict:
    """The monthly and the two ice arms. Same stream, different identities."""
    arm = prediction["arm"]
    nlat, nlon = int(prediction["nlat"]), int(prediction["nlon"])
    _, w, _ = grid(nlat)
    fields, heads = read_stream(stream_paths(run_dir), OCEAN_CODES)
    if fields["yfsst"].shape[1:] != (nlat, nlon):
        raise SystemExit(f"the stream is on {fields['yfsst'].shape[1:]} and the "
                         f"staged field on {(nlat, nlon)}")
    nrec = fields["yfsst"].shape[0]
    results: dict[str, dict] = {"terms_that_must_be_off": _terms_that_must_be_off(fields)}
    land = fields["yls"][0] >= 0.5
    ocean = ~land
    checks: list[bool] = [all(v["pass"] for v in results["terms_that_must_be_off"].values())]

    if arm == "monthly":
        written = read_sra_months(run_dir / prediction["files"]["flux"], nlat, nlon)
        if digest(written) != prediction["field"]["sha256_months_as_written"]:
            raise SystemExit(
                f"{prediction['files']['flux']} is not the file this prediction "
                "was written for. The interpolation cannot be verified against a "
                "field that changed between staging and the run.")
        scales = np.asarray(prediction["field"]["monthly_scales"], dtype=float)
        pattern = written[0] / scales[0]
        envelope = expand_to_fourteen(written)
        lo, hi = envelope.min(axis=0), envelope.max(axis=0)

        # M1 SEPARABILITY. Every month is the same pattern scaled, so a convex
        # combination of two months is that pattern scaled, and the record must
        # be one scalar times it at every cell.
        denom = float((pattern * pattern).sum())
        recovered = np.array([float((fields["yfsst"][k] * pattern).sum()) / denom
                              for k in range(nrec)])
        worst_rel, worst_at = 0.0, -1
        for k in range(nrec):
            residual = fields["yfsst"][k] - recovered[k] * pattern
            scale = abs(recovered[k]) * float(np.abs(pattern).max())
            rel = float(np.abs(residual).max() / scale) if scale else float("inf")
            if rel > worst_rel:
                worst_rel, worst_at = rel, k
        results["M1_separability"] = {
            "statement": "getflxco returns a convex combination of two months and "
                         "every month is the same spatial pattern times a scalar, "
                         "so each record must be that pattern times ONE scalar at "
                         "every cell. A per-cell mixing error, a month/space "
                         "transposition and a partially updated field all break "
                         "this and none of them breaks the envelope",
            "worst_relative_residual": worst_rel, "at_record": worst_at,
            "criterion_relative": TOL_SEPARABILITY_RELATIVE,
            "pass": bool(worst_rel <= TOL_SEPARABILITY_RELATIVE),
        }
        checks.append(results["M1_separability"]["pass"])

        # M2 THE ENVELOPE. The identity the calendar cannot influence.
        over = float(np.max(fields["yfsst"] - hi[None, :, :]))
        under = float(np.max(lo[None, :, :] - fields["yfsst"]))
        excess = max(over, under)
        results["M2_envelope"] = {
            "statement": "the month weights lie in [0, 1] whatever the calendar "
                         "says, so the interpolated value is a convex combination "
                         "of two of the fourteen months and must lie inside their "
                         "envelope at every step and every cell",
            "max_above_envelope_w_m2": over, "max_below_envelope_w_m2": under,
            "criterion_w_m2": TOL_REPRODUCTION_W_M2,
            "pass": bool(excess <= TOL_REPRODUCTION_W_M2),
        }
        checks.append(results["M2_envelope"]["pass"])

        # M3 GLOBAL INTEGRAL, at every record and not only in the mean.
        integrals = [global_integral(fields["yfsst"][k], w) for k in range(nrec)]
        results["M3_global_integral"] = {
            "statement": "every month is a scale of one pattern orthogonal to the "
                         "constant on a Gauss-Legendre quadrature, so the "
                         "prescribed advection adds no net energy at any point in "
                         "the cycle",
            "max_abs_over_records_w_m2": float(np.abs(integrals).max()),
            "criterion_w_m2": TOL_OUTPUT_INTEGRAL_W_M2,
            "pass": bool(float(np.abs(integrals).max()) <= TOL_OUTPUT_INTEGRAL_W_M2),
        }
        checks.append(results["M3_global_integral"]["pass"])

        # M4 THE SERIES ACTUALLY VARIES. Without this the arm passes vacuously.
        declared_span = float(scales.max() - scales.min())
        observed_span = float(recovered.max() - recovered.min())
        fraction = observed_span / declared_span if declared_span else float("nan")
        results["M4_series_varies"] = {
            "statement": "a model that ignored the month index and used one record "
                         "forever would pass M1, M2 and M3. An interpolation that "
                         "was never exercised is not one that was verified",
            "declared_scale_span": declared_span,
            "observed_scale_span": observed_span,
            "fraction_of_declared": fraction,
            "criterion_fraction": MONTHLY_MIN_OBSERVED_SPAN_FRACTION,
            "pass": bool(fraction >= MONTHLY_MIN_OBSERVED_SPAN_FRACTION),
        }
        checks.append(results["M4_series_varies"]["pass"])

        # M5 DELIVERY, on the ice-free ocean this arm is constructed to have.
        strict = ocean & (fields["yicec"].max(0) <= 0.0) \
            & (np.abs(fields["yiflux"]).max(0) == 0.0)
        dT = fields["ysst"][1:] - fields["ysst"][:-1]
        flux = fields["yheat"][1:] + fields["yfsst"][1:]
        delivery = _fit_delivery(dT, flux, strict, run_dir, mixed_layer_m,
                                 record_seconds)
        delivery["statement"] = (
            "the same identity the channel arm closed, on a flux that now varies "
            "in time, which conditions the fit far better than a constant one")
        results["M5_delivery"] = delivery
        checks.append(delivery["pass"])
        checks.append(delivery.get("record_seconds_pass", True))

        # REPORTED, not a bar: the stronger piecewise-linear statement needs the
        # calendar, and that is a check on the calendar port rather than on this
        # channel. The series and its extremes are here for whoever owns it.
        results["recovered_series"] = {
            "statement": "REPORTED and not judged. The breakpoints of this series "
                         "are where momint changes month, which is a statement "
                         "about the calendar and not about getflxco",
            "scalars": [float(x) for x in recovered],
            "declared_months": [float(x) for x in scales],
            "min": float(recovered.min()), "max": float(recovered.max()),
        }
    else:
        written = read_sra(run_dir / prediction["files"]["flux"], nlat, nlon)
        if digest(written) != prediction["field"]["sha256_as_written"]:
            raise SystemExit(
                f"{prediction['files']['flux']} is not the file this prediction "
                "was written for.")
        water = sea_water.constants(run_dir)
        freezing = float(water["TFREEZE"])
        iced_any = (fields["yicec"].max(0) > 0.0) & ocean
        iced_always = (fields["yicec"].min(0) > 0.0) & ocean
        free_always = (fields["yicec"].max(0) <= 0.0) & ocean
        ocean_cells = int(ocean.sum())

        # I1 THE BRANCH IS REACHED. Vacuous passing is the failure mode an ice
        # arm has that an ice-free one does not.
        reached = {
            "statement": "the arm has to grow the ice it is built to exercise, or "
                         "every identity below is the ice-free one again",
            "ocean_cells": ocean_cells,
            "iced_in_some_record": int(iced_any.sum()),
            "iced_fraction": float(iced_any.sum() / ocean_cells) if ocean_cells else float("nan"),
            "criterion_fraction": MIN_ICED_OCEAN_FRACTION,
            "pass": bool(ocean_cells and iced_any.sum() / ocean_cells >= MIN_ICED_OCEAN_FRACTION),
        }
        if arm == "ice_clim":
            ice_file = prediction["files"].get("ice_climatology")
            clim = read_sra(run_dir / ice_file, nlat, nlon) > 0.0
            overlap = iced_any & clim
            reached["climatological_ice_cells"] = int((clim & ocean).sum())
            reached["overlap_cells"] = int(overlap.sum())
            reached["overlap_fraction"] = (float(overlap.sum() / ocean_cells)
                                           if ocean_cells else float("nan"))
            reached["statement"] += (". Branch (b) needs the modelled ice and the "
                                     "climatological ice on the SAME cell, so the "
                                     "overlap is what has to clear the bar")
            reached["pass"] = bool(
                reached["pass"] and ocean_cells
                and overlap.sum() / ocean_cells >= MIN_ICED_OCEAN_FRACTION)
        results["I1_branch_reached"] = reached
        checks.append(reached["pass"])

        # I2 THE ICE-FREE CONTROL SURVIVES.
        results["I2_ice_free_control"] = {
            "statement": "branch (a) has to run in the same segment, so a delivery "
                         "failure localises to the ice branches rather than to the "
                         "slab integration they share",
            "ice_free_in_every_record": int(free_always.sum()),
            "fraction": float(free_always.sum() / ocean_cells) if ocean_cells else float("nan"),
            "criterion_fraction": MIN_ICE_FREE_OCEAN_FRACTION,
            "pass": bool(ocean_cells and free_always.sum() / ocean_cells >= MIN_ICE_FREE_OCEAN_FRACTION),
        }
        checks.append(results["I2_ice_free_control"]["pass"])

        # I3 THE THREE-TERM SLAB BUDGET.
        dT = fields["ysst"][1:] - fields["ysst"][:-1]
        flux = (SLAB_BUDGET_SIGNS["yheat"] * fields["yheat"][1:]
                + SLAB_BUDGET_SIGNS["yfsst"] * fields["yfsst"][1:]
                + SLAB_BUDGET_SIGNS["yiflux"] * fields["yiflux"][1:])
        delivery = _fit_delivery(dT, flux, ocean, run_dir, mixed_layer_m,
                                 record_seconds)
        delivery["statement"] = (
            "every flux mksst, mkiflux, addfc and mkiflx touch is either applied "
            "to the slab or charged to yiflux, and yiflux carries the NEGATIVE of "
            "what was withheld, so CRHOS*CPS*mld*dSST = (yheat + yfsst + yiflux)"
            "*dt_record on every ocean cell, iced or not. This is 'what addfc "
            "withholds equals what mkiflx adds' in the quantities the stream "
            "carries: yifluxr is not one of them")
        delivery["declared_signs"] = SLAB_BUDGET_SIGNS
        # If it fails, name the failure instead of reporting a residual. The pass
        # is on the declared signs and on nothing else.
        alternatives = {}
        for sign in (1.0, -1.0):
            for heat in (1.0, -1.0):
                alt = (heat * fields["yheat"][1:] + fields["yfsst"][1:]
                       + sign * fields["yiflux"][1:])
                alternatives[f"yheat{heat:+.0f}_yiflux{sign:+.0f}"] = \
                    _fit_delivery(dT, alt, ocean, run_dir, mixed_layer_m,
                                  record_seconds)["relative_residual"]
        delivery["residual_under_each_sign_assignment"] = alternatives
        results["I3_slab_budget_with_the_ice_term"] = delivery
        checks.append(delivery["pass"])
        checks.append(delivery.get("record_seconds_pass", True))

        # I4 THE ICE TERM IS MATERIAL.
        peak = float(np.abs(fields["yiflux"]).max())
        results["I4_ice_term_is_material"] = {
            "statement": "an identity evaluated where the new term is zero is the "
                         "identity the channel arm already closed",
            "peak_abs_yiflux_w_m2": peak,
            "criterion_w_m2": MIN_ICE_FLUX_W_M2,
            "pass": bool(peak >= MIN_ICE_FLUX_W_M2),
        }
        checks.append(results["I4_ice_term_is_material"]["pass"])

        # I5 THE FREEZING CLAMP.
        under = float(np.max(freezing - fields["ysst"][:, ocean])) if ocean.any() else float("nan")
        results["I5_freezing_clamp"] = {
            "statement": "addfc's final clamp and mkiflux both hold the modelled "
                         "sea surface at or above the freezing point and charge "
                         "the difference to the ice, so a reported value below it "
                         "is a clamp that did not fire",
            "freezing_point_k": freezing, "sea_water": water,
            "max_below_freezing_k": under,
            "criterion_k": TOL_FREEZING_CLAMP_K,
            "pass": bool(under <= TOL_FREEZING_CLAMP_K),
        }
        checks.append(results["I5_freezing_clamp"]["pass"])

        # I6 THE CORRECTION IS STILL REPORTED UNCHANGED. getflxco does not mask
        # and the field is constant in time, so the ice branches must not alter
        # what code 903 reports -- they alter what is DONE with it.
        err = float(np.abs(fields["yfsst"] - written[None, :, :]).max())
        results["I6_reported_correction_unchanged"] = {
            "statement": "yfssta accumulates yfsst2, which is what getflxco "
                         "interpolated and not what addfc applied, so the ice "
                         "branches change the delivery and not the report",
            "max_abs_error_w_m2": err,
            "criterion_w_m2": TOL_REPRODUCTION_W_M2,
            "pass": bool(err <= TOL_REPRODUCTION_W_M2),
        }
        checks.append(results["I6_reported_correction_unchanged"]["pass"])
        results["ice_cover"] = {
            "iced_in_every_record": int(iced_always.sum()),
            "iced_in_some_record": int(iced_any.sum()),
            "ice_free_in_every_record": int(free_always.sum()),
        }

    stamps = [(int(h[2]), int(h[3])) for h in heads]
    results["record_headers"] = {"first": stamps[0], "last": stamps[-1],
                                 "records": nrec}
    results["verdict"] = "PASS" if all(checks) else "FAIL"
    results["prediction"] = prediction
    return results


# ---------------------------------------------------------------------------
# selftest
#
# THE MODEL IS NOT NEEDED TO KNOW WHETHER A CRITERION DISCRIMINATES. Every
# criterion above is a statement about a stream, so a stream that satisfies the
# statements exactly can be synthesised here and each criterion perturbed in one
# named way. A criterion that passes its clean stream and fails its own
# perturbation has been shown to have a right answer; one that does not is a
# defect in this instrument rather than in the model.
#
# This does NOT verify the model. It verifies the instrument, and it is what
# stands between staging an arm and running it.
# ---------------------------------------------------------------------------

SELFTEST_RECORD_SECONDS = 86400.0
SELFTEST_MIXED_LAYER_M = 50.0


def write_service_stream(path: Path, records: list[tuple[int, np.ndarray]],
                         nlat: int, nlon: int) -> None:
    """A PlaSim service stream: the layout `read_stream` reads, and no other."""
    with open(path, "wb") as handle:
        for index, (code, field) in enumerate(records):
            day = 20260101 + index
            head = np.array([code, 0, day, 0, nlon, nlat, 0, 0], dtype=np.int32)
            raw = head.tobytes()
            handle.write(struct.pack("i", len(raw)) + raw + struct.pack("i", len(raw)))
            body = np.ascontiguousarray(field, dtype=np.float32).tobytes()
            handle.write(struct.pack("i", len(body)) + body + struct.pack("i", len(body)))


def _synthetic(run_dir: Path, prediction: dict, nrec: int = 16) -> dict:
    """Fields that satisfy the arm's identities exactly, as (nrec, nlat, nlon)."""
    arm = prediction["arm"]
    nlat, nlon = int(prediction["nlat"]), int(prediction["nlon"])
    mu, _, _ = grid(nlat)
    capacity = sea_water.slab_heat_capacity(SELFTEST_MIXED_LAYER_M, run_dir)
    freezing = float(sea_water.constants(run_dir)["TFREEZE"])
    per_w = SELFTEST_RECORD_SECONDS / capacity

    # A land mask with real land in it, so the ocean-only statements are not
    # trivially the whole grid.
    land = np.zeros((nlat, nlon), dtype=bool)
    land[:, : nlon // 6] = True
    ocean = ~land
    zeros = np.zeros((nrec, nlat, nlon))

    # COLD-START SIZED, deliberately. The delivery criterion is a fraction of
    # the largest temperature change in the window, so a fixture with small
    # fluxes puts that bar under the float32 resolution of ysst and the
    # instrument's own MIN_DELIVERY_CHANGE_K guard fires on the CLEAN stream. A
    # segment that can answer the delivery question looks like this one.
    shape = np.repeat((-500.0 * mu ** 2 + 200.0)[None, :, None], nlon, axis=2)
    heat = shape * (0.85 ** np.arange(nrec))[:, None, None]

    if arm == "channel":
        written = read_sra(run_dir / prediction["files"]["flux"], nlat, nlon)
        yfsst = np.repeat(written[None, :, :], nrec, axis=0)
        yicec = zeros.copy()
        yiflux = zeros.copy()
    elif arm == "monthly":
        written = read_sra_months(run_dir / prediction["files"]["flux"], nlat, nlon)
        scales = np.asarray(prediction["field"]["monthly_scales"], dtype=float)
        pattern = written[0] / scales[0]
        # A convex combination of two adjacent months at every record, which is
        # what getflxco returns whatever the calendar says.
        phase = np.linspace(0.0, 12.0, nrec, endpoint=False)
        recovered = np.array([
            (1.0 - (t % 1.0)) * scales[int(t) % 12] + (t % 1.0) * scales[(int(t) + 1) % 12]
            for t in phase])
        yfsst = recovered[:, None, None] * pattern[None, :, :]
        yicec = zeros.copy()
        yiflux = zeros.copy()
    else:
        written = read_sra(run_dir / prediction["files"]["flux"], nlat, nlon)
        yfsst = np.repeat(written[None, :, :], nrec, axis=0)
        cold = latitude_field(nlat, nlon, COLD_SST_POLE_K, COLD_SST_EQUATOR_K)
        iced = (cold <= freezing) & ocean
        yicec = np.repeat((iced * 1.0)[None, :, :], nrec, axis=0)
        # On an iced cell the slab is held near the freezing point, so the
        # temperature change is small and yiflux carries what was withheld.
        yiflux = zeros.copy()

    ysst = np.empty((nrec, nlat, nlon))
    if arm in ("channel", "monthly"):
        ysst[0] = np.where(ocean, CLIMATOLOGICAL_SST_K, 0.0)
        for k in range(1, nrec):
            ysst[k] = ysst[k - 1] + np.where(
                ocean, (heat[k] + yfsst[k]) * per_w, 0.0)
    else:
        start = np.where(iced, freezing + 0.5, np.maximum(cold, freezing + 1.0))
        ysst[0] = np.where(ocean, start, 0.0)
        for k in range(1, nrec):
            # Ice-free cells take the whole flux; iced cells relax toward the
            # freezing point and yiflux is whatever the identity then requires.
            step = (heat[k] + yfsst[k]) * per_w
            target = np.where(iced, freezing + 0.5 * 0.8 ** k, ysst[k - 1] + step)
            target = np.maximum(target, freezing)
            ysst[k] = np.where(ocean, target, 0.0)
            required = (ysst[k] - ysst[k - 1]) / per_w - heat[k] - yfsst[k]
            yiflux[k] = np.where(iced, required, 0.0)
        yiflux[0] = yiflux[1]

    return {"yheat": heat, "yfsst": yfsst, "yiflux": yiflux, "yicec": yicec,
            "ysst": ysst, "yls": np.repeat((land * 1.0)[None, :, :], nrec, axis=0),
            "ydsst": zeros.copy(), "yqhd": zeros.copy(), "yfldo": zeros.copy(),
            "yclsst": np.repeat(
                np.where(land, 0.0, CLIMATOLOGICAL_SST_K)[None, :, :], nrec, axis=0)}


def _write_synthetic(run_dir: Path, prediction: dict, fields: dict) -> None:
    nlat, nlon = int(prediction["nlat"]), int(prediction["nlon"])
    order = [(code, name) for code, name in sorted(OCEAN_CODES.items())]
    records: list[tuple[int, np.ndarray]] = []
    nrec = fields["yfsst"].shape[0]
    for k in range(nrec):
        for code, name in order:
            records.append((code, fields[name][k]))
    write_service_stream(run_dir / "ocean_output", records, nlat, nlon)


def selftest(run_dir: Path) -> dict:
    """Every criterion of the staged arm, against a clean stream and its own break."""
    prediction = json.loads((run_dir / PREDICTION_NAME).read_text(encoding="utf-8"))
    arm = prediction["arm"]
    clean = _synthetic(run_dir, prediction)

    def verdict(fields: dict) -> dict:
        _write_synthetic(run_dir, prediction, fields)
        return check(run_dir, SELFTEST_RECORD_SECONDS, SELFTEST_MIXED_LAYER_M)

    def broken(fn):
        f = {k: v.copy() for k, v in clean.items()}
        fn(f)
        return f

    if arm == "channel":
        cases = [
            ("a record that is not the written field",
             broken(lambda f: f["yfsst"].__setitem__(
                 (3, slice(None), slice(None)), f["yfsst"][3] * 1.1)),
             "C1_reproduction", "pass"),
            ("a correction with a net global integral",
             broken(lambda f: f.__setitem__("yfsst", f["yfsst"] + 1.0)),
             "C3_global_integral", "pass"),
            ("a correction masked to the ocean, which getflxco does not do",
             broken(lambda f: f.__setitem__(
                 "yfsst", np.where(f["yls"] >= 0.5, 0.0, f["yfsst"]))),
             "C4_land_is_not_masked", "pass"),
            ("a slab that took twice the flux on half the grid",
             broken(lambda f: f.__setitem__(
                 "ysst", f["ysst"][0:1] + (f["ysst"] - f["ysst"][0:1])
                 * np.where(np.arange(f["ysst"].shape[1])[None, :, None]
                            < f["ysst"].shape[1] // 2, 2.0, 1.0))),
             "C5_delivery", "pass"),
            ("a slab that took twice the flux everywhere, which the slope absorbs",
             broken(lambda f: f.__setitem__(
                 "ysst", f["ysst"][0:1] + 2.0 * (f["ysst"] - f["ysst"][0:1]))),
             "C5_delivery", "record_seconds_pass"),
            ("a window whose changes are too small for the bar to sit above float32",
             broken(lambda f: f.__setitem__(
                 "ysst", f["ysst"][0:1] + 0.001 * (f["ysst"] - f["ysst"][0:1]))),
             "C5_delivery", "discriminates"),
        ]
    elif arm == "monthly":
        cases = [
            ("a record that is not one scalar times the pattern",
             broken(lambda f: f["yfsst"].__setitem__(
                 (3, slice(None), slice(None)),
                 f["yfsst"][3] + 0.5 * np.sin(np.arange(f["yfsst"].shape[2]))[None, :])),
             "M1_separability"),
            ("a record outside the fourteen-month envelope",
             broken(lambda f: f["yfsst"].__setitem__(
                 (4, slice(None), slice(None)), f["yfsst"][4] * 3.0)),
             "M2_envelope"),
            ("a correction with a net global integral",
             broken(lambda f: f.__setitem__("yfsst", f["yfsst"] + 1.0)),
             "M3_global_integral"),
            ("an interpolation that never moved off one month",
             broken(lambda f: f.__setitem__(
                 "yfsst", np.repeat(f["yfsst"][0:1], f["yfsst"].shape[0], axis=0))),
             "M4_series_varies"),
            # TWO BREAKS, because the delivery test is two tests. A uniform
            # doubling is absorbed exactly by the fitted slope and leaves no
            # scatter at all; what catches it is the fitted record interval
            # against the declared one. A doubling on half the grid cannot be
            # absorbed and is what the scatter criterion is for.
            ("a slab that took twice the flux on half the grid",
             broken(lambda f: f.__setitem__(
                 "ysst", f["ysst"][0:1] + (f["ysst"] - f["ysst"][0:1])
                 * np.where(np.arange(f["ysst"].shape[1])[None, :, None]
                            < f["ysst"].shape[1] // 2, 2.0, 1.0))),
             "M5_delivery", "pass"),
            ("a slab that took twice the flux everywhere, which the slope absorbs",
             broken(lambda f: f.__setitem__(
                 "ysst", f["ysst"][0:1] + 2.0 * (f["ysst"] - f["ysst"][0:1]))),
             "M5_delivery", "record_seconds_pass"),
        ]
    else:
        cases = [
            ("no modelled ice anywhere",
             broken(lambda f: f.__setitem__("yicec", f["yicec"] * 0.0)),
             "I1_branch_reached"),
            ("ice over the whole ocean, so branch (a) never runs",
             broken(lambda f: f.__setitem__(
                 "yicec", np.where(f["yls"] >= 0.5, 0.0, 1.0))),
             "I2_ice_free_control"),
            ("the ice term entering the budget with the wrong sign",
             broken(lambda f: f.__setitem__("yiflux", -f["yiflux"])),
             "I3_slab_budget_with_the_ice_term"),
            ("an ice term that is identically zero",
             broken(lambda f: f.__setitem__("yiflux", f["yiflux"] * 0.0)),
             "I4_ice_term_is_material"),
            ("a modelled sea surface below the freezing point",
             broken(lambda f: f["ysst"].__setitem__(
                 (2, slice(None), slice(None)), f["ysst"][2] - 1.0)),
             "I5_freezing_clamp"),
            ("a reported correction the ice branches altered",
             broken(lambda f: f["yfsst"].__setitem__(
                 (5, slice(None), slice(None)), f["yfsst"][5] * 0.5)),
             "I6_reported_correction_unchanged"),
            ("a slab that took twice the flux everywhere, which the slope absorbs",
             broken(lambda f: f.__setitem__(
                 "ysst", f["ysst"][0:1] + 2.0 * (f["ysst"] - f["ysst"][0:1]))),
             "I3_slab_budget_with_the_ice_term", "record_seconds_pass"),
        ]

    results = {"arm": arm, "cases": []}
    base = verdict(clean)
    results["cases"].append({
        "case": "the synthetic stream that satisfies every identity",
        "expected": "PASS", "verdict": base["verdict"], "criterion": None,
        "pass": base["verdict"] == "PASS",
        "failing": sorted(k for k, v in base.items()
                          if isinstance(v, dict) and v.get("pass") is False)})
    for case in cases:
        label, fields, criterion = case[0], case[1], case[2]
        key = case[3] if len(case) > 3 else "pass"
        got = verdict(fields)
        entry = got.get(criterion, {})
        results["cases"].append({
            "case": label, "expected": f"{criterion}.{key} FAIL",
            "verdict": got["verdict"], "criterion": criterion, "key": key,
            "pass": entry.get(key) is False})
    (run_dir / "ocean_output").unlink(missing_ok=True)
    results["verdict"] = "PASS" if all(c["pass"] for c in results["cases"]) else "FAIL"
    return results


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="command", required=True)
    s = sub.add_parser("stage", help="write the field and the prediction")
    s.add_argument("run_dir", type=Path)
    s.add_argument("--rung", required=True, choices=sorted(rungs.RUNGS))
    s.add_argument("--arm", default="channel", choices=ARMS,
                   help="which path to isolate; each arm is its own run directory")
    c = sub.add_parser("check", help="answer every criterion against the stream")
    c.add_argument("run_dir", type=Path)
    c.add_argument("--record-seconds", type=float, default=None,
                   help="nout * solar_day / ntspd; C5's second half needs it")
    c.add_argument("--mixed-layer-m", type=float, default=None,
                   help="oceanmod_nl MLDEPTH; read from the run directory if absent")
    c.add_argument("--output", type=Path, default=None)
    t = sub.add_parser("selftest", help="every criterion of the staged arm, "
                                        "against a clean synthetic stream and "
                                        "its own named break")
    t.add_argument("run_dir", type=Path)
    a = ap.parse_args()

    if a.command == "stage":
        out = stage(a.run_dir, a.rung, a.arm)
        print(json.dumps(out, indent=2))
        return

    if a.command == "selftest":
        out = selftest(a.run_dir)
        print(json.dumps(out, indent=2))
        sys.exit(0 if out["verdict"] == "PASS" else 2)

    mld = a.mixed_layer_m
    if mld is None:
        mld = namelist_value(a.run_dir / "oceanmod_namelist", "MLDEPTH")
    report = check(a.run_dir, a.record_seconds, mld)
    path = a.output or (a.run_dir / "ocean_flux_channel_report.json")
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    trimmed = {k: v for k, v in report.items() if k != "prediction"}
    print(json.dumps(trimmed, indent=2))
    print(f"report: {path}")
    sys.exit(0 if report["verdict"] == "PASS" else 1)


if __name__ == "__main__":
    main()
