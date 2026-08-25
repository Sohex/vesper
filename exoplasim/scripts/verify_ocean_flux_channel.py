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

USAGE, in two phases, because the model runs between them:

    python exoplasim/scripts/verify_ocean_flux_channel.py stage RUNDIR --rung T42
    ...run the model in RUNDIR with nfluko = 1...
    python exoplasim/scripts/verify_ocean_flux_channel.py check RUNDIR

`stage` writes `N<nlat>_surf_0903.sra`, a companion sea surface temperature
climatology at code 169 that `oceanini` refuses `nfluko` without, and
`ocean_flux_channel_prediction.json`, which carries the field's digest, the
grid, the predicted integrals and the criteria. `check` reads that file and the
`ocean_output` stream and answers each criterion. It exits non-zero if any
criterion fails.

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

def stage(run_dir: Path, rung: str) -> dict:
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
    }
    (run_dir / PREDICTION_NAME).write_text(
        json.dumps(prediction, indent=2) + "\n", encoding="utf-8")
    return prediction


# ---------------------------------------------------------------------------
# check
# ---------------------------------------------------------------------------

def check(run_dir: Path, record_seconds: float | None,
          mixed_layer_m: float | None) -> dict:
    prediction = json.loads((run_dir / PREDICTION_NAME).read_text(encoding="utf-8"))
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
    m = np.broadcast_to(strict[None, :, :], dT.shape)
    x, y = flux[m], dT[m]
    slope = float((x * y).sum() / (x * x).sum()) if x.size and (x * x).sum() else float("nan")
    resid = y - slope * x
    scale = float(np.abs(y).max()) if y.size else float("nan")
    rel = float(np.abs(resid).max() / scale) if scale else float("nan")
    water = sea_water.constants(run_dir)
    fitted_seconds = slope * sea_water.slab_heat_capacity(
        mixed_layer_m if mixed_layer_m else float("nan"), run_dir)
    delivery = {
        "statement": "on an ice-free ocean cell mksst and addfc give "
                     "CRHOS*CPS*mld*dSST = (yheat + yfsst)*dt_record, so the "
                     "temperature change is ONE constant times the total flux "
                     "and that constant is dt_record/(CRHOS*CPS*mld)",
        "cells": int(strict.sum()), "record_pairs": int(dT.shape[0]),
        "samples": int(x.size),
        "fitted_k_per_w_m2": slope,
        "max_abs_residual_k": float(np.abs(resid).max()) if resid.size else None,
        "largest_change_k": scale,
        "relative_residual": rel,
        "criterion_relative": TOL_DELIVERY_RELATIVE,
        "pass": bool(x.size and rel <= TOL_DELIVERY_RELATIVE),
        "mixed_layer_depth_m": mixed_layer_m,
        "sea_water": water,
        "fitted_record_seconds": fitted_seconds,
        "declared_record_seconds": record_seconds,
    }
    if record_seconds and np.isfinite(fitted_seconds):
        off_by = abs(fitted_seconds / record_seconds - 1.0)
        delivery["record_seconds_relative_error"] = float(off_by)
        delivery["record_seconds_pass"] = bool(off_by <= TOL_RECORD_SECONDS_RELATIVE)
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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="command", required=True)
    s = sub.add_parser("stage", help="write the field and the prediction")
    s.add_argument("run_dir", type=Path)
    s.add_argument("--rung", required=True, choices=sorted(rungs.RUNGS))
    c = sub.add_parser("check", help="answer every criterion against the stream")
    c.add_argument("run_dir", type=Path)
    c.add_argument("--record-seconds", type=float, default=None,
                   help="nout * solar_day / ntspd; C5's second half needs it")
    c.add_argument("--mixed-layer-m", type=float, default=None,
                   help="oceanmod_nl MLDEPTH; read from the run directory if absent")
    c.add_argument("--output", type=Path, default=None)
    a = ap.parse_args()

    if a.command == "stage":
        out = stage(a.run_dir, a.rung)
        print(json.dumps(out, indent=2))
        return

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
