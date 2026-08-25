#!/usr/bin/env python3
"""Close the planetary energy budget against the model's PROGNOSTIC STATE.

Every other closure in this project compares one flux diagnostic with another.
That can only ever show the diagnostics agreeing with themselves: `ntr` equals
`rst + rlut` by construction, the radiative heating terms are the divergence of
the same flux profile the top-of-atmosphere fluxes are read off, and the pair
`rst`/`rlut` sum to `ntr` and therefore cannot be cross-checked against each
other. `docs/src/practice/failure-modes.md` class 17 is the name for that: a check with no
possible failure is not a check.

This script uses a quantity the radiation code does not supply. Energy
conservation says

    d/dt (planetary heat content)  =  net top-of-atmosphere radiation

with the heat content assembled from the prognostic state -- the ocean mixed
layer, the sea-ice and snow mass as latent heat, the soil column, the
atmosphere's enthalpy and its column vapour. That is an identity, not a
tolerance, so it has a right answer and it can fail: if a slab or sea-ice
storage term were carrying the reported imbalance, this would show it at the
full size of the imbalance.

The second check is the same shape applied to the incoming shortwave. Over a
full orbit the time mean of the inverse-square distance factor is exactly
1/sqrt(1 - e^2), so

    <rst - rsut>  =  GSOL0 / (4 sqrt(1 - e^2))

with GSOL0 and the eccentricity taken from the run's own namelist rather than
from the radiation output. That tests the disc average, the eccentricity
weighting, the zenith-angle integration and the radiation call frequency
together, against a number the radiation code never sees.

Both are reported with the spread across orbits, because a single orbit of this
model carries a seasonal storage swing of several W/m2 and a one-orbit reading
would be noise. Run it on a block of orbits written with `NLOWIO = 0`.

Usage:

    python exoplasim/scripts/close_state_energy.py <run_dir> --first 67 --last 76
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
from numpy.polynomial.legendre import leggauss

from _paths import ANALYSIS  # noqa: F401  (also puts lib/ on sys.path)

import sea_water  # noqa: E402  from lib/, via _paths

# PlaSim's own constants, from plasim/src. They are properties of the compiled
# model, not of this planet, so they are read from the source of truth that
# applies: the model. Vesper's own parameters (gravity, gas constant, stellar
# flux, eccentricity) are read from the run's namelists below, because this is a
# statement about the run being assessed and not about the current config.
ACPV = 1870.0        # plasimmod.f90: specific heat of water vapour, J/kg/K
ALV = 2.5008e6       # plasimmod.f90: latent heat of vaporisation, J/kg
ALS = 2.8345e6       # plasimmod.f90: latent heat of sublimation, J/kg
ALF = ALS - ALV      # fusion
# SEA WATER IS READ, NOT WRITTEN DOWN. Density and specific heat are icemod_nl
# keys that icemod passes to oceanini, so oceanmod does not own them and a run
# can set them. The literals that stood here said 1030 and 4180, and 4180 is
# fresh water at about 25 C, which the model left behind for sea water's value
# at S = 34.7 and its freezing point. Every mixed-layer heat content built from
# the pair was too large by their ratio, and assess_convergence.py imports these
# two names to build the slab capacity its convergence verdict rests on.
_SEA_WATER = sea_water.constants()
CRHOS = _SEA_WATER["CRHOS"]   # icemod.f90/icemod_nl: density of sea water, kg/m3
CPS = _SEA_WATER["CPS"]       # icemod.f90/icemod_nl: specific heat, J/kg/K
CRHOI = 920.0        # oceanmod.f90: density of sea ice, kg/m3
TMELT = 273.16       # icemod.f90: freezing point
SOILCAP = 2.4e6      # landmod.f90: soil heat capacity, J/m3/K
DSOILZ = np.array([0.4, 0.8, 1.6, 3.2, 6.4])   # landmod.f90: soil layer thicknesses, m
RHO_WATER = 1000.0   # snow depth is reported as metres water equivalent
AKAP_DEFAULT = 0.286     # p_earth.f90, not overridden by this project's namelist


def namelist_values(path: Path) -> dict[str, float]:
    """Every `KEY = value` in a Fortran namelist file, upper-cased."""
    values: dict[str, float] = {}
    if not path.is_file():
        return values
    for key, raw in re.findall(r"(\w+)\s*=\s*([-+0-9.eEdD]+)",
                               path.read_text(encoding="utf-8")):
        try:
            values[key.upper()] = float(raw.replace("d", "e").replace("D", "e"))
        except ValueError:
            continue
    return values


def annual_files(run_dir: Path, first: int, last: int) -> list[tuple[int, Path]]:
    files = []
    for index in range(first, last + 1):
        path = run_dir / f"MOST.{index:05d}.nc"
        if not path.is_file():
            raise SystemExit(f"missing annual output {path}")
        files.append((index, path))
    return files


def global_mean(field: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Gaussian-weighted global mean, leaving the leading (time) axis alone."""
    return np.sum(field * weights[..., None], axis=(-2, -1)) / (2.0 * field.shape[-1])


def heat_content(nc: Dataset, gravity: float, acpd: float) -> dict[str, np.ndarray]:
    """Planetary heat content per unit area, J/m2, per output bin.

    Reference level is arbitrary and cancels in the time derivative; what has to
    be right is that every reservoir with a non-negligible trend is present and
    carries its true capacity.
    """
    read = lambda name: np.asarray(nc[name][:], dtype=float)
    adv = ACPV / acpd - 1.0
    land = read("lsm")
    ocean = 1.0 - land
    ts = read("ts")
    dsigma = np.diff(read("levp"))
    surface_pressure = read("ps") * 100.0                # hPa to Pa

    # Column integral of cp T dp/g. For a hydrostatic atmosphere this is the
    # internal plus potential energy together, so no separate geopotential term
    # belongs here. Kinetic energy is order 1e6 J/m2 against 2e9 and is left out.
    atmosphere = np.sum(
        read("ta") * acpd * (1.0 + adv * read("hus"))
        * surface_pressure[:, None] * dsigma[None, :, None, None], axis=1) / gravity
    # THE MOISTURE FIXER MOVES THIS TERM BETWEEN COLUMNS, every timestep and
    # unconditionally: `miscmod.f90:fixer` removes negative specific humidity by
    # borrowing it, first within a column, then along a latitude row, then
    # globally. The global step conserves the column integral, so this reservoir
    # is intact for a GLOBAL closure and this script is a global closure. It is
    # not intact for a regional or per-column one, and in the branch where the
    # global deficit exceeds the global surplus the fixer is a sink rather than a
    # redistribution. world-bsp.
    vapour = ALV * read("prw")

    # Over an ice-covered ocean cell `ts` is the ice surface, not the mixed
    # layer, which sits at freezing. Sea ice is a small fraction here but the
    # ice-surface temperature is not a small perturbation on it.
    sea_ice = read("sic")
    mixed_layer = (CRHOS * CPS * read("mld")
                   * ((1.0 - sea_ice) * ts + sea_ice * TMELT) * ocean)
    ice = -ALF * CRHOI * read("sit") * sea_ice * ocean
    snow = -ALF * RHO_WATER * read("snd")

    # Five soil layers; the model writes layers 2-5 (codes 207, 208, 209, 170)
    # and not layer 1, whose 0.4 m is tied to the surface temperature by the
    # implicit top-layer solve. `ts` stands in for it, worth 3% of the column.
    soil_profile = np.stack([ts, read("tso2"), read("tso3"), read("tso4"),
                             read("tsod")], axis=1)
    soil = SOILCAP * np.sum(soil_profile * DSOILZ[None, :, None, None], axis=1) * land

    return {"atmosphere": atmosphere, "vapour": vapour, "mixed_layer": mixed_layer,
            "sea_ice": ice, "snow": snow, "soil": soil}


def state_energy(run_dir: Path, first: int, last: int) -> dict:
    """The state-energy closure for one window of orbits.

    Extracted so `assess_convergence.py` can read the storage term without
    computing a second version of it. The convergence criterion now PASSES on
    this quantity, so a second implementation of it would be a second answer to
    the question the criterion asks, which is how this project has previously
    ended up with three values for one number.
    """
    planet = namelist_values(run_dir / "planet_namelist")
    gravity = planet.get("GA")
    gascon = planet.get("GASCON")
    solar_constant = planet.get("GSOL0")
    eccentricity = planet.get("ECCEN")
    if None in (gravity, gascon, solar_constant, eccentricity):
        raise SystemExit(f"{run_dir/'planet_namelist'} does not declare GA, GASCON, "
                         "GSOL0 and ECCEN; this run cannot be closed against itself")
    acpd = gascon / AKAP_DEFAULT

    manifest_path = run_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    _days = manifest.get("derived_parameters", {}).get("orbital_year_earth_days")
    if _days is None:
        raise SystemExit(f"{manifest_path} lacks derived_parameters."
                         "orbital_year_earth_days; refusing a hardcoded year")
    orbit_seconds = float(_days) * 86400.0

    reservoirs = ["atmosphere", "vapour", "mixed_layer", "sea_ice", "snow", "soil"]
    fluxes = ["ntr", "hfns", "rst", "rsut", "rlut", "rss", "rls", "hfss", "hfls",
              "ts", "mld", "lsm"]
    per_orbit: list[dict[str, float]] = []
    for index, path in annual_files(run_dir, first, last):
        with Dataset(path) as nc:
            weights = leggauss(len(nc.dimensions["lat"]))[1][::-1]
            record = {"orbit": index}
            content = heat_content(nc, gravity, acpd)
            for name in reservoirs:
                record[name] = float(global_mean(content[name], weights).mean())
            for name in fluxes:
                record[name] = float(
                    global_mean(np.asarray(nc[name][:], dtype=float), weights).mean())
        per_orbit.append(record)

    orbits = np.array([r["orbit"] for r in per_orbit], dtype=float)
    if orbits.size < 3:
        raise SystemExit("a storage trend needs at least three orbits")
    series = {name: np.array([r[name] for r in per_orbit]) for name in reservoirs + fluxes}
    total = sum(series[name] for name in reservoirs)

    # Two independent estimators of the same derivative. A least-squares slope
    # uses every orbit and is what should be quoted; the endpoint difference is
    # immune to a curved approach and is reported so a disagreement between them
    # is visible rather than averaged away.
    trend = {name: float(np.polyfit(orbits, series[name], 1)[0] / orbit_seconds)
             for name in reservoirs}
    storage_fit = float(np.polyfit(orbits, total, 1)[0] / orbit_seconds)
    storage_endpoint = float((total[-1] - total[0])
                             / ((orbits[-1] - orbits[0]) * orbit_seconds))
    surface_total = sum(series[name] for name in ["mixed_layer", "sea_ice", "snow", "soil"])
    surface_storage = float(np.polyfit(orbits, surface_total, 1)[0] / orbit_seconds)
    atmosphere_total = series["atmosphere"] + series["vapour"]
    atmosphere_storage = float(np.polyfit(orbits, atmosphere_total, 1)[0] / orbit_seconds)

    mean_toa = float(series["ntr"].mean())
    mean_surface_flux = float(series["hfns"].mean())

    # The orbit-mean of (a/r)^2 is exactly 1/sqrt(1-e^2): equal areas in equal
    # times makes the time average of r^-2 equal to the reciprocal of the
    # ellipse's semi-latus-rectum product. Nothing here comes from the radiation.
    expected_insolation = solar_constant / 4.0 / np.sqrt(1.0 - eccentricity ** 2)
    reported_insolation = float((series["rst"] - series["rsut"]).mean())

    return {
        "manifest": manifest, "per_orbit": per_orbit, "series": series,
        "orbit_seconds": orbit_seconds,
        "planet_namelist": {"GA": gravity, "GASCON": gascon,
                            "GSOL0": solar_constant, "ECCEN": eccentricity},
        "storage_w_m2_least_squares": storage_fit,
        "storage_w_m2_endpoint": storage_endpoint,
        "mean_toa_w_m2": mean_toa,
        "residual_w_m2": mean_toa - storage_fit,
        "reservoir_trends_w_m2": trend,
        "surface_storage_w_m2": surface_storage,
        "mean_surface_flux_w_m2": mean_surface_flux,
        "surface_residual_w_m2": mean_surface_flux - surface_storage,
        "atmosphere_storage_w_m2": atmosphere_storage,
        "expected_insolation_w_m2": expected_insolation,
        "reported_insolation_w_m2": reported_insolation,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--first", type=int, required=True,
                        help="first orbit index of the window")
    parser.add_argument("--last", type=int, required=True,
                        help="last orbit index of the window")
    parser.add_argument("--output", type=Path, default=ANALYSIS / "energy")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()

    c = state_energy(run_dir, args.first, args.last)
    manifest = c["manifest"]
    per_orbit = c["per_orbit"]
    series = c["series"]
    orbit_seconds = c["orbit_seconds"]
    gravity = c["planet_namelist"]["GA"]
    gascon = c["planet_namelist"]["GASCON"]
    solar_constant = c["planet_namelist"]["GSOL0"]
    eccentricity = c["planet_namelist"]["ECCEN"]
    storage_fit = c["storage_w_m2_least_squares"]
    storage_endpoint = c["storage_w_m2_endpoint"]
    trend = c["reservoir_trends_w_m2"]
    surface_storage = c["surface_storage_w_m2"]
    atmosphere_storage = c["atmosphere_storage_w_m2"]
    mean_toa = c["mean_toa_w_m2"]
    mean_surface_flux = c["mean_surface_flux_w_m2"]
    expected_insolation = c["expected_insolation_w_m2"]
    reported_insolation = c["reported_insolation_w_m2"]
    orbits = np.array([r["orbit"] for r in per_orbit], dtype=float)

    report = {
        "schema_version": 1,
        "generator": "exoplasim/scripts/close_state_energy.py",
        "generator_sha256": hashlib.sha256(
            Path(__file__).read_bytes()).hexdigest(),
        "run_dir": str(run_dir),
        "run_id": manifest.get("run_id"),
        "source_build": manifest.get("source_build"),
        "config_sha256": manifest.get("config_sha256"),
        "geography": manifest.get("physical", {}).get("geography"),
        "window": {"first_orbit": args.first, "last_orbit": args.last,
                   "orbits": len(per_orbit), "orbit_seconds": orbit_seconds},
        "planet_namelist": {"GA": gravity, "GASCON": gascon,
                            "GSOL0": solar_constant, "ECCEN": eccentricity},
        "state_energy_closure": {
            "identity": "d(planetary heat content)/dt = mean net TOA radiation",
            "storage_w_m2_least_squares": storage_fit,
            "storage_w_m2_endpoint": storage_endpoint,
            "mean_toa_w_m2": mean_toa,
            "residual_w_m2": mean_toa - storage_fit,
            "reservoir_trends_w_m2": trend,
            "surface_storage_w_m2": surface_storage,
            "mean_surface_flux_w_m2": mean_surface_flux,
            "surface_residual_w_m2": mean_surface_flux - surface_storage,
            "atmosphere_storage_w_m2": atmosphere_storage,
            "atmosphere_flux_w_m2": mean_toa - mean_surface_flux,
            "atmosphere_residual_w_m2":
                (mean_toa - mean_surface_flux) - atmosphere_storage,
            # As a temperature, because that is the form the drift is visible in
            # and the form the convergence criteria are stated in. The heat
            # content is a global mean over an ocean that covers part of the
            # planet, so the divisor is the global-mean mixed-layer depth.
            "mixed_layer_drift_k_per_orbit": float(
                np.polyfit(orbits, series["mixed_layer"], 1)[0]
                / (CRHOS * CPS * series["mld"].mean())),
            "mixed_layer_drift_implied_by_toa_k_per_orbit": float(
                mean_toa * orbit_seconds / (CRHOS * CPS * series["mld"].mean())),
        },
        "insolation_closure": {
            "identity": "<rst - rsut> = GSOL0 / (4 sqrt(1 - e^2))",
            "expected_w_m2": float(expected_insolation),
            "reported_w_m2": reported_insolation,
            "residual_w_m2": reported_insolation - float(expected_insolation),
            "planetary_albedo": float(-series["rsut"].mean() / reported_insolation),
        },
        "executable_sha256": manifest.get("executable", {}).get("sha256"),
        "software": manifest.get("software"),
        "per_orbit": per_orbit,
    }

    args.output.mkdir(parents=True, exist_ok=True)
    # Named by the run AND the window. The window is a parameter of the
    # measurement, not a detail of it -- a storage trend taken across the
    # NLOWIO change is a different reading from one taken inside a regime -- and
    # a fixed filename would let the second overwrite the evidence for the first.
    path = args.output / f"{run_dir.name}_state_energy_{args.first}-{args.last}.json"
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({"state_energy_closure": report["state_energy_closure"],
                      "insolation_closure": report["insolation_closure"],
                      "report": str(path)}, indent=2))


if __name__ == "__main__":
    main()
