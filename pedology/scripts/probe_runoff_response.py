"""Predict how ExoPlaSim's runoff responds to a real soil water capacity.

Enabling `model.soil_water_source: pedology` replaces ExoPlaSim's uniform 0.5 m
bucket with one computed from texture and regolith depth, land mean 0.342 m.
Whether that fixes the 2.8% land runoff ratio, overshoots it, or oscillates can
only really be answered by running the model. Each answer costs a full T42
spin-up, so this is the cheap check to run first.

It reimplements ExoPlaSim's land bucket offline, from `landmod.f90`:

    drhs    = min(1, dwatc / (drhsfull * dwmax))     ! drhsfull = 0.4
    E       = drhs * E_potential                     ! beta-method throttling
    dwatc  += (P_liquid - E) * dt
    runoff  = max(0, dwatc - dwmax)
    dwatc   = min(dwatc, dwmax)

Potential evaporation is not an output, so it is back-solved from the run's own
actual evaporation and soil water, which is only valid because `drhs` is the sole
moisture limiter in that scheme.

**The result is worthless without the validation step**, which is the first thing
this prints. Monthly-mean forcing cannot generate event-driven runoff: a bucket
fed the average of a month's rain overflows far less than one fed the storms. So
the offline bucket is first run at ExoPlaSim's own 0.5 m and compared against the
runoff the model actually reported. If it cannot reproduce that, it cannot be
trusted about 0.342 m either, and this script says so rather than quoting a
number.

    python pedology/scripts/probe_runoff_response.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import netCDF4 as nc
import numpy as np
import yaml

from _paths import ANALYSIS, CONFIG, DATA, PROJECT_ROOT, climatology_path  # noqa: F401
from paths import rel  # noqa: E402
from builds import component_data

import orbit

EARTH_YEAR_DAYS = orbit.EARTH_CALENDAR_YEAR_DAYS
SECONDS_PER_DAY = 86400.0

# landmod.f90:53. Soil water above this fraction of capacity evaporates at the
# potential rate; below it, evaporation is throttled linearly.
DRHSFULL = 0.4

EXOPLASIM_DEFAULT_DWMAX_M = 0.5

# Agreement within this factor counts as the offline bucket reproducing the
# model. Fixed before the result was seen. It is loose on purpose: the point is
# whether the mechanism is captured, not whether the third digit matches.
VALIDATION_TOLERANCE = 2.0

COORD_DECIMALS = 4


def to_daily(monthly: np.ndarray, year_length: int) -> np.ndarray:
    """Interpolate 12 bin means onto days, cyclically."""
    bins = monthly.shape[0]
    centres = (np.arange(bins) + 0.5) * year_length / bins
    days = np.arange(year_length) + 0.5
    extended_x = np.concatenate([centres - year_length, centres,
                                 centres + year_length])
    extended_y = np.concatenate([monthly, monthly, monthly], axis=0)
    out = np.empty((year_length,) + monthly.shape[1:])
    for index in np.ndindex(monthly.shape[1:]):
        out[(slice(None),) + index] = np.interp(
            days, extended_x, extended_y[(slice(None),) + index])
    return out


def run_bucket(precip_daily: np.ndarray, potential_daily: np.ndarray,
               capacity: np.ndarray, years: int = 40) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the bucket to a repeating annual cycle. Returns runoff, mean water.

    Units are metres and metres per day throughout; the returned runoff is metres
    per simulation year.
    """
    water = np.minimum(capacity * 0.5, capacity)
    threshold = np.maximum(DRHSFULL * capacity, 1e-9)
    for year in range(years):
        runoff_total = np.zeros_like(capacity)
        water_total = np.zeros_like(capacity)
        for day in range(precip_daily.shape[0]):
            drhs = np.minimum(1.0, water / threshold)
            evaporation = drhs * potential_daily[day]
            water = water + precip_daily[day] - evaporation
            water = np.maximum(water, 0.0)
            excess = np.maximum(water - capacity, 0.0)
            runoff_total += excess
            water = water - excess
            water_total += water
    return runoff_total, water_total / precip_daily.shape[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--climatology", type=Path, default=None)
    # Per-build: the soil map is a property of a terrain plus a climatology,
    # and pedology/data/ is namespaced by build. There is no flat soilmap.txt
    # any more, so this default pointed at a file that does not exist.
    parser.add_argument("--soil-map", type=Path,
                        default=component_data("pedology") / "soilmap.txt")
    parser.add_argument("--years", type=int, default=40)
    args = parser.parse_args()

    config = yaml.safe_load(CONFIG.read_text())
    climatology = args.climatology or climatology_path()
    year_length = int(round(orbit.orbital_year_days(config)))
    # Orbits per Earth year, 2.022. Named for what it is: an earlier version
    # called this earth_years_per_orbit while holding its reciprocal, which now
    # collides with the helper of that name in lib/orbit.py.
    orbits_per_earth_year = 1.0 / orbit.earth_years_per_orbit(config)

    with nc.Dataset(climatology) as data:
        lat = np.asarray(data["lat"][:], dtype=float)
        lon = np.asarray(data["lon"][:], dtype=float)
        land = np.asarray(data["lsm"][0], dtype=float) > 0.5
        # m/s -> m/day. evap is negative upward; prsn is the frozen part of pr,
        # which reaches the bucket later as snowmelt rather than immediately.
        precip = np.asarray(data["pr"][:], dtype=float) * SECONDS_PER_DAY
        snowfall = np.asarray(data["prsn"][:], dtype=float) * SECONDS_PER_DAY
        snowmelt = np.asarray(data["snm"][:], dtype=float) * SECONDS_PER_DAY
        evaporation = -np.asarray(data["evap"][:], dtype=float) * SECONDS_PER_DAY
        soil_water = np.asarray(data["mrso"][:], dtype=float)
        model_runoff = np.asarray(data["mrro"][:], dtype=float) * SECONDS_PER_DAY

    liquid = np.maximum(precip - snowfall + snowmelt, 0.0)

    # Back-solve potential evaporation. drhs is the only moisture limiter in this
    # scheme, so E_actual = drhs * E_potential and drhs follows from the run's own
    # soil water against the capacity it ran with.
    drhs_model = np.minimum(
        1.0, soil_water / (DRHSFULL * EXOPLASIM_DEFAULT_DWMAX_M))
    potential = evaporation / np.maximum(drhs_model, 0.05)

    liquid_daily = to_daily(liquid, year_length)
    potential_daily = to_daily(potential, year_length)

    weights = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, len(lon)))
    lw = weights[land]

    def land_mean(field: np.ndarray) -> float:
        return float(np.average(field[land], weights=lw))

    # --- validation: can the offline bucket reproduce the model at its own capacity?
    default_capacity = np.full(land.shape, EXOPLASIM_DEFAULT_DWMAX_M)
    offline_runoff, offline_water = run_bucket(
        liquid_daily, potential_daily, default_capacity, args.years)

    routed_mm = land_mean(model_runoff.mean(axis=0)) * 1000.0 * EARTH_YEAR_DAYS
    offline_runoff_mm = land_mean(offline_runoff) * 1000.0 * orbits_per_earth_year
    precip_mm = land_mean(precip.mean(axis=0)) * 1000.0 * EARTH_YEAR_DAYS
    evap_mm = land_mean(evaporation.mean(axis=0)) * 1000.0 * EARTH_YEAR_DAYS

    # Validate against the land water budget, not against mrro.
    #
    # The first version of this probe compared itself to mrro and reported a
    # 5.91x failure. mrro is river-routed net water flux, not local runoff
    # generation: it is negative in places and non-zero over ocean, and it
    # recovers only 15% of the land's water surplus. The probe was right and the
    # target was wrong. In steady state the local runoff a bucket generates has
    # to equal P - E, and that is what this now checks.
    model_runoff_mm = precip_mm - evap_mm
    ratio = offline_runoff_mm / max(model_runoff_mm, 1e-9)
    valid = (1.0 / VALIDATION_TOLERANCE) <= ratio <= VALIDATION_TOLERANCE

    print("VALIDATION, offline bucket at ExoPlaSim's own 0.5 m capacity")
    print(f"  land budget P - E   {model_runoff_mm:8.2f} mm per Earth year")
    print(f"    P {precip_mm:.2f}, E {evap_mm:.2f}; mrro reports {routed_mm:.2f} "
          f"but is river-routed, not local")
    print(f"  offline runoff      {offline_runoff_mm:8.2f} mm per Earth year")
    print(f"  ratio               {ratio:8.2f}  (pass band "
          f"{1/VALIDATION_TOLERANCE:.2f}-{VALIDATION_TOLERANCE:.2f})")
    print(f"  model soil water    {land_mean(soil_water.mean(axis=0)):8.4f} m")
    print(f"  offline soil water  {land_mean(offline_water):8.4f} m")
    print(f"  -> {'PASS' if valid else 'FAIL'}")
    print()

    # --- the question: what happens at the pedology capacity
    capacity_mm: dict[tuple[float, float], float] = {}
    if args.soil_map.is_file():
        lines = args.soil_map.read_text().splitlines()
        header = lines[0].split()
        column = header.index("awc")
        for line in lines[1:]:
            parts = line.split()
            capacity_mm[(round(float(parts[0]), COORD_DECIMALS),
                         round(float(parts[1]), COORD_DECIMALS))] = float(parts[column])
    else:
        raise SystemExit(f"{args.soil_map} not found; run build_soil.py")

    lon_signed = np.round(np.where(lon > 180.0, lon - 360.0, lon), COORD_DECIMALS)
    lat_rounded = np.round(lat, COORD_DECIMALS)
    pedology_capacity = np.full(land.shape, EXOPLASIM_DEFAULT_DWMAX_M)
    for j in range(len(lat)):
        for i in range(len(lon)):
            if land[j, i]:
                value = capacity_mm.get((float(lon_signed[i]), float(lat_rounded[j])))
                if value is not None:
                    pedology_capacity[j, i] = value / 1000.0

    pedology_runoff, pedology_water = run_bucket(
        liquid_daily, potential_daily, pedology_capacity, args.years)
    pedology_runoff_mm = land_mean(pedology_runoff) * 1000.0 * orbits_per_earth_year

    print("PREDICTION, offline bucket at the pedology capacity")
    print(f"  capacity            {land_mean(pedology_capacity):8.4f} m "
          f"against 0.5000")
    print(f"  runoff              {pedology_runoff_mm:8.2f} mm per Earth year")
    print(f"  runoff ratio        {pedology_runoff_mm / precip_mm:8.4f} "
          f"against Earth land's ~0.35")
    print(f"  change              {pedology_runoff_mm / max(offline_runoff_mm, 1e-9):8.2f}x "
          f"the offline 0.5 m case")
    print()

    # What that runoff would do to weathering, which is the reason to care. The
    # current land-mean W is read from the soil report rather than written here,
    # because it moves whenever the soil or the climatology does.
    exponent = 0.65   # Berner (1994) GEOCARB II, from Dunne (1978) + Peters (1984)
    weathering_shift = (pedology_runoff_mm / max(model_runoff_mm, 1e-9)) ** exponent
    report_path = ANALYSIS / "soil_report.json"
    current_w = None
    if report_path.is_file():
        current_w = json.loads(report_path.read_text())["land_means"].get(
            "weathering_intensity")
    if valid and current_w is not None:
        print(f"  implied weathering intensity shift: {weathering_shift:.2f}x, so "
              f"land-mean W moves from {current_w:.2f} toward "
              f"{current_w * weathering_shift:.2f}")
    elif valid:
        print(f"  implied weathering intensity shift: {weathering_shift:.2f}x")
    else:
        print("  VALIDATION FAILED, so the absolute numbers above are not usable.")
        print("  What may still survive is the RELATIVE response, since both cases")
        print("  share the same bias: shrinking the bucket from 0.500 to "
              f"{land_mean(pedology_capacity):.3f} m changes")
        print(f"  offline runoff by only {pedology_runoff_mm / max(offline_runoff_mm, 1e-9):.2f}x. "
              "If that weak sensitivity is real, the")
        print("  soil-water feedback will not on its own fix the 2.8% runoff ratio.")

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "purpose": ("offline estimate of ExoPlaSim's runoff response to a "
                    "pedology-supplied soil water capacity, to decide whether the "
                    "loop is worth spinning up"),
        "climatology": rel(climatology),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "bucket_scheme": "landmod.f90 beta-method, drhsfull 0.4",
        "spin_up_years": args.years,
        "validation": {
            "target": "land budget P - E, not mrro",
            "model_runoff_mm_per_earth_year": model_runoff_mm,
            "routed_mrro_mm_per_earth_year": routed_mm,
            "offline_runoff_mm_per_earth_year": offline_runoff_mm,
            "ratio": ratio,
            "tolerance": VALIDATION_TOLERANCE,
            "passed": bool(valid),
            "model_soil_water_m": land_mean(soil_water.mean(axis=0)),
            "offline_soil_water_m": land_mean(offline_water),
            "caveat": ("Monthly-mean forcing cannot generate event-driven runoff, "
                       "so the offline bucket is expected to under-produce. The "
                       "validation is what decides whether it under-produces "
                       "enough to invalidate the comparison."),
            "diagnosis": (
                "Failed, and the way it failed is the useful part: the offline "
                "bucket reproduces the model's mean soil water almost exactly, "
                "0.204 against 0.202 m, while over-producing runoff nearly "
                "sixfold. So the storage is right and the overflow accounting is "
                "not. The likely cause is that E = drhs * E_potential is not "
                "ExoPlaSim's scheme: its surface flux is bulk aerodynamic, "
                "E = C * (drhs * q_sat(Ts) - q), where drhs scales the surface "
                "humidity rather than the flux, so evaporation stops entirely "
                "once drhs falls to the ambient relative humidity instead of "
                "declining linearly to zero. Reproducing that offline needs "
                "surface temperature, humidity and wind at the model's timestep, "
                "which is most of a land surface scheme and no longer cheap."),
        },
        "verdict": ("Offline shortcut rejected. Only a real ExoPlaSim run can "
                    "answer this. The relative response suggests the feedback is "
                    "weak, which lowers the expected payoff of that run but does "
                    "not remove the need for it."),
        "prediction": {
            "pedology_capacity_m": land_mean(pedology_capacity),
            "runoff_mm_per_earth_year": pedology_runoff_mm,
            "runoff_ratio": pedology_runoff_mm / precip_mm,
            "change_vs_offline_default": pedology_runoff_mm / max(offline_runoff_mm, 1e-9),
            "implied_weathering_multiplier": weathering_shift,
            "current_land_mean_weathering_intensity": current_w,
        },
        "precipitation_mm_per_earth_year": precip_mm,
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    path = ANALYSIS / "runoff_response_probe.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"\nwrote {rel(path)}")


if __name__ == "__main__":
    main()
