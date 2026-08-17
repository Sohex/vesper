#!/usr/bin/env python3
"""Test the lake solver's equilibrium relation against Earth's terminal lakes.

    python hydrography/scripts/validate_lake_solver.py --data-dir <downloads>

HYD-4. The design, the thresholds and the sources of scatter were registered in
`hydrography/notes/lake-solver-validation.md` before any data was fetched; this
script is the run, and the note carries the result.

## What is under test

A terminal lake grows until net evaporation over its surface consumes what its
catchment delivers, which is the relation `lake_balance.solve` iterates to a
fixed point:

    R * (A_catch - A_lake) = (E_lake - P_lake) * A_lake

    A_lake / A_catch = R / (R + E_lake - P_lake)

The prediction needs no hypsometry, so it is scale free and can be checked across
four orders of magnitude of lake size. The hypsometric inversion and the overflow
cascade are NOT tested here; Earth basins have their own hypsometry, and testing
it would test the data.

## Why these sources

The first attempt at this test failed on the literature, and the reason is in the
note: reviews state a terminal lake's balance per unit LAKE area, so it closes on
itself and cannot predict an area. The supply term has to be a depth over the
CATCHMENT, measured independently of the lake.

HydroLAKES carries `Lake_area`, `Wshd_area` and `Dis_avg` per lake under one set
of definitions, so `R = Dis_avg / (Wshd_area - Lake_area)` is exactly that. Its
discharge is modelled rather than gauged, from WaterGAP, so it is independent of
the equation under test but not an observation; what that costs is bounded by the
literature basins in the note.

GLEV supplies lake evaporation for the same 1.42 million lakes, keyed by
`Hylak_id`, so evaporation and area come from datasets built on each other rather
than reconciled by hand.

## Selection, fixed before the result was seen

- HydroBASINS level 5, `ENDO > 0`: any basin belonging to an endorheic system.
- Grouped by `MAIN_BAS`, the identifier of the whole system, taking the lake with
  the largest `Wshd_area`, which is the most downstream and therefore terminal.
- `Lake_type == 1`, a natural lake. Reservoirs and regulated lakes are excluded
  because a dam sets their area, not their evaporation.
- `Lake_area >= 10 km2`, below which the area is a few HydroLAKES polygons.
- `Dis_avg > 0` and `Wshd_area > Lake_area`, without which the relation has no
  supply term.

**These are an amendment, made after seeing which lakes the registered rule
selected and before computing any prediction.** The registered rule was
`ENDO == 1`, an endorheic *sink* basin, grouped per basin and with no lake-type
filter. It returned 54 lakes, four of the five largest were reservoirs, and it
missed every major terminal lake on Earth: Balkhash, Chad, Turkana, Urmia, Van,
Eyre and Qinghai all sit in basins HydroBASINS flags `ENDO == 2`, part of an
endorheic system rather than its sink. The amendment is recorded here and in the
note rather than quietly applied, and no predicted area had been computed when it
was made.

The Caspian Sea is lost either way: HydroBASINS classifies its shore basins as
`ENDO == 0`, treating it as a sea rather than a terminal lake, so the largest
endorheic lake on Earth is outside this test by the source's own taxonomy.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import ANALYSIS, PROJECT_ROOT  # noqa: E402

# Registered in the note before any data was fetched.
PASS_MEDIAN_ABS_LOG10 = 0.30
MARGINAL_MEDIAN_ABS_LOG10 = 0.60
MIN_LAKE_KM2 = 10.0

HYBAS_REGIONS = ("af", "ar", "as", "au", "eu", "gr", "na", "sa", "si")
OPEN_METEO = "https://archive-api.open-meteo.com/v1/archive"


def endorheic_sinks(data_dir: Path):
    """Level-5 basins flagged as endorheic sinks, all continents."""
    from pyogrio import read_dataframe
    frames = []
    for region in HYBAS_REGIONS:
        path = data_dir / f"hybas_{region}_lev05_v1c.shp"
        if not path.is_file():
            print(f"  missing {path.name}, skipping region")
            continue
        gdf = read_dataframe(path, columns=["HYBAS_ID", "MAIN_BAS", "SUB_AREA", "ENDO"])
        frames.append(gdf[gdf["ENDO"] > 0])
    if not frames:
        raise SystemExit("no HydroBASINS files found")
    import geopandas as gpd
    return gpd.GeoDataFrame(pd.concat(frames, ignore_index=True), crs=frames[0].crs)


def terminal_lakes(data_dir: Path, sinks):
    """HydroLAKES entries that are the most downstream lake of a sink basin."""
    import geopandas as gpd
    from pyogrio import read_dataframe
    cols = ["Hylak_id", "Lake_name", "Country", "Lake_type", "Lake_area",
            "Wshd_area", "Dis_avg", "Pour_long", "Pour_lat", "Elevation"]
    lakes = read_dataframe(data_dir / "HydroLAKES_polys_v10_shp"
                           / "HydroLAKES_polys_v10.shp",
                           columns=cols, read_geometry=False)
    lakes = lakes[(lakes["Lake_area"] >= MIN_LAKE_KM2)
                  & (lakes["Lake_type"] == 1)
                  & (lakes["Dis_avg"] > 0)
                  & (lakes["Wshd_area"] > lakes["Lake_area"])]
    pts = gpd.GeoDataFrame(
        lakes,
        geometry=gpd.points_from_xy(lakes["Pour_long"], lakes["Pour_lat"]),
        crs="EPSG:4326")
    joined = gpd.sjoin(pts, sinks[["HYBAS_ID", "MAIN_BAS", "SUB_AREA", "geometry"]],
                       how="inner", predicate="within")
    idx = joined.groupby("MAIN_BAS")["Wshd_area"].idxmax()
    return joined.loc[idx].drop(columns="geometry").reset_index(drop=True)


def glev_evaporation(path: Path, ids) -> pd.Series:
    """Annual mean lake evaporation, mm/yr, from GLEV monthly rates."""
    wanted = set(int(i) for i in ids)
    keep = []
    for chunk in pd.read_csv(path, chunksize=200_000):
        idcol = chunk.columns[0]
        sel = chunk[chunk[idcol].isin(wanted)]
        if len(sel):
            keep.append(sel)
    if not keep:
        raise SystemExit("no GLEV rows matched the selected lakes")
    df = pd.concat(keep, ignore_index=True)
    idcol = df.columns[0]
    months = [c for c in df.columns if c != idcol]
    # Monthly rates in mm/month; the annual total is their sum per year, and the
    # columns are a monthly time series, so the mean month times twelve is the
    # long-term annual rate.
    annual = df[months].apply(pd.to_numeric, errors="coerce").mean(axis=1) * 12.0
    return pd.Series(annual.values, index=df[idcol].astype(int).values)


def lake_precipitation(lat: float, lon: float, cache: Path,
                       tries: int = 4) -> float | None:
    """Long-term annual precipitation over the lake, mm/yr, from ERA5.

    Cached and retried for the reason `pedology/scripts/validate_against_earth.py`
    gives at length: a validation whose answer depends on which points the API
    rate-limited that minute is worse than none, because it looks like a
    measurement. A thirty-year daily series is a heavy request and the first pass
    over these 146 lakes lost 117 of them to 429s.
    """
    cache.mkdir(parents=True, exist_ok=True)
    key = cache / f"pr_{lat:.3f}_{lon:.3f}.json"
    if key.is_file():
        payload = json.loads(key.read_text(encoding="utf-8"))
    else:
        url = (f"{OPEN_METEO}?latitude={lat:.4f}&longitude={lon:.4f}"
               "&start_date=1991-01-01&end_date=2020-12-31"
               "&daily=precipitation_sum&timezone=UTC")
        payload = None
        for attempt in range(tries):
            try:
                with urllib.request.urlopen(url, timeout=120) as response:
                    payload = json.load(response)
                break
            except Exception as exc:                   # noqa: BLE001
                if attempt == tries - 1:
                    print(f"    precipitation fetch failed: {type(exc).__name__}")
                    return None
                time.sleep(15 * (attempt + 1))
        if payload is None:
            return None
        key.write_text(json.dumps(payload), encoding="utf-8")
    daily = (payload.get("daily") or {}).get("precipitation_sum")
    if not daily:
        return None
    values = np.array([v for v in daily if v is not None], dtype=float)
    return float(values.sum() / (len(values) / 365.25))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, required=True,
                    help="directory holding the HydroBASINS, HydroLAKES and "
                         "GLEV downloads")
    ap.add_argument("--glev", type=Path, default=None)
    ap.add_argument("--output", type=Path,
                    default=ANALYSIS / "lake_solver_validation.json")
    ap.add_argument("--max-lakes", type=int, default=0,
                    help="0 for all; a positive value truncates the set AFTER "
                         "selection, which is reported")
    args = ap.parse_args()
    glev_path = args.glev or (args.data_dir / "glev_evaporation_rate.csv")

    print("selecting endorheic sink basins")
    sinks = endorheic_sinks(args.data_dir)
    print(f"  {len(sinks)} endorheic basin units, "
          f"{sinks['MAIN_BAS'].nunique()} systems")

    print("selecting the terminal lake of each")
    lakes = terminal_lakes(args.data_dir, sinks)
    print(f"  {len(lakes)} terminal lakes")
    truncated = 0
    if args.max_lakes and len(lakes) > args.max_lakes:
        lakes = lakes.nlargest(args.max_lakes, "Lake_area").reset_index(drop=True)
        truncated = args.max_lakes
        print(f"  truncated to the {args.max_lakes} largest, which is reported")

    print("reading GLEV evaporation")
    evap = glev_evaporation(glev_path, lakes["Hylak_id"])
    lakes["evap_mm_yr"] = lakes["Hylak_id"].astype(int).map(evap)

    print("fetching lake precipitation")
    cache = PROJECT_ROOT / "hydrography" / "analysis" / "cache_era5"
    lakes["precip_mm_yr"] = [
        lake_precipitation(row.Pour_lat, row.Pour_long, cache)
        for row in lakes.itertuples()
    ]

    have = lakes.dropna(subset=["evap_mm_yr", "precip_mm_yr"]).copy()
    # Supply as a depth over the dry catchment, which is the whole point.
    seconds = 365.25 * 86400.0
    have["runoff_mm_yr"] = (have["Dis_avg"] * seconds
                            / ((have["Wshd_area"] - have["Lake_area"]) * 1e6)
                            * 1000.0)
    have["net_evap_mm_yr"] = have["evap_mm_yr"] - have["precip_mm_yr"]
    usable = have[have["net_evap_mm_yr"] > 0].copy()
    have["predicted_ratio"] = (have["runoff_mm_yr"]
                               / (have["runoff_mm_yr"] + have["net_evap_mm_yr"]))
    have["predicted_km2"] = have["predicted_ratio"] * have["Wshd_area"]
    have["observed_ratio"] = have["Lake_area"] / have["Wshd_area"]
    have["log10_pred_over_obs"] = np.log10(have["predicted_km2"]
                                           / have["Lake_area"])

    scored = have[np.isfinite(have["log10_pred_over_obs"])
                  & (have["net_evap_mm_yr"] > 0)].copy()
    median_abs = float(scored["log10_pred_over_obs"].abs().median())
    verdict = ("pass" if median_abs < PASS_MEDIAN_ABS_LOG10 else
               "marginal" if median_abs < MARGINAL_MEDIAN_ABS_LOG10 else "fail")

    def corr(a, b):
        m = np.isfinite(a) & np.isfinite(b)
        return float(np.corrcoef(a[m], b[m])[0, 1]) if m.sum() > 2 else float("nan")

    err = scored["log10_pred_over_obs"].to_numpy()
    correlations = {
        "with_log_lake_area": corr(err, np.log10(scored["Lake_area"].to_numpy())),
        "with_log_catchment_area": corr(err, np.log10(scored["Wshd_area"].to_numpy())),
        "with_log_runoff": corr(err, np.log10(scored["runoff_mm_yr"].to_numpy())),
        "with_net_evaporation": corr(err, scored["net_evap_mm_yr"].to_numpy()),
    }

    report = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "test": "A_lake / A_catch = R / (R + E_lake - P_lake)",
        "registered_in": "hydrography/notes/lake-solver-validation.md",
        "thresholds": {
            "pass_median_abs_log10": PASS_MEDIAN_ABS_LOG10,
            "marginal_median_abs_log10": MARGINAL_MEDIAN_ABS_LOG10,
            "min_lake_km2": MIN_LAKE_KM2,
        },
        "sources": {
            "areas_and_discharge": "HydroLAKES v1.0, Messager et al. (2016); "
                                   "Dis_avg is modelled by WaterGAP, not gauged",
            "endorheic_flag": "HydroBASINS level 5, ENDO > 0, grouped by MAIN_BAS",
            "lake_evaporation": "GLEV, Zhao et al. (2022), keyed by Hylak_id",
            "lake_precipitation": "ERA5 via the Open-Meteo archive, 1991-2020 "
                                  "at the lake pour point",
        },
        "counts": {
            "sink_basins": int(len(sinks)),
            "terminal_lakes_selected": int(len(lakes)),
            "truncated_to": truncated,
            "with_evaporation_and_precipitation": int(len(have)),
            "scored": int(len(scored)),
            "dropped_net_evaporation_not_positive":
                int((have["net_evap_mm_yr"] <= 0).sum()),
        },
        "result": {
            "median_abs_log10_pred_over_obs": median_abs,
            "verdict": verdict,
            "median_log10_pred_over_obs":
                float(scored["log10_pred_over_obs"].median()),
            "fraction_within_factor_2":
                float((scored["log10_pred_over_obs"].abs() < np.log10(2)).mean()),
            "fraction_within_factor_10":
                float((scored["log10_pred_over_obs"].abs() < 1.0).mean()),
            "error_correlations": correlations,
        },
        "lakes": json.loads(
            scored.sort_values("Lake_area", ascending=False)
            [["Hylak_id", "Lake_name", "Country", "Lake_area", "Wshd_area",
              "Dis_avg", "runoff_mm_yr", "evap_mm_yr", "precip_mm_yr",
              "predicted_km2", "log10_pred_over_obs"]]
            .round(4).to_json(orient="records")),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"\n  scored {len(scored)} lakes")
    print(f"  median |log10 predicted/observed| {median_abs:.3f}  -> {verdict}")
    print(f"  within a factor of 2: {100*report['result']['fraction_within_factor_2']:.0f}%")
    print(f"  bias, median log10: {report['result']['median_log10_pred_over_obs']:+.3f}")
    for k, v in correlations.items():
        print(f"  error {k}: {v:+.3f}")
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
