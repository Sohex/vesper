#!/usr/bin/env python3
"""Assemble Australian bore water table depths into one site table, for GW-3.

    python hydrography/scripts/build_earth_wtd_sites.py

GW-3 is the only EXTERNAL check the water table solver has. Everything else it
passes is an identity, a conservation law or a reduction, and none of those can
say the model is right about a real water table.

Fan et al. (2013)'s compilation is the usual target and its host is gone; the
route and the dead ends are in `docs/src/reference/external-data.md`. This
rebuilds the Australian part from the archive it came from, and Australia is the
part worth having: this world is dominated by arid endorheic basins, and Fan's
set is 88% United States and Canada, which are humid-to-semiarid mid-latitudes.
A calibration covering most of the wells and none of the regime tests the wrong
thing.

Input is the Australian Groundwater Explorer's per-state download, unpacked
under `--source`. Each state ships an NGIS geodatabase AND a `level_<state>.csv`
BESIDE it; the levels are not inside the geodatabase, which is documented and is
the trap that costs a download.

## The three datums are three views of one water level

`obs_point_datum` labels every reading, and the three are not interchangeable:

    DTW           depth below the GROUND SURFACE. Fan's quantity exactly.
    RSWL (mAHD)   water level ELEVATION on the Australian Height Datum.
    SWL           depth below the bore's own REFERENCE POINT, e.g. top of casing.

DTW is used as it stands and needs no elevation at all, which matters because
`LandElev` is sparse -- 10.8% of Western Australian bores carry it and none of
Tasmania's. RSWL converts as `LandElev - result`. SWL needs the casing stickup.

**The stickup is `TsRefElev - LandElev`, not `RefElev - LandElev`.** `RefElev` is
a copy of `LandElev` in this data, so using it silently applies no correction:
SWL then reads 0.570 m too deep in the median, on 91.8% of paired readings.
`TsRefElev` is the TIME SERIES reference and gives a median stickup of 0.633 m,
5th to 95th percentile 0.13 to 1.08 m, which are ordinary casing heights. It is
populated for only about 4% of bores, so SWL is dropped where it is absent
rather than assumed zero.

**Both conversions are checked against an identity that can fail.** Where one
bore reports two datums on one date, `DTW + RSWL` must equal `LandElev`, and
`SWL - stickup` must equal `DTW`. On 238,185 and 156,037 paired readings the
median residuals are +0.000 m and -0.004 m. That is what licenses the
conversions; without it they are assertions about somebody else's datum.

## What is NOT filtered, and why

**The quality flag.** Codes run A to F and are defined neither in the download's
README nor on the Explorer's metadata page, and they are per-agency: every
Western Australian reading is `quality-A` and no Australian Capital Territory or
Tasmanian reading is. Filtering on A would therefore select which STATE survives
rather than which reading, which is a selection bias wearing a quality label.

So all readings are kept and the flag is REPORTED instead: `frac_qualityA` per
site, `wtd_m_qualityA` beside `wtd_m`, and the shift between them in the report.
It comes out at a median of 0.0000 m with a 90th percentile of 0.069 m, which is
the measured answer to a question the documentation could not settle.

**Bores of unknown purpose.** `FTypeClass` is `Unknown` for 26k sites against
47k `Monitoring`, so dropping unknowns would discard a third of the set. Only
KNOWN production classes are dropped, because those carry pumping drawdown that
no model here reproduces. Fan's four columns cannot make this distinction at
all, and being able to is the reason to rebuild from the archive rather than
download her merged file.
"""

from __future__ import annotations

import argparse
import collections
import glob
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyogrio

from _paths import DATA, PROJECT_ROOT

# Purposes whose water level carries pumping drawdown. `Unknown` is NOT here:
# see the docstring. `Exploration` is not production and is kept.
PRODUCTION = {"water supply", "irrigation", "stock and domestic",
              "industrial", "dewatering", "commercial and industrial"}

# A depth outside this is a datum error or a typo, not a water table.
DEPTH_MIN_M, DEPTH_MAX_M = -50.0, 1500.0

BORE_COLS = ["HydroID", "Latitude", "Longitude", "LandElev", "TsRefElev", "FTypeClass"]


def assemble(source: Path) -> tuple[pd.DataFrame, collections.Counter]:
    frames, stats = [], collections.Counter()
    for d in sorted(glob.glob(str(source / "gdb_*"))):
        state = os.path.basename(d).replace("gdb_", "")
        levels, gdbs = glob.glob(f"{d}/level_*.csv"), glob.glob(f"{d}/*.gdb")
        if not levels or not gdbs:
            print(f"  {state:4} no level csv or geodatabase, skipped")
            continue
        lv = pd.read_csv(levels[0], usecols=["hydroid", "bore_date",
                                             "obs_point_datum", "result", "quality_flag"])
        stats["raw_readings"] += len(lv)
        bores = pyogrio.read_dataframe(gdbs[0], layer="NGIS_Bore",
                                       columns=BORE_COLS, read_geometry=False)
        lv = lv.join(bores.set_index("HydroID"), on="hydroid")

        stickup = lv.TsRefElev - lv.LandElev
        depth = pd.Series(np.nan, index=lv.index)
        m = lv.obs_point_datum == "DTW"
        depth[m] = lv.result[m]
        stats["from_dtw"] += int(m.sum())
        m = lv.obs_point_datum == "RSWL (mAHD)"
        depth[m] = lv.LandElev[m] - lv.result[m]
        stats["from_rswl"] += int(m.sum())
        m = (lv.obs_point_datum == "SWL") & stickup.notna()
        depth[m] = lv.result[m] - stickup[m]
        stats["from_swl_corrected"] += int(m.sum())
        stats["swl_dropped_no_stickup"] += int(((lv.obs_point_datum == "SWL")
                                                & stickup.isna()).sum())
        lv["depth_m"] = depth

        lv = lv[lv.depth_m.notna() & lv.Latitude.notna() & lv.Longitude.notna()]
        lv = lv[(lv.depth_m > DEPTH_MIN_M) & (lv.depth_m < DEPTH_MAX_M)]
        stats["usable_readings"] += len(lv)
        production = lv.FTypeClass.astype(str).str.strip().str.lower().isin(PRODUCTION)
        stats["production_dropped"] += int(production.sum())
        lv = lv[~production]

        lv["is_a"] = lv.quality_flag.eq("quality-A")
        site = lv.groupby("hydroid").agg(
            lat=("Latitude", "first"), lon=("Longitude", "first"),
            land_elev_m=("LandElev", "first"),
            wtd_m=("depth_m", "mean"), wtd_sd_m=("depth_m", "std"),
            n_readings=("depth_m", "size"), frac_qualityA=("is_a", "mean"),
            purpose=("FTypeClass", "first"),
            first_reading=("bore_date", "min"), last_reading=("bore_date", "max"))
        only_a = lv[lv.is_a].groupby("hydroid")["depth_m"].mean().rename("wtd_m_qualityA")
        site = site.join(only_a)
        site["state"] = state
        frames.append(site.reset_index())
        print(f"  {state:4} {len(lv):>10,} readings -> {len(site):>8,} sites")
    return pd.concat(frames, ignore_index=True), stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, required=True,
                    help="directory holding the unpacked gdb_<STATE> folders")
    ap.add_argument("--out", type=Path, default=DATA / "earth_validation" / "aus_wtd_sites.csv")
    args = ap.parse_args()

    df, stats = assemble(args.source)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    paired = df[df.wtd_m_qualityA.notna()]
    shift = paired.wtd_m_qualityA - paired.wtd_m
    report = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "source": str(args.source),
        "sites": int(len(df)),
        "readings_behind_sites": int(df.n_readings.sum()),
        "counts": dict(stats),
        "sites_by_state": df.state.value_counts().to_dict(),
        "purpose": df.purpose.astype(str).value_counts().to_dict(),
        "wtd_m": {"median": float(df.wtd_m.median()),
                  "p05": float(df.wtd_m.quantile(0.05)),
                  "p95": float(df.wtd_m.quantile(0.95)),
                  "artesian_fraction": float((df.wtd_m < 0).mean())},
        "quality_flag_sensitivity": {
            "note": "codes undefined in the download README and on the Explorer "
                    "metadata page, and per-agency, so reported rather than filtered",
            "sites_with_any_quality_a": int(len(paired)),
            "median_shift_m": float(shift.median()),
            "p90_abs_shift_m": float(shift.abs().quantile(0.9)),
            "fraction_of_sites_all_quality_a": float((df.frac_qualityA == 1).mean())},
        "sha256": hashlib.sha256(args.out.read_bytes()).hexdigest(),
        "software": {"pandas": pd.__version__, "pyogrio": pyogrio.__version__,
                     "gdal": pyogrio.__gdal_version_string__},
    }
    rp = args.out.with_suffix(".provenance.json")
    rp.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {args.out.relative_to(PROJECT_ROOT)}  ({len(df):,} sites)")
    print(f"wrote {rp.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
