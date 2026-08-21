#!/usr/bin/env python3
"""GW-22: assemble United States bore water table depths into one site table.

    python hydrography/scripts/build_us_wtd_sites.py --stage all

Worldbuilding. Vesper is an invented planet with no bores; this fetches real
measurements from a real one, because the groundwater solver's only external
check is Earth and a model nothing can contradict is not a model.

The companion to `build_earth_wtd_sites.py`, which did Australia.
`docs/src/reference/external-data.md` says take Australia FIRST, being the arid
analogue for a world of endorheic basins, and the United States second for well
count; the second leg is this. GW-21 concluded the limit on the Australian score
is the INPUT FIELDS rather than the formulation or the mesh, and that conclusion
rests on one region, so a second one is what tests it.

## Why the US set can be better than the Australian one, not merely bigger

Australia's confound is the Great Artesian Basin: a bore screened well below the
water table measures a POTENTIOMETRIC HEAD in a confined aquifer, a different
surface, and `aus_wtd_sites.csv` can only guard that with a `bore_depth_m`
heuristic. The USGS `monitoring-locations` collection carries
`aquifer_type_code` directly, so confined and unconfined are SEPARABLE here
rather than inferred, along with `well_constructed_depth` for the same
consistency check Australia uses.

Parameter code 72019 is depth to water level in feet below land surface, which
is the quantity wanted rather than a head needing a datum conversion.

`numberMatched` is absent from an ordinary response, so paging follows the
`next` link until a page comes back short. The total IS obtainable, with
`resulttype=hits` -- lower-case t, because `resultType` 400s -- and that is
worth one request at the start so progress reads as a fraction rather than a
number climbing toward nothing.
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
import time
from pathlib import Path

import os

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "hydrography" / "data" / "earth_validation"
CACHE = ROOT / "hydrography" / "data" / "earth_validation_cache"
BASE = "https://api.waterdata.usgs.gov/ogcapi/v0/collections"
BBOX = "-125,24,-66,50"          # conterminous US, matching the DEM window
PAGE = 10000                     # 5.6 s/page; 50000 is 3x worse per record
PARAM_WTD = "72019"              # depth to water level, ft below land surface

# From the API's own `aquifer-types` collection, not assumed. This is the whole
# reason the US leg can settle what Australia could not: a confined aquifer's
# bore measures a POTENTIOMETRIC HEAD, a different surface from the water table
# the solver computes, and here that is a label rather than an inference from
# bore depth.
AQUIFER_CONFINEMENT = {"U": "unconfined",   # Unconfined single aquifer
                       "N": "unconfined",   # Unconfined multiple aquifer
                       "C": "confined",     # Confined single aquifer
                       "M": "confined",     # Confined multiple aquifers
                       "X": "mixed"}        # Mixed, confined and unconfined

SITE_KEEP = ["id", "agency_code", "state_name", "site_type_code", "altitude",
             "vertical_datum", "well_constructed_depth", "hole_constructed_depth",
             "aquifer_code", "aquifer_type_code", "national_aquifer_code"]
MEAS_KEEP = ["monitoring_location_id", "value", "unit_of_measure", "time",
             "approval_status", "qualifier", "vertical_datum"]


KEY_FILE = Path.home() / ".usgs_key"


def api_key() -> str:
    """`USGS_API_KEY`, else `~/.usgs_key`, else keyless.

    The anonymous budget is 4000 requests and paging a whole collection exhausts
    it; `docs/src/reference/external-data.md` has the measured 429. The key is a
    CREDENTIAL: it is read at call time, never written to a provenance file, and
    never printed.
    """
    k = os.environ.get("USGS_API_KEY", "").strip()
    if not k and KEY_FILE.is_file():
        k = KEY_FILE.read_text().strip()
    return k


def api_headers() -> dict:
    k = api_key()
    return {"X-Api-Key": k} if k else {}


def _get(session, url, params, tries=9):
    """One request, with 429 handled as the API asks rather than as an error.

    The first version treated 429 like any exception and gave up after five
    short backoffs, losing the whole fetch. Rate limiting is not a failure, it
    is the server pacing us; honour `Retry-After` when it is sent and back off
    hard when it is not.
    """
    for att in range(tries):
        try:
            r = session.get(url, params=params, timeout=300,
                            headers=api_headers())
            if r.status_code == 429:
                wait = float(r.headers.get("Retry-After", 0) or 0) or min(30 * 2 ** att, 900)
                print(f"    429, waiting {wait:.0f}s", flush=True)
                time.sleep(wait)
                continue
            r.raise_for_status()
            return r.json()
        except Exception as e:
            if att == tries - 1:
                raise
            w = min(10 * 2 ** att, 300)
            print(f"    retry in {w}s ({type(e).__name__}: {str(e)[:80]})", flush=True)
            time.sleep(w)
    raise RuntimeError("exhausted retries")


def total_hits(session, collection: str, params: dict) -> int | None:
    """Record count for a query, or None. `resulttype=hits` is the spelling."""
    try:
        d = _get(session, f"{BASE}/{collection}/items",
                 dict(params, limit=1, resulttype="hits", f="json"), tries=3)
        return d.get("numberMatched")
    except Exception:
        return None


def paged_to_chunks(collection: str, params: dict, chunkdir: Path, keep: list[str],
                    extra=None, pace: float = 2.0) -> int:
    """Page the collection into numbered chunk files, resuming where it stopped.

    `docs/src/reference/large-data.md`: write results as they are produced and
    make re-running resume rather than restart. The state file records the query
    it belongs to and is refused when that does not match, because a checkpoint
    that has lost track of its input silently skips work never done.
    """
    chunkdir.mkdir(parents=True, exist_ok=True)
    state_p = chunkdir / "state.json"
    ident = {"collection": collection, "params": params, "keep": keep}
    url = f"{BASE}/{collection}/items"
    q = dict(params, limit=PAGE, f="json")
    state = {}
    if state_p.is_file():
        state = json.loads(state_p.read_text())
        if state.get("ident") != ident:
            raise SystemExit(
                f"{state_p} belongs to a different query. Delete {chunkdir} to refetch.")
        if state.get("done"):
            print(f"    resume: already complete, {state['records']:,} records")
            return state["records"]
        if state.get("next"):
            url, q = state["next"], {}
            print(f"    resume: {state['records']:,} records, "
                  f"chunk {state['chunk']}", flush=True)
    n = state.get("records", 0)
    chunk = state.get("chunk", 0)
    s = requests.Session()
    total = state.get("total") or total_hits(s, collection, params)
    if total:
        print(f"    {total:,} records to fetch, {-(-total // PAGE)} chunks", flush=True)
    t0 = time.time()
    while True:
        d = _get(s, url, q)
        feats = d.get("features") or []
        rows = []
        for f in feats:
            pr = f.get("properties") or {}
            rec = {k: pr.get(k) for k in keep}
            if extra:
                extra(f, pr, rec)
            for k, v in list(rec.items()):
                if isinstance(v, list):
                    rec[k] = "|".join(str(x) for x in v)
            rows.append(rec)
        if rows:
            pd.DataFrame(rows).to_csv(chunkdir / f"chunk_{chunk:05d}.csv", index=False)
            chunk += 1
        n += len(feats)
        nxt = next((l["href"] for l in d.get("links", []) if l.get("rel") == "next"), None)
        done = len(feats) < PAGE or not nxt
        state_p.write_text(json.dumps(
            {"ident": ident, "next": nxt, "records": n, "chunk": chunk,
             "done": done, "total": total}))
        pct = f"  {n/total:6.1%}" if total else ""
        print(f"    {n:,} records{pct}  {time.time()-t0:.0f}s", flush=True)
        if done:
            return n
        url, q = nxt, {}
        time.sleep(pace)


def combine(chunkdir: Path, out: Path) -> Path:
    parts = sorted(chunkdir.glob("chunk_*.csv"))
    if not parts:
        raise SystemExit(f"no chunks in {chunkdir}")
    pd.concat((pd.read_csv(p, low_memory=False) for p in parts),
              ignore_index=True).to_csv(out, index=False)
    return out


def stage_sites(quiet: bool) -> Path:
    out = CACHE / "us_gw_sites_raw.csv"
    if out.exists():
        print(f"  sites: cached ({sum(1 for _ in open(out))-1:,} rows)")
        return out
    def geom(f, pr, rec):
        g = (f.get("geometry") or {}).get("coordinates") or [None, None]
        rec["id"] = f.get("id") or pr.get("id")
        rec["lon"], rec["lat"] = g[0], g[1]
    n = paged_to_chunks("monitoring-locations", {"site_type_code": "GW", "bbox": BBOX},
                        CACHE / "us_sites_chunks", SITE_KEEP, extra=geom)
    combine(CACHE / "us_sites_chunks", out)
    print(f"  sites: {n:,} groundwater monitoring locations -> {out}")
    return out


def stage_levels(quiet: bool, full_history: bool = False) -> Path:
    """One reading per site by default, every reading with --full-history.

    The Australian table means over a site's readings. This means over one, or
    over all of them at about twenty times the paging, and the difference is
    declared rather than hidden: `latest-field-measurements` is the most recent
    reading that EXISTS, which for many sites is decades old, so it carries a
    temporal spread the Australian means do not. For a steady-state model
    scored on cell means that is a second-order concern against the 1.4 million
    sites it makes reachable, but it is the thing to revisit if a US and
    Australian score ever disagree in a way that matters.
    """
    collection = "field-measurements" if full_history else "latest-field-measurements"
    out = CACHE / (f"us_gw_levels_raw{'_full' if full_history else ''}.csv")
    if out.exists():
        print(f"  levels: cached ({collection})")
        return out
    print(f"  collection: {collection}", flush=True)
    cd = CACHE / (f"us_levels_chunks{'_full' if full_history else ''}")
    n = paged_to_chunks(collection, {"parameter_code": PARAM_WTD, "bbox": BBOX},
                        cd, MEAS_KEEP)
    combine(cd, out)
    print(f"  levels: {n:,} measurements -> {out}")
    return out


def stage_assemble(quiet: bool) -> None:
    sites = pd.read_csv(CACHE / "us_gw_sites_raw.csv", low_memory=False)
    lvp = CACHE / "us_gw_levels_raw_full.csv"
    if not lvp.is_file():
        lvp = CACHE / "us_gw_levels_raw.csv"
    print(f"  levels from {lvp.name}")
    lv = pd.read_csv(lvp, low_memory=False)
    lv["depth_m"] = pd.to_numeric(lv.value, errors="coerce") * 0.3048   # ft -> m
    lv = lv[lv.unit_of_measure.astype(str).str.strip().str.lower().eq("ft")]
    lv = lv[lv.depth_m.notna()]
    # Static readings only: a reading taken while pumping is not a water table.
    # Fan's four columns cannot support this filter; NWIS can.
    q = lv.qualifier.astype(str)
    lv["is_static"] = q.str.contains("Static", case=False, na=False)
    g = lv.groupby("monitoring_location_id").agg(
        wtd_m=("depth_m", "mean"), wtd_sd=("depth_m", "std"), n=("depth_m", "size"),
        frac_static=("is_static", "mean"),
        first=("time", "min"), last=("time", "max"))
    gs = lv[lv.is_static].groupby("monitoring_location_id")["depth_m"].mean().rename("wtd_m_static")
    g = g.join(gs).reset_index().rename(columns={"monitoring_location_id": "id"})
    df = g.merge(sites, on="id", how="inner")
    df["well_depth_m"] = pd.to_numeric(df.well_constructed_depth, errors="coerce") * 0.3048
    df["land_elev_m"] = pd.to_numeric(df.altitude, errors="coerce") * 0.3048
    # The same one-sided-filter trap build_earth_wtd_sites.py records: the test
    # is 0 < wtd <= well depth, not wtd <= well depth, or every bore shallower
    # than its own water level survives.
    df["depth_consistent"] = (df.wtd_m > 0) & (
        df.well_depth_m.isna() | (df.wtd_m <= df.well_depth_m))
    df["confinement"] = df.aquifer_type_code.map(AQUIFER_CONFINEMENT)
    df = df[df.lat.notna() & df.lon.notna() & df.wtd_m.notna()]
    p = OUT / "us_wtd_sites.csv"
    df.to_csv(p, index=False)
    print(f"\nSITES {len(df):,} from {int(df.n.sum()):,} readings -> {p}")
    print(f"  depth-consistent      {int(df.depth_consistent.sum()):,} "
          f"({df.depth_consistent.mean():.1%})")
    print(f"  with a well depth     {int(df.well_depth_m.notna().sum()):,}")
    print(f"  confinement           {df.confinement.value_counts(dropna=False).to_dict()}")
    unc = df[(df.confinement == "unconfined") & df.depth_consistent]
    print(f"  UNCONFINED and depth-consistent  {len(unc):,}  "
          f"median WTD {unc.wtd_m.median():.2f} m")
    d = df[df.depth_consistent]
    print(f"  WTD m: median {d.wtd_m.median():.2f}  5-95pct "
          f"{d.wtd_m.quantile(.05):.2f} to {d.wtd_m.quantile(.95):.2f}")
    prov = {"generated": datetime.datetime.now(datetime.timezone.utc)
            .isoformat(timespec="seconds"),
            "source": f"{BASE}/field-measurements + monitoring-locations",
            "parameter_code": PARAM_WTD, "bbox": BBOX,
            "sites": int(len(df)), "readings": int(df.n.sum())}
    (OUT / "us_wtd_sites.provenance.json").write_text(json.dumps(prov, indent=2) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stage", default="all",
                    choices=["sites", "levels", "assemble", "all"])
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--full-history", action="store_true",
                    help="every reading per site rather than the latest one; about "
                         "20x the paging, and what the Australian table's means use")
    a = ap.parse_args()
    print(f"  api key: {'yes' if api_key() else 'NO (anonymous, 4000 requests)'}", flush=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    for st in (["sites", "levels", "assemble"] if a.stage == "all" else [a.stage]):
        print(f"--- {st} ---", flush=True)
        if st == "levels":
            stage_levels(a.quiet, a.full_history)
        else:
            {"sites": stage_sites, "assemble": stage_assemble}[st](a.quiet)


if __name__ == "__main__":
    main()
