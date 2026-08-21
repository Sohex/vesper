#!/usr/bin/env python3
"""Fetch a window of the global recharge grid without downloading the globe.

    python hydrography/scripts/fetch_recharge.py --region aus
    python hydrography/scripts/fetch_recharge.py --region us-wide

Worldbuilding. Vesper is an invented planet; this fetches real recharge for the
Earth comparator that is the invented planet's only external check.

`docs/src/reference/external-data.md` records the source. What it cannot record
compactly, and what this script IS, are the exact windows and the trick that
makes them cheap: `RechargeTotal.nc` is a classic netCDF3 whose single variable
is a contiguous big-endian float64 array, so a byte RANGE over whole rows pulls
a region out of a multi-gigabyte file without fetching it. `BEGIN` is the offset
of the data section and the geometry below is the file's own grid.

Blocks are written as they arrive and a re-run resumes, because
`docs/src/reference/large-data.md` sets the test bluntly: kill it at 80% and run
it again, and if that costs more than the last chunk it is not checkpointed. The
first version of this failed that test at 85% of a 75 minute fetch, holding 1.6
GB in memory with nothing on disk.

Windows are DELIBERATELY WIDER than the region they serve where a bbox edge
would cut land. A cell with no coverage is not land, so it is ocean, so it is a
fixed head at sea level: `us-wide` reaches past Canada and Mexico for that
reason, where `aus` does not need to because Australia is an island.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import requests

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "hydrography" / "data" / "earth_validation_cache"
URL = "https://zenodo.org/records/7611675/files/RechargeTotal.nc"
BEGIN = 490224          # offset of the data section in the netCDF3 header
NX = 43200              # columns per row
ITEM = 8                # big-endian float64
LON0 = -179.99583333333334
LAT0 = 89.99583333333334
STEP = 1 / 120.0        # 30 arc-seconds

REGIONS = {
    "aus": {"rows": (12000, 16199), "cols": (34800, 40319), "out": "recharge_aus.npz"},
    "us": {"rows": (4800, 7920), "cols": (6600, 13680), "out": "recharge_us.npz"},
    # past Canada and Mexico, so neither border becomes a fake coast
    "us-wide": {"rows": (3600, 9000), "cols": (6000, 14400),
                "out": "recharge_us_wide.npz"},
}
BLOCK_ROWS = 210


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--region", required=True, choices=sorted(REGIONS))
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    R = REGIONS[a.region]
    r0, r1 = R["rows"]
    c0, c1 = R["cols"]
    out = a.out or (CACHE / R["out"])
    if out.exists():
        print(f"{out} exists; delete it to refetch")
        return
    nrow, ncol = r1 - r0 + 1, c1 - c0 + 1
    chunks = out.parent / f"{out.stem}_chunks"
    chunks.mkdir(parents=True, exist_ok=True)
    ident = {"url": URL, "rows": [r0, r1], "cols": [c0, c1], "block": BLOCK_ROWS}
    state_p = chunks / "state.json"
    if state_p.is_file() and json.loads(state_p.read_text()) != ident:
        raise SystemExit(f"{chunks} holds a different window. Delete it to refetch.")
    state_p.write_text(json.dumps(ident))
    starts = list(range(r0, r1 + 1, BLOCK_ROWS))
    print(f"{a.region}: {nrow} x {ncol}  "
          f"lat {LAT0 - r1 * STEP:.3f}..{LAT0 - r0 * STEP:.3f}  "
          f"lon {LON0 + c0 * STEP:.3f}..{LON0 + c1 * STEP:.3f}", flush=True)
    print(f"  {nrow * NX * ITEM / 1e9:.2f} GB of byte-range over {len(starts)} blocks",
          flush=True)
    s = requests.Session()
    t0 = time.time()
    for i, start in enumerate(starts):
        cp = chunks / f"block_{start:06d}.npy"
        if cp.is_file():
            print(f"    rows {start}: cached", flush=True)
            continue
        stop = min(start + BLOCK_ROWS - 1, r1)
        nb = stop - start + 1
        lo = BEGIN + start * NX * ITEM
        hi = BEGIN + (stop + 1) * NX * ITEM - 1
        for att in range(6):
            try:
                resp = s.get(URL, headers={"Range": f"bytes={lo}-{hi}"}, timeout=600)
                if resp.status_code in (206, 200) and len(resp.content) == hi - lo + 1:
                    blk = np.frombuffer(resp.content, dtype=">f8").reshape(nb, NX)
                    # written before the next request, so a kill costs one block
                    np.save(cp, blk[:, c0:c1 + 1].astype(np.float32))
                    break
                raise IOError(f"status {resp.status_code}, {len(resp.content)} bytes")
            except Exception as e:
                if att == 5:
                    raise SystemExit(f"block {start} failed: {e}")
                w = 8 * (att + 1)
                print(f"    retry {start} in {w}s ({e})", flush=True)
                time.sleep(w)
        print(f"    rows {start}-{stop}  {(i + 1) / len(starts):6.1%}  "
              f"{time.time() - t0:.0f}s", flush=True)
        time.sleep(3)
    grid = np.concatenate([np.load(chunks / f"block_{st:06d}.npy") for st in starts])
    assert grid.shape == (nrow, ncol), f"{grid.shape} != {(nrow, ncol)}"
    lat = LAT0 - np.arange(r0, r1 + 1) * STEP
    lon = LON0 + np.arange(c0, c1 + 1) * STEP
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, recharge=grid, lat=lat, lon=lon)
    f = np.isfinite(grid)
    v = grid[f]
    print(f"done in {time.time() - t0:.0f}s, finite {f.mean():.1%} -> {out}")
    print(f"recharge mm/yr: p50 {np.median(v):.2f}  p95 {np.percentile(v, 95):.1f}  "
          f"max {v.max():.1f}")
    print(f"  block cache kept at {chunks}; delete it once the npz is trusted")


if __name__ == "__main__":
    main()
