#!/usr/bin/env python3
"""GW-3: score the groundwater solver on Earth, where a water table is measured.

    python hydrography/scripts/earth_calibration.py --stage all --edge-km 15.19
    python hydrography/scripts/earth_calibration.py --stage all --edge-km 7.60

Worldbuilding. Vesper is an invented planet; this script is the one place its
groundwater solver is pointed at a REAL planet, because Vesper has no bores and
a model nothing can contradict is not a model. Everything else in `hydrography/`
runs on the invented world.

`hydrography/notes/earth-calibration-criterion.md` is the argument and declares
the thresholds BEFORE any scoring. This is that procedure made re-runnable, and
it exists because the first version of it was written inline and lost: the one
number judging this component was recorded but not reproducible (GW-20).

## What it does, and the one knob that matters

The same solver, operator and configuration as Vesper, on Earth's radius, with a
mesh built by a faithful port of Orogen's own Fibonacci generator so the
discretisation is the one GW-8 measured rather than a different one. Topography
is ETOPO, permeability is GLHYMPS, recharge is the fetched product, and the
score is against Australian bores assembled by `build_earth_wtd_sites.py`.

`--edge-km` is the knob. The original ran at 15.19 km, matched to Vesper's mesh
on purpose so the discretisation was shared. That choice is also the note's
declared confound: Fan et al.'s figures are near 1 km, and her own finding is
that terrain dominates water table depth at local scales, so a miss at 15.19 km
confounds FORMULATION with RESOLUTION. Running the same case at several edge
lengths is how that gets separated, and the fetched ETOPO is 1 arc-minute, about
1.855 km, so it supports the range without new data.

Stages cache into `hydrography/data/earth_validation_cache/<tag>/` and each is
skipped when its artifact is present. The cache is gitignored and regenerable;
the RESULT is written to `hydrography/analysis/earth_calibration.json`.
"""
from __future__ import annotations

import argparse
import datetime
import json
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(HERE))

import yaml  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402
import groundwater as gw  # noqa: E402
from orogen import LAND  # noqa: E402

CACHE = ROOT / "hydrography" / "data" / "earth_validation_cache"

# A cell with no DEM coverage has no elevation, so it is not land, so it is
# OCEAN and a fixed head at sea level. Australia is an island and its window's
# edges are real coast; the recovered harness still had to clip New Guinea for
# exactly this reason. A United States window cut at the Canadian and Mexican
# borders would put sea along both of them, so the modelled domain is WIDER than
# the scored one and `score_bbox` is what brings it back.
REGIONS = {
    "aus": {
        "dem": "etopo_aus.nc", "recharge": "recharge_aus.npz",
        "sites": "aus_wtd_sites.csv",
        "glhymps_bbox_deg": (108.0, -46.0, 158.0, -9.0),
        "dem_lat_max": -10.0,          # clip New Guinea
        "score_bbox": None,            # the window is already the region
    },
    # First pass on the narrow window, kept because the wide fetch is slow and
    # this needs no new data. Its edges cut Canada and Mexico, so those become
    # fake coast at sea level; `score_bbox` is pulled 3 degrees inside every
    # edge, about 330 km, so the scored bores sit outside that influence.
    "us-narrow": {
        "dem": "etopo_us.nc", "recharge": "recharge_us.npz",
        "sites": "us_wtd_sites.csv",
        "glhymps_bbox_deg": (-126.0, 23.0, -65.0, 51.0),
        "dem_lat_max": None,
        "score_bbox": (-122.0, 27.0, -69.0, 47.0),
    },
    "us": {
        "dem": "etopo_us_wide.nc", "recharge": "recharge_us_wide.npz",
        "sites": "us_wtd_sites.csv",
        "glhymps_bbox_deg": (-131.0, 14.0, -59.0, 61.0),
        "dem_lat_max": None,
        "score_bbox": (-125.0, 24.0, -66.0, 50.0),   # conterminous states only
    },
}
REFERENCE = ROOT / "hydrography" / "data" / "reference"
SITES = ROOT / "hydrography" / "data" / "earth_validation" / "aus_wtd_sites.csv"
EARTH_R_KM = 6371.0
VESPER_CELL_KM2 = 4 * np.pi * 7645.2 ** 2 / 2_500_001     # the 15.19 km baseline
OROGEN_SEED = 16236323
# The split seed for the benchmark, fixed so its number is a number and not a
# draw. Its predecessor was quoted in four documents and reproduced by nothing.
BENCH_SEED = 20260820
OROGEN_JITTER = 0.75
# Penman open-water potential ET, the arm groundwater.yaml names as et_max.
ET_MAX_MM_YR = 1500.0


def regions_for_edge(edge_km: float) -> int:
    """Orogen's own avgEdgeKm relation, pi*R/sqrt(N), inverted."""
    return int(round((np.pi * EARTH_R_KM / edge_km) ** 2))


def edge_for_regions(n: int) -> float:
    return float(np.pi * EARTH_R_KM / np.sqrt(n))


def make_rng(seed):
    s = (abs(int(seed * 9301 + 49297)) % 2147483646) + 1

    def r():
        nonlocal s
        s = (s * 16807) % 2147483647
        return (s - 1) / 2147483646
    return r


def orogen_fibonacci(n: int, jitter: float, seed: int):
    """Faithful port of vendor/orogen/js/sphere-mesh.js generateFibonacciSphere.

    Kept as a port rather than replaced by an equivalent-looking spiral: GW-8's
    finding is about THIS discretisation, and the operator's 12% noise floor was
    reproduced by an Earth mesh from this generator to within 2%.
    """
    rng = make_rng(seed)
    s = 3.6 / np.sqrt(n)
    dlong = np.pi * (3 - np.sqrt(5))
    dz = 2.0 / n
    lat = np.empty(n)
    lon = np.empty(n)
    z = 1 - dz / 2
    lng = 0.0
    for k in range(n):
        r = np.sqrt(max(0.0, 1 - z * z))
        lat_deg = np.degrees(np.arcsin(z))
        lon_deg = np.degrees(lng)
        if jitter > 0:
            j_lat = rng() - rng()
            j_lon = rng() - rng()
            next_z = max(-1.0, z - dz * 2 * np.pi * r / s)
            lat_deg += jitter * j_lat * (lat_deg - np.degrees(np.arcsin(next_z)))
            lon_deg += jitter * j_lon * (s / r * 180 / np.pi) if r > 0 else 0.0
        lat[k] = np.radians(lat_deg)
        lon[k] = np.radians(lon_deg)
        z -= dz
        lng += dlong
    return np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)


@dataclass
class EarthExport:
    """The seven attributes `groundwater.Geometry` and `solve` actually read."""
    radius_km: float
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    cell_area: np.ndarray
    surface_class: np.ndarray

    @property
    def n_regions(self) -> int:
        return self.x.size


def stage_mesh(tag: Path, n: int, quiet: bool) -> None:
    out = tag / "earth_xyz.npy"
    if out.exists():
        print(f"  mesh: cached ({np.load(out, mmap_mode='r').shape[1]:,} generators)")
        return
    t = time.time()
    x, y, z = orogen_fibonacci(n, OROGEN_JITTER, OROGEN_SEED)
    print(f"  mesh: {n:,} generators in {time.time()-t:.0f}s, "
          f"mean edge {edge_for_regions(n):.2f} km")
    np.save(out, np.stack([x, y, z]))


def solution_name(drain: str, et_lambda: float | None = None) -> str:
    """Every variant is its own artifact, never an overwrite of the baseline."""
    base = "earth_solution" if drain == "none" else f"earth_solution_drain-{drain}"
    return f"{base}.npz" if et_lambda is None else f"{base}_lam{et_lambda:g}.npz"


def region_dir(tag: Path, region: str) -> Path:
    """Mesh is shared per edge length; everything sampled onto it is per region."""
    d = tag / region
    d.mkdir(parents=True, exist_ok=True)
    return d


def stage_fields(tag: Path, quiet: bool, region: str) -> None:
    """ETOPO onto the mesh. Clips New Guinea: a bbox edge through a real
    landmass is a fake coast, and a fake coast is a fixed head at sea level."""
    R = REGIONS[region]
    out = region_dir(tag, region) / "earth_fields.npz"
    if out.exists():
        print("  fields: cached")
        return
    import netCDF4 as nc
    xyz = np.load(tag / "earth_xyz.npy")
    n = xyz.shape[1]
    d = nc.Dataset(str(CACHE / R["dem"]))
    z = np.asarray(d["z"][:], dtype=np.float32)
    la = np.asarray(d["latitude"][:])
    lo = np.asarray(d["longitude"][:])
    if R["dem_lat_max"] is not None:
        keep = la <= R["dem_lat_max"]
        z, la = z[keep], la[keep]
    LO, LA = np.meshgrid(np.radians(lo), np.radians(la))
    px = (np.cos(LA) * np.cos(LO)).ravel()
    py = (np.cos(LA) * np.sin(LO)).ravel()
    pz = np.sin(LA).ravel()
    _, idx = cKDTree(xyz.T).query(np.stack([px, py, pz], axis=1), workers=-1)
    zf = z.ravel().astype(np.float64)
    cnt = np.bincount(idx, minlength=n)
    s = np.bincount(idx, weights=zf, minlength=n)
    s2 = np.bincount(idx, weights=zf * zf, minlength=n)
    elev = np.full(n, np.nan)
    has = cnt > 0
    elev[has] = s[has] / cnt[has]
    # The WITHIN-CELL relief, which the mean throws away and GW-23 needs: a cell
    # drains to the gully inside it, not across itself to a mapped river.
    elev_sd = np.full(n, np.nan)
    elev_sd[has] = np.sqrt(np.maximum(s2[has] / cnt[has] - elev[has] ** 2, 0.0))
    elev_min = np.full(n, np.inf)
    np.minimum.at(elev_min, idx, zf)
    elev_min[~np.isfinite(elev_min)] = np.nan
    land = np.zeros(n, bool)
    land[has] = elev[has] > 0.0
    cell_km2 = 4 * np.pi * EARTH_R_KM ** 2 / n
    print(f"  fields: {has.sum():,} cells with DEM, {land.sum():,} land, "
          f"{land.sum()*cell_km2:,.0f} km2 (Australia is 7,688,000)")
    np.savez_compressed(out, elev=elev, land=land, cnt=cnt,
                        elev_sd=elev_sd, elev_min=elev_min)


def stage_perm(tag: Path, quiet: bool, region: str) -> None:
    out = region_dir(tag, region) / "earth_perm.npz"
    if out.exists():
        print("  permeability: cached")
        return
    import geopandas as gpd
    xyz = np.load(tag / "earth_xyz.npy")
    f = np.load(region_dir(tag, region) / "earth_fields.npz")
    land, elev = f["land"], f["elev"]
    li = np.flatnonzero(land)
    lat = np.degrees(np.arcsin(np.clip(xyz[2, li], -1, 1)))
    lon = np.degrees(np.arctan2(xyz[1, li], xyz[0, li]))
    pts = gpd.GeoDataFrame({"cell": li, "elev": elev[li]},
                           geometry=gpd.points_from_xy(lon, lat), crs="EPSG:4326")
    gdb = REFERENCE / "GLHYMPS" / "GLHYMPS.gdb"
    r = 6378137.0
    lo0, la0, lo1, la1 = REGIONS[region]["glhymps_bbox_deg"]
    bbox = (r * np.radians(lo0), r * np.sin(np.radians(la0)),
            r * np.radians(lo1), r * np.sin(np.radians(la1)))
    gdf = gpd.read_file(str(gdb), layer="Final_GLHYMPS_Polygon", bbox=bbox,
                        engine="pyogrio",
                        columns=["Permeability_no_permafrost",
                                 "Permeability_standard_deviation", "Porosity"])
    j = gpd.sjoin(pts.to_crs(gdf.crs), gdf, how="left", predicate="within")
    j = j[~j.index.duplicated(keep="first")]
    logk = j.Permeability_no_permafrost.to_numpy()
    sd = j.Permeability_standard_deviation.to_numpy()
    print(f"  permeability: {len(gdf):,} polygons, matched "
          f"{np.isfinite(logk).mean():.1%} of land cells, "
          f"median log10 k {np.nanmedian(logk):.2f}")
    np.savez_compressed(out, cell=j.cell.to_numpy(), logk=logk, logk_sd=sd)


def _csr(src, dst, n):
    """Symmetric CSR adjacency. The recovered harness used a dict of lists,
    which is fine at 1.7M cells and neither fast nor small at 7M."""
    a = np.concatenate([src, dst])
    b = np.concatenate([dst, src])
    o = np.argsort(a, kind="stable")
    a, b = a[o], b[o]
    off = np.zeros(n + 1, np.int64)
    np.add.at(off, a + 1, 1)
    np.cumsum(off, out=off)
    return off, b


def stage_drainage(tag: Path, quiet: bool, region: str) -> None:
    """Priority flood, steepest descent, accumulated area. The river network is
    DERIVED rather than assumed, the same shape hydrography does for Vesper,
    because the fixed heads GW-17 added need to know where the rivers are."""
    out = region_dir(tag, region) / "earth_drainage.npz"
    if out.exists():
        print("  drainage: cached")
        return
    import heapq
    xyz = np.load(tag / "earth_xyz.npy")
    n = xyz.shape[1]
    f = np.load(region_dir(tag, region) / "earth_fields.npz")
    elev = np.nan_to_num(f["elev"], nan=0.0)
    land = f["land"]
    ex = EarthExport(EARTH_R_KM, xyz[0], xyz[1], xyz[2],
                     np.full(n, 4 * np.pi * EARTH_R_KM ** 2 / n),
                     np.where(land, LAND, 0).astype(np.int16))
    geom = gw.Geometry(ex)
    off, nb = _csr(geom.src, geom.dst, n)
    A = geom.voronoi_area_m2

    t = time.time()
    filled = np.where(land, elev, -1e30).astype(np.float64)
    done = ~land
    h = []
    for i in np.flatnonzero(~land):
        for v in nb[off[i]:off[i + 1]]:
            if land[v] and not done[v]:
                done[v] = True
                heapq.heappush(h, (float(elev[v]), int(v)))
    while h:
        e, u = heapq.heappop(h)
        filled[u] = e
        for v in nb[off[u]:off[u + 1]]:
            if not done[v]:
                done[v] = True
                heapq.heappush(h, (max(float(elev[v]), e), int(v)))
    print(f"  drainage: priority flood {time.time()-t:.0f}s")

    # steepest descent on the filled surface, as a segment-min over CSR
    deg = np.diff(off)
    fv = filled[nb]
    fv[~land[nb]] = -1e30                       # the sea is always downhill
    best = np.full(n, -1, np.int64)
    ok = deg > 0
    idx = np.arange(n)[ok]
    starts = off[:-1][ok]
    pos = np.zeros(idx.size, np.int64)
    for j in range(idx.size):
        a, b = off[idx[j]], off[idx[j] + 1]
        pos[j] = a + int(np.argmin(fv[a:b]))
    cand = nb[pos]
    take = filled[cand] < filled[idx]
    best[idx[take]] = cand[take]

    acc = A.copy()
    for u in np.argsort(-filled):
        r = best[u]
        if land[u] and r >= 0:
            acc[r] += acc[u]
    km2 = acc / 1e6
    for thr in (1e3, 5e3, 1e4):
        m = land & (km2 >= thr)
        print(f"    accumulation >= {thr:>7.0f} km2: {m.sum():7,} cells "
              f"({m.sum()/max(land.sum(),1):5.1%} of land)")
    np.savez_compressed(out, filled=filled, recv=best, acc_km2=km2)


def stage_solve(tag: Path, quiet: bool, river_km2: float | None, region: str,
                drain: str = "none", drain_min_relief_m: float = 10.0,
                et_lambda: float | None = None) -> None:
    rd = region_dir(tag, region)
    out = rd / solution_name(drain, et_lambda)
    if out.exists():
        print("  solve: cached")
        return
    cfg = yaml.safe_load((ROOT / "hydrography" / "config" / "groundwater.yaml")
                         .read_text(encoding="utf-8"))
    D = cfg["aquifer"]["thickness_m"]
    rho = cfg["fluid"]["density_kg_m3"]
    mu = cfg["fluid"]["dynamic_viscosity_pa_s"]
    xyz = np.load(tag / "earth_xyz.npy")
    n = xyz.shape[1]
    f = np.load(rd / "earth_fields.npz")
    elev, land = f["elev"], f["land"]
    perm = np.load(rd / "earth_perm.npz")
    logk = np.full(n, np.nan)
    logk[perm["cell"]] = perm["logk"]

    rz = np.load(CACHE / REGIONS[region]["recharge"])
    rr = rz["recharge"]
    LO, LA = np.meshgrid(np.radians(rz["lon"]), np.radians(rz["lat"]))
    ok = np.isfinite(rr)
    _, idx = cKDTree(xyz.T).query(
        np.stack([(np.cos(LA) * np.cos(LO))[ok], (np.cos(LA) * np.sin(LO))[ok],
                  np.sin(LA)[ok]], axis=1), workers=-1)
    cnt = np.bincount(idx, minlength=n)
    ssum = np.bincount(idx, weights=rr[ok].astype(np.float64), minlength=n)
    rech = np.zeros(n)
    h = cnt > 0
    rech[h] = ssum[h] / cnt[h]

    cell_km2 = 4 * np.pi * EARTH_R_KM ** 2 / n
    ex = EarthExport(EARTH_R_KM, xyz[0], xyz[1], xyz[2],
                     np.full(n, cell_km2),
                     np.where(land, LAND, 0).astype(np.int16))
    t = time.time()
    geom = gw.Geometry(ex)
    ex.cell_area = geom.voronoi_area_m2 / 1e6
    geom.volume_area_m2 = geom.voronoi_area_m2
    print(f"  geometry: {time.time()-t:.0f}s, faces {geom.src.size:,}, "
          f"area closure {geom.voronoi_area_m2.sum()/(4*np.pi*(EARTH_R_KM*1e3)**2):.10f}")
    K = np.nan_to_num(gw.conductivity(10.0 ** np.nan_to_num(logk, nan=-30.0),
                                      rho, 9.80665, mu), nan=0.0)
    cond = land & np.isfinite(logk) & (rech > 0)

    # GW-15's sink and GW-17's baselevels are BOTH ON, as config/pipeline.yaml's
    # gate says: neither is optional physics. Without the sink the table pins at
    # the surface over most of the land -- the recovered pre-GW-15 harness does
    # exactly that, at 93.8% pinned -- and without local baselevels the sea is
    # the only boundary and the depth range collapses.
    et = cfg.get("evapotranspiration", {})
    YR = 365.25 * 86400.0
    et_max = ET_MAX_MM_YR / 1000.0 / YR
    et_lambda = float(et_lambda if et_lambda is not None else et.get("lambda_m", 1.0))
    # GW-23's sub-grid drain at its C -> infinity limit. A true Robin condition
    # needs a solver term; pinning the head at the cell's own valley floor is the
    # infinitely-conductive END of that family and BOUNDS what it can do. Only
    # cells with relief to vent into are drained: a flat cell has no gully.
    fixed = None
    if drain != "none":
        z_valley = (f["elev_min"] if drain == "min"
                    else np.nan_to_num(elev, nan=0.0) - 2.0 * f["elev_sd"])
        relief = np.nan_to_num(elev, nan=0.0) - z_valley
        vent = land & np.isfinite(z_valley) & (relief >= drain_min_relief_m)
        fixed = np.full(n, np.nan)
        fixed[vent] = z_valley[vent]
        print(f"  sub-grid drain ({drain}): {vent.sum():,} of {land.sum():,} land cells "
              f"vent, median relief {np.median(relief[vent]):.1f} m")
    if river_km2 is not None:
        dr = np.load(rd / "earth_drainage.npz")
        if fixed is None:
            fixed = np.full(n, np.nan)
        wet = land & (dr["acc_km2"] >= river_km2)
        fixed[wet] = np.nan_to_num(elev, nan=0.0)[wet]
        print(f"  baselevels: {wet.sum():,} river cells at >= {river_km2:.0e} km2")

    t = time.time()
    res = gw.solve(ex, geom, k0_m_s=K, thickness_m=D,
                   recharge_m_s=rech / 1000.0 / (365.25 * 86400.0),
                   surface_m=np.nan_to_num(elev, nan=0.0), conductive=cond,
                   sea_level_m=0.0, et_max_m_s=et_max, et_lambda_m=et_lambda,
                   fixed_head_m=fixed, verbose=not quiet)
    dep = res["depth_m"]
    print(f"  solve: {time.time()-t:.0f}s, conductive land {cond.sum():,}, "
          f"pinned {res['pinned'][cond].mean():.1%}, "
          f"median depth {np.median(dep[cond]):.2f} m")
    np.savez_compressed(out, depth=dep, head=res["head_m"], elev=elev, land=land,
                        cond=cond, rech=rech, K=K, pinned=res["pinned"])


def stage_score(tag: Path, edge_km: float, quiet: bool, region: str,
                confinement: str | None, drain: str = "none",
                surface: str = "cell-mean", et_lambda: float | None = None) -> dict:
    """Skill, not just residual moments. The bar and the ceiling are the note's."""
    import pandas as pd
    R = REGIONS[region]
    sol = np.load(region_dir(tag, region) / solution_name(drain, et_lambda))
    xyz = np.load(tag / "earth_xyz.npy")
    obs = pd.read_csv(ROOT / "hydrography" / "data" / "earth_validation" / R["sites"],
                      low_memory=False)
    if "depth_consistent" in obs.columns:
        obs = obs[obs.depth_consistent.astype(bool)]
    # Score only inside the region, never the margin the domain was widened by.
    if R["score_bbox"]:
        x0, y0, x1, y1 = R["score_bbox"]
        obs = obs[obs.lon.between(x0, x1) & obs.lat.between(y0, y1)]
    # A confined aquifer's bore reads a potentiometric head, not a water table.
    # Australia can only guess at this; the US labels it.
    if confinement and "confinement" in obs.columns:
        obs = obs[obs.confinement == confinement]
    if obs.empty:
        raise SystemExit(f"no sites left for region={region} confinement={confinement}")
    la = np.radians(obs.lat.to_numpy())
    lo = np.radians(obs.lon.to_numpy())
    p = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=1)
    _, cell = cKDTree(xyz.T).query(p, workers=-1)
    use = sol["cond"][cell]
    # WHICH SURFACE the depth hangs from. The model solves a head; a depth needs
    # a ground level, and the cell mean is not where the bore is. Both options
    # take elevation from the DEM, never from the observation's own survey: the
    # two disagree at sd 51.97 m and mixing them puts all of that on the residual.
    if surface == "dem-at-bore":
        fld = np.load(region_dir(tag, region) / "earth_fields.npz")
        import netCDF4 as nc
        dd = nc.Dataset(str(CACHE / R["dem"]))
        zz = np.asarray(dd["z"][:], float)
        dla = np.asarray(dd["latitude"][:]); dlo = np.asarray(dd["longitude"][:])
        blon = obs.lon.to_numpy() % 360.0 if dlo.max() > 180 else obs.lon.to_numpy()
        ri = np.clip(np.searchsorted(dla, obs.lat.to_numpy()), 0, dla.size - 1)
        if dla[0] > dla[-1]:
            ri = np.clip(dla.size - 1 - np.searchsorted(dla[::-1], obs.lat.to_numpy()),
                         0, dla.size - 1)
        ci = np.clip(np.searchsorted(dlo, blon), 0, dlo.size - 1)
        ground = zz[ri, ci].astype(float)
        # ground is per BORE, head is per CELL: index each on its own axis.
        model = ground[use] - sol["head"][cell][use]
    else:
        model = sol["depth"][cell][use]
    observed = obs.wtd_m.to_numpy()[use]
    r = model - observed

    # cell-mean skill: the model has ONE value per cell, so scoring per bore
    # charges it for within-cell variance it cannot represent. Both are reported.
    import collections
    agg = collections.defaultdict(list)
    magg = collections.defaultdict(list)
    for c, o, mv in zip(cell[use], observed, model):
        agg[c].append(o)
        magg[c].append(mv)
    cells = np.array(sorted(agg))
    obs_cell = np.array([np.mean(agg[c]) for c in cells])
    # Aggregate the MODEL the same way as the observations, over the same bores.
    # Reading `sol["depth"]` here instead silently ignored `--surface`, which is
    # the class of bug where an option appears to run and changes nothing.
    mod_cell = np.array([np.mean(magg[c]) for c in cells])
    within = sum(((np.array(agg[c]) - np.mean(agg[c])) ** 2).sum() for c in cells)
    total = ((observed - observed.mean()) ** 2).sum()
    ceiling = 1.0 - within / total

    def r2(m, o):
        return float(1.0 - ((m - o) ** 2).sum() / ((o - o.mean()) ** 2).sum())

    out = {
        "region": region,
        "drain": drain,
        "surface": surface,
        "confinement": confinement or "all",
        "edge_km": edge_km,
        "n_regions": int(xyz.shape[1]),
        "sites_scored": int(use.sum()),
        "cells_with_bores": int(cells.size),
        "observed_sd_m": float(observed.std()),
        "model_sd_m": float(model.std()),
        "residual_mean_m": float(r.mean()),
        "residual_sd_m": float(r.std()),
        "r2_per_bore": r2(model, observed),
        "r2_cell_mean": r2(mod_cell, obs_cell),
        "pearson_cell_mean": float(np.corrcoef(mod_cell, obs_cell)[0, 1]),
        "ceiling_r2_at_this_cell_size": float(ceiling),
        "within_cell_variance_fraction": float(within / total),
    }
    print(f"  score: {out['sites_scored']:,} bores in {out['cells_with_bores']:,} cells")
    print(f"    ceiling at this cell size   {out['ceiling_r2_at_this_cell_size']:.4f}")
    print(f"    R2 on cell means            {out['r2_cell_mean']:+.4f}   "
          f"(bar 0.07, statistical fit 0.2098)")
    print(f"    Pearson on cell means       {out['pearson_cell_mean']:+.4f}")
    print(f"    observed sd {out['observed_sd_m']:.2f} m   model sd {out['model_sd_m']:.2f} m")
    return out


def stage_benchmark(tag: Path, region: str, confinement: str | None) -> dict:
    """What the RESOLVABLE FIELDS support, which is the model's real target.

    This exists because the number it replaces did not. Four documents and this
    script's own score line quoted 0.28 as what a flexible fit reaches on
    held-out cells; nothing ever computed it. Reproducing the original gives
    +0.0703 linear and +0.1215 binned, and a gradient boosting fit reaches
    +0.1435 with the physical model EXCLUDED. A criterion nothing can reproduce
    is worse than none, so the criterion now ships with its computation.

    Target is `log(1 + cell-mean depth)` on cells with at least two bores, which
    is the original's choice and is kept so the numbers stay comparable: depth is
    heavily right-skewed and a raw-depth R2 is dominated by a few deep cells.
    """
    import collections
    import pandas as pd
    from sklearn.ensemble import HistGradientBoostingRegressor
    R = REGIONS[region]
    rd = region_dir(tag, region)
    # The BASELINE solution by name, not `solution_name(drain, et_lambda)`, and
    # deliberately: nothing here reads the model's head or depth. The physical
    # model is EXCLUDED from this fit -- that is what makes the number a ceiling
    # on what the resolvable fields support -- and the file is opened only for
    # the land and conductive masks, which no drain or sink moves. A variant run
    # therefore has nothing to forward here; it would fail loudly on a missing
    # file rather than score the wrong solution.
    sol = np.load(rd / "earth_solution.npz")
    fld = np.load(rd / "earth_fields.npz")
    perm = np.load(rd / "earth_perm.npz")
    dr = np.load(rd / "earth_drainage.npz")
    xyz = np.load(tag / "earth_xyz.npy")
    n = xyz.shape[1]
    elev = np.nan_to_num(sol["elev"], nan=0.0)
    logk = np.full(n, np.nan)
    logk[perm["cell"]] = perm["logk"]

    obs = pd.read_csv(ROOT / "hydrography" / "data" / "earth_validation" / R["sites"],
                      low_memory=False)
    if "depth_consistent" in obs.columns:
        obs = obs[obs.depth_consistent.astype(bool)]
    if confinement and "confinement" in obs.columns:
        obs = obs[obs.confinement == confinement]
    if R["score_bbox"]:
        x0, y0, x1, y1 = R["score_bbox"]
        obs = obs[obs.lon.between(x0, x1) & obs.lat.between(y0, y1)]
    la = np.radians(obs.lat.to_numpy())
    lo = np.radians(obs.lon.to_numpy())
    _, cid = cKDTree(xyz.T).query(
        np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=1),
        workers=-1)
    df = pd.DataFrame({"cell": cid, "wtd": obs.wtd_m.to_numpy()})
    nn = df.groupby("cell").size()
    cm = df[df.cell.isin(nn[nn >= 2].index)].groupby("cell")["wtd"].mean()
    c = cm.index.to_numpy()
    y = np.log1p(np.clip(cm.to_numpy(), 0, None))
    X = np.column_stack([
        elev[c], np.nan_to_num(fld["elev_sd"][c], nan=0.0),
        np.log10(np.maximum(dr["acc_km2"][c], 1e-3)),
        np.log1p(np.clip(sol["rech"][c], 0, None)),
        np.nan_to_num(logk[c], nan=-13.0),
        (elev - np.nan_to_num(fld["elev_min"], nan=0.0))[c],
        sol["depth"][c]])
    ok = np.all(np.isfinite(X), axis=1)
    X, y = X[ok], y[ok]
    m = len(y)
    rng = np.random.default_rng(BENCH_SEED)
    o = rng.permutation(m)
    tr, te = o[:m // 2], o[m // 2:]

    def r2(pred, truth):
        return float(1 - ((truth - pred) ** 2).sum() / ((truth - truth.mean()) ** 2).sum())

    Xs = (X - X[tr].mean(0)) / np.maximum(X[tr].std(0), 1e-12)
    A = np.column_stack([np.ones(m), Xs])
    beta, *_ = np.linalg.lstsq(A[tr], y[tr], rcond=None)
    lin = r2(A[te] @ beta, y[te])
    g = HistGradientBoostingRegressor(max_iter=400, random_state=0).fit(X[tr], y[tr])
    boost_all = r2(g.predict(X[te]), y[te])
    g5 = HistGradientBoostingRegressor(max_iter=400, random_state=0).fit(
        np.delete(X[tr], 6, 1), y[tr])
    boost_fields = r2(g5.predict(np.delete(X[te], 6, 1)), y[te])
    gm = HistGradientBoostingRegressor(max_iter=400, random_state=0).fit(
        X[tr][:, [6]], y[tr])
    model_only = r2(gm.predict(X[te][:, [6]]), y[te])
    out = {"cells": int(m), "linear_six": lin, "boosting_all": boost_all,
           "boosting_fields_only": boost_fields, "boosting_model_only": model_only}
    print(f"  benchmark on {m:,} multi-bore cells, target log(1+depth), 50/50 held out")
    print(f"    linear, all predictors            R2 = {lin:+.4f}")
    print(f"    boosting, all predictors          R2 = {boost_all:+.4f}")
    print(f"    boosting, MODEL EXCLUDED          R2 = {boost_fields:+.4f}  <- the target")
    print(f"    boosting, MODEL ALONE             R2 = {model_only:+.4f}")
    if boost_all < boost_fields:
        print("    the model is INFORMATION-NEGATIVE: the fit does better without it")
    return out


def stage_calibrate(tag: Path, region: str, confinement: str | None,
                    drain: str, et_lambda: float | None, splits: int = 20) -> dict:
    """Held-out skill after a two-parameter bias correction, over many splits.

    Direct R2 is negative because the model runs shallow, so the question worth
    asking is whether the PATTERN survives once bias and scale are removed. One
    split is not an answer: the first one tried here read +0.0748 and the mean
    over twenty is +0.068, so a single favourable draw would have reported a
    pass on the declared 0.07 bar that twenty splits do not support.
    """
    import collections
    import pandas as pd
    from scipy.stats import pearsonr
    R = REGIONS[region]
    sol = np.load(region_dir(tag, region) / solution_name(drain, et_lambda))
    xyz = np.load(tag / "earth_xyz.npy")
    obs = pd.read_csv(ROOT / "hydrography" / "data" / "earth_validation" / R["sites"],
                      low_memory=False)
    if "depth_consistent" in obs.columns:
        obs = obs[obs.depth_consistent.astype(bool)]
    if confinement and "confinement" in obs.columns:
        obs = obs[obs.confinement == confinement]
    if R["score_bbox"]:
        x0, y0, x1, y1 = R["score_bbox"]
        obs = obs[obs.lon.between(x0, x1) & obs.lat.between(y0, y1)]
    la = np.radians(obs.lat.to_numpy())
    lo = np.radians(obs.lon.to_numpy())
    _, cid = cKDTree(xyz.T).query(
        np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=1),
        workers=-1)
    keep = sol["cond"][cid]
    agg = collections.defaultdict(list)
    for cc, oo in zip(cid[keep], obs.wtd_m.to_numpy()[keep]):
        agg[cc].append(oo)
    cs = np.array(sorted(agg))
    y = np.array([np.mean(agg[c]) for c in cs])
    x = sol["depth"][cs]

    def r2(p, t):
        return float(1 - ((p - t) ** 2).sum() / ((t - t.mean()) ** 2).sum())

    tes = []
    for seed in range(splits):
        rng = np.random.default_rng(seed)
        o = rng.permutation(len(y))
        tr, te = o[:len(y) // 2], o[len(y) // 2:]
        A = np.column_stack([x[tr], np.ones(tr.size)])
        b, *_ = np.linalg.lstsq(A, y[tr], rcond=None)
        tes.append(r2(np.column_stack([x[te], np.ones(te.size)]) @ b, y[te]))
    tes = np.array(tes)
    rho = float(pearsonr(x, y).statistic)
    out = {"cells": int(len(y)), "splits": splits, "pearson": rho,
           "held_out_r2_mean": float(tes.mean()), "held_out_r2_sd": float(tes.std()),
           "held_out_r2_min": float(tes.min()), "held_out_r2_max": float(tes.max()),
           "splits_clearing_bar": int((tes >= 0.07).sum()), "bar": 0.07}
    print(f"  calibrate on {len(y):,} cells, {splits} random 50/50 splits")
    print(f"    Pearson over all cells   {rho:+.4f}  (rho^2 = {rho**2:.4f})")
    print(f"    held-out R2  mean {tes.mean():+.4f}  sd {tes.std():.4f}  "
          f"range {tes.min():+.4f} to {tes.max():+.4f}")
    print(f"    splits clearing the {0.07} bar: {(tes >= 0.07).sum()} of {splits}")
    return out


def stage_diagnostics(tag: Path, region: str, confinement: str | None,
                      river_km2: float | None = 1e4) -> dict:
    """Which term set the depth, and what out-predicts the model.

    Four tables the notes cite and nothing could reproduce. Each is cheap, and
    each answered a question that looked like it needed a solver change.

    The SINK sweep shows the steady state is the local balance
    `d = lambda ln(et_max A / supply)`: median depth is linear in lambda while
    the correlation barely moves. The TRANSMISSIVITY sweep shows the dependence
    structure IS reachable by raising `D` -- terrain correlation goes from about
    zero to +0.47 -- and that agreement does not follow. The ATTRIBUTION shows
    `delta`, where a bore sits inside its own cell, out-predicting the model on
    Pearson, which is why GW-25 did not make the bore surface the default. And
    the DRAIN scan shows every terrain-derived drain elevation anti-correlating
    with observed depth, which is why GW-23 closed.
    """
    import collections
    import pandas as pd
    import netCDF4 as nc
    from scipy.stats import pearsonr, spearmanr
    import groundwater as gw
    from orogen import LAND
    R = REGIONS[region]
    rd = region_dir(tag, region)
    xyz = np.load(tag / "earth_xyz.npy")
    n = xyz.shape[1]
    sol = np.load(rd / "earth_solution.npz")
    dr = np.load(rd / "earth_drainage.npz")
    elev = np.nan_to_num(sol["elev"], nan=0.0)
    land, cond, K, rech = sol["land"], sol["cond"], sol["K"], sol["rech"]

    obs = pd.read_csv(ROOT / "hydrography" / "data" / "earth_validation" / R["sites"],
                      low_memory=False)
    if "depth_consistent" in obs.columns:
        obs = obs[obs.depth_consistent.astype(bool)]
    if confinement and "confinement" in obs.columns:
        obs = obs[obs.confinement == confinement]
    if R["score_bbox"]:
        x0, y0, x1, y1 = R["score_bbox"]
        obs = obs[obs.lon.between(x0, x1) & obs.lat.between(y0, y1)]
    la = np.radians(obs.lat.to_numpy())
    lo = np.radians(obs.lon.to_numpy())
    _, cid = cKDTree(xyz.T).query(
        np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=1),
        workers=-1)
    keep = cond[cid]
    agg = collections.defaultdict(list)
    for c, o in zip(cid[keep], obs.wtd_m.to_numpy()[keep]):
        agg[c].append(o)
    cs = np.array(sorted(agg))
    oc = np.array([np.mean(agg[c]) for c in cs])

    ex = EarthExport(EARTH_R_KM, xyz[0], xyz[1], xyz[2],
                     np.full(n, 4 * np.pi * EARTH_R_KM ** 2 / n),
                     np.where(land, LAND, 0).astype(np.int16))
    geom = gw.Geometry(ex)
    ex.cell_area = geom.voronoi_area_m2 / 1e6
    geom.volume_area_m2 = geom.voronoi_area_m2
    YR = 365.25 * 86400.0
    fh = np.full(n, np.nan)
    # The river threshold comes from `--river-km2`, the same option `stage_solve`
    # takes, rather than from a literal here. It stood at a literal 1e4 and the
    # option therefore ran and changed nothing for this stage; the default is
    # unchanged, so every table below reads as it did. `--et-lambda`, `--drain`
    # and `--surface` are the flags this stage genuinely IGNORES, and it must:
    # the sink sweep, the drain scan and the attribution table each set those
    # for themselves across their own arms, and taking a command-line value
    # would fix one arm of a sweep whose point is to vary it.
    wet = land & (dr["acc_km2"] >= (river_km2 if river_km2 is not None else np.inf))
    fh[wet] = elev[wet]
    base = dict(k0_m_s=K, recharge_m_s=rech / 1000.0 / YR, surface_m=elev,
                conductive=cond, sea_level_m=0.0, fixed_head_m=fh, verbose=False)
    out = {"sink_sweep": [], "transmissivity_sweep": [], "attribution": {},
           "drain_scan": []}

    def run(**kw):
        r = gw.solve(ex, geom, **base, **kw)
        d = r["depth_m"]
        md = d[cs]
        return (float(np.median(d[cond])), float(d[cond].std()),
                float(pearsonr(md, oc).statistic) if md.std() > 1e-9 else float("nan"),
                float(spearmanr(md, rech[cs]).statistic),
                float(spearmanr(md, elev[cs]).statistic))

    print(f"  {len(cs):,} cells with bores\n")
    print(f"  SINK: does lambda set the depth locally?")
    print(f"    {'lambda':>8s}{'median d':>10s}{'sd':>8s}{'r(mod,obs)':>12s}{'rho rech':>10s}{'rho elev':>10s}")
    for lam in (1.0, 2.0, 5.0, 10.0):
        v = run(thickness_m=100.0, et_max_m_s=1500.0 / 1000.0 / YR, et_lambda_m=lam)
        out["sink_sweep"].append({"lambda_m": lam, "median_depth_m": v[0],
                                  "sd_m": v[1], "pearson": v[2],
                                  "rho_recharge": v[3], "rho_elevation": v[4]})
        print(f"    {lam:8.1f}{v[0]:10.2f}{v[1]:8.2f}{v[2]:+12.4f}{v[3]:+10.3f}{v[4]:+10.3f}")

    print(f"\n  TRANSMISSIVITY: does raising D restore the dependence structure?")
    print(f"    {'D (m)':>8s}{'median d':>10s}{'sd':>8s}{'r(mod,obs)':>12s}{'rho rech':>10s}{'rho elev':>10s}")
    for D in (100.0, 1000.0, 5000.0, 20000.0):
        v = run(thickness_m=D, et_max_m_s=1500.0 / 1000.0 / YR, et_lambda_m=1.0)
        out["transmissivity_sweep"].append({"thickness_m": D, "median_depth_m": v[0],
                                            "sd_m": v[1], "pearson": v[2],
                                            "rho_recharge": v[3], "rho_elevation": v[4]})
        print(f"    {D:8.0f}{v[0]:10.2f}{v[1]:8.2f}{v[2]:+12.4f}{v[3]:+10.3f}{v[4]:+10.3f}")
    print(f"    observed: rho vs recharge {spearmanr(oc, rech[cs]).statistic:+.3f}, "
          f"vs elevation {spearmanr(oc, elev[cs]).statistic:+.3f}")

    # attribution and the drain scan both need the DEM inside each cell
    d_ = nc.Dataset(str(CACHE / R["dem"]))
    z = np.asarray(d_["z"][:], np.float32)
    dla = np.asarray(d_["latitude"][:])
    dlo = np.asarray(d_["longitude"][:])
    if R["dem_lat_max"] is not None:
        k = dla <= R["dem_lat_max"]
        z, dla = z[k], dla[k]
    blon = obs.lon.to_numpy() % 360.0 if dlo.max() > 180 else obs.lon.to_numpy()
    ri = (np.clip(dla.size - 1 - np.searchsorted(dla[::-1], obs.lat.to_numpy()), 0, dla.size - 1)
          if dla[0] > dla[-1] else
          np.clip(np.searchsorted(dla, obs.lat.to_numpy()), 0, dla.size - 1))
    ci = np.clip(np.searchsorted(dlo, blon), 0, dlo.size - 1)
    zb = z[ri, ci].astype(float)[keep]
    dg = collections.defaultdict(list)
    for c, v in zip(cid[keep], zb - elev[cid[keep]]):
        dg[c].append(v)
    dc = np.array([np.mean(dg[c]) for c in cs])
    md = (elev - sol["head"])[cs]
    print(f"\n  ATTRIBUTION: what predicts observed cell-mean depth?")
    for nm, v in (("the model alone", md), ("delta alone (bore within cell)", dc),
                  ("model + delta", md + dc), ("model + 0.2 delta", md + 0.2 * dc)):
        pr, sp = pearsonr(v, oc).statistic, spearmanr(v, oc).statistic
        out["attribution"][nm] = {"pearson": float(pr), "spearman": float(sp)}
        print(f"    {nm:34s} pearson {pr:+.4f}   spearman {sp:+.4f}")

    # GW-25: the two scoring surfaces are the ends of one family, a water table
    # as a subdued replica of terrain, h(x) = h_cell + alpha (z(x) - z_cell).
    # alpha = 1 is the cell mean, depth constant; alpha = 0 is the bore surface,
    # head flat. The three measures disagree about where to stand, which is the
    # sign that alpha is not being chosen by the physics.
    print(f"\n  SUBDUED REPLICA: depth = cell_depth + (1 - alpha) delta")
    print(f"    {'alpha':>7s}{'model sd':>10s}{'pearson':>10s}{'spearman':>10s}{'R2':>9s}")
    out["subdued_replica"] = []
    dmap = collections.defaultdict(list)
    for c, v in zip(cid[keep], zb - elev[cid[keep]]):
        dmap[c].append(v)
    for alpha in (1.0, 0.9, 0.8, 0.4, 0.0):
        pred = md + (1 - alpha) * np.array([np.mean(dmap[c]) for c in cs])
        r2v = float(1 - ((pred - oc) ** 2).sum() / ((oc - oc.mean()) ** 2).sum())
        pr, sp = pearsonr(pred, oc).statistic, spearmanr(pred, oc).statistic
        out["subdued_replica"].append({"alpha": alpha, "sd_m": float(pred.std()),
                                       "pearson": float(pr), "spearman": float(sp),
                                       "r2": r2v})
        print(f"    {alpha:7.2f}{pred.std():10.2f}{pr:+10.4f}{sp:+10.4f}{r2v:+9.3f}")
    print(f"    observed sd {oc.std():.2f} m")

    print(f"\n  DRAIN SCAN: does any terrain drain elevation track observed depth?")
    LO, LA = np.meshgrid(np.radians(dlo), np.radians(dla))
    _, idx = cKDTree(xyz.T).query(
        np.stack([(np.cos(LA) * np.cos(LO)).ravel(), (np.cos(LA) * np.sin(LO)).ravel(),
                  np.sin(LA).ravel()], axis=1), workers=-1)
    zf = z.ravel().astype(np.float64)
    o = np.argsort(idx, kind="stable")
    idxs, zs = idx[o], zf[o]
    bnd = np.searchsorted(idxs, np.arange(n + 1))
    for q in (0, 5, 10, 25, 50):
        zd = np.array([zs[bnd[c]:bnd[c + 1]].min() if q == 0
                       else np.percentile(zs[bnd[c]:bnd[c + 1]], q) for c in cs])
        depth = elev[cs] - zd
        pr, sp = pearsonr(depth, oc).statistic, spearmanr(depth, oc).statistic
        lab = "minimum" if q == 0 else f"p{q}"
        out["drain_scan"].append({"percentile": lab, "median_depth_m": float(np.median(depth)),
                                  "pearson": float(pr), "spearman": float(sp)})
        print(f"    {lab:>10s}  median {np.median(depth):7.2f} m   "
              f"pearson {pr:+.4f}   spearman {sp:+.4f}")
    print(f"    observed median {np.median(oc):.2f} m -- a drain elevation can match the "
          f"SCALE and still\n    anti-correlate with the pattern, which is why GW-23 closed")
    return out



# --- GW-26: the compound topographic index on the Earth mesh, and the score --
#
# The saturated-fraction closure this scores is built for Vesper in
# `build_topographic_index.py`; everything below points the same construction at
# the Earth mesh the solver is already calibrated on, because that is the only
# place a saturated fraction can be told from a well-shaped number. The bar, the
# support, the comparator and the attribution identity are declared in
# `hydrography/config/topographic_index.yaml` and are read from there rather
# than restated.

TI_CFG = ROOT / "hydrography" / "config" / "topographic_index.yaml"


def _auc(score, label) -> float:
    """Area under the ROC curve, with ties taking the rank they earn.

    The closure caps at one, so ties are not a corner case here: every cell the
    cap binds in carries the same score and a tie-blind AUC would order them by
    accident. The Mann-Whitney form with average ranks gives such a pair the 0.5
    it is worth.
    """
    from scipy.stats import rankdata
    label = np.asarray(label, dtype=bool)
    n1 = int(label.sum())
    n0 = int(label.size - n1)
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(np.asarray(score, dtype=np.float64))
    return float((r[label].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def stage_cti(tag: Path, region: str, quiet: bool) -> None:
    """The index per mesh region, both slope arms, cached beside the solve.

    `ln(a / tan beta)` with `a` the upslope contributing area per unit contour
    width, exactly as `build_topographic_index.py` computes it for Vesper: the
    accumulation is `stage_drainage`'s own, on the filled surface and in real
    area, and the contour width is the square root of the cell's own area.

    BOTH SLOPE ARMS, and the plane fit comes from `lib/orogen.py` rather than
    from a second implementation that would look the same. The Voronoi dual is
    what supplies the neighbour structure and the exact cell areas, and it is
    the same dual `stage_drainage` accumulated over, so the index and the
    accumulation sit on one geometry.
    """
    out = region_dir(tag, region) / "earth_cti.npz"
    if out.exists():
        print("  index: cached")
        return
    from orogen import plane_fit_slope_deg

    cfg = yaml.safe_load(TI_CFG.read_text(encoding="utf-8"))
    if str(cfg["index"]["contour_width"]) != "sqrt_cell_area":
        raise SystemExit("index.contour_width is not 'sqrt_cell_area'; this "
                         "harness computes the same index build_topographic_"
                         "index.py does and has no other one")
    floor = float(cfg["index"]["min_tan_slope"])

    xyz = np.load(tag / "earth_xyz.npy")
    n = xyz.shape[1]
    f = np.load(region_dir(tag, region) / "earth_fields.npz")
    d = np.load(region_dir(tag, region) / "earth_drainage.npz")
    elev = np.nan_to_num(f["elev"], nan=0.0)
    land = f["land"]

    ex = EarthExport(EARTH_R_KM, xyz[0], xyz[1], xyz[2],
                     np.full(n, 4 * np.pi * EARTH_R_KM ** 2 / n),
                     np.where(land, LAND, 0).astype(np.int16))
    t = time.time()
    geom = gw.Geometry(ex)
    off, nb = _csr(geom.src, geom.dst, n)
    area_m2 = np.asarray(geom.voronoi_area_m2, dtype=np.float64)
    print(f"  index: Voronoi dual in {time.time()-t:.0f}s")

    pos = np.ascontiguousarray(xyz.T, dtype=np.float64)
    pos /= np.linalg.norm(pos, axis=1, keepdims=True)
    tan_plane = np.tan(np.deg2rad(
        plane_fit_slope_deg(off, nb, pos, elev, EARTH_R_KM * 1000.0)
        .astype(np.float64)))

    # The receiver-drop arm: the fall to the region's own steepest-descent
    # receiver over the distance to it. Zero inside a filled pit by
    # construction, which is one of the two reasons the index needs a floor.
    recv = d["recv"].astype(np.int64)
    ok = recv >= 0
    to = np.where(ok, recv, 0)
    dot = np.clip(np.einsum("ij,ij->i", pos, pos[to]), -1.0, 1.0)
    dist_m = EARTH_R_KM * 1000.0 * np.arccos(dot)
    tan_recv = np.where(ok & (dist_m > 0.0),
                        (elev - elev[to]) / np.maximum(dist_m, 1.0), 0.0)

    a_m = (d["acc_km2"].astype(np.float64) * 1e6) / np.sqrt(area_m2)
    arms = {name: np.log(a_m / np.maximum(t_, floor))
            for name, t_ in (("plane_fit", tan_plane),
                             ("receiver_drop", tan_recv))}
    for name, t_ in (("plane_fit", tan_plane), ("receiver_drop", tan_recv)):
        print(f"    {name:14s} land mean index {arms[name][land].mean():6.2f}   "
              f"on the slope floor over "
              f"{float((np.maximum(t_, 0.0)[land] <= floor).mean()):.1%} of land")
    np.savez_compressed(out, cti_plane_fit=arms["plane_fit"],
                        cti_receiver_drop=arms["receiver_drop"],
                        area_m2=area_m2)


def _score_boxes(xyz, nlat: int, nlon: int):
    """Equal-angle boxes at the configured grid's own row and column counts.

    THE SUPPORT IS DECLARED IN THE CONFIG and this is it in code. `f_sat_max` is
    a rank statistic over the terrain population inside a climate-grid cell, so
    scoring it needs such a population, and the Earth harness has no climate-grid
    export to take one from. Equal-angle boxes at the configured counts give a
    population of the right size; a Gaussian row spacing and an equal-angle one
    differ by far less than the cell, and nothing here depends on where a box
    edge sits.

    ONE COORDINATE SOURCE, which is what `CLAUDE.md` rule 3 is about. The boxes
    are constructed from the mesh's own generators and nothing is matched to an
    axis from anywhere else; the flat cell index is `row * nlon + col`, the
    convention `lib/gridding.py` uses.
    """
    pos = np.asarray(xyz, dtype=np.float64)
    lat = np.degrees(np.arcsin(np.clip(pos[2] / np.linalg.norm(pos, axis=0), -1, 1)))
    lon = np.degrees(np.arctan2(pos[1], pos[0])) % 360.0
    row = np.clip(((90.0 - lat) / (180.0 / nlat)).astype(np.int64), 0, nlat - 1)
    col = np.clip((lon / (360.0 / nlon)).astype(np.int64), 0, nlon - 1)
    return row * nlon + col


def stage_fsat(tag: Path, edge_km: float, region: str, confinement: str | None,
               quiet: bool) -> dict:
    """GW-26's score: does the terrain half of the closure discriminate?

    WHAT IS BEING TESTED, and it is not the closure as a whole. `f_sat =
    min(f_sat_max exp(-f_grad z_wt), 1)` is strictly monotone in the cell-mean
    depth for a fixed `f_sat_max`, so ranking cells by it with a CONSTANT
    `f_sat_max` reproduces ranking them by depth exactly. Any discrimination it
    gains over the depth alone is therefore `f_sat_max`'s and nothing else's,
    which is what makes this a test rather than a comparison. That identity is
    CHECKED here at the tolerance the config declares, because a harness where
    it fails is measuring something other than what it says.

    TWO CONDITIONS, BOTH DECLARED BEFORE ANY FRACTION WAS SCORED: the declared
    bar `bar_auc`, which is the depth's own discrimination per bore against the
    mesh region it falls in, and the cell-mean depth's AUC AT THIS SUPPORT. The
    second is there because a coarser support moves an AUC on its own, so
    clearing the first alone would not show the terrain half had done anything.
    """
    import pandas as pd

    cfg = yaml.safe_load(TI_CFG.read_text(encoding="utf-8"))
    sc = cfg["score"]
    if str(sc["support"]) != "model_grid_equal_angle_boxes":
        raise SystemExit(f"score.support {sc['support']!r} is not the support "
                         "this harness takes; it is "
                         "'model_grid_equal_angle_boxes'")
    if str(sc["metric"]) != "auc_water_table_within_1m":
        raise SystemExit(f"score.metric {sc['metric']!r} is not the metric this "
                         "harness computes")
    if str(cfg["closure"]["f_sat_max"]) != "area_share_above_cell_mean":
        raise SystemExit(
            "closure.f_sat_max is not 'area_share_above_cell_mean'. This "
            "harness scores the share of a cell's land AREA above that cell's "
            "area-weighted mean index; a region COUNT share is a different "
            "quantity and is not available here.")
    bar = float(sc["bar_auc"])
    tol = float(sc["attribution_identity_tolerance"])
    min_regions = int(sc["min_regions_per_cell"])
    f_grads = [float(v) for v in cfg["closure"]["f_grad_bracket_per_m"]]

    planet = yaml.safe_load((ROOT / "config" / "planet.yaml")
                            .read_text(encoding="utf-8"))
    nlat = int(planet["model"]["latitudes"])
    nlon = int(planet["model"]["longitudes"])
    grid = f"{planet['model']['resolution']} ({nlat}x{nlon} equal-angle boxes)"

    import gridding

    R = REGIONS[region]
    rd = region_dir(tag, region)
    xyz = np.load(tag / "earth_xyz.npy")
    fld = np.load(rd / "earth_fields.npz")
    sol = np.load(rd / "earth_solution.npz")
    idx = np.load(rd / "earth_cti.npz")
    land = fld["land"]
    cond = sol["cond"]
    depth = sol["depth"].astype(np.float64)
    area_m2 = idx["area_m2"]

    cell = _score_boxes(xyz, nlat, nlon)
    ncell = nlat * nlon

    # THE CELL-MEAN DEPTH IS THE INTENSIVE REDUCTION over the population that
    # has a depth at all, which is the conducting land the solver ran on. The
    # index's population is the LAND, which is the population `f_sat_max` is a
    # rank statistic over. Two populations, named separately, because a cell's
    # mean depth over its land and its index over its conducting land are
    # different quantities and collapsing them would hide which one moved.
    z_wt, _, z_covered = gridding.cell_mean(cell, ncell, area_m2, depth, cond)
    depth_ledger = gridding.transfer_ledger(cell, ncell, area_m2, cond)

    per_arm = {}
    for arm in sc["slope_arms"]:
        cti = idx[f"cti_{arm}"].astype(np.float64)
        mean, var, count, covered = gridding.cell_moments(
            cell, ncell, area_m2, cti, land)
        f_max, _ = gridding.cell_fraction(
            cell, ncell, area_m2, cti > mean[cell], land)
        have = covered & (count >= min_regions) & z_covered
        per_arm[arm] = {"f_sat_max": f_max, "cti_mean": mean,
                        "cti_sd": np.sqrt(var), "have": have,
                        "count": count}

    obs = pd.read_csv(ROOT / "hydrography" / "data" / "earth_validation" / R["sites"],
                      low_memory=False)
    if "depth_consistent" in obs.columns:
        obs = obs[obs.depth_consistent.astype(bool)]
    if R["score_bbox"]:
        x0, y0, x1, y1 = R["score_bbox"]
        obs = obs[obs.lon.between(x0, x1) & obs.lat.between(y0, y1)]
    if confinement and "confinement" in obs.columns:
        obs = obs[obs.confinement == confinement]
    if obs.empty:
        raise SystemExit(f"no sites left for region={region} confinement={confinement}")

    la = np.radians(obs.lat.to_numpy())
    lo = np.radians(obs.lon.to_numpy())
    p = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)],
                 axis=1)
    _, bore_region = cKDTree(xyz.T).query(p, workers=-1)
    use = cond[bore_region]
    bore_cell = cell[bore_region]
    offered = int(use.sum())

    label_all = obs.wtd_m.to_numpy() <= 1.0
    # The depth's own discrimination at the support the declared bar was
    # measured on: per bore, against the mesh region it falls in. Reported first
    # because reproducing it is what says the harness is wired to the same
    # quantity the bar came from.
    auc_region_depth = _auc(-depth[bore_region][use], label_all[use])

    out = {
        "issue": "GW-26",
        "region": region,
        "confinement": confinement or "all",
        "edge_km": edge_km,
        "support": grid,
        "min_regions_per_cell": min_regions,
        "bar_auc": bar,
        "label": "observed water table within 1 m of the surface",
        "bores_in_window": int(obs.shape[0]),
        "bores_on_conducting_land": offered,
        "base_rate_all_bores": float(label_all[use].mean()),
        "auc_depth_at_mesh_region": auc_region_depth,
        "reduction": {
            "cell_mean_depth": "gridding.cell_mean, INTENSIVE, area-weighted, population sol['cond']",
            "cell_mean_index": "gridding.cell_moments, INTENSIVE, area-weighted, population land",
            "f_sat_max": "gridding.cell_fraction, CATEGORICAL, share of the cell's land AREA above the cell mean",
            "ledger_depth_population": depth_ledger,
        },
        "arms": {},
    }
    print(f"  f_sat score: {offered:,} bores on conducting land, base rate "
          f"{out['base_rate_all_bores']:.1%}")
    print(f"    AUC of the depth at the mesh region  {auc_region_depth:.4f}   "
          f"(the bar {bar:.3f} was measured here)")

    verdicts = []
    for arm, a in per_arm.items():
        keep = use & a["have"][bore_cell]
        lab = label_all[keep]
        z = z_wt[bore_cell[keep]]
        f_max = a["f_sat_max"][bore_cell[keep]]
        auc_cell_depth = _auc(-z, lab)
        # WHAT THE SUPPORT CAN CARRY AT ALL, which is the instrument this score
        # is read on and not a second criterion. Every cell-scale predictor has
        # ONE value per cell, so the best any of them can do is order the cells
        # by their own observed wet rate; that ceiling is computed here from the
        # observations themselves. A ceiling near 0.5 would mean the cells this
        # support cuts do not separate wet bores from dry ones at all, and a
        # miss under it would say nothing about the terrain statistic. It cannot
        # move the verdict and is reported so the verdict can be read.
        rate = np.zeros(ncell)
        nb_ = np.bincount(bore_cell[keep], minlength=ncell)
        np.divide(np.bincount(bore_cell[keep], weights=lab.astype(float),
                              minlength=ncell), nb_, out=rate, where=nb_ > 0)
        auc_ceiling = _auc(rate[bore_cell[keep]], lab)
        arm_out = {
            "bores_scored": int(keep.sum()),
            "auc_ceiling_at_this_support": auc_ceiling,
            "bores_dropped_for_thin_cells": int((use & ~a["have"][bore_cell]).sum()),
            "cells_with_bores": int(np.unique(bore_cell[keep]).size),
            "base_rate": float(lab.mean()) if keep.any() else float("nan"),
            "auc_cell_mean_depth_same_support": auc_cell_depth,
            "f_sat_max_at_scored_bores": {
                str(q): float(np.percentile(f_max, q)) for q in (5, 50, 95)},
            "f_grad": {},
        }
        print(f"    slope arm {arm}: {arm_out['bores_scored']:,} bores in "
              f"{arm_out['cells_with_bores']:,} cells")
        print(f"      AUC of the cell-mean depth at this support "
              f"{auc_cell_depth:.4f}   (ceiling any cell-scale predictor has "
              f"here {auc_ceiling:.4f})")
        for fg in f_grads:
            f_sat = np.minimum(f_max * np.exp(-fg * z), 1.0)
            # THE ATTRIBUTION IDENTITY. Hold f_sat_max constant and the closure
            # is a strictly decreasing function of the depth, so its ranking is
            # the depth's ranking and the two AUCs must agree to round-off.
            ident = np.minimum(1.0 * np.exp(-fg * z), 1.0)
            auc_ident = _auc(ident, lab)
            gap = abs(auc_ident - auc_cell_depth)
            auc = _auc(f_sat, lab)
            capped = float((f_max * np.exp(-fg * z) >= 1.0).mean())
            arm_out["f_grad"][f"{fg:g}"] = {
                "auc_f_sat": auc,
                "gain_over_cell_mean_depth": auc - auc_cell_depth,
                "beats_declared_bar": bool(auc > bar),
                "beats_same_support_depth": bool(auc > auc_cell_depth),
                "attribution_identity_gap": gap,
                "attribution_identity_holds": bool(gap <= tol),
                "cap_binds_fraction_of_bores": capped,
            }
            verdicts.append(auc > bar and auc > auc_cell_depth)
            print(f"      f_grad {fg:g}/m: AUC {auc:.4f}  "
                  f"({auc - auc_cell_depth:+.4f} on the same-support depth)  "
                  f"identity gap {gap:.2e}  cap binds {capped:.1%}")
        out["arms"][arm] = arm_out

    out["passes"] = bool(verdicts) and all(verdicts)
    out["verdict"] = ("f_sat clears the declared bar and beats the cell-mean "
                      "depth at its own support on every arm"
                      if out["passes"] else
                      "f_sat does NOT clear both declared conditions on every "
                      "arm; the index adds nothing to the depth it multiplies "
                      "and the fraction is withdrawn rather than shipped with a "
                      "caveat")
    print(f"    VERDICT: {'PASS' if out['passes'] else 'MISS'} -- {out['verdict']}")
    return out

def provenance() -> dict:
    try:
        rev = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT,
                             capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        rev = None
    return {"generated": datetime.datetime.now(datetime.timezone.utc)
            .isoformat(timespec="seconds"),
            "git_rev": rev, "numpy": np.__version__,
            "sites": str(SITES.relative_to(ROOT))}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--edge-km", type=float, default=None,
                    help="mean mesh edge in km; N follows from Orogen's pi*R/sqrt(N)")
    ap.add_argument("--regions", type=int, default=None,
                    help="region count directly, overriding --edge-km. The original "
                         "baseline is 1736112, which is Vesper's cell area at Earth's "
                         "radius and is what the cached 15.19 km mesh holds")
    ap.add_argument("--et-lambda", type=float, default=None,
                    help="groundwater ET e-folding depth in m, overriding groundwater.yaml. "
                         "Its DECLARED physical bracket is 0.5 to 2.0; going outside that is "
                         "fitting rather than calibrating")
    ap.add_argument("--surface", default="cell-mean",
                    choices=["cell-mean", "dem-at-bore"],
                    help="ground level the depth hangs from; both come from the DEM")
    ap.add_argument("--drain", default="none", choices=["none", "min", "mean-2sd"],
                    help="GW-23 sub-grid drain at its C -> infinity limit: pin the head "
                         "at the cell's own valley floor, from the DEM minimum or mean-2sd")
    ap.add_argument("--drain-min-relief-m", type=float, default=10.0,
                    help="only cells with at least this much within-cell relief vent")
    ap.add_argument("--region", default="aus", choices=sorted(REGIONS),
                    help="which Earth region to score; the mesh is shared between them")
    ap.add_argument("--confinement", default=None,
                    choices=["unconfined", "confined", "mixed"],
                    help="restrict scoring to bores of this aquifer confinement, "
                         "which only the US table labels")
    ap.add_argument("--stage", default="all",
                    choices=["mesh", "fields", "perm", "drainage", "solve", "score",
                             "benchmark", "calibrate", "diagnostics",
                             "cti", "fsat", "all"])
    ap.add_argument("--river-km2", type=float, default=1e4,
                    help="upstream area above which a cell is a river and becomes a "
                         "fixed head at its own elevation (GW-17). 0 disables baselevels")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "hydrography" / "analysis" / "earth_calibration.json")
    args = ap.parse_args()

    if args.regions is None and args.edge_km is None:
        args.regions = 1_736_112          # the declared baseline, matched to Vesper
    if args.river_km2 is not None and args.river_km2 <= 0:
        args.river_km2 = None
    n = args.regions if args.regions is not None else regions_for_edge(args.edge_km)
    edge_km = edge_for_regions(n)
    tag = CACHE / f"edge{edge_km:.2f}km"
    tag.mkdir(parents=True, exist_ok=True)
    print(f"=== Earth calibration: {args.region} at {edge_km:.2f} km mean edge "
          f"({n:,} regions) ===")

    order = ["mesh", "fields", "perm", "drainage", "solve", "score"]
    run = order if args.stage == "all" else [args.stage]
    result = None
    for st in run:
        if st == "mesh":
            stage_mesh(tag, n, args.quiet)
        elif st == "fields":
            stage_fields(tag, args.quiet, args.region)
        elif st == "perm":
            stage_perm(tag, args.quiet, args.region)
        elif st == "drainage":
            stage_drainage(tag, args.quiet, args.region)
        elif st == "solve":
            stage_solve(tag, args.quiet, args.river_km2, args.region,
                        args.drain, args.drain_min_relief_m, args.et_lambda)
        elif st == "benchmark":
            result = stage_benchmark(tag, args.region, args.confinement)
        elif st == "diagnostics":
            result = stage_diagnostics(tag, args.region, args.confinement,
                                       args.river_km2)
        elif st == "calibrate":
            result = stage_calibrate(tag, args.region, args.confinement,
                                     args.drain, args.et_lambda)
        elif st == "score":
            result = stage_score(tag, edge_km, args.quiet, args.region,
                                 args.confinement, args.drain, args.surface,
                                 args.et_lambda)
        elif st == "cti":
            stage_cti(tag, args.region, args.quiet)
        elif st == "fsat":
            # GW-26's score has its own artifact rather than a key in
            # earth_calibration.json: that file is the DEPTH's certification and
            # this is a different quantity, and the license the topographic
            # index step reads has to be findable without knowing which run key
            # it was filed under.
            stage_cti(tag, args.region, args.quiet)
            fs = stage_fsat(tag, edge_km, args.region, args.confinement,
                            args.quiet)
            sp = ROOT / "hydrography" / "analysis" / "topographic_index_score.json"
            pay = json.loads(sp.read_text()) if sp.exists() else {}
            pay.setdefault("provenance", provenance())
            pay["provenance"] = provenance()
            key = f"{args.region}/{args.confinement or 'all'}/{edge_km:.2f}km"
            pay.setdefault("runs", {})[key] = fs
            pay["observation_sets_required"] = yaml.safe_load(
                TI_CFG.read_text(encoding="utf-8"))["score"]["observation_sets"]
            # THE LICENSE IS THE AND OF EVERY REQUIRED SET, never of the ones
            # that happen to be in the file. A score run on one continent and
            # filed alone would otherwise license a consumer on half the bar.
            got = {("australia_depth_consistent"
                    if r["region"] == "aus" else
                    f"{r['region']}_{r['confinement']}"): r["passes"]
                   for r in pay["runs"].values()}
            need = pay["observation_sets_required"]
            pay["scored_sets"] = sorted(got)
            pay["consumers_licensed"] = bool(
                need) and all(got.get(k, False) for k in need)
            sp.write_text(json.dumps(pay, indent=2) + "\n")
            print(f"  wrote {sp}\n  consumers_licensed: "
                  f"{pay['consumers_licensed']}")

    if result is not None:
        payload = {}
        if args.out.exists():
            payload = json.loads(args.out.read_text())
        payload.setdefault("provenance", provenance())
        # EVERY OPTION THAT MOVES THE RESULT IS IN THE KEY. `--et-lambda`
        # changes which solution is scored and was not, so two lambdas wrote
        # over each other under one name. Appended rather than inserted, and
        # only when it is set, so the keys already in the file keep their names.
        key = (f"{args.region}/{args.confinement or 'all'}/"
               f"drain-{args.drain}/{args.surface}/{edge_km:.2f}km"
               + ("" if args.et_lambda is None else f"/lam{args.et_lambda:g}"))
        payload.setdefault("runs", {})[key] = result
        args.out.write_text(json.dumps(payload, indent=2) + "\n")
        print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()
