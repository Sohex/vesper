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
          f"(bar 0.07, statistical fit 0.1435)")
    print(f"    Pearson on cell means       {out['pearson_cell_mean']:+.4f}")
    print(f"    observed sd {out['observed_sd_m']:.2f} m   model sd {out['model_sd_m']:.2f} m")
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
                    choices=["mesh", "fields", "perm", "drainage", "solve", "score", "all"])
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
        elif st == "score":
            result = stage_score(tag, edge_km, args.quiet, args.region,
                                 args.confinement, args.drain, args.surface,
                                 args.et_lambda)

    if result is not None:
        payload = {}
        if args.out.exists():
            payload = json.loads(args.out.read_text())
        payload.setdefault("provenance", provenance())
        key = (f"{args.region}/{args.confinement or 'all'}/"
               f"drain-{args.drain}/{args.surface}/{edge_km:.2f}km")
        payload.setdefault("runs", {})[key] = result
        args.out.write_text(json.dumps(payload, indent=2) + "\n")
        print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()
