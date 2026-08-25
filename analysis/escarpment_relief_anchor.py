#!/usr/bin/env python3
"""What relief does a real Great Escarpment carry over 30 km? WORLD-XGAJ.

    python analysis/escarpment_relief_anchor.py --dem-dir <dir of Copernicus GLO-90 tiles>

This measures EARTH, deliberately and only. Vesper is an invented super-Earth
whose geography comes from a fork of World Orogen, and `computeScarpPotential`
decides where a sub-grid escarpment belongs on it by putting a smoothstep on a
relief estimator. `notes/audits/orogen-resolution.md` settles which estimator
transports between region counts and leaves the two thresholds unsettled,
because choosing them by matching the areal extent the old gate happened to
produce would be calibrating against the defect the change exists to remove.
The only honest anchor is a measured escarpment, and Earth is where one can be
measured. Nothing about Vesper is simulated here.

## The source

Copernicus DEM GLO-90, the 90 m global DSM, tiles fetched from the AWS Registry
of Open Data bucket `copernicus-dem-90m` (public, no credentials). Each tile is
one degree square, EPSG:4326, heights in metres above the EGM2008 geoid.

## The estimator, fixed before anything was measured

The same one `orogen-resolution.md` measured on the two Vesper builds, in the
same units, so the number transfers to the module without a conversion:

    R = max(0, own height - mean(land height within the outer ball)) / outer radius

`R` is a rise ABOVE the neighbourhood mean per kilometre of ground, so it is
non-zero only where the ground stands high relative to what surrounds it. That
is the plateau side of a margin, which is the side the cliff is cut into, and it
is zero on the lowland below. The measure is one-sided by construction and its
median over any large area is zero.

Where an inner length is declared, the own height is itself the mean over a ball
of that radius, so both terms are fixed-length averages and neither carries the
mesh's own support. Both balls are the nearest whole number of support cells to
their declared length, and `R` is divided by the radius that whole number
REALISES rather than by the declared length, which is how the mesh
implementation rounds a hop count.

## The support, and why there are two of it

An escarpment measured at 90 m is not the quantity the module reads: Orogen
reports a mesh region tens of kilometres across. So the DEM is block-averaged
to the two mesh spacings this project has builds at, 15.19 km and 7.60 km,
before the estimator runs. The two answers bracket the support dependence that
`orogen-resolution.md` measured at 1.15x, and the bracket is reported rather
than hidden. DECLARED before the run: the adopted value is the MEAN of the two
supports, rounded to two significant figures, and the bracket is the pair.

## The two domains, fixed before the run by geography and not by relief

Both are boxes on the southern African Great Escarpment and the plateau behind
it, chosen from where the escarpment is known to run rather than from anything
this measurement reports.

  ESCARPMENT   28.80 to 29.80 E, 29.60 to 28.60 S. The Drakensberg front, from
               the Lesotho highlands down into the KwaZulu-Natal midlands;
               Cathedral Peak, Giant's Castle and Sani Pass all sit in it. The
               front occupies of order a tenth of the box's area.
  PLATEAU      25.60 to 26.60 E, 29.00 to 28.00 S. Free State interior plateau,
               about 250 km west of the front, gently undulating high ground
               with no plateau margin in it.

## The rule that turns the two distributions into the two thresholds

Declared before the run, from what the smoothstep's own comments say each
threshold means:

  SCARP_FULL_RELIEF   "at and above this, relief is not the limiting factor".
                      The 90th percentile of `R` over the ESCARPMENT box: the
                      value only the front itself reaches, since the front is
                      about a tenth of that box.
  SCARP_MIN_RELIEF    "below this there is no relief to hang a cliff on". The
                      90th percentile of `R` over the PLATEAU box: the value
                      ordinary non-escarpment ground essentially never reaches,
                      so terrain like that box scores zero through the gate.

## Gravity

The anchor is Earth relief at Earth gravity. Orogen's own heights carry a 1/g
relief scaling, `reliefScale = REFERENCE_GRAVITY_MS2 / gravityMS2`, and the
maintainable relief of the sub-grid cliff the gate is asking about scales the
same way. So the constants are declared as the Earth measurement and multiplied
by `reliefScale` where they are used, which leaves them checkable against this
script and leaves the gate gravity-aware. The Vesper-side value at this
project's gravity is reported here for reference and is NOT what the module
carries.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "lib"))


OUTPUT = PROJECT_ROOT / "analysis" / "escarpment_relief_anchor.json"

# The two lengths come from the module, so this cannot quote an anchor for a
# form the generator does not run. `analysis/scarp_relief_transport.py` is what
# chose them; both are overridable, which is how the candidates were anchored
# before one was adopted.
TERRAIN_CONFIG = (PROJECT_ROOT / "vendor" / "orogen" / "js" / "terrain-config.js")


def module_constant(name: str) -> float:
    m = re.search(rf"export const {name} = ([0-9.eE+-]+);", TERRAIN_CONFIG.read_text())
    if not m:
        raise SystemExit(f"{name} is not in {TERRAIN_CONFIG}")
    return float(m.group(1))


BASELINE_KM = module_constant("SCARP_RELIEF_BASELINE_KM")
INNER_KM = module_constant("SCARP_RELIEF_INNER_KM")
SUPPORTS_KM = (15.19, 7.60)
DOMAINS = {
    "escarpment": {"lon": (28.80, 29.80), "lat": (-29.60, -28.60),
                   "name": "Drakensberg front, KwaZulu-Natal and Lesotho"},
    "plateau": {"lon": (25.60, 26.60), "lat": (-29.00, -28.00),
                "name": "Free State interior plateau"},
}
QUANTILE_FULL = 0.90     # over the escarpment box
QUANTILE_MIN = 0.90      # over the plateau box
REPORT_QUANTILES = (0.50, 0.75, 0.90, 0.95, 0.99)

# The local equal-area projection the block averaging is done in. Southern
# African Albers is not used: an azimuthal equal-area centred on the domain
# keeps the two boxes on the same footing and needs no standard parallels.
PROJ = ("+proj=laea +lat_0=-29 +lon_0=27.7 +x_0=0 +y_0=0 "
        "+datum=WGS84 +units=m +no_defs")


def build_support(dem_dir: Path, work: Path, support_km: float) -> Path:
    """Block-average every tile to one equal-area raster at `support_km`."""
    vrt = work / "tiles.vrt"
    if not vrt.exists():
        tiles = sorted(str(p) for p in dem_dir.glob("*.tif"))
        if not tiles:
            raise SystemExit(f"no GLO-90 tiles in {dem_dir}")
        subprocess.run(["gdalbuildvrt", "-q", str(vrt), *tiles], check=True)
    out = work / f"support_{support_km:g}km.tif"
    if not out.exists():
        metres = support_km * 1000.0
        subprocess.run(
            ["gdalwarp", "-q", "-overwrite", "-t_srs", PROJ, "-r", "average",
             "-tr", str(metres), str(metres), "-dstnodata", "-9999",
             str(vrt), str(out)], check=True)
    return out


def disc_kernel(radius_cells: int) -> np.ndarray:
    r = radius_cells
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    return ((x * x + y * y) <= r * r).astype(np.float64)


def relief_field(raster: Path, support_km: float, baseline_km: float,
                 inner_km: float):
    """`R` over the whole raster, plus the realised outer radius in km."""
    import rasterio
    from scipy.signal import fftconvolve

    with rasterio.open(raster) as src:
        band = src.read(1).astype(np.float64)
        nodata = src.nodata
        transform = src.transform
        crs = src.crs
    valid = np.isfinite(band) & (band != nodata) if nodata is not None else np.isfinite(band)
    # Sea level is a value, not a gap: GLO-90 carries zeros over water, and both
    # boxes are inland, so nothing is discarded on that account.
    height = np.where(valid, band, 0.0)
    good = valid.astype(np.float64)

    def ball_mean(radius_km: float):
        cells = max(1, int(round(radius_km / support_km)))
        kernel = disc_kernel(cells)
        total = fftconvolve(height, kernel, mode="same")
        count = fftconvolve(good, kernel, mode="same")
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = np.where(count > 0.5, total / np.maximum(count, 1e-9), np.nan)
        return mean, cells * support_km

    outer, realised_km = ball_mean(baseline_km)
    if inner_km > 0:
        own, _ = ball_mean(inner_km)
    else:
        own = band
    relief = np.maximum(0.0, own - outer) / (realised_km * 1000.0)
    relief = np.where(valid & np.isfinite(outer) & np.isfinite(own), relief, np.nan)
    return relief, realised_km, transform, crs, valid


def domain_mask(shape, transform, crs, lon_range, lat_range) -> np.ndarray:
    from pyproj import Transformer

    rows, cols = shape
    j, i = np.meshgrid(np.arange(cols), np.arange(rows))
    x = transform.c + (j + 0.5) * transform.a + (i + 0.5) * transform.b
    y = transform.f + (j + 0.5) * transform.d + (i + 0.5) * transform.e
    to_wgs = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    lon, lat = to_wgs.transform(x, y)
    return ((lon >= lon_range[0]) & (lon <= lon_range[1])
            & (lat >= lat_range[0]) & (lat <= lat_range[1]))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dem-dir", type=Path, required=True,
                    help="directory of Copernicus GLO-90 GeoTIFF tiles")
    ap.add_argument("--work", type=Path, default=None,
                    help="scratch directory for the mosaic (default: --dem-dir)")
    ap.add_argument("--baseline-km", type=float, default=BASELINE_KM,
                    help="outer ball radius the rise is quoted over")
    ap.add_argument("--inner-km", type=float, default=INNER_KM,
                    help="inner ball radius the own height is averaged over; "
                         "0 takes the support cell itself")
    args = ap.parse_args()
    work = args.work or args.dem_dir
    work.mkdir(parents=True, exist_ok=True)

    # The relief scaling is read from the build rather than recomputed from a
    # reference gravity written down here: `manifest.planet.reliefScale` is the
    # factor Orogen actually applied to the heights the gate reads.
    from orogen import Export
    import builds as builds_lib
    export = Export(builds_lib.mesh_export())
    gravity = export.gravity_m_s2
    relief_scale = float(export.relief_scale)

    per_support = {}
    for support_km in SUPPORTS_KM:
        raster = build_support(args.dem_dir, work, support_km)
        relief, realised_km, transform, crs, valid = relief_field(
            raster, support_km, args.baseline_km, args.inner_km)
        entry = {"realised_ball_radius_km": realised_km, "domains": {}}
        for key, spec in DOMAINS.items():
            mask = domain_mask(relief.shape, transform, crs, spec["lon"], spec["lat"])
            values = relief[mask & np.isfinite(relief)]
            entry["domains"][key] = {
                "name": spec["name"],
                "cells": int(values.size),
                "quantiles": {f"p{int(q * 100)}": float(np.quantile(values, q))
                              for q in REPORT_QUANTILES},
                "max": float(values.max()),
            }
        per_support[f"{support_km:g}km"] = entry

    full_by_support = [per_support[f"{s:g}km"]["domains"]["escarpment"]
                       ["quantiles"][f"p{int(QUANTILE_FULL * 100)}"] for s in SUPPORTS_KM]
    min_by_support = [per_support[f"{s:g}km"]["domains"]["plateau"]
                      ["quantiles"][f"p{int(QUANTILE_MIN * 100)}"] for s in SUPPORTS_KM]

    def two_sig(x: float) -> float:
        if x == 0:
            return 0.0
        from math import floor, log10
        return round(x, -int(floor(log10(abs(x)))) + 1)

    adopted_full = two_sig(float(np.mean(full_by_support)))
    adopted_min = two_sig(float(np.mean(min_by_support)))

    out = {
        "measured": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "source": "Copernicus DEM GLO-90, AWS Open Data bucket copernicus-dem-90m",
        "baseline_km": args.baseline_km,
        "inner_km": args.inner_km,
        "supports_km": list(SUPPORTS_KM),
        "domains": {k: {"lon": v["lon"], "lat": v["lat"], "name": v["name"]}
                    for k, v in DOMAINS.items()},
        "rule": {
            "full": f"p{int(QUANTILE_FULL * 100)} of R over the escarpment box",
            "min": f"p{int(QUANTILE_MIN * 100)} of R over the plateau box",
            "adopt": "mean of the two supports, two significant figures",
        },
        "per_support": per_support,
        "earth": {
            "scarp_full_relief": adopted_full,
            "scarp_full_bracket": [min(full_by_support), max(full_by_support)],
            "scarp_min_relief": adopted_min,
            "scarp_min_bracket": [min(min_by_support), max(min_by_support)],
        },
        "vesper_reference_only": {
            "gravity_m_s2": gravity,
            "relief_scale": relief_scale,
            "scarp_full_relief": adopted_full * relief_scale,
            "scarp_min_relief": adopted_min * relief_scale,
        },
    }
    OUTPUT.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")

    for support_km in SUPPORTS_KM:
        entry = per_support[f"{support_km:g}km"]
        print(f"support {support_km:g} km, ball radius {entry['realised_ball_radius_km']:.2f} km")
        for key, dom in entry["domains"].items():
            q = dom["quantiles"]
            print(f"  {key:<11}{dom['cells']:>7} cells  " +
                  "  ".join(f"{k} {v:.5f}" for k, v in q.items()))
    print(f"Earth anchor: min {adopted_min:.4f} (bracket {min(min_by_support):.5f}"
          f" to {max(min_by_support):.5f}), full {adopted_full:.4f} "
          f"(bracket {min(full_by_support):.5f} to {max(full_by_support):.5f})")
    print(f"wrote {OUTPUT.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
