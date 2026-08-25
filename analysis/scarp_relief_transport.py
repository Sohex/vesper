#!/usr/bin/env python3
"""Which scarp relief estimator survives a change of region count? WORLD-XGAJ.

    python analysis/scarp_relief_transport.py
    python analysis/scarp_relief_transport.py --control <dir>/exoplasim-T42

Worldbuilding. Vesper is an invented super-Earth; every quantity here is a
modelled field of that planet, measured on World Orogen exports of the same
planet at the same seed and different region counts. Nothing is run: no climate
model, no soil solve, no biosphere.

`vendor/orogen/js/lithology.js:computeScarpPotential` gates a scarp marker on a
relief estimator with an absolute smoothstep, so the estimator has to report the
same landform the same way at any region count. This measures which one does,
reproduces the published shift of the one-edge form as a control on the harness,
and sizes the realisation noise floor that says which of these ratios can be
believed at all.

## The estimator family, and the definitions

For a land region, walk the LAND subgraph to a whole number of hops and take the
plain mean of the physical land height over that ball:

    R = max(0, own - mean(height over the outer ball)) / (outer hops * avgEdgeKm)

where `own` is the region's own height, or the mean over an INNER ball where an
inner length is declared. Positive only where the ground stands above its
neighbourhood, which is the plateau side of a margin and the side a cliff is cut
into; the lowland below it scores zero, so the measure is one-sided and its
median over land is exactly zero. Both balls are quoted over the run their whole
hop count REALISES, and on this project's two builds every candidate below
realises the same two radii to within 0.1 per cent, because 15.19 km is twice
7.60 km.

The land test is `surface_class`, never `land_mask` (CLAUDE.md rule 1): a dry
closed-basin floor below sea level is land, and a plateau margin standing inside
one is an escarpment like any other.

## The candidates and the bar, declared before anything was run

Outer baseline over 30, 45, 60 and 90 km; inner length 0 km, meaning the region
itself, and 15 km, which is the coarse build's own spacing and so the shortest
inner ball both meshes can express. Below about 20 km Orogen designs no terrain
(`notes/audits/orogen-resolution.md`), so a shorter baseline than 30 km would
quote a rise over ground the generator never shaped.

  BAR         a candidate TRANSPORTS when every land quantile from p50 to p99
              agrees between the 2.5M and the 10M build within 1.15x. This is
              the bar `orogen-resolution.md` fixed for this class and it is not
              re-chosen here. The measure is one-sided, so p50 is zero on both
              builds and is reported as undefined rather than as a pass.

  INSTRUMENT  changing `--regions` also moves the first plate seed, so two
              region counts are two REALISATIONS of the same world and not only
              two resolutions. The control is a third build 4 per cent away from
              the coarse one, where the resolution is to all intents identical
              and anything that moves is scatter. A support ratio inside that
              scatter says nothing, and a bar tighter than the scatter cannot be
              judged at all.

  ADOPT       the SHORTEST outer baseline that transports, preferring the inner
              length 0 at equal baseline because it declares one length rather
              than two. If nothing transports but the noise floor swallows the
              bar, that is the finding and no swap is made on this evidence.

## The control on the harness

The one-edge form -- the steepest drop to a lower land neighbour over that edge
-- is measured beside the family. `orogen-resolution.md` publishes its shift as
1.408, 1.412, 1.537, 1.610, 1.695 at p50 to p99, so reproducing those numbers is
what says this harness measures what the audit measured.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy.sparse as sp

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "lib"))

from orogen import Export, LAND                                # noqa: E402

OUTPUT = PROJECT_ROOT / "analysis" / "scarp_relief_transport.json"

BUILDS = ("precarve-craton", "precarve-craton-10m")
OUTER_KM = (30.0, 45.0, 60.0, 90.0)
INNER_KM = (0.0, 15.0)
SAMPLE = 200_000
SEED = 16236323
QUANTILES = (0.50, 0.75, 0.90, 0.95, 0.99)
TRANSPORT_BAR = 1.15
CHUNK = 50_000


class Mesh:
    """One export, with the land subgraph and a seeded land sample."""

    def __init__(self, root: Path, known: bool = True):
        self.export = Export(root, require_known_build=known)
        self.land = self.export.surface_class == LAND
        self.height = self.export.elevation_km.astype(np.float64)
        self.spacing = float(
            self.export.manifest["basins"]["resolution"]["avgEdgeKm"])
        n = self.export.n_regions
        off, lst = self.export.adjacency
        src = np.repeat(np.arange(n, dtype=np.int64),
                        np.diff(off.astype(np.int64)))
        dst = lst.astype(np.int64)
        keep = self.land[src] & self.land[dst]
        a = sp.csr_matrix((np.ones(int(keep.sum()), dtype=np.float32),
                           (src[keep], dst[keep])), shape=(n, n))
        self.step = (a + sp.identity(n, dtype=np.float32, format="csr")).tocsr()
        rng = np.random.default_rng(SEED)
        ids = np.flatnonzero(self.land)
        self.sample = rng.choice(ids, size=min(SAMPLE, ids.size), replace=False)
        self._means = {0: self.height[self.sample]}

    def hops(self, radius_km: float) -> int:
        return 0 if radius_km <= 0 else max(1, int(round(radius_km / self.spacing)))

    def ball_mean(self, hops: int) -> np.ndarray:
        if hops in self._means:
            return self._means[hops]
        out = np.empty(self.sample.size)
        for a in range(0, self.sample.size, CHUNK):
            reach = self.step[self.sample[a:a + CHUNK]]
            for _ in range(hops - 1):
                reach = reach @ self.step
                reach.data[:] = 1.0
            reach.data[:] = 1.0
            count = np.asarray(reach.sum(axis=1)).ravel()
            out[a:a + CHUNK] = (reach @ self.height) / np.maximum(count, 1.0)
        self._means[hops] = out
        return out

    def relief(self, outer_km: float, inner_km: float) -> np.ndarray:
        ho, hi = self.hops(outer_km), self.hops(inner_km)
        rise = self.ball_mean(hi) - self.ball_mean(ho)
        return np.maximum(0.0, rise) / (ho * self.spacing)

    def one_edge_gradient(self) -> np.ndarray:
        """The form the module used before, as the harness control."""
        off, lst = self.export.adjacency
        xyz = np.stack([self.export.x, self.export.y, self.export.z],
                       axis=1).astype(np.float64)
        radius = self.export.radius_km
        out = np.zeros(self.sample.size)
        for i, r in enumerate(self.sample):
            nb = lst[off[r]:off[r + 1]].astype(np.int64)
            nb = nb[self.land[nb] & (self.height[nb] < self.height[r])]
            if nb.size == 0:
                continue
            d = np.linalg.norm(xyz[nb] - xyz[r], axis=1) * radius
            out[i] = np.max((self.height[r] - self.height[nb]) / np.maximum(d, 1e-9))
        return out


def quantiles(values: np.ndarray) -> dict[str, float]:
    return {f"p{int(q * 100)}": float(np.quantile(values, q)) for q in QUANTILES}


def ratio(coarse: dict[str, float], fine: dict[str, float]):
    out = {}
    for k, a in coarse.items():
        b = fine[k]
        out[k] = None if (a == 0.0 and b == 0.0) else (float(b / a) if a else None)
    return out


def within(bar: float, ratios: dict[str, float | None]) -> bool:
    vals = [v for v in ratios.values() if v is not None]
    return bool(vals) and all(abs(np.log(v)) <= np.log(bar) for v in vals)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--control", type=Path, default=None,
                    help="a T42 export a few per cent from the coarse build's "
                         "region count, for the realisation noise floor")
    args = ap.parse_args()

    coarse = Mesh(PROJECT_ROOT / "source" / BUILDS[0] / "exoplasim-T42")
    fine = Mesh(PROJECT_ROOT / "source" / BUILDS[1] / "exoplasim-T42")
    control = Mesh(args.control, known=False) if args.control else None

    candidates = {}
    for outer in OUTER_KM:
        for inner in INNER_KM:
            key = f"outer{outer:g}_inner{inner:g}"
            qc = quantiles(coarse.relief(outer, inner))
            qf = quantiles(fine.relief(outer, inner))
            entry = {
                "outer_km": outer, "inner_km": inner,
                "realised_outer_km": {
                    "coarse": coarse.hops(outer) * coarse.spacing,
                    "fine": fine.hops(outer) * fine.spacing},
                "coarse": qc, "fine": qf,
                "support_ratio": ratio(qc, qf),
            }
            if control is not None:
                qn = quantiles(control.relief(outer, inner))
                entry["control"] = qn
                entry["noise_ratio"] = ratio(qc, qn)
            candidates[key] = entry

    edge_c, edge_f = quantiles(coarse.one_edge_gradient()), quantiles(fine.one_edge_gradient())
    harness = {"coarse": edge_c, "fine": edge_f, "support_ratio": ratio(edge_c, edge_f),
               "published": {"p50": 1.408, "p75": 1.412, "p90": 1.537,
                             "p95": 1.610, "p99": 1.695}}
    if control is not None:
        harness["noise_ratio"] = ratio(edge_c, quantiles(control.one_edge_gradient()))

    for key, entry in candidates.items():
        entry["transports"] = within(TRANSPORT_BAR, entry["support_ratio"])
        if "noise_ratio" in entry:
            floor = max(abs(np.log(v)) for v in entry["noise_ratio"].values()
                        if v is not None)
            entry["noise_floor_x"] = float(np.exp(floor))
            entry["support_clears_noise"] = not within(
                float(np.exp(floor)), entry["support_ratio"])

    adopted = None
    for outer in OUTER_KM:
        for inner in INNER_KM:
            key = f"outer{outer:g}_inner{inner:g}"
            if candidates[key]["transports"]:
                adopted = key
                break
        if adopted:
            break

    out = {
        "measured": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "bar": TRANSPORT_BAR,
        "sample": SAMPLE,
        "seed": SEED,
        "builds": {
            "coarse": {"regions": coarse.export.n_regions,
                       "terrain_hash": coarse.export.terrain_hash[:8],
                       "spacing_km": coarse.spacing},
            "fine": {"regions": fine.export.n_regions,
                     "terrain_hash": fine.export.terrain_hash[:8],
                     "spacing_km": fine.spacing},
        },
        "harness_control_one_edge": harness,
        "candidates": candidates,
        "adopted": adopted,
    }
    if control is not None:
        out["builds"]["control"] = {"regions": control.export.n_regions,
                                    "terrain_hash": control.export.terrain_hash[:8],
                                    "spacing_km": control.spacing}
    OUTPUT.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")

    def row(d):
        return "".join("      n/a" if d[f"p{int(q * 100)}"] is None
                       else f"{d[f'p{int(q * 100)}']:>9.3f}" for q in QUANTILES)

    print("  " + " " * 26 + "".join(f"{'p' + str(int(q * 100)):>9}" for q in QUANTILES))
    print(f"  {'one edge, support':<26}{row(harness['support_ratio'])}")
    if "noise_ratio" in harness:
        print(f"  {'one edge, realisation':<26}{row(harness['noise_ratio'])}")
    for key, entry in candidates.items():
        print(f"  {key + ', support':<26}{row(entry['support_ratio'])}"
              f"  {'PASS' if entry['transports'] else 'fail'}")
        if "noise_ratio" in entry:
            print(f"  {key + ', realisation':<26}{row(entry['noise_ratio'])}"
                  f"  floor {entry['noise_floor_x']:.3f}x")
    print(f"  adopted: {adopted}")
    print(f"wrote {OUTPUT.relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
