"""Render the Vesper base map into five projections.

Where a projection has a free orientation, it is chosen against the geography
rather than assumed: central meridians, the butterfly's cut meridians and the
icosahedron's attitude are all searched for the arrangement that puts the
interruption over the most water. Fuller picked his net so Earth's continents
came out whole; the same argument applied to this world gives a different net,
so the net is derived here instead of copied.

    python maps/render_projections.py [--width 2600] [--only winkel_tripel]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import binary_erosion, map_coordinates

sys.path.insert(0, str(Path(__file__).resolve().parent))
import projections as P  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BUILD = Path(__file__).resolve().parent / "build"
OUT = Path(__file__).resolve().parent

BACKGROUND = np.array([13, 16, 21.0])
GRATICULE = np.array([255, 255, 255.0])
OUTLINE = np.array([196, 204, 214.0])


# --------------------------------------------------------------------------
# Sampling the base map


class BaseMap:
    def __init__(self):
        rgb = np.load(BUILD / "basemap.npy")
        self.h, self.w = rgb.shape[:2]
        # One wrapped column at each end so bilinear sampling crosses the seam.
        self.padded = np.concatenate([rgb[:, -1:], rgb, rgb[:, :1]], axis=1).astype(np.float32)
        self.land = np.load(BUILD / "land.npy")

    def sample(self, lon, lat):
        col = np.mod((lon + 180.0) / 360.0 * self.w - 0.5, self.w) + 1.0
        row = np.clip((90.0 - lat) / 180.0 * self.h - 0.5, 0, self.h - 1)
        out = np.empty(lon.shape + (3,), np.float32)
        for c in range(3):
            out[..., c] = map_coordinates(self.padded[..., c], [row, col], order=1, mode="nearest")
        return out

    def land_at(self, lon, lat):
        col = np.mod(np.round((lon + 180.0) / 360.0 * self.w - 0.5).astype(int), self.w)
        row = np.clip(np.round((90.0 - lat) / 180.0 * self.h - 0.5).astype(int), 0, self.h - 1)
        return self.land[row, col]


def arc_land_fraction(base, a, b, n=400):
    """Land fraction along the great-circle arc between two unit vectors."""
    t = np.linspace(0.0, 1.0, n)[:, None]
    v = a[None, :] * (1 - t) + b[None, :] * t
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    lon, lat = P.to_lonlat(v)
    return base.land_at(lon, lat).mean()


def meridian_land_fraction(base, lon0, n=720):
    lat = np.linspace(-89.9, 89.9, n)
    return base.land_at(np.full(n, lon0), lat).mean()


def parallel_land_fraction(base, lat0, n=1440):
    lon = np.linspace(-180.0, 180.0, n, endpoint=False)
    return base.land_at(lon, np.full(n, lat0)).mean()


# --------------------------------------------------------------------------
# Orientation searches


def best_central_meridian(base, step=1.0):
    """Central meridian whose antimeridian crosses the least land."""
    cand = np.arange(-180.0, 180.0, step)
    cost = np.array([meridian_land_fraction(base, P.wrap180(c + 180.0)) for c in cand])
    return float(cand[np.argmin(cost)]), float(cost.min())


def best_hemisphere_meridian(base, step=1.0):
    """Centre longitude whose bounding meridian pair crosses the least land."""
    cand = np.arange(-180.0, 180.0, step)
    cost = np.array(
        [
            meridian_land_fraction(base, P.wrap180(c + 90.0))
            + meridian_land_fraction(base, P.wrap180(c - 90.0))
            for c in cand
        ]
    )
    return float(cand[np.argmin(cost)]), float(cost.min() / 2)


def best_butterfly(base, step=0.5):
    """Equatorial vertex longitude, scored on the two cut meridians only.

    The wings are continuous across their own equatorial and polar interior
    edges; only the meridian pair 90 degrees either side of the fan centres is
    actually torn.
    """
    best = None
    for a in np.arange(0.0, 90.0, step):
        for join in range(4):
            eq = [a + 90.0 * k for k in range(4)]
            cut_lons = [eq[(join + 2) % 4], eq[join]]
            # The join meridian is torn too: the wings meet only at a point.
            cost = sum(meridian_land_fraction(base, P.wrap180(m)) for m in cut_lons)
            if best is None or cost < best[0]:
                best = (cost, float(a), join)
    return best[1], best[2], best[0] / 2


def random_quaternions(n, seed=20240815):
    rng = np.random.default_rng(seed)
    q = rng.normal(size=(n, 4))
    return q / np.linalg.norm(q, axis=1, keepdims=True)


def icosahedron_edges(faces):
    edges = set()
    for f in faces:
        for i in range(3):
            edges.add(tuple(sorted((f[i], f[(i + 1) % 3]))))
    return sorted(edges)


def best_icosahedron_attitudes(base, trials=1500, keep=8):
    """Attitudes whose edge skeleton lies over the most sea.

    Every cut the net can make runs along one of the 30 edges, so this is the
    choice that decides how much land the net is able to keep whole. The
    shortlist is handed to the net search, which knows which edges actually
    get cut."""
    verts, faces = P.icosahedron()
    edges = icosahedron_edges(faces)
    scored = []
    for q in random_quaternions(trials):
        v = P.rotate(verts, q)
        cost = sum(arc_land_fraction(base, v[i], v[j], n=90) for i, j in edges)
        scored.append((cost / len(edges), q))
    scored.sort(key=lambda s: s[0])
    return scored[:keep]


def face_adjacency(faces):
    """(face, face) -> shared edge, for every pair sharing two vertices."""
    adj = {}
    for i, fi in enumerate(faces):
        for j, fj in enumerate(faces):
            if j <= i:
                continue
            shared = sorted(set(fi) & set(fj))
            if len(shared) == 2:
                adj[(i, j)] = tuple(shared)
    return adj


def _kruskal_tree(scored, n_faces):
    parent = list(range(n_faces))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    tree = {f: [] for f in range(n_faces)}
    cut = []
    for w, i, j in scored:
        ri, rj = find(i), find(j)
        if ri == rj:
            cut.append(w)
            continue
        parent[ri] = rj
        tree[i].append(j)
        tree[j].append(i)

    tree_parent = [-1] * n_faces
    seen, stack = {0}, [0]
    while stack:
        f = stack.pop()
        for g in tree[f]:
            if g not in seen:
                seen.add(g)
                tree_parent[g] = f
                stack.append(g)
    return tree_parent, cut


def _best_rotation(faces2d, steps=90):
    """Rotate the net to the angle with the tightest bounding box, landscape."""
    pts = np.concatenate(faces2d)
    best = None
    for a in np.linspace(0, np.pi, steps, endpoint=False):
        c, s = np.cos(a), np.sin(a)
        r = pts @ np.array([[c, -s], [s, c]]).T
        area = np.ptp(r[:, 0]) * np.ptp(r[:, 1])
        if best is None or area < best[0]:
            best = (area, a, np.ptp(r[:, 0]), np.ptp(r[:, 1]))
    _, angle, w, h = best
    if h > w:
        angle += np.pi / 2
    c, s = np.cos(angle), np.sin(angle)
    m = np.array([[c, -s], [s, c]])
    return [f @ m.T for f in faces2d]


def _shape(faces2d):
    """Bounding-box fill ratio and aspect of a laid-out net."""
    pts = np.concatenate(faces2d)
    w, h = np.ptp(pts[:, 0]), np.ptp(pts[:, 1])
    bbox = w * h
    def area(f):
        u, v = f[1] - f[0], f[2] - f[0]
        return abs(u[0] * v[1] - u[1] * v[0]) / 2.0

    tri = sum(area(f) for f in faces2d)
    return tri / bbox, w / h


def choose_net(base, verts, faces, candidates=240, compactness=1.0):
    """Spanning tree of the face graph: tear little land, and lie compact.

    Kruskal on the shared edges, heaviest first, where the weight is the land
    fraction along that edge; an edge kept inside the tree is an edge the map
    stays continuous across. Weight alone gives a net that sprawls diagonally,
    so candidates are generated by jittering the ordering and scored on how
    much land they tear and how much of their bounding box they fill. A net
    that self-overlaps cannot lie flat and is discarded.
    """
    adj = face_adjacency(faces)
    weights = [
        (arc_land_fraction(base, verts[a], verts[b], n=200), i, j)
        for (i, j), (a, b) in adj.items()
    ]
    best = None
    for attempt in range(candidates):
        if attempt == 0:
            order = sorted(weights, reverse=True)
        else:
            rng = np.random.default_rng(attempt)
            sigma = 0.03 + 0.30 * (attempt % 7) / 6.0
            order = sorted(weights, key=lambda s: -(s[0] + rng.normal(0, sigma)))
        tree_parent, cut = _kruskal_tree(order, len(faces))
        faces2d, faces3d = P.unfold(verts, faces, tree_parent, 0)
        if P.net_overlaps(faces2d):
            continue
        faces2d = _best_rotation(faces2d)
        torn = float(np.mean(cut))
        fill, aspect = _shape(faces2d)
        # A net that fills its box can still be a 4:1 ribbon, which is a worse
        # page than Fuller's, so penalise anything past about 2:1.
        score = (
            2.5 * torn
            + compactness * (1.0 - fill)
            + 0.6 * max(0.0, aspect / 2.0 - 1.0)
        )
        if best is None or score < best[0]:
            best = (score, faces2d, faces3d, torn, fill)
    if best is None:
        raise RuntimeError("no non-overlapping net found")
    return best[1], best[2], best[3], best[4]


# --------------------------------------------------------------------------
# Rendering


def render(proj, base, width, supersample=2, graticule=30.0):
    x0, x1, y0, y1 = proj.extent
    height = int(round(width * (y1 - y0) / (x1 - x0)))
    W, H = width * supersample, height * supersample

    xs = x0 + (np.arange(W) + 0.5) / W * (x1 - x0)
    ys = y1 - (np.arange(H) + 0.5) / H * (y1 - y0)

    acc = np.zeros((height, width, 3), np.float64)
    block = max(supersample, (1 << 22) // W // supersample * supersample)
    for r0 in range(0, H, block):
        r1 = min(r0 + block, H)
        X, Y = np.meshgrid(xs, ys[r0:r1])
        lon, lat, valid = proj.inverse(X, Y)
        rgb = base.sample(lon, lat).astype(np.float64)
        rgb[~valid] = BACKGROUND

        if graticule:
            rgb = draw_graticule(rgb, lon, lat, valid, graticule, supersample)
        rgb = draw_outline(rgb, valid)

        rows = (r1 - r0) // supersample
        acc[r0 // supersample : r0 // supersample + rows] = (
            rgb.reshape(rows, supersample, width, supersample, 3).mean(axis=(1, 3))
        )
    return np.clip(acc, 0, 255).astype(np.uint8)


def _line_mask(field, spacing, valid, width_px, wrap=None):
    """Pixels within half a line width of a multiple of `spacing`."""
    d = np.abs(np.remainder(field + spacing / 2.0, spacing) - spacing / 2.0)
    gy = np.gradient(field, axis=0)
    gx = np.gradient(field, axis=1)
    if wrap:
        gy = np.remainder(gy + wrap / 2, wrap) - wrap / 2
        gx = np.remainder(gx + wrap / 2, wrap) - wrap / 2
    grad = np.hypot(gx, gy)
    # A face boundary or the antimeridian shows up as a huge gradient; a line
    # drawn there would be a seam artefact, not a parallel.
    sane = grad < np.maximum(5 * np.median(grad[valid]) if valid.any() else 1.0, 1e-6)
    return (d < 0.5 * width_px * np.maximum(grad, 1e-9)) & valid & sane


def draw_graticule(rgb, lon, lat, valid, spacing, supersample):
    lw = 1.1 * supersample
    par = _line_mask(lat, spacing, valid, lw)
    mer = _line_mask(lon, spacing, valid, lw, wrap=360.0)
    eq = _line_mask(lat, 180.0, valid, 1.5 * supersample)
    rgb[par | mer] = rgb[par | mer] * 0.84 + GRATICULE * 0.16
    rgb[eq] = rgb[eq] * 0.74 + GRATICULE * 0.26
    return rgb


def draw_outline(rgb, valid, thickness=1):
    edge = valid & ~binary_erosion(valid, iterations=thickness, border_value=1)
    rgb[edge] = rgb[edge] * 0.45 + OUTLINE * 0.55
    return rgb


# --------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--width", type=int, default=2600)
    ap.add_argument("--supersample", type=int, default=2)
    ap.add_argument("--only", default=None)
    ap.add_argument("--trials", type=int, default=1500)
    args = ap.parse_args()

    base = BaseMap()
    print(f"base map {base.w}x{base.h}")

    prov = {}
    lon0, cut = best_central_meridian(base)
    print(f"central meridian {lon0:+.0f} (antimeridian is {cut * 100:.0f}% land)")
    prov["central_meridian"] = {"lon0": lon0, "antimeridian_land_fraction": cut}

    hemi_lon, hemi_cut = best_hemisphere_meridian(base)
    print(f"hemisphere centre {hemi_lon:+.0f} (bounding meridians {hemi_cut * 100:.0f}% land)")
    prov["hemispheres"] = {"lon0": hemi_lon, "cut_land_fraction": hemi_cut}

    # The polar plate has no orientation to search: its cut is the equator, and
    # the only free parameter is which meridian points down. It reuses the
    # central meridian so the set of maps stays consistent with itself.
    equator_land = parallel_land_fraction(base, 0.0)
    print(f"polar plate cuts the equator ({equator_land * 100:.0f}% land)")
    prov["polar_hemispheres"] = {
        "lon0": lon0,
        "cut": "equator",
        "cut_land_fraction": equator_land,
    }

    bfly_lon, bfly_join, bfly_cut = best_butterfly(base)
    print(f"butterfly vertices at {bfly_lon:+.1f}+90k, join {bfly_join} "
          f"(cut meridians {bfly_cut * 100:.0f}% land)")
    prov["butterfly"] = {
        "equatorial_vertex_lon": bfly_lon,
        "join_index": bfly_join,
        "cut_land_fraction": bfly_cut,
    }

    verts, faces = P.icosahedron()
    best = None
    for edge_land, q in best_icosahedron_attitudes(base, trials=args.trials):
        v = P.rotate(verts, q)
        f2, f3, torn, fill = choose_net(base, v, faces)
        score = 2.5 * torn + (1.0 - fill)
        if best is None or score < best[0]:
            best = (score, q, v, f2, f3, torn, fill, edge_land)
    _, quat, verts_r, faces2d, faces3d, torn, fill, edge_land = best
    print(f"icosahedron edges {edge_land * 100:.0f}% land; net tears "
          f"{torn * 100:.0f}% land per cut edge, fills {fill * 100:.0f}% of its box")
    prov["icosahedron"] = {
        "quaternion": [float(q) for q in quat],
        "mean_edge_land_fraction": float(edge_land),
        "mean_cut_edge_land_fraction": float(torn),
        "net_fill_ratio": float(fill),
        "vertices_lonlat": [
            [float(a), float(b)] for a, b in zip(*P.to_lonlat(np.asarray(verts_r)))
        ],
    }

    oct_verts = P.octahedron(bfly_lon)
    b2d, b3d = P.butterfly_net(oct_verts, join_index=bfly_join)

    todo = [
        P.WinkelTripel(lon0=lon0),
        P.EqualEarth(lon0=lon0),
        P.LambertAzimuthalHemispheres.equatorial(lon0=hemi_lon),
        P.LambertAzimuthalHemispheres.polar(lon0=lon0),
        P.PolyhedralNet(b2d, b3d, "waterman_butterfly", "Waterman Butterfly"),
        P.PolyhedralNet(faces2d, faces3d, "dymaxion", "Dymaxion"),
    ]

    from PIL import Image

    written = []
    for proj in todo:
        if args.only and proj.name != args.only:
            continue
        img = render(proj, base, args.width, args.supersample)
        path = OUT / f"vesper_{proj.name}.png"
        Image.fromarray(img).save(path, optimize=True)
        print(f"  {path.name}  {img.shape[1]}x{img.shape[0]}")
        written.append(path.name)

    src_prov = json.loads((BUILD / "basemap_provenance.json").read_text())
    prov.update(
        {
            "basemap": src_prov,
            "outputs": written,
            "width": args.width,
            "supersample": args.supersample,
            "git_commit": subprocess.run(
                ["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True
            ).stdout.strip(),
        }
    )
    (OUT / "projections_provenance.json").write_text(json.dumps(prov, indent=2) + "\n")


if __name__ == "__main__":
    main()
