"""Render the Vesper base map into a projection, or into several.

The AuthaGraph-style tetrahedral rectangle is the one drawn by default: it is
the projection that shows the whole world at once with the least distortion of
shape, so it is the picture to look at when the question is what the world looks
like. The other six are behind their own flags, and each costs its own
orientation search -- the icosahedral net alone searches 1500 attitudes -- so a
default that drew all seven charged every caller for six pictures it did not
ask for.

Where a projection has a free orientation, it is chosen against the geography
rather than assumed: central meridians, the butterfly's cut meridians and the
icosahedron's attitude are all searched for the arrangement that puts the
interruption over the most water. Fuller picked his net so Earth's continents
came out whole; the same argument applied to this world gives a different net,
so the net is derived here instead of copied.

    python maps/render_projections.py                      # the authagraph
    python maps/render_projections.py --dymaxion --winkel-tripel
    python maps/render_projections.py --all [--width 4000]

Output goes to `maps/data/<build>/<frame>/`, one UUID-named frame per state of
the world, with `maps/data/<build>/INDEX.json` saying what each frame is.
`frames.py` carries why the name is a UUID and what a frame is keyed on.
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
import frames  # noqa: E402
import projections as P  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BUILD = Path(__file__).resolve().parent / "build"

# The base map is 5760 wide by default and the mesh beneath it finer still, so
# this is a render size rather than a limit. Past the base map's own width the
# extra pixels are interpolation and not detail.
DEFAULT_WIDTH = 4000

BACKGROUND = np.array([13, 16, 21.0])
GRATICULE = np.array([255, 255, 255.0])
GRATICULE_HALO = np.array([6, 10, 16.0])
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


def distortion_fields(proj, nx=460, ny=265, h=2e-5):
    """Anisotropy and area scale, sampled over a projection's plane."""
    x0, x1, y0, y1 = proj.extent
    X, Y = np.meshgrid(
        np.linspace(x0 + 3 * h, x1 - 3 * h, nx), np.linspace(y0 + 3 * h, y1 - 3 * h, ny)
    )

    def sphere(a, b):
        lon, lat, _ = proj.inverse(a, b)
        return np.stack(P.unit(lat, lon), axis=-1)

    dx = (sphere(X + h, Y) - sphere(X - h, Y)) / (2 * h)
    dy = (sphere(X, Y + h) - sphere(X, Y - h)) / (2 * h)
    e = np.sum(dx * dx, -1)
    f = np.sum(dx * dy, -1)
    g = np.sum(dy * dy, -1)
    tr = e + g
    det = np.maximum(e * g - f * f, 1e-30)
    disc = np.sqrt(np.maximum(tr * tr / 4 - det, 0))
    aniso = np.sqrt((tr / 2 + disc) / np.maximum(tr / 2 - disc, 1e-30))
    area = np.sqrt(det)
    return X, Y, aniso, area / np.median(area)


def best_tetrahedron_attitude(base, trials=40000, samples=1600, grid=7000,
                              cut_max=0.11, gross_max=0.03, seed=11081972):
    """Attitude that keeps land off the cut and out of the distortion.

    Three things are being placed at once, and the first version of this map
    only placed one of them. The border is the map's only interruption, so land
    on it is torn. But the distortion pattern is also fixed in the rectangle --
    shapes worst at the four side midpoints where the tetrahedron's vertices
    are, and area worst there too -- so turning the tetrahedron decides which
    parts of the world land in it. Searching the border alone left exactly as
    much distortion over land as chance would give.

    Both area terms are weighted by the area a sample stands for on the globe,
    not by the room it takes up on the page: a region drawn at a third of its
    size is a third of the evidence on the page and three times the error.

    Turning the tetrahedron rotates the whole map rigidly, so the sample points
    are mapped to the sphere once and merely rotated per trial.
    """
    proj = P.TetrahedralRectangle()
    x0, x1, y0, y1 = proj.extent
    eps = 1e-4
    n = samples // 4
    t = np.linspace(0, 1, n)
    bx = np.concatenate([x0 + t * (x1 - x0), x0 + t * (x1 - x0), np.full(n, x0 + eps),
                         np.full(n, x1 - eps)])
    by = np.concatenate([np.full(n, y0 + eps), np.full(n, y1 - eps), y0 + t * (y1 - y0),
                         y0 + t * (y1 - y0)])
    lon, lat, _ = proj.inverse(bx, by)
    border = np.stack(P.unit(lat, lon), axis=-1)

    gx, gy, aniso, area = distortion_fields(proj)
    rng = np.random.default_rng(seed)
    pick = rng.choice(gx.size, size=min(grid, gx.size), replace=False)
    lon, lat, _ = proj.inverse(gx.ravel()[pick], gy.ravel()[pick])
    interior = np.stack(P.unit(lat, lon), axis=-1)
    scale = area.ravel()[pick]
    globe = 1.0 / np.maximum(scale, 1e-3)  # world area each sample stands for
    excess = np.clip(aniso.ravel()[pick] - 1.0, 0.0, None)
    # Only gross size errors are counted. A continent drawn a fifth too large
    # reads as correct; one drawn at half size does not.
    misscaled = ((scale < 0.6) | (scale > 1.67)).astype(float)

    def metrics(q):
        w, i, j, k = q
        r = np.array(
            [
                [1 - 2 * (j * j + k * k), 2 * (i * j - w * k), 2 * (i * k + w * j)],
                [2 * (i * j + w * k), 1 - 2 * (i * i + k * k), 2 * (j * k - w * i)],
                [2 * (i * k - w * j), 2 * (j * k + w * i), 1 - 2 * (i * i + j * j)],
            ]
        )
        bl, ba = P.to_lonlat(border @ r.T)
        il, ia = P.to_lonlat(interior @ r.T)
        land = globe * base.land_at(il, ia)
        total = land.sum()
        return (
            float(base.land_at(bl, ba).mean()),
            float(np.dot(land, excess)) / total,
            float(np.dot(land, misscaled)) / total,
        )

    # Random attitudes resolve SO(3) to about six degrees, too coarse to land
    # in the narrow band that satisfies all three at once, so the best few are
    # refined with a shrinking step. Refinement descends a penalised cost,
    # because an attitude usually has to pass through infeasible ground to
    # reach the good region; the winner is then picked by the standards
    # themselves, which a penalty alone would let a bad attitude buy its way
    # past.
    def penalised(m):
        cut, shape, gross = m
        return shape + 20.0 * max(0.0, cut - cut_max) + 6.0 * max(0.0, gross - gross_max)

    quats = random_quaternions(trials, seed=seed)
    order = sorted(range(trials), key=lambda n: penalised(metrics(quats[n])))

    refined = []
    rng = np.random.default_rng(seed)
    for n in order[:32]:
        q = quats[n]
        c = penalised(metrics(q))
        for sigma in (0.06, 0.03, 0.015, 0.007, 0.003):
            for _ in range(100):
                trial = q + rng.normal(0, sigma, 4)
                trial /= np.linalg.norm(trial)
                ct = penalised(metrics(trial))
                if ct < c:
                    q, c = trial, ct
        refined.append((q, metrics(q)))

    for cmax, gmax in [(cut_max, gross_max), (cut_max * 1.3, gross_max * 2), (1.0, 1.0)]:
        ok = [(q, m) for q, m in refined if m[0] <= cmax and m[2] <= gmax]
        if ok:
            break
    q, (cut, shape, gross) = min(ok, key=lambda t: t[1][1])
    return q, {"cut": cut, "land_shape_excess": shape, "land_grossly_misscaled": gross}


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
        # The graticule needs a vertical derivative, so each block is computed
        # one row proud at each end and trimmed afterwards. Without it every
        # block boundary leaves a one-pixel discontinuity along the lines.
        pad0, pad1 = (1 if r0 > 0 else 0), (1 if r1 < H else 0)
        X, Y = np.meshgrid(xs, ys[r0 - pad0 : r1 + pad1])
        lon, lat, valid = proj.inverse(X, Y)
        rgb = base.sample(lon, lat).astype(np.float64)
        rgb[~valid] = BACKGROUND

        if graticule:
            rgb = draw_graticule(rgb, lon, lat, valid, graticule, supersample)
        rgb = draw_outline(rgb, valid)

        if pad0 or pad1:
            rgb = rgb[pad0 : rgb.shape[0] - pad1]
        rows = (r1 - r0) // supersample
        acc[r0 // supersample : r0 // supersample + rows] = (
            rgb.reshape(rows, supersample, width, supersample, 3).mean(axis=(1, 3))
        )
    return np.clip(acc, 0, 255).astype(np.uint8)


def _field_gradient(field, wrap=None):
    """Degrees of the field per pixel, following the shorter way round."""
    gy = np.gradient(field, axis=0)
    gx = np.gradient(field, axis=1)
    if wrap:
        gy = np.remainder(gy + wrap / 2, wrap) - wrap / 2
        gx = np.remainder(gx + wrap / 2, wrap) - wrap / 2
    return np.hypot(gx, gy)


def _line_coverage(field, spacing, grad, half_width):
    """Antialiased coverage of the lines at multiples of `spacing`.

    Distance to the nearest line is in degrees; dividing by the local gradient
    puts it in pixels, which is what makes a line the same weight everywhere
    however hard the projection is stretching. Coverage then ramps across the
    last pixel instead of switching, so the lines do not crawl.
    """
    d = np.abs(np.remainder(field + spacing / 2.0, spacing) - spacing / 2.0)
    return np.clip(half_width + 0.5 - d / np.maximum(grad, 1e-12), 0.0, 1.0)


def draw_graticule(rgb, lon, lat, valid, spacing, supersample):
    """A dark halo under a light line, so the graticule reads over deep ocean
    and bright desert alike."""
    glat = _field_gradient(lat)
    glon = _field_gradient(lon, wrap=360.0)

    # A cut shows up as a jump of many degrees inside one pixel, where a real
    # parallel is a fraction of a degree; anything above the threshold is a
    # seam and drawing on it would be an artefact. The old rule was relative to
    # the median gradient, which also erased lines wherever the projection was
    # legitimately compressed.
    ok = valid & (glat < 20.0) & (glon < 20.0)
    ok &= binary_erosion(valid, iterations=2 * supersample, border_value=1)
    ok = ok.astype(np.float64)

    s = float(supersample)
    for field, spacing_, grad, half, colour, alpha in [
        (lat, spacing, glat, 1.7 * s, GRATICULE_HALO, 0.32),
        (lon, spacing, glon, 1.7 * s, GRATICULE_HALO, 0.32),
        (lat, spacing, glat, 0.65 * s, GRATICULE, 0.60),
        (lon, spacing, glon, 0.65 * s, GRATICULE, 0.60),
        (lat, 180.0, glat, 0.95 * s, GRATICULE, 0.80),
        (lon, 180.0, glon, 0.95 * s, GRATICULE, 0.80),
    ]:
        a = (_line_coverage(field, spacing_, grad, half) * alpha * ok)[..., None]
        rgb = rgb * (1.0 - a) + colour * a
    return rgb


def draw_outline(rgb, valid, thickness=1):
    edge = valid & ~binary_erosion(valid, iterations=thickness, border_value=1)
    rgb[edge] = rgb[edge] * 0.45 + OUTLINE * 0.55
    return rgb


# --------------------------------------------------------------------------
# The projections, and the orientation each one needs
#
# One entry per output. The searches are LAZY and memoised: the tetrahedron
# refines 40000 attitudes and the icosahedral net scores 1500, so a run that
# draws one projection must not pay for the searches of the other six. Each
# search records what it decided in `prov` as it decides it, so the frame's
# provenance holds exactly the orientations that were actually chosen.


class Orientation:
    """The orientation searches, each run at most once and only if asked for."""

    def __init__(self, base, prov, trials):
        self.base, self.prov, self.trials = base, prov, trials
        self._memo = {}

    def _once(self, key, fn):
        if key not in self._memo:
            self._memo[key] = fn()
        return self._memo[key]

    def central_meridian(self):
        def run():
            lon0, cut = best_central_meridian(self.base)
            print(f"central meridian {lon0:+.0f} "
                  f"(antimeridian is {cut * 100:.0f}% land)")
            self.prov["central_meridian"] = {
                "lon0": lon0, "antimeridian_land_fraction": cut}
            return lon0
        return self._once("central_meridian", run)

    def hemisphere_meridian(self):
        def run():
            lon0, cut = best_hemisphere_meridian(self.base)
            print(f"hemisphere centre {lon0:+.0f} "
                  f"(bounding meridians {cut * 100:.0f}% land)")
            self.prov["hemispheres"] = {"lon0": lon0, "cut_land_fraction": cut}
            return lon0
        return self._once("hemisphere_meridian", run)

    def polar_meridian(self):
        """The polar plate has no orientation to search: its cut is the equator,
        and the only free parameter is which meridian points down. It reuses the
        central meridian so the set of maps stays consistent with itself."""
        def run():
            lon0 = self.central_meridian()
            equator_land = parallel_land_fraction(self.base, 0.0)
            print(f"polar plate cuts the equator "
                  f"({equator_land * 100:.0f}% land)")
            self.prov["polar_hemispheres"] = {
                "lon0": lon0, "cut": "equator",
                "cut_land_fraction": equator_land}
            return lon0
        return self._once("polar_meridian", run)

    def butterfly(self):
        def run():
            lon, join, cut = best_butterfly(self.base)
            print(f"butterfly vertices at {lon:+.1f}+90k, join {join} "
                  f"(cut meridians {cut * 100:.0f}% land)")
            self.prov["butterfly"] = {
                "equatorial_vertex_lon": lon, "join_index": join,
                "cut_land_fraction": cut}
            return P.butterfly_net(P.octahedron(lon), join_index=join)
        return self._once("butterfly", run)

    def icosahedron(self):
        def run():
            verts, faces = P.icosahedron()
            best = None
            for edge_land, q in best_icosahedron_attitudes(
                    self.base, trials=self.trials):
                v = P.rotate(verts, q)
                f2, f3, torn, fill = choose_net(self.base, v, faces)
                score = 2.5 * torn + (1.0 - fill)
                if best is None or score < best[0]:
                    best = (score, q, v, f2, f3, torn, fill, edge_land)
            _, quat, verts_r, faces2d, faces3d, torn, fill, edge_land = best
            print(f"icosahedron edges {edge_land * 100:.0f}% land; net tears "
                  f"{torn * 100:.0f}% land per cut edge, fills "
                  f"{fill * 100:.0f}% of its box")
            self.prov["icosahedron"] = {
                "quaternion": [float(q) for q in quat],
                "mean_edge_land_fraction": float(edge_land),
                "mean_cut_edge_land_fraction": float(torn),
                "net_fill_ratio": float(fill),
                "attitude_trials": self.trials,
                "vertices_lonlat": [
                    [float(a), float(b)]
                    for a, b in zip(*P.to_lonlat(np.asarray(verts_r)))
                ],
            }
            return faces2d, faces3d
        return self._once("icosahedron", run)

    def tetrahedron(self):
        def run():
            q, m = best_tetrahedron_attitude(self.base)
            tetra = P.TetrahedralRectangle(verts=P.rotate(P.tetrahedron()[0], q))
            print(f"tetrahedral rectangle border {m['cut'] * 100:.0f}% land; "
                  f"over land, mean anisotropy "
                  f"{1 + m['land_shape_excess']:.2f} and "
                  f"{m['land_grossly_misscaled'] * 100:.0f}% of land area "
                  f"grossly mis-scaled")
            self.prov["tetrahedron"] = {
                "quaternion": [float(v) for v in q],
                "border_land_fraction": m["cut"],
                "land_mean_anisotropy": 1 + m["land_shape_excess"],
                "land_area_grossly_misscaled_fraction": m["land_grossly_misscaled"],
                "subdivision": 16,
                "area_weight": 0.5,
                "vertices_lonlat": [
                    [float(a), float(b)]
                    for a, b in zip(*P.to_lonlat(np.asarray(tetra.verts)))
                ],
            }
            return tetra
        return self._once("tetrahedron", run)


# name -> what it takes to build it. The order is the order they are drawn in.
PROJECTIONS = {
    "authagraph": lambda o: o.tetrahedron(),
    "winkel_tripel": lambda o: P.WinkelTripel(lon0=o.central_meridian()),
    "equal_earth": lambda o: P.EqualEarth(lon0=o.central_meridian()),
    "lambert_azimuthal":
        lambda o: P.LambertAzimuthalHemispheres.equatorial(
            lon0=o.hemisphere_meridian()),
    "lambert_azimuthal_polar":
        lambda o: P.LambertAzimuthalHemispheres.polar(lon0=o.polar_meridian()),
    "waterman_butterfly":
        lambda o: P.PolyhedralNet(*o.butterfly(), "waterman_butterfly",
                                  "Waterman Butterfly"),
    "dymaxion": lambda o: P.PolyhedralNet(*o.icosahedron(), "dymaxion",
                                          "Dymaxion"),
}
DEFAULT_PROJECTION = "authagraph"

# The blocks the searches above write, and the only part of a frame's sidecar
# that survives from one invocation to the next.
ORIENTATIONS = ("central_meridian", "hemispheres", "polar_hemispheres",
                "butterfly", "icosahedron", "tetrahedron")


def selected(args) -> list[str]:
    """Which projections this invocation asks for.

    No flag means the authagraph alone, and a flag REPLACES that default rather
    than adding to it: `--dymaxion` asks for the Dymaxion, and a caller who
    wants both says so. Every projection has a flag, including the default one,
    so there is one way to name each of them.
    """
    if args.all:
        return list(PROJECTIONS)
    chosen = [n for n in PROJECTIONS if getattr(args, n)]
    return chosen or [DEFAULT_PROJECTION]


# --------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(
        description="Render the base map into the AuthaGraph-style projection, "
                    "or into the ones named.")
    ap.add_argument("--width", type=int, default=DEFAULT_WIDTH,
                    help=f"render width in pixels (default {DEFAULT_WIDTH})")
    ap.add_argument("--supersample", type=int, default=2,
                    help="samples per output pixel on each axis")
    ap.add_argument("--all", action="store_true",
                    help="every projection rather than the authagraph alone")
    for name in PROJECTIONS:
        ap.add_argument(f"--{name.replace('_', '-')}", dest=name,
                        action="store_true", help=f"render {name}")
    ap.add_argument("--trials", type=int, default=1500,
                    help="icosahedron attitudes searched, for --dymaxion")
    ap.add_argument("--step", default=None,
                    help="the pipeline step this frame stands after; recorded "
                         "in INDEX.json against the frame")
    ap.add_argument("--force", action="store_true",
                    help="re-render even where the frame already holds the "
                         "projection at this size")
    args = ap.parse_args()
    names = selected(args)
    if args.step:
        frames.check_step(args.step)

    base = BaseMap()
    print(f"base map {base.w}x{base.h}")

    # The frame is keyed on what the RASTER was drawn from, taken from the base
    # map's own provenance rather than re-resolved from the config here. The two
    # differ exactly when the base map is older than the config -- which is the
    # case that matters, and recording the config's answer would file the frame
    # under a world it does not show.
    src_prov = json.loads((BUILD / "basemap_provenance.json").read_text())
    if "inputs" not in src_prov:
        raise SystemExit(
            "maps/build/basemap_provenance.json carries no `inputs` block, so "
            "the frame it belongs to cannot be identified. Rebuild it:\n"
            "    python maps/build_basemap.py")

    # REGISTERED BEFORE ANYTHING IS DRAWN, so the id and its row come into
    # existence together; see frames.py. A render is minutes long and this
    # process is not the only one that can be rendering.
    entry = frames.register(src_prov["inputs"])
    entry["climatology_caveat"] = src_prov.get("climatology_caveat")
    entry["basemap_resolution"] = src_prov.get("resolution")
    out_dir = frames.frame_dir(entry)
    out_dir.mkdir(parents=True, exist_ok=True)

    todo = names if args.force else frames.needs_render(
        entry, names, args.width, args.supersample)
    for skipped in [n for n in names if n not in todo]:
        print(f"  {skipped}: already in {entry['frame_id']} at this size")

    # Orientations ACCUMULATE across invocations, because a frame can be
    # extended with a projection whose search has not been run yet, and the one
    # already chosen for a projection beside it still describes the file on
    # disk. Everything else in the sidecar is rewritten below, so only the
    # orientation blocks are carried forward.
    prov = {}
    if (out_dir / "provenance.json").is_file():
        prov = {k: v for k, v in
                json.loads((out_dir / "provenance.json").read_text()).items()
                if k in ORIENTATIONS}
    orientation = Orientation(base, prov, args.trials)

    from PIL import Image

    for name in todo:
        proj = PROJECTIONS[name](orientation)
        img = render(proj, base, args.width, args.supersample)
        path = out_dir / f"vesper_{name}.png"
        Image.fromarray(img).save(path, optimize=True)
        print(f"  {path.relative_to(ROOT)}  {img.shape[1]}x{img.shape[0]}")
        entry.setdefault("outputs", {})[path.name] = {
            "projection": proj.title,
            "width": args.width,
            "supersample": args.supersample,
            "pixels": [int(img.shape[1]), int(img.shape[0])],
        }

    frames.note_step(entry, args.step)
    entry["git_commit"] = frames.git_commit()
    # The sizes are per output rather than per invocation: a frame can be
    # extended with another projection later, at another width, and one pair of
    # top-level width/supersample keys would then describe the last run rather
    # than the files beside it.
    prov.update({
        "frame_id": entry["frame_id"],
        "basemap": src_prov,
        "outputs": entry.get("outputs", {}),
        "git_commit": entry["git_commit"],
    })
    (out_dir / "provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    frames.record(entry)
    print(f"frame {entry['frame_id']}  ->  {out_dir.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
