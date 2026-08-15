"""Inverse map projections, written against the render grid rather than a globe.

Every projection here is expressed as an inverse: given plane coordinates,
return longitude, latitude and a validity mask. That is the direction a
raster renderer needs, and for the polyhedral projections it is also the
easy direction, because a gnomonic face inverts to a barycentric combination
of the face's vertices.

Conventions match the World Orogen export: a unit vector is
(cos lat sin lon, sin lat, cos lat cos lon), so y is north.
"""

from __future__ import annotations

import numpy as np

PHI1 = np.arccos(2.0 / np.pi)  # Winkel tripel standard parallel
EE = (1.340264, -0.081106, 0.000893, 0.003796)  # Equal Earth coefficients
EE_M = np.sqrt(3.0) / 2.0


def unit(lat_deg, lon_deg):
    la, lo = np.deg2rad(lat_deg), np.deg2rad(lon_deg)
    return np.array([np.cos(la) * np.sin(lo), np.sin(la), np.cos(la) * np.cos(lo)])


def to_lonlat(v):
    """Unit vectors (..., 3) to degrees."""
    lat = np.rad2deg(np.arcsin(np.clip(v[..., 1], -1, 1)))
    lon = np.rad2deg(np.arctan2(v[..., 0], v[..., 2]))
    return lon, lat


def wrap180(deg):
    return (deg + 180.0) % 360.0 - 180.0


# --------------------------------------------------------------------------
# Winkel tripel


def winkel_forward(lam, phi):
    half = lam / 2.0
    cosalpha = np.clip(np.cos(phi) * np.cos(half), -1, 1)
    alpha = np.arccos(cosalpha)
    # sinc(alpha) = sin(alpha)/alpha, continuous at 0
    sinc = np.where(np.abs(alpha) < 1e-9, 1.0, np.sin(alpha) / np.where(alpha == 0, 1, alpha))
    x = 0.5 * (lam * np.cos(PHI1) + 2.0 * np.cos(phi) * np.sin(half) / sinc)
    y = 0.5 * (phi + np.sin(phi) / sinc)
    return x, y


class WinkelTripel:
    """Neither equal-area nor conformal, and the usual compromise world map.

    Inverted by damped Newton on the two-variable forward, which has no closed
    form. Points off the oval simply fail to converge and are dropped.
    """

    name = "winkel_tripel"
    title = "Winkel Tripel"

    def __init__(self, lon0=0.0):
        self.lon0 = lon0
        self.extent = (-(2.0 + np.pi) / 2.0, (2.0 + np.pi) / 2.0, -np.pi / 2, np.pi / 2)

    def inverse(self, x, y):
        lam = np.clip(x / np.cos(PHI1), -np.pi, np.pi)
        phi = np.clip(y, -np.pi / 2, np.pi / 2)
        h = 1e-6
        for _ in range(40):
            fx, fy = winkel_forward(lam, phi)
            rx, ry = fx - x, fy - y
            ax, ay = winkel_forward(lam + h, phi)
            bx, by = winkel_forward(lam, phi + h)
            j11, j21 = (ax - fx) / h, (ay - fy) / h
            j12, j22 = (bx - fx) / h, (by - fy) / h
            det = j11 * j22 - j12 * j21
            det = np.where(np.abs(det) < 1e-12, 1e-12, det)
            dlam = (rx * j22 - ry * j12) / det
            dphi = (ry * j11 - rx * j21) / det
            lam = np.clip(lam - 0.9 * dlam, -np.pi, np.pi)
            phi = np.clip(phi - 0.9 * dphi, -np.pi / 2, np.pi / 2)
        fx, fy = winkel_forward(lam, phi)
        valid = np.hypot(fx - x, fy - y) < 1e-4
        return wrap180(np.rad2deg(lam) + self.lon0), np.rad2deg(phi), valid


# --------------------------------------------------------------------------
# Equal Earth


def _ee_y(theta):
    a1, a2, a3, a4 = EE
    t2 = theta * theta
    return theta * (a1 + t2 * (a2 + t2 * t2 * (a3 + a4 * t2)))


def _ee_dy(theta):
    a1, a2, a3, a4 = EE
    t2 = theta * theta
    return a1 + t2 * (3 * a2 + t2 * t2 * (7 * a3 + 9 * a4 * t2))


class EqualEarth:
    """Savric, Patterson and Jenny (2018). Equal-area, and the pseudocylinder
    that flatters high latitudes least."""

    name = "equal_earth"
    title = "Equal Earth"

    def __init__(self, lon0=0.0):
        self.lon0 = lon0
        ymax = _ee_y(np.arcsin(EE_M))
        xmax = np.pi * np.cos(np.arcsin(EE_M) * 0) / (EE_M * _ee_dy(0.0))
        self.extent = (-xmax, xmax, -ymax, ymax)

    def inverse(self, x, y):
        theta = np.clip(y / EE[0], -np.arcsin(EE_M), np.arcsin(EE_M))
        for _ in range(20):
            theta = theta - (_ee_y(theta) - y) / _ee_dy(theta)
        ok = np.abs(theta) <= np.arcsin(EE_M) + 1e-9
        theta = np.clip(theta, -np.arcsin(EE_M), np.arcsin(EE_M))
        phi = np.arcsin(np.clip(np.sin(theta) / EE_M, -1, 1))
        lam = x * EE_M * _ee_dy(theta) / np.maximum(np.cos(theta), 1e-12)
        valid = ok & (np.abs(lam) <= np.pi + 1e-9) & (np.abs(_ee_y(theta) - y) < 1e-6)
        return wrap180(np.rad2deg(np.clip(lam, -np.pi, np.pi)) + self.lon0), np.rad2deg(phi), valid


# --------------------------------------------------------------------------
# Lambert azimuthal equal-area, as a two-hemisphere plate


class LambertAzimuthalHemispheres:
    """Two equal-area hemispheres side by side.

    A single azimuthal disc can hold the whole sphere, but the outer half of
    it is unreadable, so the classic double-hemisphere plate is the honest
    layout: each disc is one hemisphere. The two centres are opposite points,
    and which pair to use is the only real choice. An equatorial pair cuts on a
    meridian, which can be steered onto water; a polar pair cuts on the
    equator, which cannot be steered at all but buys an undistorted view of
    both caps.
    """

    def __init__(
        self,
        centres=((0.0, 0.0), (0.0, 180.0)),
        gap=0.12,
        name="lambert_azimuthal",
        title="Lambert Azimuthal Equal-Area",
    ):
        (self.lat_a, self.lon_a), (self.lat_b, self.lon_b) = centres
        self.name, self.title = name, title
        self.r = np.sqrt(2.0)  # hemisphere radius for a unit sphere
        self.gap = gap * self.r
        self.extent = (-(2 * self.r + self.gap), (2 * self.r + self.gap), -self.r, self.r)

    @classmethod
    def equatorial(cls, lon0=0.0, **kw):
        return cls(centres=((0.0, lon0), (0.0, lon0 + 180.0)), **kw)

    @classmethod
    def polar(cls, lon0=0.0, **kw):
        """North on the left, south on the right.

        The north disc is drawn with `lon0` towards the bottom and the south
        disc with `lon0 + 180` towards the top, which is what makes the two
        rims carry the same longitude where they meet, so the equator reads
        continuously across the join.
        """
        kw.setdefault("name", "lambert_azimuthal_polar")
        kw.setdefault("title", "Lambert Azimuthal Equal-Area, polar aspect")
        return cls(centres=((90.0, lon0), (-90.0, lon0 + 180.0)), **kw)

    def _disc(self, x, y, lon0, lat0):
        rho = np.hypot(x, y)
        inside = rho <= self.r
        rho_s = np.where(rho < 1e-12, 1e-12, rho)
        c = 2.0 * np.arcsin(np.clip(rho_s / 2.0, -1, 1))
        la0 = np.deg2rad(lat0)
        sin_lat = np.cos(c) * np.sin(la0) + y * np.sin(c) * np.cos(la0) / rho_s
        lat = np.rad2deg(np.arcsin(np.clip(sin_lat, -1, 1)))
        lon = lon0 + np.rad2deg(
            np.arctan2(x * np.sin(c), rho_s * np.cos(c) * np.cos(la0) - y * np.sin(c) * np.sin(la0))
        )
        return wrap180(lon), lat, inside

    def inverse(self, x, y):
        left_cx = -(self.r + self.gap / 2.0)
        right_cx = self.r + self.gap / 2.0
        lo_a, la_a, in_a = self._disc(x - left_cx, y, self.lon_a, self.lat_a)
        lo_b, la_b, in_b = self._disc(x - right_cx, y, self.lon_b, self.lat_b)
        lon = np.where(in_a, lo_a, lo_b)
        lat = np.where(in_a, la_a, la_b)
        return lon, lat, in_a | in_b


# --------------------------------------------------------------------------
# Polyhedral projections: a net of gnomonic triangular faces


class PolyhedralNet:
    """A net of flat triangles, each the gnomonic image of a spherical triangle.

    Gnomonic projection sends great circles to straight lines, so the spherical
    triangle spanned by three polyhedron vertices maps exactly onto the planar
    triangle they span. Inverting is therefore just barycentric coordinates:
    the same weights that place a point inside the planar triangle place it,
    after normalising, inside the spherical one.
    """

    def __init__(self, faces2d, faces3d, name, title):
        self.faces2d = [np.asarray(f, float) for f in faces2d]
        self.faces3d = [np.asarray(f, float) for f in faces3d]
        self.name, self.title = name, title
        pts = np.concatenate(self.faces2d)
        pad = 0.02 * (pts[:, 0].max() - pts[:, 0].min())
        self.extent = (
            pts[:, 0].min() - pad,
            pts[:, 0].max() + pad,
            pts[:, 1].min() - pad,
            pts[:, 1].max() + pad,
        )

    def inverse(self, x, y):
        lon = np.zeros(x.shape)
        lat = np.zeros(x.shape)
        valid = np.zeros(x.shape, bool)
        eps = 1e-9
        for tri2, tri3 in zip(self.faces2d, self.faces3d):
            a, b, c = tri2
            det = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1])
            px, py = x - c[0], y - c[1]
            w0 = ((b[1] - c[1]) * px + (c[0] - b[0]) * py) / det
            w1 = ((c[1] - a[1]) * px + (a[0] - c[0]) * py) / det
            w2 = 1.0 - w0 - w1
            hit = (w0 >= -eps) & (w1 >= -eps) & (w2 >= -eps) & ~valid
            if not hit.any():
                continue
            v = (
                w0[hit, None] * tri3[0][None, :]
                + w1[hit, None] * tri3[1][None, :]
                + w2[hit, None] * tri3[2][None, :]
            )
            v /= np.linalg.norm(v, axis=1, keepdims=True)
            lo, la = to_lonlat(v)
            lon[hit], lat[hit] = lo, la
            valid |= hit
        return lon, lat, valid


def _place_third(p, q, other):
    """Third vertex of an equilateral triangle on edge p-q, away from `other`."""
    d = q - p
    s = np.hypot(*d)
    n = np.array([-d[1], d[0]]) / s
    m = (p + q) / 2.0
    cand = m + n * (np.sqrt(3.0) / 2.0) * s
    if other is not None and np.dot(cand - m, other - m) > 0:
        cand = m - n * (np.sqrt(3.0) / 2.0) * s
    return cand


def unfold(vertices, faces, tree_parent, root):
    """Unfold a triangulated polyhedron into the plane along a spanning tree.

    `tree_parent[f]` is the face `f` was unfolded from, or -1 for the root.
    Returns the 2D triangle of each face, vertex order matching `faces`.
    """
    v3 = np.asarray(vertices, float)
    side = float(np.linalg.norm(v3[faces[0][0]] - v3[faces[0][1]]))
    placed = {}  # face -> {vertex id: 2D point}

    a, b, c = faces[root]
    placed[root] = {
        a: np.array([0.0, 0.0]),
        b: np.array([side, 0.0]),
        c: np.array([side / 2.0, side * np.sqrt(3.0) / 2.0]),
    }

    children = {}
    for f, p in enumerate(tree_parent):
        children.setdefault(p, []).append(f)

    order = [root]
    while order:
        parent = order.pop()
        for f in children.get(parent, []):
            if f == root:
                continue
            shared = [v for v in faces[f] if v in placed[parent]]
            assert len(shared) == 2, "child face must share an edge with its parent"
            (p_id, q_id) = shared
            odd = next(v for v in faces[parent] if v not in shared)
            p, q = placed[parent][p_id], placed[parent][q_id]
            third = next(v for v in faces[f] if v not in shared)
            placed[f] = {
                p_id: p,
                q_id: q,
                third: _place_third(p, q, placed[parent][odd]),
            }
            order.append(f)

    faces2d = [np.array([placed[f][v] for v in faces[f]]) for f in range(len(faces))]
    faces3d = [np.array([v3[v] for v in faces[f]]) for f in range(len(faces))]
    return faces2d, faces3d


def triangles_overlap(t1, t2, shrink=0.94):
    """Separating-axis test on triangles pulled in from their shared edges."""

    def shrunk(t):
        c = t.mean(axis=0)
        return c + (t - c) * shrink

    a, b = shrunk(t1), shrunk(t2)
    for tri in (a, b):
        for i in range(3):
            edge = tri[(i + 1) % 3] - tri[i]
            axis = np.array([-edge[1], edge[0]])
            pa, pb = a @ axis, b @ axis
            if pa.max() < pb.min() or pb.max() < pa.min():
                return False
    return True


def net_overlaps(faces2d):
    for i in range(len(faces2d)):
        for j in range(i + 1, len(faces2d)):
            if triangles_overlap(faces2d[i], faces2d[j]):
                return True
    return False


# --------------------------------------------------------------------------
# Polyhedron geometry


def icosahedron():
    """Unit icosahedron with a vertex at each pole, and its 20 faces."""
    verts = [unit(90.0, 0.0)]
    lat = np.rad2deg(np.arctan(0.5))
    for k in range(5):
        verts.append(unit(lat, -180.0 + 72.0 * k))
    for k in range(5):
        verts.append(unit(-lat, -144.0 + 72.0 * k))
    verts.append(unit(-90.0, 0.0))
    v = np.array(verts)

    faces = []
    for k in range(5):
        faces.append((0, 1 + k, 1 + (k + 1) % 5))  # north cap
    for k in range(5):
        faces.append((1 + k, 6 + k, 1 + (k + 1) % 5))  # upper band
    for k in range(5):
        faces.append((1 + (k + 1) % 5, 6 + k, 6 + (k + 1) % 5))  # lower band
    for k in range(5):
        faces.append((11, 6 + (k + 1) % 5, 6 + k))  # south cap
    return v, faces


def octahedron(lon0=0.0):
    verts = [unit(90.0, 0.0)] + [unit(0.0, lon0 + 90.0 * k) for k in range(4)] + [unit(-90.0, 0.0)]
    return np.array(verts)


def rotate(v, quat):
    """Rotate unit vectors by a quaternion (w, x, y, z)."""
    w, x, y, z = quat
    r = np.array(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
            [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
            [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
        ]
    )
    return np.asarray(v) @ r.T


def butterfly_net(verts, join_index=1):
    """Cahill-style butterfly: two fans of four faces around opposite equatorial
    vertices, joined at a third and cut at the fourth.

    Four equilateral faces meet at each equatorial vertex and cover 240 degrees,
    so a fan unfolds flat with a 120 degree notch. Two of them, notches
    outward, is the butterfly.
    """
    n, s = 0, 5
    eq = [1, 2, 3, 4]
    join = eq[join_index % 4]
    left_c = eq[(join_index - 1) % 4]
    right_c = eq[(join_index + 1) % 4]
    cut = eq[(join_index + 2) % 4]
    side = float(np.linalg.norm(verts[n] - verts[eq[0]]))

    def fan(centre, out_sign):
        """Fan around `centre`; the join vertex points towards x = 0."""
        cx = -out_sign * side
        c2 = np.array([cx, 0.0])

        def at(angle_deg):
            a = np.deg2rad(angle_deg)
            return c2 + side * np.array([np.cos(a), np.sin(a)])

        j = 0.0 if out_sign > 0 else 180.0
        pole_n = j + out_sign * 60.0
        cut_n = j + out_sign * 120.0
        pole_s = j - out_sign * 60.0
        cut_s = j - out_sign * 120.0
        f2 = [
            np.array([c2, at(pole_n), at(j)]),
            np.array([c2, at(pole_n), at(cut_n)]),
            np.array([c2, at(pole_s), at(j)]),
            np.array([c2, at(pole_s), at(cut_s)]),
        ]
        f3 = [
            np.array([verts[centre], verts[n], verts[join]]),
            np.array([verts[centre], verts[n], verts[cut]]),
            np.array([verts[centre], verts[s], verts[join]]),
            np.array([verts[centre], verts[s], verts[cut]]),
        ]
        return f2, f3

    l2, l3 = fan(left_c, 1)
    r2, r3 = fan(right_c, -1)
    return l2 + r2, l3 + r3
