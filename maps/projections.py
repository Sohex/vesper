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


def tetrahedron():
    """Unit regular tetrahedron, and its four faces wound anticlockwise seen
    from outside."""
    v = np.array(
        [
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ]
    ) / np.sqrt(3.0)
    faces = []
    for f in [(0, 1, 2), (0, 2, 3), (0, 3, 1), (1, 3, 2)]:
        a, b, c = (v[i] for i in f)
        if np.dot(np.cross(b - a, c - a), a + b + c) < 0:
            f = (f[0], f[2], f[1])
        faces.append(f)
    return v, faces


# --------------------------------------------------------------------------
# AuthaGraph-style tetrahedral rectangle


class SnyderTriangle:
    """Snyder's equal-area mapping between a planar and a spherical triangle.

    Each triangle is cut into six sub-sectors by the rays from its centre to
    its vertices and edge midpoints. Within a sub-sector, azimuth is remapped
    so that equal fractions of the sector's area lie either side of it, and
    radius is then remapped so that equal fractions lie inside it. The result
    is exactly equal-area, and because every face performs the identical
    construction against the shared edge, it is continuous across face
    boundaries as well.

    Working inwards from the plane, the azimuth step is the only one without a
    closed form: the swept spherical area is

        A(psi) = psi - arcsin(cos(rho) sin(psi))

    for a face whose spherical inradius is rho, which Newton solves in a few
    iterations.
    """

    def __init__(self, verts3d, side=1.0):
        self.side = side
        self.r_in = side / (2.0 * np.sqrt(3.0))  # planar inradius
        a, b, c = verts3d[0], verts3d[1], verts3d[2]
        n = a + b + c
        n /= np.linalg.norm(n)
        pole = np.cross(a, b)
        pole /= np.linalg.norm(pole)
        self.rho = np.arcsin(np.clip(abs(np.dot(n, pole)), 0, 1))  # spherical inradius
        self.sector = np.pi / 3.0
        self.a_total = self._area(self.sector)

    def _area(self, psi):
        return psi - np.arcsin(np.clip(np.cos(self.rho) * np.sin(psi), -1, 1))

    def _psi(self, phi):
        """Azimuth on the sphere carrying the same area fraction as `phi`."""
        target = self.a_total * np.tan(phi) / np.tan(self.sector)
        psi = np.clip(phi, 0.0, self.sector)
        cr = np.cos(self.rho)
        for _ in range(24):
            f = psi - np.arcsin(np.clip(cr * np.sin(psi), -1, 1)) - target
            d = 1.0 - cr * np.cos(psi) / np.sqrt(np.maximum(1.0 - (cr * np.sin(psi)) ** 2, 1e-15))
            psi = np.clip(psi - f / np.maximum(d, 1e-9), 0.0, self.sector)
        return psi

    def to_sphere(self, r, phi):
        """Planar (radius, azimuth from an edge midpoint) to (z, psi)."""
        psi = self._psi(phi)
        r_max = self.r_in / np.cos(phi)
        t = np.tan(self.rho)
        cos_zmax = np.cos(psi) / np.sqrt(np.cos(psi) ** 2 + t * t)
        cos_z = 1.0 - (1.0 - cos_zmax) * np.clip(r / r_max, 0, 1) ** 2
        return np.arccos(np.clip(cos_z, -1, 1)), psi


CANONICAL2D = np.array([[0.0, 1.0 / np.sqrt(3.0)],
                        [-0.5, -0.5 / np.sqrt(3.0)],
                        [0.5, -0.5 / np.sqrt(3.0)]])


def snyder_face(bary, verts3d, snyder):
    """Snyder's map for one face, from barycentric coordinates.

    Barycentric input keeps the correspondence between the planar corners and
    the sphere's, so the face can be laid out canonically and no orientation
    bookkeeping is needed.
    """
    bary = np.asarray(bary, float)
    p = bary @ CANONICAL2D
    r = np.linalg.norm(p, axis=-1)
    ang = np.mod(np.arctan2(p[..., 1], p[..., 0]) - np.pi / 2.0, 2 * np.pi)

    sector = snyder.sector
    k = np.floor(ang / sector).astype(int)
    u = ang - k * sector
    even = (k % 2) == 0
    phi = np.where(even, sector - u, u)
    z, psi = snyder.to_sphere(r, phi)
    psi_total = np.where(even, (k + 1) * sector - psi, k * sector + psi)

    a3, b3, c3 = verts3d[..., 0, :], verts3d[..., 1, :], verts3d[..., 2, :]
    n = a3 + b3 + c3
    n = n / np.linalg.norm(n, axis=-1, keepdims=True)
    e1 = a3 - n * np.sum(a3 * n, axis=-1, keepdims=True)
    e1 = e1 / np.linalg.norm(e1, axis=-1, keepdims=True)
    e2 = np.cross(n, e1)
    d = np.cos(psi_total)[..., None] * e1 + np.sin(psi_total)[..., None] * e2
    v = np.cos(z)[..., None] * n + np.sin(z)[..., None] * d
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


class TetrahedralRectangle:
    """AuthaGraph-style world map: a tetrahedron, unfolded into a rectangle.

    A flat tetrahedron is the plane folded by half turns about the points of a
    triangular lattice, so its development tiles the plane and any fundamental
    domain of that symmetry is a complete world map. One such domain is a
    1 by sqrt(3) rectangle whose four side midpoints are the four vertices of
    the tetrahedron: the pillowcase. It is continuous everywhere inside, cut
    only along its border, and rectangular, which is the point of the exercise.

    Faces are mapped either by `SnyderTriangle`, which is exactly equal-area
    and shears badly because one face carries a quarter of the sphere, or by a
    `SubdividedFace`, which is the idea behind AuthaGraph's 96 triangles:
    cut the face up, map each piece nearly rigidly, and give up a little area
    accuracy for much better shapes.

    This is not Narukawa's projection, whose subdivision has never been
    published; it is the same idea built from parts that can be checked.
    """

    name = "authagraph"
    title = "AuthaGraph-style Tetrahedral Rectangle"

    def __init__(self, verts=None, faces=None, margin=2, face_map=None,
                 subdivision=16, area_weight=0.5):
        v, f = tetrahedron()
        self.faces = f if faces is None else faces
        self.h = np.sqrt(3.0) / 2.0
        self.extent = (0.0, np.sqrt(3.0), 0.0, 1.0)
        self.snyder = SnyderTriangle(v[list(self.faces[0])])
        if face_map is None and subdivision:
            face_map = optimised_face(subdivision, area_weight)
        self.face_map = face_map
        self._build_table(margin)
        self.verts = v if verts is None else np.asarray(verts, float)

    # -- the development ---------------------------------------------------

    def _corners(self, i, j, up):
        """Lattice cell corners, in the canonical order used by the table."""

        def p(a, b):
            return np.array([a + b * 0.5, b * self.h])

        if up:
            return [p(i, j), p(i + 1, j), p(i, j + 1)]
        return [p(i + 1, j), p(i, j + 1), p(i + 1, j + 1)]

    def _neighbours(self, i, j, up):
        """(cell, indices of the two corners it shares with this one)."""
        if up:
            return [((i, j, 0), (1, 2)), ((i, j - 1, 0), (0, 1)), ((i - 1, j, 0), (0, 2))]
        return [((i, j, 1), (0, 1)), ((i, j + 1, 1), (1, 2)), ((i + 1, j, 1), (0, 2))]

    def _build_table(self, margin):
        """Assign a tetrahedron face to every lattice cell by developing it.

        Walking across a shared edge keeps the two shared vertices and gives
        the new corner the one vertex the current face does not use. That rule
        closes consistently: three faces meet at a tetrahedron vertex and each
        contributes 60 degrees, so six cells meet round a lattice point and the
        assignment comes back to itself.
        """
        i0, i1 = -margin - 2, margin + 2
        j0, j1 = -margin, margin + 3
        seed = (0, 0, 1)
        assign = {seed: tuple(self.faces[0])}
        stack = [seed]
        while stack:
            cell = stack.pop()
            tri = assign[cell]
            i, j, up = cell
            missing = next(v for v in range(4) if v not in tri)
            corners = self._corners(i, j, up)
            for nb, (s0, s1) in self._neighbours(i, j, up):
                if not (i0 <= nb[0] <= i1 and j0 <= nb[1] <= j1):
                    continue
                nb_corners = self._corners(*nb)
                shared = [corners[s0], corners[s1]]
                mapping = []
                for pt in nb_corners:
                    hit = [k for k, s in enumerate(shared) if np.allclose(pt, s, atol=1e-9)]
                    mapping.append(tri[(s0, s1)[hit[0]]] if hit else missing)
                mapping = tuple(mapping)
                if nb in assign:
                    assert assign[nb] == mapping, "development is inconsistent"
                    continue
                assign[nb] = mapping
                stack.append(nb)

        self.i0, self.j0 = i0, j0
        shape = (i1 - i0 + 1, j1 - j0 + 1, 2)
        self.tri_id = np.zeros(shape, np.int16)
        self.swap = np.zeros(shape, bool)
        base, _ = tetrahedron()
        order = {}
        for (i, j, up), tri in assign.items():
            # Half the cells develop face-down, so their triple runs clockwise
            # seen from outside. Put every triple the same way round and record
            # the swap, which the barycentric coordinates then follow.
            a, b, c = (base[t] for t in tri)
            flip = np.dot(np.cross(b - a, c - a), a + b + c) < 0
            if flip:
                tri = (tri[0], tri[2], tri[1])
            k = (i - i0, j - j0, up)
            self.tri_id[k] = order.setdefault(tri, len(order))
            self.swap[k] = flip
        self.ordered_triples = [t for t, _ in sorted(order.items(), key=lambda kv: kv[1])]

    # -- vertex positions --------------------------------------------------

    @property
    def verts(self):
        return self._verts

    @verts.setter
    def verts(self, value):
        self._verts = np.asarray(value, float)
        self._corners3 = np.array(
            [[self._verts[a] for a in tri] for tri in self.ordered_triples]
        )  # (num triples, 3, 3)
        if self.face_map is not None:
            w = self.face_map.weights  # (V, 3)
            v = np.einsum("vw,twc->tvc", w, self._corners3)
            self._face_vertex = v / np.linalg.norm(v, axis=-1, keepdims=True)

    # -- the projection ----------------------------------------------------

    def _cell(self, x, y):
        """Lattice cell and barycentric position within it."""
        lx = y - 0.5
        ly = np.sqrt(3.0) - x
        b = ly / self.h
        a = lx - ly / np.sqrt(3.0)
        i = np.floor(a).astype(int)
        j = np.floor(b).astype(int)
        fx, fy = a - i, b - j
        up = (fx + fy) < 1.0
        bary = np.where(
            up[..., None],
            np.stack([1.0 - fx - fy, fx, fy], axis=-1),
            np.stack([1.0 - fy, 1.0 - fx, fx + fy - 1.0], axis=-1),
        )
        ii = np.clip(i - self.i0, 0, self.tri_id.shape[0] - 1)
        jj = np.clip(j - self.j0, 0, self.tri_id.shape[1] - 1)
        kk = up.astype(int)
        swap = self.swap[ii, jj, kk]
        bary = np.where(swap[..., None], bary[..., [0, 2, 1]], bary)
        return ii, jj, kk, np.clip(bary, 0.0, 1.0)

    def inverse(self, x, y):
        ii, jj, kk, bary = self._cell(x, y)
        if self.face_map is None:
            corners3 = self._corners3[self.tri_id[ii, jj, kk]]
            v = snyder_face(bary, corners3, self.snyder)
        else:
            v = self.face_map.to_sphere(bary, self._face_vertex[self.tri_id[ii, jj, kk]])
        lon, lat = to_lonlat(v)
        x0, x1, y0, y1 = self.extent
        valid = (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)
        return lon, lat, valid


# --------------------------------------------------------------------------
# Subdivided faces: the idea behind AuthaGraph's 96 triangles


def subdivision_vertices(n):
    """Barycentric integer coordinates of a regular triangular subdivision."""
    verts = [(i, j, n - i - j) for i in range(n + 1) for j in range(n + 1 - i)]
    return verts, {v: k for k, v in enumerate(verts)}


def subdivision_triangles(n, index):
    tris = []
    for i in range(n):
        for j in range(n - i):
            k = n - 1 - i - j
            tris.append((index[(i + 1, j, k)], index[(i, j + 1, k)], index[(i, j, k + 1)]))
            if k >= 1:
                tris.append(
                    (
                        index[(i + 1, j + 1, k - 1)],
                        index[(i + 1, j, k)],
                        index[(i, j + 1, k)],
                    )
                )
    return tris


def _classes(verts):
    """Group subdivision vertices into orbits of the triangle's symmetry.

    A face has the symmetry of its three corners, so only the multiset
    {i, j, k} is free. Parameterising by orbit keeps every face identical and
    every shared edge palindromic, which is what makes the four faces agree
    where they meet.
    """
    perms = [(0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)]
    reps, of, how = [], {}, {}
    for v in verts:
        rep = tuple(sorted(v, reverse=True))
        if rep not in of:
            of[rep] = len(reps)
            reps.append(rep)
        # a permutation taking the representative to this vertex
        how[v] = (of[rep], next(p for p in perms if tuple(rep[p[m]] for m in range(3)) == v))
    stab = [[p for p in perms if tuple(r[p[m]] for m in range(3)) == r] for r in reps]
    return reps, how, stab


class SubdividedFace:
    """Piecewise-affine map from a planar face to a spherical one.

    A single smooth mapping of a whole tetrahedron face has to absorb the
    curvature of a quarter of the sphere, and Snyder's does it by shearing
    everything. Cutting the face into small triangles and mapping each one
    affinely spreads the work: each piece is nearly flat, so it can be nearly
    similar to its image, and the curvature is taken up as small kinks along
    the seams instead of shear across the whole face.

    The vertex positions are the free parameters. Each is a set of weights on
    the three face corners, normalised onto the sphere, so a vertex on an edge
    stays on that edge and the corners stay pinned. They are chosen by
    minimising, over the sub-triangles,

        sum (s1/s2 + s2/s1 - 2)  +  w * sum (log(area / target))^2

    the first term being shape and the second area, which is the trade
    AuthaGraph makes when it calls itself approximately equal-area.
    """

    def __init__(self, n, weights):
        self.n = n
        self.verts, self.index = subdivision_vertices(n)
        self.weights = np.asarray(weights, float)  # (num vertices, 3)

    # -- placing the vertices ---------------------------------------------

    @staticmethod
    def _expand(theta, reps, how, stab, verts):
        """Free parameters to per-vertex corner weights."""
        w = np.exp(np.clip(theta, -30, 30)).reshape(len(reps), 3)
        w = np.where(np.array(reps) > 0, w, 0.0)  # a zero stays zero: edges hold
        sym = np.zeros_like(w)
        for r, group in enumerate(stab):
            for p in group:
                sym[r] += w[r, list(p)]
            sym[r] /= len(group)
        sym /= sym.sum(axis=1, keepdims=True)
        out = np.empty((len(verts), 3))
        for m, v in enumerate(verts):
            r, p = how[v]
            out[m] = sym[r, list(p)]
        return out

    @classmethod
    def optimise(cls, n, face3d, planar_corners=None, area_weight=6.0, initial=None):
        from scipy.optimize import minimize

        if planar_corners is None:
            planar_corners = CANONICAL2D
        verts, index = subdivision_vertices(n)
        tris = np.array(subdivision_triangles(n, index))
        reps, how, stab = _classes(verts)

        bary = np.array(verts, float) / n
        planar = bary @ planar_corners  # (V, 2), the regular mesh, held fixed
        corners3 = np.asarray(face3d)  # (3, 3)
        target = spherical_triangle_area(*corners3) / (n * n)

        # Alternate sub-triangles are wound the other way in the plane. Put them
        # all the same way round, so a negative Jacobian means a genuine fold.
        p1, p2, p3 = (planar[tris[:, m]] for m in range(3))
        flip = ((p2 - p1)[:, 0] * (p3 - p1)[:, 1] - (p2 - p1)[:, 1] * (p3 - p1)[:, 0]) < 0
        tris[flip] = tris[flip][:, [0, 2, 1]]

        p1, p2, p3 = (planar[tris[:, m]] for m in range(3))
        r2, r3 = p2 - p1, p3 - p1
        det_r = r2[:, 0] * r3[:, 1] - r2[:, 1] * r3[:, 0]

        def unpack(theta):
            w = cls._expand(theta, reps, how, stab, verts)
            v = w @ corners3
            return v / np.linalg.norm(v, axis=1, keepdims=True)

        def energy(theta):
            v = unpack(theta)
            v1, v2, v3 = (v[tris[:, m]] for m in range(3))
            e1, e2 = v2 - v1, v3 - v1
            nrm = np.cross(e1, e2)
            len_n = np.linalg.norm(nrm, axis=1)
            if np.any(len_n < 1e-12):
                return 1e6
            len_e1 = np.linalg.norm(e1, axis=1)
            u1 = e1 / len_e1[:, None]
            u2 = np.cross(nrm / len_n[:, None], u1)
            q2 = np.stack([len_e1, np.zeros_like(len_e1)], axis=1)
            q3 = np.stack([np.sum(e2 * u1, axis=1), np.sum(e2 * u2, axis=1)], axis=1)
            det_q = q2[:, 0] * q3[:, 1] - q2[:, 1] * q3[:, 0]
            det_j = det_q / det_r
            # A folded triangle is not a map, but a hard cliff here traps the
            # optimiser: once it steps across, everything nearby is equally
            # infinite and there is no gradient home. Penalise smoothly.
            fold = np.clip(1e-3 - det_j, 0.0, None)
            det_j = np.maximum(det_j, 1e-3)
            inv = 1.0 / det_r
            j11 = (q2[:, 0] * r3[:, 1] - q3[:, 0] * r2[:, 1]) * inv
            j12 = (-q2[:, 0] * r3[:, 0] + q3[:, 0] * r2[:, 0]) * inv
            j21 = (q2[:, 1] * r3[:, 1] - q3[:, 1] * r2[:, 1]) * inv
            j22 = (-q2[:, 1] * r3[:, 0] + q3[:, 1] * r2[:, 0]) * inv
            frob = j11 * j11 + j12 * j12 + j21 * j21 + j22 * j22
            shape = frob / det_j - 2.0
            area = spherical_triangle_area(v1, v2, v3)
            bad = np.clip(1e-4 - area, 0.0, None)
            area = np.maximum(area, 1e-4)
            return float(
                np.sum(shape)
                + area_weight * np.sum(np.log(area / target) ** 2)
                + 1e5 * np.sum(fold**2 + bad**2)
            )

        if initial is None:
            initial = snyder_weights(n, corners3)
        pick = {v: m for m, v in enumerate(verts)}
        theta0 = np.log(
            np.clip(np.array([initial[pick[rep]] for rep in reps]), 1e-9, None)
        ).ravel()

        res = minimize(energy, theta0, method="L-BFGS-B",
                       options={"maxiter": 4000, "maxfun": 40000, "ftol": 1e-14,
                                "gtol": 1e-10, "eps": 1e-6})
        res = minimize(energy, res.x, method="Powell",
                       options={"maxiter": 40000, "maxfev": 60000, "xtol": 1e-10,
                                "ftol": 1e-12})
        return cls(n, cls._expand(res.x, reps, how, stab, verts)), float(res.fun)


def spherical_triangle_area(v1, v2, v3):
    """Signed spherical excess, by Van Oosterom and Strackee."""
    v1, v2, v3 = np.atleast_2d(v1), np.atleast_2d(v2), np.atleast_2d(v3)
    num = np.sum(v1 * np.cross(v2, v3), axis=-1)
    den = 1.0 + np.sum(v1 * v2, -1) + np.sum(v2 * v3, -1) + np.sum(v3 * v1, -1)
    out = 2.0 * np.arctan2(num, den)
    return out if out.size > 1 else float(out[0])


def _subdivided_to_sphere(self, bary, face_vertex):
    """Locate the sub-triangle holding `bary` and interpolate inside it."""
    n = self.n
    x = np.asarray(bary, float) * n
    idx = np.floor(x).astype(int)
    idx = np.clip(idx, 0, n)
    s = idx.sum(axis=-1)

    # An exact grid point floors onto a vertex rather than a triangle; step the
    # largest coordinate back so it lands in the triangle below it.
    high = s >= n
    big = np.argmax(idx, axis=-1)
    for m in range(3):
        idx[..., m] -= (high & (big == m)).astype(int)
    s = idx.sum(axis=-1)
    frac = x - idx
    up = s == (n - 1)

    i, j, k = idx[..., 0], idx[..., 1], idx[..., 2]
    vid = self.vertex_id
    up_v = np.stack([vid[i + 1, j], vid[i, j + 1], vid[i, j]], axis=-1)
    dn_v = np.stack([vid[i + 1, j + 1], vid[i + 1, j], vid[i, j + 1]], axis=-1)
    which = np.where(up[..., None], up_v, dn_v)

    w_up = frac
    w_dn = np.stack([1.0 - frac[..., 2], 1.0 - frac[..., 1], 1.0 - frac[..., 0]], axis=-1)
    w = np.where(up[..., None], w_up, w_dn)

    take = np.take_along_axis(face_vertex, which[..., None].repeat(3, axis=-1), axis=-2)
    v = np.sum(w[..., None] * take, axis=-2)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def _subdivided_post_init(self):
    n = self.n
    self.vertex_id = np.zeros((n + 2, n + 2), np.int32)
    for m, (i, j, k) in enumerate(self.verts):
        self.vertex_id[i, j] = m


SubdividedFace.to_sphere = _subdivided_to_sphere
_orig_init = SubdividedFace.__init__


def _sub_init(self, n, weights):
    _orig_init(self, n, weights)
    _subdivided_post_init(self)


SubdividedFace.__init__ = _sub_init


def snyder_weights(n, face3d):
    """Corner weights whose Snyder images sit at the regular mesh points.

    Snyder is exactly equal-area, so every sub-triangle of the regular planar
    mesh lands on a sub-triangle of exactly the right spherical area. That
    makes it the natural starting point: the area term of the energy begins at
    zero and the optimiser spends its effort on shape.
    """
    verts, _ = subdivision_vertices(n)
    bary = np.array(verts, float) / n
    face3d = np.asarray(face3d)
    snyder = SnyderTriangle(face3d)
    v = snyder_face(bary, np.broadcast_to(face3d, (len(bary), 3, 3)), snyder)
    w = np.linalg.solve(face3d.T[None, :, :], v[:, :, None])[:, :, 0]
    w = np.clip(w, 0.0, None)
    return w / w.sum(axis=1, keepdims=True)


_FACE_CACHE = {}


def optimised_face(n=16, area_weight=4.0, cache_dir=None):
    """The subdivided face, optimised once and then cached.

    The mesh depends only on the geometry of a tetrahedron face, not on how
    the tetrahedron is turned against the world, so the search runs once and
    every attitude reuses it.
    """
    from pathlib import Path

    key = (n, float(area_weight))
    if key in _FACE_CACHE:
        return _FACE_CACHE[key]
    if cache_dir is None:
        cache_dir = Path(__file__).resolve().parent / "build"
    cache_dir = Path(cache_dir)
    path = cache_dir / f"authagraph_mesh_n{n}_w{area_weight:g}.npy"
    if path.exists():
        face = SubdividedFace(n, np.load(path))
    else:
        verts, faces = tetrahedron()
        face, _ = SubdividedFace.optimise(n, verts[list(faces[0])], area_weight=area_weight)
        cache_dir.mkdir(parents=True, exist_ok=True)
        np.save(path, face.weights)
    _FACE_CACHE[key] = face
    return face
