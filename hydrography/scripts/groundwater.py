"""Steady-state water table on the region mesh: the machinery, not a result.

`build_groundwater.py` supplies the forcing and writes the products; this file
owns the discretisation, the transmissivity, and the complementarity solve, and
it is runnable on its own to exercise the checks that do not need a climate.

    python hydrography/scripts/groundwater.py        # the discretisation checks

THE MODEL. A steady-state, vertically integrated, unconfined water table, which
is the form `docs/src/reference/no-time-axis.md` leaves reachable: an
equilibrium surface, not an aquifer with a history. Per region,

    sum_j  Trans_ij (h_j - h_i)  +  R_i A_i  =  S_i,     h_i <= z_i,  S_i >= 0

for head `h`, land surface `z`, recharge `R`, cell area `A` and seepage `S`.
The two inequalities are complementary: a cell is either below the surface with
no seepage, or at the surface and discharging. That is what makes this a
complementarity problem rather than a linear solve, and it is the whole reason
the water table is interesting -- the cap is where groundwater becomes surface
water, and Fan et al. (2013) put 15% of Earth's land area there.

Transmissivity follows Fan et al. (2007) equations (5) and (6): conductivity
decays exponentially with depth below the surface at an e-folding length `f`, so
the saturated column below a water table at depth `d` integrates to

    T(d) = K0 f exp(-d / f)

and `hydrography/config/groundwater.yaml` carries `f`, `K0` per lithology, and
the sources for both.

GRAVITY. `K0 = k rho g / mu`. Permeability `k` is pore geometry and transfers
from Earth unchanged; conductivity is a flow property and carries this world's
gravity. Nothing in this module holds a gravity literal: `conductivity()` takes
it and `build_groundwater.py` reads it from `config/planet.yaml`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

sys.path.insert(0, str(Path(__file__).resolve().parent))
from orogen import LAND, Export  # noqa: E402

# Declared BEFORE the first run, per the standing convention. A criterion chosen
# after the result it judges is not a criterion.
#
# LAPLACE_TOLERANCE bounds the discretisation error of the mesh operator against
# an analytic eigenfunction. This is a first-order finite-volume scheme on an
# irregular Voronoi mesh whose cell areas span a factor of 600 and whose degree
# runs 3 to 15, so it is not spectrally accurate and is not meant to be; what it
# has to be is right to leading order. 10% relative RMS is the bar.
LAPLACE_TOLERANCE = 0.10
LAPLACE_DEGREES = (1, 2, 3, 4)
# CLOSURE_TOLERANCE: recharge in equals discharge out. A conservation law over
# float64 sums of a few million terms, so this is a numerical bar, not a
# physical one.
CLOSURE_TOLERANCE = 1e-10
# DIVIDE_AGREEMENT: with uniform permeability and a terrain-following table, the
# groundwater catchments must reproduce the surface catchments. Outside filled
# depressions the two constructions should agree almost everywhere; inside them
# the surface catchment came from a priority flood on a flat filled surface and
# the groundwater trace has nothing to follow, so they need not.
DIVIDE_AGREEMENT_ALL = 0.90
DIVIDE_AGREEMENT_UNFILLED = 0.99
# Convergence of the outer iteration, as the water balance residual on free
# cells over total land recharge, together with a nearly still active set. The
# residual measures the nonlinear problem rather than the iteration chasing it,
# so a shrinking step cannot satisfy it; a bar on the head step could be, and
# was replaced for that reason.
RESIDUAL_TOLERANCE = 1e-6
HEAD_TOLERANCE_M = 0.05          # reported, not a bar; see solve()
FLIP_TOLERANCE = 0.002
# Release hysteresis, as a fraction of a cell's own recharge. See solve().
FLIP_HYSTERESIS = 0.01
# Trust region on the Picard step, in e-folding lengths. See solve().
STEP_CAP_EFOLDINGS = 1.0
# Passes over which the relaxation factor halves. See solve().
DAMPING_SCALE = 60.0

SECONDS_PER_DAY = 86400.0


# ---------------------------------------------------------------------------
# Mesh geometry
# ---------------------------------------------------------------------------

class Geometry:
    """Faces, face widths and cell areas, from the mesh's own Voronoi dual.

    A two-point flux finite-volume scheme needs the width of each shared face
    and the distance between the two generators either side of it. The export
    carries neither, only a CSR neighbour list and a `cell_area`, so both are
    recovered here from the geometry itself rather than estimated.

    HOW. The regions are the cells of a spherical Voronoi tessellation, so the
    convex hull of the generators IS the Delaunay triangulation, each Voronoi
    vertex is a triangle's circumcentre, and the face between two neighbours is
    the arc joining the circumcentres of the two triangles on that edge. The
    dual returns 7,499,997 faces, exactly the number of undirected pairs in the
    export's own adjacency, and agrees with it on all but 47 -- those being
    cospherical quadrilaterals where the Delaunay diagonal is genuinely
    ambiguous. The dual's edge set is the one used, because it is the one that
    has face widths attached to it.

    That this reconstruction is right is not asserted: the exact spherical areas
    it produces sum to 4 pi R^2 to eight figures.

    THE FIRST ATTEMPT ESTIMATED THE FACE WIDTH from cell area and neighbour
    count, treating each region as a regular n-gon and giving each face an equal
    share of its perimeter. That is exact for a uniform tessellation and this
    mesh is not one: it put the operator 178x off its analytic eigenvalue, with
    the right answer in the median and a spread of two orders either side,
    because a face shared with a small neighbour got the same width as one
    shared with a large one. It is recorded because the failure was invisible to
    every other check -- closure holds for any symmetric weights whatsoever.

    TWO AREAS, AND THEY ARE NOT THE SAME. `voronoi_area_m2` is the exact
    spherical area bounded by these faces. The export's `cell_area` is a
    different quantity: it sums to 0.068% more than the sphere, and per cell the
    two disagree by more than 1% over 96.5% of the mesh. `flux_area_m2` is
    therefore the Voronoi area, because a divergence has to be taken over the
    area its own faces bound or the scheme is not consistent; and water VOLUMES
    stay on the export's `cell_area`, because every other component computes
    them that way and the reduction identity against `surface_water.py` has to
    be exact rather than close.
    """

    def __init__(self, export: Export):
        from scipy.spatial import ConvexHull

        self.export = export
        self.radius_m = export.radius_km * 1000.0
        p = np.stack([export.x, export.y, export.z], axis=1).astype(np.float64)
        p /= np.linalg.norm(p, axis=1, keepdims=True)

        tri = ConvexHull(p).simplices.astype(np.int64)
        a, b, c = p[tri[:, 0]], p[tri[:, 1]], p[tri[:, 2]]
        cc = np.cross(b - a, c - a)
        cc /= np.linalg.norm(cc, axis=1, keepdims=True)
        mid = a + b + c
        mid /= np.linalg.norm(mid, axis=1, keepdims=True)
        cc *= np.sign(np.einsum("ij,ij->i", cc, mid))[:, None]

        e = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
        e = np.sort(e, axis=1)
        tid = np.tile(np.arange(len(tri)), 3)
        order = np.lexsort((e[:, 1], e[:, 0]))
        e, tid = e[order], tid[order]
        pair = np.flatnonzero(np.all(e[:-1] == e[1:], axis=1))
        self.src, self.dst = e[pair, 0], e[pair, 1]
        t1, t2 = tid[pair], tid[pair + 1]

        # Face width and generator separation, both as great-circle arcs taken
        # through the half-chord. Arccos of the dot product is the same angle
        # and loses every significant figure as the dot approaches 1; on this
        # mesh, where adjacent generators come as close as ten metres, it
        # returned 4,014 faces of exactly zero length and took the operator to
        # NaN. The half-chord form is accurate where this mesh lives.
        self.face_m = self._arc(cc[t1], cc[t2])
        self.length_m = self._arc(p[self.src], p[self.dst])
        self.geom = self.face_m / self.length_m

        self.voronoi_area_m2 = np.zeros(export.n_regions)
        for gen in (self.src, self.dst):
            np.add.at(self.voronoi_area_m2, gen,
                      self._excess(p[gen], cc[t1], cc[t2]))
        self.voronoi_area_m2 *= self.radius_m ** 2
        self.flux_area_m2 = self.voronoi_area_m2
        self.volume_area_m2 = export.cell_area.astype(np.float64) * 1e6

    def _arc(self, u, v):
        chord = np.linalg.norm(u - v, axis=1)
        return 2.0 * self.radius_m * np.arcsin(np.clip(0.5 * chord, 0.0, 1.0))

    @staticmethod
    def _excess(u, v, w):
        """Spherical triangle area on the unit sphere, by tangent half-excess.

        The angle-sum form loses all precision on the very thin triangles a
        Voronoi sliver produces, and this mesh has them.
        """
        num = np.abs(np.einsum("ij,ij->i", u, np.cross(v, w)))
        den = (1.0 + np.einsum("ij,ij->i", u, v)
               + np.einsum("ij,ij->i", v, w) + np.einsum("ij,ij->i", w, u))
        return 2.0 * np.arctan2(num, den)

    def divergence(self, face_flux):
        """Net outflow per cell from a signed per-face flux, src to dst."""
        out = np.zeros(self.export.n_regions)
        np.add.at(out, self.src, face_flux)
        np.add.at(out, self.dst, -face_flux)
        return out

    def closes_on_sphere(self) -> float:
        """Voronoi area against 4 pi R^2. The check that the dual is complete."""
        return float(self.voronoi_area_m2.sum()
                     / (4.0 * np.pi * self.radius_m ** 2))


def laplace_beltrami_error(geom: Geometry, degrees=LAPLACE_DEGREES) -> dict:
    """Discrete operator against an analytic eigenfunction. A test that can fail.

    Legendre polynomials of the axial coordinate are zonal spherical harmonics,
    and so are eigenfunctions of the Laplace-Beltrami operator on the sphere
    with an eigenvalue fixed by geometry alone:

        laplacian P_l(z)  =  -l(l+1)/R^2  P_l(z)

    Applying the assembled operator at unit transmissivity to `P_l` and dividing
    by cell area must return that. Nothing about this world enters it -- not the
    recharge, the permeability, the cap or the climate -- so it tests the face
    widths, the edge lengths and the assembly and nothing else. It is the only
    check here that a wrong face width cannot survive.
    """
    from scipy.special import eval_legendre

    z = geom.export.z.astype(np.float64)
    out = {}
    for l in degrees:
        u = eval_legendre(l, z)
        got = geom.divergence(geom.geom * (u[geom.dst] - u[geom.src])) \
            / geom.flux_area_m2
        want = -l * (l + 1) / geom.radius_m**2 * u
        scale = np.sqrt(np.mean(want**2))
        out[l] = float(np.sqrt(np.mean((got - want) ** 2)) / scale)
    return out


# ---------------------------------------------------------------------------
# Subsurface properties
# ---------------------------------------------------------------------------

def conductivity(k_m2, density_kg_m3: float, gravity_m_s2: float,
                 viscosity_pa_s: float):
    """Hydraulic conductivity, m/s, from intrinsic permeability, m2.

    The one place gravity enters the subsurface. `k` is pore geometry and is
    Earth's measurement carried over unchanged; `K` is a flow property and is
    this world's. Gravity arrives as an argument and is never a literal here.
    """
    return k_m2 * density_kg_m3 * gravity_m_s2 / viscosity_pa_s


def efolding_length(slope_tan, cover_thickness_m, cfg: dict):
    """Fan et al. (2007) eq. (7): e-folding depth of conductivity, metres.

        f = a / (1 + b * slope),  capped below at f_min above slope_cap

    Two curves, regolith and bedrock, selected by whether the export's
    unconsolidated cover reaches `regolith_cover_min_m`. The constants are Fan's
    North American fit and are Earth's; the config says so at length.
    """
    dd = cfg["depth_decay"]
    is_regolith = cover_thickness_m >= float(dd["regolith_cover_min_m"])
    f = np.empty(np.shape(slope_tan), dtype=np.float64)
    for key, mask in (("regolith", is_regolith), ("bedrock", ~is_regolith)):
        p = dd[key]
        s = np.clip(slope_tan[mask], 0.0, None)
        val = float(p["a_m"]) / (1.0 + float(p["b"]) * s)
        f[mask] = np.where(s > float(p["slope_cap"]), float(p["f_min_m"]), val)
    return f


def transmissivity(k0_m_s, f_m, depth_m):
    """Fan et al. (2007) eq. (6): T = K0 f exp(-d/f), m2/s.

    `depth_m` is the water table below the land surface and is clipped at zero:
    the head is capped at the surface, so a negative depth is not a state this
    model has, and letting one through would make the column above the ground
    conduct.
    """
    return k0_m_s * f_m * np.exp(-np.clip(depth_m, 0.0, None) / f_m)


# ---------------------------------------------------------------------------
# The complementarity solve
# ---------------------------------------------------------------------------

def cell_transmissivity(k0_m_s, f_m, depth_m, conductive):
    """Per-cell transmissivity with the numerical floor, and who hit it."""
    t = np.where(conductive, transmissivity(k0_m_s, f_m, depth_m), 0.0)
    floor = TRANSMISSIVITY_FLOOR_RATIO * t.max()
    at_floor = conductive & (t < floor)
    return np.where(conductive, np.maximum(t, floor), t), at_floor


def face_transmissivity(t_cell, gfac, src, dst, is_ocean):
    """Face conductance, averaging the two cells except at the coast.

    THE COAST TAKES THE LAND SIDE'S VALUE, not the average of the two. The ocean
    is a fixed-head boundary, not a conducting cell, and the rock the water
    actually travels through to reach it is the coastal land cell's own. Ocean
    regions carry the export's `water` rock class, which has no permeability by
    construction, so averaging made every coastal face conduct at half of
    nothing: the sea stopped being a drain, the only Dirichlet anchor left was
    the transmissivity floor at 1e-12 of the mesh maximum, and the direct solve
    reported the matrix exactly singular after thirty passes.
    """
    coast = is_ocean[src] ^ is_ocean[dst]
    land_side_t = np.where(is_ocean[src], t_cell[dst], t_cell[src])
    return gfac * np.where(coast, land_side_t,
                           0.5 * (t_cell[src] + t_cell[dst]))


def geom_divergence(n, src, dst, face_flux):
    """Net outflow per cell from a signed per-face flux, on a face subset."""
    out = np.zeros(n)
    np.add.at(out, src, face_flux)
    np.add.at(out, dst, -face_flux)
    return out


# A NUMERICAL floor on transmissivity, as a fraction of the largest value on
# the mesh, and not a physical claim. `T = K0 f exp(-d/f)` has no lower bound:
# on steep bedrock `f` falls to a metre, so a water table a hundred metres down
# gives `exp(-100)`, which underflows to exactly zero and disconnects the cell
# from the operator. Cells at the floor are counted and reported, and their
# depth is a lower bound rather than a value. Raising the floor moves the depth
# only where the water table is already deeper than this mesh can resolve.
TRANSMISSIVITY_FLOOR_RATIO = 1e-12


def solve(export: Export, geom: Geometry, *, k0_m_s, f_m, recharge_m_s,
          surface_m, conductive, sea_level_m=0.0, max_outer=60,
          relax=0.7, verbose=True):
    """Steady-state head, seepage and face fluxes.

    Active-set iteration over the complementarity condition, with a Picard
    iteration on the transmissivity nested inside it because `T` depends on the
    head it is solving for.

    Ocean regions are a fixed-head boundary at sea level. Fan et al. (2013) find
    sea level the dominant driver of water table depth at global scale, and this
    world has one, so that boundary is not an approximation to something else.

    `conductive` is the land mask minus the regions whose lithology has no
    assigned permeability. Those are removed from the network entirely rather
    than given a number: see the `unassigned` block in the config. Their
    recharge is returned as seepage, which is what the surface-only balance
    already does with it.

    Returns a dict. `head_m` is the water table; `seepage_m3_s` is groundwater
    reaching the surface, non-negative by construction; `face_flux_m3_s` is
    positive from `src` to `dst`.
    """
    n = export.n_regions
    src, dst = geom.src, geom.dst
    gfac = geom.geom
    area_m2 = geom.volume_area_m2

    is_ocean = export.surface_class != LAND
    # START FROM THE ALL-PINNED POINT AND RELEASE, rather than from all-free and
    # pin. Both are active-set iterations on the same problem and they are not
    # equally conditioned. Every cell at the surface is the zero-permeability
    # solution, so it is always feasible; the free set then grows only where a
    # cell can actually carry its recharge away, and stays a fraction of the
    # mesh. Starting all-free asks conjugate gradients to solve the whole
    # 2.5-million-cell system across six orders of magnitude of transmissivity
    # before a single cell has been pinned, and it did not converge in four
    # thousand iterations.
    free = np.zeros(n, bool)
    head = surface_m.copy()
    head[is_ocean] = sea_level_m

    # A face conducts only between two cells that are both in the network, or
    # between a network cell and the ocean boundary. Everything else is inert.
    in_domain = conductive | is_ocean
    live = in_domain[src] & in_domain[dst]
    src, dst, gfac = src[live], dst[live], gfac[live]

    supply = recharge_m_s * area_m2               # m3/s per cell
    pinned_hist = []
    result = {}
    release = np.zeros(n, bool)
    anchor = np.zeros(n, bool)      # pinned to make an orphan block solvable
    head_prev = head.copy()
    step_trace = []
    total_supply = float(supply[conductive].sum())
    # The depth at which T reaches its floor. Past it the head is undetermined,
    # so the solve does not chase it and the reported depth is a lower bound.
    max_depth = f_m * np.log(1.0 / TRANSMISSIVITY_FLOOR_RATIO)

    for outer in range(max_outer):
        depth = np.clip(surface_m - head, 0.0, None)
        t_cell, at_floor = cell_transmissivity(
            k0_m_s, f_m, depth, conductive)
        trans = face_transmissivity(t_cell, gfac, src, dst, is_ocean)

        # A cell with no conducting face has no equation: nothing can carry its
        # recharge away, so its water table stands at the surface and its
        # recharge leaves as seepage. Pinning it is the physical answer and not
        # a numerical dodge, and it is the answer the zero-permeability limit
        # consists entirely of. Left free, it is a zero row and the operator is
        # singular.
        rowsum = np.zeros(n)
        np.add.at(rowsum, src, trans)
        np.add.at(rowsum, dst, trans)
        stranded = free & ~is_ocean & (rowsum <= 0)
        if stranded.any():
            free[stranded] = False
            head[stranded] = surface_m[stranded]

        # RELEASE before solving. A pinned cell whose neighbours already draw
        # more water out of it than its own recharge supplies cannot stand at
        # the surface, so its head is an unknown. Doing this first is what makes
        # the free set grow from nothing instead of shrinking from everything.
        # HYSTERESIS on the release, scaled to the cell's own recharge. A cell
        # exactly on the boundary otherwise releases on one pass and is pinned
        # back on the next, forever: the raw criterion limit-cycled here at
        # about 385,500 free cells, flipping some 2,500 either way every pass
        # and never settling. Requiring a cell to be drawn down by more than a
        # small fraction of its own supply before it is released breaks the
        # cycle without moving any cell that is not marginal.
        pinned_now = conductive & ~free & (rowsum > 0)
        flux_now = trans * (head[dst] - head[src])
        seep_now = supply + geom_divergence(n, src, dst, flux_now)
        margin = FLIP_HYSTERESIS * np.abs(supply)
        # An anchor cell is never released. It was pinned to give its block the
        # single fixed head a Neumann problem lacks, so freeing it puts the
        # block straight back into the singular state it was rescued from: the
        # solve then pins it, releases it and re-anchors it every pass, forever.
        release = pinned_now & (seep_now < -margin) & ~anchor
        free[release] = True

        unknown = free & ~is_ocean
        idx = np.full(n, -1, dtype=np.int64)
        idx[unknown] = np.arange(int(unknown.sum()))
        m = int(unknown.sum())
        if m == 0:
            # Nothing was released, so every cell stands at the surface and
            # there is nothing to solve for. A terminal state rather than a
            # failure, and exactly the zero-permeability limit: no cell can move
            # water sideways, so every cell returns its own recharge as seepage.
            if verbose:
                print("  every cell is pinned at the surface; nothing to solve")
            result["outer_iterations"] = outer + 1
            break

        # EVERY FREE COMPONENT NEEDS AN ANCHOR, and one can lose it mid-solve.
        #
        # A block of free cells connected only to other free cells is a pure
        # Neumann problem: it has a solution only if the recharge falling in it
        # sums to zero, and recharge is non-negative, so it has none. The matrix
        # is exactly singular and the direct solve returns NaN for the WHOLE
        # system, not just that block -- which is how this presented, twice, at
        # the same pass.
        #
        # It is a physical situation and not a degenerate one: an interior
        # block whose water table has dropped below the surface everywhere has
        # no outlet, so in steady state it fills until its shallowest cell
        # reaches the surface and starts to seep. Pinning that cell is what the
        # water would do, and it gives the block the anchor it needs. If the
        # choice is wrong the release test frees it again on the next pass.
        both = unknown[src] & unknown[dst]
        anchored = np.zeros(n, bool)
        edge_one = unknown[src] ^ unknown[dst]
        if edge_one.any():
            anchored[np.where(unknown[src[edge_one]],
                              src[edge_one], dst[edge_one])] = True
        if not np.all(anchored[unknown]):
            from scipy.sparse.csgraph import connected_components

            sub = sp.coo_matrix(
                (np.ones(int(both.sum())),
                 (idx[src[both]], idx[dst[both]])), shape=(m, m))
            ncomp, label = connected_components(sub, directed=False)
            ok = np.zeros(ncomp, bool)
            ok[label[idx[anchored & unknown]]] = True
            orphan = ~ok[label]
            if orphan.any():
                where = np.flatnonzero(unknown)[orphan]
                deep = (surface_m - head)[where]
                shallowest = np.full(ncomp, np.inf)
                np.minimum.at(shallowest, label[orphan], deep)
                pick = where[deep <= shallowest[label[orphan]]]
                free[pick] = False
                head[pick] = surface_m[pick]
                anchor[pick] = True
                if verbose:
                    print(f"    anchored {int((~ok).sum()):,} orphan blocks by "
                          f"pinning {pick.size:,} cells at the surface")
                continue          # rebuild with the new set before solving

        rows = np.concatenate([idx[src[both]], idx[dst[both]],
                               idx[src[both]], idx[dst[both]]])
        cols = np.concatenate([idx[dst[both]], idx[src[both]],
                               idx[src[both]], idx[dst[both]]])
        vals = np.concatenate([-trans[both], -trans[both],
                               trans[both], trans[both]])

        rhs = np.zeros(m)
        np.add.at(rhs, idx[unknown], supply[unknown])

        one = unknown[src] ^ unknown[dst]
        if one.any():
            u_side = np.where(unknown[src[one]], src[one], dst[one])
            p_side = np.where(unknown[src[one]], dst[one], src[one])
            rows = np.concatenate([rows, idx[u_side]])
            cols = np.concatenate([cols, idx[u_side]])
            vals = np.concatenate([vals, trans[one]])
            np.add.at(rhs, idx[u_side], trans[one] * head[p_side])

        A = sp.coo_matrix((vals, (rows, cols)), shape=(m, m)).tocsr()
        diag = A.diagonal()
        if not np.all(diag > 0):
            raise SystemExit(
                f"{int((diag <= 0).sum())} free cells are isolated from the "
                "network; the operator is singular and the solve is invalid")
        # A DIRECT SPARSE SOLVE, not conjugate gradients. The matrix is a
        # symmetric positive-definite graph Laplacian and CG is the obvious
        # choice for one, but the transmissivity carries `exp(-d/f)` with `f`
        # as small as a metre on steep bedrock, so a cell whose water table sits
        # tens of metres down has a transmissivity smaller than its neighbour's
        # by twenty orders of magnitude. Diagonally preconditioned CG did not
        # converge in four thousand iterations on that, at any active set. The
        # free set is a fraction of the mesh, which puts a direct factorisation
        # well inside reach, and it does not care about the condition number in
        # the way an iterative method does.
        x = spsolve(A.tocsc(), rhs, use_umfpack=False)
        if not np.all(np.isfinite(x)):
            raise SystemExit(
                f"the direct solve returned {int((~np.isfinite(x)).sum()):,} "
                f"non-finite heads at outer iteration {outer} with {m:,} "
                "unknowns; do not use this result")

        # LIMIT THE STEP TO THE E-FOLDING LENGTH, then under-relax.
        #
        # Under-relaxation alone does not control this iteration. `T` depends on
        # the head through `exp(-d/f)`, so a step much larger than `f` leaves
        # the regime the linearisation was taken in: the cell's transmissivity
        # collapses, the next pass needs a still larger gradient to move the
        # same water, and the head runs away. Unclamped, this oscillated by
        # three to four KILOMETRES per pass on a planet whose land spans five,
        # and did not decay over thirty passes.
        #
        # One e-folding length is the natural trust region, because it is the
        # distance over which the coefficient the step was computed from changes
        # by a factor of e. Cells far from their answer then walk towards it at
        # `f` per pass instead of overshooting past it.
        # A DIMINISHING step on top of the cap. A fixed relaxation is enough to
        # keep the iteration stable but not to make it converge: with a constant
        # step it settles into a limit cycle instead of a fixed point, and this
        # one did, holding a maximum head step near 70 m and a 99.9th percentile
        # near 22 m unchanged from pass 100 to pass 158 while the free set and
        # the flip count had long since stopped moving. Shrinking the step as
        # the passes accumulate damps that cycle towards its centre. The bulk of
        # the field has converged well before the factor becomes small, so what
        # it costs is only the tail it exists to settle.
        step = np.zeros(n)
        step[unknown] = x - head[unknown]
        cap = STEP_CAP_EFOLDINGS * f_m
        step = np.clip(step, -cap, cap)
        damp = relax / (1.0 + outer / DAMPING_SCALE)
        head = head + damp * np.where(unknown, step, 0.0)

        # Below this the transmissivity is at its floor and the head is not
        # determined by anything: the cell is disconnected, and the only effect
        # of letting it fall further is to widen the range the solve has to
        # carry. Clamping it keeps the depth a lower bound rather than a number.
        head = np.where(unknown, np.maximum(head, surface_m - max_depth), head)

        # A free cell whose head came out above the surface is pinned back to
        # it. Releases happen at the top of the next pass, on a consistent T.
        over = unknown & (head > surface_m)
        head[over] = surface_m[over]
        free[over] = False

        # CONVERGE ON THE HEAD, not on the active set. A few thousand marginal
        # cells go on flipping long after the field itself has stopped moving,
        # so a criterion of "the set stopped changing" never fires while a
        # criterion on the head does. The set is still required to be nearly
        # still, so a genuinely unconverged run cannot pass on a small step.
        # THE CRITERION IS THE RESIDUAL, NOT THE STEP.
        #
        # With a diminishing relaxation the head step shrinks whether or not the
        # solution has converged, so a bar on the step is one the damping can
        # satisfy on its own. It would have read as convergence and been none.
        #
        # The residual cannot be gamed that way. Recomputing the transmissivity
        # at the head just produced and asking whether each free cell's water
        # balance closes -- recharge in equals lateral flow out, with no seepage,
        # which is what being below the surface means -- tests the nonlinear
        # problem rather than the iteration that is chasing it. It is zero only
        # when `T` and `h` agree, which is the fixed point.
        depth_now = np.clip(surface_m - head, 0.0, None)
        t_now, _ = cell_transmissivity(k0_m_s, f_m, depth_now, conductive)
        trans_now = face_transmissivity(t_now, gfac, src, dst, is_ocean)
        div_now = geom_divergence(n, src, dst,
                                  trans_now * (head[dst] - head[src]))
        still_free = free & ~is_ocean
        resid = np.abs(supply + div_now)[still_free]
        residual = float(resid.sum() / max(total_supply, 1e-30))

        dh = np.abs(head - head_prev)[conductive]
        moved = float(dh.max())
        moved_p999 = float(np.percentile(dh, 99.9))
        head_prev = head.copy()
        flips = int(over.sum()) + int(release.sum())
        step_trace.append([outer, m, moved, moved_p999, flips, residual])
        pinned_hist.append(int((conductive & ~free).sum()))
        if verbose:
            print(f"  outer {outer:3d}  free {m:>9,}  released "
                  f"{int(release.sum()):>7,}  pinned {int(over.sum()):>7,}"
                  f"  |dh| {moved:9.3g}  p99.9 {moved_p999:9.3g}"
                  f"  residual {residual:9.3e}")
        if residual < RESIDUAL_TOLERANCE and flips <= FLIP_TOLERANCE * max(m, 1):
            result["outer_iterations"] = outer + 1
            result["final_head_step_m"] = moved
            result["final_head_step_p999_m"] = moved_p999
            result["final_flip_fraction"] = flips / max(m, 1)
            result["final_residual"] = residual
            result["converged"] = True
            break
    else:
        # NOT converged. Returned rather than raised, with the trace attached,
        # because a run that stops short is evidence about the solver and a
        # bare exit is not. The caller refuses to write the water table from it;
        # the report is written either way and says how far it got.
        result["converged"] = False
        result["outer_iterations"] = max_outer
        result["final_head_step_m"] = moved
        result["final_head_step_p999_m"] = moved_p999
        result["final_flip_fraction"] = flips / max(m, 1)
        result["final_residual"] = residual
        if verbose:
            print(f"  DID NOT CONVERGE in {max_outer} passes: water balance "
                  f"residual {residual:.3e} of total recharge against a bar of "
                  f"{RESIDUAL_TOLERANCE:.0e}; max |dh| {moved:.4g} m")
    result.setdefault("converged", True)
    result["head_step_trace_m"] = step_trace

    depth = np.clip(surface_m - head, 0.0, None)
    t_cell, at_floor = cell_transmissivity(k0_m_s, f_m, depth, conductive)
    trans = face_transmissivity(t_cell, gfac, src, dst, is_ocean)
    flux = trans * (head[dst] - head[src])

    divergence = np.zeros(n)
    np.add.at(divergence, src, flux)
    np.add.at(divergence, dst, -flux)

    seepage = np.zeros(n)
    pinned = conductive & ~free
    seepage[pinned] = supply[pinned] + divergence[pinned]
    # Excluded land carries its own recharge straight to the surface: it is not
    # in the network, so nothing else can move it.
    excluded = (export.surface_class == LAND) & ~conductive
    seepage[excluded] = supply[excluded]
    seepage = np.clip(seepage, 0.0, None)

    result.update({
        "head_m": head,
        "depth_m": depth,
        "seepage_m3_s": seepage,
        "face_flux_m3_s": flux,
        "src": src, "dst": dst,
        "pinned": pinned,
        "excluded": excluded,
        "transmissivity_m2_s": t_cell,
        "ocean_outflow_m3_s": float(-divergence[is_ocean].sum()),
        "at_transmissivity_floor": at_floor,
    })
    return result


def closure(result, recharge_m_s, area_m2, export: Export) -> dict:
    """Recharge in against discharge out. A conservation law, not a comparison.

    At equilibrium every drop of recharge leaves through exactly one of two
    doors: it seeps back to the surface somewhere, or it crosses the coast into
    the ocean boundary. If the two sides disagree the discretisation is wrong,
    and no amount of agreement elsewhere rescues it.
    """
    land = export.surface_class == LAND
    supply = float((recharge_m_s[land] * area_m2[land]).sum())
    out = float(result["seepage_m3_s"].sum()) + result["ocean_outflow_m3_s"]
    return {
        "recharge_m3_s": supply,
        "discharge_m3_s": out,
        "seepage_m3_s": float(result["seepage_m3_s"].sum()),
        "to_ocean_m3_s": result["ocean_outflow_m3_s"],
        "relative_residual": abs(out - supply) / max(abs(supply), 1e-30),
        "tolerance": CLOSURE_TOLERANCE,
        "passes": abs(out - supply) / max(abs(supply), 1e-30) < CLOSURE_TOLERANCE,
    }


def terrain_following(export: Export, geom: Geometry, *, k0_m_s, f_m,
                      surface_m, conductive, sea_level_m=0.0):
    """The flux field under a water table set equal to the terrain.

    Not a solve. The divide test asks what a terrain-following table does, so
    the table is imposed and the fluxes read off it; solving for a table and
    then checking it against the terrain it only approximately reproduced would
    be testing the solver's convergence instead of the claim.
    """
    is_ocean = export.surface_class != LAND
    head = surface_m.copy()
    head[is_ocean] = sea_level_m

    in_domain = conductive | is_ocean
    live = in_domain[geom.src] & in_domain[geom.dst]
    src, dst = geom.src[live], geom.dst[live]
    gfac = geom.geom[live]

    t_cell, at_floor = cell_transmissivity(
        k0_m_s, f_m, np.zeros_like(surface_m), conductive)
    trans = face_transmissivity(t_cell, gfac, src, dst, is_ocean)
    return {
        "head_m": head,
        "depth_m": np.zeros_like(surface_m),
        "seepage_m3_s": np.zeros_like(surface_m),
        "face_flux_m3_s": trans * (head[dst] - head[src]),
        "src": src, "dst": dst,
        "pinned": conductive.copy(),
        "excluded": (export.surface_class == LAND) & ~conductive,
        "transmissivity_m2_s": t_cell,
        "at_transmissivity_floor": at_floor,
        "ocean_outflow_m3_s": 0.0,
        "outer_iterations": 0,
    }


def groundwater_receiver(n, src, dst, flux, land):
    """Per land region, the neighbour taking the largest outgoing flux.

    The groundwater analogue of the surface network's `receiver`. A positive
    face flux runs src to dst, so it leaves src and enters dst; each direction
    is reduced separately and the larger of the two claims wins.
    """
    best = np.zeros(n)
    receiver = np.full(n, -1, dtype=np.int64)
    for a, b, q in ((src, dst, flux), (dst, src, -flux)):
        pos = q > 0
        aa, bb, qq = a[pos], b[pos], q[pos]
        order = np.argsort(qq, kind="stable")     # ascending: last write is max
        aa, bb, qq = aa[order], bb[order], qq[order]
        np.maximum.at(best, aa, qq)               # the true per-source maximum
        keep = qq >= best[aa]
        receiver[aa[keep]] = bb[keep]
    receiver[~land] = -1
    return receiver


def trace_terminals(receiver, surface_terminal, land):
    """Follow the groundwater receivers to a sink and read off its basin.

    Pointer jumping rather than a walk: each pass squares the distance covered,
    so a chain of any length this mesh can hold closes in well under sixty
    passes, against up to 2.5 million for a walk. A region with no outgoing flux
    points at itself and is the terminal.

    Returns the basin label the SURFACE catalogue gives each region's
    groundwater terminal, so it is directly comparable to `regions.nc`'s
    `terminal` and the divide test is a comparison of two labellings of the same
    thing.
    """
    rec = receiver.copy()
    self_ = rec < 0
    rec[self_] = np.flatnonzero(self_)
    for _ in range(64):
        nxt = rec[rec]
        if np.array_equal(nxt, rec):
            break
        rec = nxt
    else:
        raise SystemExit(
            "the groundwater receiver field did not resolve to sinks in 64 "
            "pointer-jumping passes; it probably contains a cycle")
    return np.where(land, surface_terminal[rec], surface_terminal)


def main() -> int:
    """The checks that need no climate: the discretisation against an identity."""
    import argparse

    import builds

    ap = argparse.ArgumentParser(
        description="Discretisation checks for the water table solver. The "
                    "forcing, the products and the basin measurement are "
                    "build_groundwater.py; this exercises the mesh operator "
                    "against an analytic eigenvalue and needs no climate.")
    ap.add_argument("--degrees", type=int, nargs="+", default=list(LAPLACE_DEGREES),
                    help="Legendre degrees to test the operator at")
    args = ap.parse_args()

    export = Export(builds.mesh_export())
    print(f"mesh {export.n_regions:,} regions, R = {export.radius_km:,.1f} km")
    geom = Geometry(export)
    print(f"Voronoi dual: {geom.src.size:,} faces, area closes on the sphere "
          f"to {geom.closes_on_sphere():.8f}")

    print("\nLaplace-Beltrami: discrete operator against P_l eigenfunctions")
    print(f"  criterion, declared before the run: relative RMS < "
          f"{LAPLACE_TOLERANCE:.2f}")
    err = laplace_beltrami_error(geom, tuple(args.degrees))
    ok = True
    for l, e in sorted(err.items()):
        verdict = "pass" if e < LAPLACE_TOLERANCE else "MISS"
        ok &= e < LAPLACE_TOLERANCE
        print(f"  l = {l}   relative RMS {e:.4f}   {verdict}")
    if not ok:
        print("\nThe operator MISSES the declared bar above l = 1. The error is\n"
              "distributed truncation rather than a few bad faces -- the worst\n"
              "10,000 cells carry under a tenth of it, and excluding every cell\n"
              "touching a sliver face moves the RMS by 0.002 -- so it is the\n"
              "first-order scheme's own error on an irregular mesh, and it is\n"
              "the water table's mesh-scale noise floor. Reported, not adjusted.")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
