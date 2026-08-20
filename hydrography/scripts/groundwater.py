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
# Imported for its side effect: `_paths` is what puts `lib/` on the path, so
# `orogen` below is unimportable without it. It reads as an unused import and
# is not one -- removing it on that basis broke this module's entry point.
import _paths  # noqa: E402,F401  # noqa
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
# Convergence of the outer iteration, DECLARED BEFORE THE FIRST RUN of the
# constant-transmissivity model. The transmissivity no longer depends on the
# head, so the matrix is fixed and this is a box-constrained linear
# complementarity problem: an active-set method on one terminates exactly rather
# than approaching a fixed point, and the bar is set accordingly tight. It is
# still the water balance residual on free cells over total land recharge, for
# the reason the previous model made unavoidable -- a bar on the head step can
# be satisfied by a shrinking step rather than by convergence.
# Tighter than CLOSURE_TOLERANCE, and it has to be: a pinned cell left with a
# negative balance is clipped to zero on the way out, and that invented water is
# EXACTLY the closure error. A leak bar looser than the closure bar lets a solve
# report convergence and then fail the conservation law it just broke, which is
# what 1e-8 did -- leak 2.53e-09, closure 2.47e-09 against 1e-10.
RESIDUAL_TOLERANCE = 1e-12
# The active set must also have stopped moving. Both, not either.
FLIP_TOLERANCE = 0.0
# Complementarity is tested as the MASS a pinned cell would have to invent, over
# total recharge, against the same bar as the residual -- not as a count of
# cells with a negative balance. The count was tried first and is a proxy that
# round-off defeats: four cells out of 2.5 million held it above zero forever at
# balances of order 1e-18 m3/s, on a planet recharging at 1.3e6, so the solve
# ran to its pass limit while the closure it was protecting passed at 1.7e-16.
# The intent is unchanged and is the one that matters -- no water may be
# invented by clipping a negative seepage to zero -- and this measures it
# directly instead of through a proxy.
# CROSS-SCHEME AGREEMENT, declared before running it. With transmissivity
# independent of the head the two code paths solve the IDENTICAL linear system,
# so this is not an independent identity the way it would be for a nonlinear
# model -- it is a check that two assemblies of the same operator agree, and it
# is held to round-off rather than to discretisation error.
SCHEME_HEAD_RELATIVE = 1e-9

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


def transmissivity(k0_m_s, thickness_m: float):
    """`T = K D`, m2/s. Conductivity over a constant saturated thickness.

    THE WHOLE DEPTH MODEL, and its flatness is deliberate. `K` and `D` come from
    one source at one scale: Gleeson et al. (2011) put their permeabilities at
    5-100 km and say the lithology maps carrying them "represent the shallow
    subsurface (on the order of 100 m)". A region here is about 15 km, inside
    that range.

    WHAT THIS REPLACED. Fan et al. (2007) make conductivity decay exponentially
    with depth at an e-folding length set by terrain slope, and that was used
    here first. It has no value at this resolution: `exp(h/f)` is convex, so a
    cell mean over a water table varying within the cell by `sigma` carries
    `exp(sigma^2 / 2 f^2)`, and 100 m of sub-grid relief against the 0.95 m `f`
    that Fan's curve reaches on steep bedrock makes that `exp(5000)`. That is
    not a correction to apply but a statement that the parameterisation cannot
    be evaluated on a 15 km cell at all. The config carries the argument; the
    two solvers that failed trying are in `notes/water-table-convergence.md`.

    The cost is real and is not hidden: no slope dependence of aquifer depth,
    and no thinning of transmissivity as the water table falls. This is the
    coarser model, taken because it is the one the sources support here.
    """
    return k0_m_s * float(thickness_m)


def face_transmissivity(t_cell, gfac, src, dst, is_ocean):
    """Face conductance, averaging the two cells except at the coast.

    THE COAST TAKES THE LAND SIDE'S VALUE, not the average of the two. The ocean
    is a fixed-head boundary, not a conducting cell, and the rock the water
    actually travels through to reach it is the coastal land cell's own. Ocean
    regions carry the export's `water` rock class, which has no permeability by
    construction, so averaging made every coastal face conduct at half of
    nothing: the sea stopped being a drain, the system lost its only Dirichlet
    anchor, and the direct solve reported the matrix exactly singular.
    """
    coast = is_ocean[src] ^ is_ocean[dst]
    land_side_t = np.where(is_ocean[src], t_cell[dst], t_cell[src])
    return gfac * np.where(coast, land_side_t,
                           0.5 * (t_cell[src] + t_cell[dst]))


def geom_divergence(n, src, dst, face_flux):
    """Net INFLOW per cell from a signed per-face flux.

    The convention throughout: a positive face flux is flow from `dst` into
    `src`, so this returns what each cell gains.
    """
    out = np.zeros(n)
    np.add.at(out, src, face_flux)
    np.add.at(out, dst, -face_flux)
    return out


# ---------------------------------------------------------------------------
# The complementarity solve
# ---------------------------------------------------------------------------

def solve(export: Export, geom: Geometry, *, k0_m_s, thickness_m, recharge_m_s,
          surface_m, conductive, sea_level_m=0.0, max_outer=200,
          start_all_free=False, verbose=True):
    """Steady-state head, seepage and face fluxes.

    With `T = K D` the transmissivity does not depend on the head, so the matrix
    is FIXED and the only thing left to iterate is the active set:

        sum_j Trans_ij (h_j - h_i) + R_i A_i = S_i,   h_i <= z_i,  S_i >= 0

    That is a box-constrained linear complementarity problem, and an active-set
    method on one terminates exactly. There is no relaxation, no trust region,
    no damping and no Picard loop here, and their absence is the point: every
    one of them existed to control a fixed-point iteration on an exponential
    transmissivity that could not be evaluated on this mesh in the first place.

    Ocean regions are a fixed-head boundary at sea level. Fan et al. (2013) find
    sea level the dominant driver of water table depth at global scale, and this
    world has one, so that boundary is not an approximation to something else.

    `conductive` is the land mask minus the regions whose lithology has no
    assigned permeability. Those leave the network rather than being given a
    number, and their recharge returns as seepage, which is what the
    surface-only balance already does with it.
    """
    n = export.n_regions
    src, dst = geom.src, geom.dst
    gfac = geom.geom
    area_m2 = geom.volume_area_m2

    is_ocean = export.surface_class != LAND
    in_domain = conductive | is_ocean
    live = in_domain[src] & in_domain[dst]
    src, dst, gfac = src[live], dst[live], gfac[live]

    t_cell = np.where(conductive, transmissivity(k0_m_s, thickness_m), 0.0)
    trans = face_transmissivity(t_cell, gfac, src, dst, is_ocean)

    rowsum = np.zeros(n)
    np.add.at(rowsum, src, trans)
    np.add.at(rowsum, dst, trans)

    supply = recharge_m_s * area_m2               # m3/s per cell
    total_supply = float(supply[conductive].sum())

    # Start with every cell at the surface, which is the zero-permeability
    # solution and is always feasible, then release. Growing the free set from
    # nothing keeps it a fraction of the mesh; starting all-free asks the solver
    # to factor the whole 2.5 million cells before a single cell is pinned.
    head = surface_m.copy()
    head[is_ocean] = sea_level_m
    # `start_all_free` walks the opposite trajectory, for the uniqueness check.
    # The solution of this complementarity problem is unique, so where it lands
    # cannot depend on where it started.
    free = conductive.copy() if start_all_free else np.zeros(n, bool)
    anchor = np.zeros(n, bool)
    anchor_failed = np.zeros(n, bool)
    result, trace = {}, []

    for outer in range(max_outer):
        # A cell with no conducting face has no equation: nothing can carry its
        # recharge away, so it stands at the surface and seeps.
        stranded = free & ~is_ocean & (rowsum <= 0)
        free[stranded] = False
        head[stranded] = surface_m[stranded]

        # RELEASE before solving. A pinned cell whose neighbours already draw
        # more water out of it than its recharge supplies cannot stand at the
        # surface. An anchor is never released; see below.
        flux = trans * (head[dst] - head[src])
        seep = supply + geom_divergence(n, src, dst, flux)
        # AN ANCHOR MAY BE RELEASED, and once released it is never chosen
        # again. Anchors were permanent first, on the reasoning that freeing one
        # returns its block to the singular state it was rescued from. That is
        # true only if the block still needs an anchor, and it makes an anchor
        # that turns out infeasible unfixable: four of them held 2.5e-9 of the
        # planet's recharge as invented water, forever, because the release test
        # was forbidden to touch them.
        #
        # Retiring a failed candidate is what makes this terminate. A block with
        # no outlet must dispose of its own recharge internally, so its balances
        # sum positive and at least one of its cells genuinely seeps; each
        # release strikes one candidate off, so the search cannot cycle and
        # cannot exhaust the block.
        release = conductive & ~free & (rowsum > 0) & (seep < 0)
        free[release] = True
        anchor_failed |= release & anchor
        anchor &= ~release

        unknown = free & ~is_ocean
        m = int(unknown.sum())
        if m == 0:
            if verbose:
                print("  every cell is at the surface; nothing to solve")
            result.update(outer_iterations=outer + 1, converged=True,
                          final_residual=0.0)
            break
        idx = np.full(n, -1, dtype=np.int64)
        idx[unknown] = np.arange(m)

        # Every free block needs a fixed head somewhere. A block connected only
        # to other free cells is a pure Neumann problem with non-negative
        # recharge, so it has no solution and the direct solve returns NaN for
        # the whole system. Physically it is a block with no outlet, which fills
        # until its shallowest cell seeps, so pinning that cell is what the
        # water does. Anchors are permanent: released, the block returns to the
        # singular state it was rescued from and the solve cycles.
        both = unknown[src] & unknown[dst]
        anchored = np.zeros(n, bool)
        edge_one = unknown[src] ^ unknown[dst]
        if edge_one.any():
            anchored[np.where(unknown[src[edge_one]],
                              src[edge_one], dst[edge_one])] = True
        if not np.all(anchored[unknown]):
            from scipy.sparse.csgraph import connected_components

            sub = sp.coo_matrix((np.ones(int(both.sum())),
                                 (idx[src[both]], idx[dst[both]])), shape=(m, m))
            ncomp, label = connected_components(sub, directed=False)
            ok = np.zeros(ncomp, bool)
            ok[label[idx[anchored & unknown]]] = True
            orphan = ~ok[label]
            # A BLOCK WITH NO RECHARGE IS DRY, and needs no anchor.
            #
            # An orphan block has no conducting face to anything outside it, so
            # the only water it can hold is what falls on it. Where that is zero
            # there is nothing to dispose of, no seepage, and no water table:
            # the block is dry and belongs out of the network, exactly as an
            # unassigned lithology does.
            #
            # Anchoring one of its cells instead is not merely unnecessary, it
            # is what generated the last of the mass error. Pinning a cell at
            # the surface inside a block whose terrain is not flat drives a flux
            # between the block's own cells, and the pinned cell absorbs the
            # resulting imbalance as a negative seepage that is then clipped.
            # All four cells that held this solve short of its bar had
            # `supply` of exactly zero.
            if orphan.any():
                block_supply = np.zeros(ncomp)
                np.add.at(block_supply, label[orphan],
                          supply[np.flatnonzero(unknown)[orphan]])
                arid = orphan & (block_supply[label] <= 0.0)
                if arid.any():
                    cells = np.flatnonzero(unknown)[arid]
                    conductive = conductive.copy()
                    conductive[cells] = False
                    free[cells] = False
                    head[cells] = surface_m[cells]
                    t_cell = np.where(conductive,
                                      transmissivity(k0_m_s, thickness_m), 0.0)
                    trans = face_transmissivity(t_cell, gfac, src, dst, is_ocean)
                    rowsum = np.zeros(n)
                    np.add.at(rowsum, src, trans)
                    np.add.at(rowsum, dst, trans)
                    total_supply = float(supply[conductive].sum())
                    if verbose:
                        print(f"    {cells.size:,} cells in "
                              f"{int((block_supply <= 0).sum())} orphan blocks "
                              f"take no recharge at all; marking them dry")
                    continue

                # WHICH cell to anchor, and it is not a free choice. A block
                # with no external face must dispose of all its own recharge
                # internally, so its cells' balances sum to something positive
                # and at least one of them genuinely seeps. Anchoring THAT one
                # is the only choice that stays feasible; anchoring the cell
                # with the shallowest water table, which is what this did
                # first, picked a cell that immediately wanted releasing and
                # could not be, so four of them sat infeasible for three
                # hundred passes and the solve never converged.
                where = np.flatnonzero(unknown)[orphan]
                bal_now = (supply + geom_divergence(
                    n, src, dst, trans * (head[dst] - head[src])))[where]
                # Candidates that already failed as anchors are pushed to the
                # back rather than removed, so a block whose every cell has
                # failed still gets one and the solve reports the leak instead
                # of going singular.
                score = np.where(anchor_failed[where], -np.inf, bal_now)
                best = np.full(ncomp, -np.inf)
                np.maximum.at(best, label[orphan], score)
                pick = where[score >= best[label[orphan]]]
                free[pick] = False
                head[pick] = surface_m[pick]
                anchor[pick] = True
                if verbose:
                    print(f"    anchored {int((~ok).sum()):,} orphan blocks by "
                          f"pinning {pick.size:,} cells")
                continue

        rows = np.concatenate([idx[src[both]], idx[dst[both]],
                               idx[src[both]], idx[dst[both]]])
        cols = np.concatenate([idx[dst[both]], idx[src[both]],
                               idx[src[both]], idx[dst[both]]])
        vals = np.concatenate([-trans[both], -trans[both],
                               trans[both], trans[both]])
        rhs = np.zeros(m)
        np.add.at(rhs, idx[unknown], supply[unknown])
        if edge_one.any():
            u_side = np.where(unknown[src[edge_one]], src[edge_one], dst[edge_one])
            p_side = np.where(unknown[src[edge_one]], dst[edge_one], src[edge_one])
            rows = np.concatenate([rows, idx[u_side]])
            cols = np.concatenate([cols, idx[u_side]])
            vals = np.concatenate([vals, trans[edge_one]])
            np.add.at(rhs, idx[u_side], trans[edge_one] * head[p_side])

        A = sp.coo_matrix((vals, (rows, cols)), shape=(m, m)).tocsr()
        if not np.all(A.diagonal() > 0):
            raise SystemExit(
                f"{int((A.diagonal() <= 0).sum())} free cells are isolated; "
                "the operator is singular and the solve is invalid")
        # A direct factorisation, not conjugate gradients. Conductivity still
        # spans the 3.4 orders of magnitude between Gleeson's classes, which
        # diagonally preconditioned CG handles badly, and the free set is small
        # enough that a direct solve is comfortable.
        x = spsolve(A.tocsc(), rhs, use_umfpack=False)
        if not np.all(np.isfinite(x)):
            raise SystemExit(
                f"the direct solve returned {int((~np.isfinite(x)).sum()):,} "
                f"non-finite heads at pass {outer}; do not use this result")
        head[unknown] = x

        # No relaxation and no step limit: the operator does not depend on the
        # answer, so a full step invalidates nothing.
        over = unknown & (head > surface_m)
        head[over] = surface_m[over]
        free[over] = False

        flux = trans * (head[dst] - head[src])
        bal = supply + geom_divergence(n, src, dst, flux)
        still = free & ~is_ocean
        residual = float(np.abs(bal)[still].sum() / max(total_supply, 1e-30))

        # COMPLEMENTARITY, tested on the head just produced. A converged LCP
        # needs BOTH halves: the balance closes on free cells, AND no pinned
        # cell is being drawn below the surface it is pinned to. Only the first
        # was checked, and because the break came straight after a solve that
        # had moved the head, pinned cells left infeasible by that last solve
        # were never re-examined -- their negative seepage was clipped to zero
        # on the way out, which does not discard water, it INVENTS it. That was
        # 740 m3/s of it, and closure caught it at 5.6e-4 where the free-cell
        # residual read 3e-18.
        pinned_now = conductive & ~free
        infeasible = int((pinned_now & (bal < 0)).sum())
        leak = float(-bal[pinned_now & (bal < 0)].sum())
        leak_frac = leak / max(total_supply, 1e-30)

        flips = int(over.sum()) + int(release.sum())
        trace.append([outer, m, flips, residual, infeasible, leak_frac])
        if verbose:
            print(f"  pass {outer:3d}  free {m:>9,}  released "
                  f"{int(release.sum()):>7,}  pinned {int(over.sum()):>7,}"
                  f"  residual {residual:10.3e}  leak {leak_frac:10.3e}"
                  f"  ({infeasible:,} cells)")
        if (residual < RESIDUAL_TOLERANCE and flips <= FLIP_TOLERANCE * m
                and leak_frac < RESIDUAL_TOLERANCE):
            result.update(outer_iterations=outer + 1, converged=True,
                          final_residual=residual, seepage_leak_m3_s=leak,
                          seepage_leak_fraction=leak_frac,
                          infeasible_pinned_cells=infeasible)
            break
    else:
        result.update(outer_iterations=max_outer, converged=False,
                      final_residual=residual, seepage_leak_m3_s=leak,
                      seepage_leak_fraction=leak_frac,
                      infeasible_pinned_cells=infeasible)
        if verbose:
            print(f"  DID NOT CONVERGE in {max_outer} passes: residual "
                  f"{residual:.3e}, leak {leak_frac:.3e}, both against "
                  f"{RESIDUAL_TOLERANCE:.0e}; {infeasible:,} pinned cells "
                  f"infeasible")
            bad = np.flatnonzero(pinned_now & (bal < 0))
            for i in bad[:8]:
                print(f"    region {i}: supply {supply[i]:.6e} m3/s  "
                      f"balance {bal[i]:.6e}  anchor {bool(anchor[i])}  "
                      f"failed_anchor {bool(anchor_failed[i])}  "
                      f"rowsum {rowsum[i]:.3e}  depth {surface_m[i]-head[i]:.2f} m")
    result["residual_trace"] = trace

    depth = np.clip(surface_m - head, 0.0, None)
    flux = trans * (head[dst] - head[src])
    divergence = geom_divergence(n, src, dst, flux)
    pinned = conductive & ~free
    seepage = np.zeros(n)
    seepage[pinned] = supply[pinned] + divergence[pinned]
    excluded = (export.surface_class == LAND) & ~conductive
    seepage[excluded] = supply[excluded]
    # At convergence every pinned cell is feasible, so this clip removes
    # nothing; what it WOULD remove is recorded so that a run which stops short
    # cannot quietly balance its books against it.
    result["seepage_clipped_m3_s"] = float(-seepage[seepage < 0].sum())
    seepage = np.clip(seepage, 0.0, None)

    result.update({
        "head_m": head, "depth_m": depth, "seepage_m3_s": seepage,
        "face_flux_m3_s": flux, "src": src, "dst": dst,
        "pinned": pinned, "excluded": excluded,
        "transmissivity_m2_s": t_cell,
        "at_transmissivity_floor": np.zeros(n, bool),
        # What the ocean GAINS. `geom_divergence` returns net inflow, so this is
        # a plus; it was a minus, which made closure subtract the coastal
        # discharge instead of adding it.
        "ocean_outflow_m3_s": float(divergence[is_ocean].sum()),
    })
    return result


def closure(result, recharge_m_s, area_m2, export: Export) -> dict:
    """Recharge in against discharge out. A conservation law, not a comparison.

    At equilibrium every drop of recharge leaves through exactly one of two
    doors: it seeps back to the surface somewhere, or it crosses the coast into
    the ocean boundary. If the two disagree the discretisation is wrong, and no
    amount of agreement elsewhere rescues it.
    """
    land = export.surface_class == LAND
    supply = float((recharge_m_s[land] * area_m2[land]).sum())
    out = float(result["seepage_m3_s"].sum()) + result["ocean_outflow_m3_s"]
    rel = abs(out - supply) / max(abs(supply), 1e-30)
    return {
        "recharge_m3_s": supply,
        "discharge_m3_s": out,
        "seepage_m3_s": float(result["seepage_m3_s"].sum()),
        "to_ocean_m3_s": result["ocean_outflow_m3_s"],
        "relative_residual": rel,
        "tolerance": CLOSURE_TOLERANCE,
        "passes": rel < CLOSURE_TOLERANCE,
    }

def terrain_following(export: Export, geom: Geometry, *, k0_m_s,
                      thickness_m, surface_m, conductive, sea_level_m=0.0):
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

    t_cell = np.where(conductive, transmissivity(k0_m_s, thickness_m), 0.0)
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
        "at_transmissivity_floor": np.zeros(export.n_regions, bool),
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
