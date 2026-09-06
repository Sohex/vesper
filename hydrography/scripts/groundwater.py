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
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components as _components
from scipy.sparse.linalg import spsolve, splu

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
# GW-8. WHAT THE OPERATOR'S ERROR IS AN ERROR IN, and the criteria that decide
# it. DECLARED BEFORE `--instrument` WAS FIRST RUN.
#
# `LAPLACE_TOLERANCE` above bounds a TRUNCATION error: the assembled operator
# applied to a field that is already known. GW-8 records that error at about
# 12% relative RMS above l = 1 and calls it "the mesh-scale noise floor of any
# field solved with this operator". Read that way it puts every effect smaller
# than 12% out of reach, and GW-18 and GW-24 both propose effects of a few per
# cent, so the reading decides whether either is worth building. A number doing
# that much work is measured rather than inherited.
#
# Nothing downstream reads a truncation error. What the water table hands its
# consumers is the SOLUTION of `L h = f`, and the two errors are different
# quantities joined by `L^-1`. An elliptic inverse damps the grid-scale part of
# a residual by the square of the scale ratio, so a truncation error that is
# grid-scale noise -- which is what GW-8 already measured it to be, distributed
# rather than carried by a few faces -- need not appear in the solution at
# anything like its own size. Two claims follow, each able to fail.
#
# CLAIM 1, THE TWO ERRORS BEHAVE DIFFERENTLY UNDER REFINEMENT. A two-point flux
# scheme on a mesh whose faces are perpendicular but whose face centres are not
# the midpoints of the generator segments has a leading truncation term set by
# mesh irregularity rather than by cell size, so it does not converge; the
# solution error does, at the scheme's own order. Measured as fitted orders in
# the cell size over a resolution sweep of meshes from Orogen's own generator,
# at the jitter the active build was generated with: the truncation order must
# come out BELOW `TRUNCATION_ORDER_MAX` and the solution order ABOVE
# `SOLUTION_ORDER_MIN`. If the two orders come out the same, the claim is
# refuted and the 12% stands exactly as GW-8 wrote it.
TRUNCATION_ORDER_MAX = 0.5
SOLUTION_ORDER_MIN = 1.0
# CLAIM 2, THE ERROR IS RELATIVE TO WHATEVER IS BEING MEASURED. The truncation
# figure is a RELATIVE RMS. If the same operator is asked for the RESPONSE to a
# transmissivity perturbation of relative size `a` -- which is exactly what
# GW-18 proposes -- and returns that response with an error that is again a
# fraction of the response rather than a fraction of the field, then an effect
# of size `a` is resolved to a fraction of `a` and small effects are not
# drowned. The response is tested against an analytic answer, on a manufactured
# case with `T = 1 + a P_m(z)` where `div(T grad P_l)` is exact in closed form.
#
# Two arms. The size arm: the response's own relative error must be below
# `RESPONSE_ERROR_MAX`, of the same order as the truncation error rather than
# orders above it. The amplitude arm: the confined assembly is LINEAR in the
# transmissivity, so the response per unit amplitude cannot depend on the
# amplitude, and two amplitudes two orders apart must agree within
# `RESPONSE_AMPLITUDE_INVARIANCE`. That second arm is what GW-24's unconfined
# form does not get, because there `T` depends on the head and the assembly is
# not linear in it.
RESPONSE_ERROR_MAX = 0.30
RESPONSE_AMPLITUDE_INVARIANCE = 1.2
# THE DECISION BAR, and the reason it is this number. The smallest effect this
# component currently argues about is GW-24's unconfined correction at the 2 km
# configuration, which its own sizing puts just under half a per cent. If the
# solved head's mesh error on the mesh actually in use is below that, then the
# mesh is not what stops either GW-18 or GW-24 being measured, and the question
# moves off the discretisation. Above it, the mesh does stop them and the
# instrument statement stands.
HEAD_MESH_FLOOR = 0.005
# The linear solve has to be finished far below the discretisation error it is
# being used to measure, or the measurement is of the solver. The same direct
# factorisation `solve` uses, and this is the relative residual it must reach.
POISSON_SOLVER_RELATIVE = 1e-10
# The resolution sweep, declared with the bars. Four meshes, each a factor of
# four in region count and so a factor of two in cell size, spanning three
# refinement steps. Orogen's own generator at the active build's own jitter, so
# what is being refined is this mesh family and not a tidier one.
INSTRUMENT_SIZES = (20_000, 80_000, 320_000, 1_280_000)
# CLOSURE_TOLERANCE: recharge in equals discharge out. A conservation law over
# float64 sums of a few million terms, so this is a numerical bar, not a
# physical one.
CLOSURE_TOLERANCE = 1e-10
# DIVIDE_AGREEMENT: with uniform permeability and a terrain-following table, the
# groundwater catchments must reproduce the surface catchments. Outside filled
# depressions the two constructions should agree almost everywhere; inside them
# the surface catchment came from a priority flood on a flat filled surface and
# the groundwater trace has nothing to follow, so they need not.
# THE CATCHMENT CHECK, re-specified. Declared before the re-specified test was
# first run.
#
# The first version compared a groundwater trace on the RAW surface against
# `regions.nc:terminal`, which is a priority flood on the FILLED one, and missed
# at 73% because the two describe drainage on different surfaces. No solver
# passes that. Both sides now sit on the filled surface, so the only thing left
# differing is whether the water was routed over the ground or under it.
#
# The bar is on land OUTSIDE the pits the flood filled. Inside them the filled
# surface is flat by construction, so a flux field driven by its gradient
# carries no direction at all and the trace there is arbitrary -- a property of
# flats, not of groundwater. The share of land that is is reported rather than
# assumed small.
DIVIDE_AGREEMENT_UNFILLED = 0.95
# THE TRACE IDENTITY, declared before it was first run, and it is exact.
#
# The terrain-following comparison above is a MEASUREMENT and cannot be an
# identity, because `regions.nc:terminal` comes from a priority flood's
# discovery pointer while a flux trace is a steepest-descent rule, and
# hydrography's README says why the flood cannot use steepest descent: a filled
# pit is flat, so descent would drop whole tributaries. The two routing rules
# genuinely disagree on the same terrain and no solver reconciles them.
#
# What CAN be tested exactly is the trace machinery itself. Hand it a flux field
# whose only outgoing flux at each cell is the face to that cell's own
# `receiver`, and it must reproduce `terminal` on every land region -- not
# almost. That isolates `groundwater_receiver` and `trace_terminals` from the
# physics entirely, and it is the check that caught the receiver reading its own
# sign convention backwards.
TRACE_IDENTITY_EXACT = 1.0
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
# GW-24, DECLARED BEFORE THE FIRST UNCONFINED RUN. With `T = K (h - z_bottom)`
# the matrix depends on the head, so the complementarity problem is no longer
# LINEAR and the argument that licensed `--uniqueness-check` as an identity --
# one symmetric positive definite matrix, therefore exactly one solution --
# lapses. Two trajectories can then differ by whatever head the outer residual
# bar leaves unresolved, so the check becomes a MEASUREMENT against a declared
# relative bar rather than a bit-comparison. 1e-6 is three orders looser than
# the confined identity and still four orders below a millimetre of head, so it
# cannot be met by a solve that has not converged and cannot fail for a reason
# any consumer of the depth field could see.
UNCONFINED_HEAD_RELATIVE = 1e-6
# THE RESIDUAL BAR HAS A ROUND-OFF FLOOR, and under the unconfined form it can
# sit ABOVE `RESIDUAL_TOLERANCE`. The confined model never met it: there the
# matrix is fixed, the direct solve makes the balance zero to the linear
# solver's own accuracy, and the residual reads 3e-13 whatever the
# transmissivity. Under the unconfined form the balance is re-evaluated with the
# conductance the NEW head implies, so it is a difference of face fluxes formed
# from heads whose magnitude is the aquifer thickness while their DIFFERENCE is
# a fraction of a metre. That cancellation puts a floor on the balance of about
#
#     eps * sum_faces Trans (|h_src| + |h_dst|) / total recharge
#
# which is arithmetic about float64 and not a property of any result. The bar
# below is that bound with an ordinary safety factor; on the one-dimensional
# Dupuit case the observed plateau sat at 0.3 and 0.7 of the un-inflated bound
# at 400 m and 2 km of saturated thickness, so the bound is the right shape.
# A configuration whose floor exceeds `CLOSURE_TOLERANCE` cannot be certified at
# all, and `solve` says so rather than reporting a convergence it cannot have.
RESIDUAL_FLOOR_SAFETY = 4.0
# The Dupuit identity, declared before `--dupuit-test` was first run. Part one
# is exact: over a FLAT aquifer base the discrete unconfined face flux equals
# the Kirchhoff-linear flux `g K (u_j - u_i)` by algebra, so any difference is
# round-off. Part two is a discretisation error against the analytic Dupuit
# parabola and is not exact; the bar is what a first-order finite-volume scheme
# on a uniform one-dimensional mesh should comfortably beat.
DUPUIT_IDENTITY_RELATIVE = 1e-12
DUPUIT_ANALYTIC_RELATIVE = 1e-3

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
        """Net INFLOW per cell from a signed per-face flux. Same as
        `geom_divergence`, and the same convention: a positive face flux is flow
        from `dst` into `src`, so this returns what each cell GAINS.

        This docstring said "net outflow" and the arithmetic below has always
        said inflow. Nothing was wrong -- its one consumer is
        `laplace_beltrami_error`, which wants the Laplacian and so wants inflow,
        and that check passes at 0.0006 relative at l=1 -- but the name most
        suggests the other sign, and reading a convention off a docstring rather
        than off the code is exactly what put `groundwater_receiver` backwards
        for the whole of its life. GW-13.
        """
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
# GW-8: what the operator's error is an error IN
# ---------------------------------------------------------------------------

def orogen_mesh(n_regions: int, *, jitter: float, radius_km: float,
                seed: int = 0):
    """A mesh from Orogen's own generator, at a resolution of your choosing.

    `Geometry` needs only generator positions, a radius and a region count, so
    a resolution sweep of THIS mesh family costs no build and no export. The
    construction is `vendor/orogen/js/sphere-mesh.js:generateFibonacciSphere`
    line for line: a Fibonacci spiral, then a per-point jitter in latitude and
    longitude scaled by the local spacing, which is what gives the Voronoi
    cells their irregularity and so is the whole of what is being refined here.

    JITTER IS NOT A DEFAULT. It is a generation parameter of the build, read
    from the build's own manifest by the caller, because a sweep at a tidier
    jitter than the mesh in use would answer a question about a different mesh.

    THE COORDINATE MAP IS THE EXPORT'S, not the generator's. Orogen builds
    `r_xyz` as `(cos lat cos lon, cos lat sin lon, sin lat)` and the export
    declares `lat = asin(y)`, `lon = atan2(x, z)`, so the exported axes are a
    permutation of the generator's: the spiral axis is the export's `y`. The
    Legendre checks are taken about the export's `z`, an equatorial axis, so
    the permutation has to be carried or the sweep tests a different direction
    through the mesh than the planet arm does.

    `cell_area` is set uniform and is UNUSED by everything here: these checks
    take their areas from the Voronoi dual, which is the area the faces bound.

    THE RADIUS IS REQUIRED AND HAS NO DEFAULT. Every error this file reports is
    relative, so the radius cancels out of all of them, which is exactly why a
    literal here would never be caught by a result. It comes from the build's
    own export, and `config/planet.yaml` is the source of that.
    """
    from types import SimpleNamespace

    n = int(n_regions)
    rng = np.random.default_rng(seed)
    s = 3.6 / np.sqrt(n)
    dlong = np.pi * (3.0 - np.sqrt(5.0))
    dz = 2.0 / n
    k = np.arange(n, dtype=np.float64)
    z = 1.0 - dz / 2.0 - k * dz
    lng = k * dlong

    r = np.sqrt(np.maximum(1.0 - z * z, 0.0))
    lat_deg = np.degrees(np.arcsin(z))
    lon_deg = np.degrees(lng)
    if jitter > 0:
        j_lat = rng.random(n) - rng.random(n)
        j_lon = rng.random(n) - rng.random(n)
        next_z = np.maximum(-1.0, z - dz * 2.0 * np.pi * r / s)
        lat_deg = lat_deg + jitter * j_lat * (lat_deg - np.degrees(np.arcsin(next_z)))
        lon_deg = lon_deg + jitter * j_lon * (s / r * 180.0 / np.pi)
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    # Float32 because the export carries float32 coordinates and the face
    # widths this mesh is being used to test are differences of them.
    gx = (np.cos(lat) * np.cos(lon)).astype(np.float32)
    gy = (np.cos(lat) * np.sin(lon)).astype(np.float32)
    gz = np.sin(lat).astype(np.float32)

    return SimpleNamespace(
        n_regions=n, radius_km=float(radius_km),
        x=gy, y=gz, z=gx,
        cell_area=np.full(n, 4.0 * np.pi * radius_km ** 2 / n),
        surface_class=np.full(n, LAND, dtype=np.int64))


def _laplacian_operator(geom: Geometry):
    """`-L` as a matvec, with its diagonal. Positive semidefinite, null space
    the constants, and never formed: the sweep runs to a few million cells and
    the operator is already stored as face weights."""
    n = geom.export.n_regions
    src, dst, g = geom.src, geom.dst, geom.geom
    diag = (np.bincount(src, weights=g, minlength=n)
            + np.bincount(dst, weights=g, minlength=n))

    def matvec(u):
        u = np.asarray(u, dtype=np.float64).ravel()
        return diag * u - _gather(n, src, dst, g, u)

    return matvec, diag


def _gather(n, src, dst, g, u):
    """`sum_j g_ij u_j` per cell. `np.bincount` and not `np.add.at`: this runs
    once per conjugate-gradient iteration on a few million faces, and the
    unbuffered form is what makes the sweep unaffordable rather than slow."""
    return (np.bincount(src, weights=g * u[dst], minlength=n)
            + np.bincount(dst, weights=g * u[src], minlength=n))


def poisson_solution_error(geom: Geometry, degrees=LAPLACE_DEGREES) -> dict:
    """The error in the field the operator SOLVES for, against the same
    analytic answer `laplace_beltrami_error` scores the operator on.

    `laplace_beltrami_error` applies `L` to a known `P_l` and scores the
    residual. This does the opposite, and it is the one a consumer of the water
    table receives: it hands the operator the analytic right-hand side
    `-l(l+1)/R^2 P_l` and asks for the field, then scores that field against
    `P_l`. Nothing about this world enters it either.

    THE NULL SPACE COSTS ONE EQUATION AND NOTHING ELSE. On a closed sphere `L`
    annihilates the constants, so the system is singular and, because the
    right-hand side is made mean-zero, consistent. Deleting one equation and
    fixing that cell's head at zero is then EXACT rather than an
    approximation: `1^T L = 0`, so the deleted equation is the sum of the
    others and holds automatically. What comes back is the true solution up to
    the constant, and both sides are compared with their means removed.

    A DIRECT FACTORISATION, for the reason `solve` gives at length: this is a
    planar-graph Laplacian, minimum degree on `A + A^T` is the right ordering
    for a symmetric pattern, and diagonally preconditioned conjugate gradients
    on it did not converge in a usable number of passes when it was tried here.
    The residual is measured and returned rather than assumed.

    AND THE ORDERING IS NOT ENOUGH ON ITS OWN. SuperLU still pivots for
    stability by default, which reorders rows during the factorisation and
    throws away the fill the symmetric ordering just bought. This matrix is
    symmetric positive definite, so it needs no pivoting at all, and
    `diag_pivot_thresh=0` with `SymmetricMode` says so. Measured on the 80,000
    region mesh of the sweep, whose face weights span five orders of magnitude:
    6.8 million nonzeros in the factor with symmetric mode against a 2.5 GB
    working set and no answer in minutes without it. The setting is not an
    optimisation here, it is what makes the arm runnable at all.
    """
    from scipy.special import eval_legendre

    n = geom.export.n_regions
    z = geom.export.z.astype(np.float64)
    area = geom.flux_area_m2
    src, dst, g = geom.src, geom.dst, geom.geom
    diag = (np.bincount(src, weights=g, minlength=n)
            + np.bincount(dst, weights=g, minlength=n))
    rows = np.concatenate([np.arange(n), src, dst])
    cols = np.concatenate([np.arange(n), dst, src])
    vals = np.concatenate([diag, -g, -g])
    keep = (rows > 0) & (cols > 0)
    a = sp.coo_matrix((vals[keep], (rows[keep] - 1, cols[keep] - 1)),
                      shape=(n - 1, n - 1)).tocsc()

    from scipy.sparse.linalg import splu
    lu = splu(a, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
              options=dict(SymmetricMode=True))

    matvec, _ = _laplacian_operator(geom)
    out = {}
    for l in degrees:
        u = eval_legendre(l, z)
        rhs = area * (l * (l + 1) / geom.radius_m ** 2) * u
        rhs = rhs - rhs.mean()
        h = np.zeros(n)
        h[1:] = lu.solve(rhs[1:])
        resid = float(np.linalg.norm(matvec(h) - rhs) / np.linalg.norm(rhs))
        want = u - u.mean()
        got = h - h.mean()
        out[l] = {
            "relative_rms": float(np.sqrt(np.mean((got - want) ** 2))
                                  / np.sqrt(np.mean(want ** 2))),
            "solver_relative_residual": resid,
            "solver_converged": bool(np.all(np.isfinite(h))
                                     and resid < POISSON_SOLVER_RELATIVE),
        }
    return out


def _legendre_poly(l: int):
    from numpy.polynomial import legendre as _leg, polynomial as _poly
    c = np.zeros(l + 1)
    c[l] = 1.0
    return _poly.Polynomial(_leg.leg2poly(c))


def response_error(geom: Geometry, l: int = 2, m: int = 3,
                   amplitudes=(0.01, 0.2)) -> dict:
    """The operator's error on a RESPONSE to a transmissivity perturbation.

    GW-18 hands the solver a spatially varying multiplier on `T`. What that is
    worth depends on how well this operator resolves a CHANGE rather than a
    field, and that has an analytic answer of its own. For zonal fields on a
    sphere,

        div(T grad u)  =  (1/R^2) d/dz [ T(z) (1 - z^2) du/dz ]

    so with `T = 1 + a P_m(z)` and `u = P_l(z)` the response per unit amplitude
    is `(1/R^2) d/dz [ P_m (1 - z^2) P_l' ]`, a polynomial in `z` evaluated
    exactly. The discrete side assembles face transmissivities through
    `face_transmissivity`, the same arithmetic mean the planet arm uses, and
    subtracts the unperturbed operator applied to the same field.

    Two numbers come back per amplitude and they are the point of the check:
    the error as a fraction of the RESPONSE, and the same error as a fraction
    of the unperturbed FIELD. If the first is of order the truncation error and
    the second is smaller by the amplitude, then the operator's error scales
    with what is being measured and a small effect is not lost in it.
    """
    n = geom.export.n_regions
    z = geom.export.z.astype(np.float64)
    src, dst, gfac = geom.src, geom.dst, geom.geom
    area = geom.flux_area_m2
    no_ocean = np.zeros(n, bool)

    p_l, p_m = _legendre_poly(l), _legendre_poly(m)
    u = p_l(z)
    # q_l = (1 - z^2) P_l', so the response per unit amplitude is
    # (1/R^2) d/dz [P_m q_l], taken exactly in the polynomial basis.
    q_l = np.polynomial.polynomial.Polynomial([1.0, 0.0, -1.0]) * p_l.deriv()
    exact_unit = (p_m.deriv()(z) * q_l(z) + p_m(z) * q_l.deriv()(z)) \
        / geom.radius_m ** 2

    def apply(t_cell):
        t_face = face_transmissivity(t_cell, gfac, src, dst, no_ocean)
        return geom_divergence(n, src, dst, t_face * (u[dst] - u[src])) / area

    base = apply(np.ones(n))
    base_scale = np.sqrt(np.mean(base ** 2))
    exact_scale = np.sqrt(np.mean(exact_unit ** 2))

    out = {}
    for a in amplitudes:
        got = (apply(1.0 + a * p_m(z)) - base) / a
        err = np.sqrt(np.mean((got - exact_unit) ** 2))
        out[float(a)] = {
            "relative_to_response": float(err / exact_scale),
            "relative_to_field": float(a * err / base_scale),
        }
    ratios = [v["relative_to_response"] for v in out.values()]
    out["amplitude_ratio"] = float(max(ratios) / min(ratios))
    return out


def _fitted_order(cell_sizes, errors) -> float:
    """Slope of log(error) against log(cell size). The convergence order."""
    x = np.log(np.asarray(cell_sizes, dtype=np.float64))
    y = np.log(np.asarray(errors, dtype=np.float64))
    return float(np.polyfit(x, y, 1)[0])


def instrument_sweep(sizes=INSTRUMENT_SIZES, *, jitter: float,
                     radius_km: float, degrees=(2, 3, 4), seed: int = 0,
                     verbose: bool = True) -> dict:
    """Truncation and solution error against resolution, on Orogen's generator."""
    rows = []
    for n in sizes:
        export = orogen_mesh(n, jitter=jitter, seed=seed, radius_km=radius_km)
        geom = Geometry(export)
        trunc = laplace_beltrami_error(geom, tuple(degrees))
        soln = poisson_solution_error(geom, tuple(degrees))
        row = {"n_regions": int(n),
               "cell_size_m": float(np.sqrt(4.0 * np.pi / n) * geom.radius_m),
               "closes_on_sphere": geom.closes_on_sphere(),
               "truncation": {int(k): float(v) for k, v in trunc.items()},
               "solution": {int(k): v for k, v in soln.items()}}
        rows.append(row)
        if verbose:
            for l in degrees:
                print(f"    {n:>9,}  l = {l}   truncation "
                      f"{trunc[l]:8.5f}   solution {soln[l]['relative_rms']:10.3e}"
                      f"   {'' if soln[l]['solver_converged'] else 'SOLVER DID NOT CONVERGE'}")
        del geom, export
    h = [r["cell_size_m"] for r in rows]
    orders = {int(l): {"truncation": _fitted_order(h, [r["truncation"][l] for r in rows]),
                       "solution": _fitted_order(
                           h, [r["solution"][l]["relative_rms"] for r in rows])}
              for l in degrees}
    return {"rows": rows, "orders": orders}



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


def transmissivity(k0_m_s, thickness_m):
    """`T = K D`, m2/s. Conductivity over a saturated thickness, CONFINED.

    `thickness_m` is a scalar or a per-region array; the arithmetic is the same
    either way and the choice belongs to the config, not here. GW-18.

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
    `unconfined_transmissivity` below removes the second half of that cost and
    states what it costs in exchange.
    """
    return k0_m_s * np.asarray(thickness_m, dtype=np.float64)


def saturated_thickness(head, aquifer_base_m, min_saturated_m: float):
    """`b = max(h - z_bottom, b_min)`, m. The unconfined saturated column.

    THE FLOOR IS STRUCTURAL, NOT COSMETIC. With `b` free to reach zero a face
    stops conducting when the water table falls to the aquifer base, so which
    cells are connected to a recharge source or to the sea would depend on the
    head -- and the static dry set that `solve` computes once before the
    iteration, which is what made the uniqueness identity pass bit-identically,
    would stop being static. A strictly positive `b_min` keeps every face of the
    conductive network conducting at every head, so that connectivity argument
    survives the change unaltered. The cells sitting on it are reported, and
    their depths are a lower bound rather than a value.
    """
    return np.maximum(head - aquifer_base_m, float(min_saturated_m))


def unconfined_transmissivity(k0_m_s, head, aquifer_base_m,
                              min_saturated_m: float):
    """`T = K (h - z_bottom)`, m2/s. The textbook Dupuit-Forchheimer case.

    WHAT THIS IS FOR. `T = K D` is the CONFINED approximation: transmissivity
    does not fall as the water table drops, so deep dry ground conducts as
    freely as a full aquifer. This form is depth-selective in the direction a
    real unconfined aquifer is -- a shallow table has a thick saturated column
    and drains easily, a deep one has a thin column and drains badly -- which is
    the coupling that lets terrain rather than recharge organise the head.

    AND IT IS NOT FAN'S EXPONENTIAL, which is excluded twice over and must not
    be reached for again. `exp(h/f)` is convex, so a cell mean carries
    `exp(sigma^2/2f^2)` -- `exp(5000)` at 100 m of sub-grid relief against the
    0.95 m `f` Fan's curve reaches on steep bedrock -- and Picard on it
    limit-cycles because `T` moves by a factor of e per e-folding length. This
    form is LINEAR in the head, so a cell mean of `T` is `T` of the cell mean
    exactly and there is no Jensen term at any resolution, and its Picard map
    changes `T` by `|dh| / b` per pass rather than by `exp(|dh| / f)`.

    THE PRICE, stated rather than discovered. `T` now depends on the head, so
    the matrix is no longer fixed, the problem is no longer a LINEAR
    complementarity problem, and the uniqueness argument that licensed
    `--uniqueness-check` as an IDENTITY -- one symmetric positive definite
    matrix, therefore exactly one solution -- no longer applies. Under this form
    that check is a measurement against a declared bar, `UNCONFINED_HEAD_RELATIVE`.
    """
    return k0_m_s * saturated_thickness(head, aquifer_base_m, min_saturated_m)


def kirchhoff_potential(head, aquifer_base_m, min_saturated_m: float):
    """`u = b^2 / 2` for `b` the saturated thickness. Used only by the tests.

    WHY IT IS AN IDENTITY AND NOT A SOLVER. With the aquifer base at the SAME
    elevation on both ends of a face and the face thickness taken as the
    arithmetic mean, the discrete unconfined flux is exactly linear in `u`:

        b_face (h_j - h_i) = (b_i + b_j)(b_j - b_i) / 2 = u_j - u_i

    so the discrete problem in `u` is the same fixed linear system the confined
    model solves, and the Dupuit parabola comes back exactly. That is what
    `--dupuit-test` checks. It is NOT how `solve` runs, because a real aquifer
    base follows the terrain and is not equal across a face, and with a sloping
    base no Kirchhoff transform exists. The identity is the instrument; the
    solver iterates.
    """
    b = saturated_thickness(head, aquifer_base_m, min_saturated_m)
    return 0.5 * b * b


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

def resident_gb() -> float:
    """This process's resident set, GB, read from /proc. Zero where there is none.

    `docs/src/reference/large-data.md` asks for resident memory beside every
    progress line, and gives the reason: it is what turns "it is still going"
    into "it will not finish", and those two have different answers. The solve
    factorises a matrix whose fill is the thing most likely to exhaust the host,
    so the number belongs on the pass line rather than in a postmortem.
    """
    try:
        with open("/proc/self/status", encoding="ascii") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    return float(line.split()[1]) / (1 << 20)
    except OSError:
        pass
    return 0.0


class _Phases:
    """Wall time per named phase of the outer loop, and the pass it belongs to.

    A solve of this size has three plausible costs -- the scatter-adds over
    thirty million faces, the sparse factorisation, and the connectivity passes
    -- and which of them dominates decides what is worth changing. Measuring
    them separately is the difference between the row this instrumentation was
    added for and a guess: the file's own recorded reasoning was sized at
    400,000 unknowns on a smaller build, and a cost that grows superlinearly in
    the unknown count does not stay where it was measured.

    The clock is `perf_counter`, so this measures WALL time and is therefore
    only as meaningful as the load it was taken under. Record the load.
    """

    __slots__ = ("t", "_open", "_t0")

    def __init__(self):
        self.t = {}
        self._open = None
        self._t0 = 0.0

    def __call__(self, name):
        self.stop()
        self._open = name
        self._t0 = time.perf_counter()
        return self

    def stop(self):
        if self._open is not None:
            self.t[self._open] = self.t.get(self._open, 0.0) + (
                time.perf_counter() - self._t0)
            self._open = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.stop()
        return False

    def report(self) -> dict:
        self.stop()
        return {k: round(v, 3) for k, v in
                sorted(self.t.items(), key=lambda kv: -kv[1])}


def et_rate(et_max_m_s, et_lambda_m, surface_m, head, conductive):
    """Groundwater ET as a rate, exponential in depth below the surface.

    Shah, Nachabe and Ross (2007): "The decline of ET with DTWT is better
    simulated by an exponential decay function than the commonly used linear
    decay." Depth is clamped at zero so a cell at the surface takes the full
    rate rather than an extrapolated one.
    """
    depth = np.maximum(surface_m - head, 0.0)
    return np.where(conductive, et_max_m_s * np.exp(-depth / et_lambda_m), 0.0)


def et_balance_depth(et_max_m_s, et_lambda_m, supply_m3_s, area_m2, conductive):
    """Depth at which the sink alone consumes a cell's own recharge.

        R A = ET_max A exp(-d / lambda)   =>   d = lambda ln(ET_max A / R A)

    One equation in one unknown, so it needs no solver, and it is where a cell
    with no lateral exchange comes to rest. Zero where the sink cannot match the
    supply even at the surface, which is the case that seeps instead.
    """
    cap = et_max_m_s * area_m2
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.where((cap > supply_m3_s) & conductive,
                     et_lambda_m * np.log(cap / np.maximum(supply_m3_s, 1e-300)),
                     0.0)
    return np.clip(np.nan_to_num(d, nan=0.0, posinf=0.0), 0.0, None)


def solve(export: Export, geom: Geometry, *, k0_m_s, thickness_m, recharge_m_s,
          surface_m, conductive, sea_level_m=0.0, max_outer=200,
          start_all_free=False, verbose=True,
          et_max_m_s=None, et_lambda_m=None, fixed_head_m=None,
          aquifer_base_m=None, min_saturated_m=1.0):
    """Steady-state head, seepage and face fluxes.

    `aquifer_base_m` selects the UNCONFINED form, GW-24. Left None the model is
    the confined `T = K D` this component has always run and every array below
    is assembled once; given a per-region aquifer base elevation the
    transmissivity becomes `K (h - z_bottom)` and is reassembled from the head
    on every pass. `thickness_m` is then unused for flow and only names the
    geometry the base came from. The two forms share every other line here on
    purpose: a second solver is what this component already has two dead ones of.

    With `T = K D` the transmissivity does not depend on the head, so the matrix
    is FIXED and the only thing left to iterate is the active set:

        sum_j Trans_ij (h_j - h_i) + R_i A_i - E_i(h_i) A_i = S_i,
                                                     h_i <= z_i,  S_i >= 0

    `E` is groundwater evapotranspiration and is OFF unless `et_max_m_s` is
    given, in which case the result is bit-identical to a run without it because
    no ET code path executes at all. GW-15. Shah et al. (2007) find the decline
    with water table depth is better described by an exponential than by the
    linear form MODFLOW's EVT uses, so

        E(h) = et_max * exp(-(z - h) / et_lambda)

    It is head-dependent, so the matrix is no longer fixed and this becomes a
    Newton iteration folded into the active-set loop already here. That is
    benign in a way GW-9's exponential transmissivity was not: dE/dh is POSITIVE,
    so it lands on the diagonal with the sign that strengthens diagonal
    dominance, and the sink weakens as the table falls rather than the
    conductance collapsing as it did before.

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
    phases = _Phases()
    _t_start = time.perf_counter()
    phases("setup")
    n = export.n_regions
    src, dst = geom.src, geom.dst
    gfac = geom.geom
    area_m2 = geom.volume_area_m2

    is_ocean = export.surface_class != LAND
    # LOCAL BASELEVELS. GW-17. `fixed_head_m` is finite on cells whose head is
    # imposed rather than solved -- a river or a lake, whose surface IS the water
    # table there. Without them the ocean is the only fixed head in the problem,
    # so recharge in a continental interior has to reach a coast a thousand
    # kilometres off, no defensible transmissivity carries it, and the local sink
    # balance sets the depth everywhere at a few times lambda whatever the
    # terrain does. That is the range collapse GW-3 measured, and a baselevel
    # 15 km away rather than 1,000 km is the difference between needing 2 km of
    # aquifer and needing 600.
    #
    # They behave exactly as ocean cells do: a fixed head, not a conducting cell,
    # and outside the set being solved for. `is_ocean` stays a separate name
    # because the report attributes coastal discharge with it.
    imposed = (np.isfinite(fixed_head_m) & ~is_ocean
               if fixed_head_m is not None else np.zeros(n, bool))
    is_boundary = is_ocean | imposed
    conductive = conductive & ~imposed
    in_domain = conductive | is_boundary
    live = in_domain[src] & in_domain[dst]
    src, dst, gfac = src[live], dst[live], gfac[live]

    # GW-24. Under the unconfined form the face thickness is the ARITHMETIC mean
    # of the two cells' saturated columns and the face conductivity keeps the
    # existing rule, land-side at a boundary. That pairing is not a free choice:
    # over a flat aquifer base the arithmetic mean is what makes the discrete
    # flux exactly `g K (u_j - u_i)` in the Kirchhoff potential `u = b^2/2`, so
    # the scheme reduces to the linear one it replaces rather than to something
    # near it, and `--dupuit-test` checks that reduction as an identity.
    #
    # A cell with no aquifer base -- ocean, and anything else the base field
    # leaves non-finite -- contributes no thickness, so such a face takes the
    # other side's column, exactly as it takes the other side's conductivity.
    unconfined = aquifer_base_m is not None
    has_base = np.isfinite(aquifer_base_m) if unconfined else None

    def assemble(head_now):
        """Per-cell and per-face transmissivity at the head handed in."""
        if not unconfined:
            t = np.where(conductive, transmissivity(k0_m_s, thickness_m), 0.0)
            return t, face_transmissivity(t, gfac, src, dst, is_boundary)
        b = np.where(has_base,
                     saturated_thickness(head_now,
                                         np.where(has_base, aquifer_base_m, 0.0),
                                         min_saturated_m), 0.0)
        t = np.where(conductive, k0_m_s * b, 0.0)
        k_face = face_transmissivity(np.where(conductive, k0_m_s, 0.0),
                                     gfac, src, dst, is_boundary)
        b_face = np.where(has_base[src] & has_base[dst],
                          0.5 * (b[src] + b[dst]),
                          np.where(has_base[src], b[src], b[dst]))
        return t, k_face * b_face

    def row_sums(trans_now):
        rs = np.zeros(n)
        np.add.at(rs, src, trans_now)
        np.add.at(rs, dst, trans_now)
        return rs

    t_cell, trans = assemble(surface_m)
    rowsum = row_sums(trans)

    supply = recharge_m_s * area_m2               # m3/s per cell

    # DRY GROUND, DECIDED BEFORE THE ITERATION STARTS AND NOT DURING IT.
    #
    # Water enters this system in exactly two places: recharge falling on a
    # conductive cell, and the ocean boundary. A conductive cell that connects
    # to neither, through any chain of conducting faces, can never hold
    # groundwater at all. It is dry as a matter of terrain, lithology and
    # recharge, and that is a static property of the graph.
    #
    # It was previously discovered mid-solve instead, as blocks that happened to
    # be surrounded by free cells at some pass. Whether a block LOOKS isolated
    # depends on the active set, so which cells were found dry depended on the
    # route taken: the forward run marked 14 and the reverse run 80. The two
    # trajectories were therefore solving slightly different problems, which is
    # what the uniqueness identity was detecting when it missed by 12.76 m.
    #
    # A face conducts whenever both its ends are in the network, whatever their
    # pinned or free status, so this connectivity is the same on every pass and
    # the set it produces cannot depend on traversal order.
    both_cond = conductive[src] & conductive[dst]
    comp_graph = sp.coo_matrix(
        (np.ones(int(both_cond.sum())),
         (src[both_cond], dst[both_cond])), shape=(n, n))
    ncomp_static, comp = _components(comp_graph)
    fed = np.zeros(ncomp_static, bool)
    has_recharge = conductive & (supply > 0)
    fed[comp[has_recharge]] = True
    touches_sea = np.zeros(n, bool)
    coastal = is_boundary[src] ^ is_boundary[dst]
    if coastal.any():
        touches_sea[np.where(is_boundary[src[coastal]],
                             dst[coastal], src[coastal])] = True
    fed[comp[conductive & touches_sea]] = True
    dry = conductive & ~fed[comp]
    if dry.any():
        conductive = conductive & ~dry
        t_cell, trans = assemble(surface_m)
        rowsum = row_sums(trans)
        if verbose:
            print(f"  {int(dry.sum()):,} conductive regions reach neither "
                  f"recharge nor the sea; dry before the solve starts")

    total_supply = float(supply[conductive].sum())

    # Start with every cell at the surface, which is the zero-permeability
    # solution and is always feasible, then release. Growing the free set from
    # nothing keeps it a fraction of the mesh; starting all-free asks the solver
    # to factor the whole 2.5 million cells before a single cell is pinned.
    head = surface_m.copy()
    head[is_ocean] = sea_level_m
    if fixed_head_m is not None:
        head[imposed] = fixed_head_m[imposed]
    if et_max_m_s is not None:
        # START AT THE LOCAL SINK BALANCE, not at the surface. Every land cell
        # used to begin where the sink runs at its maximum, so each had to walk
        # down lambda*ln(ET_max/R) -- five to seven lambda for this world's
        # recharge -- with the active set churning the whole way. On Vesper that
        # did not converge in 400 passes.
        #
        # The local balance is a closed form and is where an isolated cell comes
        # to rest, so lateral exchange is left as the only thing to solve for
        # and is a perturbation rather than the whole answer. GW-12 proved this
        # solution unique, and that uniqueness is what licenses moving where the
        # iteration starts: it must land in the same place, and the identity
        # says so rather than my say-so.
        head -= et_balance_depth(et_max_m_s, et_lambda_m, supply,
                                 area_m2, conductive)
        head[is_ocean] = sea_level_m
        if fixed_head_m is not None:
            head[imposed] = fixed_head_m[imposed]
    # `start_all_free` walks the opposite trajectory, for the uniqueness check.
    # The solution of this complementarity problem is unique, so where it lands
    # cannot depend on where it started.
    # PINNED MEANS AT THE SURFACE, and the initial free set has to say so.
    # Seeding the head below the surface while leaving every cell pinned breaks
    # that invariant: `pinned` is read as "at the surface" by the seepage
    # accounting and by the report, and a cell sitting below it is neither
    # seeping nor free. The two trajectories then converge, both pass closure,
    # and land 2.8 m apart -- which the uniqueness identity caught at a bar of
    # 1e-9, and which is the whole reason that identity exists.
    free = conductive.copy() if start_all_free else (conductive & (head < surface_m))
    anchor = np.zeros(n, bool)
    anchor_failed = np.zeros(n, bool)
    result, trace = {}, []
    # Defined before the loop because the anchoring branch `continue`s past the
    # residual computation. A run that anchors on every pass used to fall out of
    # the loop and raise UnboundLocalError on these, which reports a Python
    # problem where the real one is that the solve never converged.
    residual = float("inf")
    leak = leak_frac = float("inf")
    infeasible = -1
    bar = RESIDUAL_TOLERANCE

    # THE MATRIX FINGERPRINT, per pass, and it is the measurement world-wfge
    # asks for rather than a diagnostic. The one lever the call-site comment
    # below names -- `splu` re-solves at a fraction of the cost of factorising
    # -- is reachable only on a pass whose matrix REPEATS, so how often that
    # happens is the whole question and it must be counted rather than assumed.
    # Two levels, because they have different remedies: the PATTERN repeats when
    # the free set is unchanged, and the VALUES repeat when the numbers in it
    # are unchanged too. A confined run's off-diagonals are fixed, so under the
    # confined form the second differs from the first only through the ET
    # diagonal.
    prev_free_key = None
    prev_value_key = None
    pattern_repeats = value_repeats = 0
    factors = {}
    n_factorised = n_reused = 0
    n_blocks_last = largest_block_last = 0
    for outer in range(max_outer):
        phases("active_set")
        # A cell with no conducting face has no LATERAL equation: nothing can
        # carry its recharge away sideways, so without a sink it stands at the
        # surface and seeps.
        stranded = free & ~is_boundary & (rowsum <= 0)
        free[stranded] = False
        head[stranded] = surface_m[stranded]
        if et_max_m_s is not None and stranded.any():
            # WITH a sink it has a LOCAL one, and pinning it at the surface
            # instead makes it permanently infeasible: at the surface the sink
            # runs at full rate, and where that exceeds the cell's own recharge
            # the balance is negative and no later pass can fix it, because the
            # cell has no face to release across. That is a closure failure of
            # 1.7e-4 rather than a wrong number in a corner.
            #
            # Its balance is one equation in one unknown, so solve it:
            #     R_i A_i = ET_max_i A_i exp(-(z_i - h_i) / lambda)
            # which the water table reaches by falling until the sink matches
            # the supply. Where the sink cannot match it even at the surface,
            # the surface is right and the excess seeps.
            drop = et_balance_depth(et_max_m_s, et_lambda_m, supply,
                                    area_m2, conductive)[stranded]
            head[stranded] = surface_m[stranded] - drop

        # GW-24. Under the unconfined form the conductance belongs to the head
        # it was computed at, so it is refreshed wherever the head has moved --
        # here, after the stranded cells were placed, and again below after the
        # linear solve. That is a Picard step on the saturated thickness and it
        # is the whole of the nonlinear iteration: it changes `T` by `|dh| / b`
        # per pass, where the exponential form GW-9 abandoned changed it by
        # `exp(|dh| / f)`, which is why that one limit-cycled and this one is
        # a contraction wherever the head step is smaller than the saturated
        # column. A confined run reassembles nothing and is bit-identical.
        if unconfined:
            with phases("assemble"):
                t_cell, trans = assemble(head)
                rowsum = row_sums(trans)

        phases("release")
        # RELEASE before solving. A pinned cell whose neighbours already draw
        # more water out of it than its recharge supplies cannot stand at the
        # surface. An anchor is never released; see below.
        flux = trans * (head[dst] - head[src])
        et_m3_s = (et_rate(et_max_m_s, et_lambda_m, surface_m, head, conductive)
                   * area_m2) if et_max_m_s is not None else 0.0
        seep = supply - et_m3_s + geom_divergence(n, src, dst, flux)
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

        phases("anchor")
        unknown = free & ~is_boundary
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
        if et_max_m_s is not None:
            # A CELL WITH AN ET SINK ANCHORS ITSELF. The anchoring below exists
            # because a block touching no pinned cell is a pure Neumann problem
            # and singular. GW-15's sink puts a positive term on the diagonal,
            # so such a block is no longer Neumann and no longer singular: it
            # disposes of its own recharge by evaporating it, which is what the
            # water does. Anchoring it anyway pins a cell the release test then
            # frees, and the two cycle until max_outer with no residual ever
            # computed.
            anchored |= unknown & (et_rate(et_max_m_s, et_lambda_m,
                                           surface_m, head, conductive) > 0.0)
        if not np.all(anchored[unknown]):
            sub = sp.coo_matrix((np.ones(int(both.sum())),
                                 (idx[src[both]], idx[dst[both]])), shape=(m, m))
            ncomp, label = _components(sub, directed=False)
            ok = np.zeros(ncomp, bool)
            ok[label[idx[anchored & unknown]]] = True
            orphan = ~ok[label]
            if orphan.any():
                # WHICH cell to anchor, and it is not a free choice. A block
                # with no external face must dispose of all its own recharge
                # internally, so its cells' balances sum to something positive
                # and at least one of them genuinely seeps. Anchoring THAT one
                # is the only choice that stays feasible; anchoring the cell
                # with the shallowest water table, which is what this did
                # first, picked a cell that immediately wanted releasing and
                # could not be, so four of them sat infeasible for three
                # hundred passes and the solve never converged.
                #
                # AND THAT PREMISE IS FALSE WHEN THE SINK IS ON. A block can
                # dispose of its recharge by evaporating it, so its balances
                # need not sum positive and none of its cells need seep. It is
                # also no longer singular, which is why a cell with a live sink
                # anchors itself above and never reaches here. The ET term is
                # subtracted anyway, because the ranking below is meaningless
                # against a balance that ignores a door the water leaves by.
                where = np.flatnonzero(unknown)[orphan]
                et_now = (et_rate(et_max_m_s, et_lambda_m, surface_m, head,
                                  conductive) * area_m2
                          if et_max_m_s is not None else 0.0)
                bal_now = (supply - et_now + geom_divergence(
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

        # THE FREE SET IS NOT ONE PROBLEM, AND SOLVING IT AS ONE IS THE COST.
        # GW-17 put every river and lake cell at a fixed head, which makes them
        # BOUNDARY rather than unknown, so a face touching one carries no
        # off-diagonal and the channel network CUTS the free set. What is left
        # is one block per interfluve, plus one per island, and the blocks do
        # not exchange water directly -- they exchange it through the rivers
        # between them, which are a Dirichlet condition on both sides.
        #
        # A direct factorisation's work and fill are SUPERLINEAR in the
        # unknowns, so a partition is not bookkeeping: K blocks of n/K cells
        # cost `K (n/K)^1.5`, which is `n^1.5 / sqrt(K)`, and the PEAK is one
        # block rather than the whole. That is what makes the difference between
        # a solve that fits in this host's memory and one that swaps.
        #
        # AND IT IS EXACT, not an approximation to the coupled solve. There is
        # no fill between blocks because there are no entries between them, so
        # the elimination inside a block is the same sequence of operations it
        # would be inside the whole matrix. `--factorisation-test` asserts that
        # equality bitwise on a scrambled multi-block case, which is the arm
        # that could show it failing.
        #
        # THE REUSE LIVES HERE TOO, and it is why this row's 27x lever is
        # reachable at all. The whole matrix repeats almost never, because the
        # active set churns SOMEWHERE on a planet every pass and the
        # evapotranspiration diagonal moves with the head. A single interfluve
        # is a different matter: most blocks are untouched between passes, so
        # most blocks skip the factorisation and pay the back-substitution
        # alone.
        phases("partition")
        blk_graph = sp.coo_matrix(
            (np.ones(int(both.sum()), dtype=np.int8),
             (idx[src[both]], idx[dst[both]])), shape=(m, m))
        n_blocks, blk = _components(blk_graph, directed=False)
        del blk_graph
        # STABLE, so the cells inside a block keep their relative order and the
        # block's submatrix is the one the whole matrix would have held. That is
        # the condition the bitwise equality rests on.
        order = np.argsort(blk, kind="stable")
        bounds = np.searchsorted(blk[order], np.arange(n_blocks + 1))
        n_blocks_last = n_blocks
        largest_block_last = int(np.diff(bounds).max()) if n_blocks else 0
        pos = np.empty(m, dtype=np.int64)
        pos[order] = np.arange(m)
        unknown_ids = np.flatnonzero(unknown)
        free_key = hash(unknown_ids.tobytes())
        # `idx[unknown]` was `arange(m)`, so this relabels every free cell by
        # its position in block order and leaves -1 everywhere else.
        idx[unknown] = pos

        phases("assembly")
        rows = np.concatenate([idx[src[both]], idx[dst[both]],
                               idx[src[both]], idx[dst[both]]])
        cols = np.concatenate([idx[dst[both]], idx[src[both]],
                               idx[src[both]], idx[dst[both]]])
        vals = np.concatenate([-trans[both], -trans[both],
                               trans[both], trans[both]])
        rhs = np.zeros(m)
        np.add.at(rhs, idx[unknown], supply[unknown])
        if et_max_m_s is not None:
            # Newton on E(h) about the current head: E(h) ~ E* + (E*/lambda)(h - h*).
            # The slope goes on the diagonal, where it is POSITIVE and so only
            # improves the M-matrix; the rest moves to the right-hand side.
            e_star = et_rate(et_max_m_s, et_lambda_m, surface_m, head, conductive)
            slope = e_star[unknown] * area_m2[unknown] / et_lambda_m
            rows = np.concatenate([rows, idx[unknown]])
            cols = np.concatenate([cols, idx[unknown]])
            vals = np.concatenate([vals, slope])
            # `rhs` is in BLOCK order and `slope` is in free-cell order, so
            # this scatters rather than adds elementwise. Adding them directly
            # puts every cell's sink term on a different cell.
            rhs[pos] += (slope * head[unknown]
                         - e_star[unknown] * area_m2[unknown])
        if edge_one.any():
            u_side = np.where(unknown[src[edge_one]], src[edge_one], dst[edge_one])
            p_side = np.where(unknown[src[edge_one]], dst[edge_one], src[edge_one])
            rows = np.concatenate([rows, idx[u_side]])
            cols = np.concatenate([cols, idx[u_side]])
            vals = np.concatenate([vals, trans[edge_one]])
            np.add.at(rhs, idx[u_side], trans[edge_one] * head[p_side])

        # ONE CONVERSION, NOT TWO. The operator is symmetric -- the
        # off-diagonals are `-trans` on both sides of every face, and every
        # other contribution lands on the diagonal -- so its CSR arrays ARE its
        # CSC arrays and the transpose is free. This used to build the CSR, read
        # the diagonal off it, and then convert the whole thing again for a
        # solver that wants CSC; the second conversion sorts and rebuilds arrays
        # that already held the right numbers in the right order.
        _csr = sp.coo_matrix((vals, (rows, cols)), shape=(m, m)).tocsr()
        A = sp.csc_matrix((_csr.data, _csr.indices, _csr.indptr), shape=(m, m))
        del _csr
        # Cheap, exact, and order-independent: the free set IS the pattern under
        # a fixed face list, and the assembled values are the matrix. Hashing
        # the bytes costs one pass over arrays the assembly above just built.
        value_key = hash((free_key, A.data.tobytes()))
        pattern_repeats += int(free_key == prev_free_key)
        value_repeats += int(value_key == prev_value_key)
        repeated = value_key == prev_value_key
        prev_free_key, prev_value_key = free_key, value_key
        if not np.all(A.diagonal() > 0):
            raise SystemExit(
                f"{int((A.diagonal() <= 0).sum())} free cells are isolated; "
                "the operator is singular and the solve is invalid")
        # A direct factorisation, not conjugate gradients. Conductivity still
        # spans the 3.4 orders of magnitude between Gleeson's classes, which
        # diagonally preconditioned CG handles badly, and the free set is small
        # enough that a direct solve is comfortable.
        #
        # MEASURED 2026-08-20 against three alternatives, on a planar-graph
        # Laplacian of this problem's size and coefficient span -- 399,424
        # unknowns, 1,994,592 nonzeros, conductivity drawn across Gleeson's
        # -15.2 to -11.8. SuperLU is the fastest of the four:
        #
        #     SuperLU spsolve, this line      1.39 s
        #     UMFPACK via scikit-umfpack      1.93 s   0.72x
        #     pyamg smoothed aggregation + CG 7.05 s   0.20x
        #     CG, diagonal preconditioner      20 s    DID NOT CONVERGE
        #
        # So the comment above survives testing rather than merely sounding
        # right, and neither library is worth a dependency. A sphere mesh is a
        # planar graph, where nested dissection gives a direct solve near-optimal
        # fill; algebraic multigrid earns its setup cost in 3D or at far larger
        # sizes, and neither applies here.
        #
        # THAT COMPARISON IS SIZED AT 400,000 UNKNOWNS AND THE ACTIVE BUILD IS
        # NOT. It carried the conclusion that the linear solve is about 30 s
        # inside a multi-minute run, which was true of the build it was measured
        # on. A direct factorisation's cost is superlinear in the unknowns, so
        # that conclusion does not travel: on canonical-10m-carve2 the free set
        # opens at millions of cells rather than hundreds of thousands, and one
        # confined solve ran past forty minutes at about 16 GB and was killed
        # without writing a water table. The solve IS the cost there, and what
        # the partition below is for.
        #
        # THE 27x REUSE LEVER, MEASURED RATHER THAN PROPOSED, AND IT IS A LOSS.
        # The lever named here was to hold the evapotranspiration diagonal fixed
        # across inner steps so the matrix repeats and `splu` can re-solve at a
        # fraction of the cost of factorising. On the uniqueness case, which
        # carries the sink and the imposed baselevels:
        #
        #                          passes   factorisations   reused   residual
        #     Newton, as run           10               30        0   2.25e-13
        #     diagonal held fixed     246               15      723   9.58e-13
        #
        # The reuse arrives -- the matrix repeats on 241 of 246 passes, which is
        # the 98 per cent the lever promised -- and the iteration still loses,
        # because it needs twenty-five times the passes and every pass pays an
        # assembly, a partition and two water balances over the whole mesh
        # whether or not it factorises. The head also lands 2.7e-07 m away
        # rather than bit-identical.
        #
        # AND THE REASON IS THE ONE THIS FILE ALREADY RECORDS ELSEWHERE. The
        # sink is exponential in depth at an e-folding of `et_lambda_m`, about a
        # metre, so a linearisation of `E` about the current head is good only
        # within about a metre of it, while a Newton step on this problem moves
        # the head by tens of metres. Freezing the diagonal is a fixed-point
        # iteration on an exponential outside its own e-folding length, which is
        # exactly why Fan's exponential transmissivity limit-cycled. One
        # exponential, two solvers, the same disqualification.
        # ORDERING, AND IT IS A NO-OP AT THIS SIZE. scipy's default for spsolve
        # is COLAMD, which orders for an unsymmetric pattern; this operator's
        # pattern IS symmetric, a graph Laplacian with Dirichlet rows, so
        # minimum degree on A + A^T is nominally the right ordering. MEASURED
        # with `splu` on a planar-graph Laplacian of this problem's shape, the
        # fill halves and the factor is about 3x faster at 2M unknowns:
        #
        #        n        COLAMD              MMD_AT_PLUS_A
        #      399,424     50.3 M   1.5 s      27.5 M   1.1 s
        #    1,000,000    145.1 M   6.2 s      78.5 M   3.6 s
        #    1,999,396    335.6 M  23.4 s     168.0 M   8.0 s
        #
        # THAT BENEFIT DOES NOT SURVIVE ON THIS OPERATOR AT THIS SIZE, and the
        # negative result is recorded here because the table above would
        # otherwise read as a promise. On a synthetic of `solve()`'s own shape
        # at 55,977 free cells the fill is IDENTICAL either way: SuperLU's
        # default partial pivoting reorders for stability and undoes the
        # symmetric ordering unless `SymmetricMode` is on with
        # `diag_pivot_thresh=0`, which `spsolve` cannot pass and `splu` can --
        # see the analytic check in this file, which does exactly that. The
        # keyword is kept because it is free and correct, not because it is
        # known to buy anything on the forward solve.
        #
        # WHERE IT WOULD MATTER IS THE REVERSE UNIQUENESS TRAJECTORY, which
        # starts with every land cell free: 834 M entries against 400 M at this
        # build's 4,328,736 land cells, about 10 GB against 4.8 GB before
        # SuperLU's working set. That arm needs `splu` with `SymmetricMode`
        # rather than this call, and it is the one lever WORLD-V3UR's memory
        # blocker has not been tried against.
        # `splu(A, permc_spec=...).solve(b)` IS BIT-IDENTICAL to
        # `spsolve(A, b, permc_spec=...)`: both drive the same SuperLU
        # factorisation and the same back-substitution. That equality is
        # ASSERTED by `--factorisation-test` rather than assumed, because it is
        # what licenses factorising a block on its own at all. What `splu` adds
        # is that the factor survives the call.
        #
        # `SymmetricMode` is deliberately NOT set here. It halves the fill and
        # it is measurably NOT bit-identical -- 2e-14 to 6e-14 relative on a
        # Laplacian of this operator's coefficient span -- which forfeits the one
        # thing that makes a restructuring of this solve checkable. That lever
        # belongs to world-v3ur's memory blocker, where an arm that changes the
        # answer at round-off can be declared as such.
        #
        # ONE FACTORISATION PER BLOCK, KEPT ONLY WHILE THAT BLOCK IS UNCHANGED.
        # The key is the block's own cells and its own numbers, so a block the
        # active set did not touch and whose sink did not move is recognised
        # whatever happened elsewhere on the planet. `factors` is rebuilt each
        # pass from what was actually used, which is what evicts the blocks that
        # went away: a dict that only grows would hold the fill of every free
        # set the iteration ever passed through.
        x = np.empty(m)
        kept = {}
        pass_factorised = n_factorised
        with phases("factor_solve"):
            # Sliced off the CSC arrays directly rather than through scipy's
            # indexing, because the block count runs to the number of
            # interfluves on a planet and a per-block object built the general
            # way costs more than the factorisation of a small one. The columns
            # of a block are contiguous in this ordering and every row they
            # touch is inside it, so the slice is exact.
            indptr, indices, data = A.indptr, A.indices, A.data
            for b0, b1 in zip(bounds[:-1], bounds[1:]):
                if b1 <= b0:
                    continue
                p0, p1 = indptr[b0], indptr[b1]
                sub = sp.csc_matrix(
                    (data[p0:p1], indices[p0:p1] - b0, indptr[b0:b1 + 1] - p0),
                    shape=(b1 - b0, b1 - b0))
                key = (hash(unknown_ids[order[b0:b1]].tobytes()),
                       hash(sub.data.tobytes()))
                # POPPED, not read: a block that survives MOVES from the old
                # table to the new one instead of being held in both, so the
                # transient while the pass rebuilds its factors is the blocks
                # that CHANGED and not the whole free set twice over.
                lu = factors.pop(key, None)
                if lu is None:
                    lu = splu(sub, permc_spec="MMD_AT_PLUS_A")
                    n_factorised += 1
                else:
                    n_reused += 1
                kept[key] = lu
                x[b0:b1] = lu.solve(rhs[b0:b1])
        # The previous pass's factors are dropped here, not before the loop, so
        # a block that survived is never rebuilt to be immediately discarded.
        factors = kept
        del kept, lu
        pass_factorised = n_factorised - pass_factorised
        if not np.all(np.isfinite(x)):
            raise SystemExit(
                f"the direct solve returned {int((~np.isfinite(x)).sum()):,} "
                f"non-finite heads at pass {outer}; do not use this result")
        head[unknown_ids[order]] = x

        # No relaxation and no step limit: under the confined form the operator
        # does not depend on the answer, so a full step invalidates nothing.
        over = unknown & (head > surface_m)
        head[over] = surface_m[over]
        free[over] = False

        # THE RESIDUAL IS MEASURED AGAINST THE CONDUCTANCE THE NEW HEAD IMPLIES.
        # Under the unconfined form the linear system just solved used the
        # PREVIOUS head's saturated thickness, so a balance evaluated with that
        # same `trans` would close to round-off on every pass and report a
        # converged solve the moment the linear solver worked. Reassembling here
        # makes the residual the true nonlinear one: it is small only when the
        # head has stopped moving the thickness that produced it. A confined run
        # reassembles nothing and the residual keeps its old meaning exactly.
        if unconfined:
            with phases("assemble"):
                t_cell, trans = assemble(head)
                rowsum = row_sums(trans)

        phases("residual")
        flux = trans * (head[dst] - head[src])
        et_m3_s = (et_rate(et_max_m_s, et_lambda_m, surface_m, head, conductive)
                   * area_m2) if et_max_m_s is not None else 0.0
        bal = supply - et_m3_s + geom_divergence(n, src, dst, flux)
        still = free & ~is_boundary
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

        # The round-off floor of the balance, computed rather than assumed. Zero
        # for the confined form, whose fixed matrix makes the balance exact.
        bar = RESIDUAL_TOLERANCE
        if unconfined:
            floor = (np.finfo(float).eps
                     * float((trans * (np.abs(head[src]) + np.abs(head[dst]))).sum())
                     / max(total_supply, 1e-30))
            bar = max(RESIDUAL_TOLERANCE, RESIDUAL_FLOOR_SAFETY * floor)
            result["residual_floor"] = float(floor)
            result["residual_bar"] = float(bar)

        flips = int(over.sum()) + int(release.sum())
        trace.append([outer, m, flips, residual, infeasible, leak_frac])
        if verbose:
            # FLUSHED, and that is not decoration. A solve of this size has been
            # killed under contention with its stdout block-buffered to a file,
            # which lost every pass line and left a zero-byte log: there was not
            # even a trace to say how far it got. `docs/src/reference/large-data.md`
            # asks for progress, elapsed and resident memory per unit of work,
            # and a solve whose progress is knowable only if it exits normally
            # cannot be diagnosed when it does not.
            print(f"  pass {outer:3d}  free {m:>9,}  released "
                  f"{int(release.sum()):>7,}  pinned {int(over.sum()):>7,}"
                  f"  residual {residual:10.3e}  leak {leak_frac:10.3e}"
                  f"  ({infeasible:,} cells)  {time.perf_counter() - _t_start:6.1f}s"
                  f"  {resident_gb():5.2f} GB"
                  f"  blocks {n_blocks_last:,}/{largest_block_last:,}"
                  f"  refactored {pass_factorised:,}"
                  f"{'  matrix repeated' if repeated else ''}", flush=True)
        if (residual < bar and flips <= FLIP_TOLERANCE * m
                and leak_frac < bar):
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
                  f"{bar:.3e}; {infeasible:,} pinned cells "
                  f"infeasible")
            bad = np.flatnonzero(pinned_now & (bal < 0))
            for i in bad[:8]:
                # The ET is printed because the balance contains it: a report
                # showing only supply against a balance that subtracts the sink
                # reads as arithmetic that does not add up.
                e_i = (et_m3_s[i] if et_max_m_s is not None else 0.0)
                print(f"    region {i}: supply {supply[i]:.6e} m3/s  "
                      f"et {e_i:.6e}  "
                      f"balance {bal[i]:.6e}  anchor {bool(anchor[i])}  "
                      f"failed_anchor {bool(anchor_failed[i])}  "
                      f"rowsum {rowsum[i]:.3e}  depth {surface_m[i]-head[i]:.2f} m")
    factors = {}
    phases.stop()
    result["residual_trace"] = trace
    # ON THE ARTIFACT, not in a log line, because the row this was added for
    # asks for the share of passes that reuse a factorisation to be STATED
    # rather than assumed, and a number that lives only in a terminal cannot be
    # cited later. `passes` is the count of passes that assembled a matrix, so
    # the shares below are against the number of chances there were.
    result["phase_seconds"] = phases.report()
    result["matrix_repeat"] = {
        "passes_assembled": int(len(trace)),
        "pattern_repeats": int(pattern_repeats),
        "value_repeats": int(value_repeats),
        "blocks": int(n_blocks_last),
        "largest_block_cells": int(largest_block_last),
        "block_factorisations": int(n_factorised),
        "block_factorisations_reused": int(n_reused),
    }
    if verbose:
        print(f"  phases: {result['phase_seconds']}")
        print(f"  the free set was {n_blocks_last:,} blocks, largest "
              f"{largest_block_last:,} cells; {n_factorised:,} block "
              f"factorisations against {n_reused:,} reused")
        print(f"  the WHOLE matrix repeated on {value_repeats} of "
              f"{len(trace)} passes (pattern on {pattern_repeats})", flush=True)
    # CHECKED AT THE END AND NOT PER PASS. An early pass can hold a head far
    # from the answer, so a floor computed there says nothing about the solve;
    # what matters is whether the CONVERGED configuration can be certified.
    if unconfined and result.get("residual_bar", 0.0) >= CLOSURE_TOLERANCE:
        raise SystemExit(
            f"the balance's round-off floor is "
            f"{result['residual_bar']:.3e}, at or above the closure bar "
            f"{CLOSURE_TOLERANCE:.0e}: this aquifer is thick enough that the "
            "head differences carrying the flow are lost inside the heads "
            "themselves, so no convergence this solve reports could be "
            "certified. Reduce the aquifer thickness, or measure the head "
            "against a local datum rather than sea level.")

    depth = np.clip(surface_m - head, 0.0, None)
    flux = trans * (head[dst] - head[src])
    divergence = geom_divergence(n, src, dst, flux)
    pinned = conductive & ~free
    et_final = (et_rate(et_max_m_s, et_lambda_m, surface_m, head, conductive)
                * area_m2) if et_max_m_s is not None else np.zeros(n)
    seepage = np.zeros(n)
    seepage[pinned] = supply[pinned] - et_final[pinned] + divergence[pinned]
    # `imposed` cells are a BOUNDARY, not excluded ground: they take their own
    # recharge and whatever flows in from their neighbours, and that is baseflow
    # leaving the groundwater system into the river or lake whose head they
    # carry. It is a fourth door and closure needs it by name -- without it the
    # inflow to every river cell simply vanishes, which is 1.8e-3 of the budget.
    excluded = (export.surface_class == LAND) & ~conductive & ~imposed
    seepage[excluded] = supply[excluded]
    baseflow = np.zeros(n)
    baseflow[imposed] = supply[imposed] + divergence[imposed]
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
        "unconfined": unconfined,
        # GW-24. A cell whose saturated column has reached `min_saturated_m` is
        # sitting on a numerical bound rather than on a solved thickness, so its
        # depth is a LOWER bound and not a value. Empty under the confined form,
        # where no such bound exists.
        "at_transmissivity_floor": (
            conductive & has_base
            & (head - np.where(has_base, aquifer_base_m, 0.0) <= min_saturated_m)
            if unconfined else np.zeros(n, bool)),
        # What the ocean GAINS. `geom_divergence` returns net inflow, so this is
        # a plus; it was a minus, which made closure subtract the coastal
        # discharge instead of adding it.
        "ocean_outflow_m3_s": float(divergence[is_ocean].sum()),
        # GW-15. Zero, and an array of zeros rather than a None, when the sink
        # is off, so closure can add it unconditionally and a reader cannot get
        # a different answer by forgetting it.
        "et_m3_s": et_final,
        "et_total_m3_s": float(et_final.sum()),
        # GW-17. Zero, and an array of zeros, when no head is imposed.
        "baseflow_m3_s": baseflow,
        "baseflow_total_m3_s": float(baseflow.sum()),
    })
    return result


def closure(result, recharge_m_s, area_m2, export: Export) -> dict:
    """Recharge in against discharge out. A conservation law, not a comparison.

    At equilibrium every drop of recharge leaves through exactly one of four
    doors: it seeps back to the surface, it crosses the coast into the ocean
    boundary, it evaporates from a shallow water table under GW-15's sink, or it
    discharges as baseflow into a river or lake holding a fixed head under
    GW-17. If they disagree the discretisation is wrong, and no amount of
    agreement elsewhere rescues it.

    **Every door added so far was added late and cost a closure failure first.**
    The sink was one, at 5.6e-4; the baseflow was another, at 1.8e-3, because
    the water flowing into a river cell had nowhere to be counted. Each is zero
    when its term is off, so this stays one conservation law rather than four
    behind flags -- and the next boundary added will need its own door.
    """
    land = export.surface_class == LAND
    supply = float((recharge_m_s[land] * area_m2[land]).sum())
    out = (float(result["seepage_m3_s"].sum()) + result["ocean_outflow_m3_s"]
           + result.get("et_total_m3_s", 0.0)
           + result.get("baseflow_total_m3_s", 0.0))
    rel = abs(out - supply) / max(abs(supply), 1e-30)
    return {
        "recharge_m3_s": supply,
        "discharge_m3_s": out,
        "seepage_m3_s": float(result["seepage_m3_s"].sum()),
        "to_ocean_m3_s": result["ocean_outflow_m3_s"],
        "to_groundwater_et_m3_s": result.get("et_total_m3_s", 0.0),
        "to_baseflow_m3_s": result.get("baseflow_total_m3_s", 0.0),
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
    # `is_ocean`, not GW-17's wider boundary set: this imposes a table rather
    # than solving one, and it has no fixed-head argument to widen it with.
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


def confined_test(n_cells=200, dx_m=500.0, k_m_s=1e-5, recharge_m_s=2.5e-10,
                  thickness_m=100.0):
    """The CONFINED solver against its own analytic answer. A regression guard.

    Steady one-dimensional drainage at constant transmissivity has
    `h(x) = R (x_out^2 - x^2) / 2 K D`, which the three-point scheme reproduces
    exactly at cell centres. It exists so that the unconfined path's arrival
    cannot quietly change the confined one: the two now share every line of
    `solve` except the reassembly, and a shared line is a shared failure.
    """
    from types import SimpleNamespace

    n = n_cells
    x = (np.arange(n) + 0.5) * dx_m
    src, dst = np.arange(n - 1), np.arange(1, n)
    export = SimpleNamespace(n_regions=n,
                             surface_class=np.full(n, LAND, dtype=np.int64))
    geom = SimpleNamespace(src=src, dst=dst, geom=np.full(n - 1, 1.0 / dx_m),
                           volume_area_m2=np.full(n, dx_m))
    fixed = np.full(n, np.nan)
    fixed[n - 1] = 0.0
    recharge = np.full(n, recharge_m_s)
    recharge[n - 1] = 0.0
    # `et_max_m_s` AND `aquifer_base_m` ARE OMITTED ON PURPOSE, and the reason
    # is different for each. The analytic head above is the solution of the
    # equation with no sink in it, so a run carrying one would not be the case
    # this compares against; and the aquifer base is what selects the unconfined
    # form, which is the other arm and is `dupuit_test`. Neither is the missing
    # argument of world-60x0's class: this is a control that must not carry the
    # term, not a second trajectory through the same problem.
    res = solve(export, geom, k0_m_s=np.full(n, k_m_s), thickness_m=thickness_m,
                recharge_m_s=recharge, surface_m=np.full(n, 1e7),
                conductive=np.ones(n, bool), max_outer=20,
                start_all_free=True, verbose=False, fixed_head_m=fixed)
    ana = recharge_m_s * (x[-1] ** 2 - x ** 2) / (2.0 * k_m_s * thickness_m)
    err = float(np.abs(res["head_m"][:n - 1] - ana[:n - 1]).max()
                / max(float(ana.max()), 1e-300))
    return {"converged": bool(res.get("converged", False)),
            "passes": int(res.get("outer_iterations", -1)),
            "analytic_relative_error": err}


def dupuit_test(n_cells=200, dx_m=500.0, k_m_s=1e-5, recharge_m_s=2.5e-10,
                b_out_m=20.0, base_slope=0.0, min_saturated_m=1.0,
                max_outer=200, verbose=True):
    """The unconfined solver against the analytic Dupuit parabola. GW-24.

    A test that can FAIL, on a case with a right answer, which is what
    `CLAUDE.md` requires of a check. Nothing about this world enters it: one
    dimension, uniform conductivity, uniform recharge, no-flow at one end and a
    fixed head at the other, which is the textbook steady unconfined problem.

        q(x) = R x,   q = -K b db/dx   =>   b(x)^2 = b_out^2 + R (x_out^2 - x^2)/K

    TWO THINGS ARE CHECKED AND THEY ARE DIFFERENT CLAIMS.

    The IDENTITY. Over a FLAT aquifer base the arithmetic-mean face thickness
    makes the discrete flux exactly `g K (u_j - u_i)` for `u = b^2/2`, so the
    nonlinear scheme and the linear one in `u` are the same matrix. Any
    difference is round-off, and the bar is `DUPUIT_IDENTITY_RELATIVE`.

    The DISCRETISATION. The solved head against the analytic parabola above,
    against `DUPUIT_ANALYTIC_RELATIVE`. This one is not exact in general and the
    bar is set for a first-order scheme.

    `base_slope` tilts the aquifer base, which is what a real one does and is
    the case the Kirchhoff identity does NOT cover: with the base at different
    elevations on the two ends of a face, `b_face (h_j - h_i)` is no longer a
    difference of any per-cell potential. Run with a slope the identity residual
    is REPORTED rather than judged, because it is measuring the size of a term
    the transform drops rather than an error in the code.
    """
    from types import SimpleNamespace

    n = n_cells
    x = (np.arange(n) + 0.5) * dx_m
    base = base_slope * (x[-1] - x)                      # 0 at the outlet
    k0 = np.full(n, k_m_s)
    area = np.full(n, dx_m)                              # unit face width
    src = np.arange(n - 1)
    dst = src + 1
    gfac = np.full(n - 1, 1.0 / dx_m)                    # width / separation

    # The analytic solution, and a surface set well clear of it so the box
    # constraint never binds: this is a test of the flow term, not of the
    # active set, which every other check here already exercises.
    b_ana = np.sqrt(b_out_m ** 2 + recharge_m_s * (x[-1] ** 2 - x ** 2) / k_m_s)
    surface = base + 3.0 * float(b_ana.max())

    export = SimpleNamespace(n_regions=n,
                             surface_class=np.full(n, LAND, dtype=np.int64))
    geom = SimpleNamespace(src=src, dst=dst, geom=gfac, volume_area_m2=area)

    fixed = np.full(n, np.nan)
    fixed[n - 1] = base[n - 1] + b_out_m
    recharge = np.full(n, recharge_m_s)
    recharge[n - 1] = 0.0

    # `et_max_m_s` is omitted for the same reason the box constraint is held
    # clear above: the Dupuit parabola is the solution with no sink, so the
    # control must not carry one. world-60x0's class is a re-solve that drops a
    # term the primary run has; this is a case whose right answer is known only
    # while the term is absent.
    res = solve(export, geom, k0_m_s=k0, thickness_m=1.0,
                recharge_m_s=recharge, surface_m=surface,
                conductive=np.ones(n, bool), max_outer=max_outer,
                start_all_free=True, verbose=False,
                fixed_head_m=fixed, aquifer_base_m=base,
                min_saturated_m=min_saturated_m)
    head = res["head_m"]

    b = saturated_thickness(head, base, min_saturated_m)
    u = kirchhoff_potential(head, base, min_saturated_m)
    b_face = 0.5 * (b[src] + b[dst])
    lhs = gfac * k_m_s * b_face * (head[dst] - head[src])
    rhs = gfac * k_m_s * (u[dst] - u[src])
    scale = max(float(np.abs(rhs).max()), 1e-300)
    identity = float(np.abs(lhs - rhs).max() / scale)

    free = np.arange(n - 1)
    err = float(np.abs(head[free] - (base[free] + b_ana[free])).max()
                / max(float(b_ana.max()), 1e-300))

    out = {
        "converged": bool(res.get("converged", False)),
        "passes": int(res.get("outer_iterations", -1)),
        "residual": float(res.get("final_residual", float("inf"))),
        "kirchhoff_identity_relative": identity,
        "analytic_relative_error": err,
        "base_slope": float(base_slope),
        "at_floor": int(res["at_transmissivity_floor"].sum()),
    }
    if verbose:
        flat = base_slope == 0.0
        print(f"  cells {n}, dx {dx_m:g} m, base slope {base_slope:g}")
        print(f"    converged {out['converged']} in {out['passes']} passes, "
              f"residual {out['residual']:.3e}")
        print(f"    Kirchhoff identity {identity:.3e}"
              + (f"   bar {DUPUIT_IDENTITY_RELATIVE:.0e}   "
                 f"{'pass' if identity < DUPUIT_IDENTITY_RELATIVE else 'MISS'}"
                 if flat else "   REPORTED, not judged: a sloping base is the "
                              "case the transform does not cover"))
        if flat:
            print(f"    against the analytic parabola {err:.3e}   bar "
                  f"{DUPUIT_ANALYTIC_RELATIVE:.0e}   "
                  f"{'pass' if err < DUPUIT_ANALYTIC_RELATIVE else 'MISS'}")
        else:
            # NOT an error. The flat-base parabola is not this case's solution,
            # so this measures how far a sloping base moves the answer -- which
            # is the quantity worth knowing, since a real aquifer base follows
            # the terrain and this is the term the Kirchhoff transform drops.
            print(f"    departure from the FLAT-base parabola {err:.3e}, "
                  f"which is the size of the term, not an error")
    return out


def factorisation_test(sizes=(40, 80, 160), verbose=True) -> dict:
    """Does keeping the factor change the answer? It must not, and this can fail.

    world-wfge. `solve` reuses a SuperLU factorisation across passes whose
    matrix has not moved, which is only licensed if `splu(A, ...).solve(b)` is
    the same arithmetic as the `spsolve(A, b, ...)` it replaces -- not close to
    it, the SAME. Both drive `dgstrf` and `dgstrs` at the same column ordering,
    so the equality is exact and a difference in the last bit is a defect rather
    than a tolerance to widen. The check is therefore `array_equal` and there is
    no bar to choose.

    AND THE CONTROL, which is the half that makes this a test. SuperLU's
    `SymmetricMode` with `diag_pivot_thresh=0` is the ordering lever this file's
    own comment records as halving the fill, and it MUST be rejected here: it
    changes the pivot sequence, so it changes the answer at round-off, and a
    solve whose active set turns on a cell sitting at the surface to the last
    bit can find a different set of cells dry. Recording it as a REJECTION is
    what says the equality above is a property of this particular substitution
    and not of any two ways of solving the same system.

    The operand is a planar-graph Laplacian carrying the coefficient span the
    real operator has -- Gleeson's -15.2 to -11.8 in log10 permeability, 3.4
    orders -- because a well-conditioned test matrix would agree bitwise under
    any pivoting and would therefore not be able to reject anything.
    """
    rng = np.random.default_rng(20260905)

    def lattice(nx, base=0):
        """One Dirichlet-anchored lattice block, labelled from `base`."""
        n = nx * nx
        g = np.arange(n).reshape(nx, nx)
        rows, cols, vals = [], [], []
        diag = np.zeros(n)
        for di, dj in ((1, 0), (0, 1)):
            a = g[:nx - di, :nx - dj].ravel()
            b = g[di:, dj:].ravel()
            t = 10.0 ** rng.uniform(-15.2, -11.8, a.size)
            rows += [a + base, b + base]
            cols += [b + base, a + base]
            vals += [-t, -t]
            np.add.at(diag, a, t)
            np.add.at(diag, b, t)
        # A Dirichlet edge, so the operator is non-singular for the same reason
        # the real one is: a boundary the water leaves by.
        edge = np.zeros(n)
        edge[:nx] = diag[:nx].mean()
        rows.append(np.arange(n) + base)
        cols.append(np.arange(n) + base)
        vals.append(diag + edge)
        return n, np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)

    out = {"identity": [], "control": [], "partition": [],
           "passes": True, "control_rejected": False}
    for nx in sizes:
        n, r, c, v = lattice(nx)
        A = sp.coo_matrix((v, (r, c)), shape=(n, n)).tocsc()
        rhs = rng.random(n)
        direct = spsolve(A, rhs, use_umfpack=False, permc_spec="MMD_AT_PLUS_A")
        kept = splu(A, permc_spec="MMD_AT_PLUS_A").solve(rhs)
        same = bool(np.array_equal(direct, kept))
        out["identity"].append({"unknowns": n, "bitwise": same})
        out["passes"] &= same
        sym = splu(A, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                   options=dict(SymmetricMode=True)).solve(rhs)
        differs = not np.array_equal(direct, sym)
        rel = float(np.abs(direct - sym).max() / max(np.abs(direct).max(), 1e-300))
        out["control"].append({"unknowns": n, "differs": differs,
                               "relative": rel})
        out["control_rejected"] |= differs
        if verbose:
            print(f"    n {n:>7,}   splu kept the factor: "
                  f"{'BITWISE IDENTICAL' if same else 'DIFFERS -- MISS'}"
                  f"   SymmetricMode control: "
                  f"{'differs by ' + format(rel, '.1e') if differs else 'IDENTICAL'}")
    # THE PARTITION ARM. `solve` splits the free set into the blocks GW-17's
    # river and lake baselevels cut it into and factorises each on its own, so
    # what has to be exact is that a block-diagonal system solved block by block
    # is the SAME arithmetic as the whole system solved at once. It is, because
    # there are no entries between blocks and therefore no fill between them, so
    # the elimination inside a block is the sequence it would have been inside
    # the whole matrix -- PROVIDED the cells inside a block keep their relative
    # order, which is why `solve` sorts them stably.
    #
    # The case is built to be hard rather than tidy: five lattices of different
    # sizes, then the whole global labelling scrambled, so the blocks interleave
    # and a partition that quietly reordered a block's cells would show up.
    offs, R, C, V, tot = [], [], [], [], 0
    for nx in (37, 53, 29, 61, 41):
        n, r, c, v = lattice(nx, base=tot)
        R.append(r)
        C.append(c)
        V.append(v)
        offs.append(tot)
        tot += n
    scramble = rng.permutation(tot)
    relabel = np.empty(tot, dtype=np.int64)
    relabel[scramble] = np.arange(tot)
    A = sp.coo_matrix((np.concatenate(V),
                       (relabel[np.concatenate(R)], relabel[np.concatenate(C)])),
                      shape=(tot, tot)).tocsc()
    rhs = rng.random(tot)
    whole = splu(A, permc_spec="MMD_AT_PLUS_A").solve(rhs)
    nblk, lab = _components(A, directed=False)
    order = np.argsort(lab, kind="stable")
    bounds = np.searchsorted(lab[order], np.arange(nblk + 1))
    Ap = A[order][:, order].tocsc()
    piece = np.empty(tot)
    for b0, b1 in zip(bounds[:-1], bounds[1:]):
        col = Ap[:, b0:b1]
        sub = sp.csc_matrix((col.data, col.indices - b0, col.indptr),
                            shape=(b1 - b0, b1 - b0))
        piece[b0:b1] = splu(sub, permc_spec="MMD_AT_PLUS_A").solve(rhs[order][b0:b1])
    same = bool(np.array_equal(whole[order], piece))
    out["partition"] = {"unknowns": tot, "blocks": int(nblk), "bitwise": same}
    out["passes"] &= same
    if verbose:
        print(f"    n {tot:>7,}   {nblk} interleaved blocks solved separately: "
              f"{'BITWISE IDENTICAL' if same else 'DIFFERS -- MISS'}")

    # The control has to fire on at least one size or the identity above is
    # untested: a test that cannot reject anything has no right answer.
    out["passes"] = bool(out["passes"] and out["control_rejected"])
    return out


# GW-15 and GW-17 are what made this necessary and what it must therefore
# carry. `--uniqueness-check` in `build_groundwater.py` re-solves the PLANET
# from the opposite initial active set, and that arm is the measurement; this is
# the same identity on a case small enough to run without a climatology, a
# build or a mesh, and it exists because the planet arm has never been shown
# going red. A check whose only recorded value is a pass is not evidence that it
# can fail. world-qq10.
UNIQUENESS_CONTROLS = {
    # THE DEFECT THIS ISSUE IS NAMED FOR, as a mutation of the SOLVER rather
    # than of its arguments. `PINNED MEANS AT THE SURFACE` is an invariant of
    # the forward trajectory: the head is seeded at the local sink balance, so
    # every cell that starts below the surface must start FREE. Declaring them
    # all pinned instead leaves the seepage accounting reading a head that is
    # not where it says it is, and the two trajectories then both converge, both
    # close, and land apart. That is the 2.8 m recorded in
    # `hydrography/notes/water-table-convergence.md`.
    "pinned_seed": (
        "free = conductive.copy() if start_all_free else (conductive & (head < surface_m))",
        "free = conductive.copy() if start_all_free else np.zeros(n, bool)",
    ),
    # A THIRD TRAJECTORY, and it must NOT be rejected. Half the conductive cells
    # free at random, with the pinned half put back at the surface so the
    # invariant above is kept. The solution of the complementarity problem does
    # not know which cells the iteration started from, so this has to land where
    # the other two did. It is the arm that says the identity is about
    # trajectories and not about the two particular ones the check happens to
    # run.
    "third_start": (
        "free = conductive.copy() if start_all_free else (conductive & (head < surface_m))",
        "free = conductive & (np.random.default_rng(20260825).random(n) < 0.5)\n"
        "    head = np.where(free, head, surface_m)\n"
        "    head[is_ocean] = sea_level_m\n"
        "    if fixed_head_m is not None:\n"
        "        head[imposed] = fixed_head_m[imposed]",
    ),
}


def _mutant_solve(name: str):
    """`solve` with one named line replaced, compiled against this module.

    The controls have to be mutations of the code under test rather than of a
    copy of it, so the source is read from `solve` itself and the anchor line is
    required to be present. If a refactor moves the anchor this raises instead
    of quietly running an unmutated solver and reporting a control that passed,
    which is the failure mode a hand-maintained copy has.
    """
    import inspect
    import textwrap

    old, new = UNIQUENESS_CONTROLS[name]
    src = textwrap.dedent(inspect.getsource(solve))
    if src.count(old) != 1:
        raise SystemExit(
            f"the uniqueness control {name!r} anchors on a line of `solve` that "
            f"appears {src.count(old)} times; the control cannot be built and "
            "must not be reported as passing")
    ns: dict = {}
    exec(compile(src.replace(old, new), f"<control {name}>", "exec"),
         dict(globals()), ns)
    return ns["solve"]


def uniqueness_case(nx=48, ny=48, dx_m=2000.0, seed=0):
    """A synthetic problem carrying every term the planet's solve carries.

    The point of the case is coverage, not realism: an identity is a statement
    about the solver, so what the case has to do is reach the code paths a
    trajectory could differ on. It carries

      - an ocean fixed-head boundary down one edge,
      - GW-17's imposed local baselevels, as a river of fixed heads inland,
      - GW-15's evapotranspiration sink, spatially varying, so the matrix is
        head-dependent and the solve is a Newton iteration rather than one
        factorisation,
      - conductivity across three orders of magnitude, which is Gleeson's
        within-class spread,
      - a wet, tight patch where the box constraint binds and cells stay pinned
        at the surface,
      - a block enclosed by regions the lithology excludes, which is what the
        static dry-set query exists for: whether a block LOOKS cut off was once
        decided during the iteration, and that is the trajectory dependence the
        identity caught at 12.76 m.
    """
    from types import SimpleNamespace

    rng = np.random.default_rng(seed)
    n = nx * ny
    idx = np.arange(n).reshape(ny, nx)
    x, y = np.meshgrid(np.arange(nx) * dx_m, np.arange(ny) * dx_m)

    surface = (400.0 * (x / x.max())
               + 60.0 * np.sin(2 * np.pi * y / (ny * dx_m) * 3)
               * np.cos(2 * np.pi * x / (nx * dx_m) * 2)).ravel()
    sclass = np.full(n, LAND, dtype=np.int64)
    ocean = idx[:, 0].ravel()
    sclass[ocean] = LAND + 1
    surface[ocean] = -10.0

    src = np.concatenate([idx[:, :-1].ravel(), idx[:-1, :].ravel()])
    dst = np.concatenate([idx[:, 1:].ravel(), idx[1:, :].ravel()])
    export = SimpleNamespace(n_regions=n, surface_class=sclass)
    geom = SimpleNamespace(src=src, dst=dst,
                           geom=np.full(src.size, 1.0),   # square cells: w/l = 1
                           volume_area_m2=np.full(n, dx_m * dx_m))

    land = sclass == LAND
    k0 = np.where(land, 10.0 ** rng.uniform(-7.0, -4.0, n), 0.0)
    recharge = np.where(land, 1.0e-9 * (0.2 + 1.6 * (y.ravel() / y.max())), 0.0)
    et_max = np.where(land, 2.0e-8 * (0.5 + x.ravel() / x.max()), 0.0)

    fixed = np.full(n, np.nan)
    river = idx[:, nx // 2].ravel()
    fixed[river] = surface[river] - 2.0

    conductive = land.copy()
    ring = np.zeros((ny, nx), bool)
    ring[30:39, 30:39] = True
    ring[31:38, 31:38] = False
    conductive[idx[ring]] = False

    k0[idx[2:14, 34:46]] = 1e-8
    recharge[idx[2:14, 34:46]] = 6.0e-9

    return export, geom, dict(
        k0_m_s=k0, thickness_m=np.full(n, 100.0), recharge_m_s=recharge,
        surface_m=surface, conductive=conductive, max_outer=200,
        et_max_m_s=et_max, et_lambda_m=1.0, fixed_head_m=fixed)


def _head_difference(a, b, conductive):
    d = np.abs(a["head_m"] - b["head_m"])[conductive]
    scale = max(float(np.abs(b["head_m"][conductive]).max()), 1.0)
    return float(d.max()), float(d.max() / scale)


def uniqueness_test(verbose=True) -> dict:
    """The uniqueness identity, with controls that must be rejected.

    THE BAR IS `SCHEME_HEAD_RELATIVE` AND THERE IS NO SECOND ONE. A control is
    rejected when it fails the same criterion the identity is judged by, so
    nothing here is a threshold chosen after a result was seen: the controls
    either miss the bar the check already had or they do not.

    Four arms.

    IDENTITY. Two trajectories, one argument list, all terms on. The matrix is
    fixed under the confined form, so the complementarity problem has exactly
    one solution and the two heads must agree.

    THIRD START. A third, arbitrary initial active set. It must ALSO pass, and
    it is here because an identity that only ever compares two hand-chosen
    starts can be satisfied by a solver that is merely deterministic.

    DROPPED TERM. The second trajectory solved without GW-15's sink, and again
    without GW-17's baselevels. This is world-qq10's defect exactly: the
    re-solve was given neither for as long as both existed, so the arm was
    comparing two MODELS. Both must be rejected, and the size of the miss is the
    measurement of what the broken arm was reporting on.

    MUTATED SOLVER. `pinned_seed`, a one-line change to `solve` that breaks the
    invariant tying a pinned cell's head to its surface. It must be rejected.
    This is the arm that says the identity has teeth against the code rather
    than against its arguments.
    """
    export, geom, kwargs = uniqueness_case()
    cond = kwargs["conductive"]
    bar = SCHEME_HEAD_RELATIVE
    out: dict = {"criterion_relative_head": bar, "arms": {}}

    def record(name, a, b, must_pass, note=""):
        gap_m, rel = _head_difference(a, b, cond)
        converged = bool(a.get("converged")) and bool(b.get("converged"))
        passed = converged and rel < bar
        ok = passed if must_pass else not passed
        out["arms"][name] = {"max_absolute_head_difference_m": gap_m,
                             "relative_head_difference": rel,
                             "both_converged": converged,
                             "within_bar": passed,
                             "must_pass": must_pass, "as_expected": ok}
        if verbose:
            want = "must pass" if must_pass else "must be REJECTED"
            print(f"    {name:<22} {gap_m:10.3e} m   relative {rel:9.3e}   "
                  f"{want:<16} {'OK' if ok else 'WRONG VERDICT'}"
                  + (f"   {note}" if note else ""))
        return ok

    if verbose:
        print("UNIQUENESS: the identity, and the controls that must fail it")
        print(f"  criterion, unchanged: relative head difference < {bar:.0e}, "
              "and both\n  trajectories converged. A control is rejected by "
              "that same criterion;\n  there is no separate bar for one.")
        print(f"  case: {export.n_regions:,} cells, sink on, imposed "
              f"baselevels on, an enclosed\n  block and a pinned wet patch")

    forward = solve(export, geom, **kwargs, verbose=False)
    reverse = solve(export, geom, **kwargs, start_all_free=True, verbose=False)
    ok = record("identity", reverse, forward, True)

    third = _mutant_solve("third_start")(export, geom, **kwargs, verbose=False)
    ok &= record("third_start", third, forward, True,
                 "a third arbitrary active set")

    no_sink = dict(kwargs, et_max_m_s=None, et_lambda_m=None)
    ok &= record("dropped_sink", solve(export, geom, **no_sink,
                                       start_all_free=True, verbose=False),
                 forward, False, "GW-15 missing from the re-solve")

    no_base = dict(kwargs, fixed_head_m=None)
    ok &= record("dropped_baselevels", solve(export, geom, **no_base,
                                             start_all_free=True, verbose=False),
                 forward, False, "GW-17 missing from the re-solve")

    mutant = _mutant_solve("pinned_seed")
    ok &= record("mutated_solver", mutant(export, geom, **kwargs,
                                          start_all_free=True, verbose=False),
                 mutant(export, geom, **kwargs, verbose=False), False,
                 "pinned cells seeded below their surface")

    # STATED RATHER THAN LEFT TO BE FOUND. The explicit anchoring path, which
    # pins one cell of a free block that reaches no fixed head, is not exercised
    # here and no control is offered on it: a cell carrying GW-15's sink anchors
    # itself, so with the sink on -- the default -- a block only reaches that
    # code where `et_max` is zero. What a control on the anchor RELEASE would
    # take is tracked separately.
    out["passes"] = bool(ok)
    return out


def instrument_report(mesh_arm: str = "auto") -> int:
    """GW-8, printed. Every criterion is stated before the number that meets it.

    The verdict this returns is on the two claims, not on the 12% itself: the
    truncation figure is not in dispute and is reprinted here unchanged. What
    is under test is whether it is the noise floor of a SOLVED field, which is
    the reading that would put GW-18 and GW-24 out of reach before either is
    built.
    """
    import json

    import builds

    export_dir = builds.mesh_export()
    jitter = float(json.loads((export_dir / "manifest.json")
                              .read_text(encoding="utf-8"))["params"]["jitter"])
    radius_km = Export(export_dir).radius_km

    print("INSTRUMENT: what the operator's error is an error in. GW-8.")
    print("  criteria, declared before the run:")
    print(f"    truncation order in the cell size  < {TRUNCATION_ORDER_MAX}")
    print(f"    solution   order in the cell size  > {SOLUTION_ORDER_MIN}")
    print(f"    response error, relative to the response  < {RESPONSE_ERROR_MAX}")
    print(f"    response error, at two amplitudes, within a factor of "
          f"{RESPONSE_AMPLITUDE_INVARIANCE}")
    print(f"    solved head's mesh error on the mesh in use  < {HEAD_MESH_FLOOR}")
    print(f"    linear solve finished to a relative residual of "
          f"{POISSON_SOLVER_RELATIVE:.0e}")

    print(f"\n  the sweep, on Orogen's generator at the active build's own "
          f"jitter of {jitter}\n  and its own radius:")
    sweep = instrument_sweep(jitter=jitter, radius_km=radius_km)
    ok = True
    print("\n  fitted orders in the cell size:")
    for l, o in sorted(sweep["orders"].items()):
        good = (o["truncation"] < TRUNCATION_ORDER_MAX
                and o["solution"] > SOLUTION_ORDER_MIN)
        ok &= good
        print(f"    l = {l}   truncation {o['truncation']:6.3f}   "
              f"solution {o['solution']:6.3f}   "
              f"{'as claimed' if good else 'CLAIM 1 REFUTED'}")

    print("\n  the response to a transmissivity perturbation, against its "
          "analytic answer,\n  on the largest mesh of the sweep:")
    geom = Geometry(orogen_mesh(INSTRUMENT_SIZES[-1], jitter=jitter,
                                radius_km=radius_km))
    resp = response_error(geom)
    for a in sorted(k for k in resp if k != "amplitude_ratio"):
        r = resp[a]
        print(f"    amplitude {a:6.3f}   error relative to the response "
              f"{r['relative_to_response']:8.5f}   to the field "
              f"{r['relative_to_field']:9.2e}")
    size_ok = all(resp[a]["relative_to_response"] < RESPONSE_ERROR_MAX
                  for a in resp if a != "amplitude_ratio")
    amp_ok = resp["amplitude_ratio"] < RESPONSE_AMPLITUDE_INVARIANCE
    ok &= size_ok and amp_ok
    print(f"    amplitude ratio {resp['amplitude_ratio']:.6f}   "
          f"{'size and linearity as claimed' if (size_ok and amp_ok) else 'CLAIM 2 REFUTED'}")
    del geom

    if mesh_arm == "skip":
        print("\n  the real-mesh arm was skipped, so the sweep's own "
              "extrapolation is all there is\n  and the decision bar is not "
              "answered on the mesh in use.")
        return 0 if ok else 1

    root = export_dir if mesh_arm == "auto" else Path(mesh_arm)
    if not (root / "manifest.json").is_file():
        root = builds.mesh_export_of(root)
    print(f"\n  a real export, {root.parent.name}/{root.name}. THE ARM'S COST IS "
          f"A DIRECT\n  FACTORISATION OF THE WHOLE MESH, not of a free set: "
          f"name a smaller build\n  than the active one where the fill will "
          f"not fit.")
    export = Export(root)
    geom = Geometry(export)
    print(f"    {export.n_regions:,} regions, Voronoi area closes on the "
          f"sphere to {geom.closes_on_sphere():.8f}")
    trunc = laplace_beltrami_error(geom, (2, 3, 4))
    soln = poisson_solution_error(geom, (2, 3, 4))
    for l in (2, 3, 4):
        r = soln[l]
        good = r["solver_converged"] and r["relative_rms"] < HEAD_MESH_FLOOR
        ok &= good
        note = ("below the decision bar" if good
                else "ABOVE THE DECISION BAR"
                if r["relative_rms"] >= HEAD_MESH_FLOOR
                else "below the bar, but the solve MISSES its declared residual")
        print(f"    l = {l}   truncation {trunc[l]:8.5f}   solved head "
              f"{r['relative_rms']:10.3e}   residual "
              f"{r['solver_relative_residual']:8.2e}   {note}")

    # The bar is about the mesh the world is solved on. Where the arm ran on a
    # smaller build than the configured one, the sweep's own fitted order
    # carries the number across, and the extrapolation is labelled as one.
    configured = Export(export_dir).n_regions
    if configured != export.n_regions:
        ratio = np.sqrt(export.n_regions / configured)
        print(f"    the configured build carries {configured:,} regions, so its "
              f"cells are {ratio:.2f}x the size\n    of this arm's. Carried "
              f"across at the sweep's own fitted order, the solved head\n    "
              f"there is")
        for l in (2, 3, 4):
            o = sweep["orders"][l]["solution"]
            print(f"      l = {l}   {soln[l]['relative_rms'] * ratio ** o:10.3e}"
                  f"   (extrapolated, order {o:.2f})")
    return 0 if ok else 1


def groundwater_receiver(n, src, dst, flux, land):
    """Per land region, the neighbour taking the largest outgoing flux.

    The groundwater analogue of the surface network's `receiver`.

    THE SIGN IS THE WHOLE OF THIS FUNCTION. `geom_divergence` and everything
    built on it take a positive face flux to be flow from `dst` INTO `src`, so
    what LEAVES `src` across that face is `-flux`. This read the sign the other
    way and traced every cell to the neighbour it receives most water FROM,
    which is to say it followed the water uphill. The catchment check sat at
    71.6% against a 95% bar and survived a re-specification of the test aimed at
    the wrong cause, because a backwards trace and a mis-specified comparison
    look identical from the outside.
    """
    best = np.zeros(n)
    receiver = np.full(n, -1, dtype=np.int64)
    for a, b, q in ((src, dst, -flux), (dst, src, flux)):
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
    ap.add_argument("--dupuit-test", action="store_true",
                    help="GW-24: the unconfined solver against the analytic "
                         "Dupuit parabola, on a synthetic one-dimensional "
                         "aquifer. Needs no mesh and no climate.")
    ap.add_argument("--instrument", action="store_true",
                    help="GW-8: what the operator's truncation error is an "
                         "error IN. A resolution sweep on meshes from Orogen's "
                         "own generator, separating truncation from solution "
                         "error, and the operator's error on a transmissivity "
                         "RESPONSE against an analytic answer.")
    ap.add_argument("--instrument-mesh", default="auto",
                    help="the real-mesh arm of --instrument: 'auto' for the "
                         "configured build, 'skip', or an export directory")
    ap.add_argument("--factorisation-test", action="store_true",
                    help="world-wfge: that keeping a SuperLU factorisation "
                         "across passes is the same arithmetic as solving "
                         "afresh, bitwise, with the ordering lever that is NOT "
                         "as the control. Needs no mesh and no climate.")
    ap.add_argument("--uniqueness-test", action="store_true",
                    help="world-qq10: the uniqueness identity on a synthetic "
                         "case carrying GW-15's sink and GW-17's baselevels, "
                         "with the controls it must reject. The planet arm is "
                         "build_groundwater.py --uniqueness-check; this is the "
                         "one that can be shown going red, and needs no mesh "
                         "and no climate.")
    args = ap.parse_args()

    if args.factorisation_test:
        print("FACTORISATION REUSE: that keeping the factor is the same "
              "arithmetic, bitwise")
        print("  criterion, and there is no bar to choose: `splu(A).solve(b)` "
              "equals\n  `spsolve(A, b)` at the same ordering EXACTLY. The "
              "control is SymmetricMode,\n  which must be rejected by that "
              "same criterion.")
        r = factorisation_test()
        print("  " + ("PASS" if r["passes"] else "MISS"))
        return 0 if r["passes"] else 1

    if args.uniqueness_test:
        return 0 if uniqueness_test()["passes"] else 1

    if args.instrument:
        return instrument_report(args.instrument_mesh)


    if args.dupuit_test:
        print("DUPUIT: the unconfined transmissivity against a case with an "
              "analytic answer")
        print(f"  criteria, declared before the run: Kirchhoff identity < "
              f"{DUPUIT_IDENTITY_RELATIVE:.0e} on a flat base, head against "
              f"the analytic parabola < {DUPUIT_ANALYTIC_RELATIVE:.0e}")
        ok = True
        print("  the CONFINED path first, against its own analytic answer, so "
              "that the unconfined\n  arrival cannot quietly change it:")
        for D in (100.0, 2000.0, 20000.0):
            c = confined_test(n_cells=100, dx_m=1000.0, thickness_m=D)
            good = c["converged"] and c["analytic_relative_error"] < DUPUIT_ANALYTIC_RELATIVE
            ok &= good
            print(f"    D {D:8.0f} m   {c['passes']} pass, error "
                  f"{c['analytic_relative_error']:.2e}   "
                  f"{'pass' if good else 'MISS'}")
        print()
        for cells in (100, 200):
            r = dupuit_test(n_cells=cells, dx_m=100000.0 / cells)
            ok &= (r["converged"]
                   and r["kirchhoff_identity_relative"] < DUPUIT_IDENTITY_RELATIVE
                   and r["analytic_relative_error"] < DUPUIT_ANALYTIC_RELATIVE)

        # WHAT THE PICARD ITERATION COSTS, as a function of the one number that
        # governs it. The saturated thickness contrast across the domain sets
        # both how strongly the transmissivity varies -- which is the whole
        # point of the unconfined form -- and how fast the iteration converges.
        # They are the same number, so this table is the price list.
        print("\n  passes against the saturated-thickness contrast, which is "
              "the same number\n  that sets how much the unconfined form "
              "changes the answer at all:")
        print(f"    {'b at outlet':>12} {'b_max/b_out':>12} {'passes':>7} "
              f"{'converged':>10} {'head error':>11}")
        for b_out in (20.0, 100.0, 400.0, 1000.0, 2000.0):
            b_max = np.sqrt(b_out ** 2 + 2.5e-10 * (99500.0 ** 2 - 500.0 ** 2) / 1e-5)
            try:
                r = dupuit_test(n_cells=100, dx_m=1000.0, b_out_m=b_out,
                                max_outer=60, verbose=False)
            except SystemExit:
                print(f"    {b_out:12.0f} {b_max / b_out:12.2f} "
                      f"{'--':>7} {'refused':>10} "
                      f"{'round-off floor':>15}")
                continue
            print(f"    {b_out:12.0f} {b_max / b_out:12.2f} {r['passes']:7d} "
                  f"{str(r['converged']):>10} {r['analytic_relative_error']:11.2e}")

        print("\n  and the case the transform does NOT cover, for size:")
        dupuit_test(n_cells=200, dx_m=500.0, base_slope=1e-3)

        # THE FLOOR GUARD IS ITSELF A CHECK THAT CAN FAIL. This strip carries a
        # few tens of microlitres a second under half a kilometre of head, so
        # its balance is lost in its own heads and the solve must refuse rather
        # than report a convergence. A real mesh is nowhere near it -- there the
        # denominator is a planet's recharge -- and that is the point of showing
        # where the boundary is.
        print("\n  the round-off guard, on a case built to trip it:")
        try:
            dupuit_test(n_cells=400, dx_m=250.0, verbose=False)
            print("    NOT TRIPPED, and it should have been: the balance's "
                  "round-off floor no longer stops an uncertifiable solve")
            ok = False
        except SystemExit as exc:
            print(f"    refused, correctly: {str(exc).splitlines()[0][:96]}...")
        return 0 if ok else 1

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
