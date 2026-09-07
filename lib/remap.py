"""Move a field between two grids that are not the same grid, and say what it cost.

`lib/gridding.py` owns the mesh-to-grid reduction and the column convention, and
its vocabulary -- extensive, intensive, categorical, moments, expectation -- is
the vocabulary here too. This module adds the one thing that vocabulary has no
operator for: a crossing between two DIFFERENT grids, which is what the offline
ocean needs and what nothing in the tree could do.

The crossing this exists for is ExoPlaSim's Gaussian grid to GOLDSTEIN's, and
the reason it is not simply averaging is that the two are different in KIND.
A Gaussian row is not an equal-area row, so a remap that treats either grid as
regular is wrong at the poles first and everywhere second.

## Why an offline crossing has to conserve

A synchronous coupling's interpolation error is transient and partly
self-cancelling. This coupling is offline and passes climatologies, so a
non-conservative crossing is a FIXED spurious source term in an ocean that is
then integrated for of order ten thousand model years: a permanent, spatially
structured bias with unlimited time to express itself. Conservation is therefore
a property of the weight matrix rather than an error budget to be tracked, and
`scripts/build_ocean_remap.py` checks it against an integral with a known
answer.

## Three field semantics, two normalisations, no default

Which normalisation a field takes is a property of the FIELD and is stated by
the contract, not chosen at the call site, so `apply` takes it as a required
argument out of a closed set:

    EXTENSIVE_TOTAL  an already-integrated quantity, watts or kilograms per
                     CELL. Weight is intersection over SOURCE area, so the
                     global sum is preserved exactly.
    FLUX_DENSITY     a rate or a density whose integral is a real budget --
                     net heat flux, freshwater flux, wind stress. Weight is
                     intersection over DESTINATION area, so `sum(F*A)` is
                     preserved exactly. A field of ones does NOT remap to ones
                     under this, and that is correct: what comes back is the
                     coverage, which is why coverage travels with the result.
    INTENSIVE        a state with no integral to conserve -- sea surface
                     temperature, an ice fraction, an exported surface
                     velocity. Weight is intersection over COVERED area, so a
                     field of ones remaps to ones exactly and an uncovered
                     destination contributes nothing rather than a zero.

Routing an extensive quantity through `INTENSIVE`, or a state through either of
the others, is the error `lib/gridding.py`'s two words exist to prevent. Nothing
here can detect it -- a temperature and a heat flux are both float64 arrays --
so the semantics is required, closed, and stamped into the ledger, and the
contract that names each field is `config/pipeline.yaml`'s OCN row rather than
this module.

## The masks will not agree, and the residual is coastline-shaped

The ocean's wet mask and the atmosphere's come from the same Orogen mesh at
different resolutions, so they cannot agree at the coast: a cell that is 40 per
cent land at the atmosphere's support sits inside an ocean cell that is either
wet or dry. Two things follow and both are stated rather than absorbed.

**A destination cell with no valid source over it is an ERROR.** It is a mask
disagreement, not sparse coverage, and backfilling it from the nearest source --
which `lib/gridding.py:land_weighted` does for a grid finer than the mesh, where
the substitution is defensible -- would put an invented flux into an ocean cell
that will integrate it for ten thousand years.

**A valid source cell with no valid destination under it goes to the nearest
valid destination cell**, by great-circle distance between cell centres, and the
ledger records how much area moved and how far the furthest went. This is the
rule `references/climber-x/src/geo/coast_cells.f90` takes for river discharge,
narrowed to a single neighbour; the difference here is that it has to carry heat
and salt as well as water, which is why it is the CONSERVING matrix that carries
the move and the intensive one that does not. An intensive state in an orphaned
source cell has nothing to conserve, so it is dropped and counted.

## A vector is not two scalars

The ocean returns a surface VELOCITY, so this module carries a vector semantics
beside the three scalar ones. `Crossing.apply_vector` lifts the two components
into the sphere's own three cartesian components at the SOURCE cells' frames,
remaps those three with the same weights, and projects onto the DESTINATION
cells' frames. Remapping "east" and "north" as independent scalars instead
averages components in frames that are not the same frame; the error vanishes at
the equator, grows with latitude, and is worst exactly where an ocean grid and a
Gaussian grid disagree most, which is where sea ice is.

THE SIZE OF THAT IS SET BY THE DESTINATION CELL'S LONGITUDE SPAN and not by
latitude, because the frame rotates by the span whatever row the cell is in.
Measured against a rigid rotation it is a median 0.50 per cent of the local
speed at T21 onto 36 x 36, 0.17 at T42 onto the same, and 0.070 at T85 onto
72 x 72, and it is not smaller at the equator.

THE CHECK THAT DISCRIMINATES IS NOT A FULL-SPHERE INTEGRAL, which is the part
that is not obvious. A rigid rotation has an analytically zero cartesian
integral over the sphere -- the rotation vector crossed into the integral of
position, and that integral is the origin -- and the crossing must return it.
So does the COMPONENTWISE remap, to round-off: both grids are uniform in
longitude, so the frame error is a pure phase and it cancels around every row.
Taken over a region that is not zonally symmetric it does not cancel and the two
separate by four orders. `ocean/scripts/build_ocean_grid.py` runs both forms and
`ocean/notes/vector-crossing.md` carries the measurement.

## Where the radius went

Every weight is a ratio of areas, so the planetary radius cancels and nothing
here reads `config/planet.yaml`. Areas are shares of the sphere. A caller
wanting square metres multiplies by `GridSpec.cell_area(radius)` at the end.
"""

from __future__ import annotations

import numpy as np

from gridding import GridSpec

EXTENSIVE_TOTAL = "extensive_total"
FLUX_DENSITY = "flux_density"
INTENSIVE = "intensive"
SEMANTICS = (EXTENSIVE_TOTAL, FLUX_DENSITY, INTENSIVE)
CONSERVING = (EXTENSIVE_TOTAL, FLUX_DENSITY)

# The acceptance tolerances, FIXED BEFORE ANY RESULT WAS SEEN, in the units the
# quantity is reported in. Each is set from the round-off floor of the operation
# it bounds rather than from what the implementation happened to return.
#
# CLOSURE_TOLERANCE bounds `|sum(F_dst*A_dst) - sum(F_src*A_src)| / sum(|F_src|*A_src)`.
# float64 carries 2.2e-16; numpy sums pairwise, so a sum over the 32768 cells of
# a T85 grid carries about eps*log2(N) = 3.4e-15, and each weight is a
# difference of sines good to a few ulp. 1e-12 is therefore about two and a half
# orders above the floor -- loose enough not to be a round-off alarm -- and
# three orders TIGHTER than the 1e-9 relative that `references/esmf/` holds
# itself to on analytic fields, so it is a bar this can miss.
CLOSURE_TOLERANCE = 1e-12

# CONSTANT_TOLERANCE bounds `max|remap(ones) - 1|` under INTENSIVE. This is a
# ratio of two sums of the same weights, so the floor is a few ulp and there is
# no cancellation; 1e-13 is three orders above it.
CONSTANT_TOLERANCE = 1e-13

# AREA_TOLERANCE bounds `|sum(cell_area_fraction) - 1|`, a pairwise sum of
# positive terms. Floor about 3.4e-15; 1e-13 is one and a half orders above it.
AREA_TOLERANCE = 1e-13


def _periodic_overlap(dst_edges: np.ndarray, src_edges: np.ndarray,
                      period: float) -> np.ndarray:
    """Overlap of two partitions of a circle, shaped (n_dst, n_src).

    Both edge arrays increase and each spans exactly one period, so they are two
    cuts of the same circle offset from one another. The source is first shifted
    by a whole number of periods into the destination's window, after which the
    two windows overlap within one period and only the shifts -1, 0 and +1 can
    contribute anything; the sum over those three is exact.

    The offset is REAL and is the reason this is not an index join: GOLDSTEIN's
    first column edge sits at `phi0`, whose shipped value is 260 degrees west,
    while the export's sits at 180 west. Nothing here compares one grid's
    longitude LABEL with the other's -- both edge arrays come from their own
    constructor in `lib/gridding.py` -- which is the form `CLAUDE.md` rule 3
    takes on a boundary where the identity does not hold.
    """
    dst_edges = np.asarray(dst_edges, dtype=np.float64)
    src = np.asarray(src_edges, dtype=np.float64)
    src = src - period * np.floor((src[0] - dst_edges[0]) / period)
    lo_d, hi_d = dst_edges[:-1, None], dst_edges[1:, None]
    lo_s, hi_s = src[None, :-1], src[None, 1:]
    out = np.zeros((dst_edges.size - 1, src.size - 1), dtype=np.float64)
    for shift in (-1.0, 0.0, 1.0):
        out += np.clip(np.minimum(hi_d, hi_s + shift * period)
                       - np.maximum(lo_d, lo_s + shift * period), 0.0, None)
    return out


def _interval_overlap(dst_edges: np.ndarray, src_edges: np.ndarray) -> np.ndarray:
    """Overlap of two partitions of an interval, shaped (n_dst, n_src).

    Monotonic in either direction, because a `GridSpec` carries its rows in its
    own order and this project's grids run north to south while GOLDSTEIN's `j`
    runs south to north. The bounds are taken per row rather than assumed
    ascending, so neither orientation is privileged and neither needs flipping
    before it gets here.
    """
    d = np.asarray(dst_edges, dtype=np.float64)
    s = np.asarray(src_edges, dtype=np.float64)
    lo_d = np.minimum(d[:-1], d[1:])[:, None]
    hi_d = np.maximum(d[:-1], d[1:])[:, None]
    lo_s = np.minimum(s[:-1], s[1:])[None, :]
    hi_s = np.maximum(s[:-1], s[1:])[None, :]
    return np.clip(np.minimum(hi_d, hi_s) - np.maximum(lo_d, lo_s), 0.0, None)


class Crossing:
    """The weight matrix between two grids, with its coverage and its ledger.

    Built from two `GridSpec`s and nothing else, so the geometry is exact and
    the construction is cheap: because neither grid's row boundaries depend on
    longitude and neither grid's column boundaries depend on latitude, the
    overlap AREA factorises into a longitude overlap times a sine-of-latitude
    overlap. Two one-dimensional matrices, no polygon intersection.

    A mask does NOT factorise -- it is genuinely two-dimensional -- so the
    factorisation is used to find the sparsity pattern and the masked matrix is
    then assembled sparse over it. At a T85 atmosphere against a 72 x 72 ocean
    that pattern holds of order 1e5 nonzeros, which is why the exact operator is
    cheaper than the approximate ones it replaces.
    """

    def __init__(self, src: GridSpec, dst: GridSpec):
        if not isinstance(src, GridSpec) or not isinstance(dst, GridSpec):
            raise TypeError(
                "a crossing is built from two GridSpecs, which carry cell "
                "BOUNDARIES and the constructor that made them. Centre arrays "
                "read off a file are labels, not a coordinate correspondence.")
        self.src, self.dst = src, dst
        self.lat_overlap = _interval_overlap(dst.sin_edges, src.sin_edges)
        self.lon_overlap = _periodic_overlap(dst.lon_edges, src.lon_edges,
                                             period=float(src.lon_edges[-1] - src.lon_edges[0]))
        self.src_area = src.cell_area_fraction().ravel()
        self.dst_area = dst.cell_area_fraction().ravel()
        self._overlap = None          # intensive matrix: raw restricted overlap
        self._conserving = None       # conserving matrix: orphans moved, columns closed
        self.src_valid = None
        self.dst_valid = None
        self.orphans: dict = {}

    # -- the geometry, before any mask ---------------------------------------

    def _sparse_overlap(self):
        """The unmasked intersection areas as a sparse (ncell_dst, ncell_src).

        Areas are shares of the sphere: a cell between two sines of latitude and
        two longitudes covers `dsin * dlon / 2 / full turn`, so the intersection
        of two such cells is the product of the two one-dimensional overlaps
        under the same normalisation.
        """
        from scipy import sparse
        turn = float(self.src.lon_edges[-1] - self.src.lon_edges[0])
        jd, js = np.nonzero(self.lat_overlap)
        idd, iss = np.nonzero(self.lon_overlap)
        lat_v = self.lat_overlap[jd, js]
        lon_v = self.lon_overlap[idd, iss] / (2.0 * turn)
        rows = (jd[:, None] * self.dst.nlon + idd[None, :]).ravel()
        cols = (js[:, None] * self.src.nlon + iss[None, :]).ravel()
        vals = (lat_v[:, None] * lon_v[None, :]).ravel()
        return sparse.coo_matrix((vals, (rows, cols)),
                                 shape=(self.dst.ncell, self.src.ncell)).tocsr()

    # -- the mask, and the coast ---------------------------------------------

    def restrict(self, src_valid=None, dst_valid=None) -> "Crossing":
        """Fix the two masks and build both weight matrices. Returns self.

        `src_valid` and `dst_valid` are boolean arrays shaped like their grids.
        `None` means every cell, which is the geometric crossing and is what the
        acceptance checks run on.

        Raises on a valid destination cell with no valid source over it, because
        that is a mask disagreement and an invented flux there would be
        integrated for the life of the ocean run. Moves a valid source cell with
        no valid destination under it to the nearest valid destination cell in
        the CONSERVING matrix only, and records the move.
        """
        overlap = self._sparse_overlap()
        sv = (np.ones(self.src.ncell, dtype=bool) if src_valid is None
              else np.asarray(src_valid, dtype=bool).ravel())
        dv = (np.ones(self.dst.ncell, dtype=bool) if dst_valid is None
              else np.asarray(dst_valid, dtype=bool).ravel())
        if sv.size != self.src.ncell or dv.size != self.dst.ncell:
            raise ValueError(
                f"masks must match their grids: source {sv.size} against "
                f"{self.src.ncell}, destination {dv.size} against {self.dst.ncell}")
        self.src_valid, self.dst_valid = sv, dv

        from scipy import sparse
        overlap = (sparse.diags(dv.astype(np.float64)) @ overlap
                   @ sparse.diags(sv.astype(np.float64))).tocsr()
        overlap.eliminate_zeros()

        landed = np.asarray(overlap.sum(axis=0)).ravel()      # per source cell
        received = np.asarray(overlap.sum(axis=1)).ravel()    # per destination cell
        blind = dv & (received <= 0)
        if blind.any():
            rows, cols = np.divmod(np.flatnonzero(blind), self.dst.nlon)
            raise SystemExit(
                f"{int(blind.sum())} destination cells are valid and have no "
                f"valid source over them, first at (row {rows[0]}, col {cols[0]}). "
                "That is a mask disagreement between two supports of the same "
                "mesh, not sparse coverage, and it is an error rather than a "
                "backfill: an invented flux in an ocean cell is integrated for "
                "the life of the run. Reconcile the masks, or declare the cell "
                "invalid.")

        orphan = sv & (landed <= 0)
        conserving = overlap.copy()
        moved_area = 0.0
        moved_km = 0.0
        if orphan.any():
            targets, distance = self._nearest_valid(np.flatnonzero(orphan), dv)
            extra = sparse.coo_matrix(
                (self.src_area[orphan], (targets, np.flatnonzero(orphan))),
                shape=conserving.shape)
            conserving = (conserving + extra).tocsr()
            landed = np.asarray(conserving.sum(axis=0)).ravel()
            moved_area = float(self.src_area[orphan].sum())
            moved_km = float(distance.max())
        self.orphans = {
            "source_cells_orphaned": int(orphan.sum()),
            "source_area_fraction_moved": moved_area,
            "furthest_move_degrees": moved_km,
        }

        # Close every valid source column: all of a valid source cell's content
        # lands on valid destinations. Without this a coastal source cell leaks
        # the part of itself that lies over destination land, and the leak is a
        # coastline-shaped field rather than noise.
        scale = np.zeros(self.src.ncell)
        np.divide(self.src_area, landed, out=scale, where=(landed > 0))
        self._conserving = (conserving @ sparse.diags(scale)).tocsr()
        self._overlap = overlap
        return self

    def _nearest_valid(self, orphan_cells, dst_valid):
        """Nearest valid destination cell to each orphaned source cell.

        Great-circle nearest by cell CENTRE, which is the single-neighbour form
        of what `coast_cells.f90` does with an expanding stencil. It is a rule
        and not an approximation: the flux has to arrive somewhere named, and
        the ledger reports how far the furthest one travelled so a reader can
        see when the answer has stopped being local.
        """
        from scipy.spatial import cKDTree

        def unit(spec):
            lat, lon = spec.cell_centres()
            la = np.deg2rad(lat)[:, None] * np.ones((1, spec.nlon))
            lo = np.deg2rad(lon)[None, :] * np.ones((spec.nlat, 1))
            return np.stack([(np.cos(la) * np.cos(lo)).ravel(),
                             (np.cos(la) * np.sin(lo)).ravel(),
                             np.sin(la).ravel()], axis=1)

        dst_pts = unit(self.dst)
        src_pts = unit(self.src)
        keep = np.flatnonzero(dst_valid)
        chord, which = cKDTree(dst_pts[keep]).query(src_pts[orphan_cells])
        arc = np.rad2deg(2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0)))
        return keep[which], arc

    # -- applying it ---------------------------------------------------------

    def apply(self, field, semantics: str):
        """Remap `field` onto the destination grid. Returns (field, coverage).

        `semantics` is required and closed. `coverage` is the share of each
        destination cell covered by valid source area, in the destination's
        shape, and it travels with the result because a discarded coverage
        fraction is a conservation identity that can no longer be evaluated.
        """
        if semantics not in SEMANTICS:
            raise ValueError(
                f"{semantics!r} is not a field semantics. It is one of "
                f"{', '.join(SEMANTICS)}, and there is no default because "
                "which one a field takes is a property of the field.")
        if self._overlap is None:
            self.restrict()
        f = np.asarray(field, dtype=np.float64)
        if f.shape != self.src.shape:
            raise ValueError(
                f"field is {f.shape} and the source grid is {self.src.shape}")
        flat = f.ravel()
        coverage = (np.asarray(self._overlap.sum(axis=1)).ravel() / self.dst_area)
        if semantics == INTENSIVE:
            received = self._overlap @ flat
            weight = np.asarray(self._overlap.sum(axis=1)).ravel()
            out = np.zeros(self.dst.ncell)
            np.divide(received, weight, out=out, where=(weight > 0))
        elif semantics == EXTENSIVE_TOTAL:
            # The weight is intersection over SOURCE area, so each source cell's
            # own total is split among the destination cells it lands in and the
            # column sums to one rather than to an area. Dividing the field by
            # the source area first is that normalisation: `_conserving` closes
            # its columns on the AREA, because that is what a density needs.
            out = self._conserving @ (flat / self.src_area)
        else:
            out = (self._conserving @ flat) / self.dst_area
        return out.reshape(self.dst.shape), coverage.reshape(self.dst.shape)

    def apply_vector(self, east, north, semantics: str):
        """Remap a tangent VECTOR field. Returns (east, north, coverage, radial).

        A vector does not remap like a scalar, and remapping the east and north
        components as two independent scalars is the error this exists to
        prevent: "east" at one longitude is not "east" at another, so the two
        numbers being averaged are components in different frames. What sets the
        size of it is the destination cell's LONGITUDE SPAN rather than its
        latitude, and a full-sphere integral cannot see it at all; the module
        docstring above has the measurement and the check that can.

        So direction is carried by the geometry rather than by the labels. The
        field is lifted into the sphere's own three cartesian components at the
        SOURCE cells' frames, each of the three is remapped by the SAME weights
        under the SAME semantics as an ordinary scalar -- the lift is linear, so
        whatever `apply` conserves it conserves componentwise -- and the result
        is projected onto the DESTINATION cells' frames. Both frames come from
        `GridSpec.cell_frames`, so neither grid is ever asked what the other
        calls a direction.

        `radial` is what the projection discarded: the part of the mapped vector
        that no longer lies in the destination cell's tangent plane, in the same
        units as the components. It is not an error to be corrected away and it
        is not renormalised, because scaling the horizontal part back up to the
        original magnitude would destroy the conservation the lift just bought.
        It is the honest measure of how far the two grids' tangent planes are
        apart over a cell, it is zero when the two grids coincide, and it belongs
        in the ledger beside the coverage. It is also what lets a caller rebuild
        the mapped cartesian vector exactly, which is what the acceptance
        integral is taken over.
        """
        e_s, n_s, _ = self.src.cell_frames()
        e_d, n_d, u_d = self.dst.cell_frames()
        u = np.asarray(east, dtype=np.float64)
        v = np.asarray(north, dtype=np.float64)
        if u.shape != self.src.shape or v.shape != self.src.shape:
            raise ValueError(
                f"the two components are {u.shape} and {v.shape} and the source "
                f"grid is {self.src.shape}")
        cart = e_s * u.ravel()[:, None] + n_s * v.ravel()[:, None]
        out = np.empty((self.dst.ncell, 3))
        coverage = None
        for k in range(3):
            got, coverage = self.apply(cart[:, k].reshape(self.src.shape), semantics)
            out[:, k] = got.ravel()
        return ((out * e_d).sum(axis=1).reshape(self.dst.shape),
                (out * n_d).sum(axis=1).reshape(self.dst.shape),
                coverage,
                (out * u_d).sum(axis=1).reshape(self.dst.shape))

    # -- what it cost --------------------------------------------------------

    def ledger(self, field=None, semantics: str | None = None) -> dict:
        """What the crossing carries and what it does not.

        The same shape as `lib/gridding.py:transfer_ledger`, which is the ledger
        every conversion in this project emits, plus the two facts a grid-to-grid
        crossing has that a mesh-to-grid reduction does not: the coast moves, and
        the closure residual against a field whose integral is known.

        For a conserving semantics the residual is
        `|sum(F_dst*A_dst) - sum(F_src*A_src)| / sum(|F_src|*A_src)` and its bar
        is `CLOSURE_TOLERANCE`. For `INTENSIVE` there is no integral to close, so
        what is reported instead is `max|remap(ones) - 1|` against
        `CONSTANT_TOLERANCE`: the identity that operator does have.
        """
        if self._overlap is None:
            self.restrict()
        sv, dv = self.src_valid, self.dst_valid
        coverage = (np.asarray(self._overlap.sum(axis=1)).ravel() / self.dst_area)
        out = {
            "source_grid": self.src.name,
            "source": self.src.source,
            "destination_grid": self.dst.name,
            "destination": self.dst.source,
            "source_cells_valid": int(sv.sum()),
            "destination_cells_valid": int(dv.sum()),
            "source_area_fraction_valid": float(self.src_area[sv].sum()),
            "destination_area_fraction_valid": float(self.dst_area[dv].sum()),
            "destination_coverage_min": float(coverage[dv].min()) if dv.any() else 0.0,
            "destination_coverage_mean": float(
                (coverage[dv] * self.dst_area[dv]).sum() / self.dst_area[dv].sum())
            if dv.any() else 0.0,
            "weights_nonzero": int(self._conserving.nnz),
            **self.orphans,
        }
        if field is not None:
            if semantics not in SEMANTICS:
                raise ValueError("a ledger about a field needs that field's semantics")
            out["semantics"] = semantics
            f = np.asarray(field, dtype=np.float64)
            if semantics == INTENSIVE:
                ones, _ = self.apply(np.ones(self.src.shape), INTENSIVE)
                out["constant_residual"] = float(np.abs(ones.ravel()[dv] - 1.0).max())
                out["constant_tolerance"] = CONSTANT_TOLERANCE
                out["passes"] = bool(out["constant_residual"] <= CONSTANT_TOLERANCE)
            else:
                dstf, _ = self.apply(f, semantics)
                if semantics == EXTENSIVE_TOTAL:
                    got = float(dstf.ravel()[dv].sum())
                    want = float(f.ravel()[sv].sum())
                    scale = float(np.abs(f.ravel()[sv]).sum())
                else:
                    got = float((dstf.ravel() * self.dst_area)[dv].sum())
                    want = float((f.ravel() * self.src_area)[sv].sum())
                    scale = float((np.abs(f.ravel()) * self.src_area)[sv].sum())
                out["destination_integral"] = got
                out["source_integral"] = want
                out["closure_residual_relative"] = (abs(got - want) / scale
                                                    if scale > 0 else 0.0)
                out["closure_tolerance"] = CLOSURE_TOLERANCE
                out["passes"] = bool(out["closure_residual_relative"] <= CLOSURE_TOLERANCE)
        return out
