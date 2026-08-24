"""How a model field crosses a change of horizontal resolution.

Worldbuilding frame: these operate on the saved fields of the Vesper climate
model. Nothing here is about the real world.

TWO OPERATORS, because the model holds its state two ways.

SPECTRAL. A triangular truncation is a set of (m,n) modes, and the same mode
means the same thing at every truncation. Going up, the modes above the source
truncation are exactly zero; going down, they are discarded and their norm is
what the conversion cost. What must NOT be done is a prefix copy of the packed
array: the offset of a later `m` block moves with the truncation even though
the shared coefficients do not.

GRIDPOINT. Longitude is uniform and periodic; latitude is Gaussian, and the
quadrature weights ARE the cell areas in the sine-latitude measure, so the cell
edges are the cumulative sum of the weights. The remap is separable, and the
operator that conserves a global integral and the operator that area-averages
an intensity are the same operator -- the difference is only what gets
reported.

LONGITUDE IS NEVER IN DEGREES HERE. Cell positions are fractions of a turn,
taken from the cell index alone. Both grids put a cell centre at index zero, so
index zero means the same meridian at every resolution and no degree
convention has to be agreed on to make them line up. CLAUDE.md rule 3, and the
one-grid-convention check in `scripts/smoke_test.py`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from gauss_weight_reference import gauss_legendre


# ---------------------------------------------------------------------------
# Spectral
# ---------------------------------------------------------------------------

def triangular_modes(ntru: int) -> list[tuple[int, int]]:
    """The (m,n) modes of a triangular truncation, in the model's packed order.

    `m = 0...NTRU`, `n = m...NTRU`, which is the order `mpputsp` writes and
    `legmod` indexes. The two real components of each mode are adjacent, so the
    flat offset of component c of the k-th mode is `2*k + c`.
    """
    return [(m, n) for m in range(ntru + 1) for n in range(m, ntru + 1)]


def project_spectral(values: np.ndarray, ntru_src: int, ntru_tgt: int):
    """Project one level of a packed spectral field onto the target triangle.

    Returns the target level and the (L2 norm, largest magnitude) of whatever
    was discarded, both zero when the target truncation is the higher.
    """
    src_modes = triangular_modes(ntru_src)
    tgt_modes = triangular_modes(ntru_tgt)
    if values.size != 2 * len(src_modes):
        raise ValueError(
            f"a T{ntru_src} level holds {2 * len(src_modes)} reals and this "
            f"one holds {values.size}")
    src_at = {mode: k for k, mode in enumerate(src_modes)}
    out = np.zeros(2 * len(tgt_modes), dtype=np.float64)
    kept = set()
    for k, mode in enumerate(tgt_modes):
        j = src_at.get(mode)
        if j is None:
            continue                       # a mode the source never resolved
        out[2 * k:2 * k + 2] = values[2 * j:2 * j + 2]
        kept.add(mode)
    dropped = [k for k, mode in enumerate(src_modes) if mode not in kept]
    if not dropped:
        return out, (0.0, 0.0)
    lost = np.concatenate([values[2 * k:2 * k + 2] for k in dropped])
    return out, (float(np.sqrt(np.sum(lost * lost))), float(np.max(np.abs(lost))))


# ---------------------------------------------------------------------------
# Gaussian grid
# ---------------------------------------------------------------------------

def latitude_edges(nlat: int) -> np.ndarray:
    """Cell edges in sine latitude, from +1 at the north to -1 at the south.

    The Gaussian weights sum to two and each one IS its row's extent in this
    measure, so the cumulative sum gives the edges exactly. Computed in
    extended precision by `gauss_weight_reference.py`, which exists because
    three double-precision implementations of these weights disagree at 1e-10
    and could not settle it between themselves.
    """
    _, weights = gauss_legendre(nlat)
    edges = np.empty(nlat + 1, dtype=np.longdouble)
    edges[0] = np.longdouble(1.0)
    edges[1:] = np.longdouble(1.0) - np.cumsum(weights)
    # The last edge is -1 up to the accumulated rounding of nlat additions.
    edges[-1] = np.longdouble(-1.0)
    return edges


def _overlaps_1d(edges_src: np.ndarray, edges_tgt: np.ndarray):
    """(target index, source index, overlap) for two monotone edge arrays."""
    ti, si, wt = [], [], []
    for b in range(len(edges_tgt) - 1):
        hi_t, lo_t = edges_tgt[b], edges_tgt[b + 1]
        for a in range(len(edges_src) - 1):
            hi_s, lo_s = edges_src[a], edges_src[a + 1]
            ov = min(hi_t, hi_s) - max(lo_t, lo_s)
            if ov > 0:
                ti.append(b)
                si.append(a)
                wt.append(float(ov))
    return np.array(ti), np.array(si), np.array(wt)


def _circle_overlaps(nlon_src: int, nlon_tgt: int):
    """(target column, source column, overlap) on the periodic longitude circle.

    Positions are fractions of a turn taken from the cell index, so nothing
    here has a degree in it and the two grids share their origin by
    construction: cell zero is centred on the same meridian at every
    resolution. Cell i of an n-column grid spans (i-0.5)/n to (i+0.5)/n.

    The two grids start half a cell apart and neither start is the other's, so
    the source is tiled over the turn before and the turn after and clipped to
    the target's own turn. Without that tiling the cells at the seam find no
    partner: they are the ones that wrap.
    """
    # Measure from the target's first edge, so the target covers exactly [0,1].
    origin = -0.5 / nlon_tgt
    edges_s = (np.arange(nlon_src + 1, dtype=np.float64) - 0.5) / nlon_src - origin
    edges_t = (np.arange(nlon_tgt + 1, dtype=np.float64) - 0.5) / nlon_tgt - origin
    ti, si, wt = [], [], []
    for b in range(nlon_tgt):
        lo_t, hi_t = edges_t[b], edges_t[b + 1]
        for turn in (-1.0, 0.0, 1.0):
            for a in range(nlon_src):
                ov = (min(hi_t, edges_s[a + 1] + turn)
                      - max(lo_t, edges_s[a] + turn))
                if ov > 0.0:
                    ti.append(b)
                    si.append(a)
                    wt.append(ov)
    return np.array(ti), np.array(si), np.array(wt)


@dataclass(frozen=True)
class RemapWeights:
    """Source-cell to target-cell overlap areas, built once and reused.

    Sparse by construction: a target cell overlaps only the few source cells
    its edges cross, so the triple is a few times the target cell count however
    far apart the two resolutions are.
    """

    tgt: np.ndarray          # flat target cell index
    src: np.ndarray          # flat source cell index
    area: np.ndarray         # overlap area, in the same measure as cell area
    nugp_src: int
    nugp_tgt: int
    shape_src: tuple
    shape_tgt: tuple

    @property
    def identity(self) -> bool:
        return self.shape_src == self.shape_tgt


def build_weights(nlat_src: int, nlat_tgt: int) -> RemapWeights:
    """Overlap areas between two Gaussian grids, as a sparse triple."""
    nlon_src, nlon_tgt = 2 * nlat_src, 2 * nlat_tgt
    lat_t, lat_s, lat_w = _overlaps_1d(latitude_edges(nlat_src),
                                       latitude_edges(nlat_tgt))
    lon_t, lon_s, lon_w = _circle_overlaps(nlon_src, nlon_tgt)
    # The outer product of the two one-dimensional overlaps is the area
    # overlap, because the grids are separable in latitude and longitude.
    tgt = (lat_t[:, None] * nlon_tgt + lon_t[None, :]).ravel()
    src = (lat_s[:, None] * nlon_src + lon_s[None, :]).ravel()
    area = (lat_w[:, None] * lon_w[None, :]).ravel()
    return RemapWeights(tgt=tgt, src=src, area=area,
                        nugp_src=nlat_src * nlon_src,
                        nugp_tgt=nlat_tgt * nlon_tgt,
                        shape_src=(nlat_src, nlon_src),
                        shape_tgt=(nlat_tgt, nlon_tgt))


def remap(field: np.ndarray, weights: RemapWeights, *,
          mask: np.ndarray | None = None):
    """Area-weighted remap of one level, optionally restricted to a class.

    The same operator serves an intensity and a per-area stock: dividing by the
    target cell's own overlapped area gives the area-weighted mean, and that
    mean times the target area sums to the source integral exactly. What
    differs between the two is only what the caller then reports.

    `mask` selects the source cells eligible to contribute -- land under the
    source land mask, ocean under the ocean one -- so a coastline change cannot
    average a soil temperature into a sea surface temperature. Target cells
    with no eligible source overlap are returned as NaN and named, for the
    caller's declared fallback to fill.
    """
    if field.shape != (weights.nugp_src,):
        raise ValueError(
            f"a level of this source grid holds {weights.nugp_src} cells and "
            f"this one holds {field.shape}")
    eligible = np.ones(weights.src.shape, dtype=bool) if mask is None \
        else mask[weights.src]
    w = np.where(eligible, weights.area, 0.0)
    den = np.bincount(weights.tgt, weights=w, minlength=weights.nugp_tgt)
    hit = den > 0
    # Normalise the weights before they touch the field rather than dividing
    # the weighted sum afterwards. A target cell fed by exactly one source cell
    # then carries a weight of exactly 1.0 and its value survives unchanged,
    # which is what makes a same-resolution conversion bit-identical instead of
    # merely close.
    share = np.zeros_like(w)
    live = den[weights.tgt] > 0
    share[live] = w[live] / den[weights.tgt][live]
    out = np.full(weights.nugp_tgt, np.nan, dtype=np.float64)
    total = np.bincount(weights.tgt, weights=share * field[weights.src],
                        minlength=weights.nugp_tgt)
    out[hit] = total[hit]
    return out, ~hit


def cell_area(nlat: int) -> np.ndarray:
    """Relative area of every cell of a Gaussian grid, flattened row-major.

    In the same measure the overlaps use, so an inventory computed with it
    before and after a remap is comparable without a planet radius entering.
    """
    edges = latitude_edges(nlat)
    rows = np.asarray(edges[:-1] - edges[1:], dtype=np.float64)
    return np.repeat(rows / (2 * nlat), 2 * nlat)
