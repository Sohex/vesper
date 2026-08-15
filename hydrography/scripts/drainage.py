"""Resolve drainage over the native mesh.

World Orogen exports `drain_to` as raw steepest descent and flags every chain
that ends in neither the ocean nor a preserved basin sink with
`drainage_terminal == -2`. On this planet that is 63% of the land, draining into
220,649 unpreserved pits, most of them a single mesh cell of noise. That is not
an oversight: routing water is a hydrology decision and the exporter leaves it
here, the same way it leaves lake levels alone.

This module resolves it with a priority flood (Barnes, Lehman & Mulla 2014),
which fills those pits while keeping the preserved basins as genuine terminals.
The result is a complete assignment: every land region drains either to the
world ocean or to exactly one preserved basin sink.
"""

from __future__ import annotations

from dataclasses import dataclass
import heapq

import numpy as np

from orogen import Export, OCEAN, LAND

TERMINAL_OCEAN = -1
TERMINAL_NONE = -2


@dataclass
class Drainage:
    """Outcome of the priority flood.

    terminal:  per region, the basin index it drains to, TERMINAL_OCEAN for the
               world ocean, or TERMINAL_NONE for regions that are not land.
    filled_km: the surface a droplet must rise to in order to leave the region,
               i.e. the depression-filled elevation. Equal to `elevation_km`
               wherever the terrain already drained.
    """

    terminal: np.ndarray
    filled_km: np.ndarray

    def catchment_areas(self, export: Export, n_basins: int) -> np.ndarray:
        """Total land area draining to each basin, in km2."""
        area = export.cell_area
        land = export.surface_class == LAND
        out = np.zeros(n_basins, dtype=np.float64)
        sel = land & (self.terminal >= 0)
        np.add.at(out, self.terminal[sel], area[sel])
        return out


def resolve(export: Export, *, progress=None) -> Drainage:
    """Priority-flood the land surface, seeding ocean margins and basin sinks.

    Seeding both fronts at their own elevation is what keeps a preserved basin
    intact: its sink sits below its rim, so the basin front claims the interior
    long before the ocean front crosses the rim from outside. Where the two
    fronts meet is the spill point, which falls out of the algorithm rather than
    having to be imposed.
    """
    elev = export.elevation_km.astype(np.float64)
    sc = export.surface_class
    off, adj = export.adjacency
    n = export.n_regions

    land = sc == LAND
    terminal = np.full(n, TERMINAL_NONE, dtype=np.int32)
    filled = np.full(n, np.nan, dtype=np.float64)
    heap: list[tuple[float, int]] = []

    # Land regions touching the world ocean drain straight out. Done over the
    # flat CSR edge list rather than per region; the loop form costs minutes.
    src = np.repeat(np.arange(n, dtype=np.int64), np.diff(off).astype(np.int64))
    ocean_adjacent = np.zeros(n, dtype=bool)
    np.logical_or.at(ocean_adjacent, src, sc[adj] == OCEAN)
    ocean_adjacent &= land
    # Heap key carries a rank after the level so that a tie goes to the ocean.
    # Ties are common and not incidental: a divide where both fronts arrive at
    # the same level is exactly a saddle. Resolving it toward the ocean is the
    # conservative reading, because a basin that fills to its spill overflows
    # to the ocean anyway, so the ocean is the honest default base level.
    for r in np.flatnonzero(ocean_adjacent):
        terminal[r] = TERMINAL_OCEAN
        filled[r] = elev[r]
        heapq.heappush(heap, (elev[r], 0, int(r)))

    # Preserved basin sinks are terminals in their own right.
    for b in export.basins:
        s = b.sink
        if terminal[s] != TERMINAL_NONE:
            # A sink on the ocean margin would leak; keep the basin, note it.
            continue
        terminal[s] = b.index
        filled[s] = elev[s]
        heapq.heappush(heap, (elev[s], 1, s))

    seeded = len(heap)
    done = 0
    while heap:
        level, rank, r = heapq.heappop(heap)
        if level > filled[r]:
            continue  # stale heap entry
        done += 1
        if progress and done % 200_000 == 0:
            progress(done, seeded)
        t = terminal[r]
        for nb in adj[off[r]:off[r + 1]]:
            if not land[nb] or terminal[nb] != TERMINAL_NONE:
                continue
            # Water leaving nb must rise to at least this front's level.
            lvl = elev[nb] if elev[nb] > level else level
            terminal[nb] = t
            filled[nb] = lvl
            heapq.heappush(heap, (lvl, rank, int(nb)))

    return Drainage(terminal=terminal, filled_km=filled)


def spill_levels(export: Export, drainage: Drainage, n_basins: int) -> np.ndarray:
    """Level at which each basin first touches terrain draining somewhere else.

    Measured on the finished terrain, so it supersedes the catalogue's natural
    spill elevation wherever erosion has since cut the rim down. Returns the
    level, the terminal each basin overflows into, and how many basins were
    merged away to break spill cycles.
    """
    off, adj = export.adjacency
    n = export.n_regions
    term = drainage.terminal
    filled = drainage.filled_km

    src = np.repeat(np.arange(n, dtype=np.int64), np.diff(off).astype(np.int64))
    dst = adj.astype(np.int64)
    # Edges leaving a basin: source in a basin, destination somewhere else.
    inside = term[src] >= 0
    foreign = term[dst] != term[src]
    sel = inside & foreign
    if not sel.any():
        return np.full(n_basins, np.inf), np.full(n_basins, TERMINAL_OCEAN, np.int32), 0

    a, b = filled[src[sel]], filled[dst[sel]]
    edge = np.where(np.isfinite(b) & (b > a), b, a)
    owner = term[src[sel]]
    partner = term[dst[sel]]

    # Lowest exit per (basin, where it leads). Keeping the whole table rather
    # than only the single lowest exit is what makes the next step possible.
    key = owner.astype(np.int64) * (n_basins + 1) + (partner.astype(np.int64) + 1)
    uniq, inv = np.unique(key, return_inverse=True)
    pair_level = np.full(uniq.size, np.inf)
    np.minimum.at(pair_level, inv, edge)
    pair_owner = (uniq // (n_basins + 1)).astype(np.int32)
    pair_target = (uniq % (n_basins + 1)).astype(np.int32) - 1

    order = np.lexsort((pair_level, pair_owner))
    pair_owner, pair_target, pair_level = pair_owner[order], pair_target[order], pair_level[order]
    first = np.searchsorted(pair_owner, np.arange(n_basins + 1))

    # Basins that meet at a shared saddle each name the other as their outlet,
    # which is a cycle the overflow cascade cannot resolve. Physically they are
    # not spilling into each other, they are merging into one lake once all of
    # them reach the saddle. Collapse each cycle onto its lowest-sunk member and
    # give that member its next exit outside the cycle.
    #
    # Residual approximation: the survivor's hypsometry still describes only its
    # own depression, so a merged group filled past the shared saddle holds more
    # water than the curve says. Recorded in the report as `merged_basins`.
    sink_elev = np.array([b.sink_elevation_km for b in export.basins])
    excluded: dict[int, set[int]] = {}
    fixed = np.zeros(n_basins, dtype=bool)

    def pick(b_i: int) -> tuple[float, int]:
        skip = excluded.get(b_i, ())
        for k in range(first[b_i], first[b_i + 1]):
            if pair_target[k] not in skip:
                return float(pair_level[k]), int(pair_target[k])
        # Every neighbour was in the cycle we just collapsed. Keep the saddle
        # height, which is real, and send the overflow to the ocean rather than
        # leaving an infinite spill that would poison the hypsometry.
        if first[b_i] < first[b_i + 1]:
            return float(pair_level[first[b_i]]), TERMINAL_OCEAN
        return np.inf, TERMINAL_OCEAN

    out = np.full(n_basins, np.inf)
    target = np.full(n_basins, TERMINAL_OCEAN, dtype=np.int32)
    for b_i in range(n_basins):
        out[b_i], target[b_i] = pick(b_i)

    merged = 0
    for _ in range(n_basins):
        colour = np.zeros(n_basins, np.int8)
        cycles: list[list[int]] = []
        for s in range(n_basins):
            if colour[s]:
                continue
            path, u = [], s
            while u >= 0 and colour[u] == 0:
                colour[u] = 1
                path.append(u)
                u = int(target[u])
            if u >= 0 and colour[u] == 1:
                cycles.append(path[path.index(u):])
            for v in path:
                colour[v] = 2
        if not cycles:
            break
        for cyc in cycles:
            primary = min(cyc, key=lambda m: sink_elev[m])
            for m in cyc:
                if m != primary:
                    target[m] = primary
                    fixed[m] = True
                    merged += 1
            excluded.setdefault(primary, set()).update(cyc)
            out[primary], target[primary] = pick(primary)
    else:
        raise RuntimeError("spill graph cycle breaking did not terminate")

    return out, target, merged
