# Two things the carve verdict gets wrong, and neither is measured yet

Recorded 2026-08-16. Both concern `hydrography/scripts/export_carve_list.py`,
which turns the overflow test into the per-basin retain fraction Orogen consumes.
They were carried in `world_state.json` as curated judgement, which is a state
file and the wrong place for an argument, so they are written down here with
their evidence and tracked as HYD-1 and HYD-2 in `TASKS.md`.

Neither has been re-measured against the active build, because there is no
climatology on it. The counts below are from the verdict that was current when
each was noticed, and they are quoted as history rather than as current values.

## 1. A basin with no catchment can still fill from its own surface

The retain fraction is built from the fractional margin

    margin = (E_penman - E*) / E_penman,   E* = P + critical * runoff

which says how far open-water evaporation would have to fall before the basin
starts overflowing. With `runoff` zero the expression loses its second term and
reduces to a comparison of evaporation against precipitation alone. That sent 228
dry salt pans to be carved wide open, and the guard added at
`hydrography/scripts/export_carve_list.py:180` pins every zero-runoff basin to
retain 1 instead.

The guard's justification is that a basin with no catchment runoff receives
nothing and can never overflow. That is true of the catchment and not of the
basin. The lake surface receives precipitation directly, so where P exceeds
open-water E the water body grows, and since the imbalance does not close as the
surface expands it grows until it spills. The 228 are exactly the set where that
inequality held, which is why they were carved before the guard existed.

So both treatments are wrong at the same 228 basins, in opposite directions: the
old one carved them on an expression that had lost its runoff term, the current
one preserves them on an argument that only covers the catchment. What the set
needs is the marginal band, with a retain fraction built from the direct
precipitation balance over the lake surface rather than from the catchment one.

The size of the error is bounded and small in area -- these are pans, and 1,248
basins had zero catchment runoff in the same verdict -- but it is not bounded in
kind, because a basin that spills is a basin that stops existing.

## 2. Retain measures distance from a threshold, not the ability to cut

`retain = min(1, margin / TOLERANCE)` scales the notch Orogen cuts by how close
the basin sits to its own overflow threshold, in units of the evaporation
uncertainty we actually have. TOLERANCE is 0.25, set from the biosphere
assumption rather than from anything about the basin.

What incises an outlet is discharge. Stream power goes with the flux over the
sill, so a basin that barely trickles over and one that pours the runoff of a
large catchment through the same saddle receive the same treatment from this
mapping, when the second should cut orders of magnitude faster. The verdict is
correct in kind -- both overflow, both eventually drain -- and wrong in degree,
which matters because the retain fraction is precisely a statement of degree.

The consequence is a systematic over-cut on marginal trickle basins and an
under-cut on the large spillers, which is the opposite of the ordering the
landscape would produce. Making retain a function of overflow discharge needs the
discharge, which is the catchment integral the verdict already computes, so the
input exists and only the mapping is missing.

Not quantified: no one has compared the two orderings on a real verdict.
