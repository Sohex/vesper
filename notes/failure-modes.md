# How this project goes wrong

Every bug in this list was found late, by comparing two things that were supposed
to agree. None announced itself: array shapes matched, fields looked plausible,
runs converged, and the answer was wrong. They are recorded by *class* rather
than as a changelog, because the classes recur and the instances do not.

Run `python scripts/check_consistency.py` before an expensive run and after
changing `source_build`. It mechanises the checks that come out of the classes
below, and it exists because prose in a document does not stop anything.

---

## 1. One quantity, several consumers, a fix that reaches some of them

The most expensive class here, twice over.

**`mrro`.** Four scripts integrated catchment runoff. `pedology` established that
`mrro` is river-routed net divergence rather than local generation and switched
to P − E, documenting why. `surface_water.py` reached the same conclusion
independently and also switched. `export_carve_list.py` and `carve_verdict.py`
were never updated, and the first of those is the one whose output leaves the
project and changes the terrain. Two independent corrections of one
misunderstanding, neither of which propagated to the consumer that mattered most.
Cost: an entire terrain build.

**`field_lon`.** The longitude-convention fix landed in `basin_means` and two of
its three callers. The third was `export_carve_list.py`. Again the one that
changes the terrain.

**The rule.** When a shared quantity is redefined, enumerate its consumers before
declaring the fix done — `grep` for the symbol, not for the bug. A fix that lands
where the bug was noticed is half a fix. Prefer one canonical accessor over four
call sites that each decide for themselves.

## 2. An optional argument whose absence means "do the wrong thing"

`basin_means(..., field_lon=None)` defaulted to the unremapped path. So the fix
for the antipodal bug reproduced the antipodal bug in any caller that had not
been updated — the original error wearing the shape of its own correction.

**The rule.** If a function cannot do the right thing without a parameter, make
it required. A `TypeError` at the call site is free; a wrong answer is not. The
same applies to config keys: `--runoff-source` is now explicit rather than
implied by whichever field the code happens to read.

## 3. Artifacts from different builds, paired silently

Drainage is a property of a terrain, so `basins.nc`, `regions.nc`,
`coupling_*.nc` and `surface_water.nc` are per-build. Four scripts were written
against a flat `<component>/data/` and would read whichever build was there:

- `world_state.py` reported 2,107 basins and 1,522 carves for a build with 2,540
  and 1,089.
- `carve_verdict.py` raised an `IndexError` **only because the basin counts
  differed**. Had they matched it would have paired row *i* of one terrain with
  row *i* of another and returned a number.
- `surface_water.py` and `export_carve_list.py` had the same default.

**The rule.** Anything derived from a terrain carries `terrain_hash`, and
anything consuming two such artifacts checks they agree. Use
`builds.component_data()` rather than a hardcoded path. The consistency checker
verifies this across the tree.

## 4. Physical inputs missing from an artifact's identity

`run_id` named resolution, flux, CO2, rotation, obliquity, eccentricity, physics
switches and a digest of every surface file — but not the stellar spectrum, and
it rounded flux to hundredths. Two near-misses followed: a `k25v` re-baseline
would have been written into the completed `k2` run's directory, where
`finalize()` takes `sorted(glob("MOST*"))[-1]` and copies out the last match; and
`0.945` and `0.94` both resolved to `s094`. Only the geography digest separated
them, which was luck rather than a guard.

**The rule.** Anything physical that is not in the name is a collision waiting for
the run that changes it. When adding an input that changes the answer, add it to
the identity in the same commit.

## 5. Provenance written once and never recomputed

The carve list header reported 749 carved for a file containing 1,838, and
described "0.96 S-Earth, 292.88 K" two re-baselines after both had moved. A
year length of `189.6145` sat as a literal while the baseline ran at a different
flux. `climatology_s096` was described as current long after it was not.

**The rule.** Derive provenance from the artifact at write time. If a header
states a number, compute that number in the same function that writes it.

## 6. A claim recorded as settled that measurement contradicts

- The energy residual was recorded as "a fixed offset, −0.455 ± 0.010, which
  rules out every state-dependent candidate". Measured the same way across four
  converged runs it spans −0.400 to −0.489, nine times wider, putting the
  ruled-out candidates back in play.
- `lake-representation.md` prescribed a large `dwmax` because "a large full
  bucket evaporates at open water's rate". Wetness reaches 1 above 40% *of*
  `dwmax`, so a deeper bucket needs proportionally more water; and routed river
  water never re-enters the evaporating bucket at all, so no setting sustains a
  lake.
- A stellar-spectrum correction was estimated at 0.4 to 0.7 K and measured at
  0.04 W/m², because the estimate was computed on a world with far more snow.

**The rule.** A recorded conclusion is evidence about what was true of the cases
it was measured on. Before relying on one, check whether the current
configuration is inside that envelope. Re-measure rather than inherit.

## 7. Unguarded arithmetic at a physical boundary

Switching catchment runoff to P − E made it negative for 57% of basins, which is
physically meaningful — those catchments evaporate more than they receive — but
the lake-area solve `A = R·C/(E − P + R)` was not written for it and reported a
lake covering −6.7 × 10¹⁵ percent of the planet. Earlier, a spill-level fallback
returned `inf` and produced an infinite basin capacity.

**The rule.** Clamp at the physical bound where the quantity is defined, not
where it is used. A catchment delivers zero or more, never less.

---

## The pattern behind the pattern

Almost every entry here was found by comparing two artifacts that should have
agreed, and almost none by reading code. The productive habit is to look for
quantities computed two ways and check them against each other: mesh against
grid, our routing against the exporter's, Penman against the model over ocean
cells, one build's composition against another's. Where a cross-check exists,
these bugs surface in minutes. Where none exists, they survive until something
downstream looks strange.

So when adding a component, add the cross-check with it, and prefer the check
that would have caught the last bug.
