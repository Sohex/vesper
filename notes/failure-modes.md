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

---

## Things already checked and disproved: do not re-derive these

Each cost a real investigation and each is the kind of plausible-sounding claim
that gets raised again by the next reviewer. Recorded so the answer is cheaper
than the check was.

**LPJ-GUESS does not have a daylength irradiance bias from the 30-hour day.**
`canexch.cpp:718` divides daily PAR by daylength, so energy and daylength both
carry the day-length factor and the ratio is exactly 1. The mirror-image
suggestion, scaling energy up to compensate, would introduce the bug it thinks it
is fixing.

**The stellar cycle does not cross PFT survival thresholds.** `tcmin_surv` tests
`mtemp_min20`, a twenty-year running mean of coldest-month means. A cycle shorter
than that is invisible to it. Any claim that N% of land sits within a cycle's
swing of a threshold has to reckon with the averaging first.

**`mrro` is not defective.** It is river-routed net divergence rather than local
generation, which explains the negative values, the nonzero ocean values and the
99.1% global conservation. It is the wrong field for catchment runoff, which is a
different statement, and the one that matters.

**Bedrock water fraction is worth 6-9% of AET, not 20-41%.** The larger figure
came from a test that changed a formula and a parameter together.

**Per-gridcell available water capacity does reach LPJ-GUESS.** Nothing is
spelled "AWC" or "water-holding capacity", which is why a grep for those finds
nothing and concludes the loop is open. It is not: `soilinput.cpp:41` takes the
texture path for any file with more than three columns, `soilinput.cpp:318-334`
runs a Cosby pedotransfer per gridcell on the sand and clay columns
`build_soil.py` writes, and `iforganicsoilproperties` blends organic retention in
and fails loudly without a SoilC column rather than silently doing nothing.
`vesperinput.cpp:293` then scales capacity by regolith depth per layer.

Two traps in that code. LPJ-GUESS carries soil water as a *fraction* of layer
capacity, so a layer scaled to exactly zero divides by zero and the resulting NaN
kills the gridcell **silently** -- zero LAI, zero AET, no crash. And because
water is fractional, thinner soil reads as relatively *wetter* for the same
absolute water, so a single cell can move either way through PFT competition.
Trust a controlled sweep, never one cell.

## One quantity, two meanings, three times the value

"The endorheic share of land" names two different measurements that differ by a
factor of about 3.5:

- **12.4% of land** is inside a preserved basin (`is_endorheic` cell area). This
  is the lithology and thermostat figure.
- **43.0% of land** *drains* to a closed basin. This is the hydrography figure
  and the one the carve verdict is about.

Both are correct. A basin's catchment is far larger than its floor, which is the
whole reason a small area of fill can decouple a large share of weathering. They
are named apart in `world_state.json`; quote the name, never "the endorheic
share".
