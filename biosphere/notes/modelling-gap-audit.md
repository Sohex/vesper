# Biosphere modelling gap audit

**Recorded:** 2026-08-21  
**Scope:** the Vesper input module, LPJ-GUESS run harness, soil feedback,
productivity scoring, and the return of vegetation state to ExoPlaSim.  This was
a source and pipeline audit only.  No LPJ-GUESS build or simulation was run.

BIO-5 through BIO-10 are the activation chain for the CNP fork.  They do not
replace the C-N acceptance criteria or close the vegetation-climate loop.  The
findings below are independent of whether phosphorus limitation is enabled.

## 1. LPJ-GUESS does not share the surface's definition of habitable land

`build_lpj_driver.py` selects every cell for which the climatology has
`lsm > 0.5`.  That includes solved lakes, salt crust and playa mud because those
surfaces deliberately remain ExoPlaSim land cells.  The driver carries no
rootable-area fraction and LPJ-GUESS receives a valid soil for every such cell.

`build_surface_albedo.py` already records the consequence: LPJ-GUESS is not told
which ground is barren and will grow on it.  The script multiplies tree and grass
cover by `1 - barren_fraction` when painting albedo, but this is only a graphical
correction to the feedback field.  The unmasked vegetation still contributes to
`anpp.out`, `cpool.out`, the productivity score, and therefore the soil-carbon
iteration.  Persistent lake area is not included even in that late mask.

The required interface is one provenance-stamped **rootable surface fraction**
per climate cell, derived from the same surface classification and lake solution
used by ExoPlaSim.  A fully non-rootable cell need not be simulated.  A partially
rootable cell may still use one natural stand, but extensive quantities must be
weighted by its rootable fraction.  The same fraction must govern:

- NPP and carbon totals;
- soil-carbon feedback;
- foliar and forest cover returned to the climate;
- aerodynamic roughness; and
- prediction denominators and barren-area diagnostics.

This is BIO-11.  Applying the contract to the surface compositing itself is
BIO-17, so defining the interface does not get mixed with one consumer's code.

## 2. The last output year is not an equilibrium statistic

Three downstream consumers independently select the highest year number in an
LPJ-GUESS output file:

- `build_surface_albedo.py:read_foliar_cover`;
- `score_prediction.py:read_table`; and
- `build_soil.py:read_soil_carbon`.

The model enables stochastic establishment, stochastic mortality, generic
disturbance and fire, while the run harness defaults to only five patches.  The
registered nitrogen sweep has already measured a non-monotonic 12% NPP movement
and attributes it to patch stochasticity at `npatch 5`.  With a multi-year
stellar forcing, the last year is additionally one arbitrary phase of the
forcing cycle.

There must be one shared output-window rule.  It should reject a trending end
window and otherwise average a fixed number of complete forcing cycles; for a
one-year climatology that reduces to a fixed multi-year window.  The reducer must
be used by feedback, scoring and soil carbon, and must report temporal spread.
Patch-count or seed sensitivity belongs in the uncertainty report rather than
being hidden by a single-year value.  This is BIO-12.

## 3. The wet-days setting is not implemented by VesperInput

`global.ins` declares `ifrainonwetdaysonly 1`, and the normal monthly INTERP
paths call `prdaily` when that flag is true.  `VesperInput::interpolate` instead
calls `interp_monthly_totals_conserve`, which turns each monthly total into a
smooth sequence of positive mean-daily values, and never calls `prdaily`.
Consequently the instruction file says wet days only while Vesper effectively
rains a little almost every day.

The correction must not import an Earth observational GWGEN field.  The climate
archive or snapshot stream should provide a Vesper-native daily sequence or at
least monthly wet-day counts, and the binary driver must carry the selected
statistic explicitly.  Monthly totals must still close exactly.  Smooth and
event-distributed precipitation should be compared before the latter is adopted,
because the affected mechanisms include interception, drought stress, soil
water and GLOBFIRM fire.  This is BIO-13.

## 4. Process success is not a scientific acceptance gate

`run_lpj_guess.py` treats a zero MPI return code as success.  During merging it
uses whichever rank pieces exist, skips an output entirely if no piece exists,
and records row counts without checking them against expected ranks, cells and
years.  Downstream tools can therefore consume a spatially incomplete run.

The retained output list also omits the stock `tot_runoff`, `nmass`, `npool` and
`nflux` products.  Existing C outputs are sufficient to start a carbon check,
but the current artifact cannot close even the baseline nitrogen and water
budgets.  BIO-10 adds phosphorus diagnostics later; it needs a general
acceptance mechanism underneath it rather than a P-only exception.

The run harness must write a machine-readable pass/fail assessment covering
exact output/rank/cell/year coverage, finite and physically valid values,
end-window stability, and C, N and water residuals.  No albedo, soil or scoring
step should accept an unassessed or failed run.  This is BIO-14.

## 5. The declared soil-biosphere convergence criteria have no assessor

`pedogenesis.yaml` fixes a 2% relative tolerance on land-mean soil carbon, a 5%
limit on cells whose clay fraction moves by more than 0.02, and a maximum of six
iterations.  No code reads those keys.  The README accordingly records that the
loop has never been iterated to its own criteria.

A convergence step must compare successive provenance-compatible soil and LPJ
summaries, apply both criteria, write the measured residuals, and fail rather
than silently continuing beyond the maximum iteration.  It depends on BIO-12's
output statistic and BIO-14's accepted-run artifact.  This is BIO-15.

## 6. Modelled vegetation does not reach aerodynamic roughness

The coupling documentation says vegetation returns to ExoPlaSim through albedo,
forest fraction and roughness.  `build_surface_roughness.py` explicitly excludes
the `modelled` mode from its forest-fraction mapping and falls back to one canopy
roughness over all non-barren land.  Its pipeline step has no dependency on
`lpj_run`.

This leaves the intended aerodynamic feedback open.  Roughness must consume the
same temporally aggregated tree/grass cover and rootable fraction as albedo,
while retaining the existing orographic contribution and lake treatment.  The
pipeline must record the LPJ back-edge.  This is BIO-16.

## 7. Modelled cover repaints solved lakes

The albedo builder says lakes are applied last, but they are composited into the
substrate before the `mode == "modelled"` block.  That later block blends tree and
grass albedo over every ExoPlaSim land cell without consulting `lake_mask`, and
the modelled forest field likewise writes tree cover over all land.  The
non-modelled forest branch correctly excludes lakes, demonstrating the intended
semantics.

Modelled albedo and forest fraction must apply the shared rootable fraction and
must leave persistent water at water albedo with zero forest.  A regression
fixture should cover a partial lake cell, a fully barren cell and an ordinary
rootable cell.  This is BIO-17.

## 8. Modelled vegetation discards the two-band spectral correction

SPEC-5 derived a K-star vegetation endmember and the two-band pair `[0.075,
0.225]`.  The modelled path instead uses scalar tree and grass albedos, after
which the writer copies the same broadband grid to fields 174, 175 and 176.  It
therefore reintroduces the assumption SPEC-5 removed: equal vegetation
reflectance on both sides of 0.75 micrometres.  The archived SPEC-5 result
already calls the tree and grass constants its residual.

Tree and grass need separately justified broadband levels and band pairs whose
stellar-flux-weighted recombination reproduces those levels.  Modelled mode must
then write distinct 175 and 176 fields and report the recombination residual.
This is BIO-18.

## 9. The pipeline cannot prove soil and LPJ used the same climate

The `soil` step declares a dependency on `bootstrap_climatology`, while
`build_soil.py` defaults to the configured baseline climatology and `lpj_driver`
depends on `baseline_climatology` plus that soil artifact.  The graph therefore
describes a bootstrap-derived soil feeding a baseline-driven biosphere even when
manual practice has rebuilt the soil from the baseline.  `check_consistency.py`
checks the soil report's terrain hash but not its climatology hash against the
climatology selected for LPJ.

The bootstrap soil and the iterative baseline soil need distinguishable graph
states, or an equivalent mode-bearing artifact contract.  Consistency checking
must reject a soil map whose recorded climatology does not match the LPJ driver's
climatology.  This is BIO-19.

## 10. The project adapter has no cheap regression boundary

The vendored project carries upstream unit tests, but there are no project-level
fixtures for the Vesper binary contract and its consumers.  The highest-risk
interfaces are consequently exercised only by a full LPJ run.

A no-simulation fixture suite should cover driver serialization and parsing,
coordinate bijection, one- and multi-year indexing, exact monthly precipitation
conservation, missing-rank/output rejection, shared end-window reduction, and
rootable/lake/two-band feedback compositing.  It must remain cheap enough to run
without launching LPJ-GUESS.  This is BIO-20.

## Recorded structural limits, not reopened tasks

BIO-4 already established that ExoPlaSim lacks stomatal, LAI, interception and
rooting-depth control of evaporation; enabling SIMBA does not supply those
mechanisms.  That is a limitation of the climate model rather than an unfinished
LPJ adapter feature, so this audit does not reopen it under a new id.

The 30-hour-day frost limitation is likewise already recorded in
`biosphere/README.md`: LPJ-GUESS receives the diurnal range but uses it only in
the disabled BVOC path.  Adding a daily-minimum mortality mechanism would be a
new ecological model, not completion of the current port.

Earth-analogue PFT limits and the static annual return of vegetation properties
remain interpretation boundaries.  BIO-12 will expose stochastic and patch
uncertainty, but seasonal vegetation-climate coupling would require ExoPlaSim to
accept time-varying boundary fields and is not implied by the present loop.
