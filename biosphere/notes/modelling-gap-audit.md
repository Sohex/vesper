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
paths call `prdaily` when that flag is true.  `VesperInput` never calls
`prdaily`.  As measured, `VesperInput::interpolate` called
`interp_monthly_totals_conserve`, which turned each monthly total into a smooth
sequence of positive mean-daily values; `VesperInput::integrate_year` now takes
each absolute day's duration-weighted mean over the forcing intervals covering
it, which stops manufacturing the smooth curve but does not produce wet days
either, because an interval mean rate spread over its own days rains a little
every day just as the smooth curve did.  Either way the instruction file says
wet days only while Vesper rains a little almost every day, and the missing
thing is a source with events in it rather than a better reconstruction.

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

## 8. The two-band surface albedo, and the two rules it establishes

**Settled 2026-08-24.**  Every material this project writes to codes 174, 175
and 176 now carries a band pair, `--mode modelled` included, and two rules come
out of getting there that bind anything else that writes one.

The blocker was never the canopy.  A modelled cell is
`cover*canopy + (1-cover)*substrate`, and the substrate had no split at all:
the rock table carries one broadband number per class.
`analysis/rock_albedo_bands.py` supplies `band2/band1` for every Orogen rock
class from ECOSTRESS directional-hemispherical reflectance under this star.  It
derives the SHAPE and never the level, because laboratory preparation moves a
level by a factor of 2.2 to 5.1 on the same rock and cancels from a ratio; the
bracket is measured rather than argued, since the library carries the same rock
at four preparations and the ratio's sensitivity to preparation can be read off
it directly.

**A band pair is anchored on the model's own band weights.**  `radmod.f90`
forms `zsolars(1)*dsalb(1) + zsolars(2)*dsalb(2)`, and `lib/stellar` reproduces
those weights including the `minwavel` cut and the band-edge interval.  Anchor
on anything else and the pair does not return the broadband level the
derivation set: normalising instead on the flux share of the 0.35-2.5
micrometre range leaf spectra are measured over, which is 0.399 against the
model's 0.382, left vegetated ground about 0.0025 too bright, one-signed, and
`analysis/vegetation_albedo.json` reports that as
`band_recombination_residual_naive_anchor`.  Anchoring on the model's weights
removes it by construction, and `build_surface_albedo.py` asserts
`z1*175 + z2*176 == 174` cell by cell before writing, at twice the `.sra`
format's own quantum.  The assertion is what catches a repaint that moves a
level without moving the material's ratio with it.

**Whether tree and grass need separate band ratios is not resolved, and the
test says so.**  The broadband correction takes the population leaf ratio,
because the classes differ there by 0.001 in albedo under a grass sample of
four spectra.  The band ratio is a much larger difference -- tree 3.15 against
grass 2.56, measured 2026-08-24 -- and the pre-fixed test was whether the grass
sample lies entirely below the trees' tenth percentile.  It misses, 2.713
against 2.710.  So both endmembers take the population ratio and carry the
per-class ratio as the other end of the bracket, and their pairs are distinct
because their levels are.  What settles it is more grass spectra rather than
more argument: the trees' own ratios run 1.78 to 4.92, so a class mean here is
only as good as its sample.  BIO-18's own note argues the same class dependence
from leaf transmittance, and the leaf library cannot settle either half of it,
because it carries reflectance and not transmittance.

BIO-18 and `world-36g` are closed.

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
