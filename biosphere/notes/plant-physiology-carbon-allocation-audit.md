# Plant physiology and carbon allocation audit

**Scope:** the active LPJ-GUESS-CNP path from absorbed radiation through
photosynthesis, stomatal coupling, autotrophic respiration, tissue allocation,
nutrient storage and reproduction. This is a static source and primary-literature
audit. No build, LPJ-GUESS execution or other material computation was run.

The neighbouring boundaries are intentionally not re-audited here. Atmospheric
pressure and humidity enter through BIO-23; photon supply through BIO-25;
hydraulic supply and groundwater through PLHY; gravity-dependent support through
GRAV-7; establishment and mortality through DEMO; and litter fate through SDEC.
This note identifies the physiology and allocation contracts those tasks must
meet.

## Active model path

The Vesper instruction file enables cohort vegetation, nitrogen limitation,
SLA and minimum leaf C:N calculated from leaf longevity, vegetation carbon debt,
and the fork's acclimated-respiration option. Phosphorus limitation and the
Walker joint N--P photosynthesis option remain off. The configured natural PFTs
are Earth C3/C4 functional types with fixed photosynthetic temperature limits,
leaf longevity, tissue turnover, wood density, allometry, leaf:root ratio,
respiration, nutrient storage and a 10% reproductive fraction.

Daily carbon enters through the Haxeltine--Prentice analytical light-use model in
`modules/canexch.cpp`. Leaf respiration is proportional to Rubisco capacity;
sapwood and root maintenance respiration scale with tissue N and temperature;
growth respiration is 25% of carbon left after maintenance. NPP accumulates
through the model orbit and `modules/growth.cpp::growth()` allocates it once at
the orbit boundary. Leaf:root allocation follows the most limiting of annual
water, N and P scalars, while trees additionally satisfy the pipe-model and
height--diameter allometries. A retrospective sapwood loan (`cmass_debt`) is the
only carbon buffer. Positive NPP loses a fixed 10% to reproduction before tissue
allocation; that carbon is subsequently reported as nitrogen-free reproductive
carbon returned to the atmosphere.

This is a coherent Earth DGVM baseline, but it is not yet a coherent physiology
contract for a 30-hour rotating, 181-Earth-day-orbit, high-gravity K-star world.

## Findings

### 1. The photosynthesis day is internally fixed at 24 hours

`modules/driver.cpp::daylengthinsoleet()` converts the illuminated angular
fraction to `24 * hh / pi`. The Haxeltine--Prentice implementation independently
uses 24 hours in the optimized Vmax expression, Rubisco-limited hourly rate and
daytime respiration fraction (`vmax*()`, `photosynthesis()`). Its subdaily
contract also says that `daylength` must be 24 hours to obtain daily units.

This is not merely the already-declared error in astronomical daylength. The
Vesper adapter deliberately retains 24 hours of absolute time per step so that
per-day turnover and respiration constants keep their calibration, but one such
step spans only 0.8 of the 30-hour rotation. Conversely, a complete light/dark
cycle lasts 30 hours and cannot be represented by changing an equinoctial
photoperiod from 12 to 15 hours while continuing to call a 24-hour radiation
integral a complete daily integral. The analytical model's optimization,
light-saturation and dark-respiration terms all assume that its integration
interval is one complete solar day.

Haxeltine and Prentice (1996) make the assumption explicit: total daily absorbed
PAR is integrated under a sinusoidal daylight curve and respiration covers
`24 - daylength`. Their derivation supports the algorithm on its own time base,
not the mixed time base presently supplied to it.

This is a correctness blocker. The solution needs an explicit absolute-time
physiology interval and rotation phase: integrate incident/APAR and light versus
dark respiration across each 24-hour forcing interval, carrying phase across
steps, or use a subdaily forcing path that does so. It must conserve the ExoPlaSim
radiative energy integral and must not rescale every biological daily rate to a
30-hour day. BIO-20 can test constant-light, equinox, solstice and phase-wrapping
fixtures without launching LPJ-GUESS.

### 2. Photon correction alone does not close the K-star physiology contract

`VESPER_FRADPAR` converts surface shortwave into energy inside the declared
0.40--0.75 micrometre photosystem window. BIO-25 correctly owns replacement of
the fixed 550 nm `CQ` energy-to-photon conversion. Downstream, however,
`canexch.h` retains fixed Earth C3/C4 quantum efficiencies (`ALPHA_C3`,
`ALPHA_C4`), curvature, leaf-respiration-to-Vmax ratios and the globally tuned
canopy scalar `ALPHAA_NLIM`.

Kiang, Marosvölgyi--van Gorkom and Lehmer justify the declared spectral window
and photon calculation, not invariant action efficiency at every wavelength.
Walker et al. (2014) also show that the balance between electron transport and
carboxylation carries a photoprotection trade-off. Red/far-red acclimation,
blue-light stomatal response, within-canopy spectral quality and dissipation of
excess excitation are therefore omitted model form, not parameters determined by
the existing spectrum calculation.

The central case should remain the declared transplanted Earth oxygenic
photosystem with spectrum-weighted photon count. Register a bounded action-
spectrum/photoprotection sensitivity around it and keep BIO-26's canopy scalar
separate. Do not invent alien pigments or tune quantum efficiency to the desired
productivity.

### 3. The enabled acclimated-respiration routine has no acclimation state

`respiration_acclimated()` documents its `Tacc` argument as a running average of
air or soil temperature and applies the Thum et al. response
`10^[-0.008 (Tacc - 10.15)]`. Its caller instead supplies the current day's air
temperature and current 25 cm soil temperature. No running-average state or
process-specific memory is stored. The same instantaneous temperature therefore
drives both the acute Lloyd--Taylor response and the nominal acclimation
multiplier.

Thum et al. (2019) explicitly separates instantaneous meteorology from
process-specific memory and acclimation to prevailing growth temperature.
Gifford (2003) distinguishes the short-term temperature response from an
acclimated response that can change within about a week, while Atkin et al.
(2014) shows that including acclimation materially changes large-scale carbon
exchange. The present call does not implement the cited model form.

Add explicit, restart-safe air and root-zone acclimation-temperature states with
a declared absolute-time response constant, initialize them deterministically,
and test step, constant and restart-continuation cases. Until then, the named
option should fail closed or the standard respiration path should remain the
declared baseline.

### 4. Autotrophic respiration and construction cost are over-compressed

The model relates leaf respiration to Vmax and root/sapwood maintenance to
tissue N. It uses one 25 cm soil temperature for all fine roots, treats all
sapwood as respiring, gives tissue P no respiration role, and assigns a uniform
25% growth-respiration cost to all positive carbon remaining after maintenance.
It does not distinguish live sapwood fraction, tissue age, root depth, tissue
biochemical composition, phloem transport, nutrient uptake/transformation or
defence and repair costs.

Gifford (2003) shows why growth and maintenance coefficients are operational
constructs whose values vary by organ, age and composition, and why nutrient
uptake and transport costs fall awkwardly between them. Atkin et al. (2014)
warns against treating respiration as an isolated tax rather than a source--sink
process; light inhibition also makes daytime leaf respiration different from
dark respiration. A more detailed mechanism is not uniquely supported, but the
fixed residual fraction is not uncertainty-free.

Retain the present formulation as the reduced baseline, then register an organ-
and-depth-resolved respiration/construction-cost bracket. It must use PLHY-3's
root layers, distinguish living from structural support tissue, and propagate
GRAV-7's added structural carbon through construction and maintenance without
simply multiplying every metabolic cost by gravity.

### 5. The active leaf economics and tissue stoichiometry are one Earth trait chain

With `ifcalcsla` and `ifcalccton` enabled, leaf longevity determines SLA and
minimum leaf C:N using regressions from Reich et al. (1992). That paper found a
co-dependent syndrome: short-lived leaves tend to have high SLA, leaf N and
mass-based photosynthesis, while long-lived leaves trade instantaneous return
for persistence. It does not establish those regression coefficients for an
alien flora or license varying each trait independently.

The CNP extension is less secure. `Pft::init_ctop_min()` derives leaf C:P from
SLA; root and sapwood ratios then reuse Friend's N multipliers. The source itself
labels the root/sapwood C:N range as "picked out thin air", asks whether using
the same range for P "MAKES SENSE?", and gives P the same unsupported narrow
range. Fixed 50% N and P resorption and PFT storage fractions complete the live-
plant stoichiometric assumptions. Dantas de Paula et al. (2025) documents this
extension but does not make the inherited ratios planetary constants.

Create one provenance-stamped live-plant trait registry that preserves covariance
among leaf longevity, SLA, N, P, photosynthetic capacity, respiration, root
strategy, wood density, turnover, resorption and storage. Use coherent strategy
ensembles or named endmembers, not independent perturbations. SDEC-5 should
consume the resulting tissue-specific litter chemistry and GRAV-7/PLHY should
consume the same wood/root strategies.

### 6. The two phosphorus--photosynthesis options have different empirical domains

The disabled default P path uses a linear relation labelled to Hidaka and
Kitayama (2013), with a positive photosynthetic-capacity intercept even at zero
active P, and combines separately calculated N- and P-limited capacities by a
minimum. The alternative Walker relation models Vcmax jointly from leaf N and P
but is disabled.

Hidaka and Kitayama measured ten tropical montane tree species and showed that
P-use efficiency depends on allocation among foliar P fractions, particularly
metabolic versus lipid P. Total P alone is therefore an incomplete state, and
their measurements do not establish a universal Vcmax slope for every PFT.
Walker et al. (2014) compiled 356 species from 24 studies and found that P changes
the sensitivity of Vcmax to N; it is a stronger basis for a joint model, but its
authors still call for broader temperate and boreal coverage and identify
mesophyll-conductance and temperature uncertainties.

Before BIO-10 activates P limitation, reproduce both equations and their units in
no-simulation fixtures, establish domain guards including zero and low P, trace
the active/non-active P pool and `P0` provenance, and pre-register the independent
minimum versus joint Walker formulation as a structural sensitivity. The choice
must not be made by which produces the registered BIO-9 answer.

### 7. Allocation combines one-orbit memory with incomplete resource feedback

Carbon is accumulated daily but allocated once per 181-day orbit. The leaf:root
ratio is `min(mean water scalar, N scalar, P scalar) * ltor_max`, after which the
tree pipe model and height--diameter allometry determine tissue increments. The
water-supply source contains an explicit TODO to make supply depend on root
allocation, so the model can allocate more carbon below ground in response to
drought without gaining the corresponding water-acquisition capacity in the
same process representation.

Poorter et al. (2012) finds real allocation plasticity, strongest under nutrient
limitation and severe drought, but also shows that ontogeny and morphology often
change more than biomass fractions. Franklin et al. (2012) shows that fixed
fractions, allometry, functional balance and optimization are competing model
classes and recommends comparing alternative hypotheses rather than declaring a
fit by one flexible model. The current minimum-scalar rule is one such
functional-balance hypothesis.

Extend BIO-22's time registry to the allocation event, `wscal_mean`, carbon-debt
repayment and nutrient-storage turnover. Then register fixed/current functional
balance versus a bounded source--sink allocation alternative. Root investment
must feed PLHY-3 water and nutrient access; support investment must consume
GRAV-7's gravity-aware structural cost; and allocation should not be adjusted
twice for the same limitation.

### 8. Carbon debt is not nonstructural carbon storage

With `ifcdebt 1`, a tree may retrospectively borrow enough sapwood carbon to
satisfy allocation constraints, bounded by 80% of the deficit and 20% of
unencumbered sapwood, then repays 20% of the debt at each orbit allocation. This
is a numerical loan attached to structural carbon. There is no positive sugar or
starch pool, diel carbohydrate balance, reserve draw for leaf-out, storage sink,
source--sink growth limitation, osmotic reserve floor or disturbance recovery.

Dietze et al. (2014) establishes that nonstructural carbon supports growth,
respiration, transport, osmoregulation, cold tolerance and recovery, that plants
can access old reserves under stress, and that active, quasi-active and passive
storage may all occur. Thum et al. (2019) gives a relevant reduced architecture:
separate labile and reserve pools, source and sink limitation, and process-
specific memory. Neither supplies a unique Vesper parameterization.

Keep debt only as a documented numerical baseline. Add an explicit C pool bracket
with mass-conserving labile/reserve stocks, minimum functional reserve,
phenological pull, maintenance priority and storage push. Couple its depletion to
PLHY drought and DEMO mortality only through pre-registered rules; do not infer
carbon starvation from low structural growth alone.

### 9. Reproduction is carbon-only and routed directly to the atmosphere

Every positive annual NPP value loses `reprfrac = 0.1` before tissue allocation;
negative NPP has no reproductive cost. The stand-level carbon helps scale
establishment, but the patch-level amount is immediately accumulated as
"nitrogen-free reproduction litter" and emitted through `REPRC`. There are no
N or P construction costs, propagule/fruit pool, seed bank, reproductive
turnover, litter transfer or consumer pathway.

The carbon ledger can close because `REPRC` is an external flux, but the C--N--P
plant economy does not: reproduction competes for carbon without competing for
the nutrients required to build it. It can therefore bias tissue allocation and
nutrient limitation before the demography code ever uses the reproductive
signal.

Make reproductive allocation explicitly C--N--P and mass conserving. Treat the
fixed 10% as a baseline bracket, route propagule and accessory tissue to named
stocks and fates, and make DEMO establishment consume the same pool. If a
fraction represents metabolic construction cost or consumption, report it under
that name rather than calling it nutrient-free litter.

### 10. Accepted outputs cannot distinguish the competing explanations

The run harness retains aggregate GPP, NPP, biomass, LAI and cover, but not the
physiological state needed to tell whether similar vegetation arose through
photon supply, biochemical capacity, stomatal limitation, respiration,
allocation, storage or reproduction. Aggregate productivity could therefore
look plausible while compensating errors cancel.

Extend BIO-14 with absorbed energy and photons, day/night interval, Vcmax and its
N/P limits, stomatal/water limitation, leaf/root/sapwood maintenance and growth
respiration, organ C--N--P increments and stocks, resorption/storage fluxes,
carbon debt or NSC, reproductive allocation and source--sink limitation. Require
daily and orbit C--N--P closure, and make the accepted equilibrium window retain
both seasonal phase and allocation-event diagnostics.

### 11. Fine-root and sapwood C:N and C:P windows are anchored on the tissue mean, so the applied proportion is the sourced one

`Pft::init_cton_limits()` and `Pft::init_ctop_limits()` set the fine-root and
sapwood window MEAN from the leaf window mean and build the window outward from
it:

    <tissue>_avr = <leaf>_avr * frac
    <tissue>_max = <tissue>_avr * (1 + m) / (2 * m)
    <tissue>_min = m * <tissue>_max

with `m` = `frac_maxtomin` / `PFRAC_MAXTOMIN` = 0.9. That is the same inversion
`Pft::init_ctop_min()` performs on the leaf window and for the same reason:
`avg_cton()` and `avg_ctop()` are harmonic means, so a mean has to be inverted
through them to come back out. `m` therefore sets the WIDTH of the tight window
and nothing else, and the four proportion constants set where its mean sits.

What reaches the simulated plant is a ratio of window means. `canexch.cpp:1402`
and `:1406` compute fine-root and sapwood nitrogen demand against
`cton_leaf_opt * <tissue>_avr / cton_leaf_avr`, and `:1574` and `:1578` do the
same for phosphorus. The applied proportion is therefore exactly the declared
constant, for any leaf window width `f` and any `m`.

| element | tissue | phenology | constant | declared | applied | applied / declared | max-anchored applied | source quantity | source value |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| N | sapwood | non-crop | `frac_leaftosap`, `guess.h:2417` | 6.9 | 6.900 | 1.000 | 12.355 | Friend 1/X_C:N(f/p) | 6.897 |
| N | fine root | non-crop | `frac_leaftoroot`, `guess.h:2414` | 1.16 | 1.160 | 1.000 | 2.077 | Friend 1/X_C:N(f/r) | 1.163 |
| P | sapwood | non-crop | `PFRAC_LEAFTOSAP`, `guess.h:302` | 6.9 | 6.900 | 1.000 | 15.296 | none; the constant is nitrogen's | -- |
| P | fine root | non-crop | `PFRAC_LEAFTOROOT`, `guess.h:261` | 1.16 | 1.160 | 1.000 | 2.571 | Yuan contrast against the applied N root proportion | 1.0 |
| N | sapwood | CROPGREEN | `frac_leaftosap` at f = 5.0 | 6.9 | 6.900 | 1.000 | 19.611 | Friend 1/X_C:N(f/p) | 6.897 |
| N | fine root | CROPGREEN | `frac_leaftoroot` at f = 5.0 | 1.16 | 1.160 | 1.000 | 3.297 | Friend 1/X_C:N(f/r) | 1.163 |
| P | sapwood | CROPGREEN | `PFRAC_LEAFTOSAP` at f = 7.77 | 6.9 | 6.900 | 1.000 | 28.664 | none | -- |
| P | fine root | CROPGREEN | `PFRAC_LEAFTOROOT` at f = 7.77 | 1.16 | 1.160 | 1.000 | 4.819 | Yuan contrast | 1.0 |

That is the complete extent of the construction. The crop stem pair is the only
other tissue in the same demand construction and it never needed the repair:
`cton_stem_avr` and `ctop_stem_avr` are set directly as means at
`guess.h:2457-2458` and `:2515-2516` rather than derived from a `_max`, and
`guess.cpp:2034` and `:2045` then use them in exactly the form the four
constants are read in. That is also the internal corroboration of which quantity
this construction wants. The leaf window widths are not part of it: `cton_leaf_min`
is a genuine endpoint from Reich et al. (1992) and `ctop_leaf_min` is already
inverted so that `ctop_leaf_avr` reproduces its regression's average.

**What the max-anchored form did.** It set `<tissue>_max = <leaf>_max * frac` and
`<tissue>_min = m * <tissue>_max`, so the constants were ratios of window
ENDPOINTS while the model applied a ratio of window MEANS. The leaf mean sits
well below the leaf maximum, which put a fixed factor between the declared
constant and the applied proportion:

    applied / declared = [2m / (1 + m)] * (1 + f) / 2

That is 1.7905 on the non-crop nitrogen window, 2.2168 on the non-crop
phosphorus window, and 2.8421 and 4.1542 on the CROPGREEN arm, and it is what
produced the max-anchored column above. The check it failed is the acceptance
test for this finding: after `init_cton_limits()` and `init_ctop_limits()`,
`<tissue>_avr / <leaf>_avr` must equal the constant the source measured, for
every pair in both phenology arms.

**The applied number is the one that has to equal the source, and the endpoints
are the free parameters.** Friend, Stevens, Knox and Cannell (1997) allocate
nitrogen so as to hold the relative C:N ratios of foliage, sapwood plus bark and
fine roots FIXED. From their Eqs. (42) to (44), the sapwood and foliage nitrogen
fractions are `d = X_C:N(f/p) * c * Cp/Cf` and `c`, so
`(Cp/Np)/(Cf/Nf) = 1/X_C:N(f/p) = 6.897` at every nitrogen status, with no window
anywhere in the formulation. `X_C:N(f/p)` is a mean across nine species of
measured tissue ratios. The construction
`cton_leaf_opt * cton_sap_avr / cton_leaf_avr` is that fixed relative ratio
applied to a simulated individual's own current leaf ratio, so
`cton_sap_avr / cton_leaf_avr` is precisely the quantity `1/X_C:N(f/p)` names.
Heineman, Turner and Dalling (2016) regress log species-mean wood on log
species-mean leaf concentration, so its proportion is likewise between tissue
means. Friend's own sensitivity analysis, Table 8 p. 280, ranks `X_C:N(f/r)`
ninth of 76 parameters at index 0.665 and `X_C:N(f/p)` twenty-sixth at 0.337, so
a factor of 1.79 on them is not a detail.

**What it moved, and by how much.** All arithmetic below is closed form on the
model source and the two papers. Nothing here is execution-verified: LPJ-GUESS
does not build on this tree and no run was made.

The repair divides every fine-root and sapwood C:N by 1.7905, on both tissues by
the same factor, and leaves the leaf window untouched. Tissue nitrogen demand per
unit tissue carbon therefore rises by 1.7905 in both. Whole-plant structural
nitrogen requirement rises less, because leaf demand does not move: at leaf C:N 30
and carbon 0.3 leaf, 0.3 fine root, the requirement rises by a factor of 1.257
with no sapwood, 1.371 at 1.5 kgC/m2 sapwood, 1.446 at 3.0 and 1.536 at 6.0.
Those carbon masses are stated as an illustration of composition dependence, not
as a prediction of this world's. On the phosphorus side the same factor is
2.2168, and it moves nothing that limits simulated growth: `ifplim 0`, and
`parameters.cpp` refuses `ifplim 1` on separate grounds.

Maintenance respiration depends on which respiration path is configured, and the
two differ completely. `respiration()` takes `respcoeff`, and
`guess.h:2450-2456` divides `respcoeff` by
`cton_root/(cton_root_avr + cton_root_min) + cton_sap/(cton_sap_avr + cton_sap_min)`.
Under the mean anchor `<tissue>_min` is `(1 + m) / 2` of `<tissue>_avr`, so that
sum is 1.95 times each tissue's own mean, exactly as it was under the max anchor.
Scaling both windows by the same factor therefore scales that sum inversely and
leaves `respcoeff / cton` unchanged, so on that path sapwood and fine-root
maintenance respiration is invariant by construction. The configured path is not
that one. `build_vesper_pfts.py` defaults its source instruction file to
`global.ins`, which sets `acclimated_respiration 1`, and
`respiration_acclimated()` substitutes a temperature function for `respcoeff` and
never reads it (finding 3), so on the configured path there is no compensation
and maintenance respiration scales as `1/cton_sap` and `1/cton_root` directly.
`global_p.ins` sets the switch to 0, so moving the run to the C-N-P instruction
file would restore the `respcoeff` compensation and remove this effect entirely.

The size of that on the configured path is bounded rather than known.
`Individual::cton_sap()` returns
`max(cton_sap_avr * cton_leaf_min / cton_leaf_avr, cmass_sap / nmass_sap)`, and
the floor moves by the same 1.7905, so a simulated individual whose nitrogen
demand is met sees the full factor while one that is nitrogen-starved sees less.
In the nitrogen-replete limit, with `resp_growth` a fixed quarter of what is left
after maintenance, NPP is `0.75 * (A - Rm)` and multiplying `Rm` by 1.7905 leaves
NPP at 0.80, 0.66 and 0.47 of its former value for `Rm/A` of 0.20, 0.30 and 0.40.
That ratio is BRACKETED 0.20 to 0.40 rather than measured; this project has no
accepted respiration diagnostic, which is finding 10.

So the C-N configuration moves by a first-order amount in a known direction on
every quantity the biosphere returns, and only a run settles the magnitude. Every
number produced under the max anchor is worthless in rule 7's sense rather than
stale, and the biosphere's re-commissioning is what replaces them. The canonical
climatology lineage does not exist yet, so a defect in the model was never weighed
against the cost of the output it invalidates.

**`PFRAC_LEAFTOROOT` is now delivered, and is decoupled from `PFRAC_MINTOMAX`.**
That constant's derivation is a CONTRAST and not a level: a C:P proportion
divided by its C:N counterpart is the fine-root-to-leaf N:P ratio, and Yuan, Chen
and Reich (2011) cannot separate live-root N:P from leaf N:P, so the contrast is
1.0 and the constant keeps nitrogen's 1.16. Under the mean anchor the applied
contrast is that ratio of constants, exactly 1.0, whatever widths the two leaf
windows carry. Under the max anchor it was `(1 + PFRAC_MINTOMAX) / (1 + 2.78)` =
4.68 / 3.78 = 1.238, because BIO-34 widened the leaf C:P window to 3.68 while the
leaf C:N window stayed at 2.78, and the point estimates Yuan's test could not
separate give a contrast of 0.879 to 1.164, so the derivation was not reaching
the simulated plant. The coupling BIO-34 introduced is removed rather than
re-tuned around: the alternative that keeps the max anchor needs
1.16 * 3.78 / 4.68 = 0.937 here and needs redoing whenever either leaf window
moves.

**What the repair does not settle.** `PFRAC_LEAFTOSAP` remains refused on its own
grounds. Its 6.9 is a NITROGEN ratio for foliage against bark plus sapwood over
nine temperate species, and Heineman et al. (2016) rejects the proportional FORM
for phosphorus while failing to reject it for nitrogen. The quantity question is
answered and the element and form questions are not: the model now applies
exactly 6.9 where a forced scalar from that one paired dataset would be BRACKETED
10.1 to 15.5, so the applied value sits below the bracket instead of at the top
of it by cancellation. `ifplim 1` keeps refusing and keeps naming it. BIO-34.

## Recommended order

1. Correct the mixed rotation/physiology time base and respiration-acclimation
   state before measuring any response.
2. Re-commission the biosphere on the mean-anchored tissue windows of finding 11.
   They decide what fine-root and sapwood nutrient demand and, on the configured
   acclimated path, maintenance respiration are, so no C-N number produced under
   the max anchor survives them.
3. Finish BIO-23, BIO-25 and the live-plant trait registry so pressure, photons
   and coherent Earth-derived priors are explicit.
4. Close the P-photosynthesis domain before BIO-10, and connect allocation to the
   same PLHY/GRAV/SDEC traits.
5. Add accepted diagnostics, then compare reduced model-form brackets for
   respiration, allocation, NSC and reproduction. No item here authorizes an
   LPJ-GUESS run.

## Sources read

- Haxeltine & Prentice (1996), *A general model for the light-use efficiency of
  primary production*, `10.2307/2390165`.
- Reich, Walters & Ellsworth (1992), *Leaf Life-Span in Relation to Leaf, Plant,
  and Stand Characteristics among Diverse Ecosystems*, `10.2307/2937116`.
- Gifford (2003), *Plant respiration in productivity models: conceptualisation,
  representation and issues for global terrestrial carbon-cycle research*,
  `10.1071/FP02083`.
- Franklin et al. (2012), *Modeling carbon allocation in trees: a search for
  principles*, `10.1093/treephys/tpr138`.
- Poorter et al. (2012), *Biomass allocation to leaves, stems and roots:
  meta-analyses of interspecific variation and environmental control*,
  `10.1111/j.1469-8137.2011.03952.x`.
- Hidaka & Kitayama (2013), *Relationship between photosynthetic phosphorus-use
  efficiency and foliar phosphorus fractions in tropical tree species*,
  `10.1002/ece3.861`.
- Atkin et al. (2014), *Improving representation of leaf respiration in
  large-scale predictive climate-vegetation models*, `10.1111/nph.12686`.
- Dietze et al. (2014), *Nonstructural Carbon in Woody Plants*,
  `10.1146/annurev-arplant-050213-040054`.
- Walker et al. (2014), *The relationship of leaf photosynthetic traits -- Vcmax
  and Jmax -- to leaf nitrogen, leaf phosphorus, and specific leaf area: a
  meta-analysis and modeling study*, `10.1002/ece3.1173`.
- Thum et al. (2019), *A new model of the coupled carbon, nitrogen, and
  phosphorus cycles in the terrestrial biosphere (QUINCY v1.0; revision 1996)*,
  `10.5194/gmd-12-4781-2019`.
- Friend, Stevens, Knox & Cannell (1997), *A process-based, terrestrial biosphere
  model of ecosystem dynamics (Hybrid v3.0)*, `10.1016/S0304-3800(96)00034-8`.
  Table 4 p. 254, Eqs. (42) to (44) and text p. 274 are the source of both
  nitrogen tissue proportions; Table 8 p. 280 is their sensitivity ranking.
- Heineman, Turner & Dalling (2016), *Variation in wood nutrients along a
  tropical soil fertility gradient*, `10.1111/nph.13904`. Table 3 is the wood
  against leaf scaling that rejects a fixed phosphorus proportion and fails to
  reject one for nitrogen.
- Yuan, Chen & Reich (2011), *Global-scale latitudinal patterns of plant
  fine-root nitrogen and phosphorus*, `10.1038/ncomms1346`.

The Dantas de Paula et al. (2025) CNP-fork paper and the previously read K-star,
hydraulics, groundwater, gravity, demography and soil papers listed in
`references/INDEX.md` were also checked at the interfaces above.
