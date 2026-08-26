# Tuned values across the tree

Audited 2026-08-25. The question: which constants in this project are held in
place by nothing except that they made a comparison come out right?

A TUNED VALUE, for this audit, is a constant whose only justification is that it
was fitted, calibrated or adjusted until a comparison landed. Its derivation
cannot be inspected, cannot be checked, and cannot be transferred to a different
planet or star. That is a narrower class than either of the two audits beside
it, and the three do not overlap:

| class | what is missing | the audit |
| --- | --- | --- |
| tuned | any derivation at all; the number IS the residual of a fit | this one |
| opaque | the derivation exists somewhere a reader of the code cannot reach | `opaque-constants.md` |
| implicit-Earth | the derivation is sound and is for the wrong planet | `ocean-tier-implicit-earth.md`, `inherited-earth-constants.md` |

An Earth measurement has provenance. A published fit has provenance: the paper
can be opened and the fit's domain checked. A tuning has none, and that is the
whole of the difference. A constant is NOT tuned here if it is a measured
material property with a citation, a design decision this project took
deliberately and recorded, a numerical parameter following from a stability or
resolution argument, or a threshold fixed in advance as a criterion.

## What a change to one of these costs, which is nothing

**The canonical climatology lineage does not exist.** Nothing has been
commissioned, so no run, climatology or build is protected, and replacing a
tuned value is not a re-commissioning cost. It is simply what the next cycle
implements. Counting the invalidated output as a price invents a sunk cost and
biases every decision in the table below toward keeping the tuned value, which
is the one outcome the rule exists to prevent. CLAUDE.md rule 7 already says
this; it is restated here because these are exactly the decisions where it gets
forgotten.

## The queue

Ranked by reach. The first six are in different units and the ranking across
them is a judgement rather than an arithmetic; within each the magnitude is
stated in whatever unit the value actually lives in, and in kelvin through
`lib/sensitivity.py` wherever the value is radiative.

---

### 1. `baseline_flux_earth = 0.945`, and the cap that was solved backwards

`config/planet.yaml:135`. The planet's insolation as a fraction of Earth's. It
fixes the semi-major axis and therefore the orbital period, which is COMPILED
INTO the biosphere through `vesper.h`, and it is read by every climate run, by
`lib/sensitivity.py`, by `scripts/error_budget.py` and by
`hydrography/scripts/export_carve_list.py`.

**The evidence is the project's own pipeline file.** `config/pipeline.yaml`'s
`design_flux` gate, on the artifact that was supposed to justify the number:

> the purged artifact inferred it at 0.0642, which reproduced 0.945 from an
> unconstrained winner of 0.8725, and a threshold fitted to its own answer is
> what CLIM-30 exists to undo

`inherited-earth-constants.md` finding 5 states the same thing and adds that the
original derivation has no script, no analysis product, no pipeline row and no
recorded threshold, and does not name the runs it was projected across; those
runs are deleted, so it is not reproducible even in principle. The file itself
carries the number as PROVISIONAL.

This is the shape in its purest form: a free parameter, the extreme-cold land
fraction cap, solved backwards so that a criterion reproduces a number that was
already there.

**Magnitude.** 0.945 against the unconstrained 0.8725 is a flux ratio of 0.0725,
98.7 W/m2 of TOA insolation and about 17 W/m2 of absorbed flux at a planetary
albedo near 0.3. At `SLOPE_K_PER_FLUX_RATIO` = 202.0 that is **14.6 K of
global-mean surface temperature**. The largest single entry in
`scripts/error_budget.py` is about 1.2 W/m2, so this one value is an order of
magnitude larger than the whole of the rest of the budget, and it is the only
member of it that is not inspectable.

**Disposition: REPLACEABLE, BLOCKED.** The replacement is a pipeline step that
already exists, `exoplasim/scripts/derive_design_flux.py`. It is blocked on two
things and only two: a cap declared in advance, which is `clim-30`, and a
baseline climatology, which this build does not have. No paper is needed.

---

### 2. `ALPHAA_NLIM = 0.6`, fitted to Earth's realised carbon cycle

`vendor/lpj-guess/modules/canexch.h:75`, selected at `canexch.cpp:503` under
`ifnlim 1`, which is the live branch. It scales PAR absorption from leaf to
plant projective area for every simulated individual on every day, so
everything downstream of gross primary production inherits it.

**The source says what it is, in the source.** `canexch.h:75`:

> Same as ALPHAA above but chosen to give pools and flux values that agree with
> published estimates when Nitrogen limitation is switched on.

and on `ALPHAA` at `:70`, "value chosen to give global carbon pool and flux
values that agree with published estimates".

**Disposition: IRREDUCIBLE.** There is no sourced form, and the reason is
structural rather than a gap in the literature. It is not a property of a leaf
or a canopy; it is an ecosystem-level scalar whose target is Earth's global
carbon budget, and this world has no measured carbon budget to put in its place.
`biosphere/notes/implicit-earth-assumptions.md:148-156` reaches the same verdict.
What can be done about an irreducible tuning is bound its influence, which is
`bio-26`.

**Magnitude.** Dimensionless, and the four declared variants span 0.5 to 0.9
with nothing bounding it beyond [0, 1]. A move of 0.1 is 17 per cent on absorbed
PAR at the plant level.

---

### 3. The regolith depth level: `erosion_weight` and `erosion_reference_relief_m`

`pedology/config/pedogenesis.yaml:557` and `:560`, consumed at
`pedology/scripts/build_soil.py:313-318`:

    erosion = erodibility * max(relief_m,0)/erosion_reference_relief_m * moisture
    depth   = maximum_depth_m * production / (production + erosion_weight * erosion)

**The file confesses the first and says nothing about the second.**
`pedogenesis.yaml:545-549`: "the one free scale parameter in the depth model.
There is no measurement of this world's soils to fit, so it is calibrated
against a declared Earth-analogue target". `pedology/README.md:575` repeats it.

**The two are exactly degenerate and only one is watched.** They enter the depth
solely through the ratio `erosion_weight / erosion_reference_relief_m`, so
halving 200 is identical to doubling 1.5. One degree of freedom is carried in
two keys, of which one is labelled free and the other carries a comment
describing a different quantity ("Normalised so Earth-like mean land conditions
give E near 1") with no source, no bracket and no note anywhere in the tree.

**What breaks.** Regolith depth times `volumetric_capacity_by_texture` is the
land-mean plant-available water capacity, and that field IS the model's `dwmax`
bucket whose overflow is the runoff, which drives the weathering law that
produced the depth. `pedology/README.md:686` says so. This is the loop-closing
field, and its level rests on a calibration.

**Magnitude.** At the land mean the erosion term carries 0.792 of the
denominator, so a 10 per cent move in either key moves land-mean regolith depth
and land-mean PAWC by 7.9 per cent in the opposite direction.

**Disposition: REPLACEABLE, SOURCE NEEDED, and one of the two papers is
already read.** Heimsath et al. (1997), `10.1038/41056`, gives the production
function and is read. Portenga and Bierman (2011), *Understanding Earth's
eroding surface with 10Be*, GSA Today 21(8), 4-10, `10.1130/G111A.1`, gives
global cosmogenic denudation by lithology and is HELD, which by
`references/INDEX.md`'s own rule is an open exposure rather than a resource: its
numbers have only ever been used secondhand here. An absolute
production-to-erosion ratio is constructible from the pair rather than fitted,
and reading the second is the first step.

---

### 4. `maximum_depth_m = 5.0`, the asymptote nobody labelled

`pedology/config/pedogenesis.yaml:543`. The whole comment is "Asymptote, not a
clip: the depth to which weathering can reach at all." The block above argues at
length for the saturating FORM over Heimsath's exponential and concedes that it
is a parameterisation rather than a derivation; neither says where 5 m comes
from. No source, no bracket.

**Magnitude.** A pure prefactor: a 10 per cent move changes every cell's
regolith depth and therefore land-mean PAWC by exactly 10 per cent. Per unit
fractional change it is a larger lever than `erosion_weight`, and unlike
`erosion_weight` it is not labelled free anywhere.

**Disposition: REPLACEABLE, SOURCE NEEDED.** Same pair of papers as row 3; a
weathering-front depth is what a production function integrates to.

---

### 5. The carve incision coefficient, and the size floor that sets its target

`hydrography/scripts/export_carve_list.py:319` (`calibrate_coefficient`), with
the target at `:142-144`: `EARTH_STANDING_BASINS = 15`,
`EARTH_BAND_LAND_MKM2 = 78.9`, `CALIBRATION_LATITUDE = 35.0`.

`C` is the leading coefficient in `cut = C * erodibility * slope^n * Q^m`, and
therefore decides `retain` for every basin in the catalogue, and therefore what
Orogen carves.

**Admitted twice.** `export_carve_list.py:324-330`: "**Calibrated here rather
than declared, because the target moves.** The coefficient is whatever leaves as
many overflowing-but-still-standing basins per unit land as Earth has".
`hydrography/notes/retain-fraction.md:47-50`: "It is calibrated against Earth
rather than declared, and the reason is that the relaxation window is undefined
rather than unmeasured."

**This is the borderline case and it lands on the tuned side.** It is re-solved
against this world's own basin population every run, which is better than a
carried-over literal. What it is solved AGAINST is an Earth observable
transplanted whole, and its only justification is that the comparison comes out.

**The larger lever is the size floor, and it is not bracketed.**
`export_carve_list.py:137-139` says the floor "decides the answer: counting from
10 km2 instead gives a density 33 times higher and reads as 'carve almost
nothing'". The floor is 1,000 km2, defended as what the mesh can resolve, but
the mesh cell is 285 km2 and the median overflowing basin at spill is 7,338 km2,
so 1,000 sits between the two and is chosen rather than derived.

**Magnitude.** `retain-fraction.md:310-312`: across the Poisson bracket on
Earth's 15 the marginal class moves by a factor of 1.9. A single verdict change
moved the solved `C` from 161 to 235.9, outside its own published 128 to 218
bracket, which is the strongest evidence that this number is a residual rather
than a property.

**Disposition: IRREDUCIBLE for `C`; REPLACEABLE NOW for the floor.** `C` absorbs
a relaxation window that `docs/src/reference/no-time-axis.md` establishes cannot
exist here, so no dataset can supply it. The size floor is a different matter: a
sensitivity across the floor, reported beside the verdict, converts a hidden
33-fold lever into a declared bracket, and that costs one re-run of a script.

---

### 6. `eddy_wind_m_s = 5.94`, measured on a run that was deleted

`config/planet.yaml:731`. The single free quantity in the hyperdiffusion rule
`tau_vorticity = pi * radius / (NTRU * eddy_wind)`, dividing into every entry of
the `timescales_days` table at `:747-752` for all five rungs and four fields,
which reach the model as TDISSD, TDISSZ, TDISST and TDISSQ.

**The comment is a measurement claim**: "Measured, not assumed: eddy RMS wind
... from run_4182235e9781 orbit 39." That run is not on this tree. It survives
as a `restart_from_run` string in seven manifests and one INDEX row in
`archive/runs/_bulk_2026-08-24/INDEX_AT_DELETION.json`; the archive's README
records that those directories were removed by an `rm -rf` whose argument came
from an unguarded command substitution that expanded to nothing. No manifest, no
series, no NetCDF.

So this fails the audit's criterion on "cannot be inspected, cannot be checked"
rather than on "fitted". The distinction matters for the fix. Note also that the
block's own corroboration is on the RULE and not on the value: ECHAM5's
empirically tuned table implies a constant eddy wind of 12.6 to 17.5 m/s, and
5.94 is a factor of two to three below that band.

**Magnitude.** Linear: a factor of two in `eddy_wind` halves every damping
timescale in the table.

**Disposition: REPLACEABLE NOW ONCE A BASELINE EXISTS.** No paper is needed. The
recipe is fully stated in the config block, the quantity is well defined, and it
needs one converged climatology, which loop A produces anyway.

---

### 7. `acllwr = 0.100`, the longwave cloud absorption coefficient. FIXED

`vendor/exoplasim/exoplasim/plasim/src/radmod.f90:77`, used at `:3361`.

It carried no comment, no unit and no source, standing in an expression the
fork's own reading identifies as Kiehl et al. (1998) Eqs. 12-13 written out.
Eq. 14 gives `k_l` = 0.090361 m2/g for liquid cloud water and Eq. 15's ice
coefficient is smaller than `k_l` at every effective radius, so no phase-weighted
`k_abs` can exceed `k_l`; this model has no phase split, so the liquid-only limit
is both the consistent value and the ceiling. **Upstream's 0.100 was above that
ceiling**, which is the check that could have failed and did not: it cannot be
read as a rounded weighted `k_abs`.

**Magnitude.** Confined to thin cloud. On the layer cloud water paths
`cloud-water-reference.md` tabulates, the modelled longwave cloud emissivity
moves by 0.0175 at 17.6 g/m2, 0.0371 at 6.9 g/m2 and 0.0249 at 2.2 g/m2, peaking
at 0.0373 near 6.35 g/m2, and by nothing at all in the two lowest layers, which
are saturated at either coefficient. Those are the cold high layers, so the sign
on outgoing longwave is set by the modelled high cloud fraction, which only a
baseline can give.

**Disposition: REPLACEABLE NOW, and done.** The value is 0.090361 with the
derivation in-code. `acllwr` leaves PHYS-11's task list, which carried it because
it was untraced.

---

### 8. `th2oc = 0.024`, the only water vapour continuum in the model

`radmod.f90:72`, applied at `:3463` as
`zah2o = min(zah2o + (1 - exp(-th2oc*zsumwv)), 1)`.

The comment is "absorption coefficient h2o continuum (lwr)". No unit, no source,
no derivation anywhere in the tree. It is the fourth member of the per-truncation
tuning table `radmod.f90:976-989` describes upstream carrying, alongside `tswr1`,
`tswr2` and `tswr3`, selected through a flag that never fired at any truncation.
The three siblings are named as tunings in their own declarations; this one is
not, and it is the only one of the four this project has never re-weighted or
bracketed.

**What it does.** Sasamori's broadband absorptances have no window absorption at
all, so this term is the entire continuum contribution to the modelled longwave.
`config-rationale.md:776-778` reaches it, and only to establish that there is no
shortwave counterpart; it does not ask where 0.024 came from.

**Magnitude.** The whole term is `1 - exp(-0.024 w)` on the pressure-weighted
column path `w` in g/cm2, so setting `th2oc` to zero, which the surrounding
`if(th2oc > 0.)` supports, removes between 0.012 and 0.058 of broadband
absorptance over a column path of 0.5 to 2.5 g/cm2. Against a clear-sky
greenhouse trapping of order 150 W/m2 that brackets the term at roughly 2 to
9 W/m2, hence 2 to 7 K by `lib/sensitivity.py`. The bracket is wide because the
path is not measurable without a baseline climatology; what is not in doubt is
that the term is first-order and the coefficient is a residual.

**Disposition: REPLACEABLE, SOURCE NEEDED, and the source is the one this
project already knows it lacks.** `references/INDEX.md` records that Mlawer et
al. (2012), `10.1098/rsta.2011.0295`, read, "does not supply an evaluable
continuum": the coefficients ship as data with LBLRTM. So the route is a
correlated-k or line-by-line calculation with MT_CKD on this path, which is the
same bundle `exoplasim/notes/corrk-cross-check.md` says the project does not
have and which also blocks `h2o_sw_level`'s bracket. The two open together.

---

### 9. `EXOPLASIM_DZ0LAND_M = 2.0`, the anchor a derived field is rescaled onto

`exoplasim/scripts/build_surface_roughness.py:137`, used as the bisection target
at `:302-313`. `aeolian/config/dust.yaml:152-153` confirms it from the consumer
side: "its land mean is anchored to 2.0 m so the value ExoPlaSim was tuned
against does not move."

**The decision is inspectable and the anchor is not.** The builder derives a
per-cell roughness field from this world's own lithology and subgrid relief,
spanning 0.025 to 11.2 m, and then solves a free coefficient so that the
area-weighted land mean of that field equals `dz0land`. The argument for doing
so is recorded (`config-rationale.md:764-766`, "redistributes roughness without
moving the global value the model was tuned against") and is a legitimate design
decision. What it anchors ON is a compiled default whose only documentation in
the model source is the comment "roughness length land", and nothing in this
tree establishes that 2.0 m is anything but an unexamined number.

So the tuned value here is not the field and not the decision. It is the target,
and the step that discards a physically sourced land mean in favour of it.

**Magnitude.** `opaque-constants.md` finding 7 prices the neighbouring error at
8.4 times too much land exchange; the anchor itself sets the level of land
sensible heat and land evaporation, and therefore the carve criterion's
numerator. `analysis/spatial_reduction_gap.py:376-395` already exists to price
the ladder half of it: the solved coefficient rises by 1.98 from T21 to T170 to
hold the mean, so two rungs built with the defaults differ by their terrain AND
by their calibration.

**Disposition: REPLACEABLE NOW.** The field's own area-weighted land mean, taken
from the lithology roughness map with no rescaling, is a sourced number this
project already computes. Reporting it beside the anchored one, and running the
free arm, converts the anchor from an assumption into a bracket. What it needs
is a decision, not a paper.

---

### 10. `vdiff_lamm = 160.0` and the three Louis fives

`vendor/exoplasim/exoplasim/plasim/src/fluxmod.f90:47` and `:56-58`.

`vdiff_lamm` is the asymptotic mixing length, the single number that sets
free-tropospheric vertical diffusivity through `zmixm = lamm*k*z/(lamm + k*z)`
and, through `zlamh`, its heat counterpart. The declaration says what it is:
"160 m is ECHAM's and is a fit to Earth's free troposphere. It is not obviously
transferable, and it is left at ECHAM's value because nothing in this project has
measured a replacement." The three stability-function coefficients are cited to
ECHAM REPORT 218, which is "a report and not a derivation", and the declaration
says none of the three is separately justified there or here: "they are one
fitted set and are used as one."

**Disposition: IRREDUCIBLE as stated, REPLACEABLE in principle.** Louis (1979),
*A parametric model of vertical eddy fluxes in the atmosphere*, Boundary-Layer
Meteorology 17, 187-202, `10.1007/BF00117978`, is the primary and is not in
`references/`; fetching it would establish whether the three fives are Louis's
own fitted set or ECHAM's re-fit, which is a different question from whether
they transfer. The mixing
length is the harder half: 160 m is a length scale of Earth's free troposphere
and the corresponding scale here follows from the atmosphere's own scale height,
which is 0.766 of Earth's. That is a derivation this project can do, and doing it
is a physics claim rather than a transcription.

---

### 11. `frac_labile_carbon = 0.5`, in no register at all

`vendor/lpj-guess/data/ins/global_soiln.ins:6`, consumed at
`vendor/lpj-guess/modules/somdynam.cpp:1316`.

It sets the fraction of total microbial respiration declared labile. That pool
is split by water-filled pore space and becomes the Michaelis-Menten substrate
term that multiplies BOTH denitrification steps
(`vendor/lpj-guess/modules/ntransform.cpp:412`). It is the only substrate control
on denitrification in the model.

**The evidence is an absence, and it is a complete one.** The instruction file
gives it no comment. It is not in `biosphere/config/ntransform.yaml`'s
`instruction_parameters.values`, not in `parser_bounds`, not in the
`calibration.entries` block, and not in `somdynam.yaml`. A grep for it across
`biosphere/` returns nothing. `parameters.cpp:55` defaults it to 1.0 and the
instruction file halves that with no note.

This is the finding the biosphere sweep exists for: `ntransform.yaml` is an
honest register that names seven soil-nitrogen constants as `unsourced` and
refuses on them under `--strict`, and this one falls outside it.

**Magnitude.** Dimensionless on [0, 1]. A factor of two against the compiled
default, sitting inside a saturation, so between one and two times on the
denitrification substrate and hence on the mineral nitrogen the simulated plants
can reach.

**Disposition: REPLACEABLE, SOURCE NEEDED.** The labile fraction traces to Li et
al. (1992), which `ntransform.yaml:743` records as not held here: *A model of
nitrous oxide evolution from soil driven by rainfall events: 1. Model structure
and sensitivity*, J. Geophys. Res. 97(D9), 9759-9776, `10.1029/92JD00509`,
confirmed against Crossref. Fetching it settles this and the
`michaelis_menten_divisor` entry together. Registering it in `ntransform.yaml`
is the immediate move and does not wait on the paper.

---

### 12. `gamma = 0.01`, the precipitation re-evaporation fraction

`vendor/exoplasim/exoplasim/plasim/src/rainmod.f90:78`, applied at four sites
around `:2475`. The fraction of the sub-saturation deficit falling precipitation
evaporates per timestep, which is what sets how much precipitation reaches the
ground, which is P minus E over land directly.

The rung-keyed branch is gone and the declaration now says why, but it also says
the thing this audit is about: "with no derivation on either side". Upstream had
0.01 and 0.007 and neither carries one.

**Magnitude.** The two upstream values differ by 43 per cent, which is the only
scale anyone has ever put on it.

**Disposition: REPLACEABLE, SOURCE NEEDED.** A re-evaporation fraction per
timestep is not a published quantity; what is published is a fall-speed and a
ventilated evaporation rate for a drop-size distribution, from which a per-step
fraction follows given the layer depth and the step. Kessler (1969), *On the
Distribution and Continuity of Water Substance in Atmospheric Circulations*,
`10.1007/978-1-935704-36-2`, is the canonical form and is not in `references/`.
This is a derivation to do, not a number to copy.

---

### 13. `zcca = 0.245`, `zccb = 0.125`, and the `rcrit` floor of 0.85

`rainmod.f90:2086` and `:207`. The convective cloud fraction as
`zcca + zccb*log(convective rain rate)`, and the cell-mean relative humidity at
which stratiform cloud starts.

Both declarations state the defect and neither states a source. On the pair:
"The rate it is fitted against is a CELL MEAN, so the same simulated storm
spread over a smaller cell gives a larger rate and more cover: the cloud
fraction of a convecting region is rung dependent through this pair alone.
Anchored to T21." On `rcrit`: "it encodes an assumed distribution of humidity
WITHIN the cell: a smaller cell holds a narrower distribution and should start
later ... the horizontal assumption has no NLAT term anywhere. Anchored to T21."

**Magnitude.** `rcrit` enters as `((rh - rcrit)/(1 - rcrit))^2`, so at the eight
interior levels `1/(1-rcrit)^2` is 44 and a 0.01 move in the floor moves cloud
fraction by 0.01 to 0.03 absolute wherever the modelled relative humidity sits
near threshold. Cloud fraction is the largest single lever on planetary albedo,
so this is a first-order radiative quantity resting on two anchored-to-T21 fits.

**Disposition: IRREDUCIBLE as scalars; REPLACEABLE as a form.** Neither number
can be sourced, because what they encode is a subgrid distribution and the
subgrid scale is a property of the mesh rather than of the world. What CAN be
done is make the resolution dependence EXPLICIT rather than declared, which is
what `rcritmod` and `rcritslope` exist for and neither appears anywhere outside
`vendor/`. A tuning that is a declared function of the grid can be tested
against the grid; one anchored to a truncation with a comment cannot.
`world-khn` asked the resolution question and closed it by declaring the anchor
with no numeric changed, which is the right answer to that question and leaves
this one.

---

### 14. `a1 = 0.75`, the one fitted coefficient in the sigma quartic

`vendor/exoplasim/exoplasim/plasim/src/plasim.f90:2025`. The declaration already
separates the determined part from the fitted part, which is exactly the right
treatment and is why this ranks low: two of the three coefficients follow from
`sigmah(1) = 1` and a vanishing derivative at the surface, and "The ONE fitted
number is a1, the slope at the model top, and it sets how much of the column the
upper half spans. It is upstream's and carries no derivation there."

**Disposition: IRREDUCIBLE.** A vertical grid is a discretisation choice, not a
property of the atmosphere. What a sourced form would look like is a convergence
argument across `a1`, which is a measurement this project can make and no paper
can supply.

---

### 15. The pedology blocks that carry no bracket

The pedology sweep's own headline: the tuned values in that component are
concentrated in `pedology/config/pedogenesis.yaml`'s `regolith:` and `texture:`
sections, and those are the only blocks in the file that do not ship brackets.
Every other declared-without-a-source value in either that component or
hydrography (`catena`, `water`, `andisol`, `bedrock_water`, the loess threshold,
`f_grad`, `lambda_m`) ships a bracket and an instruction not to tune within it.

| where | value | what is missing | magnitude |
| --- | --- | --- | --- |
| `pedogenesis.yaml:236` | `clay_conversion = 0.55` | no comment on the value, no source, no bracket; numerically identical to `clay_yield` at `:298`, which WAS fitted to a SoilGrids slope of 0.543 | near-linear at this world's weathering intensity: a 10 per cent move is +9.2 per cent converted clay |
| `pedogenesis.yaml:298` | `clay_yield = 0.55` | fully admitted as an Earth fit over sixteen type localities, and applied BELOW the range it was calibrated in, where the file records the bias is largest | linear multiplier; documented residual bias of +0.084 clay fraction at this world's mean |
| `pedogenesis.yaml:565` | `dry_erosion_baseline = 0.15` | the mechanism is argued and the value is not attached to it; no source, no bracket | 21 per cent of the moisture term at the land mean, 100 per cent of it over the arid 18 per cent of land, where depth goes as its reciprocal |
| `pedogenesis.yaml:619-633` | the pH block: six parent values, `leaching_slope` at `:629`, `endorheic_alkalinity_bonus = 0.8` | no comment on ANY value, no citation, no bracket, nothing in README, notes or references | 0.9 moves land-mean soil pH by 0.04 per 10 per cent and wet-cell pH by 0.13; the bonus is a flat 0.08 over the closed-basin fraction that the brine and duricrust rules key on |
| `pedogenesis.yaml:301` | `sand_to_silt_loss_ratio = 2.0` | no source; and the mechanism sentence beside it has the surface-area argument backwards, which suggests the number did not come from the reasoning printed next to it | redistributes 2 to 3 per cent between sand and silt at fixed clay; under 1 per cent on PAWC |

**Disposition: REPLACEABLE NOW, as brackets.** The cheapest correct move here is
not new physics. It is the bracket discipline the rest of the same file already
enforces on itself, applied to the five rows above, so that what is unknown is
declared as unknown and gets swept rather than quoted. SoilGrids `phh2o` over the
sixteen sites `pedology/scripts/validate_against_earth.py` already joins to would
bracket `leaching_slope` directly, with the caveat the file itself raises that
those sites have absorbed one fit already.

---

### 16. The LPJ-GUESS constants whose own source says they were invented

| where | value | the source's own words |
| --- | --- | --- |
| `vendor/lpj-guess/framework/guess.h:2552`, comment at `:2546` | `frac_maxtomin = 0.9` | `// Tighter C:N ratio range for roots and sapwood: picked out thin air.` |
| `vendor/lpj-guess/modules/soil.cpp:111` | `sompool[PASSIVESOM].ntoc = 1.0/9.0` | `// passive has a fixed value (why? passive SOM should also vary.)` |
| `vendor/lpj-guess/modules/somdynam.cpp:50` | `TAU_LITTER = 2.85` | `// Thonicke, Sitch, pers comm, 26/11/01` |
| `vendor/lpj-guess/modules/vegdynam.cpp:1435` | `fireprob` floor `0.001` | `// c.f. LPJF` |

The first two are live under `ifnlim 1`. `frac_maxtomin` sets the width of the
fine-root and sapwood C:N windows, and this project has already repaired their
PLACEMENT from Friend et al. (1997) Table 4 under `world-xms4` while leaving the
width unsourced for both elements; `guess.h:313` names that as a whole-model gap.
`ntoc` sets how much nitrogen is locked in a pool with a turnover of order 1400
years, so it sets steady-state mineral nitrogen after `equilsom`.

`TAU_LITTER` is on the `ifcentury 0` path and this project sets `ifcentury 1`, so
it is dead. The fire floor is not: GLOBFIRM is the live fire model, and the floor
puts a minimum burn of 0.1 per cent per year on every patch of the planet forever
on the authority of "the other model does it".

**Disposition.** `ntoc` is REPLACEABLE NOW: Parton, Stewart and Cole (1988) is
held and read, and it is the paper this project already used to settle the
phosphorus analogue, `PMASS_SAT` and `SURFHUMUS`. The nitrogen side was simply
not given the same treatment. `frac_maxtomin` is REPLACEABLE, SOURCE NEEDED and
the source is a C:N range rather than a ratio, so it is a literature question
rather than a fetch. The fire floor is REPLACEABLE NOW by deletion: a floor whose
only justification is another model's floor is not a physics term, and removing
it is a fix.

---

### 17. `ntransform.yaml`'s seven `unsourced` entries, already registered

`biosphere/config/ntransform.yaml` lines 687, 702, 732, 815, 834, 966 and 1054.
Recorded here so the survey is complete and so nobody re-derives them: this
project already knows, `ntransform_gate.py --strict` already refuses on them,
and the register carries the evidence per entry. The two worst by the register's
own account are the water-filled-pore-space partition, whose midpoint has two
competing numbers in the code's own comment and whose shape parameter has none at
all, and the Michaelis-Menten divisor, whose two readings differ by 6.7 to 30
times.

`f_denitri_max` is the clearest single instance in the set:

> Table 9 eqn 3 is DNmax * ftemp * NO3/(NO3 + Kn) with no maximum-rate constant
> beside it, and table 11 lists none. The port multiplies that expression by
> 0.33, a ceiling of a third of the anaerobic NO3 pool per day, which neither
> paper states.

**Disposition: correctly held.** The register is the model this audit recommends
for every other component.

---

### 18. Tuning knobs that are currently inert, and the traps in them

A tuning that no run reaches costs nothing today and costs everything on the day
someone throws its switch. Each of these is recorded with the switch that arms
it, which is the form `opaque-constants.md` finding 10 established after the
ocean radius was found under a scheduled measurement.

| where | value | armed by |
| --- | --- | --- |
| `radmod.f90:76` | `tpofmt = 1.00`, self-declared "tuning of point of mean transmittance" | nothing; it is at its identity and no run sets it. It is a live multiplier at `:3591` with no source, so any run that moves it moves an unsourced knob |
| `carbonmod.f90:97-98` | `tune1 = 5.41`, `tune2 = 2.20`, commented "Tuning adjustment to make global average match" | `NCARBON = 1`. The module header already says it is fitted to Earth by construction |
| `hurricanemod.f90:50-114` | the pLCL empirical parameters and `CL`, which the header calls an admitted fudge | the diagnostic switch. `run_exoplasim.py:2629` already refuses it on these grounds |
| `simba.f90` | Earth-fitted throughout, with turnover times scaled by the ORBITAL period | `NVEG = 1` |
| `vendor/lpj-guess/modules/blaze.cpp:59-66` | four litter combustion factors spanning 300 times, selected by ABSOLUTE EARTH LATITUDE at `:1003-1017`, under the comment `// Latitude depending tuning values mortality` | `firemodel "BLAZE"`. `run_lpj_guess.py:148` overrides to GLOBFIRM, and `global.ins:97` still says BLAZE |
| `vendor/lpj-guess/modules/soil.h:211,299` | `CH4toCO2_inundated = 0.027`, commented "to match global emissions"; `PEATLAND_WETLAND_LATITUDE_LIMIT = 40.0` selecting physics on Earth latitude | `ifmethane 1`. `biosphere/config/wetlands.yaml:138-142` already names the 40.0 and forbids latitude as a regime selector |
| `vendor/lpj-guess/modules/soilinput.cpp:350-352` | `kplab`, `spmax`, `pwtr` at Amazon-FACE site values, marked `// FIXED FOR AMAZON FACE AT THE MOMENT`, applied planet-wide | `ifplim 1` with the pedology contract absent. The live path at `:407-411` reads the contract |
| `vendor/exoplasim/exoplasim/pyburn.py:2003, 2843` | `tstar[tstar<255] = 0.5*(255+tstar)`, ECMWF's cold-surface guard, calibrated on Earth's surface temperature distribution | anything reading `psl`. Nothing does yet, and at 32 degrees of obliquity the guard fires over a large share of the land for much of the orbit |

---

### 19. The cgenie tier: a calibration subsystem, not a set of constants

`vendor/cgenie/genie-rokgem/src/fortran/rokgem_lib.f90:207-285` declares
`opt_calibrate_T_0D`, `_R_0D`, `_P_0D`, `calibrate_weath`, four
`calibrate_weather_*` factors and `opt_calibrate_T_2D`, `_R_2D`, `_P_2D`, with
`par_data_T_0D`, `par_data_R_0D` and `par_data_P_0D` as the Earth global means to
calibrate to and `par_ref_T0_2D`, `par_ref_R0_2D`, `par_ref_P0_2D` as Earth
pattern files to reproduce.

This is the class in its most concentrated form anywhere in the tree, and it is
not a list of numbers. It is machinery whose PURPOSE is to rescale modelled
temperature, runoff and productivity fields until they match Earth observations,
with the reference patterns shipped as data files. On this world every one of
those options must be off and every calibration factor must be 1.0, and nothing
in this project says so.

Beside it, `references/INDEX.md`'s entry on Holden et al. (2016) records that the
same family's ocean coupling carries three more: `scf`, which the paper states is
"a TUNED ensemble parameter, not a regridding artifact and not derivable from a
stress-product ratio"; an Atlantic-Pacific moisture flux adjustment set against
Talley (2008) basin budgets; and energy flux corrections diagnosed against
observed present-day sea-ice thickness.

**Disposition: IRREDUCIBLE, and that is the argument for the boundary rather than
for the values.** cgenie is the candidate under OCN-3, nothing reads it, and it
does not build where it stands. The correct output of this row is not a fix but a
precondition: an adoption decision for any cgenie component states which
calibration switches are off and which reference files are absent, before the
component is wired to anything. `ocean-tier-implicit-earth.md` does not reach
rokgem.

---

## Checked and clean

Recorded so it is not re-audited, and because what a project does right is the
argument for the rule.

**The aeolian component refuses the textbook tuning explicitly.** Dust emission
schemes conventionally carry a global constant scaling total emission to an
observed dust load. There is none here, and `aeolian/config/dust.yaml:26-29`
gives the reason: Kok (2014) was chosen over Marticorena-Bergametti "because K14
derives the emitted size distribution from fragmentation physics rather than
fitting it, so it carries one fewer unconstrained knob. On a synthetic planet
with no aerosol optical depth observations to tune against, a globally tuned
constant is worth nothing." Lines 8 to 12 pre-commit against tuning toward the
reopening test. `aerocore.f90:1270` repeats the argument in the model source.
Every unsourced value in that file and in `sea_salt.yaml` carries a bracket and
the model reports the spread.

**The shortwave band weights are derived, not fitted.** `h2o_sw_weight`,
`co2_sw_weight`, `cloud_absorption_scale` and `h2o_sw_level` each have a named
generator, a recorded route and a note carrying the argument. The last is the
closest to a fit and is the most carefully argued key in `config/planet.yaml`:
what is derived is the arithmetic, both factors it multiplies are declared, and
the bracket ships as arms rather than as an error bar.

**`earth_calibration.py` is a validation and not a knob.** It fits a bias and
scale inside a split loop and discards them on every iteration; what it returns
is a held-out score against a bar declared in advance in
`hydrography/notes/earth-calibration-criterion.md`. Its output artifact is read
by no code, and the solver's own parameters are sourced independently in
`groundwater.yaml`, which explicitly refuses to prefer the arm that matched best.

**`ntransform_gate.py` discusses tuning and contains none.** Its 37 occurrences of
the vocabulary are all enforcement machinery and fixture names.

**`biosphere/scripts/build_vesper_header.py`'s solstice offset is fitted and is
not tuned.** It is a grid search against the climatology's own declination, it
reports rms and maximum residuals into the header and the provenance JSON, and it
re-derives per planet. That is the distinction this audit turns on, in one
function.

**`scripts/error_budget.py:108`'s `DEFAULT_ATTENUATION = 0.5`** is the one that
looked like a fudge and is not: "A factor of 0.5 is used here as a conservative
default and both columns are reported, because the honest claim is a factor of two
rather than a number."

**Also swept and clear**, in less detail because the pattern is the same: the
`minerals/` prospectivity rules, which carry a declared four-value grounding
vocabulary per rule; `analysis/`, all 24 scripts; `lib/`, where `sensitivity.py`'s
slope is from two runs still in `archive/runs/` and corroborated on a second
terrain; `config/pipeline.yaml`, which carries no bare numeric parameters at all
and whose tuning vocabulary is entirely prohibitions; `hydrography/config/`,
where every Gleeson permeability is read from the paper's Table 1 with its own
sigma; `pedology/config/surface_classes.yaml`, `outgassing.yaml`,
`land_column_properties.yaml` and `weathering_schemes.yaml`;
`weathergen.cpp`, whose hundreds of literals are Numerical Recipes rational
approximations and which is dead at `weathergenerator "INTERP"`.

## What this audit did not cover

Stated so the coverage is not overread. Orogen's own generation constants were
not swept. `maps/` was swept and its constants are rendering choices with no
physical claim riding on them, so they are excluded by construction rather than
cleared. `vendor/lpjml/` was not swept; nothing reads it. The per-PFT
bioclimatic limits and BVOC emission capacities in the LPJ-GUESS instruction
files are Earth species-distribution fits by construction, roughly 100 and 36
vectors respectively, and are recorded as a class rather than enumerated: no
sourced replacement exists for any of them, keeping them is a declared
worldbuilding decision, and the BVOC set is already bracketed as floors and
dormant. `tswr1`, `tswr2` and `tswr3` are excluded because they were being
settled separately while this ran.

## Tasks

Tracked in the `bd` issue tracker under the `tuned-value` label, not restated
here. Two of the rows above were already tracked before this audit ran and no
new bead was filed for either: `clim-30` for the design flux's declared cap, and
`bio-26` for a bounded sensitivity on `ALPHAA_NLIM`. Row 7 was fixed rather than
filed.

| row | id |
| --- | --- |
| 1. `baseline_flux_earth` and the cap solved backwards | `clim-30` |
| 2. `ALPHAA_NLIM` | `bio-26` |
| 3, 4. the regolith depth level and its asymptote | `world-qs63` |
| 5. the carve size floor | `world-7vj6` |
| 6. `eddy_wind_m_s` | `world-9y07` |
| 7. `acllwr` | fixed |
| 8. `th2oc` | `world-2esd` |
| 9. the `dz0land` anchor | `world-u8ds` |
| 10. `vdiff_lamm` and the Louis fives | `world-x6q8` |
| 11. `frac_labile_carbon` | `world-vyvn` |
| 12. `gamma` | `world-trs3` |
| 13. `zcca`, `zccb`, `rcrit` | `world-o12h`; `world-khn` closed the resolution half by declaring it |
| 14. `a1` in the sigma quartic | none: irreducible, and a convergence sweep across it is a measurement rather than a fix |
| 15. the unbracketed pedogenesis values | `world-9ctm` |
| 16. `ntoc`, `frac_maxtomin`, the fire floor | `world-xfif`, `world-v5j1` |
| 17. `ntransform.yaml`'s seven `unsourced` entries | already registered; the gate refuses on them |
| 18. the dormant knobs | `world-9g8p` |
| 19. the cgenie tier | `world-u9kg` |

Three papers are named as needed and none is in `references/`. They are the
whole of what this audit's dispositions wait on that cannot be done from inside
the tree: Li et al. (1992) `10.1029/92JD00509` for `frac_labile_carbon`, Louis
(1979) `10.1007/BF00117978` for the three fives, and Kessler (1969)
`10.1007/978-1-935704-36-2` for the re-evaporation form. Portenga and Bierman
(2011) `10.1130/G111A.1` is already on disk and marked *held*, so reading it is
the fourth and cheapest.

