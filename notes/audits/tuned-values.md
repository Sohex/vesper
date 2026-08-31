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

The three are project vocabulary and `docs/src/reference/vocabulary.md` is
canonical for them; this table is the restatement an auditor works from. The
one that is read backwards is the third: IMPLICIT-EARTH IS A DIAGNOSIS AND NEVER
AN ENDORSEMENT. Proving a constant is Earth's own number is what establishes it
is wrong on this planet, and a reclassification out of `tuned` is progress in
knowing rather than a repair. Row 9 is the worked example.

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
baseline climatology, which this build does not have. No paper is needed. The
script now REFUSES to run with the cap undeclared, before it reads a field, so
the block is enforced rather than described.

#### The criteria that replacement rests on, audited one at a time

Re-examined 2026-08-25 under world-bnu6, in the window before the next
commissioning cycle produces the runs the step reads: a threshold re-examined
after those runs exist is not a threshold. The trigger was that the shortwave
cloud optics moved by up to 20.6 K under `world-f9ig` and `world-jgen`, and the
comfort band is only 58 K wide.

**A radiative change does not reach a threshold stated in degrees.** The band is
a set of temperatures; changing the radiation changes which stellar flux
delivers a given temperature field, so it moves the search's RANGE and leaves
the band where it is. `exoplasim/notes/trace-gas-absorbers.md` section 4 already
recorded the same mechanism from the other side, of a radiation scheme short of
greenhouse forcing: it "reaches those thresholds at a higher stellar flux, and
the search finds it". That is why only one entry in the table below moved.

| criterion | class | disposition |
| --- | --- | --- |
| `cold_floor_c` = -25 C | declared design preference, fixed ahead of its runs | KEPT, and declared with a swept bracket |
| `warm_extreme_c` = 38 C | declared design preference | KEPT |
| `cold_extreme_c` = -40 C | declared design preference | KEPT |
| `warm_ceiling_c` = 33 C | REVERSE-ENGINEERED from where the superseded physics put the tropics | RE-DECLARED at the same value, with a swept bracket |
| `cold_extreme_cap` | REVERSE-ENGINEERED; solved backwards from the answer | refused until declared: `clim-30`, above |
| candidate range | derived under the superseded optics, and could no longer reach its own answer | re-derived, and an edge winner now refuses |

The warm ceiling is the entry that moved, and not because 33 C stopped being a
reasonable preference. Its recorded justification was that the anchor flux "puts
the tropics near +33 C in their warmest month rather than +36", which is a
number read off the old answer -- the same shape as the cap. The corrected
optics delete that justification outright rather than shifting it, because at
the anchor flux the modelled tropics are now several kelvin colder and the
sentence no longer picks out 33. So the value survives only if it can be stated
as a preference in its own right, which is what the script now does. A
re-declaration at the old number is where a tuning would hide, so it takes the
third of the four dispositions rather than the bare first: declared WITH A
BRACKET THAT GETS SWEPT, `THRESHOLD_BRACKET` in the script, with the design flux
each bracket point returns carried in the artifact. That costs no run -- the
projection is already computed at every candidate -- and it replaces an
assertion that the answer does not depend on the preference with a measurement
of how much it does.

The other three thresholds are the opposite case and are worth stating rather
than leaving inherited. They were fixed in the script's docstring ahead of the
data, nothing in the model or in any run determines them, and they say what kind
of world this is meant to be: land is comfortable if its warmest binned month
averages at most the ceiling and its coldest at least the floor, harsh out to
the extreme thresholds, uninhabitable beyond. `conventions.md` allows a
threshold fixed in advance to be a judgement precisely because it was fixed in
advance, and a design preference does not move when the physics does.

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

`pedology/config/pedogenesis.yaml`'s `regolith:` block, consumed in
`pedology/scripts/build_soil.py:regolith_depth`. The two keys are now ONE,
`erosion_coefficient_per_relief_m`, and the paragraphs below record what they
were and why they collapsed:

    erosion = erosion_coefficient_per_relief_m * erodibility * max(relief_m,0) * moisture
    depth   = maximum_depth_m * production / (production + erosion)

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

**Disposition: DEGENERACY REMOVED, LEVEL BRACKETED FROM THE PAIR. Both papers
are now read.** Portenga and Bierman (2011) was read on 2026-08-25 and it does
not contain what this row expected of it, which is recorded here rather than
left standing: it publishes **no relief-to-erosion relation**. Mean basin SLOPE
is the dominant global regressor for drainage-basin denudation and the top one
in nearly every climatic, lithologic and tectonic subpopulation; relief is a
secondary bivariate correlate and is "unimportant for most categories of
outcrops"; and the article gives no regression equation for either. Its
by-lithology result is for OUTCROPS, which are bare rock and carry no regolith,
and for drainage basins lithology does not separate the rates at all.

What the pair does supply is two separate things, and separating them is the
point.

**The FORM is sourced.** `relief_m` is the height spread across a cell's
neighbours, so at fixed mesh spacing it is a slope, and erosion rising with it
is Portenga and Bierman's leading result rather than an assumption. That also
retires the unsourced normalisation comment: what the term is proportional to is
now stated by a source, and how it is scaled is one number.

**The LEVEL is a bracket and not a value, and the pair does not close.** A
soil-mantled hillslope at steady state has production equal to erosion, so
Heimsath's measured function inverts to `h = 0.435 m * ln(77 / E)` with `E` in
m/Myr, and Portenga and Bierman supply `E` on a standardised global footing.
Their drainage-basin median of 54 m/Myr gives 0.15 m. Their outcrop mean of 12
gives 0.81 m and their outcrop median of 5.4 gives 1.16 m. Their slowest
subpopulation, polar outcrops at 3.9, gives 1.30 m. **Their drainage-basin mean
of 218 and their arid-basin mean of 100 both exceed Heimsath's maximum
production rate of 77 m/Myr, where the balance has no solution and no soil is
possible at all.** So the absolute production-to-erosion ratio this row hoped
for is constructible and is not single-valued: it brackets the land-mean
regolith at 0.15 to 1.30 m and contains a regime with no steady state. The width
is a property of the sources, not of the arithmetic -- Heimsath's 77 m/Myr is
one greywacke site and Portenga and Bierman's basins cover 2.3 per cent of
Earth's land with a stated accessibility bias -- and neither paper narrows it.

**The two degenerate keys are one.** `erosion_weight` and
`erosion_reference_relief_m` are replaced by
`erosion_coefficient_per_relief_m = 1.5/200 = 0.0075`, with the law and every
number it produces unchanged. One degree of freedom is now carried in one
watched key, inside a declared bracket, with `regolith_depth_bracket_m` beside
it as what a sweep runs over. The declared value puts the land mean at 1.04 m,
inside the bracket and in its upper third, which is the branch reached by
outcrop denudation rather than basin denudation.

**What is left open** is converting the depth bracket into a coefficient
bracket, which is one re-run of `pedology/scripts/build_soil.py` per endpoint on
a build with a climatology and is not a model run.

---

### 4. `maximum_depth_m`, the asymptote nobody labelled

`pedology/config/pedogenesis.yaml`'s `regolith:` block. The whole comment used
to be "Asymptote, not a clip: the depth to which weathering can reach at all."
The block above it argues at length for the saturating FORM over Heimsath's
exponential and concedes that it is a parameterisation rather than a
derivation; neither said where the number came from. No source, no bracket.

**Magnitude.** A pure prefactor: a 10 per cent move changes every cell's
regolith depth and therefore land-mean PAWC by exactly 10 per cent. Per unit
fractional change it is a larger lever than the erosion coefficient, and unlike
that coefficient it was not labelled free anywhere.

**This row's first disposition was WRONG and is recorded rather than deleted.**
It said "same pair of papers as row 3; a weathering-front depth is what a
production function integrates to". Having read both, that is false. A
production function integrates to a STEADY-STATE depth against an erosion rate,
which is row 3's bracket; it does not integrate to a maximum reach. Heimsath's
form has no asymptote at all, running to infinity as erosion vanishes -- that
divergence is the reason this file replaced it. Portenga and Bierman measure
denudation and not the depth a weathering front reaches. Neither bears on this
number and no amount of re-reading either would have changed that.

**Disposition: SOURCED AS A BRACKET, from Shangguan, Hengl, Mendes de Jesus,
Yuan and Dai (2017), *Mapping the global depth to bedrock for land surface
modeling*, JAMES 9, 65-88, `10.1002/2016MS000686`.** Read 2026-08-26. Their
Table 1, world row, absolute depth to bedrock over 1,590,464 soil-profile and
borehole observations: median 670 cm, mean 1,309.3 cm, maximum 312,541 cm.

The measurement is not the modelled quantity, and it misses it in two
directions the paper itself names. Both become an end of the bracket, which is
the same shape row 3's disposition took and for the same reason.

- **The median bounds it from below.** The observed population is
  overwhelmingly eroding land, and an eroding profile sits under the asymptote
  it would reach without erosion. A second reason runs the same way: a soil
  survey logs intact weathered rock as R horizon, so a profile-derived map
  reads as deeper than what it prints.
- **The mean bounds it from above.** Depth to bedrock counts transported fill
  as readily as in-situ weathering product, and this model's regolith is
  in-situ product only. The mean is 1.95 times the median and the maximum is
  3.1 km, which is the size of that contamination.

So the bracket is 6.70 to 13.09 m, and the declared value is the low end
because every correction that separates the observation from the modelled
quantity runs upward from the median and downward from the mean, and the
depositional contamination is the larger of the two. The observation set is
dominated by North American and European boreholes, where glacial drift and
alluvium are most of what the drill passes through.

**The PDF is not held and the reason is worth recording, because it is a state
between "fetched" and "closed access".** The article is gold open access --
Unpaywall, OpenAlex and Semantic Scholar all name the publisher's `pdfdirect`
URL as the OA copy -- and that URL returns HTTP 403 to `paperfetch`, to plain
curl, to a browser-impersonating curl carrying a session cookie, and to a
reader proxy. There is no repository copy and no Wayback snapshot. The full
text was read through the publisher's own HTML article page instead, and the
table was transcribed twice independently before being used.

**What this row also exposed, which no key's own bracket says.** This
prefactor and `erosion_coefficient_per_relief_m` are JOINTLY constrained: the
prefactor scales and the coefficient carries the level, and the land-mean depth
the pair produces has to land inside `regolith_depth_bracket_m`, which brackets
the LEVEL and not either key. Moving one alone can put the level outside a
bracket neither key mentions. A sweep runs the two together, and nothing in the
tree checks it, because no gate reads `pedogenesis.yaml`'s brackets at all.

---

### 5. The carve incision coefficient, and the size floor that sets its target

`hydrography/scripts/export_carve_list.py`, `calibrate_coefficient` and
`sweep_size_floor`, against `EARTH_SIZE_FLOOR_KM2`, `EARTH_STANDING_BY_FLOOR`,
`EARTH_BAND_LAND_MKM2` and `CALIBRATION_LATITUDE`.

`C` is the leading coefficient in `cut = C * erodibility * slope^n * Q^m`, and
therefore decides `retain` for every basin in the catalogue, and therefore what
Orogen carves.

**`C` is IRREDUCIBLE and it is declared as such.** It absorbs a relaxation
window that `docs/src/reference/no-time-axis.md` establishes cannot exist in a
generator with no time axis, so no dataset can supply it, and re-solving it
against this world's own basin population every run is the right treatment of a
quantity of that shape: a literal would go stale the moment the verdict moved,
silently, and it did once already. What is durable is the Earth measurement and
the expected-value argument, not the number.

**The size floor on the Earth sample WAS the undeclared lever, and it is now
declared and swept.** It sets which lakes the Earth density is a density of, so
it sets the calibration target, so it sets `C`. It used to be a single value
defended in a comment as what the mesh resolves, with the file admitting in the
same breath that counting from 10 km2 instead "gives a density 33 times higher
and reads as 'carve almost nothing'".

`sweep_size_floor` re-solves at every rung of a measured ladder, over a span the
build derives rather than one chosen -- the median land mesh cell area at the
bottom, below which a counted Earth lake has no representable counterpart, and
the median area at spill of the overflowing basins at the top, above which most
of the population being solved for is outside the class the Earth sample stands
for. The result goes into the carve list's sidecar under
`method.calibration.size_floor_sensitivity`, beside the verdict, with the whole
ladder including the rungs the span excludes.

**Magnitude, and it is well clear of the noise.** Measured on
`precarve-craton-10m`; `hydrography/notes/retain-fraction.md` carries the table.
Over the derived span the solved `C` moves by a factor of 8.2 and the marginal
class by a factor of 9.2, against 1.8 for the Poisson error on the Earth count
at the floor in force. The floor is the larger lever by a factor of about five,
which is why leaving it in a comment was the defect and reporting it is the fix.

**Disposition: IRREDUCIBLE for `C`, with the floor DECLARED AND SWEPT.** Neither
is now a tuned value in the sense this audit uses: `C` has the argument for why
no source can exist, and the floor has a bracket that is reported wherever the
number it produces is.

---

### 6. `eddy_wind_m_s`, re-measured on runs that exist. RESOLVED

`config/planet.yaml`, `model.hyperdiffusion`. The single free quantity in the
hyperdiffusion rule `tau_vorticity = pi * radius / (NTRU * eddy_wind)`, dividing
into every entry of the `timescales_days` table for all five rungs and four
fields, which reach the model as TDISSD, TDISSZ, TDISST and TDISSQ.

**What was wrong.** The value was 5.94 m/s and the comment beside it was a
measurement claim on `run_4182235e9781` orbit 39. That run is not on this tree:
it survives as a `restart_from_run` string in seven manifests and one INDEX row
in `archive/runs/_bulk_2026-08-24/INDEX_AT_DELETION.json`, and the archive's
README records that the directories were removed by an `rm -rf` whose argument
came from an unguarded command substitution that expanded to nothing. No
manifest, no series, no NetCDF. So the row failed this audit's criterion on
"cannot be inspected, cannot be checked" rather than on "fitted".

**What it is now.** Re-measured 2026-08-26 at 10.09 m/s on
`run_57cecaa8391b`, 25 orbits and 800 instantaneous snapshots, by
`exoplasim/scripts/measure_eddy_wind.py`, which states the definition in full and
writes its verdict under `exoplasim/analysis/eddy-wind/` so the derivation
outlives the run directory. Five T21 runs -- four links of one 110-orbit restart
chain and one cold start independent of it -- give 10.05 to 10.18, a spread of
1.3 per cent; the within-run standard error is 0.2 per cent against a 5 per cent
criterion fixed before the measurement. Declared with a +/- 5 per cent bracket.

**The old value was a different quantity, and this is checkable.** Reducing the
TIME-MEAN output stream by the same recipe returns 5.82 to 5.97 on the same five
runs, and 5.94 sits inside that spread; reducing the INSTANTANEOUS stream returns
10.05 to 10.18. A time mean carries the stationary eddies and drops the
transients, at a ratio of 0.58 on every run measured. That is not a choice
between two defensible readings: nothing is advected by a time average, so an
advective time built from a time-mean field is not an advective time. The
correction also removes the gap the block's own corroboration could not explain,
since 5.94 sat a factor of 2.1 to 2.9 below the 12.6 to 17.5 m/s ECHAM5's tuned
table implies for Earth and 10.09 sits 1.25 to 1.7 below it.

**Magnitude of the change.** Linear. Every damping timescale in the table
shortens by 1.699x at every rung. The operator is applied implicitly in the
damped variable (`plasim.f90:5148`), so no stability ceiling and no binary
changes, and the keys are namelist keys, so no rebuild.

**The one term left, and what sweeps it.** These runs were damped with the
previous value, so 10.09 is a first iterate rather than a measured fixed point.
The next run reads the new value and the same instrument re-measures it at no
marginal cost; the criterion is the declared bracket, within 5 per cent the value
stands and outside it the value is replaced and the run repeated.
`exoplasim/notes/resolution-tuned-parameters.md` carries the bound on how far
that iterate can move and the argument for the recipe.

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
column path `w` in g/cm2. THE PATH RANGE THIS ROW USED WAS TOO NARROW: it said
0.5 to 2.5 g/cm2, and the model's own paths run 0.04 to 6.08 with a mean of
2.23, so the term is larger than this row stated at the wet end. Setting `th2oc`
to zero, which the surrounding `if(th2oc > 0.)` supports, removes 0.052 of
broadband absorptance at the mean path and 0.136 at the maximum, which through
`lib/sensitivity.py` is roughly 5 to 12 K end to end -- above the 2 to 7 K this
row first carried.

**AND THE INHERITED VALUE IS HALF ITS OWN CEILING.** The term stands in for the
window, so it cannot add more absorptance than the window's Planck share, which
over 250 to 300 K at the largest path caps `th2oc` at 0.038 to 0.050. At 0.024
it already claims more than half the whole window. That ceiling is derived from
what the term IS rather than fitted, which is what makes the bracket below a
construction rather than a guess.

**Disposition: DECLARED BRACKET, SWEPT.** `[0.0, 0.038]`, with the arms run end
to end and the static conversion's high lean stated beside them. Sourcing it is
not available and the evidence says so: `references/INDEX.md` records that
Mlawer et al. (2012), `10.1098/rsta.2011.0295`, read, "does not supply an
evaluable continuum" -- the coefficients ship as data with LBLRTM. The route to
a VALUE is a correlated-k or line-by-line calculation with MT_CKD on this path,
which is the same bundle `exoplasim/notes/corrk-cross-check.md` says the project
does not have and which also blocks `h2o_sw_level`'s bracket; the two open
together. Until then a swept bracket is the honest disposition, and it is the
third of the four this project allows.

---

### 9. The land roughness field's orographic term. DERIVED, no free coefficient

`exoplasim/scripts/build_surface_roughness.py`. Two constants were in question
here and neither survives as a tuning: `EXOPLASIM_DZ0LAND_M`, the model's
namelist fallback over land, and the coefficient that used to be bisected so
this world's land mean landed on it.

**The orographic term is now derived from this world's own subgrid slope.** The
scheme is the one large-scale models use for turbulent orographic form drag, and
every constant in it is read out of the papers that established it:

- Wood and Mason (1993) Eq. (33). The effective roughness of hilly ground is the
  length whose drag coefficient at the pressure scale height equals the SUM of
  the hills' pressure drag coefficient and the ground's skin friction. The terms
  add in the drag coefficient. `simba.f90:474` combines its vegetation and
  topography terms in quadrature in the LENGTH, which is the model's convention
  and not the one the literature derives, and this field follows the derivation.
- Beljaars, Brown and Wood (2004) Eq. (6), which is WM93's surface pressure drag
  with the hills' base area set equal to the domain area for mountainous
  terrain: `Ca = 2 alpha beta Cmd theta_bar^2`. The form drag is QUADRATIC IN
  SLOPE, not linear in relief amplitude. Their Eqs. (4) and (5) give the pressure
  scale height, the inner-layer depth and `alpha`, so none of the three is
  assumed; `beta = 1` is WM93's Table 1 for the parameter this takes, since
  packed three-dimensional hills give effective roughnesses of the same order as
  two-dimensional ones at the same `A/Sd`.
- Mason (1988). The area-average roughness is the one that reproduces the correct
  area-mean surface stress, and it is obtained by averaging DRAG COEFFICIENTS at
  the BLENDING HEIGHT, his Eq. (14) from the scale over which the cover varies.
- Lettau (1969) as an independent geometric cross-check.

**The old coefficient was not merely unsourced; the relation it sat in had the
wrong dependence.** `z0_oro = coefficient * stdev(elevation)` takes an amplitude
where the drag takes a slope, so the coefficient had to carry the missing
horizontal scale. That is why it could not be sourced, and why it drifted with
the support. Once the slope is used directly, and the mesh supplies it from its
own plane fit, there is nothing left to solve.

**The scheme is held to numbers its own papers print, and refuses when it misses
them.** Beljaars simplifies Eq. (4) to `h_m = 0.1 lambda` with `Cmd = 0.005` and
fixes `alpha = 12`; solving Eq. (4) at the wavelength those constants belong to
and evaluating Eq. (5) there returns 12.0 to three figures and `h_m/lambda`
within a twentieth of the paper's round number. Lettau's own worked orographic
example, Colorado's peaks at `h* = 1000 m`, `s = 5e6 m^2` and `S = 2e8 m^2`,
comes back at the 12.5 m he prints. WM93 Eq. (33) with no hills returns the cover
roughness exactly. A solver that has not converged, or an equation transcribed
with the wrong power, misses all of these by orders of magnitude.

**The declared bracket is on the scheme's two attributions and it is narrow.**
The wavelength the mesh-resolved slope belongs to, bracketed at one to four mesh
spacings around a central two, which enters only through logarithms; and the
turbulence closure, which WM93 measure at two levels for the same pressure force,
`5.9 theta^2` for the mixing-length closure against `3.48 theta^2` for the
second-order one. Measured on `canonical-10m-base` on 2026-08-26, the two
together move the land mean by about 4 per cent at T21 and 4 per cent at T42.
Both ends are computed on every build whichever is written.

**`EXOPLASIM_DZ0LAND_M` is implicit-Earth, and that is a reason to replace it
rather than to keep it.** The vendored model ships PlaSim's own boundary dataset
for the model's Earth configuration at `N032` and `N064`, carrying roughness in
three fields: code 172 the land mask, code 173 the total `dz0clim`, code 1730 its
topography-only part `dz0climo`. Reducing code 173 over that dataset's land in
`ce` returns the effective land roughness the model carries when it is given a
real map, and at `N032` -- the resolution PlaSim's own configuration defaults to
-- that reproduces the namelist fallback to within a tenth. So the fallback has a
derivation and it lives in a boundary dataset rather than in the source comment,
which reads "roughness length land". `earth_reference_land_z0()` computes it at
build time so a change under `vendor/exoplasim` moves it, and nothing about it is
a literal. What it is NOT is a target: a sound derivation for the wrong planet is
a diagnosis. Solving this world's land mean onto it would report the distance
between the two planets as zero by construction.

**The distance, measured rather than asserted.** On `canonical-10m-base` on
2026-08-26 the derived land mean is 0.40 m at T21 and 0.41 m at T42, against an
Earth reference of 2.20 m and 1.21 m at the same two rungs: the derived field is
0.18 and 0.34 of Earth's. The land-mean exchange coefficient is 2.1 times smaller
at T21 and 1.6 times smaller at T42 than the Earth-anchored field was, and `ce`
is what land sensible heat and land evaporation are linear in, so the carve
criterion's numerator moves by that factor.

**Where the whole of that gap is, and why it is a limit rather than a bracket.**
The two planets' LAND COVER roughness is not what separates them: the Earth
dataset's own cover term, recovered from the quadrature, and this world's cover
term agree to within a tenth at both shipped rungs, and neither is derived from
the other. The gap is entirely in the relief term. `Ca` is quadratic in slope and
slope variance is `integral k^2 F(k) dk`, so it is dominated by the shortest
wavelength present and the large scales contribute almost nothing. Orogen has a
measured terrain-information floor near 20 km and produces no relief below it;
Beljaars' scheme is written for the band BELOW 5 km and reaches it by
extrapolating a power-law orographic spectrum fitted to United States topography.
Supplying that band here would be importing Earth's small-scale terrain under a
derivation's name, so it is not supplied, and the derived roughness is what this
world's terrain supports.

**The reduction height moved, and Mason is why.** The cell average used to be
taken at the lowest model level. Mason's Eq. (14) puts the blending height an
order of magnitude below that on this mesh, and his whole result is that the
average gives extra weight to high roughnesses occupying small fractions of area
-- which is exactly the closed-basin case this field exists for, run the other
way. Measured on `canonical-10m-base` on 2026-08-26, moving the average to the
blending height raises the cover-only land mean by about a tenth. His procedure
is insensitive to the exact height and needs only its order, so one height is
solved for the grid and the per-cell spread is reported.

**Lettau does not agree at this scale, and that is a finding rather than a
failure.** His `z0 = 0.5 h* s/S` shares no constant with WM93's route and reaches
a roughness by a different argument, so it is a real check. Applied with one mesh
region as one roughness element it runs about two orders above the WM93 orographic
enhancement on this world's land. It grows linearly with the size of the elements
at fixed slope, where the WM93 form saturates through the logarithmic profile it
is built on, and the two are within a factor of a few at the kilometre scale the
operational schemes were built for. The divergence is at the ten-kilometre scale
this world's terrain floor puts the resolved band at, which is Beljaars' own
reason for assigning scales above 5 km to gravity-wave and blocking schemes
rather than to a roughness at all. The builder computes both and reports the
ratio; it is not gated, and no threshold on it was fixed in advance.

**The ladder half closed, and it closed by the term becoming a property of the
land.** The old solve was per grid, so two rungs differed by their terrain and by
their calibration at once and SPAT-8 had to pass one coefficient to every rung.
The derived orographic term takes the mesh's slope and the mesh's spacing, and
the grid rung enters only through which regions fall in which cell. Measured on
`canonical-10m-base` on 2026-08-26, the field's own land mean runs 0.403, 0.410,
0.419 and 0.432 m at T21, T42, T85 and T170, against a factor of 1.8 between the
two rungs the Earth-anchored arm could be built at at all. Most of that 7 per
cent is in the COVER term, whose own land mean runs 0.375 to 0.399 m over the
same span: the reduction is an expectation of `ce` and a finer cell holds fewer
surface classes, so the Jensen suppression falls as the ladder is climbed.
`analysis/spatial_reduction_gap.py`'s recalibration arm isolates the orographic
part by holding the cover reduction linear, and there the land mean spreads by
0.9 per cent from T21 to T170 where the solved coefficient it replaced moved by
a factor of 1.98. A derived land mean that moved with the rung would say the
scheme was still carrying the support inside it.

**And it builds on rungs the anchored arm could not.** The model ships a
boundary dataset at `N032` and `N064` only, so the Earth-anchored default
refused above T42 and a ladder comparison had to pass it a fixed coefficient.
The derived arm needs nothing from the vendored dataset and builds at every rung
the export carries.

**A consequence for SPAT-7 that has to travel with this.** The reduction-order
correction the roughness field carries -- averaging `ce` over a cell's surfaces
rather than averaging the lengths -- was marginal while the land mean was
anchored high, because a large orographic term swamped the two-decade contrast
between barren ground and canopy. It is not marginal against a derived land mean.
Measured on `canonical-10m-base` on 2026-08-26 at T42, the land-mean `ce` differs
between the two reduction orders by about a tenth, and the land area on which the
per-cell difference exceeds the step's own reference-height bracket is about three
quarters, against a few per cent before. The operator was already the right one;
what changed is that it now decides a first-order quantity rather than a tail.

**Disposition: DERIVED, with a swept bracket on two attributions and no free
coefficient. Done.** `--orographic-arm derived` is the default and what a build
writes. `--orographic-arm reference` rebuilds the Earth-anchored field as a
DIAGNOSTIC for pricing the change and no climatology lineage may rest on it;
`--orographic-arm none` writes the cover term alone. Every build reports all
three land means, both ends of the declared bracket, the Earth comparison and the
Lettau cross-check. Three refusals stand where there were none: a bisection that
ends on a bound of its search bracket has not solved anything, the `ce` inversion
is checked as the identity it is, and every implicitly defined length is checked
against its own defining equation rather than trusted to a fixed iteration count.

Two things the derived field carries that a reader should see before running on
it. Its maximum is far above the old field's: WM93's effective roughness reaches
hundreds of times the cover roughness over steep ground, their own Table 1 runs
to 300 for periodic two-dimensional ridges, and ECMWF caps at 100 m for the same
reason. This model's lowest level is lower than that cap, so a handful of cells
carry a roughness that is a large fraction of the height the model applies it at
and have no surface layer left to put a logarithmic profile in. `fluxmod.f90`
forms `ln(z/z0 + 1)`, so nothing diverges. Measured on `canonical-10m-base` on
2026-08-26 those cells are under a millionth of land area and carry under two
parts in a hundred thousand of the land-mean exchange coefficient at T42, and
none at all at T21. The builder reports both fractions and the share of land-mean
`ce` they carry, and does not gate on them.

Elvidge et al. (2019), *Uncertainty in the Representation of Orography in Weather
and Climate Models and Implications for Parameterized Drag*, JAMES 11,
`10.1029/2019MS001661`, bounds what any of this can buy: zonal-mean orographic
surface stress differs by a factor of four between operational models of
comparable resolution. The declared bracket here is well inside that, which says
the remaining uncertainty is not in the scheme's constants but in the terrain the
world does not have.

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

**Disposition: IRREDUCIBLE for the three fives, ANSWERED against the primary.
DERIVED for the mixing length.** Louis (1979) was read on 2026-08-25 and settles
the fives negatively; Blackadar (1962) was fetched and read on 2026-08-26 and
supplies the mixing length's transfer rule.

**The fives are ECHAM's re-fit and not Louis's**, which was the question this
row asked. Louis's unstable fit is `F = 1 - b Ri/(1 + c sqrt(|Ri|))` and his
stable one `F = 1/(1 + b' Ri)^2`, with `b = 2b' = 9.4` chosen so that `dF/dRi`
is continuous across neutrality -- a constraint rather than a fit, and the only
derived relation in the set. His `c` is not a constant: it is
`C* a^2 b sqrt(z/z0)` from a free-convection dimensional argument, with `C*`
7.4 for momentum and 5.3 for heat and moisture and `a^2` the neutral drag
coefficient. The model's `b = c = d = 5` with 2b and 3b folded into the
expressions corresponds to Louis's `b` near 10 and 15 against his own 9.4.
Reading further does not help: Louis says of his own values that they "are
rather uncertain because of the large scatter in the observations", so both sets
are fits and neither transfers on its authority.

**`vdiff_c` IS NOT A DROPPED DEPENDENCE, and the claim that it was is withdrawn
against the expression.** `world-awm5` asked whether the surface limb should
carry `sqrt(z/z0)` explicitly. It already does. `fluxmod.f90`'s `mktcoe` forms
`zdenom = 1 + 3*vdiff_c*vdiff_b*sqrt(-zri*(zbz0+1))*zkblnz2` with
`zbz0 = znl/dz0`, so `zbz0+1` IS `z/z0` from this world's own per-cell roughness
field, and `zkblnz2 = (vonkarman/ln(z/z0+1))^2` IS `a^2`. Factoring out the
Richardson number leaves `3*vdiff_c*vdiff_b*a^2*sqrt(z/z0)`, which is Louis's
eq. 20 term for term. The single constant fixes `C*` alone, at `3*vdiff_c` = 15
against his 7.4 and 5.3.

**The magnitude argument was wrong in the direction that decides it.**
`sqrt(z/z0)` never appears without the `a^2` that accompanies it in eq. 20, and
the two run opposite, because `a^2` carries `ln(z/z0)` in its denominator.
Measured across this world's land roughness span, read from
`exoplasim/inputs/t21/roughness_t21_report.json` rather than restated here, the
product `a^2*sqrt(z/z0+1)` moves by a factor of 5.2 at a lowest-level height of
331 m, by 4.1 at 150 m and by 6.1 at 600 m, and it is not monotonic in `z0` --
against the factor of 34 that `sqrt(z/z0)` alone suggests. The quantity that
varies by 5.2 is one the model already integrates at every gridpoint.

  RE-MEASURED 2026-08-26 on the DERIVED roughness field, span 0.003 to 3.59 m at
  T21. The figures on the superseded Earth-anchored field were 2.7, 1.9 and 3.3
  against 21, and reproducing them on that span is what checked this
  calculation. The ratios roughly doubled because the derived field's span is
  wider in logarithm, and the CONCLUSION did not move: a factor of five is still
  a quantity the model integrates rather than one it flattens. A magnitude
  argument that rests on a regenerated field has to be re-taken when the field
  is, which is why the span is now cited and not copied.

**The momentum and heat limbs are separated too, in the numerator.** `zrifm`
carries `2*vdiff_b` and `zrifh` `3*vdiff_b` over a shared denominator, where
Louis carries one `b` over denominators differing through `C*`. In the strongly
unstable limit the heat-to-momentum enhancement ratio is 3/2 here against
Louis's 7.4/5.3 = 1.40, the same distinction by the same amount to within the
scatter Louis reports on his own coefficients. Splitting `C*` would import one
fit's numbers into another fit's algebra to move that ratio from 1.50 to 1.40
with nothing to say which is right, so no numeric and no form change follows.
Recorded as irreducible with the argument in the declaration.

**The mixing length: the scale-height derivation stays unsupported, and
Blackadar supplies a different one.** 160 m is not Louis's either; his eq. 22 is
Blackadar's `l = kz/(1 + kz/lambda)`, he calls `lambda` "an adjustable
parameter", and he took it as 100 m. Nothing in Blackadar's or Louis's form makes
`lambda` proportional to an atmospheric scale height, so scaling 160 m by 0.766
would be a guess wearing a derivation's clothes and it is not done.

Blackadar (1962) eq. 25 gives `lambda = 0.00027 G/f`, and his argument for the
form is dimensional: `z0` is "ruled out as a factor affecting characteristics of
the free atmosphere", leaving `G/f` as the only length the neutral problem
supplies. **The value does not transfer and the scaling does.** The 0.00027 is
fixed by matching one observed surface wind deflection, 32 degrees at Brookhaven
at `z0` = 1 m and `G` = 10 m/s, and he names `lambda` proportional to `u*/f` as
"a priori just as acceptable". His relation returns about 26 m at Earth
mid-latitudes against the 160 m here, because his `lambda` is the neutral
boundary layer's asymptote where ECHAM's is the value that carries
free-tropospheric mixing through the whole column: importing 26 m would be a
category error.

What `fluxmod.f90` now integrates is ECHAM's anchor on Blackadar's rotation
scaling. `f = 2*Omega*sin(lat)`, so at fixed latitude and geostrophic wind
`lambda` goes as `1/Omega`, and this world turns in 30 hours against Earth's
23.93: `vdiff_lamm = 160*(OMEGA_EARTH/ww)` = 200.5 m, derived in `fluxini` and
overridable by a positive namelist value. Holding 160 m fixed is not the neutral
choice; it asserts a rotation-independent length, which is what Blackadar's
relation denies. **The bracket is `G`**, held equal to Earth's because no run on
this tree has produced a circulation to read one from. `lambda` is linear in `G`,
so a mid-latitude geostrophic wind 20 per cent above Earth's carries
`vdiff_lamm` to 241 m, and that is the sweep this number wants.

---

### 11. `frac_labile_carbon`, the only substrate control on denitrification. RESOLVED

`vendor/lpj-guess/data/ins/global_soiln.ins`, consumed in
`vendor/lpj-guess/modules/somdynam.cpp` at the end of `somfluxes`.

It sets the labile carbon the soil nitrogen transformation operator denitrifies
on. That pool is split by water-filled pore space and becomes the
Michaelis-Menten substrate term multiplying BOTH denitrification steps in
`vendor/lpj-guess/modules/ntransform.cpp`. It is the only substrate control on
denitrification in the model.

**The evidence was an absence, and it was a complete one.** The instruction file
gave it no comment. It appeared in none of
`biosphere/config/ntransform.yaml`'s `instruction_parameters.values`, its
`parser_bounds` or its `calibration.entries`, nor in `somdynam.yaml`, and a grep
for it across `biosphere/` returned nothing.
`vendor/lpj-guess/framework/parameters.cpp` defaults it to 1.0 and the
instruction file halved that with no note.

This is the finding the biosphere sweep exists for. `ntransform.yaml` is an
honest register that names its soil-nitrogen constants `unsourced` where they
are and refuses on them under `--strict`, and a constant this load-bearing sat
outside it entirely: the register's completeness was never checked against the
instruction file it declares, so absence from the register looked exactly like
absence of a constant. The gate now walks both ways and refuses on any numeric
key in `global_soiln.ins` that the register does not carry.

**It was never TUNED, and the class is what settled the repair.** 0.5 is not the
residual of a fit; it is exactly DNDC's humads-path ratio, 0.2 over 0.4, and the
derivation is reachable and has been read, so it was not OPAQUE either. It is
Earth-calibrated, as everything in DNDC is, but IMPLICIT-EARTH was not the
defect that bit first and that class's repair -- find this world's number -- was
neither available nor to the point. The defect was SCOPE.

Li, Frolking and Frolking (1992) p. 9765 states the soluble carbon pool per
path: 60 per cent of the carbon leaving microbial biomass plus 20 per cent of
that leaving humads, in both cases the share recycled into biomass, against CO2
shares of 20 and 40 per cent on the same two paths. So the soluble-to-respired
ratio DNDC implies is 3.0 on the biomass path and 0.5 on the humads path, and
the declared 0.5 was the humads path's ratio applied to all respiration.

**Disposition: SOURCED, and the form is now the paper's.** `somfluxes` computes
labile carbon as one coefficient per pool class on that class's decomposition:
0.6 on `SURFMICRO` and `SOILMICRO`, 0.2 on `SURFHUMUS`, `SLOWSOM` and
`PASSIVESOM`, and nothing on the six litter pools, which is DNDC's own structure
rather than an omission. Both numbers are the paper's own and both calibration
entries now read `agrees`.

Three things about the repair are worth carrying, because each is a shape other
rows in this audit will meet:

- **The ratio was the part that did not cross.** DNDC respires 20 per cent of
  the carbon leaving microbial biomass; CENTURY respires 0.6 on every
  `SURFMICRO` call and `0.85 - 0.68*(clay + silt)` on `SOILMICRO`. Applying
  DNDC's 3.0 to CENTURY's respiration makes the substrate larger than the carbon
  that produced it over most of the texture simplex. Two models that state the
  same process differently share their FRACTIONS, not the ratios between them.
- **The mapping needed a bracket in one place and not in the others.**
  `PASSIVESOM` has no DNDC counterpart at all, and that loose correspondence
  needs no bracket because the pool carries 0.03 to 0.35 per cent of
  decomposition. `SLOWSOM` returns 0.447 of its decomposition to the microbial
  pool where DNDC's humads returns 0.2, and that is a real disagreement about
  one quantity: bracket 0.2 to 0.447, declared on
  `instruction:frac_labile_carbon_humus`.
- **The repair went the other way from the row that opened it.** "Six times low
  on one of two paths" reads as an increase. The litter pools carry 48 to 76 per
  cent of respiration in this model and DNDC gives them nothing directly, so the
  paired form gives 0.70 of mainline's labile carbon on a leaf and root litter
  mix and 0.26 on coarse woody debris alone. A magnitude estimated from the
  coefficients alone had the sign wrong.

It is a declared divergence in two registers, `labile_carbon_paths` in
`biosphere/config/somdynam.yaml` for the source form and
`labile_carbon_microbial_share` in `ntransform.yaml` for the instruction value,
and it is NOT execution-verified: no LPJ-GUESS run this project has made reaches
acceptance (`world-qcse`). `biosphere/notes/soil-nitrogen-transformation-parameterisation.md` carries the correspondence, the bracket, the steady-state sweep and what an accepted-pair arm would settle.

Reading the paper settled a second thing the row named. Li's table 7 gives `Kc`
= 0.017 kg C/m3 and `Kn` = 0.083 kg N/m3 and attributes BOTH to Shah and Coulman
(1978) rather than measuring them, so the chain to a measurement is one paper
longer than this project recorded -- and it states the units exactly as Xu-Ri
does, per cubic metre of an unnamed volume. The `michaelis_menten_divisor`
question is therefore NOT closed by it. Shah and Coulman (1978) is what settles
it, in favour of the soil-solution reading the operator already runs, and the
magnitude argument from DNDC's own soluble carbon pool that this row once
offered for the per-cubic-metre-of-SOIL branch is withdrawn: it rested on an
assumed bulk density and on a coincidence of numerals between two different
volumes.

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

**Disposition: DERIVED, and the derived object is a FORM rather than a
constant.** Kessler (1969) was read on 2026-08-25 and the derivation this row
called for closes. His table 4 gives rain evaporation for a Marshall-Palmer
distribution as `dM/dt = k3 N0^(7/20) m M^(13/20)` in g m-3 s-1, with the
evaporation coefficient `k3` = 1.93e-6, `N0` = 1e7 m-4 so `N0^(7/20)` = 282,
`m` the saturation deficit and `M` the rain water content, both in g/m3; his
mean volume-weighted fall speed is `38.8 N0^(-1/8) M^(1/8)`, i.e.
`5.17 M^0.125` m/s. Integrated over the layer the four re-evaporation sites act
on, **the layer depth cancels** and

    gamma = 5.44e-4 * M^0.65 * deltsec2,   M = P / V

with `P` the layer's precipitation flux. Three consequences, none of which a
constant can carry. `gamma` is **proportional to the timestep**, which the
declaration already inferred from the form and which now has a derivation under
it. It goes as **P^0.578**, so it is not one number per planet, let alone per
rung. And it carries **ga^(-0.289)** through the drop terminal speed, which is
0.926 of its Earth value here.

**The magnitude disagrees with the declared value and that is the finding.** At
this project's step, `deltsec2` = 3600 s, the derived value is 0.039 at 0.5
mm/day, 0.10 at 3 mm/day and 0.22 at 10 mm/day. **0.01 is below the whole of
that span**: it is the Kessler value at 0.048 mm/day, sixty times below a
global-mean precipitation rate. Using a grid-mean flux understates `M` and so
understates `gamma`, which makes those figures floors rather than estimates.

The derivation is in the declaration. The value is left at 0.01 because
replacing it is a change of FORM in `rainmod.f90` -- the four blocks compute
`gamma` per level from `zprl`, `zprc`, `zprsl` and `zprsc`, which they already
hold -- and that is `world-trs3`'s remaining half. Kessler's own caveats travel
with it: the single-drop fit is accurate to about 40 per cent, a constant `N0`
misrepresents evaporation because the process depletes small drops
preferentially, and the rate is for standard air density with no altitude
variation.

---

### 13. `zcca = 0.245`, `zccb = 0.125`, and the `rcrit` floor of 0.85

`rainmod.f90`'s `mkclouds` and `rainini`. The convective cloud fraction as
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

**Disposition: IRREDUCIBLE as scalars; DERIVED as a form for `rcrit`, and
IRREDUCIBLE as a form for the convective pair.** Neither number can be sourced,
because what they encode is a subgrid distribution and the subgrid scale is a
property of the mesh rather than of the world. `world-khn` asked the resolution
question and closed it by declaring the anchor with no numeric changed, which is
the right answer to that question and leaves this one. The two halves separate
under it, and they separate for different reasons.

**`rcrit` takes a derived NLAT term, and it goes on the WIDTH.** What `rcrit`
encodes is the width of the subgrid humidity distribution: cloud begins when the
moist tail reaches saturation, which puts onset a distance proportional to the
width below saturation, so the quantity with a resolution dependence is
`(1 - rcrit)` and never `rcrit`. Specific humidity is a passive scalar, and in
the inertial-convective subrange Kolmogorov-Obukhov-Corrsin gives its variance
across a separation `L` as `L^(2/3)`; the standard deviation therefore goes as
`L^(1/3)`, and a cell's width is its grid spacing, which goes as `1/NLAT`. So
`(1 - rcrit)` is proportional to `NLAT^(-1/3)`, integrated as
`1 - rcrit(jlev) = (1 - rcrit_T21(jlev))*(32/NLAT)^(1/3)` through a new
`rcritwidth` key that derives from NLAT unless given a positive value.

`rcritmod` is the wrong operator for it and that is worth recording: it
multiplies `rcrit`, and at the top and bottom levels the sigma limb already puts
`rcrit` near 0.95, so the factor giving the right floor at T85 drives those
levels past 1.0 where `(rh-rcrit)/(1-rcrit)` divides by zero. Scaling the width
is exact at every level and bounded below 1 by construction. `rcritmod` and
`rcritslope` keep their meanings and apply on top, as the sweep handles.

**T21 does not move and the higher rungs do.** The anchor is NLAT 32, so the
factor is exactly 1 at T21 and every T21 answer is unchanged. On the 0.85 floor
it is 1.024 at T31, 1.036 at T42, 1.054 at T63 and 1.065 at T85, taking `rcrit`
from 0.85 to 0.869, 0.881, 0.896 and 0.906, and `1/(1-rcrit)^2` from 44 to 58,
71, 92 and 112. Cloud starts later in a smaller cell, which is the direction the
argument requires.

**The convective pair does NOT take that scaling, and the reason is that they
are not the same kind of quantity.** `rcrit` encodes the width of a subgrid
scalar distribution, which the Corrsin argument gives an exponent for; `zcca`
and `zccb` encode how a convecting AREA dilutes into a cell, and the exponent
there is set by the unresolved convective area fraction. A storm small against
the cell dilutes as `L^-2` and one that fills the cell does not dilute at all, so
the correction that would hold `zcc` fixed,
`zcca -> zcca + 2*zccb*ln(NLAT_ref/NLAT)`, is exact only in the first limit and
wrong by that whole term in the second. The model carries no convective area
fraction, so nothing in this scheme can tell the two limits apart and no exponent
between them is derivable from what it holds. Writing one in would be a fit
wearing a derivation's clothes. What would settle it is a convective area
fraction carried alongside the rate, which is a scheme change and not a
constant.

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

The pedology sweep's own headline: the tuned values in that component were
concentrated in `pedology/config/pedogenesis.yaml`'s `regolith:`, `texture:`
and `ph:` sections, and those were the only blocks in the file that did not
ship brackets. Every other declared-without-a-source value in either that
component or hydrography (`catena`, `water`, `andisol`, `bedrock_water`, the
loess threshold, `f_grad`, `lambda_m`) ships a bracket and an instruction not
to tune within it.

| where | what was missing | magnitude | disposition |
| --- | --- | --- | --- |
| `texture.clay_conversion` | no comment on the value, no source, no bracket; numerically identical to `clay_yield`, which WAS fitted to a SoilGrids slope of 0.543 | near-linear at this world's weathering intensity: a 10 per cent move is +9.2 per cent converted clay | BRACKETED on the e-folding intensity, 0.29 to 1.00, and the near-degeneracy with `clay_yield` recorded |
| `texture.clay_yield` | fully admitted as an Earth fit over sixteen type localities, and applied BELOW the range it was calibrated in, where the file records the bias is largest | linear multiplier; documented residual bias of +0.084 clay fraction at this world's mean | UNCHANGED. It is disclosed, and disclosure is the standard this row asks for |
| `regolith.dry_erosion_baseline` | the mechanism was argued and the value was not attached to it; no source, no bracket | 21 per cent of the moisture term at the land mean, 100 per cent of it over the arid fraction of land, where depth goes as its reciprocal | SOURCED from Portenga and Bierman's arid basins, and bracketed from their stated uncertainties |
| the `ph:` block: six parent values, `leaching_slope`, `endorheic_alkalinity_bonus` | no comment on ANY value, no citation, no bracket, nothing in README, notes or references | the slope moves land-mean soil pH by 0.04 per 10 per cent and wet-cell pH by 0.13; the bonus is a flat 0.08 over the closed-basin fraction that the brine and duricrust rules key on | REBUILT. The block's FORM was refuted by its own source and replaced; BOTH BUFFERS ARE NOW DERIVED at this world's pCO2, the lithology contrast is SOURCED as a base-cation supply and no longer stated as a pH at all, and the closed basin is now the SUMP's buffer inside that relaxation rather than a flat offset on the whole basin floor, and what remains declared is `gibbsite_buffer_ph` (implicit-Earth, bracketed), `soda_buffer_ph` (implicit-Earth, bracketed, emitted so the source is stated once) and `leaching_slope` (bracketed one-sided) |
| `texture.sand_to_silt_loss_ratio` | no source; and the mechanism sentence beside it had the surface-area argument backwards, which is why the number did not follow from the reasoning printed next to it | redistributes 2 to 3 per cent between sand and silt at fixed clay; under 1 per cent on PAWC | BRACKETED over both mechanisms, sentence corrected, and the value moved to the one point in the bracket that carries a statement |

**Disposition: RESOLVED, and one of the five turned out not to be a bracketing
job at all.** The arithmetic behind every end is in
`pedology/notes/pedogenesis-value-provenance.md`; what follows is what changed
and what reading the sources found.

**`dry_erosion_baseline` was sourced, not bracketed, and it was low by a factor
of about six.** The moisture term is `runoff / reference_runoff +
dry_erosion_baseline`, so the key IS the ratio of denudation at zero runoff to
denudation at Earth's reference runoff, as `b / (1 + b)` -- and that ratio is
something the compilation already read for row 3 reports. Portenga and Bierman
put arid drainage basins at 100 +/- 17.3 m/Myr, the slowest of their climate
zones, against a global drainage-basin mean of 218 +/- 35 m/Myr. The ratio is
0.459 and inverts to 0.85, with the two stated uncertainties giving a bracket
of 0.49 to 1.79. The old 0.15 was a suppression of dry denudation to 13 per
cent of the reference rate with nothing behind it. **The bracket reaching a
nearly flat moisture term is the same paper's doing**: its basin denudation has
no significant bivariate correlation with mean annual precipitation at all,
mean basin slope being the regressor. This is the second time that compilation
has answered a question this audit did not fetch it for.

**The pH block's alkaline end is DERIVED, at this world's atmosphere rather
than Earth's.** Slessarev, Lin, Bingham, Johnson, Dai, Schimel and Chadwick
(2016), *Water balance creates a threshold in soil pH at the global scale*,
Nature 540, 567-569, `10.1038/nature20139`, give the calcite equilibrium in
their Methods as a quartic in hydrogen ion activity at a stated CO2 partial
pressure. That is a derivation to do rather than a number to copy, which is the
shape Kessler took in row 12. Solved at their laboratory pressure it returns
8.236, and at the pre-1977 pressure they quote a check for it returns 8.248;
both round to what the paper prints, so the transcription is theirs and not
something plausible assembled beside it. Solved at Vesper's declared pCO2 it
gives 8.16 for the calcite buffer and 5.59 for water carrying no alkalinity.
The first is now TWO keys with two roles, `parent_by_category` and
`calcite_buffer_ph`: the pH of a solution on carbonate rock, and the buffer a
soil that exports nothing reaches whatever its rock was. The two values together
bound the fresh solution each declared base-cation supply implies, which is the
one pH statement the block still makes about a rock. The derivation travels to
another planet, which is what separates it from adopting 8.2.

**And the same paper refuted the block's FORM, which is more than this row
asked for, and the form was then replaced.** Slessarev sample 20,000 subsoil pH
measurements out of 60,291 and find the distribution BIMODAL: a calcite mode
near 8.2 and a gibbsite mode near 5.1, with an abrupt step where mean annual
precipitation begins to exceed potential evapotranspiration, and the neutral
range uncommon. NEITHER MODE CARRIES THE ROCK. `build_soil.py:soil_ph` used to
start every parent at its own value and take one `leaching_slope` down from it,
a line with no asymptote at either end, so no slope landed the wet cells on one
buffer and no slope held the dry ones on the other. It now runs the leaching
index BETWEEN the two buffers,

    ph = G + (D - G) * exp(-leaching_slope * L / u)

with `D` the buffer the cell reaches at zero export and `u` the parent's
base-cation supply as a fraction of a calcite-saturated soil's. Both ends are
parent-independent, which is the observation, and the lithology contrast
survives at finite leaching in `u`, which is where the paper's own carbonate
deviation lives.

**The closed basin is the third buffer, and it is inside that relaxation.** A
closed basin evaporates past calcite saturation and the brine that does so ends
up in the SUMP -- Orogen zones the basin and puts the salt crust in the lowest
ground, giving the margins the clastic load -- so `D` is raised toward
`soda_buffer_ph` by the salt crust's share of the cell. It was an additive
offset on the endorheic fraction, eleven times that area, of which 91.75 per
cent drains at a median 11.55 mm/yr, and unconditional in `L`, so a fully
leached floor cell landed on the gibbsite buffer and then had 0.8 added to it.
A solute stays where the drainage does not export it, which is the same
statement `exp(-sL/u)` already makes, so the closed basin decays the same way.
`soda_buffer_ph` is Helvaci (2019)'s closed-basin lake water, low end declared
because the top of that range is the most evaporated brine in a hydrothermally
fed basin; IMPLICIT-EARTH, because a soda brine's pH is set by its own
alkalinity at its own pCO2 and the repair is a trona saturation solve this
project cannot do today.

**`u` IS SOURCED, and it is not a pH.** It sits in the denominator of an
exponent whose numerator is a leaching index, so it is a ratio of a supply flux
to an export flux; encoding it as `(parent - G)/(C - G)` from a declared parent
pH was an affine map between quantities a logarithm separates, and it could not
express a supply below 0.159 at all. Meybeck (1987) Table 2C, already extracted
in this tree, and GEM-CO2, Amiotte Suchet and Probst (1995) as
`vendor/cgenie/genie-rokgem` implements it, both measure the quantity per
lithology on a common basis and agree to better than a factor of two on every
category. Four supplies are declared from them, IMPLICIT-EARTH by the same test
that convicts the gibbsite buffer -- the mineral stoichiometry travels, the
temperature and soil CO2 the measurements were made at do not -- with brackets
that are the spread between the two Earth compilations and not the transfer.

The bracket on `leaching_slope` is one requirement evaluated at three points of
the wettest quartile, and it is ONE-SIDED: all three ends are readings of an
upper bound, and the finding that would supply a lower one stops being a bound
under a two-buffer form rather than failing. Re-run against the sourced supply
set, the three readings are unchanged and the lower bound is still absent.

**`sand_to_silt_loss_ratio` is not derivable and the bracket says so.** The two
mechanisms point opposite ways: specific surface area goes as one over grain
diameter, so silt dissolves faster per unit mass and the USDA class limits put
that at a ratio of 0.032; while weathering as a size cascade has sand leaving
its own class faster than dissolution implies with silt replenished out of it,
which is unbounded above. The bracket is that whole span and the value now sits
at the one point inside it carrying a statement, loss in proportion to
abundance. That a bracket three orders of magnitude wide still moves PAWC by
under one per cent is the result, not a weakness of the bracket.

**What is left open, and it is not a value.** Row 4's joint constraint between
`maximum_depth_m` and `erosion_coefficient_per_relief_m` is a relation between
two keys rather than a value inside an interval, so no bracket expresses it and
nothing checks it. The single-key half of this is closed:
`scripts/smoke_test.py:DECLARED_BRACKETS` now registers every bracket in the
file and refuses one it does not name, which is what makes a bracket added here
a bracket something reads.

---

### 16. The LPJ-GUESS constants whose own source says they were invented

| where | value | the source's own words |
| --- | --- | --- |
| `vendor/lpj-guess/framework/guess.h:2552`, comment at `:2546` | `frac_maxtomin = 0.9` | `// Tighter C:N ratio range for roots and sapwood: picked out thin air.` |
| `vendor/lpj-guess/modules/soil.cpp:111` | `sompool[PASSIVESOM].ntoc = 1.0/9.0` | `// passive has a fixed value (why? passive SOM should also vary.)` |
| `vendor/lpj-guess/modules/somdynam.cpp:50` | `TAU_LITTER = 2.85` | `// Thonicke, Sitch, pers comm, 26/11/01` |

The first two are live under `ifnlim 1`. `frac_maxtomin` sets the width of the
fine-root and sapwood C:N windows, and this project has already repaired their
PLACEMENT from Friend et al. (1997) Table 4 under `world-xms4` while leaving the
width unsourced for both elements; `guess.h:313` names that as a whole-model gap.
`ntoc` sets how much nitrogen is locked in a pool with a turnover of order 1400
years, so it sets steady-state mineral nitrogen after `equilsom`.

`TAU_LITTER` is on the `ifcentury 0` path and this project sets `ifcentury 1`, so
it is dead.

**Disposition.** `ntoc` is REPLACEABLE NOW: Parton, Stewart and Cole (1988) is
held and read, and it is the paper this project already used to settle the
phosphorus analogue, `PMASS_SAT` and `SURFHUMUS`. The nitrogen side was simply
not given the same treatment. `frac_maxtomin` is REPLACEABLE, SOURCE NEEDED and
the source is a C:N range rather than a ratio, so it is a literature question
rather than a fetch.

**The GLOBFIRM burn-probability floor was audited here and is NOT tuned.** It
belongs in this section's history because `// c.f. LPJF` reads exactly like the
three comments above it, and the reading was wrong. The floor is IMPLICIT-EARTH,
which is a different defect with a different repair, and misfiling it here
taught the wrong lesson about the cleanest worked example of the other class
this project has.

Three findings settle it, and each rules out one class.

- NOT TUNED. A fit residual does not arrive at exactly 1e-3 in two independent
  codebases. LPJmL's `src/soil/fire_prob.c`, vendored at `vendor/lpjml`, returns
  0.001 for any `fire_frac` below 0.001, and no comment in either implementation
  claims a fit.
- NOT OPAQUE. The repair the opaque class prescribes is to walk the chain, and
  the chain is walked to its end. `// c.f. LPJF` points at the FORTRAN LPJ, whose
  fire module is Thonicke, Venevsky, Sitch and Cramer (2001). That paper states
  no minimum burned fraction anywhere. It states the opposite: its Eqn 5, `A(s) =
  s * f(s)`, is introduced with "A is zero when fire conditions were absent
  during the year", and its own reporting stops at "more than 900 years" for
  regions unsuitable to carry fire. Sitch et al. (2003), the LPJ-DGVM
  description, documents the fire terms as empirical global-Earth choices and
  states no floor either. Two implementations and two papers; no derivation at
  the end of the chain.
- IMPLICIT-EARTH, and the model's own source is what says so. `commonoutput.cpp`
  writes `firert_gridcell += 1000.0` under `// Set a limit of 1000 years` for any
  patch whose `fireprob` is below 0.001, so 0.001 per year is a REPORTING CAP on
  the fire return interval -- the same diagnostic ceiling Thonicke's Fig. 6
  reports at, one round cutoff on where a diagnostic stops resolving. The code
  applied that cap to the model STATE. And the cap is one per EARTH year: LPJF
  and LPJmL both run `NDAYYEAR 365` against this model's 183-day year, declared
  as `VESPER_YEAR_LENGTH_DAYS` in `framework/vesper.h`, so carrying the numeral
  across delivered a whole Earth year of the floor every 183 days. The
  neighbouring quantity of the same class IS converted for exactly this reason:
  `distinterval` reaches the model as 199.594 simulation years against an Earth
  calibration of 100.

The repair is not the rescaling the implicit-Earth diagnosis usually implies.
Rescaling would repair the unit of a number with no derivation, which is a
tidier copy of the wrong thing, and the quantity being rescaled is a diagnostic
cap rather than a mechanism. The line is deleted, under `world-v5j1`. The
deletion is registered as a declared divergence from mainline LPJ-GUESS 4.1.1 in
`biosphere/config/fire.yaml`, `biosphere/scripts/fire_gate.py` enforces the
absence, and the argument is `biosphere/notes/fire-model-audit.md` section 8.

What deleting the floor does NOT dispose of is Thonicke's Eqn 8 itself. Its four
coefficients -- 0.45, 2.83, 2.96, 1.04 -- are a least-mean-square non-linear
regression against observed area burnt in Portugal, southern California and
Kakadu National Park, which is a published fit with inspectable provenance and
an Earth domain rather than a tuning. The curve is the operator this project
runs, and it is a separate question from the floor that was sitting on top of
it.

---

### 17. `ntransform.yaml`'s six `unsourced` entries, already registered

`biosphere/config/ntransform.yaml`, the `calibration.entries` whose verdict is
`unsourced`: `instruction:f_denitri_gas_max`, `instruction:f_denitri_max`,
`function:wet_fraction`, `form:wet_fraction_shape`,
`function:nitrification_activity` and `form:denitrification_wfps_threshold`.
Named rather than cited by line, because the register is appended to and every
line number in it goes stale silently. `ntransform_gate.py` prints the set on
every invocation. Recorded here so the survey is complete and so nobody re-derives them: this
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

**`scripts/error_budget.py`'s `DEFAULT_ATTENUATION`** is the one that looked like
a fudge and was not: it was a declared 0.5 with the file saying the honest claim
was a factor of two rather than a number. It is now measured on four paired arms
at 0.38 with a span of 0.29 to 0.48, so it has left this audit's subject matter
entirely; `notes/audits/albedo-attenuation.md` is the finding.

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
| 3, 4. the regolith depth level and its asymptote | `world-qs63`; the degeneracy is gone and the level is bracketed, and `maximum_depth_m` needs a different source than this row named |
| 5. the carve size floor | `world-7vj6`; the floor is swept and reported beside the verdict, and `C` is declared irreducible |
| 6. `eddy_wind_m_s` | `world-9y07` |
| 7. `acllwr` | fixed |
| 8. `th2oc` | `world-2esd` |
| 9. the `dz0land` anchor | `world-u8ds` |
| 10. `vdiff_lamm` and the Louis fives | `world-x6q8`; the fives are answered and irreducible, the mixing length waits on Blackadar (1962), and `vdiff_c`'s dropped roughness dependence is `world-awm5` |
| 11. `frac_labile_carbon` | `world-vyvn`; the two-path form the primary states is implemented and declared, both coefficients are the paper's, and the humified path carries a bracket. Execution verification waits on `world-qcse` |
| 12. `gamma` | `world-trs3`; the form is derived and in the declaration, and implementing it in `rainmod.f90` is what is left |
| 13. `zcca`, `zccb`, `rcrit` | `world-o12h`; `world-khn` closed the resolution half by declaring it |
| 14. `a1` in the sigma quartic | none: irreducible, and a convergence sweep across it is a measurement rather than a fix |
| 15. the unbracketed pedogenesis values | `world-9ctm` |
| 16. `ntoc`, `frac_maxtomin` | `world-xfif`. The GLOBFIRM burn-probability floor left this class rather than being repaired in it: it is implicit-Earth, `world-v5j1` deleted it, and `biosphere/config/fire.yaml` registers the deletion |
| 17. `ntransform.yaml`'s six `unsourced` entries | already registered; the gate refuses on them |
| 18. the dormant knobs | `world-9g8p` |
| 19. the cgenie tier | `world-u9kg` |

## The five papers, and what reading them changed

All five are read. Four of the five dispositions they were fetched for moved,
and **four of the five expectations recorded against them were wrong** -- which
is what an expectation written against an unread paper is worth, and the reason
this section says which.

| paper | what it was expected to give | what it gave |
| --- | --- | --- |
| Portenga and Bierman (2011) `10.1130/G111A.1` | global cosmogenic denudation BY LITHOLOGY, and with Heimsath an absolute production-to-erosion ratio for rows 3 and 4 | by-lithology denudation for OUTCROPS only, which carry no regolith, and none at all for drainage basins; no relief-to-erosion relation of any kind, mean basin slope being the regressor. It sourced the erosion term's FORM and, with Heimsath, bracketed the level at 0.15 to 1.30 m with a regime in it where no soil is possible. Row 4's expectation was wrong outright: it measures denudation, not a weathering-front reach |
| Li et al. (1992) `10.1029/92JD00509` | the value of `frac_labile_carbon`, and the `michaelis_menten_divisor` volume, together | the first, by supplying a FORM where the row expected a constant. DNDC's soluble carbon is per-path, 0.6 of biomass decomposition and 0.2 of humads, so a single fraction of total respiration is not a quantity the paper forms and the operator now carries both coefficients. On the divisor it repeats Xu-Ri's silence and pushes the chain one paper further, to Shah and Coulman (1978) |
| Louis (1979) `10.1007/BF00117978` | whether the three fives are Louis's set or ECHAM's re-fit | ECHAM's re-fit, definitively: Louis's own set is `b = 2b' = 9.4` with `c` varying as `sqrt(z/z0)`, and he calls his values uncertain. It also removed a proposed derivation, since his `lambda` is 100 m and declared adjustable. The dropped roughness dependence it was read as exposing in `vdiff_c` is not there -- see row 10; the reading of the paper was right and the reading of the code beside it was not |
| Kessler (1969) `10.1007/978-1-935704-36-2` | the re-evaporation FORM, as a derivation to do rather than a number to copy | exactly that, and the derivation closes: `gamma = 5.44e-4 M^0.65 deltsec2` with the layer depth cancelling, proportional to the step, going as `P^0.578` and carrying `ga^(-0.289)`. The declared 0.01 is below the whole span the derived form reaches |
| Blackadar (1962) `10.1029/JZ067i008p03095` | a rule for the asymptotic mixing length that carries to another rotation rate | exactly that: eq. 25 is `lambda = 0.00027 G/f`, dimensionally argued with `z0` ruled out of the free atmosphere. The VALUE does not transfer -- 0.00027 is matched to one observed wind deflection, he offers `u*/f` as equally acceptable, and his 26 m is a boundary-layer asymptote against ECHAM's whole-column 160 m -- but the `1/Omega` SCALING does, taking `vdiff_lamm` to 200.5 m on a 30-hour rotator with `G` left as the bracket |

Both papers named as needed by the readings above are now HELD AND READ, and
the shape is the ordinary one for this project: a primary names its own primary,
and the chain is walked rather than stopped at the secondary.

Shah and Coulman (1978) `10.1002/bit.260200105` settled the Michaelis-Menten
volume. It was closed access on every automated route and was SUPPLIED BY THE
USER. Its constants are measured in a continuous stirred tank reactor on
dissolved substrates with no soil in the experiment at all, so the volume is the
SOLUTION'S and `ntransform.cpp` dividing by the depth of available water is
correct. The factor of 6.7 to 30 on denitrification closed in favour of the
branch already running.

Shangguan et al. (2017) `10.1002/2016MS000686` sourced `maximum_depth_m`. It is
read and NOT HELD AS A PDF, which is a third state worth distinguishing from
both "fetched" and "closed": the article is gold open access and the publisher's
own OA URL returns HTTP 403 to every automated route tried, while the full text
is reachable through the HTML article page. The numbers were transcribed from it
and cross-checked against a second independent read. `references/INDEX.md`
records that state on the row rather than implying a file exists.

