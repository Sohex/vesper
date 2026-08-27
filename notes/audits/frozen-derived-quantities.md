# Frozen derived quantities across the tree

Audited 2026-08-27. The question: which numbers in this project were COMPUTED by
something, then written down somewhere that has no way to learn they moved?

A sweep of the whole tree against the corollary `docs/src/pipeline/loops.md`
states as "a field update has to be able to carry". It reports what is frozen and
where. It repairs nothing: the disposition of each instance is a separate
decision and several are contentious.

## The class, and the test

A FROZEN DERIVED QUANTITY is a value that some other artifact, script or
measurement computes, copied into a place with no path back to the computation.

The test is one question, asked of any number: **WHAT RE-RUNS WHEN THIS
CHANGES?** If the answer is nothing, the quantity is frozen. It will drift from
the thing it describes without anything objecting, because a stale number is a
perfectly ordinary number and the artifact carrying it is a real artifact.

The test is about the PATH BACK, not about correctness today. Most instances
below currently agree with what produces them. That agreement is not evidence of
anything: it is what a number looks like the day before it drifts, and the ones
already stale looked exactly like it the day before they did.

## The boundary, which is most of the value

Four things are NOT in the class, and getting each right is what keeps the class
narrow enough to act on.

**A DECISION is not a derived quantity.** A design preference, a declared
threshold, a bracket fixed ahead of the runs it judges. These are INPUTS: they
are supposed to be written down and to stay put. `cold_extreme_cap`,
`warm_ceiling_c`, the convergence thresholds, `filter_kappa`, `filter_power`,
`lib/run_lengths.py`'s `PRODUCTION_SPAN_TAU_MULTIPLE` and `SETTLING_RESIDUAL_K`,
`lib/gridding.py`'s `COUPLING_OCEAN_FRACTION_LIMIT`, the 32 MB per-die target,
`aeolian/scripts/dust_intensity_levers.py`'s two tabulation criteria and
`maps/render_projections.py`'s `cut_max` and `gross_max` are decisions. So is
`scripts/error_budget.py`'s `DEFAULT_ATTENUATION`, which declares itself a
conservative factor of two rather than a measurement.

A measured figure USED AS a criterion is still a decision, provided it was fixed
before the results it judges. `dust_intensity_levers.py` states that in its own
header and is the model for how to say it.

**A PHYSICAL CONSTANT or a paper-sourced value is not in the class.** Nothing in
this tree computes it, so it cannot drift here. `lib/lapse.py`'s molar masses and
specific heats, `lib/orbit.py`'s Earth year lengths, `lib/snow.py`'s
`RHOICE_F2021` and `ADOPTED_TEMPERATURE_K`, `analysis/soil_thermal_inertia.py`'s
Horai and Farouki coefficients, `config/planet.yaml`'s `ozone_scale` and
`cloud_water_reference_kg_m3`, every `phosphorus_ppm` row in
`pedology/config/pedogenesis.yaml`, the whole emission block of
`aeolian/config/dust.yaml`, every grounded rule in `minerals/config/`, and
`maps/README.md`'s polyhedral geometry table are out.

**A value already inside a loop with a check is the target shape, not a
defect.** `baseline_flux_earth` is re-derived on every new terrain. The land
column contract re-derives its medians and refuses when the declaration
disagrees. Both are listed below, because the next person needs the pattern more
than the list.

**A dated MEASUREMENT in a note under `notes/` is a historical record and is
SUPPOSED to be frozen.** `docs/src/practice/conventions.md` says so. Such a note
enters this class only when something CONSUMES its number as a live input, or
when a document elsewhere restates it AS CURRENT.

## The three neighbouring classes, and how this one differs

This class is orthogonal to the three the vocabulary already carries. A frozen
value's derivation is sound, inspectable, and on this planet. What is missing is
only the path back to it.

| class | what is missing | the audit |
| --- | --- | --- |
| tuned | any derivation at all; the number IS the residual of a fit | `tuned-values.md` |
| opaque | the derivation exists somewhere a reader of the code cannot reach | `opaque-constants.md` |
| implicit-Earth | the derivation is sound and is for the WRONG PLANET | `inherited-earth-constants.md`, `ocean-tier-implicit-earth.md` |
| frozen | nothing; the derivation is sound and reachable, and the COPY cannot learn it moved | this one |

A value can sit in two at once. Nothing here is a reason to reopen a value's
classification in another.

## The two dispositions

Either the consumer READS THE EMITTED VALUE, so there is nothing to freeze; or
the declaration lives INSIDE A LOOP that re-derives it, with a check that fires
when the declaration and the derivation disagree. A declared constant with
neither is a drift with a timer on it.

This tree already contains good instances of both, and they are the argument that
the repair is cheap rather than structural.

### Disposition one: state no number at all

`config/planet.yaml`'s `precip_reevaporation_gamma`, `asymptotic_mixing_length_m`
and `cloud_fraction_subgrid_width` are set to the string `derived`. The model
reads each key on its sign, computes the value from what it already holds, and
prints what it chose; `run_exoplasim.py` reads that print back onto the run
manifest under `derived_model_constants`, and `check_consistency.py` refuses a
run that staged a sentinel and recorded no derived value. The config's own
comment states the reasoning: a number written there would be a second statement
of a derivation the model owns, silently staleable against it and against a
change of rung.

`aeolian/config/volcanic_sulfate.yaml` and
`minerals/config/downstream_prospectivity.yaml` do the same in configuration
rather than in code, each refusing to copy a value another config owns and saying
so. `biosphere/config/snow_thermal.yaml` refuses to restate a polynomial's
coefficients in the same terms.

This is the cheapest disposition available and it is the one to reach for first.

### Disposition two: declare, re-derive, and refuse on disagreement

| the declaration | what re-derives it | where the check fires |
| --- | --- | --- |
| `pedology/config/land_column_properties.yaml` `thermal.saturation_mapping` and `thermal.endpoints` | `pedology/scripts/land_column_properties.py`, from the states the contract itself emits | the same script, at the declared `tolerance` |
| `lib/sensitivity.py` `SLOPE_K_PER_FLUX_RATIO` and `SLOPE_BRACKET_RUNS` | `sensitivity.verify()`, reading the named runs' asymptotes out of `exoplasim/runs/INDEX.json` | `check_consistency.py`, and `scripts/error_budget.py` |
| `config/planet.yaml` `h2o_sw_weight`, `co2_sw_weight`, `h2o_sw_level`, `cloud_absorption_scale` | `exoplasim/scripts/shortwave_band_weights.py` and `cloud_band_weight.py`, into named JSON | `check_consistency.py`, at the tolerance the config's own written precision implies |
| `config/planet.yaml` `model.hyperdiffusion.timescales_days` | the rule `pi*a/(NTRU*eddy_wind)` | `check_consistency.py` |
| `lib/rungs.py` `STABILITY_CEILING_MINUTES` | `check_stability_ceilings()`, from `exoplasim/analysis/stability_probe.json` | `check_consistency.py` and `smoke_test.py` |
| `lib/rungs.py` `RUNGS`, restated in CMake, in the installed package and in a Fortran comment | `check_restatements()` | `smoke_test.py` |
| `lib/rungs.py` `ESCALATION_ROUTE`, restated in `sequencing.md` section D | `check_timestep_restatements()` | `smoke_test.py` |
| `lib/stellar.py` `SOLAR_PARTITION` | reproduced from the spectrum in the same module | `check_consistency.py` |
| `exoplasim/config/ocean_tier.yaml`'s per-constant `declaration:` lines | `exoplasim/scripts/ocean_tier_gate.py`, against the Fortran source | the gate, with its own fixtures |
| `vendor/exoplasim/.../landmod.f90` `sicediff`, `sicecap` | `analysis/ice_properties.py`, which holds the compiled literals to what it computes | that script, on every run |
| `analysis/ice_albedo.py` `RECORDED` | reproduced from the vendored source | that script, at 1e-9 |
| `biosphere/config/wetlands.yaml` `production_ratio` | `biosphere/scripts/wetland_gate.py`, re-extracting both constants from `modules/soil.h` | the gate, at one per cent |
| `biosphere/config/mineral_reactivity.yaml` `literals` | `mineral_reactivity_gate.py`, against `somdynam.cpp`, with the closure evaluated on the live soil map | the gate |
| `pedology/config/outgassing.yaml` `source_partition` fractions | `outgassing_gate.py`, from each source's flux over the assembled total | the gate, at a declared tolerance |
| `biosphere/scripts/acclimation_gate.py` and `wetland_gate.py` constants | both parse `vendor/lpj-guess` source at run time rather than copying | the gates |

The `check_consistency.py` band-weight table is the most reusable of these: four
rows of `(config key, artifact, node path, generator)`, the generator named so a
failure says what to re-run, and a tolerance that is arithmetic on the config's
own written precision rather than a number chosen to make today's values pass.
Its own comment states the general problem it was built for. Most of the
config-side instances below are one row in that table away from being closed.

## What the sweep found

Sixty-four instances. No area came back empty. Within areas, these returned
nothing: `lib/climatology.py`, `lib/remap.py`, `lib/gridding.py`, `lib/orbit.py`,
`lib/surface_classes.py`, `lib/sea_water.py` and `lib/provenance.py`;
`vendor/lpjml/` and `vendor/cgenie/`, because nothing in this project reads
either yet; and `vendor/orogen/js/lithology.js` and `basins.js` apart from one
ambiguous sentence recorded below.

| area | instances | on the path to a baseline run |
| --- | --- | --- |
| `config/planet.yaml` | 14 | 12 |
| `exoplasim/` scripts and config | 9 | 5 |
| `lib/` | 6 | 5 |
| `scripts/` | 5 | 1 |
| `analysis/` | 4 | 0 |
| `pedology/` | 6 | 3 |
| `hydrography/` | 7 | 5 |
| `biosphere/` | 3 | 1 |
| `aeolian/` | 6 | 0 |
| `minerals/` | 4 | 0 |
| `maps/` | 4 | 0 |
| `vendor/` fork comments | 6 | 0 |
| `docs/src/` | 8 | 0 |

Eleven are DEMONSTRABLY STALE against the artifact that produces them today. The
rest currently agree, which the test says nothing about. A twelfth, the land
column contract's saturation mapping, was stale when the sweep began and is the
one instance that was repaired while it ran -- because it is the one that has a
check.

## Ranked by blast radius

Blast radius is what reads the quantity and how far an error travels. The tiers
are ordered by that, not by how wrong anything is. Tiers 1 to 3 are on the path
to a baseline run; tiers 4 and 5 are not.

### Tier 1: it sizes every run on the ladder

These decide how many orbits get bought and when a run is declared finished. An
error is multiplied by every rung of the resolution ladder.

| quantity | where | what computes it | check |
| --- | --- | --- | --- |
| `NOMINAL_TAU_ORBITS`, `NOMINAL_ORBIT_SCATTER_K` | `exoplasim/scripts/assess_convergence.py` | `lib/autocorrelation.py` over a run's per-orbit mean surface temperature; each run's convergence artifact carries the result | none |
| `NOMINAL_RELAXATION_ORBITS` | `exoplasim/scripts/assess_convergence.py` | `relaxation_orbits()` in the same file, from the run's slab capacity and `lib/sensitivity.py`'s feedback | none |
| `TAU_MEMORY_ORBITS_BRACKET` | `lib/run_lengths.py` | the same autocorrelation reading; a second statement of `NOMINAL_TAU_ORBITS` at a different precision | `smoke_test.py` checks only that it is an ordered pair of positive times |
| `TAU_RELAXATION_ORBITS_BRACKET` | `lib/run_lengths.py` | `exoplasim/scripts/check_relaxation_ceiling.py`, which writes `relaxation_ceiling.json` beside the convergence reports | none |
| `COMMISSIONING_EVIDENCE` orbit counts and verdicts | `lib/rungs.py` | the runs themselves, through `exoplasim/runs/INDEX.json` and `archive/runs/*/INDEX_ENTRY.json` | `_check_ceilings()` and `_check_route()` validate the rungs and the ordering; nothing re-reads an orbit count or a verdict |
| `FIT_TAIL_FRACTION` | `exoplasim/scripts/check_relaxation_ceiling.py`, restated twice more in `assess_convergence.py` | one algorithm parameter, stated three times, the restatement declaring itself as one | none |

**The convergence constants and the window they size are circular.** The declared
scatter and memory time produce `DEFAULT_WINDOW_ORBITS`; a run is assessed over
that window; the assessment reports the scatter and memory time IT measured; and
nothing carries that back. The artifact makes the loop visible:

| quantity | declared in `assess_convergence.py` | `run_432e5e46adef_convergence.json` reports |
| --- | --- | --- |
| orbit scatter, K | 0.091 | 0.08333 |
| memory time, orbits | 2.2 | 1.7146 |
| window, orbits | 35, derived from the two above | 35, inherited from that same derivation |

Both differences run in the direction the constants' own comment predicts, since
both were taken as upper bounds on a series that still drifts. That is the point:
a number can be right, be labelled a bound, and still have no way to learn it
moved.

**`TAU_RELAXATION_ORBITS_BRACKET` is the cleanest instance in the audit.** The
producer exists, the consumer exists, the comment names the producer, and
`exoplasim/analysis/convergence/relaxation_ceiling.json` is not on disk. The wire
between them is a human.

**This tier is why the class is worth a name.** The predecessors of the
convergence constants were wrong by about a factor of two, and by a factor of
five in the production span they set, for as long as they stood. Nothing in the
tree objected.

### Tier 2: it reaches the model's namelist or a staged surface field

Each is written into a run's namelist or into a `.sra` the model integrates. An
error is in the physics from the first timestep.

| quantity | where | what computes it | check | stale |
| --- | --- | --- | --- | --- |
| `eddy_wind_m_s` | `config/planet.yaml` | `exoplasim/scripts/measure_eddy_wind.py`, into `exoplasim/analysis/eddy-wind/*.json` | the `timescales_days` table derived FROM it is held to the rule; the wind itself is compared with nothing | no |
| `cold_start_profile.surface_temperature_k` | `config/planet.yaml` | the fitted asymptote of `run_524fbed77a9a`, through `exoplasim/runs/INDEX.json` | none on this copy | no |
| `cold_start_profile.lapse_rate_k_per_m` | `config/planet.yaml` | `lib/lapse.py:dry_adiabat_k_per_km` from the declared composition and gravity, times the declared fraction of neutrality | none | no |
| `cold_start_profile.tropopause_height_m` | `config/planet.yaml` | the same gas properties, and the surface temperature two lines above it | none | no |
| `semi_implicit_reference_temperature_k` | `config/planet.yaml` | its own comment states the rule: the `dsigma`-weighted mean of the `setzt` profile over the model's ten sigma levels | none, and no artifact emits it | no |
| `land_longwave_emissivity` | `config/planet.yaml` | `analysis/rock_emissivity.py` per class, area-weighted over the build's land by `analysis/emissivity_contrast.py` | none | YES |
| `vegetation_albedo`, its bracket and bands, `tree_albedo`, `grass_albedo`, their brackets | `config/planet.yaml` | `analysis/vegetation_albedo.py`, into `analysis/vegetation_albedo.json`, under the same key names | none | no |
| `lithology_albedo_overrides.playa_clastic.albedo` | `config/planet.yaml` | `analysis/playa_albedo.py` | the sibling `replaces` key guards the INPUT; nothing guards the output | no |
| `ozone_height_m`, `ozone_spread_m` | `config/planet.yaml` | upstream's two lengths times `gascon/ga` over Earth's, from `lib/lapse.py:gas_properties` and the declared gravity | none | no |
| `ozone_visible_weight` | `config/planet.yaml` | unresolved; see the ambiguities | none | unknown |
| `snow_fusion_j_kg` | `config/planet.yaml` | `analysis/ice_properties.py` from IAPWS-06, into `analysis/ice_properties.json` | none on the config copy; the script does hold the COMPILED Fortran literals to what it computes | no |
| `surface.land_water_column.layer_capacity_fraction` | `config/planet.yaml` | `layer_thickness_m` on the line below it, over the column base | `run_exoplasim.py` checks length and positivity only | no |
| `CORRK_PATH_CM`, `CORRK_RATIO_TO_EQ21` | `exoplasim/scripts/shortwave_band_weights.py` | the baseline climatology's water path, and `corrk_cross_check.py`'s absorptance at it | none, and the gate above them is circular | no |
| `DEFAULT_CO2_PLANET`, `DEFAULT_CO2_EARTH`, `DEFAULT_WATER_CM` | `exoplasim/scripts/corrk_cross_check.py` | the paths `shortwave_band_weights.py` evaluates at; the comment says to recompute them there rather than trust these | none | no |
| `SIGMA_LOWEST` | `lib/lapse.py`, restated in `hydrography/scripts/carve_verdict.py` | the model's level construction at ten layers; a climatology's `lev` axis carries it | none | no |
| `thermal.saturation_mapping.sr_at_wilting_point`, `sr_at_field_capacity` | `pedology/config/land_column_properties.yaml`, restated in `config/planet.yaml` and in `landmod.f90` | `pedology/scripts/land_column_properties.py` | THE CHECK EXISTS, and is what caught this one | no, and the drift it caught is what this class is named after |
| `thermal.constants`: the two glacier-ice and two snow thermal values | `pedology/config/land_column_properties.yaml` | `analysis/ice_properties.py`, from `rhoglac` and `rhosnow` | the script holds the COMPILED Fortran literals and never opens this yaml, so this is a fourth restatement outside the loop | no |

**`land_longwave_emissivity` was computed on a build the registry refuses.** The
declared value reproduces the mesh land mean of `precarve-craton-10m` exactly.
The active build gives a different one, and nothing said so:

| build | terrain hash | mesh land-mean emissivity |
| --- | --- | --- |
| `precarve-craton-10m`, refused, and where the value was taken | `ab0d679b` | 0.9357 |
| `canonical-10m-base`, the configured build | `20046729` | 0.9350 |

The drift is far inside what the surface energy balance can distinguish: the
audit that sized this key measured the whole per-cell field as worth under half a
watt against a criterion of 1.4. The physics has not moved. What has moved is a
number declared to four decimals, by more than ten times what its own written
precision allows. The blast radius is not this drift but the next one: the key is
DEFINED as an area weighting over the build's land, so every terrain change moves
it. `analysis/emissivity_contrast.json` also still records the declared scalar as
the blackbody value this key replaced, so the artifact predates the key it would
be checked against.

**The shortwave gate is circular with respect to a frozen quantity upstream of
it, and this failure mode belongs in the class definition.**
`check_consistency.py` holds `h2o_sw_level` against
`exoplasim/analysis/h2o_sw_level.json`. That artifact is written by
`shortwave_band_weights.py --level` from `CORRK_RATIO_TO_EQ21` and
`H2O_CONTINUUM_FRACTION` and from nothing else -- no climatology, no spectrum, no
network, which is exactly why the row is regenerable on any tree. So the gate
catches a retyping error into config and cannot catch the climatology moving
under `CORRK_PATH_CM`, the path `CORRK_RATIO_TO_EQ21` was evaluated at. A check
that regenerates its own reference from the frozen value is not a check on that
value. `H2O_CONTINUUM_FRACTION` beside it is correctly out of the class: declared
from published numbers, papers cited, arithmetic in the note.

**`SIGMA_LOWEST` is declared twice, one import apart.**
`hydrography/scripts/carve_verdict.py` imports `reference_height_m` from
`lib/lapse.py` and then passes its own copy of the constant into it. Both agree
with the model's `lev` axis today. Both are properties of the layer count, and
this project compiles one executable per layer count.

### Tier 3: it decides the carve list

Loop A's verdict leaves the project and changes the terrain, so it is the one
output that cannot be revised by re-running something.

| quantity | where | what computes it | check | stale |
| --- | --- | --- | --- | --- |
| `HYDROLOGICAL_RESPONSE_PER_KELVIN` | `scripts/error_budget.py` | a secant between two converged fluxes sharing a surface on this build, over the last ten orbits of each | none; unlike `SLOPE_BRACKET_RUNS` the runs are not even named | unknown, the runs being unnamed |
| `GASCON` | `hydrography/scripts/carve_verdict.py`, restated in `land_water_ledger.py` | the model, from the declared composition; runs record it on their energy artifacts | none | no |
| `WATER_ALBEDO` | `hydrography/scripts/carve_verdict.py` | the export's own water rock class | none | no |
| `EARTH_BAND_LAND_MKM2` | `hydrography/scripts/export_carve_list.py` | nothing in this tree; the companion numerator HAS a re-derivation function beside it and this denominator has none | none | unknown |
| `VESPER_CELL_KM2` | `hydrography/scripts/earth_calibration.py` | planet radius and region count, both carried in `hydrography/data/*/hydrography_report.json` | none | no |
| `score.bar_auc` | `hydrography/config/topographic_index.yaml` | `hydrography/scripts/earth_calibration.py`, as the AUC of the model depth per bore against its mesh region | re-derived and PRINTED beside the declared value; the gate then compares a different statistic | YES |
| `ph.parent_by_category.carbonate`, both `parent_bracket_*` ends | `pedology/config/pedogenesis.yaml` | Slessarev's carbonate-system quartic at this world's `pCO2_bar`; the derivation exists only as prose in `pedology/notes/pedogenesis-value-provenance.md` | none | no |
| `endorheic_alkalinity_bonus_bracket` | `pedology/config/pedogenesis.yaml` | the carbonate pH above, subtracted from a published range | none | follows the row above |
| the land-mean regolith depth that justifies `erosion_coefficient_per_relief_m` | `pedology/config/pedogenesis.yaml` | `pedology/scripts/build_soil.py`, into `pedology/analysis/soil_report.json` | none | YES |
| `clay_conversion_bracket` low end | `pedology/config/pedogenesis.yaml` | one over the maximum weathering intensity over the sixteen type localities, from `pedology/analysis/earth_validation.json` | none | no |
| the theta_fc/theta_s statistics and the cell count they are taken over | `hydrography/config/land_water_ledger.yaml`, restated three times in `biosphere/config/ntransform.yaml` | `pedology/scripts/land_column_properties.py` over the build's soil map | none; `ntransform_gate.py` never opens a soil map, while `mineral_reactivity_gate.py` beside it does | YES |

**The carbonate pH is the sharpest of these because the file says the opposite.**
`pedogenesis.yaml` states that `config/planet.yaml` owns `pCO2_bar` and that
nothing there restates it. The declared pH IS that restatement, one function
application removed, and `build_soil.py` consumes it on every run. No
implementation of the quartic exists anywhere in the tree, so the derivation
cannot be re-run even by hand without transcribing it out of the note.

**`bar_auc` has already drifted and the project noticed.**
`hydrography/notes/subgrid-water-table.md` records the re-derived Australian
figure beside the declared bar. `earth_calibration.py` recomputes that statistic
and prints it next to the bar with the words "the bar was measured here", and
then compares a DIFFERENT statistic against the bar at the gate. This is the good
disposition with the refusal missing, which is the most instructive shape in the
audit: the re-derivation is already paid for.

**The regolith depth is stale in its ARGUMENT rather than in its value.**

| quantity | quoted in `pedogenesis.yaml` as the justification | `pedology/analysis/soil_report.json` `land_means` |
| --- | --- | --- |
| land-mean regolith depth, m | 1.04 | 0.6979 |
| the declared bracket, m | [0.15, 1.30] | unchanged |

The coefficient may well still be right. Its stated reason is not: the mean sits
in the bracket's lower half rather than its upper third, and the branch
attribution that follows no longer holds. A frozen quantity can invalidate an
argument without invalidating the number the argument defends, and that is the
harder case to see.

**`ntransform.yaml` describes the wrong soil map by a factor of four.** Three of
its comments justify a parameterisation choice on the distribution over "the
current soil map's" cells:

| quantity | quoted in `ntransform.yaml` | the map it describes | the configured build's map |
| --- | --- | --- | --- |
| gridcells | 4105 | `precarve-craton` at T42 | 1019, `canonical-10m-base` at T21 |

Every median and count derived from it is on that support. The same figures are
also carried in two dated notes under `biosphere/notes/`, which is correct: those
say when they were measured. The config comments say "current", which is what
puts them in the class.

### Tier 4: it prices, ranks, schedules, or feeds a component nothing yet reads

Nothing here changes a field the baseline run integrates. What it changes is
which work looks worth doing, or what an offline component computes.

| quantity | where | what computes it | check | stale |
| --- | --- | --- | --- | --- |
| `FORCING_ITEMS` dust, sea salt and volcanic sulfate rows | `scripts/error_budget.py` | `exoplasim/scripts/dust_forcing.py`, `aeolian/scripts/build_sea_salt.py`, `aeolian/scripts/build_volcanic_sulfate.py` | none, AND the three artifacts each row names are not on disk | unknowable |
| the roughness distribution row in `OTHER_ITEMS` | `scripts/error_budget.py` | `exoplasim/scripts/build_surface_roughness.py`, into `exoplasim/inputs/t21/roughness_t21_report.json` | none | YES |
| the two lithology land shares inside the albedo items | `scripts/error_budget.py` | the export, regenerated into `world_state.json` as `lithology_land_fractions` | none, and `analysis/rock_albedo_bands.py` states the rule for this exact number in the opposite direction | no |
| `aeolian_z0_by_class_m.playa_clastic` `z0` and `bracket`, `sigma_g` and its bracket | `aeolian/config/dust.yaml` | `aeolian/scripts/playa_roughness_mix.py`, into `aeolian/analysis/playa_roughness_mix.json`, which carries all four | none; the script names the config in comments and never opens it | no |
| `BAND_SHARES`, `CONFIGURED_BUILD` | `aeolian/scripts/playa_roughness_mix.py` | the area-weighted share of `playa_clastic` per landform band on the raw mesh | none | YES: it names a build the registry refuses and has no column for the configured one |
| `stress_reference` | `minerals/config/prospectivity.yaml` | Orogen's `elevation.js`, as a declared quantile of propagated stress; the config carries the identity that would verify it | none | unknown |
| `deep_exhumation_delta` | `minerals/config/prospectivity.yaml` | the p75 of granite exhumation on the export; no script in the tree emits it | none | unknown |
| `potassium_ueq_l` | `minerals/config/downstream_prospectivity.yaml` | the maximum solute concentration over the checked-in table, crossed with the build's rock map | one-sided: `prospectivity_scale.py` raises on an under-estimate and is silent on an over-estimate | no |
| `DECLARED_BOUNDS` | `biosphere/scripts/build_vesper_pfts.py` | the `declareitem` calls in `vendor/lpj-guess/framework/parameters.cpp` | none, though two gates in the same directory parse that source at run time | no |
| the resolution written into the basemap caveat | `maps/build_basemap.py` | the climatology's own `lat`/`lon`, which the same file already opens | none | YES |
| `T42_SECONDS_PER_ORBIT_AT_DT45`, `RUNG_FACTOR`, `OVERHEAD_S` | `exoplasim/scripts/filter_timestep_matrix.py` | this host, by the sweep the file runs, or `bench_ab.py` | `RUNG_FACTOR` entries are replaced in-flight within one invocation; the anchor never is | unknown |
| `SUPPORTS_KM`, and the mesh spacing inside `LAG_KM` | `analysis/escarpment_relief_anchor.py`, `analysis/orogen_resolution.py` | the export manifest's own mean edge, which `analysis/subgrid_slope_support.py` and `analysis/scarp_relief_transport.py` both READ | none | no |
| the three snow-albedo zenith rows | `analysis/snow_albedo_zenith.py` | `analysis/ice_albedo.py`'s `RECORDED` block, which IS reproduced from the vendored source and refuses at 1e-9 | none on the copy | no |
| `DUNNE_SCATTER` | `analysis/spatial_reduction_gap.py` | the antilog of the log-units scatter declared in `pedology/config/pedogenesis.yaml`, which already carries the same antilog | none | no |
| `NLON`, `NLAT`, `NLEV`, `NTRU`, `NSTEPS` | `exoplasim/scripts/bench_pyburn_read.py` | the configured rung, layer count and output interval | none | no |
| `SURFACES` ocean and salt-crust albedos | `exoplasim/scripts/dust_forcing.py` | `exoplasim/config/ocean_tier.yaml` and `analysis/rock_albedo_bands.json`; the other two entries in the same dict ARE read from config | none | no |
| `MEASURED_FLAT_ABOVE_UM` | `exoplasim/scripts/dust_indices.py` | the maximum wavelength of the digitised dataset file the script already reads | none | no |
| `atmosphere_fraction` and `ocean_fraction` | `pedology/config/outgassing.yaml` | the sibling `source_partition` fractions, which ARE re-derived and refused | the partition is checked; the identity linking these two to it is not, and `atmosphere_fraction` has no consumer at all | no |

**The error budget's forcing table cites three artifacts that do not exist.**
`analysis/dust_forcing.json`, `aeolian/analysis/sea_salt_baseline.json` and
`aeolian/analysis/volcanic_sulfate.json` are absent; the generators are present.
The literals are now the only copies of those numbers, so there is nothing left
to compare them against, and the budget that ranks the project's uncertainties
rests on them. `scripts/error_budget.py` is the one half-converted file in the
tree: it reads the albedo report, a run manifest and the climatology for the
water balance, and hardcodes the rest. Its own albedo items were converted to
read from an artifact after one went stale, and the forcing table was not.

**The playa roughness pair graduates to tier 2 the moment dust emission is
enabled.** The class's roughness drives emission as roughly the inverse cube of a
logarithm, `aeolian/config/dust.yaml` is what consumers read, and
`playa_roughness_mix.py` writes the derivation into a JSON nothing opens.
`aeolian/README.md` states the arrangement as the design rather than as a gap.
The mix is also computed with `CONFIGURED_BUILD` naming
`precarve-craton-10m` -- a build `lib/orogen.py` carries a refusal for -- and the
table has no column for `canonical-10m-base` at all.

**`minerals/`'s two frozen thresholds are the audit's clearest case of a
deliberate freeze needing a stated disposition.** `stress_reference` and
`deep_exhumation_delta` are both quantiles of exported fields, both measured on
`precarve-craton`, and both frozen on purpose so that the same number means the
same thing on every build. That intent is defensible. What is missing is that
neither carries a re-derivation gate keyed on the terrain hash, and
`deep_exhumation_delta` names neither a closed form nor a verifying identity, so
nothing could re-take it without a hand measurement.

**`maps/build_basemap.py` writes a wrong resolution INTO an artifact.** The
caveat string is emitted into `basemap_provenance.json` and inherited by every
frame; the same function derives the climatology label from the file it read and
hardcodes the resolution beside it, and the file's own later comments already
reason correctly about the active grid.

### Tier 5: it appears only in prose or in a comment

Nothing reads these. They are in the class because the project's own convention
puts them there: a number in a document earns its place as a decision, a
threshold, an identity, or a magnitude an argument fails without, and a current
measurement is none of those.

| quantity | where | what computes it | stale |
| --- | --- | --- | --- |
| the two stellar-cycle temperature ranges | `config/planet.yaml` `stellar_cycle` notes | the declared amplitudes and damping factors through `lib/sensitivity.py`'s slope | YES |
| the inland-water and internally-draining land shares | `docs/src/pipeline/loops.md` | `hydrography/scripts/build_hydrography.py`, through `world_state.json` | no |
| the endorheic share and the unpreserved pit count | `hydrography/scripts/drainage.py` module docstring | the same producer | YES |
| the erodible and playa land shares | `aeolian/README.md` | `aeolian/scripts/build_dust.py` into `dust_baseline.json`, and `world_state.json`'s `lithology_land_fractions` | YES |
| the Weibull shape figures | `aeolian/README.md` | `build_dust.py:weibull_shape_from_samples`, into `dust_baseline.json` | YES |
| the aeolian roughness bracket low end and the emission factor it is worth | `aeolian/README.md` | `aeolian/config/dust.yaml`, and `aeolian/analysis/dust_intensity_levers.json` | YES |
| the median land runoff | `aeolian/config/dust.yaml`, `build_dust.py`'s emitted caveat, `aeolian/README.md` | the hydrography runoff field; no artifact carries this statistic | unknown |
| the median roughness over dust source cells | `aeolian/README.md` | nothing; the roughness reports carry a land mean, not this median, and `dust.yaml` instructs readers to read the report instead | unknown |
| the transport-convergence consequences | `aeolian/README.md` | the declared cap and floor beside them; the inputs are decisions and these are their outputs | unknown |
| the stripped-arc share | `minerals/README.md`, `minerals/config/prospectivity.yaml` | the lithology census over the export | YES, by inference: the granodiorite share it rests on has moved |
| the mesh mean edge and region count | `minerals/README.md`, `prospectivity.yaml`, `build_prospectivity.py`, `maps/README.md` twice, `maps/build_basemap.py` | the export manifest, through `lib/orogen.py`'s registry | YES |
| the projection measurement tables | `maps/README.md` | `maps/render_projections.py`, which emits every one of them per frame into a provenance JSON | unknown; the only frame index on disk is a build old |
| the lapse rate quoted for the orographic correction | `maps/README.md` | `lib/lapse.py:environmental_lapse_k_per_km`, which the code calls and whose comment names the number the README states as the one it stopped using | YES |
| the land-mean bucket depth | `docs/src/reference/config-rationale.md` | `exoplasim/scripts/build_surface_soil_water.py` | unknown; no artifact carries it |
| the band-1 flux share, twice | `docs/src/reference/config-rationale.md` | `lib/stellar.py`, which reproduces `radmod.f90:solarini` | no |
| the mesh mean edge, five more restatements | `docs/src/pipeline/components.md`, `config-rationale.md`, `sequencing.md`, `economic-minerals.md` twice | `analysis/orogen_resolution.py`, from the region count and the radius | no |
| the per-rung refusal ceilings, twice | `docs/src/pipeline/sequencing.md` | `exoplasim/analysis/stability_probe.json`, which `lib/rungs.py` DOES re-derive its own copy from | no |
| the land fraction | `docs/src/reference/config-rationale.md` | the export, through `world_state.json`, which already carries a second gridded value differing in the third digit | no |
| the grid-cell width used to argue about mountain peaks | `docs/src/pipeline/state.md` | the rung and the planet radius | YES, for the declared operating support |
| the land-column capacity split distribution, three restatements | `vendor/exoplasim/.../landmod.f90` twice, `exoplasim/scripts/build_surface_soil_water.py` | `pedology/scripts/land_column_properties.py` | no |
| the build's land-cell count at T21 | `vendor/exoplasim/.../landmod.f90`, `lib/builds.py` | the T21 land mask, materialised as the row count of the soil map `lib/builds.py` itself resolves the path to | no |
| the re-evaporation magnitudes, and the timestep they were taken at | `vendor/exoplasim/.../rainmod.f90` | the form the same file evaluates, which is linear in the model timestep | YES |
| the two shortwave aerosol ratios | `vendor/exoplasim/.../aeromod.f90` | the dust effective radii, through `aeolian/config/dust.yaml` and `analysis/dust_optics.json` | unknown |
| the K-dwarf band-gap shares | `vendor/exoplasim/.../radmod.f90` | the `k25v` spectrum against the albedo grid; the only one of its group citing no producer | unknown |
| the arc land share and the convergent-boundary ratio | `vendor/orogen/js/terrain-config.js` | the generator itself; `world_state.json` carries `lithology_land_fractions` | YES |
| the mesh cell area arithmetic around the basin floor | `vendor/orogen/js/terrain-config.js` | the region count and the radius, which `analysis/orogen_resolution.py` recomputes independently | no |
| the thermal-inertia arithmetic in prose, and the bulk-density p5 | `pedology/config/land_column_properties.yaml` | `land_column_properties.py`; the endpoints they descend from ARE re-derived and refused, so most are protected transitively; the p5 is not | no |
| the roughness span, transfer-coefficient ratio and runoff shares in the carve verdict's docstring | `hydrography/scripts/carve_verdict.py` | the roughness field, and the run's own climatology | unknown |
| the land means quoted as findings | `pedology/README.md` | `pedology/analysis/soil_report.json` | YES |

**The stellar-cycle ranges are the most wrong thing found, and the arithmetic
dates exactly when they froze.** Each is amplitude times damping times the
flux-to-kelvin slope:

| component | declared amplitude | declared damping | at the superseded slope 150.2 | at the live slope 202.0 | the note says |
| --- | --- | --- | --- | --- | --- |
| medium | 0.025 | 0.83 | 3.12 K | 4.19 K | 3.1 K |
| long | 0.035 | 0.99 | 5.20 K | 7.00 K | 5.2 K |

`lib/sensitivity.py` records the superseded slope, why it was wrong -- its
measurement bracket did not contain the baseline flux -- and that it is 26 per
cent low. `sensitivity.verify()` re-derives the live slope from the run index on
every `check_consistency.py` invocation and currently passes. These two sentences
learned nothing when the slope moved, and they are what a reader quotes for what
this world's stellar cycle is worth over a life.

**`loops.md` carries two frozen measurements about two hundred lines above its own
statement of the rule.** Both are rounded `world_state.json` values. The paragraph
they sit in argues about where the thin soil is, and it survives without either.
They also demonstrate the general hazard: `world_state.json` already holds two
endorheic shares differing in the fourth digit, so a reader copying "the" value
has to pick one.

**`drainage.py`'s docstring describes no build in the tree.**

| quantity | the module docstring | `hydrography/data/canonical-10m-base/hydrography_report.json` |
| --- | --- | --- |
| endorheic share of land | 63 per cent | 0.7444 |
| unpreserved pits | 220,649 | 583,684 |

The other two builds read 0.7622 and 0.8103, so the docstring's figure has never
described any build the tree holds.

**`rainmod.f90` states its magnitudes at a timestep the project no longer runs.**
The comment computes them at twice a thirty-minute step and names that as this
project's step; `config/planet.yaml` declares forty-five minutes and the
escalation route runs it at every rung. The derived quantity is linear in the
timestep, so all three magnitudes are two thirds of what the model integrates.
The conclusion the comment draws is unaffected, which is why it went unnoticed.

**`aeolian/README.md` is the densest cluster in tier 5, and it argues against
itself.** It states, ten lines above the first stale figure, that a table of the
same numbers there would go stale the next time the roughness mosaic or the wind
sample moved, and that it had twice. Four of its figures have since drifted. The
component also contains the best instance of the opposite discipline:
`build_dust_source_fields.py:measured_weibull_shape` reads the shape out of the
artifact and RAISES rather than falling back.

**The roughness span in `fluxmod.f90` is the repaired shape, and worth reading as
the model.** The span is still written out, but it names the build it was measured
on and says the report is where the field is read from. Both the comment and
`notes/audits/tuned-values.md` cite `exoplasim/inputs/t21/roughness_t21_report.json`,
which carries `land_min_m` and `land_max_m` and currently reproduces the quoted
span. The remaining weakness is that citing an artifact is not reading one:
nothing compares the two, and the repair rests on a reader following the citation.

## Where a frozen quantity could not be told from a decision

Eight cases where the boundary genuinely does not resolve from the text. Each is
a finding in its own right, because an unstated disposition is what lets a value
be read either way by whoever needs it to be.

1. **`minerals/config/prospectivity.yaml`'s `stress_reference` and
   `deep_exhumation_delta`.** Both are quantiles of exported fields, so derived.
   Both are frozen ON PURPOSE, so that a score means the same thing on every
   build, which is a decision. The file argues both sides. The freezing is the
   decision; the value is still derived, and nothing would notice if the plate
   model moved under it.
2. **`config/planet.yaml`'s `luminosity_solar`.** Seventeen significant figures
   is a computed number, and no producer exists in the tree.
   `docs/src/reference/config-rationale.md` reasons that surface gravity need not
   be declared BECAUSE mass, luminosity and effective temperature already fix it,
   which treats this as an input. If it was itself derived off-tree, it is
   frozen.
3. **`config/planet.yaml`'s `ozone_visible_weight`.** Its sibling
   `ozone_uv_weight` is cleanly paper-sourced with the arithmetic shown in
   `exoplasim/notes/ozone.md`. The visible weight has no equivalent derivation in
   that note or in any script, which does not rule out its being computed from
   the `k25v` spectrum that `build_stellar_spectrum.py` rewrites in place.
4. **`hydrography/scripts/earth_calibration.py`'s `ET_MAX_MM_YR`.**
   `hydrography/config/groundwater.yaml` declares the corresponding Vesper-side
   value as a COMPUTATION, which `build_groundwater.py` honours by reading an
   evaporation grid. On the Earth harness there is no field to read and a flat
   rate stands in. That reads as a decision that impersonates a derived quantity,
   and it is separately duplicated as a bare literal at two of its three uses.
5. **`biosphere/scripts/build_lpj_driver.py`'s `SOIL_CODE_BY_ROCK`.** Documented
   as the judgement fallback for when no soil map is given, which makes it a
   decision. But `pedology/config/pedogenesis.yaml` carries a per-lithology parent
   texture table and `build_soil.py` computes an evolved texture per cell, so the
   same physical claim is encoded twice in two currencies with nothing tying
   them. Duplicated judgement rather than frozen derivation, and it fails the
   same way.
6. **`pedology/config/land_column_properties.yaml`'s `column_base_m`.** Declared
   three times across two components, and the pedology comment asserts that it
   matches the ledger's copy. Nothing verifies that assertion. A duplicated
   decision, not a frozen derivation, with the same failure mode.
7. **`vendor/orogen/js/basins.js`'s comparison land share.** The sentence says "on
   a real planet", which is ambiguous between Earth (a paper-sourced comparison,
   out of the class) and this build's finished terrain (a measurement, in it).
   Nothing else in the file disambiguates it.
8. **`docs/src/pipeline/state.md`'s two glacier magnitudes.** Both cite
   `notes/glacier-rough-pass.md`, which is the mitigated form the roughness-span
   case was repaired into, and both are magnitudes the paragraph's argument needs.
   They read as exempt. The grid-cell width in the same section does not: it is
   derived geometry restated in prose, and it is stated for a rung that is not
   the declared operating support.

A ninth, recorded separately because it is a different question: `maps/` and
`docs/src/reference/economic-minerals.md` both carry dated A/B measurement tables
that read correctly as historical records but live outside `notes/`, where the
convention puts dated records. That is a filing question, not a freezing one, and
it is settled by whoever owns the convention rather than by this audit.

## What is true, in one paragraph

The tree's discipline is real and uneven. It contains fifteen worked instances of
a declaration held to its own derivation by a check that fires, several of them
written specifically to close this class of defect, and three configuration
blocks that refuse to state a number at all. It also contains sixty-four places
where a computed number was written down with no path back, eleven of which have
already drifted, one of them by a third in reader-facing prose about the planet's
climate. Every frozen instance sits within a few lines of code of the machinery
that would close it, and in five cases within the same file. The one drift that
was caught and repaired rather than integrated is the one whose declaration had a
check on it, which is the whole of the argument for the second disposition.
