# Aeolian

Mineral dust: emission, transport and deposition, computed offline from a
climatology rather than inside the climate model. `notes/dust.md` carries the
decision and why it is not in the GCM; this component is the implementation.

```bash
python aeolian/scripts/build_dust.py                          # the baseline source map
python aeolian/scripts/build_dust.py --variant arid_bare_ground   # the bracket on vegetation
python aeolian/scripts/build_dust_source_fields.py            # the same map as .sra boundary fields
python aeolian/scripts/build_dust_source_fields.py --self-test  # the identity, and that it can fail
```

Products are `analysis/dust_<variant>.json` and a matching `.nc` carrying annual
mean optical depth, deposition flux, emission and the erodible fraction.

`build_dust_source_fields.py` is the bridge to the in-model scheme: it writes the
same per-cell source map as ExoPlaSim surface codes 1801, 1802 and 1803, plus the
`aero_nl` group that goes with them. It is a separate script rather than a flag
on `build_dust.py` because it produces a boundary condition rather than an
answer, and because it must be regenerated whenever the terrain, the lake
solution or the soil moves, which is a different cadence from the offline run.

## Wind input health

Runs predating the low-I/O patch carry a corrupted first output bin in `spd`,
`ua` and `va`, inflated up to 7.5x at the bottom level (failure-modes class
14; `exoplasim/notes/first-output-bin.md`). Emission goes as roughly u* cubed
above a threshold, so one bad bin can be the whole annual total.
`flag_anomalous_bins` detects it and the run excludes it loudly rather than
consuming it.

**Transport now converges, and the reason it did not is worth recording.** It
was neither a CFL violation nor a cycle nor a missing sink. The explicit step is
CFL-limited by the convergence of the longitude grid at the poles, so the polar
rows set the step for the whole globe, while the finest bin relaxes on a
timescale of years where nothing rains: the integration was being stopped a few
per cent of the way and reporting not-converged, correctly. A declared polar
`cos(lat)` floor relaxes the step and the iteration cap was raised; all cases
now report converged. The floor, the cap and the arithmetic behind both are in
`config/dust.yaml` under `transport`, and the cell width and step they imply
move with the rung, so they are not restated here.

### What it says now, with the wind tail measured and gravity in the threshold

`aeolian/analysis/dust_baseline.json` carries emission, land-mean and global
optical depth, deposition and the roughness each arm actually ran at, for all
three ends of the bracket at once, and it records the Weibull shape it used with
the sample count and the file it was fitted from. Read that; a table of the same
numbers here would go stale the next time the roughness mosaic or the wind
sample does, and it has twice.

What the artifact says, in shape rather than in figures: Earth emits about 2000
Tg per Earth year at a land-mean dust optical depth near 0.03, and every arm of
this world's bracket except the roughest is above both, the central one by
several times on emission and by more than ten on optical depth. This is a dusty
world.

**Both of the large corrections went in opposite directions and neither cancelled
the other.** Measuring the wind tail raised emission by more than an order of
magnitude, because a Weibull fitted to snapshots days apart is too narrow and
biases the shape high: snapshots resolve synoptic variance and not sub-daily
variance, so a shape fitted to them is an upper bound and the emission a lower
one. The shape actually used, the sample count and the file it was fitted from
are all in `dust_baseline.json` under `subgrid_wind`, which is where
`build_dust_source_fields.py:measured_weibull_shape` reads it and refuses if it
is absent. Putting this world's gravity into the saltation threshold then cut
emission back, a few per cent on the threshold amplified by the same u*
nonlinearity working the other way. The threshold correction is the fourth root of the gravity ratio and
not the square root; `aeolian/config/dust.yaml` carries the derivation.

The roughness bracket is the only large uncertainty left and it no longer spans
the answer: the rough end returns a real emission rather than the zero it once
read, and the ends are now the per-lithology mosaic rather than one scalar.

### The reopening test, stated without tuning

`notes/dust.md` reopens the in-model question at a land-mean optical depth above
0.10, and the baseline artifact's `reopening_test` block carries that threshold,
the range it is held against and the verdict. Only the roughest end of the
bracket is below it, and that is the treatment that tabulates the roughness by
lithology: resolving the roughness WITHIN a lithology raises the emission at
every arm, and at every measured within-class spread but the narrowest it lifts
the roughest arm across the threshold too. The crossing is therefore a floor
rather than an estimate. `notes/audits/dust-intensity-levers.md` carries the
measurement.

That was a pre-committed threshold, fixed before the answer was known, and what
it commits to is that a prescribed field is no longer defensible and the emission
scheme belongs inside the model. That is DUST-3's problem now, not a conclusion
this file should quietly draw.


## What the component does get right

**The source map, which is the thing the built-in GCM scheme gets wrong.**
`fcoeff * land_mask` emits as much from forest as from salt pan. This one starts
from `substrate_class`, which is consolidated lithology plus closed-basin fill,
so the only unconsolidated material is the fill; weights the two barren classes
separately, because a cemented salt crust is not a silicate soil and Kok's
fragmentation theory does not describe halite cement; removes standing water
from the solved lake extent; and removes snow. The bare erodible share of land
it gives is well below the `playa_clastic` share of land it starts from, and the
difference is lakes and crust rather than an assumption. Both are per build:
`dust_baseline.json` carries `source_map.erodible_fraction_of_land` and
`world_state.json` carries `lithology_land_fractions`.

**The physics is grounded rather than recited.** Every constant traces to a
fetched primary source, listed in `config/dust.yaml` with what it is worth:
Kok et al. (2014) for emission, Kok (2011) for the emitted size distribution,
Fecan et al. (1999) for the soil-moisture threshold, Marticorena and Bergametti
(1995) for the drag partition, Sportisse (2007) for below-cloud scavenging.

**Settling uses this world's gravity.** At 12.81 m/s2 a given particle falls
1.31x faster than terrestrial intuition, so dust lifetime here is shorter than
Earth analogues suggest. That is in `settling_velocity` rather than in a comment.

## What is declared rather than derived, and where the brackets are

`config/dust.yaml` carries all of it. The widest are, in order:

1. **The subgrid wind shape**, now measured but only to an upper bound. There is
   no default gust source and there deliberately is not one: `build_dust.py`
   refuses to run unless it is given `--gust-samples` or told in as many words to
   fit from the snapshot climatology instead. The two answers differ by a factor
   of 40, and an optional argument whose absence means do the wrong thing is
   `docs/src/practice/failure-modes.md` class 2. DUST-15.
2. **The aeolian roughness of the erodible surface**, bracketed per lithology
   in `config/dust.yaml` under `aeolian_z0_by_class_m` and worth the emission
   factor `dust_intensity_levers.json` reports as
   `pattern_vs_intensity.emission_bracket_factor`. Neither the ends nor the
   factor is restated here: the ends move whenever the mixture behind a class is
   re-derived, and the factor is measured against whichever ends are declared.
   The grid-cell roughness field is deliberately NOT used: its median over
   source cells is orders of magnitude above a roughness element, because it is
   a kilometric orographic variance, and feeding that to a scheme built for
   centimetre-scale elements returns zero emission everywhere. The roughness
   reports under `exoplasim/inputs/` carry that field's own statistics. Sheltering of a patch by
   the terrain around it is therefore not represented, which biases emission up.
   **The bracket belongs to ONE class.** Holding `evaporite` at its centre and
   sweeping `playa_clastic` across its own range reproduces almost the whole
   bracket, and doing the reverse moves the emission by under a percent, which
   follows from the erodible weights. Narrowing the intensity therefore means
   narrowing `playa_clastic`'s roughness and nothing else, and neither end of
   that class's bracket is a measurement of a clastic playa: the smooth end is a
   modelling convention for an active dust source and the rough one is the
   sand-desert boundary stepped down an order.
   **And the roughness within a class is not one value.** That spread is
   measured, in `config/dust.yaml` under `within_class_z0`, from the repeat
   entries of Prigent et al. (2005) Table 1. Collapsing it to the tabulated
   point, which is what the arms do, biases the emission LOW at every arm, so
   the intensity this component reports is a floor on that axis. It is reported
   as a correction rather than folded into the arms because each class's bracket
   is the same kind of range and counting it twice would be double counting.
   world-03x; `notes/audits/dust-intensity-levers.md`.
3. **The evaporite erodible weight**, 0.1 with a bracket of 0.0 to 0.3, for
   crust cementation. Declared suppression, not measured efficiency.
4. **Vegetation cover**, which is not modelled at all. Non-barren land is
   assumed to carry a canopy and not emit, following the project's existing
   declared position, which is generous on a world whose median land runoff is
   a few mm per Earth year -- `config/dust.yaml` carries that figure with the
   build it was measured on, and no artifact in the tree re-derives it. `--variant arid_bare_ground` brackets it, and that
   variant turns out to be worth under a fifth on emission at every roughness
   arm, one-signed. That is why the in-model arm's inability to express it is a
   declared gap and not a defect; world-4qem.
5. **Which refractive indices this world's dust has**, `optics.indices`, chosen
   as the measured datasets in the shortwave and OPAC in the thermal infrared,
   with OPAC kept as the absorbing bracket. This is the one entry in the file
   that something OUTSIDE this component also reads: the aerofile ExoPlaSim's
   radiation loads and the offline forcing are both built from it, through
   `exoplasim/scripts/dust_indices.py`, so that the model and the pricing
   describe one particle. The config carries the argument; DUST-12.

## The two size distributions, and which one to read

This component carries two dust size distributions. They are not alternatives
and neither is derivable from the other here; they describe different
populations at different points in a particle's life, and both are right where
they stand.

| | declared at | what it is | read it for |
| --- | --- | --- | --- |
| **emitted** | `aeolian/config/dust.yaml:size_distribution` | Kok (2011) eq. 6: median **diameter** 3.4 um, geometric sigma 3.0 | emission, the bin edges, the coarse mode near a source |
| **optics** | `aeolian/config/dust.yaml:optics_size_distribution` | Balkanski et al. (2007): number median **radius** 0.295 um, sigma_g 2.0, 2.6 g/cm3 | optical depth, single-scattering albedo, asymmetry, the transported accumulation mode |

**Check the units before comparing them.** One is a diameter and the other a
radius, so the naive comparison is wrong by about a factor of six before sigma
is reached: 3.4 um diameter is 1.7 um radius against 0.295. The remaining gap
is real and is the point -- the coarse mode has largely settled out by the time
the optics matter, so the population the model's optical properties describe is
not the population that left the ground.

The optics distribution is DECLARED in `dust.yaml` and RESTATED as module
constants in `exoplasim/scripts/dust_optics.py`, which is what integrates the
Mie calculation over it. It is restated rather than read because reading would
make `dust.yaml` an input of `analysis/dust_optics.json` and force a Mie
regeneration for no change in any number; `scripts/smoke_test.py` asserts the
two agree instead. `analysis/dust_optics.json` echoes the distribution it used,
so an artifact always says which population its numbers describe.

## What it does not do

No dust-climate feedback: the climatology is an input and does not respond.
No vertical structure: a well-mixed column of declared scale height advected by
a single steering wind. No inter-bin microphysics, which is correct for mineral
dust because it neither coagulates nor grows appreciably.

**No nocturnal-inversion low-level jet, and that one is a floor rather than a
bracket.** The Bodele mechanism -- a sharp nocturnal inversion decoupling the
surface, the jet above it mixing down after sunrise -- supplies a large share of
Earth's dust from one dry lake bed, and `playa_clastic` is declared as that
class's analogue here. The model has two levels below a kilometre and one below
half a kilometre, so the inversion and the jet core both fall inside its single
lowest layer; the wind that reaches this component is that layer's mean, and
sampling it more often cannot recover a structure the model has no levels to
form. Basin emission is therefore understated one-sidedly, at EVERY end of the
roughness and vegetation brackets rather than at one of them, which is why it is
a floor. Declared rather than parameterised: a jet parameterisation fitted to
nothing this world can measure would be precision theatre. DUST-16;
`notes/external-review-abl-obliquity-dust-carbon.md` point 1 carries the
measurement and what would revisit it.

## One-off tools

- `scripts/dust_runoff_sensitivity.py` -- one-off: converts a precipitation change into a basin count, which is what prices the catchment half of the dust question in `docs/src/pipeline/sequencing.md` A4. Registered under `one_offs` in `config/pipeline.yaml`.
- `scripts/extract_high_cadence_wind.py` -- one-off: pulls instantaneous near-surface winds out of a high-cadence run segment. This is what DUST-5 measured and what `build_dust.py --gust-samples` must be given; the 32-sample snapshot climatology fits a Weibull shape about twice too high and costs a factor of 40 on emission. Registered under `one_offs` in `config/pipeline.yaml`.
- `scripts/dust_intensity_levers.py` -- one-off: which lever moves the intensity, as RATIOS of emission totals against `dust_baseline.json`. It prices the within-class roughness collapse, the between-class collapse, whether the roughness bracket moves the pattern or the total, what the soil texture is worth, what the scalar in-model roughness costs, and what the vegetation bracket costs. Transport is not re-run and the artifact measures the transfer that licenses that. Registered under `one_offs` in `config/pipeline.yaml`; `notes/audits/dust-intensity-levers.md` reads it.
- `scripts/playa_roughness_mix.py` -- one-off: derives `playa_clastic`'s aerodynamic roughness and its bracket from Greeley et al. (1997)'s wind-profile measurements, mixed geometrically over the landform bands `notes/audits/orogen-lithology.md` measures the class to be made of. The class is a residual of two geometric conditions rather than a clay plain, so a single measured playa is the wrong endmember for it. What the script writes is a derivation; the numbers consumers read are the ones carried in `config/dust.yaml`. Registered under `one_offs` in `config/pipeline.yaml`.
- `scripts/roughness_mixing_order.py` -- one-off: whether the ORDER of that mixture changes the answer. MB95's drag partition is affine in `ln z0`, so the geometric mean reproduces the area-mean drag efficiency exactly and the first arm is that identity rather than a comparison; the friction velocity the same value divides is not affine in `ln z0`, and the second arm sizes what the order costs there against the level bracket each class declares in `config/dust.yaml`. SPAT-7; `notes/audits/nonlinear-spatial-reductions.md` section 8. Registered under `one_offs` in `config/pipeline.yaml`.

## The in-model port, and why this component survives it

`aeolian/notes/in-model-dust.md` is the design for putting emission, deposition
and scavenging inside ExoPlaSim (DUST-3), with each piece's predicted effect and
what result would falsify it. Its patches are RESIDENT in the `vendor/exoplasim`
subtree; the files under `exoplasim/patches/` are the record of what each
changed and why, and each names its base shas, which is where the stacking order
is kept. Three of them carry the emission chain:

- `exoplasim/patches/exoplasim-3.4.2-aerocore-defects.patch`, seven latent
  defects in `aerocore.f90` and `aeromod.f90`, unconditional.
- `exoplasim/patches/exoplasim-3.4.2-aerosol-deposition.patch`, a dry deposition
  velocity behind `ldepvel` and Sportisse below-cloud scavenging behind
  `lwetdep`, both defaulting to off and both reading their coefficients from
  `aeolian/config/dust.yaml` rather than carrying Fortran defaults.
- `exoplasim/patches/exoplasim-3.4.2-dust-emission.patch`, the three boundary
  fields, the gathers, the `aero_nl` calibration group and Kok (2014) equation
  18 in the source term, behind `ldustemit` and defaulting to off.

**Three fields and no fourth, and that is now a priced decision rather than a
gap.** The roughness enters the in-model scheme twice: through the drag
partition, which field 1802 carries per cell, and through `ln(zref/z0)` in the
friction velocity, which reads the scalar namelist `DUSTZ0` because the model
computes `zref` from each cell's own temperature and the ratio is not a static
prefactor. `DUSTZ0` is the roughness mosaic's own erodible-area weighted
geometric mean at the arm being written, computed from the fields as they are
built rather than looked up in the config; with it the in-model emission sits
within a fraction of a percent of the offline arm at every arm, so a fourth
boundary field carrying the roughness per cell would buy nothing. It used to
carry the class-blind scalar, whose bracket ends predate the per-lithology
values and sit outside them, and that put the two arms as much as four times
apart. world-h24h; `notes/audits/dust-intensity-levers.md`.

The run side is `model.dust_emission` in `config/planet.yaml`, which is absent
and therefore `none` by default. Setting it makes `run_exoplasim.py` stage the
three fields, write the whole `aero_nl` group from the provenance file beside
them, and set `L_AERO = 1` and `l_source = 2`. `model.dust_dry_deposition` and
`model.dust_wet_scavenging` switch on the two removal terms independently, which
is what lets each be its own A/B arm off one binary.

The emitted dust is RADIATIVELY INERT by default and
`model.dust_emission_radiative` is the switch; the driver writes `l_aerorad`
from it. The two reasons that setting was chosen for have both since closed:
`aero_ini` populates `radmod`'s `apart` from the one `aeromod` declares, and the
transported aerosol has its thermal-IR absorption through `aeroqlw`, which
`radini` aborts without and which `build_surface_dust.py` now derives beside
`DUSTQLW` so something writes it. So turning the radiation on is a decision
about which aerosols the climate carries rather than a wait on missing physics,
and it is taken in DUST-13. The exclusion that framing used to carry is gone
with `CLIM-39`: the transported tracer is one species of the array at index
`ndustrad+1`, so a prescribed field and a transported one no longer share a
slot and both can be on at once.

`aeolian/notes/multi-species-aerosol.md` is the companion design, for the other
end of the same interface: the radiation carries ONE aerosol with one global set
of optical properties, so dust, sea salt and volcanic sulfate cannot be in it at
once. It sizes what carrying N species takes and why the two-species form is not
the cheap version. `CLIM-39` and `CLIM-40`.

**This component does not retire when that lands.** ExoPlaSim's aerosol is one
tracer with one radius and one density fixed at compile time, so the in-model
chain gets the burden and the optical depth and cannot provide a
composition- and size-resolved soil nutrient input. The offline artifact emits
total mineral deposition, but `phosphorus_budget.py` does not currently read it
and it contains no elemental composition or bioavailable fraction. ANUT-3 owns
that missing soil-facing transformation. This is the declared cost of DUST-8
and the note carries the radiation-side numbers.

## Sea salt

The second aerosol, added 2026-08-19 for CLIM-27. `notes/audits/unpriced-terms.md`
finding 2 says why: this project priced mineral dust carefully and had never
asked what any other aerosol was worth, on a world whose ocean is more than half
the surface.

```bash
python aeolian/scripts/sea_salt_source.py     # the Earth check, and nothing else
python aeolian/scripts/sea_salt_optics.py     # step sea_salt_optics
python aeolian/scripts/build_sea_salt.py      # step sea_salt
python aeolian/scripts/build_sea_salt.py --bracket   # every declared end
```

`sea_salt_source.py` is the source function, Grythe et al. (2014) equation 7,
shared by the other two because one needs its mass and the other its shape.
`sea_salt_optics.py` writes `analysis/sea_salt_optics.json`, band-averaged
optics per dry size bin and per relative humidity. `build_sea_salt.py` writes
`aeolian/analysis/sea_salt_baseline.{nc,json}`: emission, burden, optical depth
and the top-of-atmosphere shortwave forcing, at both ends of every declared
bracket. Constants are in `aeolian/config/sea_salt.yaml`.

**Three things about it are not obvious from the dust component beside it.**

The particle is WET and the mass budget is DRY. A sea-salt particle at ambient
humidity is a solution droplet, twice the dry diameter at 80% relative humidity
and a quarter of the dry mass fraction, so its extinction is computed on the wet
particle and divided by the dry mass. The growth curve is OPAC's, in
`exoplasim/data/sea_salt/`, and the optics script checks that OPAC's growth
factor and its wet density say the same thing.

Wet removal is NOT dust's. Sea salt is the most cloud-condensation-active
aerosol there is and is removed in cloud rather than below it; using dust's
below-cloud coefficient gives the accumulation mode a three-week lifetime
against a measured day. The level is anchored to Jaegle et al. (2011).

The spume mode is reported apart. Grythe's third lognormal is centred at 30 um
and its lower tail dominates the emitted mass below the 10 um cut, so the
product carries an `all_modes` and a `no_spume` variant and any comparison with
an Earth number is against the second. Neither is a correction to the other.

## Deposition mass, and why it is its own artifact

`aeolian/scripts/aerosol_deposition.py` is the elemental split and the writer
shared by sea salt and volcanic sulfate, and
`--self-test` is its identities: that each composition table closes on one and
that a mistranscribed one is rejected, that the element split conserves mass and
that a monatomic ion carries its own solute mass exactly, that the writer's
land mean reproduces the analytic mean of a uniform field, and that the mass is
EXTENSIVE: the written file's area integral reproduces the input's, the bin and
element fields partition the total exactly, and the integral is invariant under
an area-weighted change of support while the unweighted aggregation is shown to
move it, so the check is seen to be able to fail.

```bash
python aeolian/scripts/aerosol_deposition.py --self-test
```

`build_sea_salt.py` writes `aeolian/analysis/sea_salt_deposition.{nc,json}` and
`build_volcanic_sulfate.py` writes
`aeolian/analysis/volcanic_sulfate_deposition.{nc,json}`, both declared against
their own steps in `config/pipeline.yaml`. Wet and dry removal per size bin, in
kg per m2 per second, plus a field per element; the JSON reads the same fields as
land and ocean means in mg per m2 per Earth year, because the magnitude anchor is
an Earth measurement.

**They are separate files rather than variables in the baseline artifacts, and
that is not tidiness.** `biosphere/config/abiotic_nutrients.yaml` lists
`aeolian/analysis/volcanic_sulfate.nc` and
`aeolian/analysis/sea_salt_optics.json` in `screen.forbidden_carriers` and
`abiotic_nutrient_ledger.py:check_screen` rejects any candidate whose declared
carrier contains one of those strings. The refusal is by path and it is right:
those files hold optical depth, an optical depth folds in refractive index, size
distribution and water uptake, and no elemental mass comes back out of it. A
carrier file therefore has to hold mass and nothing else.

**The deposition is a PARTITION, not a second calculation.** At steady state each
size bin loses `loss * m` per unit area per unit time and `loss` is already the
sum of a settling frequency and a wet scavenging rate, so the two fields are that
one rate split in two. `build_sea_salt.py` and `build_volcanic_sulfate.py` each
assert that they sum back to the removal their own mass balance computed, to
within rounding, which catches a sign error, a dropped bin, or a scale height
edited in one of its two spellings -- none of which the mass residual can see.

**The composition is Earth's Reference Composition, declared rather than
derived.** `aeolian/config/sea_salt.yaml` carries Table D.3 of the TEOS-10 Manual
in full with its argument; nothing in this project models ocean chemistry, so the
proportions are an assumption, and it is the same assumption the source function
already makes by carrying no salinity weight. Phosphorus is deliberately absent:
Chadwick et al. state the seawater contribution of P to marine aerosol is
vanishingly small, and booking one would invent it.

## Volcanic sulfate

The third aerosol, added 2026-08-19 for CLIM-28.

```bash
python aeolian/scripts/build_volcanic_sulfate.py             # step volcanic_sulfate
python aeolian/scripts/build_volcanic_sulfate.py --bracket   # every declared end
```

Writes `aeolian/analysis/volcanic_sulfate.{nc,json}`. Constants are in
`aeolian/config/volcanic_sulfate.yaml`, the OPAC sulfate tables in
`exoplasim/data/sulfate/`.

**The source is a coupling, not a scaling, and that is the reason it exists.**
Dust and sea salt are both a wind acting on a vast surface. Volcanic sulfur is
set by how much magma this world degasses, and `pedology/analysis/weathering_fluxes.json`
already carries that: at steady state the atmosphere's CO2 is set by outgassing
balancing silicate weathering, so `implied_outgassing_over_earth` is a statement
about magma. Carn et al.'s satellite-measured passive volcanic SO2 flux, scaled
by it and distributed over the arc classes the fork places about the volcanic
front, is the emission. The step therefore NEEDS `weathering_fluxes` and refuses
to run without it rather than falling back to a literal.

**It comes out one to two orders of magnitude below the other two**, and the
answer is worth having for what it bounds rather than for what it adds. The
sulfate pathway on this world is weak: a four-day aerosol from a source far
below the wind-driven ones. Reaching sea salt's forcing would need a sulfur flux
more than ten times Earth's, which is the bound that also covers the two sulfur
sources this project cannot compute -- explosive eruptions, which need an
eruption history `docs/src/reference/no-time-axis.md` says does not exist, and marine
biogenic sulfur, which needs a marine biosphere this project does not have at
all.

**That bound does not reach carbonaceous aerosol**, but carbonaceous is not one
optical species.  Primary fire black/brown carbon can absorb strongly and stays
under CLIM-29.  Secondary organic aerosol is usually represented as mostly
scattering, with pathway- and wavelength-dependent brown-carbon absorption, and
its sign cannot be declared before chemistry, aging and optics exist.  The LPJ
BVOC module is present but off; it supplies only a source inventory, not the
oxidation, partitioning, aerosol burden, optics or cloud response.  Those stages
are audited in `biosphere/notes/bvoc-soa-atmospheric-coupling-audit.md` and
tracked as BVOC-1 through BVOC-10.
