# Vegetation demography, disturbance and dispersal audit

This note records the first of the follow-up biosphere investigations: the
active LPJ-GUESS-CNP cohort-demography path, excluding fire, hydraulics and soil
biogeochemistry except where the demography code directly changes their pools.
It distinguishes a defensible potential-natural-vegetation experiment from a
claim about transient biogeography.  No LPJ-GUESS build or simulation was run.

## Scope and sources read

The source audit covered `modules/vegdynam.cpp`, the relevant framework state
and common outputs, `data/ins/global.ins`, and Vesper's PFT generator and run
harness.  The primary sources below were fetched into `references/`, read, and
recorded in `references/INDEX.md`:

- Fulton (1991), juvenile growth and adult recruitment;
- Smith et al. (2001), the original LPJ-GUESS vegetation-dynamics formulation;
- Sitch et al. (2003), the LPJ population-mode formulation and global
  evaluation;
- Pacala, Canham and Silander (1993), field-calibrated spatial forest
  demography;
- Smith et al. (2014), the active C-N LPJ-GUESS configuration;
- Hickler et al. (2012), transient and equilibrium European vegetation;
- Fisher et al. (2018), vegetation demography in Earth-system models; and
- Bugmann et al. (2019), structural uncertainty among tree-mortality models.

The Pacala PDF is an image scan; it was read after local OCR.  The findings do
not transport any fitted Earth parameter to Vesper.

## Active Vesper configuration

`global.ins` selects cohort mode, 0.1 ha nominal patches, five-orbit cohort
establishment, stochastic establishment and mortality, background
establishment, the within-stand “spatial mass effect”, and generic disturbance.
The run harness overrides the source file's 15 replicate patches with five.
`build_vesper_pfts.py` converts the 100-Earth-year generic-disturbance interval
and the nutrient-free spin-up duration into model orbits; it does not otherwise
change the demographic model.  Fire remains a separate active mechanism and is
outside this note.

## Findings

### 1. Nutrient activation is implemented as a hidden stand-replacement event

At `date.year == freenyears`, cohort mode calls `disturbance(patch, 1.0)` for
every patch whenever CENTURY soil and N or P limitation are enabled.  The
comment says this is done to obtain the “right pft composition” under nutrient
limitation faster.  The call kills all vegetation, transfers its C, N and P to
litter through `Individual::kill()`, resets patch age, and skips ordinary
mortality and establishment for that orbit.

This is not merely a numerical switch.  It creates a synchronous planet-wide
stand-replacement event and a live-biomass-to-litter nutrient pulse.  Smith et
al. (2014) documents a vegetation removal for a separate chronosequence
experiment, but the paper's global spin-up description does not establish this
forced reset as part of the physical model.  Vesper's generator deliberately
rescales `freenyears`, so this source behavior is active rather than an
unreachable remnant.

The reset must be treated as an initialization accelerator, not an ecological
disturbance.  It should be removed from the accepted physical trajectory unless
matched runs establish that it is restart-equivalent after the declared
equilibrium window, including vegetation, litter, soil C-N-P and patch-age
distributions.  If retained solely to initialize a state, the pre-reset state
and reset flux must be excluded from ecological budgets and identified in the
manifest.  This is DEMO-1.

### 2. Generic disturbance is a complete-kill hazard, not a model of a Vesper process

The normal generic-disturbance path draws one Bernoulli event per patch and
orbit with probability `1 / distinterval`.  An event kills all vegetation,
leaves litter and soil intact, resets age to zero, and carries no size,
severity, climate, wind, topography or spatial-correlation driver.  The active
rescaling gives it an approximately 100-Earth-year mean return interval.

Smith et al. (2014) describes this LPJ-GUESS mechanism as a generic
patch-destroying disturbance representing, for example, windstorms or
landslides.  Those examples do not share the code's material destination:
windthrow can leave live biomass as litter without removing soil, whereas a
landslide can export vegetation, litter and soil.  Pacala et al. (1993) also
shows that disturbance size and severity can change long-term composition; the
active path has neither.  Fire is applied separately, so this hazard must not be
interpreted as an additional fire model.

The hazard supplies a successional mosaic that can materially prevent
late-successional dominance, but no Vesper evidence identifies its central
rate or physical mixture.  The honest first result is therefore a registered
model-form bracket: no generic disturbance versus the current complete-kill
hazard under identical forcing, with the latter labelled rather than silently
treated as physical truth.  A future windthrow or mass-wasting model must name
its driver, severity and C-N-P destinations.  This is DEMO-2.

### 3. “Spatial mass effect” is not geographic dispersal

With `ifsme 1`, establishment uses reproductive carbon pooled over replicate
patches in the same stand.  The patches are close enough to share propagules in
the Smith et al. (2001) conceptual model.  No seed moves between grid cells,
and `kest_bg` permits a climatically eligible PFT to appear without a local
parent.  Fisher et al. (2018) identifies absent inter-grid-cell dispersal as a
general vegetation-demography limitation.  Hickler et al. (2012) states the
consequence directly for LPJ-GUESS: propagules are assumed available wherever
climate is suitable, so projected vegetation shifts can require migration
faster than historical estimates.

This is acceptable for a long-equilibrium *potential natural vegetation*
experiment conditional on an available propagule pool.  It is not a prediction
of colonization time, range migration, isolation, refugia or island
biogeography.  The current results and consumers must carry that label.  Any
transient or geographically constrained use must fail closed until it has a
declared initial distribution and cross-cell accessibility/dispersal contract;
renaming `ifsme` in project-facing metadata avoids presenting it as that
contract.  This is DEMO-3.

### 4. Establishment and mortality retain Earth-calibrated model form, not only Earth traits

Tree establishment combines potential forest-floor productivity, Earth PFT
values for maximum establishment, shade response and propagule terms, and a
hard-coded `3 / nwoodypfts_estab` multiplier.  The multiplier prevents adding
PFTs from raising total recruitment, but makes three eligible woody PFTs the
unstated reference richness.  Seedlings are not resolved: a Poisson number of
saplings enters a cohort with initial biomass proportional to potential
forest-floor assimilation.  Fulton (1991) supports a nonlinear relationship
between juvenile growth and recruitment under explicit assumptions; it does
not validate the active coefficients or collapse of the seedling stage on
Vesper.

Mortality combines an age-dependent background hazard with a smoothed
five-orbit growth-efficiency response.  The original Smith et al. (2001)
formulation used constant background mortality tied to longevity and a
threshold growth-efficiency penalty.  The current source instead chooses a
quadratic age hazard whose integrated survival at the stated longevity is
0.001, plus a smooth response with exponent five described in code as a global
validation choice.  PFT `greff_min` values are calibrated for that smooth
option; one is explicitly documented as improving an Earth boreal PFT balance.

Bugmann et al. (2019) shows why these are structural uncertainties rather than
minor coefficients: mortality submodels with comparable historical behavior
can diverge strongly over long projections.  The accepted Vesper result needs a
pre-registered establishment/mortality model-form sensitivity, not a retuning
to the desired biomass or PFT map.  This is DEMO-4 and builds on BIO-22's time
unit correction.

### 5. The stochastic experiment has neither an exposed seed nor independent spatial streams

Every new `Stand` initializes its random-number seed to the same literal
`12345678`.  The Vesper harness exposes and records patch count but no root seed.
Different climate-dependent call paths can make streams diverge, but identical
stands begin synchronized, and there is no supported way to generate the seed
ensemble BIO-12 already requires.  MPI rank or traversal changes must not be
allowed to redefine the scientific random realization implicitly.

The harness needs one manifest-recorded root seed and stable, independently
derived streams keyed by grid cell, stand and process/replicate as appropriate.
No stream may depend on MPI rank or iteration order.  No-simulation fixtures can
test determinism, key separation and rank/order invariance before any ensemble
is run.  This is DEMO-5.

### 6. Five patches and 0.1 ha are part of the model form, not just an error bar

The harness uses five replicate patches where the source instruction uses 15.
Poisson establishment scales with `patcharea`, and stochastic cohort mortality
rounds `density * patcharea` to an integer before drawing deaths.  Nominal
patch area therefore changes demographic granularity and extinction behavior;
increasing only `npatch` reduces Monte Carlo noise but does not test the same
thing.  Fisher et al. (2018) likewise cautions that patch meaning depends on the
disturbance and demographic scale being represented.

Current retained output gives grid-cell/PFT mean density, biomass, cover and
LAI; mean species height is available in the output module but disabled by the
harness, and cohort age structure is sent only to an interactive plotting path.
There are no retained recruitment, mortality-cause, generic-disturbance or
patch-age distributions with which to explain an equilibrium mean.

After DEMO-5 supplies controlled streams, a small pre-registered convergence
design must vary root seed, patch count and patch area separately and assess
the full equilibrium window.  The accepted artifact must retain enough
age/size, patch-age, recruitment, mortality-cause and disturbance flux
diagnostics to attribute the spread, while C-N-P closure remains under BIO-14.
This is DEMO-6; the simulations require explicit permission.

## Interpretation boundary and order

The present configuration can support a labelled equilibrium potential-natural-
vegetation experiment after the forced reset and time-base defects are closed.
It cannot support claims about transient migration or the physical frequency of
non-fire stand replacement.  DEMO-1 and DEMO-5 are source/provenance work that
can begin without a model run.  DEMO-2 through DEMO-4 define model-form
brackets.  DEMO-6 supplies the diagnostics and measures stochastic convergence
only after the correctness gates and with explicit permission to execute LPJ.
