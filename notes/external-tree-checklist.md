# What to check first in an external model tree

Vesper is a fictional planet and this project simulates it, which is why every
external Earth-system tree acquired here has to be audited for the planet welded
into it. This is the checklist that audit should run, derived from
`notes/external-model-survey.md`, which read some twenty independent trees and
hit the same classes of defect in most of them. Every item below cites the
sections that evidence it and names the trees it was observed in. Nothing here
is a rule about how to model; it is a list of where to look, in what order, and
what finishing an item looks like.

## How to read the support column

An item observed in four independently authored trees is a class. An item
observed in one is an anecdote that has not yet been contradicted. The `models`
column names them and the count is the number of INDEPENDENTLY AUTHORED trees,
counted by authorship rather than by directory: cGENIE counts once even where
GOLDSTEIN, BIOGEM, EMBM and ENTS each supply an instance, but SICOPOLIS and
PALADYN count separately from CLIMBER-X because they are separate upstreams
vendored into it. Items with a count of one are quarantined in the last section
rather than mixed into the list.

The `clean negative` column is what completing the item looks like. An item
without one cannot be closed, only abandoned, and the survey's own practice was
to record negatives so the check was not repeated
(`notes/external-model-survey.md` sections 9d, 25, 47d, 53b).

## The order, and why it is this order

The audit is not uniform over the tree and should not be run as though it were.
`notes/external-model-survey.md` section 40c is the finding that sets the order:
in cGENIE the biogeochemical parameters are namelist-exposed and the physical
scaling constants are compile-time `PARAMETER`s, so an audit that samples evenly
spends its effort where the model is already configurable. The same split
appears in CLIMBER-X, where the gas box models take namelists and
`src/main/constants.f90` welds radius, rotation and gravity at compile time
(section 7a). Go to the physics core and the unit conversions first.

Tier A is grep. It costs minutes, needs no understanding of the tree, and finds
the highest-yield class in the survey -- the planetary constant stated more than
once. Tier B is one file open per item. Tier C requires reading a use site, and
it is where the survey's wrong readings were caught: a declaration
mischaracterises its own value often enough that Tier A and B results are leads
rather than findings. Tier D requires reading a subsystem and is worth doing
only once a tree is a live candidate.

Two orderings inside the work follow from the same document and are worth
stating before any of it starts. Correctness-preserving work against a
regression suite must come BEFORE any change to planetary parameters, because
the suite's reference values are computed at Earth's parameters and moving them
destroys the comparison as its first act (section 13c). And a rule that reads
only what an open task row already names cannot find a blind spot, because a
blind spot is by definition a thing no row names (section 9); reserve part of
the reading for the tree's own shape.

## Tier A: grep, before opening anything

| # | check | clean negative | models | sections |
| --- | --- | --- | --- | --- |
| A1 | Count the DISTINCT statements of the year and of the day, in every language in the tree -- Fortran, config, XML, namelist, shell, MATLAB. Enumerate value, file, line and whether each is compile-time. | The year is stated once and the day is derived from a rotation period rather than a literal; or N statements are enumerated with their values and each disagreement is named. | 5: cGENIE, CLIMBER-X, FATES, CAABA/MECCA, LPJmL | 40a, 40b, 46b, 42e, 54b, 48e, 51c |
| A2 | Locate the planetary constant block and classify every member: compile-time `parameter`, runtime namelist, or derived. | Every planetary constant is located and classified, and each dependent member is shown to be DERIVED after the read rather than independently settable. | 8: CLIMBER-X, cGENIE, Isca, SOCRATES, ClimaParams.jl, CaMa-Flood, SICOPOLIS, FATES | 7a, 10e, 9d, 21, 5, 52b, 45d, 54d |
| A3 | Grep the Earth numeric literals directly, not the names: `9.81`, `9.80665`, `6.37e6`, `6371`, `7.2921e-5`, `5.15e18`, `1361`, `1368`, `101325`, `1013.25`. Constants escape a named block by being retyped. | Each hit is either absent, or accounted for in A2's classification with its second statement identified. | 6: cGENIE, CLIMBER-X, SICOPOLIS, CaMa-Flood, FATES, CAABA/MECCA | 10e, 40b, 43c, 45b, 52b, 54d, 48d |
| A4 | Grep for irradiance-to-photon and band-split scalars: `4.6`, `PAR`, `umol`, `quanta`, `550`, `400`, `700`, half-saturations in `W m-2`, e-folding depths in metres. | Every such scalar is located with its integration window AND its weighting spectrum named; or the quantity is supplied by the driver rather than held as a constant. | 7: LPJ-GUESS, cGENIE, FATES, ClimaLand, CLIMBER-X, SPEEDY, SOCRATES | 24, 40d, 54d, 12c, 43c, 8a, 44c |
| A5 | Grep for latitude constants and named region ids in BRANCH position, not just in data: `lat >=`, `lat <`, `40.0`, `30`, `60`, and any region-id enumeration read from a mask. | No latitude branch and no region id sits in a path the model reads at runtime; or each is enumerated with what it selects. | 6: LPJ-GUESS, FATES, CLIMBER-X, SICOPOLIS, BIG-MITgcm, LPJmL | 25, 34b, 54b, 43d, 45d, 4, 51c |
| A6 | Grep for compile-time grid and topology bounds: `#define`, `parameter` dimensions, `NLAT`, maximum-island macros, hardcoded cell-boundary indices. | Every bound is located, and each is shown either to be a namelist quantity or to be a rebuild whose trigger condition is stated. | 3: cGENIE, CLIMBER-X, ExoPlaSim | 10d, 30d, 43d |
| A7 | Grep for the tuning and validation surface: cost functions, bias corrections, observational target files, fitted-coefficient tables, known-good suites. | The scoring targets are named and each is stated as inert here, or as a separable build target that can be declined. | 6: cGENIE, CLIMBER-X, BIG-MITgcm, BLAZE, ArcSDM, Global NEWS 2 | 46e, 13c, 3f, 4, 51d, 49d, 52c |

## Tier B: one file open per item

| # | check | clean negative | models | sections |
| --- | --- | --- | --- | --- |
| B1 | Open the non-dimensionalisation or unit-conversion module -- the file that turns physical quantities into the model's internal ones. This is where the planet is welded even in trees whose parameters are otherwise exposed. | Every scale is located and each is either a planetary quantity that is configurable, or a modelling choice with no planetary content, and which is which is stated per scale. | 3: cGENIE, CLIMBER-X, ESMF | 10e, 43c, 55b |
| B2 | List every coefficient whose UNITS contain a length, a pressure or a head in the denominator, and recompute each into its SI material form. A unit can absorb a planetary constant and then hide it: head in metres absorbs `rho*g`, bar per metre absorbs `rho*g`, an atmospheric thickness in metres absorbs the whole atmospheric mole inventory. | Each recomputed coefficient matches a published material property with no residual `rho*g`; or the absorbed factor is identified and the corrected value stated. | 4: SICOPOLIS, CLIMBER-X, ClimaLand, cGENIE | 45b, 43c, 41b, 46b |
| B3 | Classify every per-day and per-year rate into three kinds: physical and kinetic, where the day is 86400 s and must NOT be rescaled; biological and entrained, where the value was calibrated under a 24-hour light cycle and whether it moves is a modelling decision; and not temporal at all, where the day sits inside a unit that is really spectral or spatial. | Every rate is assigned to one of the three, and every entrained one is recorded as an open choice rather than silently converted. | 4: cGENIE, FATES, LPJ-GUESS, LPJmL | 40d, 54c, 51c |
| B4 | Read the comment beside every constant found in Tier A, then check it against the value. The comment is usually the only place a value's provenance appears, and it is the place the provenance is wrong. | Each comment that asserts a property has been checked against the value and the use site, and each mismatch is recorded. | 5: ExoPlaSim, CLIMBER-X, cGENIE, Isca, ClimaLand | 11c, 42e, 42c, 10a, 9d, 41b |
| B5 | List every clamp, floor, cap and limiter, and evaluate where each binds at this planet's parameters rather than at Earth's. A dimensional limiter sitting on a quantity that scales steeply with gravity binds far more often and truncates the effect being studied. | Each bound is located and its binding regime here is stated, with any that bind in normal operation named. | 5: SPEEDY, SICOPOLIS, CLIMBER-X, ExoPlaSim, cGENIE | 8a, 45d, 43c, 29c, 32a, 20a |

## Tier C: read a use site

Tier A and Tier B produce leads. This tier is what turns them into findings, and
it is the tier the survey most often needed: characterising a value from its
declaration rather than from its use gave a wrong account repeatedly, in trees
that were otherwise being read carefully.

| # | check | clean negative | models | sections |
| --- | --- | --- | --- | --- |
| C1 | For every constant from Tier A, open one use site and confirm what it does. Ask specifically whether a later stage OVERRIDES the declaration, and whether the parameter is reachable at all under the default switch. | Every carried constant has one use site read, and any declaration overridden downstream or unreachable by default is marked as such. | 4: cGENIE, ExoPlaSim, CLIMBER-X, rokgem | 31a, 31b, 11c, 29b, 46a |
| C2 | Test every name that carries a condition. A predicate that reads as a property of its subject can be a property of the subject AND something else; a component named after an algorithm can implement its opposite. A call site cannot be audited for a dependency by looking at the call site. | Every predicate and every named component on the audit path has been opened rather than read by name. | 5: LPJ-GUESS, Landlab, cGENIE, CLIMBER-X, ClimaLand | 25, 53a, 27, 19, 7a |
| C3 | Separate DECLARED capability from DEMONSTRATED capability: an annotated parameter range with nothing shipped near its top, a code path with no namelist and no test entry, an abstract type with one implementation, a documented class hierarchy the code does not have. | For each capability the tree is being acquired for, either a shipped configuration or a test exercises it, or it is recorded as unexercised and the build check is named. | 6: cGENIE, CAABA/MECCA, ClimaLand, Landlab, ExoPlaSim, LPJmL | 10d, 48c, 46c, 37a, 53e, 6, 51a |
| C4 | Find the fail-open branches. A backwards-compatibility branch that reverts to an Earth constant, a list-directed read that desynchronises on a short record, a fill that invents gradient, a fallback with no threshold at which it refuses, an extrapolation past a declared validity range. | Each fallback found either refuses by default, or is recorded with the condition that triggers it and what it silently substitutes. | 6: cGENIE, CLIMBER-X, SOCRATES, SICOPOLIS, Landlab, ESMF | 10e, 31d, 47a, 44d, 45d, 53b, 55c |
| C5 | Look for duplicate live state and second constant sets: two copies of one state array, a dead module carrying its own inconsistent constants, two gravities in one model, two densities for one substance. | Every duplicated quantity is located and the authoritative copy named, or the duplicate is confirmed dead with its call site checked. | 3: cGENIE, CLIMBER-X, CaMa-Flood | 35a, 10a, 45d, 43a, 52b |

## Tier D: read a subsystem, once the tree is a live candidate

| # | check | clean negative | models | sections |
| --- | --- | --- | --- | --- |
| D1 | Read the geography preprocessing that runs before the model sees terrain. Depression filling, lake deletion, minimum-depth flooring, elevation clamping and boundary treatment all destroy features silently, and closed basins are the feature this project's terrain produces in quantity. | Every preprocessing step that alters the mask, the bathymetry or the elevation field is enumerated, with what it removes and whether it can be declined. | 3: BIG-MITgcm, CLIMBER-X, Landlab | 4, 47a, 47c, 53a |
| D2 | Read the coupling boundary field by field: units, index base, array orientation, grid convention, and the CALENDAR DIMENSION of the exchange arrays. A seasonal array dimensioned by day-of-year is a calendar assumption in the coupling rather than in a namelist. | Every exchanged field has its units, index base and dimension checked against both sides, and any calendar-dimensioned array is sized against this world's orbit. | 4: cGENIE, SICOPOLIS, BIG-MITgcm, CLIMBER-X | 10c, 30b, 30d, 45d, 4, 7a |
| D3 | For every accelerator, ask what it holds fixed, what it declares itself invalid for, and whether any single number measures the drift it introduces. A cycling accelerator that returns to the full model periodically is valid under weaker assumptions than a single long jump. | Each accelerator's validity statement is located, and either its drift scalar and guard thresholds are named or the accelerator is recorded as a fixed schedule with a declared structural cost. | 3: CLIMBER-X, cGENIE, ExoPlaSim | 3c, 26, 35, 46c |
| D4 | Run the conservation identity, which is the one class of check here that has a right answer in advance: flow accumulation against upstream area, remap integrals across a grid pair, a closed water budget. Register the tolerance before running it. | The identity closes at a tolerance registered in advance, and the result is recorded so the check is not repeated. | 3: CLIMBER-X, Landlab, ESMF | 47d, 53b, 55a, 47b |
| D5 | Check whether the tree ships an offline driver or a fake-component pattern. A component that runs against prescribed forcing without its host can be TRIED before anything is adopted, and that asymmetry decides sequencing independently of merit. | The tree either ships a standalone driver or a per-component file-reading stub, or it is recorded as requiring its host, with what adopting the host costs. | 3: CLIMBER-X, MARBL, Isca | 7a, 16b, 9c |
| D6 | Where a structure has spectral bands, check that the bands hold different numbers. A two-band array initialised from a scalar carries a star's spectral shape nowhere, so a spectrum fix upstream does nothing for it; a single attenuation coefficient cannot express that a redder star's photons are removed faster. | Every band array is confirmed to hold distinct per-band values anchored to a stated broadband target, or the flat ones are enumerated with which are used as written. | 4: ExoPlaSim, cGENIE, CLIMBER-X, ExoCAM | 11a, 11c, 12c, 43c, 22, 29b |

## The method notes that generalise

These are search techniques, and each cost a wrong claim or nearly did.

**Characterise a value from its use site, never from its declaration.** A
positional table with no header is read correctly only from the Fortran that
consumes it, and reading the archival data files instead produced a conclusion
that did not fire (section 31). A declaration can be overwritten downstream, so
a constant that looks binding contributes nothing to the surface it names
(sections 11c, 29b). A namelist parameter can be silently ignored under the
default scheme string while the run log prints it as though in force (section
46a). Tier A and Tier B produce leads; only Tier C produces findings.

**A comment is evidence about intent and never about the value.** A model can
call 87658 seconds the actual seconds per day (section 42e). Two arrays can be
commented "spectral weighted" while holding one number in both bands (section
11c). A namelist can print the published range beside a tuned value that sits
outside it at both ends (section 42c). A comment can concede an unresolved
out-of-bounds and hedge it (section 10a). Read comments as the fastest available
index to provenance, then check every one against the value.

**A name can hide a condition.** A predicate reading as a property of the stand
is a property of the stand AND its latitude, so a call site cannot be audited for
latitude dependence by looking at the call site; the predicate has to be opened
(section 25). A component named for a priority-flood router fills every
depression where this project's implementation of the same name preserves them
(section 53a). A variable in a block named for relaxation is a numerical
under-relaxation and not a physical restoring (section 27). A file named for a
free surface diagnoses sea surface height and does not replace the rigid lid
(section 7a).

**Search the whole tree in every language it contains.** The authoritative
statement of cGENIE's year length is a bash literal in a run script, appended to
the generated configuration and overriding the config file it is built from
(section 46b). A search confined to the model's own language misses it.

**Count the uses, not the declarations.** FATES carries a compile-time day count
used 146 times and a runtime host value used 21 times, and the ratio is what
says which one governs -- and the two meet inside one physical process, so the
defect is a mass-balance error that cannot appear on the calendar it was written
for (section 54b).

**A unit is a place a constant can hide.** Gravity enters SICOPOLIS's
Clausius-Clapeyron gradient with no `9.81`, no density and no `G` on the line,
because the coefficient is per metre of ice rather than per pascal (section 45b);
the same form appears in an ocean equation of state as a bar-per-metre literal
(section 43c) and in plant hydraulics as a potential stored in metres of head
(section 41b). No search for constants finds any of them. B2 exists because the
grep cannot.

**Where several independent trees make the same choice, the departure is what
needs the argument.** Three codes take the same stellar spectrum library and the
same three parameters, which is a convergence rather than a finding and means
nothing has to change (section 9a). Two independent implementations estimate the
water table from column water content rather than from a lateral solve, which
makes the column estimator normal practice and the lateral solve the departure
(section 36b). Two of three land models omit snow compaction, so its absence is
conventional at that complexity level rather than anomalous -- and the separate
question of whether the Earth-magnitude judgement behind the omission survives
this planet's gravity is not answered by the convergence (section 23d).

**Disagreement between published schemes is a free bracket.** Where a tree ships
several parameterisations of one process with their constants, the spread across
them is a structural-uncertainty estimate that costs nothing to compute and does
not need the schemes to be adopted (sections 17, 31c, 29c, 32a). Two of those
turned out to be worth more as a magnitude than the adoption question was worth.

**Record negatives.** A cheap check that comes back clean is worth writing down
so it is not run twice: a heating-rate constant checked and found to use the
model's own gravity (section 9d), a constant scan whose results already live
elsewhere (section 25), a flow accumulation confirmed correct in two independent
implementations after a third was found wrong (sections 47d, 53b), two modules
that resolve to nothing (section 46c).

**Purpose does not exempt a tree from the audit.** An exoplanet-capable code
parameterised its planetary constants at runtime and still let a derived
radiative constant through as a compile-time product of Earth gravity, Earth
specific heat and Earth day length, with the author's uncertainty left on the
record (section 9d). And an exoplanet-built code is not automatically the better
fit: what makes this world unusual is its star rather than its atmosphere, so an
Earth-built code with a substitutable stellar spectrum will usually fit better
than one aimed at exotic compositions (section 15c).

## Single-model observations, carried but not promoted

These appear once each in the survey. They are recorded so a second instance
promotes them rather than being read as new, and none is strong enough to spend
audit time on by itself.

| observation | model | section |
| --- | --- | --- |
| An accelerator whose entry is guarded on a RATE and whose extension is guarded on an accumulated error integral, with a hard cap by default | cGENIE | 35b |
| Land and ocean fraction threshold made latitude-dependent because a coarse cell's ocean fraction has different consequences where straits carry through-flow | CLIMBER-X | 47f |
| A depth QUANTILE rather than a mean handed to the ocean model per coarse cell | CLIMBER-X | 47f |
| Runoff delivered to ocean cells with embayment-first spreading over an ordered neighbour list | CLIMBER-X | 47e |
| Topography filter width computed in cells from a distance in km via Earth's radius, silently changing effective resolution at another radius | CLIMBER-X | 47f |
| A polynomial emulator anchored to absolute mole inventory rather than mixing ratio, evaluated outside its training range with no truncation | cGENIE | 46a |
| Second-order remapping silently degrading to first order where a source cell has too few unmasked neighbours, which is what a coastline looks like | ESMF | 55e |
| A prognostic field absent from the restart, so a restarted run reinitialises it | CLIMBER-X | 47b |
| A diagnostic-only block carrying the largest concentration of Earth astronomy in its module, deletable without touching a result | LPJmL | 51c |
| Normalisation by a single cell's maximum, making a field incomparable across builds | ArcSDM comparison | 49d |
