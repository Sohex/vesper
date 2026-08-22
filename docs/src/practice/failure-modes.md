# How this project goes wrong

Almost every bug in this list was found late, by comparing two things that
were supposed to agree: array shapes matched, fields looked plausible, runs
converged, and the answer was wrong. The few that did announce themselves
(classes 3, 11, 21) did so by luck, not by design. They are recorded by *class* rather
than as a changelog, because the classes recur and the instances do not.

Run `python scripts/check_consistency.py` before an expensive run and after
changing `source_build`. It mechanises the checks that come out of the classes
below.

---

## 1. One quantity, several consumers, a fix that reaches some of them

The most expensive class here, twice over.

**`mrro`.** Four scripts integrated catchment runoff. `pedology` established that
`mrro` is river-routed net divergence rather than local generation and switched
to P - E, documenting why. `surface_water.py` reached the same conclusion
independently and also switched. `export_carve_list.py` and `carve_verdict.py`
were never updated, and the first of those is the one whose output leaves the
project and changes the terrain. Two independent corrections of one
misunderstanding, neither of which propagated to the consumer that mattered most.
Cost: an entire terrain build.

**`field_lon`.** The longitude-convention fix landed in `basin_means` and two of
its three callers. The third was `export_carve_list.py`. Again the one that
changes the terrain.

**The rule.** When a shared quantity is redefined, enumerate its consumers before
declaring the fix done -- `grep` for the symbol, not for the bug. A fix that lands
where the bug was noticed is half a fix. Prefer one canonical accessor over four
call sites that each decide for themselves.

## 2. An optional argument whose absence means "do the wrong thing"

`basin_means(..., field_lon=None)` defaulted to the unremapped path. So the fix
for the antipodal bug reproduced the antipodal bug in any caller that had not
been updated -- the original error wearing the shape of its own correction.

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
switches and a digest of every surface file -- but not the stellar spectrum, and
it rounded flux to hundredths. Two near-misses followed: a `k25v` baseline re-run
would have been written into the completed `k2` run's directory, where
`finalize()` takes `sorted(glob("MOST*"))[-1]` and copies out the last match; and
`0.945` and `0.94` both resolved to `s094`. Only the geography digest separated
them, which was luck rather than a guard.

**The rule.** Anything physical that is not in the name is a collision waiting for
the run that changes it. When adding an input that changes the answer, add it to
the identity in the same commit.

## 5. Provenance written once and never recomputed

The carve list header reported 749 carved for a file containing 1,838, and
described "0.96 S-Earth, 292.88 K" two baseline re-runs after both had moved. A
year length of `189.6145` sat as a literal while the baseline ran at a different
flux. `climatology_s096` was described as current long after it was not.

**The rule.** Derive provenance from the artifact at write time. If a header
states a number, compute that number in the same function that writes it.

## 6. A claim recorded as settled that measurement contradicts

- The energy residual was recorded as "a fixed offset, -0.455 +/- 0.010, which
  rules out every state-dependent candidate". Measured the same way across four
  converged runs it spans -0.400 to -0.489, a band some four and a half times
  wider than the recorded +/-0.010, putting the ruled-out candidates back in
  play.
- `exoplasim/notes/lake-representation.md` prescribed a large `dwmax` because "a large full
  bucket evaporates at open water's rate". Wetness reaches 1 above 40% *of*
  `dwmax`, so a deeper bucket needs proportionally more water; and routed river
  water never re-enters the evaporating bucket at all, so no setting sustains a
  lake.
- A stellar-spectrum correction was estimated at 0.4 to 0.7 K and measured at
  0.04 W/m2, because the estimate was computed on a world with far more snow.

**The rule.** A recorded conclusion is evidence about what was true of the cases
it was measured on. Before relying on one, check whether the current
configuration is inside that envelope. Re-measure rather than inherit.

## 7. Unguarded arithmetic at a physical boundary

Switching catchment runoff to P - E made it negative for 57% of basins, which is
physically meaningful -- those catchments evaporate more than they receive -- but
the lake-area solve `A = R*C/(E - P + R)` was not written for it and reported a
lake covering -6.7 e15 percent of the planet. Earlier, a spill-level fallback
returned `inf` and produced an infinite basin capacity.

**The rule.** Clamp at the physical bound where the quantity is defined, not
where it is used. A catchment delivers zero or more, never less.

---

## The pattern behind the pattern

Almost every entry here was found by comparing two artifacts that should have
agreed, and almost none by reading code. The productive habit is to look for
quantities computed two ways and check them against each other: mesh against
grid, our routing against the exporter's, the Penman closure's implied
surface temperature against the model's own SST (the evaporation ratio itself
can only differ; class 17), one build's composition against another's. Where a cross-check exists,
these bugs surface in minutes. Where none exists, they survive until something
downstream looks strange.

So when adding a component, add the cross-check with it, and prefer the check
that would have caught the last bug.

## 8. Reading an intermediate configuration as the world

Three conclusions in this project have been right in mechanism and wrong in
magnitude. The dust figure was taken from the build sitting in front of me
rather than from the configuration the pipeline is converging on; the
phosphorus conclusion came from following one transport pathway where two run
in opposite directions.

- **The carbonate-silicate thermostat** was described as structurally weak
  because endorheic drainage withholds alkalinity from the ocean. The mechanism
  is real; measured, the efficiency is 0.906 and rising with each carve. The
  error was using an area share where a weathering-weighted share was needed --
  endorheic land is dry, which is *why* it is endorheic, so it weathers less per
  unit area than its extent implies.
- **The phosphorus gradient** was described as P-poor uplands against P-rich
  basin floors. That is the hydrological leg only; aeolian transport from dry
  lake beds runs the other way and partly refills the uplands, which on Earth is
  how the Amazon is supplied.
- **The dust source area** was taken as 26.6% of land, which is the *uncarved*
  build -- the hyper-arid limit in which no basin has had its outlet cut. Carved
  and excluding fill that sits under a lake, it is about 12%. The bound fell from
  -4.6 K to -2 K.

**The rule.** A pre-carve build is a *limit*, not a state. So is any figure taken
before the lakes are solved, before a verdict is applied, or before the loop this
project is built around has run. Before quoting a geography number, ask which
iteration it belongs to and whether the pipeline expects to move it. The
uncarved numbers are the ones physically present in `source/` at the start of
every cycle; presence there is availability, not currency.

Corollary: this world is unusual enough in its drainage that reasoning from one
transport pathway is reliably insufficient. Where water does something
surprising here, check whether wind or sediment does the opposite.

## 9. A number taken from a citation rather than from the paper

Parent-rock phosphorus was entered from a second-hand table that cited Hartmann
et al. 2014. The values were wrong: they came from that paper's "P content in %
relative to Cat + SiO2 release" row, a release-normalised ratio, not from its
rock-content row. Every value was low by roughly half, the units were
misidentified as P2O5 weight percent, and the ppm conversions applied on top of
that were meaningless.

Two things made it look sound. The class-to-class *contrast* was nearly right --
2.33x reported against 1.94x actual -- because a ratio row and a content row rank
the classes the same way. And it cross-checked plausibly against an independent
compilation, which agreed on the contrast for the same reason.

**The rule.** A contrast surviving a cross-check is not evidence that the
magnitudes are right; ratios are preserved by exactly the transformations that
corrupt values. When a number will be integrated rather than compared, read the
source table, not a citation of it. The tell here was available and ignored: the
same paper prints P2O5 and elemental P as adjacent rows differing by 0.4364, so
"which row is this" was a question the data itself was asking.

**And the correction was itself half wrong**, which is the more useful part. The
values were replaced with the elemental rock-content row on the grounds that the
release row is "not a rock content". True, but the intended use was to multiply
by a weathering intensity, and Hartmann states plainly that P release is computed
as the release-relative content times the cation-and-silica flux -- so the
release row was the right one for that purpose all along. The agent had the right
numbers with the wrong description; the correction had the right description of
the wrong numbers.

The deeper error was never asking what the number was *for*. Both rows are
correct data; which one is correct depends entirely on what multiplies it. The
config now carries both, labelled by the equation each belongs to.

Corollary for delegated work. The agent that produced these flagged its own
uncertainty accurately -- it said the units were "likely, not confirmed" -- and
that flag was correct and was not acted on. A caveat carried forward verbatim is
not the same as a caveat resolved.

## 10. A bound asserted in the wrong direction

The Lacis-Hansen ultraviolet weight was computed from a 4965 K blackbody as
0.469 and recorded as a **floor**, on the reasoning that a real star adds
chromospheric ultraviolet that a photospheric model lacks, so the true value
could only be higher. The measured value, from observed IUE spectra of the star
in question, is **0.335** -- lower, not higher.

The reasoning omitted the larger term. Ultraviolet **line blanketing** in a real
stellar atmosphere removes far more flux than the chromosphere puts back, so a
blackbody *overestimates* a cool star's ultraviolet rather than bounding it from
below. Two effects, opposite signs, and only the smaller one was considered.

**The rule.** Before calling an estimate a bound, enumerate what could push it
the other way. A bound asserted from one mechanism is a guess wearing a stronger
word, and it is worse than an unqualified estimate because it discourages the
check that would correct it. If both directions cannot be enumerated, say
"estimate" and leave it undefended.

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
zeroes the gridcell's vegetation **silently** -- zero LAI, zero AET, no crash. And because
water is fractional, thinner soil reads as relatively *wetter* for the same
absolute water, so a single cell can move either way through PFT competition.
Trust a controlled sweep, never one cell.

## One quantity, two meanings, several times the value

"The endorheic share of land" names two different measurements that differ by
a large factor:

- the share of land *inside* a preserved basin (`is_endorheic` cell area) --
  the lithology and thermostat figure;
- the share of land that *drains* to a closed basin -- the hydrography figure
  and the one the carve verdict is about.

Both are correct, both move with the carve, and their current values live in
`world_state.json` under separate names. A basin's catchment is far larger than its floor, which is the
whole reason a small area of fill can decouple a large share of weathering. They
are named apart in `world_state.json`; quote the name, never "the endorheic
share".

## 11. A patched source, and one binary per configuration

ExoPlaSim compiles a separate executable for every (resolution, layers, ranks)
triple: `most_plasim_t42_l10_p16.x` and `most_plasim_t21_l10_p8.x` are different
files built from the same source at different times. **Patching the source and
rebuilding rebuilds only the configuration you are running.** Every other binary
on disk keeps the old code, silently, until something asks for it.

Found moving from the T21 flux bracket to a T42 bootstrap. The ozone
band-weight patch adds `o3uvw` and `o3visw` to `radmod_nl`. The T21 8-rank binary
had been rebuilt and carried them; the T42 16-rank binary was three days older
and did not, so it died in `radini_` with

    Fortran runtime error: Cannot match namelist object name o3visw

Of the six executables present, exactly one was newer than the patched
`radmod.f90`. The other five would all have failed the same way.

**Why this one is survivable, and what would not be.** The failure is loud
because the patch added a *namelist key*: an old binary cannot parse a key it was
not compiled with, so it aborts. A patch that changed only the *value* or
*meaning* of an existing quantity would not abort. It would run to convergence
and produce a result from the unpatched physics, and nothing in the run manifest
records which binary produced it -- `run_id` names the geography, spectrum and
flux, all of which would be identical.

**How to apply.** After patching the model source, treat every executable older
than the patched file as stale. Before an expensive run at a new resolution or
rank count, compare the binary's mtime against the source it should contain.
Deleting a stale binary is cheap and it is rebuilt on demand; discovering the
staleness after a converged run is not.

The general form is the one this file keeps returning to: a build product whose
identity does not record what went into it. The same reasoning put the geography
digest, the spectrum and the flux into `run_id`. The producing binary does NOT
join them, because rule 6 says an id encodes nothing: it is recorded in
`run_manifest.json`, once for what the run was prepared with and again on every
segment, since a run long enough to span a rebuild has segments no single
executable produced. `binary_manifest.json` records the compiler and the flags
beside the source shas for the same reason, so an executable whose sha has moved
is attributed rather than merely noticed.

**A bed carries the binaries it was STARTED with, and `--verify` cannot see it.**
Two questions look alike and are not: "are the installed binaries current", which
`rebuild_binaries.py --verify` answers, and "is this run using them", which it
cannot. A working directory copied out of `exoplasim/runs/` holds executables
frozen at the moment that run began, so a bed built from a recent run measures a
model source nobody named -- it runs, it is self-consistent, and it is about the
wrong code. Hash the executables a bed actually contains against
`binary_manifest.json` before measuring anything with it;
`reproducibility_matrix.py:check_binaries` is the guard, written after this cost
a full measurement pass.

**A namelist variable that is not broadcast is absent everywhere but the root,
and nothing says so.** `radmod_nl` and its siblings are read on NROOT only and
every variable is then broadcast by hand. Add one to the namelist, forget the
`mpbcr`, and it keeps its default on every other rank -- so a term whose default
means "off" runs on one rank's latitude rows and nowhere else. CLIM-42 lost most
of a day to this: the model does not warn, the run completes, `plasim_diag`
echoes the value because that too is printed on NROOT, and the symptom is a term
that comes out weak rather than absent, which reads as a physics problem and was
chased as one through two wrong hypotheses and a fetched paper.

**What finds it is a PER-RANK diagnostic, and getting one is harder than it
looks: only NROOT's writes to `nud` reach `plasim_diag`, so the obvious
`if (mypid == NROOT)` samples rank 0's latitude rows and nothing else -- all
polar at T42 on 8 ranks. Write to `70+mypid` instead. The tell is an exact
zero.** A
global mean cannot distinguish a term that is weak from a term that is off over
most of the domain; `+0.0000` on seven ranks of eight and a real number on the
eighth is not something physics produces. The same class is
`notes/audits/nlowio-collective-deadlock.md`, a collective placed behind an
unbroadcast `nlowio`, and PHYS-9, a key applied at prepare and not per segment.

**Changing the rank count changes the answer.** Reproducibility at a fixed rank
count does not survive changing it: 8 ranks and 16 ranks integrate one restart to
restarts differing in 83 of 199 records, at round-off and growing, because the
decomposition changes which partial sums are formed in which order. So an A/B
must hold the rank count fixed, and a bracket run at one rank count is not
comparable to a control at another. `notes/audits/model-reproducibility.md` has
the measurement.

## 12. A continuation that re-derives the physics from config

`run_exoplasim.py --flux-ratio 0.91` prepares a run at 0.91 and stamps it in the
manifest. `continue_exoplasim.py` then took its flux from
`config/planet.yaml:orbit.baseline_flux_earth`, so continuing that run integrated
it at the baseline instead, and the manifest went on saying 0.91.

Found on the first flux bracket of the baseline re-run. Nothing failed:
the run converged cleanly, passed all seven convergence criteria, and reported a
temperature. It was caught because two runs whose fluxes differed by 0.035
converged to the same temperature, which no amount of internal consistency can
make physical.

**The drift check could not have caught it, and that is the general shape.**
`config_drift` compares the config against the manifest's `source_config`, and
both said 0.945. The flux that made this run different from the baseline never
lived in the config at all: it arrived as a command-line argument at prepare
time. So a value that is part of a run's identity was being re-derived from a
source that had never held it.

**The rule.** A continuation continues; it does not re-decide. Anything that
defines the run has to be read back from the run, and an argument that disagrees
with what is on disk is an error rather than an override. `continue_exoplasim.py`
now takes the flux from `physical.flux_ratio` and refuses a `--flux-ratio` that
does not match it.

**The same shape one level down: a guard that compares the reference and not the
referent.** The resume guard's allowlist calls `star.spectral_type` inert, and
that is CORRECT about the config key -- nothing passes the spectral type to the
model. It is wrong about the artifact, because `build_stellar_spectrum.py`
writes `k25v.dat` from it, in place and under a name that never moves, so the
whole radiative input can be replaced while every recorded value is unchanged.
Tightening the allowlist would not have helped: the spectrum can be regenerated
with no config edit at all. The fix is to compare the FILE, which the manifest
now records by sha256 the way `biosphere/generated/vesper_provenance.json` (generated by
`build_vesper_header.py`)
already did. Recorded as CONS-3. The general form: when a guard
compares a NAME, ask what the name points at and whether that can move under it.

**What it cost, and what to check.** Forty-five orbits, about two hours, and a
run whose manifest had to be corrected rather than trusted. When a segment
finishes, the namelist
in the run directory is the ground truth for what was integrated; the manifest is
a claim about it.

## 13. State that lives in the checkpoint rather than in the inputs

`run_exoplasim.py --restart-from` seeds a new run from another run's restart
file, and it is refused whenever the surface fields differ. The reason is that
`landmod`'s `landini` takes `dwmax`, `dz0clim` and all three `dalbcl` bands from
the restart when `restart > 0`. The `.sra` files are read only on a cold start.

So a seeded run with new soil water, new roughness or a new albedo would ignore
every one of them and silently reproduce its parent, while the run directory,
the manifest and the surface-field hashes all recorded the new values. Nothing
would look wrong.

Hit again, seeding the baseline from the bootstrap to skip 86 orbits
of spin-up: the guard refused it on the surface codes differing by 229. The guard
existed because this had been found once already; what did not exist was the rule
written anywhere a person planning a run would read it, which is why it was
attempted a second time.

**The rule, and it is not "restarts are useless".** They are valid across runs
exactly when the SURFACE is unchanged and the FORCING differs, because flux and
CO2 are namelist parameters rather than restart state. That is precisely the flux
bracket, and the 0.91 point was cold-started for want of noticing, which cost
about 45 orbits. They are invalid the moment any boundary field moves, which is
every step of loop A, because each one adds a surface field.

**The general shape.** A checkpoint is a second source of truth for state that
also has an input file, and the two silently disagree in favour of the
checkpoint. Ask of any resume mechanism which fields it restores rather than
re-reads. The same question applies to `configure()`, which clears surface fields
when handed a landmap, and to anything else that decides between a file and a
saved state.

## 14. A bad bin in a shared artifact, inherited silently by everything

Every orbit of `run_b014469b8091` writes a corrupted first output bin. In `spd`,
`ua` and `va` the lower troposphere is inflated, worsening monotonically
downward: 1.02x at the model top, 2.5x at sigma 0.30, and **7.5x at the bottom
level**. It is systematic and not a restart shock -- orbits 86, 87, 88 and 90
give 7.50, 7.57, 7.54 and 7.58 -- so it is a property of how the first output
interval of each model year is accumulated. The snapshot climatology built from
the same run is clean, which localises it to the regular binned output.

Found while chasing a dust emission flux that came out 1700x Earth's.
Emission goes as roughly u* cubed above a threshold, so the bad bin was
**100.00%** of the annual total and every other bin together was a rounding
error. The physics and the unit conversions were correct throughout; six hours
of suspicion pointed at them anyway.

**The reason this is a class and not an incident is who else reads it.**
`carve_verdict.py`, `surface_water.py` and `export_carve_list.py` all take
`np.asarray(ds["spd"][:]).mean(axis=0)[-1]` for the Penman wind. That average
includes the bad bin and is inflated **1.55x**, which inflates open-water
evaporation, which biases the carve verdict toward less overflow and therefore
toward UNDER-carving. The lake solution already committed on this build carries
it too.

**What made it invisible.** The defect is in one of twelve bins of one field of
a shared product, and every consumer reduces that field to a single time mean
before using it. A 55% error in a mean looks like a climate, not a bug. Nothing
in the chain compares two things that should agree, so nothing had a chance to
notice.

**The rule.** A field that is reduced over an axis should be checked ALONG that
axis at least once, at the point where it enters the project. A per-bin global
mean printed beside its own median would have shown this the day the climatology
was built. `aeolian/scripts/build_dust.py:flag_anomalous_bins` does that now for
wind, and it is one function that any consumer can call.

The general form is the one this file keeps arriving at from different
directions: a derived product carries no evidence of its own health, and the
consumer that reduces it destroys the evidence before looking.

## 15. Two errors that nearly cancel, and the correction that would break it

A follow-on to class 14, and the more useful half. Having found that bin 0 of
the binned climatology carries corrupted winds, the obvious fix is to drop it
from the Penman wind the carve verdict uses. **That fix would have made the
number worse**, and only measuring the whole chain showed why.

Bottom level, global area-weighted, on the baseline climatology:

| quantity | m/s |
| --- | ---: |
| binned `spd`, all 12 bins, what the consumers used | 8.411 |
| binned `spd`, excluding the corrupt bin 0 | 5.033 |
| `sqrt(ua^2+va^2)` from the binned components | 3.922 |
| snapshot `spd`, 32 instantaneous samples an orbit | 7.356 |

Dropping the bad bin turns a 14% error into a 32% one, in the other direction.

**The mechanism is in the model's output accumulation, and I got it wrong twice
before getting it right.** In the snapshots `spd` is exactly `sqrt(ua^2+va^2)`,
ratio 1.000 and correlation 1.0000 at every level on every orbit, so the field is
what it claims and the snapshot value is the mean of instantaneous speeds. The
binned `spd` is not the same average of the same thing: it lands between the
speed of the time-mean vector and the mean of the speeds, because the records
pyburn averages have already been accumulated over an output interval by the
model, and a vector accumulated over an interval loses whatever reverses inside
it. With a 30-hour day and a ~1-day output interval, that is most of a diurnal
cycle.

**Locating an averaging error means finding the earliest artifact that has it.**
My second explanation put the vector averaging in `build_climatology.py`, on the
grounds that it averages five orbits and `ua`, `va` combine as vectors there
while `spd` combines as a scalar. That is true and it is not the cause: the gap
is already 5.002 against 7.342 within a SINGLE orbit's pyburn output, which
`build_climatology.py` has not touched. One measurement on the upstream artifact
would have shown that, and I reasoned about the code instead.

**The wrong mechanism, recorded because it was nearly paved into a component.**
I concluded instead that `ua` and `va` on sigma levels carried a cos(phi) factor
and were not physical winds, from `RevCosPhi` appearing in `burn7.cpp`'s
pressure-level path at line 4542, supported by a ratio that grew with latitude.
It is refuted pointwise on the raw artifact: on `MOST_SNAP.00086` the largest
difference between `spd` and `sqrt(ua^2 + va^2)` anywhere in the field is
**1.9e-06 m/s**, and the ratio is 1.0000 at 82 degrees where cos(phi) is 0.134.
A cos(phi) factor cannot hide at 82 degrees.

The latitude dependence that convinced me is real and has the other cause:
wind direction varies more at high latitude, so vector averaging cancels more
there, which mimics a 1/cos signature closely enough to pass a casual test. My
own test returned 0.88 to 1.08 in mid-latitudes and 0.41 to 0.62 at the poles,
and I read the polar failure as the approximation breaking down rather than as
the refutation it was.

**A turbulent flux wants the mean of the speed**, because the aerodynamic term
is essentially linear in wind speed -- and the lasting lesson is that I then
mistook where that quantity had gone. The binned field was NOT structurally
incapable of being a mean of the speed; it was corrupted into something else by
the same output path as the first-bin defect. With `NLOWIO = 0` the binned wind
is the mean of the speed to 0.2%, the snapshot workaround is unnecessary, and the
1.55x correction written against it would be a 55% error. Both were removed
rather than kept.

The value in use was high, and it was high by an amount that depended on the run,
because a large positive error from the corrupt bin and a large negative one from
the interval accumulation partly cancelled. On the bootstrap climatology the
residual was 3.1% and the cancellation looked almost complete, which is what made
the whole thing look ignorable; on the baseline it is 14%. **Do not carry a
cancellation from one run to another.** Neither error is a property of the world,
so nothing constrains them to keep the same ratio -- and neither, it turned out,
was a property of the diagnostic. Both were the output path.

**Two rules, and the second is the one that nearly did damage.**

When a bug is found in one term of a chain, measure the chain before correcting
the term. A compensating error is not a reason to leave a bug in place, but it
decides what the correction has to be: here the fix is to read the snapshot
product, which is right for the right reason, rather than to patch the binned
one, which would be wrong in a new direction.

**And a mechanism that explains the right number can still be the wrong
mechanism.** The cos(phi) story predicted the observed 1.5x, was consistent with
a real line of source code, and was wrong. Agreement with the number you are
trying to explain is the weakest possible evidence for a mechanism, because the
number is what you selected the mechanism to reproduce. The test that settles it
is the one the mechanism makes a DIFFERENT prediction for -- here, the poles,
where cos(phi) is 0.134 and the story demanded a ratio of 7.5. That test was
available, I ran a version of it, and I explained away the disagreement instead
of accepting it.

This project's standing rule covers it and I did not apply it: check the claim
against the artifact, not against the source that suggested it.

## 16. Tuning physics to a metric, and calling it validation

Found in my own hands, and corrected by the person reading over my
shoulder rather than by anything in this repo.

A stability correction was added to the Penman transfer coefficient, because a
lake's surface layer is stratified and stratification changes turbulent
exchange. It moved the comparison against the model's own ocean evaporation from
1.0415 to 1.0845. **I removed the correction on that basis.** That was wrong, and
the reasoning is worth naming because it is seductive.

**Physics is not a knob.** A process is in the model because it exists, not
because including it improves a number. Adding and removing terms according to
whether a comparison tightens is fitting, and fitting to a metric that is not
truth is worse than not fitting at all -- the model's ocean scheme is a different
parameterisation, not the right answer, so agreement with it was never the
target.

**"Correct physics made the comparison worse" is information, not a verdict.** It
means one of three things and none of them is "remove the physics":

- the implementation is wrong, which it was: `fluxmod.f90:248` branches on
  `dls < 1` to a Miller et al. (1992) free-convection form over WATER, and the
  first attempt used the land branch on a lake;
- the comparison measures something else, which it also did: Penman is a
  combination equation and the model runs a bulk formula off a prognostic surface
  temperature, and those differ by several percent whatever else is true;
- or a second error was cancelling the first, which is class 15 wearing a
  different hat.

**What a real check looks like.** The useful test was not the ratio. It was
comparing the surface temperature my closure implies against the model's actual
sea surface temperature, where the model HAS one: median error +0.24 K, and the
model's ocean is genuinely unstable at +1.51 K, confirming the branch. That test
could have failed and it did not. The ratio could only ever have told me that two
formulations differ.

Ask what would falsify the implementation, not what would flatter it.

## 17. A check that cannot fail is not a check

The consequence of class 16, and general enough to stand alone. Three things get
called validation here and only one of them is.

**Cannot fail at all.** `penman_ocean_validation` returned 3.736 against 3.672
for a ratio of 1.017 as literal constants, whatever the inputs were. It was
quoted for months, including as the check that a dust perturbation had not broken
Penman -- a thing it was structurally incapable of detecting. The computed figure
is 1.0415, so every "validated to within 1.7%" in this repo was also wrong by a
factor of two and a half.

**Can only differ.** Penman against the model's own bulk evaporation scheme. Two
legitimate formulations of the same quantity, disagreeing by several percent
because they are different formulations, with neither being the right answer. A
change in that disagreement is not evidence about either of them, and treating it
as evidence is how a correct physical term came to be deleted.

**Can fail.** A comparison with a right answer, where one outcome would mean
"wrong". All of these are from a single day and every one of them could have gone
the other way:

| check | why it can fail | what it caught |
| --- | --- | --- |
| term 7 against `-dshfl` | algebraically the same expression | out by 1.6 W/m2: the unaccumulated snapshot |
| snapshot `spd` against `sqrt(ua^2+va^2)` | definitional | agreed to 1.000: settled what code 259 is |
| `planck()` integrated against sigma T^4 | physical constraint | off by pi: a double-counted 2*pi*h*c^2 |
| implied surface temperature against the model's SST | the model HAS the answer | +0.24 K: the closure holds |
| binned `spd` against the snapshot mean at NLOWIO = 0 | same quantity, same regime | 1.002: the vector cancellation was an artifact, not a property |
| land P - E against the model's global water budget | conservation | 0.03%: the fields are sound and `mrro` was misread |

**Test the implementation against something that can fail, not the outcome
against something that can only differ.**

To find one, look for an identity, a definition, a conservation law, or a
quantity the other side already knows. Those are the four that work. If you
cannot say in advance what result would mean the thing is wrong, you have not
built a test -- you have built a number that will later be quoted as though you
had.

## 18. A transfer that moves the wrong element count, while its counter survives

This is the mechanism behind class 14 rather than another
instance of it. The evidence is in `exoplasim/notes/first-output-bin.md`.

PlaSim saves six SPECTRAL output accumulators through `mpputgp`/`mpgetgp`, which
declare their dummy argument `pp(NHOR,klev)`. Fortran reinterprets the
`(NESP,NLEV)` array against that shape and moves `NHOR*NLEV` elements, so at
T42/L10/p16 5,120 of 19,040 elements cross the restart and everything below model
level three comes back as zero. The counter that divides them, `naccuout`, is a
scalar and survives INTACT.

That combination is what makes it a class of its own:

- **The counter being right is what hides it.** Class 13's lesson was that state
  living in the checkpoint drifts from state living in the inputs, so the counter
  is where suspicion goes first. Here the counter is correct and the DATA is
  short, which is the same symptom with the opposite cause. This project flagged
  `naccuout` as the obvious candidate and was wrong to.
- **The vertical structure is the tell.** A counter error scales a field
  uniformly. A truncated transfer damages the elements past the cut and leaves
  the rest exact, so the error has a boundary in it -- here, clean above model
  level three and wrong below. Whenever a defect has structure a scalar cannot
  produce, the scalar is not the cause.
- **No language error is available to catch it.** The call compiles, runs and
  reports nothing; only the shapes disagree, and Fortran is content to
  reinterpret them. Nothing short of reading the routine's declaration finds it.

The general form: **when a quantity is restored across a boundary together with
the count of what it accumulates, the two can disagree, and the one that looks
authoritative is the one that survived unharmed.** Check the transferred SIZE
against the declared shape at every such call, not the value at the far end.

What made it provable rather than plausible was that the deficit is a pure ratio
and predicts numbers in advance: the surviving coefficients must show no anomaly
and the truncated ones must all show the same one, a different run length must
give a different constant, and the wind error must be a level-independent
`cos(phi)` solid body whose amplitude follows from `plavor`, the rotation rate
and the radius with nothing fitted. All four were stated first and then measured.
Class 17 is what turned a diagnosis into a mechanism.

## 19. A defect offered as a decision

Found in this project's own working method rather than in its code,
and recorded here because it wastes the scarcest thing in the loop: the user's
attention on the questions that actually need it.

Two items were put up as decisions with options and trade-offs. Neither was a
decision. `config/planet.yaml` declares a K2.5V star, and the Lacis and Hansen
shortwave absorptances are stated as fractions of SOLAR flux, so leaving them
unweighted is running a scheme built for a different star; "physics is not a
knob" settles that already. The config also NAMES the stellar spectrum, so the
model failing to read it was a defect against the config rather than a choice
about the world. In both cases the menu was assembled from what the CODE could
be made to say, not from what the project had already declared to be true.

**The tell: one of the options was "leave the known-wrong thing as it is".**
Where that appears, the question has been mis-framed. A defect has a fix, and
offering it as an alternative launders a bug into a preference.

The damage is not only wasted attention. A menu implies the project has no
position, so it quietly discards the position the project actually holds -- and
those positions, in `CLAUDE.md` and the pipeline chapters, are the accumulated result of
having been wrong before.

The distinction that works, applied in the same conversation and correctly that
time: a question is a DECISION when the declared truths CONFLICT, or when it
needs a threshold or preference nothing has fixed. Replacing the convergence
criterion qualified, because the threshold convention's ban on retrofitting a criterion
(`docs/src/practice/conventions.md`) and
conservation's verdict that the criterion measured the wrong quantity pointed
opposite ways, and the new threshold had to be argued from the estimator rather
than looked up. Everything else that day was work wearing a decision's clothes.

## 20. Reconciling what an upstream change has already made worthless

From the same root as class 19: treating derived
state as though it were an asset.

A radiation change landed, which means every artifact below the climatology
describes a world that no longer exists -- the verdict, the carve list, the dust
chain, the soil, the derived surface classes, the prospectivity, the budget's
basin currency. The instinct was to ask which of them needed updating, which
pairs were now mismatched, and whether a comparison had been contaminated by a
rebuild. All of that is wasted: the answer to every one of those questions is
that the artifact is worthless and will be regenerated.

Three shapes this took, all in one day:

- **Staleness tracked as work.** A task row saying "rebuild these fields before
  the next run" is state, and state is what `check_consistency.py` reports. It
  was closed as not-a-task.
- **Numbers archived as findings.** A basin count and an optical depth were
  written into the permanent record beside the mechanisms that explain them.
  Only the mechanisms survive an iteration.
- **Cost attributed to regeneration.** Deferring a rebuild to avoid pairing a
  corrected field with an uncorrected climatology, when both were disposable and
  the correct move was to regenerate both.

`CLAUDE.md` rule 7 is the positive statement. The negative one is here: after
an upstream change, do not ask what needs updating. Ask what is now worthless,
and then stop thinking about it.

## 21. A wait condition that matches itself

Trivial mechanically, and expensive out of all proportion
to its content, because the failure mode is silence: nothing errors, nothing
progresses, and nothing prompts a look.

Waiting for a build with

    until ! pgrep -f build_star_cycle; do sleep 20; done

never exits. The loop's own command line contains the pattern, so `pgrep` finds
the waiter and reports the build as running forever. Worse, the command that was
supposed to START that build opened with `until ! pgrep -f continue_exoplasim`
and self-matched the same way, so the build never began -- and five successive
waiters then reported a non-existent process as busy for the best part of an
hour, each keeping the pattern alive for the others.

This is class 17 in different clothes: **a check that cannot fail.** The
condition was not testing whether the build was running; it was testing whether
the test existed. Nothing it could observe would have ended the wait.

Two rules, and the second is the durable one:

- Never `pgrep` for a string the polling command itself contains. If a process
  test is unavoidable, match on the binary or use `pgrep -f -- "$pat"` with the
  pattern built so the waiter cannot contain it.
- **The disguised variant is the one to watch for**, because it reads as sound:

      until [ ! -e /proc/$(pgrep -f build_star_cycle.sh | head -1) ]; do ...

  Testing whether a PID's `/proc` entry still exists IS a good way to wait on a
  process. The rot is one level down, in the `pgrep` that supplies the PID: it
  returns the waiter's own, whose `/proc` entry necessarily exists. Three of
  these were found running at 24, 35 and 61 minutes, and this was the oldest --
  the two obvious ones were spotted first precisely because they looked wrong.
- **Wait on the ARTIFACT, not the process.** Every long step here writes
  something: `rebuild_binaries.py` writes `binary_manifest.json`, the cycle build
  writes `cycle_binary_manifest.json`, a segment writes `MOST.NNNNN.nc`. Waiting
  for the file to appear or its mtime to move is immune to this by construction,
  and it tests the thing actually wanted -- the output -- rather than a proxy for
  it. The binary's mtime was visible throughout and would have settled the
  question in one command.


## 22. Two copies of what a segment must re-apply, and only one of them maintained

`configure()` re-copies the shipped namelists over the configured ones every time
the model is constructed on an existing run directory. So a continuation has to
REAPPLY every configured key, and `continue_exoplasim.py` did -- from its own
copy of the list.

`run_exoplasim.py` staged five shortwave gas keys. `continue_exoplasim.py`
restated the tuple with four, under a comment reading "all four", which was true
on the day it was written. PHYS-9 then added `h2o_sw_level` to one copy. The
result: `config/planet.yaml`'s 1.127 applied to the orbits the prepare script
ran and reverted to radmod.f90's 1.0 on every segment after.

Found on `run_78c22fb1a1bd`. Orbits 0-6 carried the level and 7-59 did
not, so the physics changed PARTWAY THROUGH A RUN and the convergence window and
climatology sat entirely in the wrong half. The bundle prices that term at
+1.25 K, its largest. A second instance of the same shape was found in the same
file within the hour: `configure()`'s `otherargs` was also duplicated, and
CLIM-17's `TFREEZE` was in the prepare copy only, so the run's `icemod_namelist`
carried `NICE` and nothing else.

**This is not class 12.** There a continuation RE-DECIDED a value it should have
read back from the run. Here it correctly re-applies -- from a list that had
quietly stopped being the list. The two point opposite ways and the fixes differ:
12 wants the value read from the run, 22 wants there to be only one statement of
what gets applied.

**Why nothing caught it.** `config_drift` compares the config against the
manifest's `source_config`, and both were right; the config said 1.127
throughout. `check_consistency` checks artifacts against the config, and the
namelist is not one of the artifacts it reads. `assess_convergence` reads the
trailing window, which was uniformly in the wrong regime and therefore perfectly
self-consistent. Every guard passed because each compared two things that agreed,
and the disagreement was between the config and a file neither of them opened.

**The tell, and it is cheap.** Class 12 already recorded that the namelist in the
run directory is the ground truth for what was integrated and the manifest is a
claim about it. Read the namelist. `H2OSWL` was absent from the exact slot
between `CO2SWW` and `H2OSWW` where the insertion order puts it, which is a
one-command check and was the whole of the evidence.

**The rule.** What a segment must re-apply is ONE definition, imported, never a
second list that has to agree. Two lists that agree today are a defect with a
date on it: the second copy is not maintained by whoever extends the first, and
nothing in the tooling connects them. The same bug had already hit the stellar
spectrum on this file once, and the fix then was to add the missing call rather
than to remove the duplication, which is why it recurred twice.

## 23. Settling a conceptual question by measuring it

A thing that should not exist is not refuted by measuring how badly it performs,
and measuring it is worse than saying nothing, because the measurement records
the wrong reason and the wrong reason expires.

**What happened.** `pipeline.py` read `TASKS.md` to print a "carve
gate": the open tasks whose `[step: <id>]` marker named a step upstream of
`carve_list`. Two fixes were attempted, in order, and both were answers to
questions nobody had asked.

The first moved the parsing into a new `scripts/carve_gate.py`, so that the
graph tool read no tracker. That relocated the coupling and left the machinery
standing, and the coupling was never the defect: the gate was DEFINED over
tracker prose, so any correct implementation of it must read `TASKS.md`. A
defect that cannot be removed without removing the feature it serves is a fact
about the feature, and this one was pointing at the gate the whole time. The
commit recorded "output is byte-identical" as though that were a virtue, which
is the second tell -- preserving behaviour settles nothing until the behaviour
is known to be wanted, and here it was the evidence that nothing had been
examined. The question of whether the gate should exist was never asked,
because the gate already existed.

The second deleted the script, and justified it by counting: 58% of the graph
was upstream of the carve, so the filter selected 10 of 11 open rows. The
deletion was right. **The justification was not, and it was the more damaging
half**, because it framed a category error as a performance problem. A
justification by ratio says the thing would be fine at a better ratio, so it
licenses the return of exactly what it removed, and the count that supports it
goes stale within an iteration -- see the first convention in `CLAUDE.md`.

**The actual argument, which does not contain a number.** The gate asks whether
any outstanding work still moves an artifact the verdict is computed from. That
is a fact about what closing a row would CHANGE, and it lives in the row's
prose. A `[step: <id>]` marker records where the work is FILED. A task can name
an upstream step and move nothing -- it is deferred, or it is a bound, or it is
blocked on something that happens after the carve. So the marker does not
determine the answer, at any count, and a traversal over markers cannot compute
the gate however well it appears to score.

**The tell.** You reach for a measurement to decide whether something should
exist. Before running it, ask what the number would have to be to change the
conclusion. If no value would -- if you would delete it at 10 of 11 and also at
3 of 11 -- then the question was never empirical, and the measurement is
decoration that will later be quoted as the reason.

**Related but distinct.** Class 16 is tuning physics to a metric. Class 17 is a
check that cannot fail. This is the inverse of 17: a check that CAN fail,
performed on a question where failure and success were both irrelevant.

## 24. A derived summary standing in for its primitives

A count, a partition, a pairing, or an exception claim asserted where the
primitive facts it summarizes are absent, wrong, or contradicted: "three of
the six pair naturally and one cannot" over pairings that shared a member;
"seven surface fields" over an enumeration that omitted the conditional ones;
"none of them announced itself" at the head of this very catalog, against
three of its own loud classes. The summary reads as rigor and cannot be
checked without the primitives, so it survives every review that does not
recompute it. The fix is mechanical: state what gates what, or what defines
what, and keep the derived claim only if it still checks out and changes a
decision. The instances are in git.

## 25. Byproduct kept as content

The documentation accumulates what writing it produced rather than what
reading it needs: the story of how a section came to be right, dates on
arguments that cannot go stale, notes about the note, ceremony standing where
a status label already carries the weight. Each reads as diligence, and each
is noise a future reader pays for on every load -- the cost is not the disk,
it is the context. The test is the future decision a passage could change; if
none, it goes, however true it is. Git already records the how; keeping it in
the document as well is the failure. The one-line cost history that
calibrates a rule, the measured-on date, the pre-registered threshold and the
status label are not this class: each changes how a reader acts.

## 26. A poison seed

A line that is true, brief, and harmless on its face, but that reliably
steers a reader who has it in context toward a documented antipattern -- the
attractor is the problem, not the content. The proven mechanisms: a
prohibition that prints the thing it forbids, which re-supplies the material
on every load (the vocabulary rule's first version did exactly this); a rule
system that announces its own tensions, which converts work into adjudication
of the rules instead of use of them; an exception stated more vividly than
its rule, which teaches the exception; a prediction that the reader will
fail, which reads as permission to; an untracked "worth checking someday",
which recruits sessions into work nothing asked for; and tool output whose
form outruns its content -- a formatted summary that reads as a verdict. The
test is not the sentence's truth but its expected effect over many loads.
The defusals are mechanical: state the positive form and put specifics behind
a pointer; replace announced tension with the decision procedure that
resolves it; fence an exception with the condition that licenses it; convert
an open loop to a task row or delete it; make tool output state what it is
and is not.

## 27. A comparison whose baseline was never established

Two builds were compared on a bed that starts cold with `KICK` non-zero and no
`SEED`. `plasim.f90` seeds that white noise from the CLOCK when `seed(1)` is
zero, so one binary run twice on that bed does not agree with itself, and every
comparison built on it is a comparison of two clock readings. It reported 94 of
199 records wrong at ONE step against a correct transform rewrite, in the loud
unambiguous form that invites belief, and the rewrite was most of the way to
being reverted before the bed was suspected.

The tell was in the output and was read past: the worst record was `seed`
itself, at 4.3e+170 relative. When the RNG state or another metadata record is
among the differences, the runs are not comparable at all and no finding about
the physics can be read out of them.

**Before comparing A with B, compare A with A.** If one binary run twice does
not give one answer, the comparison has no resolution and its verdict is noise,
whichever way it comes out. The self arm is cheap, it runs first, and its
failure is diagnostic rather than confusing. `verify_shared_determinism.sh`
runs it as its first arm for this reason, and `_bed_guard.sh` refuses the
specific bed that caused this, which is the mechanical half of the same lesson.

The general form is wider than beds and wider than this model. A check has
inputs it does not control -- an unfixed seed, a wall-clock, a directory it
shares, an environment variable set outside it -- and any of them can supply
the whole of the difference it reports. The distinguishing question is not "is
the difference large" but "does the null case come out null".

Related: class 17 is a check that CANNOT fail. This is its mirror, a check that
fails for a reason outside what it checks, and it is the more expensive of the
two because it produces work rather than merely permitting it.

## 28. One component measured against the whole's budget

A constraint is stated over a total -- a working set, an error budget, a mass
balance -- and then a change to one term is checked against the total's ceiling
as though the term were the whole. It reads as rigor because a real number is
compared with a real limit, and both numbers are right; only the comparison is
meaningless.

Three times in one afternoon, on one change. The Legendre weight factorisation
was reported as bringing the per-die working set "to about 36 MB, under CCD1's
32 MB" -- which is not under it, and was not caught until someone read the two
numbers side by side. Corrected to 30.82 MB by sharing what had been duplicated,
it was then reported as fitting the die, when 30.82 MB is the WEIGHTS and the
die also has to hold the Fourier fields, the spectral state and the reduction
partials. The total is about 260 MB. The weights had gone from four times the
die to about one, which is a large win on the dominant term and is not a fit.

The tell is a sentence that names one array and one cache in the same breath.
The fix is to make the total computable and cite it -- `cache_budget.py` for
this one -- so that a component's number cannot be quoted as the budget without
the budget being visible beside it.

It is not class 24. That is a summary asserted without its primitives; this is a
primitive asserted in place of the summary, and it fails in the opposite
direction: the part is measured correctly and compared to the wrong thing.

## 29. A check that simplifies away the configuration it certifies

A verification driver stands the component up outside the model, which is the
whole reason it is sharp: it can drive dense fields through both sides, isolate
one mode, and carry controls the model could never run. Standing it up means
choosing values for everything the model would have supplied, and each of those
choices is a chance to certify a configuration that does not exist.

`verify_shtns_equivalence.f90` reported nine arms at 5e-14 while the model it
certified disagreed with itself by 100% at the first step. It set `plavor = 0`
so the planetary vorticity correction would not enter the comparison, and it
never read a namelist, so `nfilter` kept its default of none. Those are exactly
the two things the SHTns wrappers were missing: legmod takes the planetary
vorticity back out of the wind, and `legini` folds `skspgp(n+1)` into `fsp`,
`fmu` and `fmv`, so every conversion it performs is filtered. Both terms are
identically one in the configuration the check chose, and neither is one in any
configuration the model runs.

The simplification is not the error. Zeroing a term to isolate another is
correct practice, and a driver that had SAID it certifies the unfiltered
non-rotating case would have been honest and useful. The error is the scope
claimed: nine green arms were read as "the wrappers compute this model".

The tell is a default taken rather than set. A term the driver never mentions
is a term whose value came from a module initialiser, and a module initialiser
is not the model's configuration. Grep the driver for every quantity the
operator under test depends on and ask which of them the model sets from a
namelist.

Two fixes, and the second is the one that generalises. The driver now sets the
beds' `nfilter` and a rotating planet, and a control drops `fsp` from every
wrapper and must fail. But an array comparison certifies the TRANSFORM, and the
model reaches that transform through a namelist, a build configuration and a
call site the array comparison never touches -- so the gate that closes it is a
second check that runs the MODEL, `verify_shtns_model.sh`, one binary with the
switch flipped. A component check and a model check answer different questions
and neither substitutes for the other.

It is not class 17. There the check has no right answer; here it has one and
computes it exactly, for a configuration nobody runs.

## 30. A stochastic property tested by a pair

A check asks whether two runs agree, they do, and the answer is recorded as
"reproducible". But agreement is the question, not the observation: if the
property fails only sometimes, a sample of two is a coin landing heads twice
and reports the same thing whether the underlying rate is one in ten or never.

`verify_shtns_model.sh` asked for bit identity twice and got it, more than once,
across several days of changes -- while the model actually had three or four
distinct outcomes, in clusters of two to four identical runs. The clustering is
what makes this worse than an ordinary sample-size error: consecutive runs are
positively correlated, because whatever varies between them (machine timing,
which algorithm a self-tuning library picked at startup) persists for a while.
So a pair is not even two independent draws.

The tell is a check whose subject is a property of a DISTRIBUTION -- run to run
identity, convergence, a race not firing -- tested with the smallest sample that
can technically show disagreement. Ask what rate the check could miss: two runs
cannot see one-in-three, and this one did not.

There is no universal number. The bar is that the sample be large enough to
catch a rate that would matter, stated where the check declares its other
thresholds, and larger where consecutive draws are correlated. Four runs here,
which catches one-in-three about four times in five; the cost is four T21 runs.

It is not class 17. There the check has no answer that could mean "wrong"; here
it has one and asks too few times to see it.
