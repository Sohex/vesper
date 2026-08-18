# How this project goes wrong

Every bug in this list was found late, by comparing two things that were supposed
to agree. None announced itself: array shapes matched, fields looked plausible,
runs converged, and the answer was wrong. They are recorded by *class* rather
than as a changelog, because the classes recur and the instances do not.

Run `python scripts/check_consistency.py` before an expensive run and after
changing `source_build`. It mechanises the checks that come out of the classes
below, and it exists because prose in a document does not stop anything.

---

## 1. One quantity, several consumers, a fix that reaches some of them

The most expensive class here, twice over.

**`mrro`.** Four scripts integrated catchment runoff. `pedology` established that
`mrro` is river-routed net divergence rather than local generation and switched
to P − E, documenting why. `surface_water.py` reached the same conclusion
independently and also switched. `export_carve_list.py` and `carve_verdict.py`
were never updated, and the first of those is the one whose output leaves the
project and changes the terrain. Two independent corrections of one
misunderstanding, neither of which propagated to the consumer that mattered most.
Cost: an entire terrain build.

**`field_lon`.** The longitude-convention fix landed in `basin_means` and two of
its three callers. The third was `export_carve_list.py`. Again the one that
changes the terrain.

**The rule.** When a shared quantity is redefined, enumerate its consumers before
declaring the fix done — `grep` for the symbol, not for the bug. A fix that lands
where the bug was noticed is half a fix. Prefer one canonical accessor over four
call sites that each decide for themselves.

## 2. An optional argument whose absence means "do the wrong thing"

`basin_means(..., field_lon=None)` defaulted to the unremapped path. So the fix
for the antipodal bug reproduced the antipodal bug in any caller that had not
been updated — the original error wearing the shape of its own correction.

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
switches and a digest of every surface file — but not the stellar spectrum, and
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

- The energy residual was recorded as "a fixed offset, −0.455 ± 0.010, which
  rules out every state-dependent candidate". Measured the same way across four
  converged runs it spans −0.400 to −0.489, nine times wider, putting the
  ruled-out candidates back in play.
- `lake-representation.md` prescribed a large `dwmax` because "a large full
  bucket evaporates at open water's rate". Wetness reaches 1 above 40% *of*
  `dwmax`, so a deeper bucket needs proportionally more water; and routed river
  water never re-enters the evaporating bucket at all, so no setting sustains a
  lake.
- A stellar-spectrum correction was estimated at 0.4 to 0.7 K and measured at
  0.04 W/m², because the estimate was computed on a world with far more snow.

**The rule.** A recorded conclusion is evidence about what was true of the cases
it was measured on. Before relying on one, check whether the current
configuration is inside that envelope. Re-measure rather than inherit.

## 7. Unguarded arithmetic at a physical boundary

Switching catchment runoff to P − E made it negative for 57% of basins, which is
physically meaningful — those catchments evaporate more than they receive — but
the lake-area solve `A = R·C/(E − P + R)` was not written for it and reported a
lake covering −6.7 × 10¹⁵ percent of the planet. Earlier, a spill-level fallback
returned `inf` and produced an infinite basin capacity.

**The rule.** Clamp at the physical bound where the quantity is defined, not
where it is used. A catchment delivers zero or more, never less.

---

## The pattern behind the pattern

Almost every entry here was found by comparing two artifacts that should have
agreed, and almost none by reading code. The productive habit is to look for
quantities computed two ways and check them against each other: mesh against
grid, our routing against the exporter's, Penman against the model over ocean
cells, one build's composition against another's. Where a cross-check exists,
these bugs surface in minutes. Where none exists, they survive until something
downstream looks strange.

So when adding a component, add the cross-check with it, and prefer the check
that would have caught the last bug.

## 8. Reading an intermediate configuration as the world

Three conclusions in this project have been right in mechanism and wrong in
magnitude, and the last two were wrong the same way: a figure was taken from the
build sitting in front of me rather than from the configuration the pipeline is
converging on.

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
every cycle, which is exactly why they keep getting picked up.

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
kills the gridcell **silently** -- zero LAI, zero AET, no crash. And because
water is fractional, thinner soil reads as relatively *wetter* for the same
absolute water, so a single cell can move either way through PFT competition.
Trust a controlled sweep, never one cell.

## One quantity, two meanings, three times the value

"The endorheic share of land" names two different measurements that differ by a
factor of about 3.5:

- **12.4% of land** is inside a preserved basin (`is_endorheic` cell area). This
  is the lithology and thermostat figure.
- **43.0% of land** *drains* to a closed basin. This is the hydrography figure
  and the one the carve verdict is about.

Both are correct. A basin's catchment is far larger than its floor, which is the
whole reason a small area of fill can decouple a large share of weathering. They
are named apart in `world_state.json`; quote the name, never "the endorheic
share".

## 11. A patched source, and one binary per configuration

ExoPlaSim compiles a separate executable for every (resolution, layers, ranks)
triple: `most_plasim_t42_l10_p16.x` and `most_plasim_t21_l10_p8.x` are different
files built from the same source at different times. **Patching the source and
rebuilding rebuilds only the configuration you are running.** Every other binary
on disk keeps the old code, silently, until something asks for it.

Found 2026-08-16, moving from the T21 flux bracket to a T42 bootstrap. The ozone
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
digest, the spectrum and the flux into `run_id`, and the binary is the input that
is still missing from it.

## 12. A continuation that re-derives the physics from config

`run_exoplasim.py --flux-ratio 0.91` prepares a run at 0.91 and stamps it in the
manifest. `continue_exoplasim.py` then took its flux from
`config/planet.yaml:orbit.baseline_flux_earth`, so continuing that run integrated
it at the baseline instead, and the manifest went on saying 0.91.

Found 2026-08-17, on the first flux bracket of the baseline re-run. Nothing failed:
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

**What it cost, and what to check.** Forty-five orbits, about two hours, and a
run whose manifest had to be corrected rather than trusted -- the correction is
recorded in the manifest itself, against the namelist on disk, which is the only
artifact that recorded what actually ran. When a segment finishes, the namelist
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

Hit again 2026-08-17, seeding the baseline from the bootstrap to skip 86 orbits
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

Found 2026-08-17 while chasing a dust emission flux that came out 1700x Earth's.
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
