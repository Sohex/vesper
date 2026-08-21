# The CH4 and N2O band: the design, and the tests declared before the work

Worldbuilding. Vesper is an invented planet and this note is about the
simulation of it: a toy climate model's longwave scheme and a term being added
to it. Every quantity named here is a modelled field or a laboratory
measurement of a gas.

CLIM-42. `exoplasim/notes/trace-gas-absorbers.md` prices what the missing term
is worth -- 1.56 to 1.98 W/m2, the largest single item in
`analysis/error_budget.json` -- and CLIM-43 settled the mixing ratios it has to
carry. This note is the design of the term itself, written before the code, on
CLIM-39's precedent: the tests are declared here so the implementation cannot
be judged against a bar chosen after it runs.

**Nothing here is a case for tuning.** The gas is in this world's atmosphere at
a modelled abundance that a photochemical model computed from fixed biogenic
fluxes and this star's ultraviolet. That is what justifies the band. A term
added because it improves an agreement is
`docs/src/practice/failure-modes.md` class 16, and the bound in the pricing
note says the term is LARGE, not that any particular value of it is right.

## 1. Why the existing scheme has no key to set

`lwr` in `radmod.f90` is Sasamori (1968): a broadband scheme in which each
absorber contributes an ABSORPTIVITY as a function of its own accumulated
amount, and the clear-sky layer transmissivity is one minus their sum,

    ztaucs = 1 - a(h2o) - a(o3) - a(co2) * t(h2o)

with `t(h2o)` the Boer et al. (1984) water vapour transmissivity in the CO2
overlap region. Sasamori fits H2O, CO2 and O3 and nothing else, so there is no
coefficient for CH4 or N2O to take -- which is what `config/planet.yaml` means
by "adding either means adding a band, not a key".

## 2. The band model, and why this one

Donner and Ramanathan (1980), *Methane and Nitrous Oxide: Their Effects on the
Terrestrial Climate*, J. Atmos. Sci. 37(1), 119-124. It is the paper that did
exactly this job: simplified band models for the CH4 and N2O longwave bands,
built to be dropped into a broadband radiation scheme.

It adopts the Cess and Ramanathan (1972) band model as modified by Ramanathan
(1976), which gives the total band absorptance `A` in cm-1 as

    A(U, beta) = 2 A0 ln[ 1 + U / sqrt(4 + U (1 + 1/beta)) ]     (1)
    U    = S W / A0                                              (2)
    beta = beta0 (P / P0)                                        (3)

`W` is the absorber amount in cm atm, `S` the band intensity in
cm-1 (cm atm)-1, `P` the broadening pressure against `P0` = 1 atm, `A0` the
bandwidth parameter and `beta0` the line shape parameter at `P0`. Table 1's
parameters, with their temperature scaling:

| | CH4 1306 cm-1 | N2O 1285 cm-1 | N2O 589 cm-1 |
| --- | --- | --- | --- |
| `A0`, cm-1 | 52 (T/300)^1/2 | 20.4 (T/300)^1/2 | 23 (T/300)^1/2 |
| `beta0` | 0.17 (300/T)^1/2 | 1.12 (300/T)^1/2 | 1.08 (300/T)^1/2 |

**Three reasons this is the right source rather than a convenient one.**

It is a BAND ABSORPTANCE, which is the quantity Sasamori's scheme is written
in. The alternative that was to hand -- Byrne and Goldblatt's forcing fits,
already used to price the omission -- is a top-of-atmosphere forcing for one
atmosphere, and putting a global-mean forcing inside a layer-by-layer solver
would be a category error however well it reproduces the global mean.

It carries the OVERLAP the same way the existing code does. Donner and
Ramanathan handle the water vapour overlap in the 1306 and 1285 cm-1 regions
"by multiplying the band absorptances of CH4 and N2O with the mean
transmissivity of water vapor in the region of overlap", which is line for line
what `lwr` already does to CO2 with `zth2o`. The 589 cm-1 N2O band overlaps CO2
as well, and that is a term the scheme has the pieces for.

It has PRESSURE and TEMPERATURE dependence built in, at `beta = beta0 P/P0`.
That matters here more than on Earth: this world runs 1 bar over 12.81 m/s2,
so its pressure at a given absorber amount is not Earth's, and a
parameterisation with the broadening folded into a constant would carry Earth's
column silently. `npbroaden` already switches the same physics for the existing
absorbers.

## 3. The band intensities, and how CH4's was recovered

Donner and Ramanathan give `A0` and `beta0` in Table 1 but take the band
intensities from elsewhere -- Cess and Chen for CH4, McClatchey et al. (1973)
for N2O -- and do not restate them as numbers.

**For CH4 the paper contains its own answer.** Table 2 lists eleven
(total pressure, absorber amount, Eq. (1) absorptance) triples for the 1306
cm-1 band. With `A0` and `beta0` fixed by Table 1, Eq. (1) has exactly one free
parameter, so `S` is determined by that table rather than chosen:

    S(CH4, 1306 cm-1) = 187.7 cm-1 (cm atm)-1

which reproduces ten of the eleven rows to a maximum of 0.23 cm-1, 1.6%
relative, rms 0.107. The eleventh row, `P` = 1.0 and `W` = 0.505, misses by
2.4 cm-1 in a table whose neighbours fit to a fifth of a wavenumber, and its
Eq. (1) entry of 41.6 is the same number printed diagonally below it in the
adjacent column. It is read as a typesetting slip in the 1980 table and is
excluded, WITH the fit reported both ways: including it gives S = 194.3 and an
rms of 0.83, so the choice moves S by 3.5% and nothing else.

**For N2O both intensities come from McClatchey et al. (1973), which is the
source the paper names.** Its Table 13 gives band SYSTEM intensities, and a
system is the right quantity: Donner and Ramanathan's 1285 cm-1 analysis
"includes the fundamental and the first hot band and, furthermore, accounts for
contribution from four isotopes", which is what a system sums. Converted by
Loschmidt's number, since McClatchey quotes per molecule cm-2 and `W` is cm atm
at STP:

| band | McClatchey Table 13 | `S`, cm-1 (cm atm)-1 |
| --- | --- | ---: |
| N2O 1285 cm-1 | 996 +- 40, times 1e-20 | 267.6 |
| N2O 589 cm-1 | 118 +- 9, times 1e-20 | 31.7 |

**The three parameters are a MATCHED TRIPLE, and that is what decides the
sourcing.** `A0` and `beta0` are not independent measurements: Donner and
Ramanathan obtained them by fitting Eq. (1), at a particular `S`, to laboratory
absorptance. Substituting a different compilation's `S` into their `A0` and
`beta0` breaks the fit rather than modernising it. So CH4's comes from their own
Table 2 and N2O's from the compilation they cite, and neither is swapped for the
other's.

That is visible in the one band both sources carry. McClatchey's Table 18 puts
CH4's 1306 cm-1 intensity at 5.87e-18 per molecule cm-2, which converts to 158
against the 188 that Donner's Table 2 requires -- a 19% gap between two
published compilations of one band. **It is used as a check on the conversion
and not as a value.** Any unit error in that arithmetic would be a factor of
1e19, 100 or 10; landing at 1.19 rules one out, and it does not validate the
number, which is why the CH4 intensity still comes from Table 2.

### The route that failed, recorded so it is not retried

Fig. 2 plots Eq. (1) for the 1285 cm-1 band at two pressures, which is two
independent constraints on one unknown, so an extraction from it self-checks.
It failed: tracing curve A by continuity off a 600 dpi render and fitting over
the same absorber range gives S = 302 from the 0.5 atm panel against 348 from
the 0.1 atm panel, each fitting its own trace to an rms of 0.2 cm-1. A 15%
systematic gap between panels that individually fit that tightly is a
calibration error in the extraction, and a number from either panel alone would
have carried it invisibly.

McClatchey's 267.6 sits 8% below the 0.5 atm trace and 15% below the 0.1 atm
one, so it is bracketed by the two panels' own disagreement rather than
contradicted by it. That is corroboration at the level the trace can support,
which is weak, and it is reported as such.

Ramanathan (1976) was also fetched, on the chance it tabulated the N2O
parameters, and does not: it carries no N2O at all. Its row in
`references/INDEX.md` says so, so that is not re-checked.

## 4. What lands in the code

Per layer, alongside `zqco2` and by the same construction:

- absorber amounts `zqch4` and `zqn2o` in cm atm, scaled from the model's
  mixing ratios by the same `zsfac` path, with their own `zfch4`/`zfn2o` STP
  density factors;
- accumulated `zsumch4` and `zsumn2o` down the column, as `zsumco2` already is;
- Eq. (1) evaluated per band and divided by the broadband flux the band sits in
  to give a fractional absorptivity, which is the conversion that puts a cm-1
  band absorptance into Sasamori's currency;
- multiplied by the water vapour overlap transmissivity `zth2o` the CO2 term
  already computes, and for the 589 cm-1 N2O band by a CO2 overlap as well;
- subtracted from `ztaucs` alongside the existing three.

Mixing ratios come from `config/planet.yaml` through a namelist key, on the
pattern `dqco2` already sets, with the values CLIM-43 measured. They are an
ASSUMPTION of this world in the same sense 450 ppm of CO2 is, and the config
should say so where it declares them.

## 4b. The CO2 overlap at 589 cm-1, which is not optional

**The 589 cm-1 band is worth carrying and cannot be carried bare.** Weighted by
the Planck function at 255 K against the band absorptance at this world's
column, the three bands split about 45 / 26 / 29 percent, CH4 1306 / N2O 1285 /
N2O 589. The 589 band is the weakest of the three -- 5.0 cm-1 of absorptance
against 22.4 -- and it lands at 17 um where the Planck function is five times
larger, which is what makes it comparable.

It also sits inside CO2's 15 um band, 78 cm-1 from the 667 cm-1 fundamental. So
including it WITHOUT the CO2 overlap credits N2O with absorption CO2 already
provides: the error is an overstatement, not the conservative understatement
that omitting the whole band would be. Once the band is in, the overlap is part
of it.

**The procedure is Ramanathan (1976) Appendix A**, which Donner and Ramanathan
name. Its shape is the Goody (1964) relation `T = exp(-A_bar)` with
`A_bar = A / (2 A0)`, evaluated for CO2 in the 589 cm-1 region using an
EFFECTIVE intensity: Edwards and Menard's line intensity distribution shifts a
band centred at `w1` into the region of a band centred at `w2` as
`S_eff = S exp(-|w2 - w1| / A0)`. The absorptance itself is the same Eq. (1)
already being implemented, so the overlap reuses the band function rather than
adding a scheme.

**What that costs is a CO2 band model the host scheme does not have.** PlaSim's
`lwr` carries Sasamori's broadband CO2 absorptivity, not a band model, so
Ramanathan's ten-band 15 um treatment has to come in. It is contained: it
computes a transmissivity MULTIPLIER for one N2O band and never touches CO2's
own contribution to the flux, which stays Sasamori's. Two CO2 representations
in one routine is still worth saying out loud, and the containment is the
reason it is acceptable.

**The parameters are now assembled**, and where each comes from:

| quantity | source |
| --- | --- |
| CO2 15 um band strengths, B1 to B7 | Dickinson (1972) Table 3, already in cm-1 (cm atm STP)-1 |
| CO2 band centres, B1 to B10 | Ramanathan (1976) Table 3, after Goody (1964) |
| bandwidth parameter A0 | Cess and Ramanathan (1972) |
| effective-intensity shift | Edwards and Menard (1964) |
| hot and isotopic band summation | Edwards (1965), via Ramanathan Eq. (12) |

### Two attributions that do not hold up, checked against the page

Ramanathan (1976) p. 1333 says, verbatim: "The parameter `q_i` is the ratio of
the abundance of the individual isotopes to the total abundance of CO2. The
`q_i`'s are obtained from Goody (1964). The values of `D_i` and `S_i` are taken
from Dickinson (1972), `A_0` is taken from Cess and Ramanathan (1972), and
`v_0` = 0.064 cm-1 atm-1." Read off the page image, because the scan's OCR
cannot be trusted on subscripts.

**`S_i` is there; `D_i` is not.** Dickinson (1972) Table 3 tabulates band
STRENGTHS and nothing else. Dickinson's own treatment is statistical -- "mean
line strength of group" against "number of lines in group" -- so a mean spacing
is DERIVABLE from it but is not tabulated, and lifting `D_i` from the cited
source is not possible as the sentence implies.

**The internal cross-reference in the same sentence is also wrong.** It sends
the reader to Table 2 for the band list; Ramanathan's Table 2 is the
sensitivity of global surface temperature. The band list is Table 3, in
Appendix A, and it carries centres only.

Two loose references in one sentence is a reason to prefer the route that does
not depend on either.

### Which is why the correlated-k route becomes the primary

**It needs neither `D_i` nor `q_i`.** The LMD tables are HITRAN 2020 line by
line at natural isotopic abundance, so the line spacing and the isotopes are
already integrated into the k-distribution rather than being parameters someone
has to supply. The two quantities that cannot be cleanly sourced are exactly
the two it does not ask for.

So the CO2 transmissivity in the N2O 589 cm-1 region is taken from
`N2-CO2var_2026` at this world's CO2 amount, with the far-wing and CIA
companions applied, and Ramanathan's band sum becomes the CROSS-CHECK. Run over
B1 to B7 -- the fundamental and all six hot bands, which is the bulk of the 15
um absorption -- it needs only Dickinson's strengths, `A_0` from Cess and
Ramanathan, and a line spacing that has to be stated as an assumption rather
than cited. B8 to B10, the minor isotopic bands, need `q_i` and are left out of
the cross-check with their omission declared.

For `q_i`, if it is ever wanted: McClatchey (1973) is already on disk and is a
line parameters compilation, so its own isotopic abundances are the internally
consistent choice beside the N2O intensities taken from it -- better than a 1964
textbook for a number that has been remeasured since.

**A second, independent route is available and is a CHECK rather than a
substitute.** The LMD Generic PCM bundle on this host carries
`N2-CO2var_2026`, a HITRAN 2020 correlated-k table with CO2 as the variable
gas, whose IR band 28 spans 546 to 630 cm-1 against this band's `2 A0` of 46.
Its far-wing and CIA companions are present too --
`CO2-N2_line-far-wings_50-1000K_2026.dat` and the CO2-CO2 CIA -- which closes
the 25 cm-1 line cutoff those tables are built with. A band-mean CO2
transmissivity computed from it is a modern line-list answer to the same
question Appendix A asks, and the two disagreeing would be information about
a 1972 band model rather than a defect.

Note why substituting here would have been legitimate where substituting `S`
was not: the overlap is applied when the band is USED in an atmosphere, not
during the fit that produced `A0` and `beta0`, so it is separable in a way the
matched triple is not. The reason to do both anyway is that neither costs
anything now that the papers are on disk.

## 4c. What the cross-validation returned

Run 2026-08-20, `exoplasim/scripts/co2_overlap_589.py`, writing
`exoplasim/analysis/co2_overlap_589.json`. Adopted route: the correlated-k band
mean, corrected onto the window N2O occupies. Cross-check: Ramanathan (1976)
Appendix A over Dickinson's seven bands.

**The comparison that means something is ABSORPTANCE, not transmissivity, and
that is the finding worth keeping.** The two conventions average over different
spectral widths -- the corrk band is 83.8 cm-1, Ramanathan's `2 A0` is 34.6 --
so their transmissivities are not the same quantity. Compared that way they
disagree by up to **6.2x** and the disagreement is pure bookkeeping. Compared as
absorptance in cm-1, which carries no width convention, they agree to a factor
of 2.02 across three decades of CO2 amount.

That trap is worth stating plainly because both quantities are numbers between
zero and one that fall with CO2 amount, so taking one where the other belongs
would have put a factor of several into the model with nothing to catch it.

| CO2, atm cm | T corrk | T Ramanathan | A corrk | A Ramanathan | adopted multiplier |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.980 | 0.976 | 1.67 | 0.83 | 0.988 |
| 10 | 0.896 | 0.805 | 8.75 | 7.50 | 0.935 |
| 100 | 0.665 | 0.292 | 28.1 | 42.5 | 0.780 |
| 272, the whole column | 0.542 | 0.132 | 38.4 | 70.1 | **0.689** |

**The declared bar was 2.0 and the result is 2.02. It is recorded as a MISS and
the bar is not moved.** The worst point is the optically THIN end rather than
the thick one: at 1 atm cm the band model gives half the line list's
absorptance, because Ramanathan's effective intensity assumes a pure
exponential falloff from 667 cm-1 and at 589 the real band has structure that
falloff does not describe. In the thin limit absorptance is proportional to that
intensity, so the error arrives undiluted. At the thick end they diverge the
other way, the corrk absorptance saturating toward its own band width while the
band model's grows logarithmically without bound.

What the check was for is gross error -- an order of magnitude, a unit slip, the
wrong band -- and there is none. Agreement inside a factor of two across three
decades, against a 1976 band model whose line spacing had to be ASSUMED rather
than cited, supports the correlated-k route as the adopted one. It does not
validate the band model, and it is not evidence about the corrk route's own
accuracy, which rests on HITRAN 2020.

**The window correction is 0.6097 and is derived, not fitted.** The corrk band
spans 546-630 cm-1 while the N2O band occupies 566-612, and the corrk band
reaches further toward the 667 cm-1 fundamental where CO2 absorbs far harder.
Weighting both windows by the same exponential falloff the band model uses, the
corrk band's mean CO2 absorption is 1.64x the absorption over the window N2O
actually occupies, so the band-mean optical depth is scaled by 0.6097 before it
becomes the multiplier. Uncorrected it would understate the transmissivity N2O
sees and so understate the 589 band.

**The far wings are applied and are negligible here.** `k` at 589 cm-1 and 250 K
is 1.95e-05, which moves the whole-column transmissivity from 0.5439 to 0.5424.
Carried anyway, because the corrk tables are built with a 25 cm-1 line cutoff
and leaving their companion out would be an omission rather than a
simplification.

**Sensitivity to the assumed line spacing**, which is the cross-check's one free
parameter: halving `D` to 0.78 cm-1 takes the band model's whole-column
transmissivity to 0.040 and doubling it to 3.12 takes it to 0.142, against
0.076 at the assumed 1.56. So the cross-check is good to about a factor of two
on that alone, which is the same size as the disagreement being measured. It
cannot be sharpened without a line spacing someone can cite.

## 4d. What the implementation returned

Built and measured 2026-08-20 at T42 L10 on 8 ranks, 56 timesteps -- a segment
that crosses two writes, so the comparison can fail. CLIM-44 settled that the
model is bit-reproducible at a fixed rank count, so a bit-identity claim is
valid here in a way it was not for CLIM-39.

**All four declared tests pass.**

**The reduction identity.** With both mixing ratios at zero, a binary carrying
this change is BIT-IDENTICAL to one built without it, on `plasim_status`,
`plasim_output` and `plasim_snapshot` alike. Zero is the model's default, so
the change is inert until a namelist turns it on.

**Sign and monotonicity.** Outgoing longwave falls when either gas is added,
everywhere, at every abundance tried.

**Bounded transmissivity.** A run at 100 ppmv of each -- the top of Byrne's
range -- produced no NaN and no abort.

**Magnitude.** At 1.600 ppmv of CH4 and 0.300 of N2O against zero, the global
mean outgoing longwave falls by **0.799 W/m2**, and no cell is unaffected.

### Reading the magnitude, and how far that reading goes

The offline pricing is 1.56 to 1.98 W/m2 above a 100 ppbv floor, so a
0-to-1.6 ppmv change should price a little above 2. The model applies 0.799,
about 40 percent of it.

**The scheme under-delivers on CO2 by a similar factor**, which is why 0.799 is
not read as a defect in this term. Measured on the same bed: the whole CO2
greenhouse is 12.28 W/m2, from turning CO2 off entirely, and a doubling gives
1.576. **The Earth figures those are held against -- roughly 25 to 30 W/m2 and
3.7 -- are quoted from general knowledge and have NOT been sourced or checked
here**, so "40 to 45 percent of a line-by-line answer" is an impression and not
a measurement. It is enough to say the two shortfalls are the same size and not
enough to say either is right.

**What the shortfall is NOT is this world's thinner column.** 1 bar over 12.81
m/s2 gives 0.756 of Earth's column for the same mixing ratio, and the obvious
reading is that a thinner column is a weaker greenhouse. That is directionally
true and quantitatively small. Measured by running the model at 595 ppm, which
is the CO2 amount that gives Earth's column at this gravity:

| CO2 | greenhouse |
| --- | ---: |
| 450 ppm, this world's column | 12.284 W/m2 |
| 595 ppm, Earth's column | 12.926 W/m2 |

**0.642 W/m2**, against a gap of order fifteen. The column reduction accounts
for about four percent of it. Note also that gravity is what does this, not
radius: the column over a square metre is `P/g` and carries no radius at all. A
larger planet holds more atmosphere in total and the same amount above each
point.

So the remaining shortfall is unattributed. It is consistent with a broadband
scheme against line-by-line, and this note does not establish that.

## 4e. The bug that hid all of this, and it was not physics

**For most of a day this term measured 0.0175 W/m2 instead of 0.799, and the
cause was a missing MPI broadcast.**

`radmod_nl` is read on NROOT only. Every namelist variable in `radini` is then
broadcast explicitly -- `call mpbcr(co2)` and forty others. `ch4` and `n2o`
were added to the namelist and NOT to that block, so on every rank but the root
they kept their default of 0.0, which means absent. At T42 on 8 ranks that left
the band running on rank 0's eight latitude rows and nowhere else.

**It is invisible in every way that matters.** The model does not warn. The run
completes. The namelist echo in `plasim_diag` shows the values, because that is
printed on NROOT too. The reduction identity still passes, because it tests the
gases OFF. And the symptom -- a term about eight times too weak, further
diluted because rank 0's rows are polar and small in area -- looks exactly like
a physics problem.

It was chased as one. The record of that is in this note's own history: a
Doppler-broadening hypothesis that sent `Cess (1973)` to be fetched and read,
a shortwave hypothesis, a claim that the scheme could not convert absorptivity
into flux, and a "40x contradiction" between two perturbations that were in
fact identical. None of those were true. What made it findable in the end was
dumping the flux PER RANK rather than as a global mean: `+0.0000` on ranks 1
through 7 and `+1.0456` on rank 0 is not a number physics produces.

**This is a class this project has recorded twice before.**
`notes/audits/nlowio-collective-deadlock.md` is a collective placed behind an
unbroadcast `nlowio`; PHYS-9 is a namelist key applied at prepare and not
reapplied per segment. Both are the same shape: a value that exists on one rank
or in one code path and silently defaults everywhere else.

The lesson that generalises, and the reason the diagnostic mattered: **a global
mean cannot distinguish a term that is weak from a term that is off in most of
the domain.** Any per-rank or per-region quantity that comes out exactly zero
is worth more attention than one that comes out small.

## 4f. What the false trails were worth

Recorded because two of them produced results that stand on their own.

**The band model was validated against a line list**, which is section 4c, and
that stands: Eq. (1) reproduces HITRAN 2020 to better than 7 percent at every
level. That was gathered to explain a shortfall that turned out not to exist,
and it is now the evidence that the band model is right.

**Cess (1973) was read and does not apply**, which is worth keeping so it is
not tried again. Its Doppler branch comes out below Lorentz at every level this
model has -- 0.31 of it at the top layer's 0.025 atm, 0.04 at the surface --
and the two would cross near 2 mbar against a model lid of 50. The Doppler
regime lies entirely above the simulated atmosphere. The paper says why in its
own words: a Lorentz line strong at low altitude stays strong as altitude
rises, and the top-layer line-centre optical depth measured here is about 53.

**The scheme's own CO2 calibration** -- 12.28 W/m2 for the whole CO2 greenhouse
and 1.576 for a doubling -- was measured to bound the shortfall and is now what
the 0.799 is read against. It is the more useful number of the two.

The Doppler branch is NOT implemented. Adding it would be dead code at this
model's pressures, and the reason is recorded here rather than in the source.

## 5. The tests, declared before the work

**Band model against its own source.** Eq. (1) with `A0` = 52 (T/300)^1/2,
`beta0` = 0.17 (300/T)^1/2 and `S` = 187.7 must reproduce the ten consistent
rows of Donner and Ramanathan's Table 2 to better than 2% relative. This is an
identity against printed numbers, not an agreement, and it fails if the band
model is transcribed wrongly in any of its three equations. It is checkable
without the model and belongs beside the implementation as a script.

**The zero-abundance identity.** With the CH4 and N2O mixing ratios set to
zero, the model must be BIT-IDENTICAL to a binary built without the term. Not
close: identical. That is the reduction identity CLIM-39 used, and it is the
only test that can catch the new absorptivity leaking into the existing three.
CLIM-44 governs how far it can be claimed -- the horizon has to be established
with a same-binary control first, and that is unfinished.

**Sign and monotonicity.** Raising either mixing ratio must lower outgoing
longwave radiation at the top of the atmosphere and must not raise the
clear-sky layer transmissivity anywhere. A grey absorber that warms in one cell
and cools in another is a coding error, and this catches a sign slip that the
global mean would average away.

**Magnitude against the offline estimate.** At Earth's own abundances and
Earth's gravity, the term must give a top-of-atmosphere forcing of order the
2.02 W/m2 that `analysis/trace_gas_forcing.py` computes independently by Byrne
and Goldblatt's fits. This is the weakest of the four and is stated as such:
the two calculations share no code but do share an atmosphere, so agreement is
evidence and disagreement is a question rather than a verdict. **It is not a
calibration target.** If the band model and the forcing fits disagree by a
factor, that is information about two published parameterisations, and the band
model is not to be adjusted to close it.

**Bounded transmissivity.** `ztaucs` is already clamped to `[zero, 1-zero]`
after the existing three absorbers. Adding a fourth and fifth absorptivity to a
sum that is subtracted from one makes saturation reachable at high abundance,
so the clamp has to be shown to be doing the work rather than assumed: a run at
100 ppmv, the top of Byrne's range, must produce no NaN and no negative
transmissivity.

## 6. What has to happen around it

CLAUDE.md rule 4: a change under `vendor/exoplasim` makes every binary stale,
and ExoPlaSim compiles one executable per (resolution, layers, ranks) triple.
Every one is rebuilt and `--verify` run before any result from this is quoted.

The file is `radmod.f90`, which CLIM-39 has already changed and landed, so
there is no longer a sequencing conflict with the aerosol work.
