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

## 4d. What the implementation returned, including the test it misses

Built and measured 2026-08-20 at T42 L10 on 8 ranks, 56 timesteps -- a segment
that crosses two writes, so the comparison can fail. CLIM-44 settled that the
model is bit-reproducible at a fixed rank count, so a bit-identity claim is
valid here in a way it was not for CLIM-39.

**The reduction identity PASSES.** With both mixing ratios at zero, a binary
carrying this change is BIT-IDENTICAL to one built without it, on
`plasim_status`, `plasim_output` and `plasim_snapshot` alike. Zero is the
model's default, so the change is inert until a namelist turns it on.

**Sign and saturation PASS.** Adding the gases lowers outgoing longwave:
238.759 W/m2 with them off, 238.742 at this world's abundances, 238.551 at 100
ppmv of each. The 100 ppmv run -- the top of Byrne's range -- produced no NaN
and no abort, which is the bounded-transmissivity test.

**THE MAGNITUDE TEST MISSES, BY A FACTOR OF ABOUT NINETY, AND IT IS NOT A
CODING ERROR.** The offline estimate is 1.56 to 1.98 W/m2 and the model returns
0.0175. What rules out a coding error is that the same model, on the same bed,
gives 1.576 W/m2 for a CO2 doubling, and that the per-layer absorptivities are
individually right: at the full column CH4 contributes 0.00475, N2O 1285
contributes 0.00269 and N2O 589 contributes 0.00445, a total of 0.0119 against
0.0115 for the CO2 doubling. Two changes of the same size in the same variable,
and the fluxes respond ninety times differently.

**The difference is where in the column they act.** Instrumented per path from
the top:

| path | CO2 amount | my added absorptivity | a CO2 doubling's |
| --- | ---: | ---: | ---: |
| top layer only | 13.7 | 0.00076 | 0.0164 |
| top three | 65.3 | 0.0032 | 0.0164 |
| whole column | 273.7 | 0.0119 | 0.0115 |

CO2's absorptivity is logarithmic, so a doubling adds the same 0.0164 whatever
the path. This band model is linear in amount where the path is thin, so aloft
it adds almost nothing -- and top-of-atmosphere forcing is made in the upper
troposphere and above, where the emitting temperature differs most from the
surface.

**The suspected cause is that the band model is Lorentz-only.** At the top
layer the path pressure is 0.025 atm, and there Eq. (1) returns 2.66 cm-1 of
absorptance where the weak-line limit `S W` would give 9.1 -- the line shape
parameter `beta = beta0 P/P0` is suppressing it by 3.4x. At that pressure
Doppler broadening is not negligible and a pressure-broadened formulation
understates the absorptance. Donner and Ramanathan applied their model to the
troposphere, where that assumption holds.

The candidate fix is named in Ramanathan's own bibliography and is now on disk:
Cess (1973), *A band absorptance formulation for Doppler broadening*, which
extends this band model to exactly this regime. It has not been read or
implemented, and until it is **the term must not be switched on**: it would
put a known ninety-fold understatement into a run while the config declares
mixing ratios that say otherwise.

So the state is: the machinery is in, inert, and its three structural tests
pass. What it computes in the troposphere matches an independent band-model
calculation. What it does not yet do is produce the forcing the pricing note
says these gases are worth, and the reason is a documented limitation of the
formulation aloft rather than of the code.

## 4e. Cess (1973) was read, and it does not apply

Recorded so it is not tried again. Section 4d suspected the Lorentz-only band
model of understating absorptance aloft, and named Cess (1973) as the fix. That
was WRONG, and reading the paper is what shows it.

Cess gives the Doppler analogue in the same variables:

    A_D = A0 u (1 - 0.18 u/delta)                    for u/delta <= 1.5
    A_D = 0.753 A0 delta {[ln(u/delta)]^1.5 + 1.21}  for u/delta >= 1.5

with `u = S P H / A0` as before and `delta = sqrt(pi) gamma_D / d`, `gamma_D`
the Doppler half-width and `d` the mean line spacing. The prescription for
combining the two is Goody and Belton's, which Cess adopts: take whichever
mechanism gives the LARGER absorptance.

`d` does not need a new source. The model's own definition of the line-structure
parameter is `beta = 4 gamma_L P / d`, so `d = 4 gamma_L / beta0` recovers it
from the `beta0` already in Donner and Ramanathan's Table 1, with `gamma_L` the
mean Lorentz half-width from McClatchey -- 0.055 cm-1 atm-1 for CH4, which it
adopts explicitly for all CH4 lines, and 0.082 for N2O from Toth's J-resolved
table. That gives 1.29 cm-1 for CH4 and 0.293 for N2O.

**Doppler absorptance is BELOW Lorentz at every level this model has:**

| band | 0.025 atm | 0.119 atm | 0.497 atm |
| --- | ---: | ---: | ---: |
| CH4 1306 | 0.31 | 0.10 | 0.04 |
| N2O 1285 | 0.32 | 0.13 | 0.07 |
| N2O 589 | 0.50 | 0.22 | 0.10 |

So Goody and Belton's rule keeps the Lorentz value everywhere, and adding the
Doppler branch would be dead code. Extrapolating the top row, the two would
cross near 2 mbar; this model's lid is `PTOP` = 50 mbar, so the Doppler regime
lies entirely ABOVE the atmosphere being simulated.

The paper also explains why, in its own words: the volumetric absorption
coefficient at the centre of a Lorentz line is independent of pressure, so "if
a Lorentz line is strong at low altitudes, it will remain strong as altitude is
increased". Checked here at the top layer, the line-centre optical depth is
about 53. The lines are saturated, which is the regime Eq. (1) is for.

## 4f. Where the shortfall is not, and what would find it

The scheme was calibrated against its own CO2 to see whether it is uniformly
weak. Turning CO2 off entirely takes OLR from 238.759 to 251.043, so **this
model's whole CO2 greenhouse effect is 12.28 W/m2**, against roughly 25 to 30
for Earth's. The same ratio shows in the doubling, 1.576 against about 3.7. So
the scheme runs at 40 to 45 percent of reality for CO2 -- weak by a factor of
two, and that is a property of a broadband scheme rather than a defect.

Against that calibration the trace-gas term is still about a HUNDRED times too
weak, and Doppler broadening is no longer a candidate. What is established:

- the absorber amounts are right, checked against an independent calculation;
- the band absorptances are what Donner and Ramanathan's Eq. (1) gives, and
  Eq. (1) is the correct branch at these pressures per Cess;
- the conversion to Sasamori's currency is consistent with Sasamori's own
  normalisation to within about 1.6x, checked by running CO2's band absorptance
  through the same conversion and comparing with `zaco2`;
- the flux machinery responds correctly to a change of the same size in the
  same variable, which is what the CO2 doubling shows.

**The next test is the one that worked for the CO2 overlap: put the band model
against a line list.** The LMD bundle on this host carries CH4 correlated-k
tables -- the `early_earth_CO2_*_CH4_*` sets -- so the band absorptance can be
compared with a HITRAN 2020 band mean at the same amount, pressure and
temperature, exactly as section 4c did for CO2. That separates a band model
that is wrong from a scheme that cannot carry the term, and it needs no new
data. It has not been done.

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
