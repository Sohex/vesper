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

**For N2O the two intensities are NOT recovered, and both routes have been
tried.** The paper's evidence for the 1285 cm-1 band is Fig. 2 rather than a
table, and the 589 cm-1 band has neither.

*The figure route was tried and its own check refuted it.* Fig. 2 plots Eq. (1)
at two pressures, which is two independent constraints on one unknown, so it
self-checks. It fails: tracing curve A by continuity off a 600 dpi render and
fitting over the SAME absorber range gives S = 302 from the 0.5 atm panel and
S = 348 from the 0.1 atm panel, while each fits its own trace to an rms of 0.2
cm-1. A 15% systematic disagreement between panels that individually fit that
tightly is a calibration error in the extraction, not noise, and a number taken
from either panel alone would carry it invisibly. The extraction is discarded.
This is what the two-panel design was for; it did its job by failing.

*The primary source is not obtainable.* Donner and Ramanathan never state the
N2O intensities and cite McClatchey et al. (1973), *AFCRL Atmospheric Absorption
Line Parameters Compilation*, AFCRL-TR-73-0096 -- a 1973 technical report with
no DOI, which `paperfetch` cannot identify. Ramanathan (1976) was fetched on the
chance it tabulated them and does not: it carries no N2O at all, being H2O, CO2
and O3, and what it contributes is the band model itself.

**So this note cannot support an N2O implementation**, and a value guessed to
fill the gap would be the failure this project calls precision theatre. CH4
alone is the larger half of the term and is fully specified, so the sensible
shape is CH4 first with the N2O slot left explicitly empty. The open routes,
neither taken: obtain AFCRL-TR-73-0096 from DTIC, or sum HITRAN line intensities
over the two bands, which is a different source from the one the band model was
fitted against and would need saying so.

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
  already computes, and for the 589 cm-1 N2O band by the CO2 overlap as well;
- subtracted from `ztaucs` alongside the existing three.

Mixing ratios come from `config/planet.yaml` through a namelist key, on the
pattern `dqco2` already sets, with the values CLIM-43 measured. They are an
ASSUMPTION of this world in the same sense 450 ppm of CO2 is, and the config
should say so where it declares them.

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
