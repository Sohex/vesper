# CH4 and N2O: what the offline calculation says they are worth

Worldbuilding. Vesper is an invented planet and this note is about the
simulation of it: a toy climate model's radiation scheme and two trace gases in
it. Every quantity named here is a modelled field.

Priced 2026-08-20 for CLIM-41; the mixing ratios measured 2026-08-20 for
CLIM-43. This is the OFFLINE number, from line-by-line fits and independent of
any band in the model. `radmod.f90` has carried both gases since CLIM-42, and
`exoplasim/notes/trace-gas-band.md` measures that band against the number
below: it returns 40 per cent of it, which is what this broadband scheme
returns for CO2 as well. So this note is the reference, and the difference is
the residual the scheme still misses.

**It is the largest single item in the error budget.** The derivation is
`analysis/trace_gas_forcing.py`, which writes `analysis/trace_gas_forcing.json`
and carries the checks; this note is the argument.

## 1. This star's mixing ratios, not Earth's

The gases are photochemical products of biogenic surface fluxes, so their
steady-state abundance is set by the host star's ultraviolet rather than by any
choice made here. Segura et al. (2003) and Rugheimer et al. (2013) both run
that calculation with the surface fluxes of H2, CH4, N2O, CO and CH3Cl held at
modern Earth values throughout: same biology, different star.

**Rugheimer's grid is the one that answers this world**, and its 5000 K point
is not an interpolation. It runs 4250 K to 7000 K at 250 K spacing, drives the
photochemistry with IUE-OBSERVED stellar ultraviolet rather than a modelled
quiet spectrum, and **epsilon Eridani is one of its grid stars**, at 5090 K and
assigned to the 5000 K point -- 35 K from this host, and the star
`config/planet.yaml` names as what this world's activity is modelled on.

Its Fig. 7 is a figure, but the paper is a vector PDF, so the profiles come out
as path coordinates rather than as a reading off a log axis. At the 5000 K grid
point:

| gas | this host, 5000 K | Rugheimer's Sun case | ratio |
| --- | --- | --- | --- |
| CH4 | 1.600 ppmv | 1.6 ppmv | 1.00 |
| N2O | 0.300 ppmv | 0.3 ppmv | 1.00 |

**Neither gas is enhanced at this star.** Both are well mixed through this
world's troposphere in that model -- CH4 falls 1.600 to 1.583 ppmv over the
lowest 20 km, N2O is flat to 15 km -- so the surface value is the column value,
which is what the forcing wants.

### Segura and Rugheimer do not disagree; they are at different temperatures

This was carried as a bracket for a day, on the reading that Segura's Fig. 12a
puts a K2V planet's CH4 at 4000 to 8000 ppb and its N2O at 600 to 1500 ppb
against the Sun's 1600 and 300, while Rugheimer says of the same method that
"up to about 20 km, there is no significant difference between stellar types in
N2O concentration". Both readings are correct. The enhancement is real, and it
turns on sharply at the cold end of Rugheimer's grid:

| grid star | CH4, ppmv | N2O, ppmv |
| --- | --- | --- |
| 5250 K | 2.464 | 0.342 |
| 5000 K | 1.600 | 0.300 |
| 4750 K | 1.600 | 0.300 |
| 4500 K | 37.27 | 1.654 |
| 4250 K | 124.5 | 0.840 |

**Segura's K2V is 4620 K, inside that transition**, which is exactly why it
reads a few ppmv of CH4. This host is 4965 K, above it. The two sources are
consistent and the apparent conflict was a comparison between stars 345 K
apart.

That also closes the caveat the earlier estimate could not. Segura's K2V is a
quiet MODELLED spectrum while this star is declared active, of the kind epsilon
Eridani represents, and the two pull opposite ways on lifetime. Rugheimer
settles it by using the observed ultraviolet of epsilon Eridani itself at the
grid point this world sits on. Note what is NOT thereby settled:
`star.surface_uv_relative_to_earth` still takes 0.4 from the same quiet-star
calculation, so that value rests on an assumption this one no longer needs.

### What the width is now

The model's own scatter over the 250 K either side of this host, taken as the
envelope of the 4750 / 5000 / 5250 K points. The two points that actually
bracket 4965 K are degenerate, so a bracket from those alone would report zero
width -- a statement about the grid rather than about the world. The cold-star
regime below 4750 K is EXCLUDED rather than bounded, because the activity an
extrapolation toward it would represent is already represented at 5000 K by
epsilon Eridani.

## 2. Gravity, which cuts the other way and must not be dropped

A mixing ratio is not a column. At 1 bar over 12.81 m/s2 against Earth's 1.013
bar over 9.81, the column for a given mixing ratio here is **0.756** of Earth's.
That is the same 23% reduction in atmospheric column mass that
`notes/audits/absent-and-inherited-physics.md` confirms is carried correctly
through every other absorber, and it is why 450 ppm of CO2 here really is a
smaller CO2 column than 450 ppm on Earth.

Every number below is gravity-corrected: this world's mixing ratio times 0.756,
which is the Earth-equivalent column the Earth-calibrated expressions expect.
It is the whole reason this world's term comes out BELOW Earth's despite an
identical mixing ratio.

## 3. What it is worth

Byrne and Goldblatt (2014) recalculate greenhouse forcing over **100 ppbv to
100 ppmv** for both gases and give fits over that whole span, split at 2.5 ppmv
because CH4 and N2O cross from the square-root regime into the logarithmic one
there. That range covers this world at both ends, where Etminan et al. (2016)
does not: Etminan is valid only to 3500 ppb of CH4 and 525 ppb of N2O.

Four checks that the fits are being applied correctly, all of which pass and
all of which `trace_gas_forcing.py` raises on:

- They are CONTINUOUS at the 2.5 ppmv junction. The square-root branch and the
  logarithmic branch both give 0.824 W/m2 for CH4 and 4.182 for N2O, which is
  what confirms the constants against a table that does not survive text
  extraction cleanly. The corrected version of Table 2 is the one used; the
  original was published with several pairs of square brackets missing.
- They reproduce the paper's own stated maxima at 100 ppmv, 6.50 against 6.66
  W/m2 for CH4 and 22.69 against 22.3 for N2O.
- They reproduce Earth's known forcings: 0.564 W/m2 for CH4 from 715 to 1800
  ppb, and 0.193 for N2O from 270 to 324 ppb.
- The extracted Sun curve reproduces Rugheimer's separately STATED Sun-case
  mixing ratios, to 1% on N2O and 8% on CH4. That is the figure extraction
  checking itself against a number printed in the same paper, and the 8% is
  the honest precision of the CH4 panel rather than a calibration error: the
  N2O panel is read the same way and lands at 1%.

**And where Etminan is valid the two methods agree to 3.5%**, which is the
cross-check that matters, because the two are independent line-by-line
calculations fitted by different groups.

Each gas is quoted above Byrne's 100 ppbv floor, so each still excludes the
residual from zero to 100 ppbv that no fit here covers. The same calculation
for Earth's own atmosphere at Earth's gravity gives 2.02 W/m2, which is the
figure to hold these against:

| | CH4 | N2O | total, W/m2 |
| --- | --- | --- | ---: |
| value, the 5000 K grid point | 1.209 ppmv | 0.227 ppmv | **1.56** |
| the 250 K envelope | 1.209 to 1.863 ppmv | 0.227 to 0.258 ppmv | **1.56 to 1.98** |

Both columns are Earth-equivalent, gravity applied.

**Two omissions with known signs, and they oppose.** The sub-100-ppbv residual
holds these numbers low. The CH4-N2O band overlap holds them high: it is about
-0.02 W/m2 for Earth's own atmosphere, and it comes out zero here not because
it is negligible but because this world's Earth-equivalent N2O column of 0.227
ppmv is below the 0.270 ppmv reference the overlap expression is written
against, so the expression has nothing to say. An overlap can only reduce a
total, so the sign is known even where the magnitude is not.

## 4. Why this ranks first, and what it does to the design flux

At 1.56 to 1.98 W/m2 this exceeds every priced item in
`analysis/error_budget.json`. Mineral dust is +0.34 to +0.61 and stands second
in the whole budget; sea salt is -0.16 to -0.89; the energy closure residual is
-0.40 to -0.49.

The budget's stopping rule -- refine an input when its plausible RANGE exceeds
the effect of the thing you last refined -- now reads differently from the
magnitude. The range is 0.42 W/m2, which is no longer the widest thing in the
budget. **What outranks everything is the term itself, not the uncertainty in
it**, and the response to a large term that is known to a few percent is to
carry it, not to measure it again. That is CLIM-42.

**`derive_design_flux.py` is absorbing it, and by construction.** The design
flux is SEARCHED to satisfy thresholds declared in advance on the modelled
climate. A radiation scheme missing 1.6 W/m2 of greenhouse forcing reaches
those thresholds at a higher stellar flux, and the search finds it. At the
recorded 202 K per unit flux ratio the compensation is small in flux-ratio
terms and exact in temperature: the world arrives at the right temperature by
too much starlight and too little greenhouse.

That does not make the climate wrong. It makes the ATTRIBUTION wrong, and it
moves the error somewhere the temperature does not show it. Everything keyed on
insolation rather than on temperature inherits the bias with the same sign:
photosynthetically active flux in the biosphere, the shortwave band split, the
aerosol forcings that scale with incident flux, and the surface UV that
`star.surface_uv_relative_to_earth` reports. `config/planet.yaml` already asks
that 450 ppm of CO2 be read as an assumption rather than a result; the derived
flux now needs the same label, and for a reason that has nothing to do with CO2.

## 4b. Checked and does not apply: the shortwave offset

Recorded so it is not re-derived. CH4 absorbs in the near-infrared as well as
the thermal infrared, and this star puts more of its flux there than the Sun
does, so shortwave absorption offsetting the greenhouse warming is the obvious
way this estimate could be too large. Byrne and Goldblatt (2015) is the paper on
exactly that effect, and it puts the diminishment at CH4 abundances **above
1e-3**, which is 1000 ppmv. This world is near 1.2 ppmv, nearly three orders of
magnitude below it, and the star's near-infrared share is larger than the Sun's
by a factor of order one rather than of order a thousand. The effect is real, it
is directionally relevant here, and it is negligible at these concentrations.

## 5. What would settle it

Adding the band. `config/planet.yaml`'s own sentence is the right one: "Adding
either means adding a band, not a key." Sasamori (1968) is a broadband
parameterisation and CH4's 7.7 um band sits inside the window the scheme already
treats, so this is a radiation-scheme change of the same kind as CLIM-39's
aerosol work rather than a namelist setting, and it lands on the same file.

Nothing here should be read as a case for tuning. The bound says the term is
large, not that any particular value is right, and a term added because it
improves an agreement is `docs/src/practice/failure-modes.md` class 16. What
justifies the band is that the gas is there.
