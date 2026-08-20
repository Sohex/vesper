# CH4 and N2O: what the radiation cannot carry, and what it is worth

Worldbuilding. Vesper is an invented planet and this note is about the
simulation of it: a toy climate model's radiation scheme and the trace gases it
has no term for. Every quantity named here is a modelled field.

Priced 2026-08-20 for CLIM-41. `config/planet.yaml` already declares that
PlaSim's longwave is Sasamori (1968) over water vapour, CO2 and ozone, that
`radmod.f90` has no CH4 and no N2O term, and that the direction is one-signed.
What it says instead of a magnitude is "by whatever they would have been worth".
This is that number, as far as it can be bounded.

**It comes out as the largest single item in the error budget**, and the reason
it is larger here than the same omission would be on Earth is the star.

## 1. This star's mixing ratios, not Earth's

Segura et al. (2003) is already on disk and read; it supplies
`star.surface_uv_relative_to_earth` and `model.ozone_scale`. The same paper runs
55 chemical species over 217 reactions in a photochemical model coupled to a
radiative-convective one, for a K2V host beside the Sun and an F2V, and carries
CH4 and N2O through to its spectra.

Its result, p. 701, verbatim: "CH4 and N2O levels are systematically higher for
the K2V planet and lower for the F2V planet because of the difference in stellar
UV fluxes. Note that trace gas concentrations are different from Earth at *all*
O2 levels (including 1 PAL) because we have assumed the same surface fluxes as
for modern Earth."

That is the whole mechanism: same biology, less near-UV, slower photolysis,
longer lifetimes, higher steady-state concentrations. The surface fluxes of H2,
CH4, N2O, CO and CH3Cl are held at modern Earth values throughout.

**The values are read off a log-scale figure, so they are a bracket.** Fig. 12a
is the K2V case and Fig. 4 is the Sun case on the same axes shifted one decade.
At 1 PAL of O2, which is this world's declared oxygen:

| gas | Sun, Fig. 4 | K2V, Fig. 12a | ratio |
| --- | --- | --- | --- |
| CH4 | 1600 ppb | 4000 to 8000 ppb | 2.5 to 5 |
| N2O | 300 ppb | 600 to 1500 ppb | 2 to 5 |

**Three checks that the figure reading is sound.** Fig. 4's Sun values at 1 PAL
reproduce modern Earth's actual CH4 and N2O, which is what that panel is for and
what calibrates the axis. The text gives one K2V number independently of the
figure -- at 0.1 PAL the K2V CH4 is "already of the order of 20 ppm" -- which is
where reading Fig. 12a's peak puts it. And Rugheimer et al. (2013), a separate
model on the same method, states its Sun-case surface mixing ratios as numbers
rather than a figure: cCH4 = 1.6e-6 and cN2O = 3.0e-7, which are exactly the
values read here.

**Two caveats, and they pull opposite ways.** Segura's K2V is 4620 K against this
host's 4965 K, so this star is the warmer end of the class and the enhancement
here is the milder one. Against that, this star is DECLARED active, of the kind
epsilon Eridani represents, and Segura's is a quiet modelled K2V; chromospheric
emission puts back near-UV, shortens the lifetimes and cuts the enhancement.
That second caveat is inherited rather than new: `surface_uv_relative_to_earth`
already takes 0.4 from the same quiet-star calculation, so both config values
rest on one assumption about this star's UV, and it is the activity declaration
that is in tension with it.

**Rugheimer et al. (2013) is the route out of both caveats and has not been
fully exploited here.** It runs the same fixed-biogenic-flux method over the
whole FGK main sequence at 250 K intervals, so 4965 K falls essentially on its
5000 K grid point rather than needing interpolation from 4620 K; it drives the
photochemistry with IUE-observed stellar ultraviolet rather than a modelled
quiet spectrum; and **epsilon Eridani is one of its grid stars**, at 5090 K and
assigned to that same 5000 K point -- the star `config/planet.yaml` names as
what this world's activity is modelled on. What is taken from it here is its
Sun-case anchor and its statement about N2O; its per-grid-star surface mixing
ratios are in its Fig. 7 and would replace the bracket above with a value.

## 2. Gravity, which cuts the other way and must not be dropped

A mixing ratio is not a column. At 1 bar over 12.81 m/s2 against Earth's 1.013
bar over 9.81, the column for a given mixing ratio here is **0.756** of Earth's.
That is the same 23% reduction in atmospheric column mass that
`notes/audits/absent-and-inherited-physics.md` confirms is carried correctly
through every other absorber, and it is why 450 ppm of CO2 here really is a
smaller CO2 column than 450 ppm on Earth.

Every number below is gravity-corrected: this world's mixing ratio times 0.756,
which is the Earth-equivalent column the Earth-calibrated expressions expect.

## 3. What it is worth

Byrne and Goldblatt (2014) recalculate greenhouse forcing over **100 ppbv to
100 ppmv** for both gases and give fits over that whole span, split at 2.5 ppmv
because CH4 and N2O cross from the square-root regime into the logarithmic one
there. That range covers this world at both ends, where Etminan et al. (2016)
does not: Etminan is valid only to 3500 ppb of CH4 and 525 ppb of N2O.

Three checks that the fits are being applied correctly, all of which pass:

- They are CONTINUOUS at the 2.5 ppmv junction. The square-root branch and the
  logarithmic branch both give 0.824 W/m2 for CH4 and 4.182 for N2O, which is
  what confirms the constants against a table that does not survive text
  extraction cleanly. The corrected version of Table 2 is the one used; the
  original was published with several pairs of square brackets missing.
- They reproduce the paper's own stated maxima at 100 ppmv, 6.50 against 6.66
  W/m2 for CH4 and 22.69 against 22.3 for N2O.
- They reproduce Earth's known forcings: 0.564 W/m2 for CH4 from 715 to 1800
  ppb, and 0.193 for N2O from 270 to 324 ppb.

**And where Etminan is valid the two methods agree to 3.5%**, which is the
cross-check that matters, because the two are independent line-by-line
calculations fitted by different groups.

Every number below is gravity-corrected and is quoted above Byrne's 100 ppbv
floor, so each still excludes the residual from zero to 100 ppbv that no fit
here covers. The same calculation for Earth's own atmosphere at Earth's gravity
gives 2.02 W/m2, which is the figure to hold these against:

| | CH4 | N2O | total, W/m2 |
| --- | --- | --- | ---: |
| A. N2O at the Sun value, per Rugheimer | 3.02 to 6.05 ppmv | 0.23 ppmv | **2.25 to 2.95** |
| B. N2O enhanced, per Segura Fig. 12 | 3.02 to 6.05 ppmv | 0.45 to 1.13 ppmv | **3.02 to 5.26** |

### The two sources disagree about N2O, and the disagreement is the bracket

They agree about CH4. Rugheimer: "CH4 abundance increases with decreasing
stellar temperature, dominated by the effects of decreasing stellar UV."

They do not agree about N2O. Rugheimer, on the same grid: "Up to about 20 km,
there is no significant difference between stellar types in N2O concentration.
Above ~20 km, Fig. 7 shows a decrease in N2O concentration for atmosphere models
around hot compared to cool grid stars ... Below 20 km N2O is shielded from
photolysis by the O3 layer." If the enhancement is confined above 20 km it
barely reaches the forcing, because the column and the forcing are both
troposphere-weighted.

Segura's Fig. 12 is captioned surface concentrations and reads higher for the
K2V planet at 1 PAL. That is a value read off a log axis at the edge of a plot,
which is the least reliable place to read one, and it is the weaker of the two
claims. **The shielding argument is a mechanism and the figure reading is not**,
so case A is the more likely of the two and case B is retained as the upper end
rather than as an equal alternative.

This is now where the width comes from. The forcing half is settled to a few
percent by two independent methods; the spread above is one physical
disagreement about tropospheric N2O, and resolving it is CLIM-43.

## 4. Why this ranks first, and what it does to the design flux

At 2 W/m2 and up this exceeds every priced item in `analysis/error_budget.json`.
Mineral dust is +0.34 to +0.61 and stood second in the whole budget; sea salt is
-0.16 to -0.89; the energy closure residual is -0.40 to -0.49. On the budget's
own stopping rule -- refine an input when its plausible range exceeds the effect
of the thing you last refined -- this outranks the aerosol work.

**`derive_design_flux.py` is absorbing it, and by construction.** The design flux
is SEARCHED to satisfy thresholds declared in advance on the modelled climate. A
radiation scheme missing 2 W/m2 of greenhouse forcing reaches those thresholds
at a higher stellar flux, and the search finds it. At the recorded 202 K per unit
flux ratio the compensation is small in flux-ratio terms and exact in
temperature: the world arrives at the right temperature by too much starlight
and too little greenhouse.

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
1e-3**, which is 1000 ppmv. This world is near 4 ppmv, more than two orders of
magnitude below it, and the star's near-infrared share is larger than the Sun's
by a factor of order one rather than of order a hundred. The effect is real, it
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
