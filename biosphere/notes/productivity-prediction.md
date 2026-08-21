# Pre-registered productivity prediction

Written before LPJ-GUESS has run a single Vesper gridcell, so that the model's
answer can be scored rather than rationalised. Every coefficient below is fixed
as of this document. If a prediction misses, the miss is the result; nothing here
gets retuned to fit.

## Rules

1. The **method and its coefficients** are what is registered, not only the
   numbers. The numbers below are that method evaluated on the only described
   climatology that exists today, which is not the one LPJ-GUESS will be given.
   When the real climatology lands, re-run the same arithmetic on it, unchanged,
   and score against that. Re-deriving is allowed; re-tuning is not.
2. Scoring compares against LPJ-GUESS's own outputs: `anpp.out`, `lai.out`,
   `fpc.out`, `cmass.out`.
3. A prediction that misses is recorded as a miss in this file, with the actual
   value, and the mechanism that broke is written down. It is not deleted.

## What this predicts against

The climatology LPJ-GUESS should consume does not exist yet. It needs to be
`carved-zoned` terrain under the corrected `k25v` spectrum. What exists is
`analysis/climatology_s096`, which is `precarve-unzoned` at 292.97 K under the
wrong `k2` spectrum, and which `world_state.json` already labels "never as
Vesper's climate".

The numbers below are computed on that superseded climatology because it is the
only one there is. Two known shifts are already in flight and both move the
answer:

- The carved-zoned baseline converged at **295.18 K**, 2.2 K warmer, with land
  albedo 0.1722 against 0.2240. Warmer and darker, so more evaporative demand.
  Expect the water-limited fraction to rise and the temperature-limited fraction
  to fall.
- The `k25v` spectrum brightens snow and ice by 0.10 to 0.17, worth roughly
  -0.4 to -0.7 K. Small, and opposite in sign to the carve.

So treat the absolute numbers as provisional and the **ratios and the mechanism**
as the real prediction.

## Assumed LPJ-GUESS configuration

The prediction assumes the porting decisions in
[`lpj-guess-porting-audit.md`](lpj-guess-porting-audit.md):

| item | value |
| --- | --- |
| calendar | 24-hour steps, 181-day year |
| `FRADPAR` | 0.4624, from the k25v spectrum over a 400-750 nm window |
| PFT set | Earth's, as an Earth-analogue biosphere |
| photosystem window | 400-750 nm, K2V-adapted (see the amendment below) |
| `nfix_a` | 0.234, LPJ-GUESS's central value; range 0.102-0.367 |
| `gdd5min` | rescaled by 180.655 / 365.2569 = 0.4946 |
| CO2 | 450 ppm |
| nitrogen | deposition 0.5 kgN/ha/yr declared; fixation and mineralisation dominate it, see the amendment |

If any of those change, the prediction is void rather than wrong.

## A unit trap that will look like a factor-of-two error

LPJ-GUESS reports fluxes per **simulation year**, which under the registered
calendar is 181 days, a Vesper year. Every number in this document is per
**Earth** year. Convert before comparing:

    NPP_earth_year = NPP_lpj_output x 365.2569 / 180.655 = NPP_lpj_output x 2.022

A raw comparison of LPJ-GUESS output against these figures will appear to miss
low by almost exactly two. That is the calendar, not the biosphere.

## Earth reference, pinned

So the comparison cannot drift as the conversation moves on:

| quantity | value used |
| --- | --- |
| land area | 149e6 km2 |
| terrestrial NPP | 55-60 PgC/yr |
| land-mean NPP | 370-400 gC/m2/yr |
| land-mean annual temperature | 8.6 C |
| land-mean annual precipitation | 750 mm/yr |
| land-mean surface downward shortwave | 180-190 W/m2 |
| PAR fraction of shortwave | 0.50 |

## The supply side

**Light.** The star is dimmer and redder, but the sky is clearer, and the second
nearly cancels the first.

| | Vesper (k25v) | Sun | ratio |
| --- | --- | --- | --- |
| PAR energy fraction of shortwave | 0.3093 | 0.3903 | 0.793 |
| PAR photons per J of shortwave | 0.1735 | 0.2132 | 0.814 |
| mean PAR wavelength | 561.0 nm | 546.3 nm | |
| land-mean surface downward SW | 195.7 W/m2 | 180-190 | 1.03-1.09 |
| land-mean cloud fraction | 0.525 | ~0.67 | |
| **PAR photons per m2 of land** | | | **0.84-0.89** |

Photon flux is the right currency, not energy: photosynthesis counts quanta.
Because this star's PAR is redder, each joule carries more photons, which
recovers about 3% of the energy deficit. Small, but it is the correct quantity
and it is registered as such.

**Climate.** Area-weighted over land, on the superseded climatology:

| quantity | Vesper | Earth |
| --- | --- | --- |
| annual mean temperature | 17.34 C | 8.6 C |
| annual precipitation | 892 mm per Earth year | 750 mm |
| land area | 314.7e6 km2 | 149e6 km2 |
| land area ratio | **2.11x** | |

The land is 43.17% of a sphere 1.44 times Earth's area. That factor, not any
per-hectare advantage, is where most of the aggregate productivity comes from.

**Empirical productivity.** The Miami model (Lieth 1975), NPP as the minimum of a
temperature and a precipitation limit, converted at 0.45 gC per g dry matter.
Applied at land means for both worlds, so the model's nonlinearity cannot
corrupt the comparison:

    Vesper (17.34 C, 892 mm)  603 gC/m2/yr   precipitation-limited
    Earth  ( 8.60 C, 750 mm)  530 gC/m2/yr   temperature-limited
    ratio 1.14

Applied per gridcell and area-averaged, Vesper gives 443 gC/m2/yr. The gap
between 603 and 443 is spatial heterogeneity against a concave function; Earth
would take a similar penalty, which is why the like-for-like ratio is the figure
carried forward.

## The predictions

Combining: per-area 1.14 for climate, times 0.86 to 0.95 for light. The light
correction is not applied in full because Miami is Earth-calibrated and already
assumes solar illumination, and because 60% of this land is water-limited, where
less light does not translate into less growth.

| # | prediction | central | scored a hit if | a clear miss beyond |
| --- | --- | --- | --- | --- |
| 1 | land-mean NPP relative to Earth's | 1.0x | 0.7-1.4x | 0.5-2.0x |
| 2 | land-mean NPP, absolute | 400 gC/m2/yr | 280-560 | 200-800 |
| 3 | total terrestrial NPP relative to Earth | 2.1x | 1.6-2.8x | 1.2-3.5x |
| 4 | total terrestrial NPP, absolute | 125 PgC per Earth year | 90-175 | 60-240 |
| 5 | tree cover, fraction of land with tree FPC > 0.1 | 55% | 40-70% | below 25 or above 80 |
| 6 | effectively barren land, LAI < 0.5 | 15% | 8-25% | above 40 |
| 7 | water-limited rather than temperature-limited land | 60% | 50-70% | below 35 |
| 8 | grass PFTs exceed tree PFTs in land area covered | yes | | |
| 9 | C4 grasses are a substantial presence, above 10% of land | yes | | |
| 10 | boreal needleleaf stays marginal, below 10% of land | yes | | |

Predictions 8 to 10 follow from the land being warm (52.3% above 20 C annual
mean, 70.4% with a warmest month above 22 C) and dry over most of its area,
while only 17.5% has a coldest month below -15 C and the short year halves the
degree-day budget that boreal trees need.

Supporting structure, from the same climatology, registered so the shape of the
prediction is checkable and not just its totals:

| land structure | fraction |
| --- | --- |
| hyper-arid, P < 250 mm | 14.2% |
| semi-arid, 250-500 mm | 15.7% |
| subhumid, 500-1000 mm | 32.5% |
| humid, P > 1000 mm | 37.6% |
| frost season, coldest month < 0 C | 34.8% |
| Miami NPP > 300 gC/m2/yr | 63.0% |

## Known biases in this prediction

Stated now so they cannot be discovered afterwards as excuses.

- **Miami is coarse and Earth-calibrated.** It knows nothing about CO2, nutrient
  limitation, fire, or competition. It implicitly assumes solar illumination,
  so multiplying by a light correction risks double-counting. That is the single
  largest reason the tolerances are wide.
- **The climatology is the wrong one**, in two known ways, described above.
- **The evaporite subtraction is not made.** 12.65% of carved-zoned land is
  evaporite or playa clastics where nothing roots regardless of what the climate
  permits. Miami will assign those cells whatever their aridity allows, and the
  dry ones will already score low, but salt crust is barren even when wet. Expect
  the true figure to sit 5-10% below predictions 2 and 4 for this reason alone.
- **Nitrogen is a declared constant.** On Earth, N limitation is a first-order
  control on NPP. A world with no biosphere history has no basis for a deposition
  field, so whatever is chosen is an assumption, and prediction 2 inherits it.
- **The climatology assumes the answer.** It was run with
  `land_albedo_source: vegetated`, 0.15 on anything that can carry a canopy. If
  LPJ-GUESS returns substantially less canopy than that assumed, the climate that
  produced these numbers is not the climate that biosphere would sustain, and the
  loop has to turn again. Prediction 5 is therefore the most decision-relevant
  line in this table, not predictions 1 to 4.

## Not part of the prediction: the photosystem window

Earth's 400-700 nm PAR window is an accident of our star, and every prediction
above assumes it transplants unchanged. The literature says it probably would
not, and says so with unusual specificity for this exact spectral type.

**Kiang et al. 2007** ("Spectral signatures of photosynthesis II", Astrobiology
7:252-274) work through pigment coevolution with the host star and state that
pigments around F2V stars may peak in the blue, **K2V in the red-orange**, and M
stars in the near infrared. Vesper's host is K2.5V, so this is nearly a direct
hit rather than an extrapolation.

**Lehmer et al. 2021** (Frontiers in Astronomy and Space Sciences 8:689441) put
numbers on it, applying the Marosvolgyi and van Gorkom power-gain optimisation
model, which solves for the absorbance spectrum maximising photon energy captured
against thermal emission losses. Predicted peak absorbance wavelengths:

| star | predicted peaks |
| --- | --- |
| F2V | 468, 476 nm |
| G2V (Sun) | 644, 672 nm |
| **K2V** | **675, 711, 746 nm** |
| M1V | 753 nm |
| M5V | 987, 1050 nm |

So the expected shift for this star is modest, about 30 to 75 nm redward, and it
straddles rather than abandons the 700 nm boundary. A window of roughly 400-750
nm is what the literature points at.

**That is very nearly the window that restores parity.** Recomputing photon flux
on the k25v spectrum:

| assumed window | photon flux vs Earth's 400-700 nm |
| --- | --- |
| 400-700 nm (Earth's, what is registered above) | 0.81x |
| **400-750 nm (Lehmer K2V optimum)** | **0.99x** |
| 400-800 nm (chlorophyll f limit) | 1.16x |
| 400-900 nm | 1.47x |
| 400-1100 nm | 2.06x |

A 50 nm extension almost exactly cancels the dimmer, redder star. That is not a
coincidence so much as the point of the optimisation argument: pigments track the
photon supply.

**The oxygenic constraint bounds this, and it comes from our own config.**
`planet.yaml` declares `pO2_bar: 0.21`, which requires oxygenic primary
production. Kiang et al.'s more dramatic result, that M-star planets could exceed
Earth's productivity if useful photons extend to 1.1 um, is explicitly for
*anoxygenic* photosynthesis, which produces no oxygen and cannot sustain that
atmosphere. Oxygenic photosynthesis is capped near 800 nm on known biochemistry:
chlorophyll f absorbs at 706 nm in vitro and was isolated from stromatolite
cyanobacteria cultured under 720 nm light (Chen, Schliep, Willows, Cai, Neilan, Scheer (2010). *A Red-Shifted Chlorophyll.* Science 329(5997), 1318-1319. `10.1126/science.1191127`),
and with chlorophyll d it lets some cyanobacteria work at 700-800 nm via
far-red photoacclimation. The further claim that such pigments contribute over
20% of gross photosynthesis in natural biofilms is DECLARED here pending its
own source -- it is real literature this project has not yet fetched, and none
of the four window references contains it (REF-9). Nothing oxygenic reaches
1.1 um without a three- or four-photon scheme nobody has observed. So the bottom two rows of that table are unavailable
to a world with a 21% oxygen atmosphere, and the honest range is 0.81x to 1.16x,
not 0.81x to 2.06x.

**Caveats that keep this out of the registered prediction.** Optimal is not the
same as realised: Earth's chlorophyll a sits where it does partly because of the
energetics of splitting water and partly because of three billion years of
contingency, not because it is optimal for a G2V. Lehmer et al. also exclude
canopy and leaf structure, which is exactly where a land-plant biosphere does
most of its light harvesting. And LPJ-GUESS has no mechanism to represent a
different photosystem at all beyond `FRADPAR`.

**What this means for the prediction.** Running Earth PFTs picks the first row,
and the registered numbers therefore sit at the conservative end by about 23% in
photon supply. If a K2V-adapted photosystem is adopted later, raise `FRADPAR`
from 0.396 to about 0.48 and re-register; this document becomes void rather than
wrong. The swing is roughly +20% on light-limited productivity, which is inside
the tolerance bands above, so adopting it would not by itself falsify anything
here.

## Amendment: the photosystem window was widened, and the prediction stands

Registered against Earth's 400-700 nm window. That window has since been changed
to 400-750 nm as a deliberate worldbuilding decision, on Lehmer et al. 2021's
predicted K2V peaks of 675, 711 and 746 nm. `FRADPAR` moved from 0.3963 to
0.4624, a rise of 16.8%.

Under the rules at the top of this document that makes the registration void
rather than wrong, and it is re-registered here rather than quietly edited. The
numbers themselves do not move: the document already said a redward shift was
worth roughly +20% on light-limited productivity, which sits inside every
tolerance band in the table, and 16.8% duly does. **No prediction changes.** What
changes is that the registered numbers now sit nearer the centre of their bands
than the conservative edge.

Setting the window exposed a real bug in `build_vesper_header.py`, which had been
scaling the star's PAR fraction against the Sun's measured over the *same*
window. Widening the window widened both and cancelled most of the effect,
giving 0.4076 instead of 0.4624. Earth's 0.5 is anchored to Earth's own
400-700 nm window, so the solar reference has to stay there whatever window this
world's biosphere is given.

## Amendment: nitrogen, measured rather than assumed

This document listed nitrogen as "a declared constant" and the biases section
called deposition the largest unbracketed assumption in the biosphere. Measured
over 18 cells spanning 63 S to 57 N, that was wrong.

At the declared 0.5 kgN/ha/yr, deposition supplies 0.50 kgN/ha/yr against 4.36
from biological fixation and 25.4 from mineralisation: **1.7% of the nitrogen
plants actually receive.** Sweeping deposition over 150-fold, 0.1 to 15
kgN/ha/yr, moves NPP by 12% non-monotonically, which at `npatch 5` is patch
stochasticity rather than signal.

The nitrogen lever is `nfix_a`, the Cleveland fixation-against-evapotranspiration
fit, and LPJ-GUESS brackets it in `global.ins` itself:

| `nfix_a` | fixation | available N | NPP | vs central |
| --- | --- | --- | --- | --- |
| 0.102 conservative | 2.60 | 23.74 | 0.2356 | -10.8% |
| 0.234 central | 4.36 | 30.30 | 0.2642 | - |
| 0.367 upper | 5.20 | 29.85 | 0.2785 | +5.4% |

A monotone 18.2% span in NPP, all units kgN/ha/yr and kgC/m2 per simulation
year. That is inside the tolerance bands here, so again no prediction changes,
but productivity figures should be quoted with the bracket rather than from the
central value alone.

One thing that does not need correcting: the fit is a flux per unit
evapotranspiration, so it scales correctly with this world's shorter year
without any correction.

**Time-base correction, 2026-08-21:** that statement holds only for the
AET-proportional term. The Cleveland fit has a non-zero intercept applied once
per model orbit, and its five-year mean is five 181-day orbits rather than five
Earth years. The deposition input labelled per year is likewise delivered once
per orbit. BIO-24 corrects those semantics and BIO-2 must recompute this bracket;
the smoke-run values above remain provenance for the earlier experiment, not a
valid final nitrogen bracket.

## Result

Not yet run. To be filled in with actual values, hit or miss per line, and the
mechanism behind any miss.
