# The two shortwave absorptances, checked against a modern line list

Measured 2026-08-18, against the LMD Generic PCM's correlated-k tables. Both
terms hold. `h2osww` is confirmed to 1.4% and `co2sww` to 0.03%; the CO2
ABSOLUTE is confirmed to 8%, which is inside the check's own uncertainty. What
does not hold is the per-band attribution inside the CO2 term, and that is the
part of this worth reading.

`exoplasim/notes/shortwave-water-vapour.md` and `exoplasim/notes/shortwave-co2.md`
are the derivations being checked. `exoplasim/scripts/corrk_cross_check.py`
re-runs everything here; `--checks` runs the falsifiable checks alone.

## Why this is a test and not a second opinion

Both derivations come from Howard, Burch and Williams (1956), which is laboratory
band data from the 1950s reached through two closed-form fits, and nothing had
checked them against a line-by-line treatment. The tables here are HITRAN2020
through SpeCT and exo-k (Chaverot et al. 2025), premixed for N2 + CO2 + variable
H2O, and they bracket this planet directly at 376 and 1000 ppm CO2.

**The comparison is like-for-like because both sides compute the same defined
quantity**: the fraction of a star's TOTAL incident flux absorbed by one gas in a
HOMOGENEOUS path at 760 mm Hg holding a stated absorber amount. That is
Yamamoto's (1962) definition, it is what Lacis and Hansen Eq. 21 fits, and it is
what `shortwave_band_weights.py` reconstructs. The absorber amounts are the
derivation's own effective paths, so the pressure reduction and the `zbetta`
magnification are already inside both numbers and cannot separate them.

**The spectra are the same objects.** `corrk_cross_check.py` imports `blend` and
`Spectrum` from `shortwave_band_weights.py`, so the BT-Settl 4965 K and 5772 K
blends are the ones the derivation used. A systematic in the stellar model
therefore cannot appear here as a disagreement about absorption, and this check
is independent of the separate audit of the spectrum itself.

Only the absorption data differs: Howard's total band absorption over 8 CO2 and
9 H2O bands, divided by the band width for his Eq. 11 band-average, against
`1 - sum_g w_g exp(-k_g u)` over the correlated-k bands from 10 to 30000 cm-1,
76 of them once the IR and VI sets are joined -- see the join section below.

## What would have meant "wrong", stated before the numbers were computed

**PHYS-1, `h2osww`.** The weight is an absorption-weighted mean of the per-band
star-over-Sun flux boost, and that boost runs 0.99 at 0.72 um to 1.55 at 6.3 um,
so any correct answer lies inside that. The note's own prediction table gives
43.1 W/m2 of atmospheric shortwave absorption, 8.3 W/m2 at the top of the
atmosphere and 7.3 K per unit weight. Against that:

- within 0.035 of 1.346 confirms, because it moves the mean by under 0.25 K and
  lies inside the note's own 1.301 to 1.363 bracket;
- 0.035 to 0.10 is a real disagreement, 0.25 to 0.7 K, comparable to the
  "not included" convective term the note already flags at -0.2 to -0.5 K;
- beyond 0.10 the derivation is wrong by more than any other term in the
  shortwave budget and the flux re-derivation has to wait for it.

A direction was predicted too. HITRAN2020 resolves the 0.72 and 0.81 um bands
that Howard never measured and Yamamoto estimated from Fowle, and those sit where
the boost is LEAST. If the derivation's guess at them understates their strength,
the correlated-k weight comes out LOWER. Above 1.363 the note's bracket would not
have been a bracket.

**PHYS-6, `co2sww`.** The WEIGHT is nearly unfalsifiable in the sense that
matters: at 2.69 W/m2 the term carries 1.78 W/m2 of atmospheric absorption per
unit weight, 0.34 W/m2 at the top of the atmosphere, 0.29 K. A 0.15 error in
`co2sww` is 0.04 K. The thing to test is the ABSOLUTE, so:

- within a factor 1.3 of the derivation's absorptance confirms it;
- beyond a factor 2 the term is wrong by more than 1.35 W/m2 of atmospheric
  absorption and about 0.2 K, which is the size of effect the note treats as
  material elsewhere.

A direction was predicted here as well. Howard's 1956 set omits weak hot bands
and overtones a modern line list carries, and the derivation caps each band at
its own width, so the correlated-k absorptance should come out HIGHER. Coming out
substantially LOWER would have pointed at the derivation's construction, not at
the line list.

## The checks on the check

Three things had to hold before any of it counted, and they are what
`--checks` runs.

**The file layout is verified, not assumed.** CO2 is a trace N2-broadened gas at
both 376 and 1000 ppm, so its k must scale EXACTLY with the mixing ratio. In the
bands where CO2 dominates the ratio is 2.6595 against an expected 2.6596, over
356 points. That pins the array shape, the Fortran ordering, the units of k and
the meaning of the pressure grid all at once, and it would have failed loudly on
any of them.

**The H2O isolation is verified.** No table holds water without CO2, so H2O is
isolated by dividing the mixture transmission by the dry-slice transmission,
which is the random-overlap assumption and the same one the derivation makes when
it charges CO2 with what water vapour leaves. It has to return the same H2O
absorptance from the 376 ppm and the 1000 ppm table, because the water is
identical and only the CO2 differs by 2.66x. It does, to 0.05 to 0.20%.

**The spectral coverage is accounted for.** The tables span 10 to 30000 cm-1,
which holds 0.972 of the solar flux and 0.992 of this star's. The remainder is
below 0.33 um, where neither gas absorbs, and the denominator is the whole flux
on both sides, so nothing is lost. All of that 0.972 is inside a band: the
per-band fractions and the one-piece integral telescope exactly, which they can
only do for a partition with no hole in it, and `run_checks` raises on the
difference. The join that makes it a partition is the section below.

## The result

Correlated-k at 290 K and 1013.25 mbar, at the derivation's own effective paths:
225.56 atmos-cm of CO2 for this planet, 298.55 at Earth's gravity, 2.7891
precipitable cm of water.

| | derivation | correlated-k | difference |
| --- | ---: | ---: | ---: |
| `h2osww`, k25v over Sun | 1.3456 | 1.3271 | -1.4% |
| `h2osww` against the 4965 K blackbody | 1.2016 | 1.1943 | -0.6% |
| `co2sww`, k25v over Sun | 1.5098 | 1.5094 | -0.03% |
| `co2sww` against the 4965 K blackbody | 1.3640 | 1.3287 | -2.6% |
| H2O absorptance, solar-weighted | 0.14952 | 0.14914 | -0.25% |
| CO2 absorptance before the H2O overlap, planet path | 0.012190 | 0.013214 | +8.4% |
| CO2 absorptance after the H2O overlap, planet path | 0.005536 | 0.005098 | -7.9% |
| the Earth check, W/m2 | 2.05 | 1.92 | -6% |

**Both weights confirm, and `co2sww` confirms to four figures.** 1.3271 is 0.0185
from 1.346, well inside the 0.035 line and inside the note's own bracket. The
predicted direction was right: the correlated-k weight is the LOWER one, which is
what resolving the 0.72 and 0.81 um bands from a real line list rather than from
a scaled guess does to it.

**The Earth check passes on the correlated-k side too**, at 1.92 W/m2 against the
measured 1.5 to 2.5 W/m2, so the check the derivation set itself is not one that
only its own construction can pass.

The correlated-k CO2 absorptance carries about 10% of its own uncertainty, and it
is worth saying where from. The table's temperature nodes give 0.010360 at 230 K,
0.013214 at 290 K and 0.011887 at 350 K for the same column: the series is not
monotone, so the 290 K point sits on a wobble of roughly the same size as the
disagreement being measured. 290 K is the like-for-like point because Howard's
laboratory was at room temperature, but the 8% is not resolved by this check and
should not be quoted as though it were.

## The finding that is not a confirmation

The CO2 totals agree to 8% and the per-band attribution inside them does not.
Charging each CO2 band with the fraction of its interval water vapour leaves is
the right idea evaluated on the wrong resolution: Howard's band-mean absorptance,
spread uniformly across a 1000 cm-1 interval, smears saturation out of the band
cores and into the wings.

| CO2 band | derivation's clear fraction | correlated-k | contribution, derivation | correlated-k |
| --- | ---: | ---: | ---: | ---: |
| 15 um | 1.000 | 0.672 | 0.000168 | 0.000098 |
| 4.8 um | 0.632 | 0.561 | 0.000094 | 0.000129 |
| 4.3 um | 1.000 | 0.913 | 0.001784 | 0.001355 |
| 2.7 um | 0.171 | 0.003 | 0.000841 | 0.000016 |
| 2.0 um | 0.479 | 0.718 | 0.001931 | 0.002578 |
| 1.6 um | 0.941 | 0.922 | 0.000525 | 0.000233 |
| 1.4 um | 0.348 | 0.139 | 0.000190 | 0.000085 |

**Two errors of opposite sign, and they nearly cancel.** At 2.7 um the water
vapour band is black at this path and the correlated-k clear fraction is 0.003,
so the CO2 band there contributes essentially nothing. That is Yamamoto's own
verdict, and `shortwave-co2.md` decided against it. At 2.0 um the derivation
suppresses too hard, because Howard's 1.87 um water band average is applied
across an interval that is mostly window, and the correlated-k contribution is a
third larger. The two are worth -0.000825 and +0.000646, and the total survives.

**The note's third reason for keeping the 2.7 um band was sound and its first was
not.** "The decision does not decide the check" is correct: the check passes at
1.75 and at 2.05 W/m2 and the correlated-k answer of 1.92 sits between them, so
keeping the band is indeed not what makes the answer look right. But "the overlap
is COMPUTED band by band rather than argued about" claims a precision the band-
mean construction does not have. It is computed, and at 2.7 um it is computed
wrong by two orders of magnitude, and only the offsetting error at 2.0 um keeps
the total honest. This is `docs/src/practice/failure-modes.md` class 15 in a new place, and
the standing rule from it applies: do not correct one band without measuring the
chain.

The 15 um band behaves exactly as the note predicted it would. Its clear fraction
is 1.0 by omission, because Howard's near-infrared water set has no pure rotation
band, and the note says "the direction is to overstate it". Correlated-k puts it
at 0.672, overstated, by 0.000070 out of 0.005536.

Two things the derivation misses entirely and that partly refund the difference:
the correlated-k total over ALL bands is 0.005098 against 0.004494 over Howard's
eight intervals, so 13% of the CO2 shortwave absorption falls outside every band
Howard measured.

## The independent confirmation nobody asked for

`shortwave-water-vapour.md` reports that its reconstruction absorbs 11 to 16%
more than Lacis and Hansen Eq. 21 over the range that matters, median 13.4%, and
attributes it to Yamamoto having computed his curves for a real column with a
Curtis-Godson effective pressure while Eq. 21 is refitted as though it held at
standard pressure. That attribution was an argument, not a measurement.

It is now a measurement. Correlated-k against Eq. 21, on the same homogeneous
760 mm Hg path, runs 1.10 to 1.15 with 1.127 at the operating path:

| water path, cm | correlated-k over Eq. 21 |
| --- | ---: |
| 0.01 | 1.153 |
| 0.1 | 1.105 |
| 1.0 | 1.124 |
| 2.7891 | 1.127 |
| 5.0 | 1.132 |
| 10.0 | 1.142 |

**A modern line list absorbs 12 to 13% more than Eq. 21 at the same amount, which
is the derivation's excess to within a percent.** So the excess is a real deficit
in Eq. 21's pressure treatment rather than an error in the reconstruction, and
this could have failed: a ratio near 1.0 would have meant the reconstruction was
13% too absorbing and the whole H2O derivation suspect.

**That deficit is now corrected by `h2oswl` (PHYS-9), whose value 1.127 is
this table's operating-path ratio; it is not what `h2osww` is for.**
`h2osww` is a ratio and a level offset divides out of it. What the ratio says is
that the model's absolute clear-sky water vapour shortwave absorption is low by
about 12% relative to a modern line list, before the missing continuum, and that
is a separate defect from the one PHYS-1 fixes. It is not being fixed here.

## Making 1.127 two-sided: the bars, fixed before the sources were read

The ratio above is one-signed because the correlated-k side carries no water
vapour continuum, so it is a FLOOR on Eq. 21's deficit and the ceiling is
missing. Five sources are held for closing it and `references/INDEX.md` says
what each answers. These are the criteria the answer will be judged by, written
down before any of them was opened.

**What a unit of `h2oswl` is worth, so the bars can be in kelvin.**
`h2oswl` multiplies the same water vapour absorptance `h2osww` scales, so it
inherits that term's sensitivity exactly.
`exoplasim/notes/shortwave-water-vapour.md` measures 43.2 W/m2 of water vapour
shortwave absorption on the baseline climatology at `h2osww` = `h2oswl` = 1, and
prices `h2osww` = 1.346 at +15.0 W/m2 of atmospheric shortwave absorption, +2.9
W/m2 at the top of the atmosphere and +2.5 K. So d(atmospheric shortwave)/dL is
1.346 x 43.2 = 58.2 W/m2 per unit `h2oswl`; the top-of-atmosphere share is that
note's own 0.193, bracket 0.133 to 0.30; and the temperature follows at its
0.861 K per W/m2. **0.01 of `h2oswl` is 0.10 K, bracket 0.07 to 0.15 K.**

**BAR 1, does the key move.** `h2oswl` changes from 1.127 if the
continuum-inclusive CENTRAL estimate differs from it by 0.02 or more, which is
0.2 K: the size of effect this note already treats as material for the CO2
absolute, and the size of the convective term `shortwave-water-vapour.md` lists
as not included. Under 0.02 the two-sided bound is recorded and the key keeps
1.127, because a config key change makes every existing run unresumable and
0.1 K does not buy that.

**BAR 2, is the bracket carried.** Same conversion. The bound is quoted as a
bracket on the key rather than as a point value if its half-width reaches 0.02,
0.2 K. Narrower than that and it is a number with a stated uncertainty; wider,
and every consumer has to run both arms.

**BAR 3, does `shortwave_band_weights.py` keep raising.** That file gates on the
Howard reconstruction lying inside the envelope of Lacis and Hansen Eqs. 21, 22
and 23 widened by Howard's +/-3%, and the gate currently fires. The gate asks
whether the reconstruction agrees with every published determination of the
quantity. It is FALSIFIED AS A GATE -- not widened, which would be a criterion
chosen after the run -- if an independent MODERN reference for the same defined
quantity falls outside the same widened envelope at the same water amounts,
because a bar that rejects a known-good answer is not a bar on the
reconstruction. The correlated-k absorptance in the table above is that
reference and the test is arithmetic against the three formulas. If correlated-k
lands INSIDE the envelope where the reconstruction lands outside, the gate
stands and the reconstruction is what is at fault.

**BAR 4, are the two checks one quantity.** They are the same quantity if and
only if numerator and denominator are the same defined thing on both sides:
the absorptance of a homogeneous 760 mm Hg path holding w precipitable cm of
water, solar-weighted, over Eq. 21 at the same w. If they are, they must agree
at the operating path to within the +/-3% Howard state for the band absorptions,
that being the only stated accuracy either side carries; a disagreement wider
than 3% at 2.7891 cm is a defect in one of them and not a difference of method.

## What the sources gave, and the number that comes out

**The formulation does not ship with its paper.** `mlawer2012` gives the MT_CKD
line shape and nothing that can be evaluated: the continuum coefficients are the
sum of that shape over every water vapour line, so they are a data file that
travels with LBLRTM, not a table in the article. Its Table 3 and Sect. 4 are
about REVISIONS to those coefficients rather than the coefficients themselves.
So the term is bounded from `shine2012`, which reports its size, and
`mlawer2012` is read for how well it is known.

**Both sides of the ratio are missing the same part of the continuum, and it is
the part in the windows.** The correlated-k tables carry line centres to +/- 25
cm-1 with the plinth removed, which is the MT_CKD definition, so the whole
MT_CKD-convention continuum is absent from them. Eq. 21 is a fit to Yamamoto,
who summed Howard's laboratory BAND absorptions, so it has no absorption at all
between the bands. The two therefore lack the same thing, and the agreement in
the section below is evidence about which part is missing: within the bands the
continuum is secondary because the lines are already strong (`shine2012` p. 536),
and Howard's cell measurements carry what there is of it, which is why a Howard
reconstruction and a modern line list land 0.26% apart at the operating path.
What neither has is the window continuum, and that is exactly the term
`shine2012` measures.

**The size.** `shine2012` p. 548: the CAVIAR laboratory continuum "adds about
1.25 W m-2 (or about 1.5%) to the global-mean clear-sky shortwave absorption,
relative to that calculated using the MT-CKD model". That fixes clear-sky
atmospheric shortwave absorption at 83 W/m2, and water vapour is 60% of it by
`shine2012` p. 536 and 72% of it on this world's own baseline
(`exoplasim/notes/shortwave-water-vapour.md`), so the CAVIAR increment is 2.1 to
2.5% of water vapour's shortwave absorption.

**The total is what is missing, not the increment**, because the correlated-k
side has no continuum of any strength. Write f for the factor by which the
laboratory determination exceeds MT_CKD in the near-infrared windows. Then the
increment is (f - 1) times MT_CKD's own contribution and the total is f/(f - 1)
times the increment, which is the quantity that survives not knowing f:
`shine2012` p. 535 and p. 547 put f at "up to around an order of magnitude",
and `mlawer2012` p. 2551 gives the laboratory factors one window at a time, 6-12
at 4600 cm-1, 2-4 at 6140 cm-1, 2-3 at 9466 cm-1 and 2-4 at 2400-2640 cm-1. Over
f from 2 to 12 the total runs 2.3 to 5.0% of water vapour shortwave absorption,
and the middle of the laboratory factors, f = 4, gives 3.1%.

| f, laboratory over MT_CKD in the near-infrared windows | MT_CKD's own share | the total continuum |
| --- | ---: | ---: |
| 2 | 2.1-2.5% | 4.2-5.0% |
| 4 | 0.70-0.83% | 2.8-3.3% |
| 12 | 0.19-0.23% | 2.3-2.7% |

**Transferring a global mean to this path is the weak link and is bracketed
rather than argued away.** The quantity here is one homogeneous 760 mm Hg path
at 290 K holding 2.7891 cm at vmr 1e-2, so a water partial pressure of 10.1 mb;
an Earth column's water-mass-weighted mean water partial pressure is half its
surface value for any power-law profile, 6 to 12 mb, and its mass-weighted total
pressure is about 810 mb against this path's 1013. The self continuum follows the
water partial pressure and the foreign continuum the total, so this path collects
MORE continuum per unit water than a column does; the band absorption it is
divided by rises as well, as the square root of pressure for a strong band, and
this path is 10 K warmer than a column's water-weighted mean, which cuts the self
continuum back. Net, the fraction transfers within about 0.8 to 1.4 and the
direction is toward more.

**So the missing term is 0.15 to 7.0% of the correlated-k absorptance, and 3.1%
is the estimate.** The floor is MT_CKD taken as correct at the largest f, which
is the position that disbelieves the CAVIAR laboratory work entirely; the ceiling
is the laboratory determination at the smallest f with the path transfer at its
largest.

### The base is corroborated by a fourth method, and the family is one-signed

The three remaining sources are about what a broadband parameterization of this
kind is scored against, and one of them turns out to score THIS one.

**`ramaswamy1992` puts Lacis and Hansen (1974) against a line-by-line
calculation directly.** Its Table I, p. 11,491, gives the solar flux absorbed by
water vapour in a mid-latitude summer column at surface albedo 0.2, in W/m2:

| | 30 deg, direct | 30 deg, total | 75 deg, direct | 75 deg, total |
| --- | ---: | ---: | ---: | ---: |
| line-by-line | 168.4 | 178.1 | 68.2 | 71.4 |
| Lacis and Hansen | 154.3 | 162.3 | 62.5 | 63.6 |

which is line-by-line over Lacis and Hansen of 1.091, 1.097, 1.091 and 1.123, and
the paper's own verdict on p. 11,510 is that "both methods yield an underestimate
in absorption". **That is a fourth construction of the same defect, and its
reference has no continuum either**: the word does not appear in the paper and the
line data is the AFGL 1982 catalogue. So four methods now put Eq. 21 low by a
continuum-free 9 to 13%, and they were built from different absorption data on
different decades' understanding.

It is NOT folded into the bound, because it is a different quantity: a real
column with a slant path and a reflecting surface rather than the homogeneous
760 mm Hg path this note and `shortwave_band_weights.py` both define. It
corroborates the sign and the size and it is not a fifth measurement of the
number.

**`oreopoulos2012` says the whole family is one-signed against a reference that
DOES carry the continuum.** CIRC Phase I scores 13 shortwave codes against
LBLRTM/CHARTS with MT_CKD_2.0, and p. 6 reports that "Only Models 2 and 5 (an
LBL code) achieve absorption errors within 2.5% for all the pristine cases"; the
conclusions on p. 15 say "Previously found underestimates of SW absorption by
less spectrally detailed models [Ackerman et al., 2003] were confirmed again
here". Every
code's mean column-absorption error in its Figure 4 is negative. So a 12 to 16%
under-absorption for a two-band 1974 scheme is at the outer edge of the
documented behaviour of this family and not outside it.

**`freidenreich1999` sets the floor on what a good scheme achieves.** Its
72-interval successor is accurate to within 2% on clear-sky column absorption,
Table 7 p. 31,404, and the sign never flips across nine profile and zenith-angle
combinations: -1 to -2% everywhere. It does NOT score original Lacis and Hansen
-- its "old" water vapour scheme is `ramaswamy1992`'s already-corrected modified
LH -- so nothing here is a second measurement of Eq. 21.

### BAR 1: the key moves, to 1.163

    correlated-k over Eq. 21, no continuum          1.127
    with the continuum, estimate                   1.163
    two-sided bound                        1.129 to 1.206

1.163 - 1.127 = 0.036, which is 0.35 K and clears the 0.02 the bar was set at,
so `h2o_sw_level` becomes **1.163** and stops being one-signed. The third digit is
carried so the arithmetic above reproduces, not because it is resolved; the
bracket is what says how well the number is known. The direction is
the one that was declared: the continuum ADDS absorption to the line-by-line
side, so it moves Eq. 21 further below a modern calculation, and 1.127 was a
floor.

The cost side of that bar is small in fact rather than in principle. A config key
change makes every existing run unresumable, but no run has yet used the k25v
spectrum or turned `H2OSWL` on -- `radmod.f90`'s default is 1.0 and every
existing run carries it -- so what the change invalidates is the SETTING and not
a result. Under rule 7 nothing is made worthless that was not already disposable:
there is no canonical climatology lineage.

### BAR 2: it is carried as a bracket

The half-width is (1.206 - 1.129)/2 = 0.039, which is 0.38 K and clears the 0.02
the bar was set at, so the key is quoted with `h2o_sw_level_bracket` and the two
ends are arms to run rather than a stated uncertainty on a settled number. The
width is dominated by not knowing f and by the transfer to this path, in that
order, and neither is resolved by anything this project can compute: closing it
needs a line-by-line calculation with an MT_CKD or CAVIAR continuum on this
path, which is a correlated-k bundle this project does not have.

**What the continuum does to `h2osww` is separately negligible and stays at
1.346.** The continuum absorbs in the windows, where this star's flux boost is
largest, so it raises the star-over-Sun ratio as well: at 3.1% of the absorption
and a local boost of 1.45 the correlated-k weight goes from 1.3271 to about
1.331, which is 0.004 and 0.03 K. That is inside the 0.035 line PHYS-1 fixed for
itself and far inside the 1.301 to 1.363 bracket, so the weight is untouched.

### BAR 3: the envelope gate in `shortwave_band_weights.py` is falsified, and stops raising

The test was whether an independent modern reference for the same defined
quantity also falls outside the envelope of Eqs. 21, 22 and 23 widened by
Howard's +/-3%. It does, at every amount but the driest:

| water path, cm | correlated-k over Eq. 21 | envelope top over Eq. 21 | correlated-k against the envelope |
| --- | ---: | ---: | --- |
| 0.01 | 1.153 | 1.506 | inside |
| 0.1 | 1.105 | 1.039 | outside by 6.4% |
| 1.0 | 1.124 | 1.030 | outside by 9.1% |
| 2.7891 | 1.127 | 1.030 | outside by 9.4% |
| 5.0 | 1.132 | 1.031 | outside by 9.8% |
| 10.0 | 1.142 | 1.076 | outside by 6.1% |

That is the same pattern the Howard reconstruction shows, inside only at 0.01 cm
and outside by 6 to 10% elsewhere, and it is a FLOOR because the correlated-k
side has no continuum: adding it pushes HITRAN2020 further out still. **A bar
that rejects a modern line list is not a bar on the reconstruction.** The
envelope is three determinations from 1915, 1956 and 1962, two of which trace to
the same laboratory data, and what has been established is that all three are low
rather than that the reconstruction is high.

So the check stops raising and becomes what it can support: a record of how far
every construction of this quantity, old and new, sits from the published
envelope. It is NOT widened -- widening a bound to admit the thing that missed it
is the move `world-zvk2` refused and this does not do it. The gate that replaces
it is BAR 4's, which has a right answer that can fail.

### BAR 4: the two checks are one quantity, and are now tied together

`shortwave_band_weights.py`'s `ratio_to_eq21` and this note's correlated-k ratio
have the same numerator definition and the same denominator: the solar-weighted
fraction of total incident flux absorbed by water vapour in a homogeneous 760 mm
Hg path holding w precipitable cm, over Eq. 21 at the same w. Only the absorption
data differs, Howard's nine bands against 76 correlated-k bands on HITRAN2020.
They are the SAME QUANTITY measured two ways, and at the operating path they
agree to 0.26%, far inside the +/-3% the bar was set at:

| water path, cm | Howard reconstruction over Eq. 21 | correlated-k over Eq. 21 | apart |
| --- | ---: | ---: | ---: |
| 0.01 | 1.319 | 1.153 | 14.4% |
| 0.1 | 1.109 | 1.105 | 0.4% |
| 1.0 | 1.126 | 1.124 | 0.2% |
| 2.7891 | 1.130 | 1.127 | 0.26% |
| 5.0 | 1.137 | 1.132 | 0.4% |
| 10.0 | 1.159 | 1.142 | 1.5% |

**The dry end is where they part, and that is a real finding rather than
scatter.** At 0.01 cm the reconstruction is 14% above the correlated-k answer,
which is the small-w Elsasser extrapolation and the 0.72 and 0.81 um bands that
Howard never measured, carried at the 0.94 um band's shape under a scale factor
no source supplies. It is also the one amount at which both land inside the
envelope, so the two failures do not overlap. The model does not operate there --
a T42 column spans roughly 0.3 to 5 cm -- and the weight is a ratio in which the
dry end largely divides out, so this is recorded and tracked and does not gate.

The two numbers were free to drift apart because nothing compared them.
`shortwave_band_weights.py` now carries the correlated-k ratio at the operating
path as a declared constant and raises if its own reconstruction disagrees with
it by more than Howard's +/-3%, which is the only stated accuracy either side
has. A change to the band set, to the spectra or to this note's table now fires
in one place instead of leaving two numbers for one quantity.

### The one thing the model has no term for at all

`radmod.f90` carries a water vapour continuum coefficient, `th2oc`, and it is in
`lwr` and only in `lwr`. `swr` has ozone in band 1 and water vapour in band 2 and
nothing else, so there is no shortwave continuum term to set and the continuum
reaches the shortwave only through `h2oswl`. That is why the correction is a
level on Eq. 21 rather than a new absorber: a term with the continuum's spectral
shape would absorb in the windows, which this two-band scheme cannot resolve.

## What the disagreements would move if they are real

Using the sensitivities `shortwave-water-vapour.md` states for itself:

| | shift | atmospheric SW | top of atmosphere | mean surface temperature |
| --- | ---: | ---: | ---: | ---: |
| `h2osww` 1.3456 to 1.3271 | -0.0185 | -0.80 W/m2 | -0.15 W/m2 | -0.135 K |
| CO2 absorptance, -7.9% | | -0.21 W/m2 | -0.04 W/m2 | -0.035 K |
| both | | -1.01 W/m2 | -0.19 W/m2 | **-0.17 K** |

Against a combined prediction of +2.94 K. That is smaller than the note's own
top-of-atmosphere bracket, smaller than the -0.2 to -0.5 K convective term it
lists as not included, and far smaller than the 3.7 to 7.1 K endmember spread the
flux window is derived inside. **Nothing here changes the flux re-derivation, and
nothing here is a reason to change either term.** Physics is not a knob; a 1.4%
disagreement between two independent absorption datasets is not evidence that
either is wrong.

## Every choice this comparison made, and what it was worth

The point of listing these is that each is a place the comparison could have
become two different quantities disagreeing about neither term.

| choice | made as | sensitivity |
| --- | --- | --- |
| temperature | 290 K, because Howard's laboratory was at room temperature | `h2osww` 1.3242 at 230 K to 1.3281 at 320 K, so 0.3%. CO2 absorptance +28% from 230 to 290 K and the nodes are not monotone, so 10% |
| water broadening fraction | vmr 1e-2, the mid-column value | `h2osww` 1.3244 at 1e-1 to 1.3277 at 1e-3, so 0.25% |
| water path | 2.7891 cm, the derivation's own magnified path | this is the LARGE one: `h2osww` 1.4435 at 0.01 cm to 1.2937 at 10 cm. The derivation reports the same shape, 1.416 to 1.313, so both agree a constant weight is worth about 7% across the range the model spans |
| pressure | 1013.25 mbar homogeneous, both sides | not scanned; it is the definition of the quantity, not a free parameter |
| isolating H2O from a premixed table | transmission ratio, random overlap | 0.05 to 0.20%, measured against the 1000 ppm table |
| CO2 mixing ratio of the table | 376 ppm against the config's 450 | none: the absorber amount is set independently and k is verified to scale exactly with the mixing ratio |

**And the one thing that could not be made comparable.** The tables carry line
centres only, plus or minus 25 cm-1, with the plinth removed on the MT_CKD
convention, and the bundle's `continuum/far_wing_data` has files for CO2-CO2,
CO2-H2O and CO2-N2 but none for water with itself or with N2. So the correlated-k
water vapour absorptance has NO self or foreign continuum in it. The continuum
absorbs in the windows BETWEEN the bands, which is where this star's flux boost is
largest, so including it would raise the correlated-k weight rather than lower it:
if it added 5% of the absorption at a local boost of 1.45 the weight would go to
about 1.333. **The -0.0185 is therefore an upper bound on the disagreement, not a
central estimate**, and the true gap is smaller than the number in the table.

## Where the weight comes from

Absorption-weighted, so this is the decomposition of 1.3271 rather than of the
flux. It is here because the derivation's equivalent table is per Howard band and
this one is per wavelength, and the two say the same thing.

| wavelength | share of the solar-weighted H2O absorptance | local star-over-Sun boost |
| --- | ---: | ---: |
| 0.30-0.75 um | 3.3% | 0.958 |
| 0.75-1.00 um | 17.9% | 1.098 |
| 1.00-1.30 um | 16.5% | 1.210 |
| 1.30-1.60 um | 27.9% | 1.352 |
| 1.60-2.20 um | 16.3% | 1.518 |
| 2.20-3.50 um | 14.7% | 1.515 |
| 3.50-6.00 um | 1.9% | 1.513 |
| beyond 6 um | 1.4% | 1.546 |

The reason `h2osww` exceeds the 1.204 flux-share ratio is visible in one column:
the absorption sits at 1.3 to 3.5 um, where the boost is 1.35 to 1.52, and only a
fifth of it sits below 1 um where the boost is near unity. The derivation reaches
the same conclusion from Howard's bands and the agreement is on the mechanism, not
only on the number.

## The join between the IR and VI sets, and the hole it used to leave

Re-measured 2026-08-24; the finding is world-olt. Everything above this section
was measured on 2026-08-18 against a 75-band join and its digits are that day's.

**The bundle does not ship one band set, and it does not ship a gap either.**
The `40x38` grid is two sets that OVERLAP: `narrowbands_IR.in` runs 10 to 3000
cm-1 in 40 bands and `narrowbands_VI.in` runs 2000 to 30000 in 38, each
contiguous inside itself to the bit, and they share 2000 to 3000. The Generic
PCM never concatenates them -- `rad_correlatedk_read_opacity_tables.F90` sets
`IR_VI_wnlimit = 3000.` and hands `WNOI` and `WNOV` to two independent solvers,
one thermal and one stellar. The two files are byte-identical between the 376
and 1000 ppm tables, so this is a property of the grid and not of a mixture.

**One flux-weighted partition is what THIS check needs, so the join is this
script's own decision, and the first version of it opened a hole.** Cutting the
IR set at a round 2000 cm-1 and dropping every band that reached past it also
dropped 1974.95 to 2000, because the IR grid has no edge at 2000: its nearest
edge below is 1974.952011. That 25.05 cm-1 window, at 5.0 to 5.06 um where the
shortwave and thermal halves hand over, was inside the table span and inside no
band. It carries 1.7e-4 of the solar flux and 2.5e-4 of this star's, and every
band-weighted total priced it as transparent in both gas sets.

**The join is now the VI set's own first edge, and the IR band that straddles it
is truncated rather than dropped.** The two halves meet to the bit, the
truncated band keeps its own k-distribution -- the only opacity the bundle
carries over that sliver -- and no edge moves outward. That last part is the
constraint: the per-band flux fractions are differences of one cumulative
integral, so they telescope exactly for a contiguous partition and cannot for
anything else, and widening a band to make the sum close would have removed the
only instrument that can see a hole. `run_checks` still counts holes and still
raises on one, and now reports the join and the truncation by name.

**What closing it moved**, at the same paths and the same spectra as the table
above:

| | 75-band join | 76-band join |
| --- | ---: | ---: |
| flux fraction inside the span, in bands | 0.97187 | 0.97204 |
| flux fraction inside the span, in no band | 0.00017 | 0.00000 |
| `h2osww`, k25v over Sun | 1.3271 | 1.3272 |
| `co2sww`, k25v over Sun | 1.5094 | 1.5093 |
| CO2 absorptance after the H2O overlap, planet path | 0.005098 | 0.005123 |
| the same over Howard's eight intervals only | 0.004494 | 0.004520 |

The two weights are unmoved at the digit that matters, because both are RATIOS
and the sliver enters numerator and denominator alike. The absorptance moves by
+0.5%, worth +0.01 W/m2 of Earth-mean insolation.

**The refit was re-landed.** `--fit` reads `radmod.f90`'s four coefficients back
out of the source and says whether they still match, so the join change showed
up there as a disagreement rather than as nothing. The coefficients follow their
derivation: `zca1` 3.1020E-4, `zcb1` 19.857, `zca2` 3.8291E-3, `zcb2` 3.9587E-3,
against 3.0658E-4, 20.376, 3.8193E-3 and 3.9672E-3. The absorptance at this
planet's path goes from 0.005027 to 0.005051, +0.5% and +0.008 W/m2, against a
fit residual of 3.3% rms -- so this is bookkeeping following a derivation, and
not a result worth chasing. The fit quality is unchanged to the digit that
matters: 8.8% at u = 1e4 and 4.5% over the 100 to 1000 atmos-cm a T42 column
occupies, against 8.9% and 4.6%.

**What a consumer has to carry forward.** A band-resolved scheme built on this
bundle should NOT inherit this script's join. The two sets are meant to be used
separately, over their own spans, by two solvers; joining them is an instrument
for comparing one broadband number against another, and the 2000 to 3000 cm-1
overlap the join throws away is real data that a thermal scheme wants.

## The data

`~/git/generic_pcm/LMDZ.GENERIC/datagcm/corrk_data/N2-0.000376CO2-H2Ovar_2026` and
`N2-0.001CO2-H2Ovar_2026`, `40x38` bands. The format was read from
`LMDZ.GENERIC/libf/phygeneric/rad_correlatedk_read_opacity_tables.F90` and
`rad_correlatedk_opacities_stellar.F90` rather than guessed. That tree is CeCILL
licensed and this project is not: it is a design reference, no code crossed, and
`docs/src/reference/external-data.md` records the route.

Chaverot et al. (2025), A&A, https://doi.org/10.1051/0004-6361/202555762 for the tables and the SpeCT line-by-line code behind them.
