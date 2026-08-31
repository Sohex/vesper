# The phosphorus cycle's parameters, and which of them are phosphorus

This is a worldbuilding project. Vesper is a simulated super-Earth around a
mid-K dwarf, and its biosphere is LPJ-GUESS-CNP, a vegetation model written for
Earth and ported to this world's calendar. Everything below is about that
model's phosphorus constants and the simulated soil and plant processes they
drive.

The C-N-P fork's phosphorus side was built by copying its nitrogen side. Each
copy is defensible on its own and together they mean the modelled phosphorus
cycle has little parameterisation of its own: uptake competition ranks the way
nitrogen competition ranks, leaf N:P is pinned, and one of the two soil
saturation thresholds is nitrogen's. The other is not a copy at all, and
reading the paper it comes from moved the defect from the constant to the pool
the constant is applied to. This document says, for every one of those
constants, what it does to the modelled result, what the measured phosphorus
value is or what its bracket is, and which of them are still nitrogen's.

Nothing here is verified by execution. LPJ-GUESS builds and runs on this tree,
and no run it has made has passed acceptance, so no simulated phosphorus stands
behind anything below. Every figure is hand arithmetic off the source and the
cited measurements.

## The register

| constant | where it lives | read at | value | unit | source or bracket | moved |
| --- | --- | --- | --- | --- | --- | --- |
| `PUPS_UPPER_ADV` | `guess.h` | `Pft::init_pupscoeff`, then `canexch.cpp:phosphorus_uptake_strength` | 3.18 | dimensionless | Jobbagy and Jackson (2001) profile contrast on the model's own nitrogen value; BRACKET 2.9 to 3.5 | yes, from 2.0 |
| `PFRAC_MINTOMAX` | `guess.h` | `Pft::init_ctop_limits` and `Pft::init_ctop_min` | 3.68 | dimensionless | McGroddy et al. (2004) foliar dispersion contrast; BRACKET 3.68 to 4.41 | yes, from 2.78 |
| `PFRAC_MINTOMAX_CROPGREEN` | `guess.h` | `Pft::init_ctop_limits` | 7.77 | dimensionless | the same transfer applied to its nitrogen counterpart 5.0; cropland is refused under BIO-27 | yes, from 5.0 |
| `PFRAC_LEAFTOROOT` | `guess.h` | `Pft::init_ctop_limits`, then `canexch.cpp` fine-root P demand | 1.16 | dimensionless | Yuan et al. (2011); root N:P is not separable from leaf N:P, so the nitrogen proportion carries; BRACKET 1.02 to 1.35. The tissue window is mean-anchored, so the contrast the model runs on is the derived 1.0 for any leaf window width | no, and now for a reason |
| `PFRAC_LEAFTOSAP` | `guess.h` | `Pft::init_ctop_limits`, then `canexch.cpp` sapwood P demand | 6.9 | dimensionless | UNDERIVED, on the LEVEL and not the form. Nitrogen's, for sapwood PLUS BARK, from Friend et al. (1997) Table 4 p. 254. The proportional form stands: Yan et al. (2016) put the phosphorus exponent at 1.58 tropical, 0.97 temperate and 0.80 boreal over 335 species and 12 sites, so it crosses 1 inside the climate range a global model spans. The 6.9 sits below the 10.1 to 15.5 the one paired leaf-and-bole-wood phosphorus dataset brackets a scalar at | no |
| `PFRAC_MAXTOMIN` | `guess.h` | `Pft::init_ctop_limits` | 0.9 | dimensionless | nothing to transfer: the nitrogen original is declared arbitrary in `Pft::init_cton_limits`. It sets the tissue window's WIDTH only, so it no longer moves the applied proportion | no |
| `PMASS_SAT` | `somdynam.cpp` | `somfluxes` through `setptoc`, and the P-limitation-off pin | 0.002 * 6.6 | kgP/m2 labile P | Parton, Stewart and Cole (1988) Fig. 3 p. 115, whose labile-P axis saturates at 2.0 gP/m2, CONVERTED into this fork's Hedley-labile currency by 6.6, the low end of the 6.6 to 11.3 bracket | yes, and it is a declared divergence from the vendored 0.002 |
| `PCONC_SAT` | was `somdynam.cpp` | nowhere | removed | -- | the ramp it was the threshold of is removed: a decomposer community's biomass C:P is homeostatic with respect to its resource's phosphorus content, so the surface microbial pool's C:P does not vary with litter P. The surface microbial pool holds the C:P `soil.cpp` initialises it to, which is `world-634q` | deleted |
| `USORB` | `somdynam.cpp` | `somfluxes` | 0.0067 / `VESPER_EARTH_YEAR_DAYS` | per absolute day | Wang et al. (2010) Appendix D | divisor, under the time-base contract |
| `USSORB` | `somdynam.cpp` | `somfluxes` | 0.0067 / `VESPER_EARTH_YEAR_DAYS` | per absolute day | the same, and equal to `USORB` in that source | divisor, under the time-base contract |
| `UOCC` | was `somdynam.cpp` | nowhere | removed | -- | no citation anywhere in the tree | deleted |
| `Soiltype::pwtr` | the soil input | `somdynam.cpp:soilpadd` | per BIO-5 | kgP/m2 per EARTH year | BIO-5 emits it against `biosphere/notes/time-base-unit-contract.md` | divisor, under that contract |

The three declared divergences from the vendored CNP fork are registered in
`biosphere/config/somdynam.yaml` and checked by
`biosphere/scripts/somdynam_gate.py`, which refuses a form the source stops
recording, a form the source goes back to running, a changed line the source
loses, a ramp whose `setptoc` call moves, and the phosphorus-limitation-off pin
being split off the ramp threshold. Its reference point is the vendored subtree
commit rather than a release, because LPJ-GUESS 4.1.1 has no phosphorus.

Every phosphorus quantity in `somdynam.cpp` is reached only through
`som_dynamics_century`, the live `ifcentury 1` path. `som_dynamics_lpj`, the
`ifcentury 0` path, carries no phosphorus at all, so the two decay paths do not
have to be told apart for any constant in this document.

## The uptake profile: the contrast transfers, the level does not

`init_pupscoeff` builds `pupscoeff = rootdist_upper * PUPS_UPPER_ADV +
rootdist_lower`, where the upper layer is 0 to 500 mm and the lower is 500 to
1500 mm. Its nitrogen counterpart is 2.0, attributed to Franzluebbers et al.
(2009) as an approximate factor of two.

Jobbagy and Jackson (2001) measure the fraction of the top metre's content held
in the top 20 cm across the USDA National Soil Characterization Database: a
median of 0.43 for extractable P against 0.36 for total N, with phosphorus the
shallower of the two in 61 percent of paired profiles. Extractable P is the
closest measured analogue this model's labile P pool has, which is what makes
the pair usable at all.

The absolute construction on those factors is refuted, not untested. Fitting a
density profile that declines exponentially with depth to each measured factor
and taking the mean density in 0 to 500 mm against 500 to 1500 mm returns 3.45
for nitrogen, where the model carries 2.0. So the level does not transfer. The
same construction returns 5.47 for phosphorus, a contrast of 1.588, and 2.0
times that contrast is the 3.18 adopted.

The bracket is the depth Jobbagy could not measure. The data stop at 1 m and the
model's lower layer reaches 1500 mm, so the contrast is recomputed with the
profile truncated at 1000 mm, giving 1.446, and extended to 3000 mm, giving
1.740. That is `PUPS_UPPER_ADV` between 2.89 and 3.48.

What it is worth: `pupscoeff` enters only through
`phosphorus_uptake_strength`, and `pcompete` distributes the available labile P
in proportion to that strength, so a uniform scaling cancels and only the
contrast between plant functional types survives. `global.ins` gives grasses
0.18 per layer in the upper five layers, so `rootdist_upper` is 0.9, and trees
0.12, so theirs is 0.6. Strength goes as `pupscoeff` to the two-thirds power.
The shallow-rooted type's share of labile P per unit fine-root mass therefore
rises from 12.1 percent above the deep-rooted type's to 18.1 percent above it.
That is the whole of the effect, and it is a competition effect and not a
productivity one.

## The leaf C:P window: phosphorus moves further than nitrogen

`PFRAC_MINTOMAX` is the ratio of maximum to minimum leaf C:P a plant functional
type may occupy. Its nitrogen counterpart is 2.78 from White et al. (2000) and
Reich et al. (1992), which measured nitrogen.

McGroddy et al. (2004) report the coefficient of variation of foliar ratios
across forests worldwide over the same 55 to 59 stands: 79 percent for C:P
against 59 percent for C:N. Treating each ratio as lognormal, the log standard
deviations are sqrt(ln(1 + 0.79^2)) = 0.696 and sqrt(ln(1 + 0.59^2)) = 0.547,
in the ratio 1.274. A window of width 2.78 widened by that exponent is
2.78^1.274 = 3.68.

The bracket's other end applies the same dispersion contrast to the window's
relative half-width instead of to its logarithm. The half-width (f - 1)/(f + 1)
is 0.4709 at f = 2.78; scaled by 0.79/0.59 it is 0.6305, which is f = 4.41. The
lower end is adopted.

Both ends rest on one assumption, stated here because it is not measured: that
within-type physiological plasticity in leaf C:P scales like across-forest
dispersion. The direction does not depend on it. Leaf N:P is not fixed across
species or across soil P supply, which is the effect a C-N-P model exists to
resolve and which a copied 2.78 removes.

`PFRAC_MINTOMAX_CROPGREEN` is its nitrogen counterpart 5.0 through the same
exponent, 7.77. Cropland is inert on the natural-vegetation baseline and fails
closed under BIO-27, so it is converted for consistency and not for a result.

## The mean convention in `init_ctop_min`, which is a defect on its own

The regression that sets leaf C:P is `exp(8.63342 + log(sla) * -0.80936)`, and
Dantas de Paula et al. (2025) Eq. A3 describes its output as the AVERAGE leaf
C:P, which is what plant P demand is computed from.

The code divided that average by `(f + 1) / 2` to obtain the minimum, which is
the arithmetic-mean relation. But `avg_ctop` returns the HARMONIC mean,
`2 / (1/min + 1/max)`, which is the C:P whose P:C is the arithmetic mean of the
two bounding P:C ratios and is the right average for a tissue concentration. The
two conventions do not agree, so `ctop_leaf_avr` came out at
`4f / (1 + f)^2` of the regression's own output: 0.778 of it at f = 2.78.

The nitrogen side has no such defect, because its regression returns the minimum
leaf C:N directly and no divisor is applied.

The fix is to multiply by `(1 + f) / (2f)` instead, which makes
`ctop_leaf_avr` equal the regression's output exactly for any f. That is a check
that can fail, and it is the one the arithmetic form failed.

What it is worth, alone and independently of every other change here: leaf
phosphorus demand per unit leaf carbon falls by 22.2 percent, for every plant
functional type, exactly. It is the only change in this document whose size does
not depend on where a run lands.

Taken together with the wider window, and writing A for the regression's output:
`ctop_leaf_min` goes from 0.529 A to 0.636 A, so the ceiling on leaf P content
falls 16.8 percent; `ctop_leaf_max` goes from 1.471 A to 2.340 A, so the floor on
leaf P content falls 37.1 percent; and `ctop_root_avr` and `ctop_sap_avr` are
proportional to `ctop_leaf_avr`, which is A itself, so they follow the
regression's own output and not either endpoint.

## Fine roots: the copy is right, and now it is a result

`PFRAC_LEAFTOROOT` multiplies leaf C:P to give fine-root C:P. Dividing a C:P
proportion by its C:N counterpart gives exactly the ratio of root N:P to leaf
N:P, so that ratio is the quantity to measure.

Yuan et al. (2011) compile 631 live fine-root N:P samples and test them against
the two global green-leaf compilations. They cannot separate them: p = 0.271
against Reich and Oleksyn (2004) and p = 0.120 against Wright et al. (2004).
Fine roots are P-poor and N-poor against leaves in the same proportion, which is
the paper's own headline. So the contrast is 1.0 and the nitrogen value carries
unchanged.

The bracket is the point estimates the test could not separate: live-root N:P of
16.0 against leaf N:P of 13.8 and 18.2, giving contrasts of 1.16 and 0.879, so
`PFRAC_LEAFTOROOT` between 1.02 and 1.35. The value in the source is inside it.

**The model delivers the contrast, and it is decoupled from `PFRAC_MINTOMAX`.**
`Pft::init_ctop_limits` sets the fine-root window's MEAN from the leaf window's
mean and builds the window outward from it, and `Pft::init_cton_limits` does the
same on the nitrogen side, so what reaches simulated fine-root phosphorus demand
is `ctop_root_avr / ctop_leaf_avr` = 1.16 exactly and its nitrogen counterpart is
1.16 exactly. The contrast is the ratio of the two constants, 1.0, whatever
widths the two leaf windows carry. That is why widening the leaf C:P window to
3.68 while nitrogen's stayed at 2.78 does not touch this constant.

The earlier max-anchored form applied the constants to a window ENDPOINT and
delivered 2.57 against nitrogen's 2.08, a contrast of 1.238, outside the 0.879 to
1.164 the point estimates above support. Keeping that anchor would have needed
1.16 * 3.78 / 4.68 = 0.937 here and would have needed redoing every time either
leaf window moved. The repair was taken on both elements together as world-xms4;
the complete table across both elements and all four tissue proportions, with
what it cost the established C-N configuration, is
`plant-physiology-carbon-allocation-audit.md` finding 11.

## Sapwood: the constant is nitrogen's, and it is not the quantity the model applies

Three claims are stacked inside `PFRAC_LEAFTOSAP`, and having the papers
separates them. The element was never phosphorus and the tissue is not the
model's sapwood; the proportional form the constant asserts looked refuted on one
dataset and is not refuted on two.

### What 6.9 is in Friend et al. (1997)

Table 4, p. 254, gives `X_C:N(f/p)` = 0.145, "relative C:N ratio between foliage
and bark plus sapwood", in kgC kgN^-1 per kgC kgN^-1. Its reciprocal is 6.897.
The companion `X_C:N(f/r)` = 0.86 has reciprocal 1.163, which is
`PFRAC_LEAFTOROOT` and the nitrogen side's `frac_leaftoroot`. Two constants
matching two reciprocals to three figures identifies the source beyond doubt.

The text on p. 274 says what was fitted. `X_C:N(f/p)` is a mean across NINE
species from Turner (1980) and Turner and Lambert (1981), calculated with a
13 percent weighting for bark; those data are Eucalyptus spp. and planted
Douglas-fir. `X_C:N(f/r)` is from Nambiar and Fife (1991) for Pinus radiata
SEEDLINGS. Hybrid v3.0 applies both means to every one of its plant types, which
is what Table 4 is for, and uses them in Eqs. 42 to 44 to hold the relative C:N
ratios of foliage, sapwood plus bark, and fine roots FIXED while nitrogen is
reallocated each year.

Friend et al. (1997) contains no phosphorus. The word occurs once in the paper,
in the title of the Turner (1980) reference. So 6.9 is a nitrogen figure carried
across, it is a nitrogen figure for a tissue that includes bark, and it was
fitted on nine temperate plantation and eucalypt species.

### What the model applies is the constant itself, and that question is settled

What reaches the simulated plant is `canexch.cpp:1578`, where sapwood phosphorus
demand is computed against a target C:P of
`ctop_leaf_opt * ctop_sap_avr / ctop_leaf_avr`, and `ctop_leaf_opt` is that
individual's own leaf C:P. So the multiplier `ctop_sap_avr / ctop_leaf_avr` is a
proportion between two tissues' C:P, which is exactly the quantity both papers
measure.

`init_ctop_limits` sets `ctop_sap_avr = ctop_leaf_avr * PFRAC_LEAFTOSAP` and then
builds the window outward from that mean, `ctop_sap_max` at
`(1 + PFRAC_MAXTOMIN) / (2 * PFRAC_MAXTOMIN)` of it and `ctop_sap_min` at
`PFRAC_MAXTOMIN` of the maximum. `avg_ctop` is the harmonic mean, and that
inversion is the one that makes it return the mean again, so
`ctop_sap_avr / ctop_leaf_avr` = 6.9 exactly, for any `PFRAC_MINTOMAX` and any
`PFRAC_MAXTOMIN`. The check that has to pass, and it passes in closed form on the
source: after `init_ctop_limits()`, `ctop_sap_avr / ctop_leaf_avr` equals the
constant the source measured. It is not execution-verified, because no
LPJ-GUESS run this project has made has passed acceptance.

It did not always pass. The earlier form anchored the maximum instead, setting
`ctop_sap_max` to `ctop_leaf_max * PFRAC_LEAFTOSAP`, so the constant was a ratio
of window ENDPOINTS while the model applied a ratio of window MEANS. Because the
two windows had very different
widths, `PFRAC_MINTOMAX` = 3.68 for leaves against 1/0.9 for sapwood, the applied
proportion was

    [2 * 0.9 / 1.9] * PFRAC_LEAFTOSAP * (1 + PFRAC_MINTOMAX) / 2
        = 0.94737 * 6.9 * 2.34
        = 15.30

The nitrogen side was the same construction on its own window, 2.78, applying
12.36 against Friend's own 6.897, so it was a whole-model form question and not a
phosphorus one and could not be settled here alone. It was settled on both
elements together as world-xms4; the four-constant table and what the repair cost
the established C-N configuration are in
`plant-physiology-carbon-allocation-audit.md` finding 11.

### Heineman et al. (2016) tests the form on the right element and nearly the right tissue

Type II major axis regression of log species mean wood on log species mean leaf
concentration, n = 58, Table 3. Phosphorus: slope 2.10, 95 percent CI 1.52 to
3.17, r2 = 0.36; on phylogenetically independent contrasts, 2.22, CI 1.85 to
2.67, r2 = 0.54. Nitrogen: slope 1.25, CI 0.83 to 1.95, r2 = 0.31; contrasts
1.01, CI 0.84 to 1.23, r2 = 0.49. Calcium 1.93, potassium 1.97, magnesium 2.93.

"Significantly above 1" is the paper's own claim and not an inference from the
intervals: "For Ca, K, Mg and P, the slope of the wood vs leaf relationship (b)
was significantly > 1", against "the slope ... did not differ from 1 for N ...,
indicating that N scales isometrically between wood and leaf tissues". Nitrogen
is the only one of the five elements that is isometric, and it is the internal
control that makes the phosphorus result mean something: same 58 species, same
regression, same plots, and a method that rejected isometry for everything would
only be reporting its own sensitivity.

One hazard for anyone re-deriving off the paper. Table 3's caption prints the
model as `log(leaf) ~ log(a) + log(wood) x b`, which is the reverse of what was
fitted. The discussion says the fitted direction outright, "the wood-leaf scaling
exponent was ~2, meaning that, for example, a 10% increase in foliar P
corresponds to a 20% increase in wood P", and the intercepts settle it: at the
Table 2 species mean wood concentrations, wood-as-response returns leaf P
1117 ug/g and leaf N 1.90 percent, both ordinary tropical foliage, while
leaf-as-response returns leaf P 0.87 ug/g. Major axis regression is exactly
reciprocal-symmetric, so reading the caption literally does not relabel the
exponent, it inverts it to 0.48 and reverses the finding.

The scope is one regional gradient. The 58 paired species are all at Fortuna in
western Panama: six lower montane plots, 700 to 1500 m, mean annual temperature
19 to 23 C, annual rainfall 4000 to 9000 mm, on rhyolitic tuff, andesite and
porphyritic dacite with soil pH 3.6 to 5.6. The wider wood-only set is 106
species over 10 plots including four lowland Canal watershed sites, and no leaf
data. Wood is the outermost 5 cm annulus at breast height; leaves are three fully
expanded SHADE leaves from three individuals per species, ground with petioles
and rachii, collected in July 2010 against February 2011 for the Fortuna cores.
Shade leaves were used because 43 percent of the species never reach the canopy,
and the justification offered is that rank order of species mean N and P is
preserved between sun and shade leaves, which supports an exponent and not a
level. Wood P falls 35 percent from the outer 5 cm to the adjacent 5 to 10 cm
annulus, in 88 of 110 trees and in 14 of 18 species tested individually; wood N
has no consistent radial direction.

So it licenses the rejection of the fixed proportion for phosphorus AT THAT
SITE, with nitrogen as its own control. It does not license a global exponent:
these are tropical lower montane angiosperms, with no needleleaf, temperate or
boreal type in the sample, and the section below is what happens when a dataset
that does span those is read. It does not license a whole-sapwood tissue ratio,
because the outer annulus is the phosphorus-rich end of a 35 percent radial
gradient. It does
not license a sun-leaf level. And it does not license a scalar of any value,
because it is the scalar form that it rejects.

### The numbers, and the nitrogen check that could have failed

At the Table 2 species mean wood concentrations, 111 ug/g P and 2557 ug/g N, the
fitted lines put leaf P at 1117 ug/g and leaf N at 19027 ug/g. The tissue
proportions those imply, as leaf concentration over wood concentration, are 7.44
for nitrogen and 10.06 for phosphorus. A concentration proportion equals the
sapwood-over-leaf C:N or C:P proportion whenever the two tissues' carbon
fractions are alike, and species variation in tissue carbon is small against the
fourfold nitrogen and thirtyfold phosphorus variation here, so it moves the level
a little and leaves the exponent alone.

The nitrogen figure is the check. Friend's nine temperate species with bark
included give 6.90 and Heineman's 58 tropical species on the outer sapwood
annulus give 7.44: an 8 percent difference across two continents, two tissue
definitions and two methods. The construction therefore reproduces the nitrogen
constant, so the phosphorus number from the same construction is not an artefact
of the construction.

The phosphorus figure is not a constant. Because the fitted exponent is 2.10 the
proportion falls as leaf phosphorus rises, and over the observed span of species
mean wood P, 19 to 668 ug/g, it runs monotonically from 25.4 down to 3.9, a
factor of 6.4. The same computation on nitrogen over its observed 1300 to
5800 ug/g runs 8.5 to 6.3, a factor of 1.35, and the slope it comes from is not
separable from 1 anyway. A single number for phosphorus is not a value that is
merely uncertain; it is a value that does not exist.

The radial gradient moves the level and not the shape. `cmass_sap` is the whole
simulated sapwood and Heineman measured its outer 5 cm, so applying the 35
percent outer-to-inner decline to the whole pool raises the proportion from 10.06
to 15.5. That is the bracket on a forced scalar: 10.1 to 15.5, against the
nitrogen 6.9 the constant carries.

### What any of this is worth to the modelled result

The level first. The proportion the model applies is the constant itself, 6.9,
which sits below that 10.1 to 15.5 bracket, so simulated sapwood phosphorus
demand per unit sapwood carbon is 1.46 to 2.25 times what the tropical
measurement supports. The max-anchored form applied 15.30 and landed at the top
of the bracket by cancellation, a nitrogen constant too low for phosphorus
multiplied by a window factor of 2.22 that should not have been there; neither
half was a derivation and their product was not one either. Now that the applied
proportion is the declared one, moving the level is a change to the constant and
nothing else: 10.06 needs `PFRAC_LEAFTOSAP` = 10.06 and lowers sapwood phosphorus
demand per unit sapwood carbon by 31 percent, and 15.5 needs 15.5 and lowers it
by 55 percent.

The shape is worth more, and it is the part that is refuted. Across simulated
plant types leaf C:P varies only through `sla`, as `sla^-0.80936`. The woody
types' calculated `sla` runs from 9.300 for the needleleaf evergreens at three
years' leaf longevity to 26.03 for the broadleaf types at half a year, so their
leaf C:P spans a factor of 2.30, and under a fixed proportion their sapwood C:P
spans the same 2.30. Under the measured exponent it would span 2.30^2.10 = 5.75,
or 6.35 on the contrast-corrected slope.

Anchored on the tropical broadleaf evergreen type, which is what Heineman
measured, and at the observed 2.10: the boreal needleleaf evergreen type's
sapwood C:P would be 56 percent higher than the fixed proportion gives it, so its
sapwood phosphorus demand per unit sapwood carbon would be 36 percent lower; the
broadleaf summergreen types' sapwood C:P would be 37 percent lower, so their
demand would be 60 percent higher. Between those two the fixed form and the
measured form differ by a factor of 2.50 in RELATIVE sapwood phosphorus demand,
2.76 on the contrast-corrected slope. That is a competition effect among
simulated plant types on a large carbon pool, and it is the whole of what the
proportional form removes.

### Yan et al. (2016) runs the same test across 32 degrees of latitude, and the refutation does not survive

The same regression, on the same pair of tissues in the sense that matters --
woody stem against leaf, log-log, reduced major axis -- for 335 woody species in
198 genera and 73 families at 12 forest sites across eastern China, 18.7 to
50.9 degrees north, mean annual temperature -5.7 to 25.3 C and annual
precipitation 423 to 2031 mm, boreal coniferous forest through tropical
rainforest. The wood is the terminal 10 to 20 cm of the twig stems supporting the
sampled leaves; the leaves are fully expanded SUN leaves.

**The phosphorus exponent crosses 1 inside the climate range a global model
spans.** By biome it is 1.58 in tropical forest, 0.97 in temperate and 0.80 in
boreal. By site it runs 1.36 at Mt. Dinghu, 23.2 N, down to 0.71 at Mt. Genhe,
50.9 N. By functional group it is 1.26 for evergreen broad-leaved, 0.96 for
deciduous broad-leaved and 0.70 for coniferous plants, and 1.86 for legumes
against 0.88 for non-legumes. It correlates significantly with mean annual
temperature and not with precipitation or with soil total N or P.

The nitrogen exponent behaves the same way and less steeply: 1.30, 0.97 and 0.89
by biome, 1.45 to 0.74 by site, 1.20, 1.19 and 0.95 by group. So the phosphorus
exponent exceeds the nitrogen one in the tropics, which is Heineman's finding,
and falls below it towards the pole, which is the part one tropical site cannot
show.

**What that does to the form claim.** Heineman's 2.10 is the tropical extreme of
a gradient rather than a universal rejection of isometry, and a proportional form
is the unbiased default over the whole range a global model has to cover. Its
residual is one-signed within a biome and changes sign between them: at a fixed
scalar the model under-supplies the simulated tropical types' sapwood P demand
relative to leaf and over-supplies the boreal ones'. That is a statement about
which way the error goes, which the one-site reading could not make, and it is
strictly better than adopting a tropical exponent everywhere.

**What it does NOT do is supply the level, and it says why with a number.** The
leaf-to-wood phosphorus ratio is a property of which wood. Read off Yan's
Fig. 4a against its own axis, leaf P over twig stem P runs 0.77 at 18.7 N to
2.63 at 50.9 N -- a factor of 3.4 across the gradient FOR ONE TISSUE -- while
Heineman's outer 5 cm bole annulus brackets a forced scalar at 10.1 to 15.5. The
nitrogen counterparts sit the same way round: Yan's leaf N over twig stem N runs
1.8 to 3.0 where Friend's bole bark-plus-sapwood figure is 6.9. This model's
sapwood pool is the whole living sapwood of the simulated individual, whose mass
is overwhelmingly bole and large branches rather than terminal twigs, so its
value belongs at the bole end of that factor of twenty -- which is where the
phosphorus bracket sits and where the applied 6.9 does not.

### What is still open

ONE thing, and it is the level. The form is settled: the proportional relation
the constant asserts is not refuted, and the reading that said it was rested on
one site. What the model applies is the number the constant names, which is
world-xms4 and is also settled. What is left is that the constant's own source
measures a different element in a different tissue on nine temperate species,
and no accessible measurement gives the right element in this model's tissue.

Nothing above derives a value and nothing above should: a proportion kept
because it is convenient is a knob. Adopting Heineman's 10.1 to 15.5 would import
one Panamanian lower montane site's climate into every simulated plant type, and
Yan's 3.4-fold latitudinal span of the same ratio is the measure of what that
costs.

What would settle it is one measurement: paired leaf and WHOLE-SAPWOOD phosphorus
concentrations over more than one region. It is not in the accessible literature
and it is not a decision anyone can make from what is here.

A second finding falls out of Yan and is NOT this row: the leaf-to-wood
phosphorus ratio, not just its exponent, varies systematically with mean annual
temperature and by plant functional type, along the same axis the model's plant
types already distinguish. A single scalar applied to every type is wrong in a
direction the model could resolve. `world-3e5n` owns it.

`ifplim 1` keeps refusing and keeps naming `PFRAC_LEAFTOSAP` until the level is
answered. BIO-34.

## The saturation pair: one is its source's value, and the other's ramp is gone

`setptoc` ramps a soil organic matter pool's C:P from its maximum down to its
minimum, linearly, as a driving quantity rises from `fmin` to `fmax`. `PMASS_SAT`
and `PCONC_SAT` were the two `fmax` values, and both ramps were inert in
opposite directions. The reason was a different one for each, and so is the
repair: `PMASS_SAT` is its source's own number in the wrong currency and is
converted, and `PCONC_SAT` is the threshold of a ramp that should not exist and
is removed with it.

### `PMASS_SAT` is Parton, Stewart and Cole (1988) Fig. 3, read line for line

Figure 3, on page 115, plots the C:P ratio of the passive, slow and active soil
organic matter against soil labile P on an axis running 0 to 2.0 gP/m2. All
three are straight lines, all three start at their maximum C:P at zero labile P,
and all three reach their minimum at the right-hand end of the axis. 2.0 gP/m2
is 0.002 kgP/m2, which is `PMASS_SAT` exactly.

The three `(ctop_max, ctop_min)` pairs the fork passes to `setptoc` are that
figure's three lines: slow 200 to 90, passive 200 to 20, soil microbial 80 to 30,
against the paper's own text on the same page giving the slow fraction 90 to 200,
the passive 20 to 200 and the active 30 to 80. `setptoc`'s linear interpolation
is the figure's straight lines, `fmin = 0` is its origin, and the ratio it sets
is the P:C of the pool RECEIVING carbon, which is the paper's own construction:
"The organic P flows are calculated by multiplying the carbon flow rate times the
C:P ratio for the state variable receiving the carbon" (p. 115). The fork's
soil microbial pool is CENTURY's active pool, so the mapping is complete.

So the constant is derived, the functional form is the source's, and the
resemblance to the 0.002 that `NMASS_SAT` is built from is a coincidence of two
unrelated readings. The audit that opened this row read the coincidence as a
copy, and that was wrong.

One inconsistency inside the paper, recorded because it bounds the active pool's
lower end: the conclusions (p. 128) give the active range as 20 to 80 where the
model description and Fig. 3 give 30 to 80. The code follows the figure.

The conversion needs no time base and no planetary correction. It is a mass per
unit area against a mass per unit area, and the relation is a statement about
microbial and humic stoichiometry rather than about climate, orbit or day
length: phosphatase mineralisation is what makes soil organic matter C:P rise
when labile P is scarce, and microbial C:P is what makes it fall when labile P
is ample (McGill and Cole 1981, as the paper cites it). That transfers to this
world's simulated soils exactly as far as the rest of CENTURY does, which is a
declared Earth-analogue assumption and not a new one.

### What does not transfer is the pool the threshold reads

Parton's labile P is defined on the same page: orthophosphate that is
isotopically exchangeable or extractable with anion exchange resin, in the
0 to 20 cm soil the whole paper works in. This fork's `soil.pmass_labile` is a
different operational pool. It is the Hedley-labile pool of Wang et al. (2010)
and Yang et al. (2013), labile inorganic P plus labile organic P, and Dantas de
Paula et al. (2025) adopt that wider definition deliberately: they justify the
absence of biomineralisation processes in the model by the larger
plant-available pool the definition gives them.

Their own global figures separate the two definitions cleanly. Simulated labile
P is 2.11 PgP (Table 1), Yang's Hedley-labile estimate is 3.6 PgP over 0 to
50 cm, and Olsen-extractable P, the narrow operational pool, is 0.319 PgP over
0 to 20 cm. Over 1.3e14 m2 of ice-free land those are 16, 28 and 2.4 gP/m2.
Parton's 2.0 gP/m2 saturation point sits on the last of them, which is what it
should do, and the pool the fork feeds `setptoc` is eight times larger.

The fork's own sorption parameters say the same thing from the inside. `kplab`,
the Langmuir half-saturation for labile P, is 10 to 78 gP/m2 by soil order
(Wang et al. 2010 as tabulated in Dantas de Paula et al. 2025 Table A1), so
`PMASS_SAT` is 5 to 39 times below the labile P at which the fork's own isotherm
expects sorption to be halfway to saturating. Two parameter sets, from two
papers, on two operational definitions of "labile P", wired to one state
variable.

Yang and Post (2011) is the same finding from the other side: Hedley-labile P
exceeds vegetation demand even in strongly weathered soils, so it is not
plant-available P and must not be read straight into a model's labile pool.

Unconverted, `fac` therefore exceeds `fmax` nearly everywhere and the slow,
passive and soil microbial pools sit at their MINIMUM C:P always, which is their
most phosphorus-rich end. That is not a wrong constant. It is a right constant
reading a pool its source did not define.

What `setptoc` writes has exactly one live reader: `transferdecomp`'s `pinc`,
the phosphorus that rides carbon into a receiving pool. So this sets the
phosphorus content of every organic transfer in the model, and it is what the
reported soil and litter phosphorus stock is made of.

### The threshold is converted, and the conversion goes on the threshold

The two candidate changes read as a pair of equivalent framings and they are
not, because `PMASS_SAT` does a second job. `somfluxes` ends by pinning
`soil.pmass_labile` to `PMASS_SAT` whenever `!ifplim`, and the only things that
move it between that pin and the next call's `setptoc` are one day of uptake,
deposition, weathering and leaching. So under the configuration this project
runs, `fac` arrives at `fmax` to within a day's net phosphorus flux, which is a
fraction of a per cent of the axis.

- Converting the THRESHOLD into this fork's labile-P currency moves the pin with
  it, because the pin reads the same symbol and reads it for the same reason.
  `fac` tracks `fmax`, and nothing changes under `ifplim 0`. It bites only under
  `ifplim 1`.
- Driving `setptoc` with a resin-equivalent FRACTION of `soil.pmass_labile`
  leaves the pin alone, puts `fac` well below `fmax`, and moves the three pools
  partway up the ramp in the current configuration, where the reported stock is.

THE THRESHOLD IS THE ONE THAT IS CONVERTED. The pin is not a coincidental second
reader of the same symbol; it is DEFINED as holding labile P at the value where
the C:P ramp stops responding, exactly as the nitrogen side pins `NH4_mass` to
`NMASS_SAT`. Moving the threshold and letting the pin follow is therefore one
definition propagating rather than a side effect, and splitting the two into
separate constants would invent a second number with no independent derivation
and let the pair drift apart.

THE SCALAR IS 6.6. The ratio between this fork's `pmass_labile` and the
resin-extractable orthophosphate Parton's figure was drawn against is BRACKETED
at 6.6 to 11.3 by the fork's own published numbers, and the two ends are not two
estimates of one quantity:

- 6.6 is the model's SIMULATED labile P, 2.11 PgP, over Olsen-extractable
  0.319 PgP over 0 to 20 cm.
- 11.3 is the OBSERVATIONAL Hedley-labile estimate, 3.6 PgP, over the same Olsen
  figure.

`setptoc` reads the simulated pool, not the observation, so 6.6 is the ratio
between the two quantities that are actually wired together. The model
under-predicts its own observational target by 41 per cent, and 11.3 carries
that error as well as the definitional conversion, so using it would count the
under-prediction twice. `PMASS_SAT` is therefore `0.002 * 6.6` = 0.0132 kgP/m2,
which is 13.2 gP/m2 against a `kplab` of 10 to 78 and against the fork's own
simulated 16 gP/m2.

The bracket cannot be tightened from those numbers, and this is a second,
independent obstacle that the choice of end does not remove. Parton's is a 0 to
20 cm quantity; LPJ-GUESS-CNP simulates soil organic matter as a bulk pool with
no explicit depth, which Dantas de Paula et al. (2025) state as a known
limitation. The two sides of the ratio do not share a support, so there is no
depth on the receiving side to convert to, and narrowing the bracket needs a run
of this fork rather than more arithmetic.

### The conversion is a declared divergence, and what it moves

Vendored LPJ-GUESS-CNP has `PMASS_SAT = 0.002`, Parton's axis maximum carried
across unconverted, and this fork now diverges from it. Both the mainline value
and the applied one are recorded at the constant.

Under `ifplim 0`, the configuration this project runs, the divergence is INERT
for the carbon and phosphorus flows: the pin moves with the threshold, `fac`
still arrives at `fmax`, and every soil organic C:P ratio is where it was. What
does change is the reported stock. `commonoutput.cpp` writes the pinned
`pmass_labile` into `PO4_mass` and `availp`, so the reported labile P moves from
2.0 to 13.2 gP/m2, which is the same order as the fork's own simulated
16 gP/m2 instead of an order below it. It remains a constant meaning "not
limiting" and not a simulated stock; what the conversion removes is the reader's
false impression that the model is reporting a phosphorus-poor world.

Under `ifplim 1` the divergence bites for real. The emergent labile P now has a
threshold it can sit below, so the slow, passive and soil microbial pools ramp
instead of saturating, and the section below is the size of that.

### What the saturated ramp was worth, in the model's own reported stocks

The fork's published global run reports litter plus soil C of 1474.1 PgC and
litter plus soil P of 51.9 PgP, a bulk organic C:P of 28.4 by mass. That is
below the minimum C:P of every pool in the ramp except the passive one, so the
passive pool's carbon dominates the sum, as CENTURY's structure says it should.

Moving each pool from its minimum to the midpoint of its own range multiplies
its C:P by 5.5 for the passive pool, 1.61 for the slow and 1.83 for the soil
microbial, so the model's largest phosphorus stock would fall by a factor
between 1.6 and 5.5 depending on how the carbon is distributed among them, and
the reported bulk ratio says the high end governs. 51.9 PgP over that range is
9.4 to 32 PgP, against Yang et al. (2013)'s organic soil P of 8.6 plus or minus
6 PgP for 0 to 50 cm. The paper's own comparison column puts 51.9 PgP against
40.6 to 89 PgP, which are TOTAL soil P estimates including inorganic, secondary
and occluded fractions, not organic ones.

This is the size of the effect and not an argument for a value. A correct
threshold that worsened that comparison would still be the correct threshold.
What the arithmetic establishes is that the saturated ramp is worth up to a
factor of five on the largest phosphorus stock this model reports, so it is not
a tidiness question. That is the scale of what the converted threshold buys back
under `ifplim 1`, and it is not a claim about where the pools land, which needs
a run.

### The phosphorus-limitation-off pin is not a second defect

`PMASS_SAT` is also the value `pmass_labile` is pinned to whenever phosphorus
limitation is off, which is the configuration this project runs. That is the
right constant for the second job as well as the first: "saturated" here means
exactly "at the value where `setptoc` stops responding", the nitrogen side pins
`NH4_mass` to `NMASS_SAT` for the same reason, and the two uses move together
when the threshold moves, as they did under the conversion above.

What follows from it is a reporting hazard rather than a modelling one. Under
`ifplim 0` the labile P written to `PO4_mass` and `availp` is a CONSTANT meaning
"not limiting", identical in every simulated cell, and comparing it against a
measured labile P or against this fork's own phosphorus-limited run is
meaningless in either direction. The conversion moves it onto the same order as
the fork's own simulated pool, which makes it less misleading to read and no
more meaningful to compare. Read it as a flag.

Before the conversion, the soil organic C:P ratios were at the same place in
both configurations: under `ifplim 0` because `fac` was held at `fmax` by the
pin, under `ifplim 1` because the emergent labile P was far above it. That
coincidence is what made the two candidate changes behave differently, and it is
argued above.

### `PCONC_SAT` is removed, because the quantity its ramp modelled does not ramp

`PCONC_SAT` was compared against `litter_pmass / (litter_cmass * 2)`, the
phosphorus fraction of litter dry mass, and set the surface microbial pool's C:P
between 80 and 30.

**The measurement that was missing has been found, and it refutes the FORM.** A
decomposer community's biomass C:P is homeostatic with respect to its resource's
phosphorus content. Mooshammer et al. (2014) Table 2, recalculated from Xu et
al. (2013) over n = 405: microbial biomass C:P regressed on soil C:P has a slope
of 0.015 with R = 0.078, R2 = 0.006 and P = 0.118, and the log-log form gives
R = 0.000 and P = 0.992; the fitted microbial C:P moves from 66.5 to 67.3 while
the soil C:P it is regressed on moves from 156 to 1611. The same paper reports
the litter case directly: in decomposing litter, resource C:N and C:P are
strongly negatively correlated with the gross N and P MINERALISATION FLUXES
while the microbial communities are homeostatic in those element ratios
(Mooshammer et al. 2012), and Achat et al. (2010) find relatively constant
microbial biomass C:P in forest soils with the C:P of the mineralisation flux
varying strongly. What varies with a resource's phosphorus content is the flux
out of the decomposer, not the stoichiometry of the decomposer.

**The control that could have failed and did not** is the nitrogen row of the
same table. Microbial C:N against soil C:N has P = 0.044 and its log-log form
P < 0.001, where the phosphorus row is indistinguishable from flat. The table
separates the two elements rather than being too noisy to show anything, and the
nitrogen ramp beside the deleted one keeps its own direct source.

**What `PCONC_SAT` was is now exactly nameable.** `NCONC_SAT`, the constant it
copied, is sourced to the character: Parton et al. (1993) p. 791 says the C:N of
newly formed surface microbial biomass "increases from 10 to 20 as the N content
decreases from 2.0% to 0.01%", which is the `setntoc` call's pair and its `fmax`
as a mass fraction of litter dry mass. So `PCONC_SAT` was a sound derivation
belonging to the other element, and no phosphorus reading could have rescued it:
Parton, Stewart and Cole (1988) has no surface microbial pool and no C:P ramp
driven by a litter concentration at all, its litter P being a fixed structural
C:P of 500 with the remainder to the metabolic pool (p. 115).

**What the removal is worth is nearly nothing, which is why it is safe as well
as right.** Senesced-litter C:P is 660 to 1596 by mass across forest biomes
(McGroddy et al. 2004, Table 1), a litter phosphorus fraction of 7.6e-4 down to
3.1e-4. Against an `fmax` of 0.02 that drove the ramp over 1.6 to 3.8 per cent
of its declared span, so it returned a C:P of 79.2 to 78.1 against a declared 80
to 30, and the pool was already at the 80 `soil.cpp` initialises it to. The
removal moves the surface microbial pool's C:P by at most 2.4 per cent and makes
what the model runs visible where it is set.

**What is not settled by this.** The 80 the pool now holds for the whole of a
run is Fig. 3's ACTIVE SOIL line's `ctop_max` end applied to a pool that paper
does not have, which is the standing its three neighbours in `soil.cpp`'s
initialiser share. That is not what `PCONC_SAT` was and the removal did not
create it; `world-634q` owns it. The obvious substitution is refused there:
measured decomposer biomass C:P is 66.5 by MOLE, which is 25.8 by mass and the
model's ratios are mass ratios, and CENTURY's microbial pools are conceptual SOM
pools rather than measured biomass -- the same offset sits in the nitrogen side,
where the pool's sourced C:N of 10 to 20 stands against a measured microbial
biomass C:N of about 7 by mass.

### The same-relative-position transfer is refuted for one and unnecessary for both

The transfer recorded earlier gave `PCONC_SAT` = 9.3e-4, bracketed 7.1e-4 to
1.7e-3, by placing the phosphorus threshold at the same fraction of observed
mean litter concentration that the nitrogen threshold sits at. For `PMASS_SAT`
it offered 0.03 to 0.05 kgP/m2, fifteen to twenty-five times the current value,
by placing the threshold at the upper end of the observed labile-P distribution.

Both are now refuted, for different reasons, and NEITHER IS ADOPTED. The
`PMASS_SAT` construction is refuted outright: the paper gives the value
directly, it is 0.002, and a threshold fifteen times higher would have been
wrong against its own source. The `PCONC_SAT` construction is refuted by the
same measurement that removed the constant: it would have placed a threshold on
a ramp whose quantity does not ramp, so the position it was calibrating did not
exist to be found.

### The fork's published methods still cite the wrong Parton, and the right one is implemented anyway

Dantas de Paula et al. (2025) section 2.2 attributes the C:P ratios of the slow,
passive and active pools and their variation with the labile P pool to "the
CENTURY P submodel (Parton et al., 2010)", the ForCent paper, which contains no
phosphorus at all. Its own reference list carries Parton, Stewart and Cole (1988)
separately, for the leached-organic-P treatment.

The citation is wrong and the implementation is right: every number in those
three `setptoc` calls is the 1988 paper's, and the fork reproduces its
functional form. What the wrong citation cost was three years of nobody being
able to check the fourth argument.

### Two defects the same reading exposed, and they were one

Neither is a constant and both are in the phosphorus path. They sit in the same
branch of `somfluxes`, and closing the second closes the first.

`SURFHUMUS` was in the nitrogen ramp and not the phosphorus one, and the
phosphorus immobilisation branch scaled `sompool[SURFHUMUS].ptoc` down by
`ptoc_reduction` alongside `SLOWSOM` and `SOILMICRO`. Those two are re-derived
by `setptoc` at the top of every call to `somfluxes`; `SURFHUMUS` was not,
because its `setptoc` line was commented out. So its P:C ratcheted downward
without bound over a run, from the 1/150 `soil.cpp` initialised it to, and the
surface humus pool asymptotically received carbon carrying no phosphorus. The
invariant to keep is that a pool whose P:C is flexed down in that branch has to
be one `setptoc` re-derives. It was kept first by taking `SURFHUMUS` out of the
list, and it is kept now from the other side: the pool has a phosphorus ramp,
so the set `setptoc` covers and the set that branch flexes are one set again,
which is what WORLD-16PB asked for. `PASSIVESOM` is the harmless other
direction: re-derived, never flexed. What gave `SURFHUMUS` its ramp is the
section below.

### The surface humus pool's C:P is the slow pool's, and it ramps

The fork initialised `sompool[SURFHUMUS].ptoc` to 1/150 and never moved it, so
150 was the surface humus pool's C:P for the whole of a run and the phosphorus
content of every transfer into it out of `SURFSTRUCT`, `SURFFWD`, `SURFCWD` and
`SURFMICRO`. `transferdecomp` reads `sompool[receiver].ptoc` with no `ifplim`
guard, so that was live in both configurations.

150 has no source for a humus pool. The block cited Fig. 2 of Parton, Stewart
and Cole (1988), which is the P submodel's flow diagram and carries no C:P
values; the citation is corrected, and 150 appears in that paper only as the
lower bound on the C:P of new plant material for wheat (p. 112) and as the
structural litter C:N. Its three neighbours in the block are Fig. 3's `ctop_max`
ends and are overwritten on the first call to `somfluxes`, so their provenance
barely matters; 150 was the one that stood.

The objection that kept the `setptoc` line commented out was that the paper has
no humus pool, so its 200 and 90 would import `SLOWSOM`'s pair without a source.
They are not imported. The PAIR is Fig. 3's slow line. The IDENTIFICATION of
surface humus with the slow pool is the model's own, and this fork states it
twice:

- `setntoc` gives `SURFHUMUS` exactly `SLOWSOM`'s `(30, 15)` off exactly
  `SLOWSOM`'s driver, the mineral nitrogen pool. So the nitrogen side already
  ramps this pool as the slow pool at the surface, from the soil's mineral
  nutrient, and the phosphorus line is that treatment carried to the other
  element.
- `soil.cpp`'s own alternative P initialisation, commented out beside the live
  one, sets `SURFHUMUS` and `SLOWSOM` to the same 1/90, which is that line's
  phosphorus-rich end.

The alternative was a measured surface humus C:P, and it is the wrong KIND of
quantity, not merely a missing one. `setptoc` sets the P:C at which a pool
RECEIVES carbon: a stoichiometric target that moves with labile P by the
phosphatase mechanism of McGill and Cole (1981) that the paper builds on. A
measured forest-floor C:P is an emergent bulk ratio, and a fixed number is the
wrong shape for this argument however well sourced. That is what the fork's
1/150 was.

So `setptoc(soil, pmin_mass, SURFHUMUS, 200.0, 90.0, 0.0, PMASS_SAT)` runs, and
the initialisation moves to 1/200, the slow line's `ctop_max` end, which puts it
back on the block's own pattern and makes it not load-bearing again.

Both changes are declared divergences from the vendored CNP fork, with its own
lines recorded verbatim beside them. What moves: under `ifplim 0` the pin holds
`fac` at `fmax`, so the surface humus pool's C:P goes from a fixed 150 to a
fixed 90, and the phosphorus riding carbon into it rises by a factor of 1.67.
Under `ifplim 1` it ramps between 200 and 90 with labile P like the other three.
Nothing here is execution-verified.

The branch it sits in is the nitrogen branch above with `n` substituted for `p`,
down to every comment inside it still saying nitrogen, and the copy dropped
`|| !ifnlim` from the first condition. That omission decided which setting the
ratchet lived in, and it is the opposite of what it looks like. Under
`ifplim 1` the `else if (!ifplim)` arm is false and control goes to
`reduce_decay_rates`, so the ratchet was never reachable under phosphorus
limitation. Under `ifplim 0` it was the live arm. What kept it rare there is the
second job `PMASS_SAT` does: `somfluxes` ends by pinning `soil.pmass_labile` to
`PMASS_SAT` whenever `!ifplim`, so `pmin_mass` is exactly 2 gP/m2 at the start
of every call and the arm fires only on a day whose net phosphorus
immobilisation exceeds that. The whole daily carbon throughput of the soil
organic matter pools, carried at the richest C:P the model has, is an order of
magnitude short of it. So the ratchet was latent, not active, and it was latent
for a reason unrelated to the code that contained it.

`|| !ifplim` is restored, which makes `ifplim 0` mean what `ifnlim 0` already
means and leaves `ifplim 1` untouched. Both arms are now unreachable, and both
are kept, because upstream's form is what the pair is legible against.

## The strongly sorbed pool was a drain, not a pool

The audit that opened BIO-35 found that `UOCC` was declared and never used and
`Soil::pmass_occluded` never assigned, and concluded that the phosphorus cycle
has no terminal sink. Half of that is right and the consequence is the opposite
of what it looks like, because `Soil::pmass_strongly_sorbed` was never assigned
either.

It was declared in `guess.h`, initialised to zero in `soil.cpp`, serialised into
the restart, reported by `commonoutput.cpp`, and read at exactly one place:

    delta_strongly_sorbed = USORB * soil.pmass_sorbed - USSORB * soil.pmass_strongly_sorbed;
    pmass_add(soil, -delta_strongly_sorbed);
    patch.fluxes.report_flux(Fluxes::P_SOIL, delta_strongly_sorbed);

With the pool pinned at zero the back term vanishes, so this is not Wang et al.
(2010) Eq. D10, a transfer between two pools that equilibrates because `USORB`
and `USSORB` are equal. It is a first-order drain of the sorbed pool that can
never shut off, and the drained phosphorus is booked as an ecosystem loss.

The magnitude: Wang et al. (2010) Table A1, which this fork adopts, gives
`Spmax` of 77 to 145 gP/m2 by soil order, so the drain runs at 0.0067 times that,
0.52 to 0.97 gP/m2 per Earth year. The fork's default `soiltype.pwtr` is
3e-6 kgP/m2 per Earth year, 0.003 gP/m2 per Earth year. The spurious loss was
therefore 170 to 320 times the weathering supply, and under `ifplim 1` the soil
phosphorus system could not have reached a steady state at all.

It did not show up as a conservation failure for two reasons. `Patch::pcont`
excluded both the strongly sorbed and the occluded pool, so the phosphorus had
genuinely left the accounted system and the books balanced. And
`MassBalance::check_patch_P` is declared in `guess.h` and called from nowhere in
this tree, so the carbon, nitrogen and phosphorus balance checks do not run.

The fix is to assign the pool, which makes both terms of Eq. D10 real. The
strongly sorbed pool then fills to the size of the sorbed pool and the net flux
goes to zero, which is what the cited equation describes. The transfer is
internal, so it is no longer reported as a soil P loss, and `Patch::pcont` counts
the pool instead. Under `ifplim 0`, which pins `pmass_sorbed` to `spmax` every
step, the strongly sorbed pool now fills and the reported soil P loss falls to
zero once it has.

## Occlusion is absent by decision

Wang et al. (2010) section 2.3 states that the flux from strongly sorbed to
occluded P is not represented, on the argument that including it would add
computation with little impact at decade to century scales. The fork's other
phosphorus source disagrees: Parton, Stewart and Cole (1988) Fig. 2 carries an
occluded P box, and p. 117 gives the rate constant that fills it, K3 = 1.0e-6
per month applied to secondary P and multiplied by the same combined
moisture-temperature factor the decay rates use. At that factor's maximum of
one, K3 is 1.2e-5 per year. The vendored fork carried `UOCC = 1.0e-5` per year,
used nowhere, and that number is within 20 percent of the paper's -- so it was
undocumented in this tree rather than unsourced, and the audit that called it
uncited was reading the tree and not the source.

The decision is that Vesper's simulated phosphorus cycle declares the absence
rather than closing with an occlusion loss, and `UOCC` stays deleted rather than
wired up. Three reasons, and none of them is that omitting it improves an
agreement:

- Having the rate does not make the sink safe to add. Occluded phosphorus is
  terminal, so the rate directly sets this world's long-run soil phosphorus
  stock, and Parton's K3 acts on his secondary P pool over 0 to 20 cm in a
  monthly model, not on this fork's strongly sorbed pool at a daily step. That
  is the same pool-definition question the `PMASS_SAT` conversion answers for
  the C:P ramp, on a flux nothing can undo.
- `equilsom` spins the soil organic matter pools for 40000 model years to solve
  their equilibrium analytically. It is a numerical device and not 40000 years
  of this world's history, and `pmass_sorbed` and `pmass_strongly_sorbed` are
  not saved and restored across it the way `pmass_labile` is. A terminal sink
  inside that loop would drain the soil phosphorus of every cell over a span the
  run does not represent.
- A soil's pedogenic age on this world is set by hydrography and pedology, not
  by the length of a spin-up, so an occlusion loss belongs with the initial
  phosphorus stocks that ANUT-1 through ANUT-10 own rather than as a flux
  inside the vegetation model.

The same page carries the rest of that model's inorganic phosphorus rates, for
the component that owns weathering rather than for this one: K1 = 0.05 per month
for the formation of secondary P from labile P, K2 = 0.0022 per month for its
solubilisation back, and K4 = 1.0e-4 per month for the weathering of primary P,
all on the same moisture-temperature factor. The fork replaces K4 with Hartmann
and Moosdorf's climate-driven chemical weathering model, so the form is not
transferable, but a first-order decay of a primary P pool is what the cited
CENTURY submodel does and it is the alternative ANUT-1 through ANUT-10 are
choosing against.

What the decision costs, stated because it is a real absence: the modelled
phosphorus cycle's only losses are leaching, fire and harvest. It cannot
reproduce the Walker and Syers depletion of phosphorus over long pedogenesis,
and any simulated soil this project asks to be old must have that depletion
supplied in its initial stocks rather than developed by the model.

`Patch::pcont` counts `pmass_occluded` even though it is structurally zero, so
that adding an occlusion flux later is a change to one file rather than a silent
break in a conservation sum.

## What is still open

- `PFRAC_LEAFTOSAP`'s LEVEL. The form is settled and the applied proportion
  equals the declared constant, so this is a question about one number, and it
  needs paired leaf and whole-sapwood phosphorus concentrations over more than
  one region, which the accessible literature does not have.
- The surface microbial pool's fixed C:P of 80, as WORLD-634Q. It is Fig. 3's
  active soil line's `ctop_max` end applied to a pool that paper does not have,
  which is what its three neighbours in the initialiser also are, and removing
  the ramp made it load-bearing for the whole of a run rather than for the
  fraction of a per cent the dead ramp left it.
- Whether to represent terminal occlusion after all, now that the cited CENTURY
  submodel is held and does carry it, as WORLD-2LCW.
- Every derived value above is BRACKETED. A run that uses them has to say which
  end of each bracket it is on, and `parameters.cpp` refuses `ifplim 1`
  meanwhile.
- `MassBalance::check_patch_C`, `check_patch_N` and `check_patch_P` are declared
  and never called. A conservation check that does not run is why the strongly
  sorbed drain could stand.
