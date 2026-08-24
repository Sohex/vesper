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

Nothing here is verified by execution. LPJ-GUESS does not build on this tree:
`framework/vesper.h` is generated, and the chain to it runs through a baseline
run and a baseline climatology that do not exist. Every figure below is hand
arithmetic off the source and the cited measurements.

## The register

| constant | where it lives | read at | value | unit | source or bracket | moved |
| --- | --- | --- | --- | --- | --- | --- |
| `PUPS_UPPER_ADV` | `guess.h` | `Pft::init_pupscoeff`, then `canexch.cpp:phosphorus_uptake_strength` | 3.18 | dimensionless | Jobbagy and Jackson (2001) profile contrast on the model's own nitrogen value; BRACKET 2.9 to 3.5 | yes, from 2.0 |
| `PFRAC_MINTOMAX` | `guess.h` | `Pft::init_ctop_limits` and `Pft::init_ctop_min` | 3.68 | dimensionless | McGroddy et al. (2004) foliar dispersion contrast; BRACKET 3.68 to 4.41 | yes, from 2.78 |
| `PFRAC_MINTOMAX_CROPGREEN` | `guess.h` | `Pft::init_ctop_limits` | 7.77 | dimensionless | the same transfer applied to its nitrogen counterpart 5.0; cropland is refused under BIO-27 | yes, from 5.0 |
| `PFRAC_LEAFTOROOT` | `guess.h` | `Pft::init_ctop_limits` | 1.16 | dimensionless | Yuan et al. (2011); root N:P is not separable from leaf N:P, so the nitrogen proportion carries; BRACKET 1.02 to 1.35 | no, and now for a reason |
| `PFRAC_LEAFTOSAP` | `guess.h` | `Pft::init_ctop_limits` | 6.9 | dimensionless | UNDERIVED. Heineman et al. (2016) reject the fixed-proportion form for phosphorus | no |
| `PFRAC_MAXTOMIN` | `guess.h` | `Pft::init_ctop_limits` | 0.9 | dimensionless | nothing to transfer: the nitrogen original is declared arbitrary in `Pft::init_cton_limits` | no |
| `PMASS_SAT` | `somdynam.cpp` | `somfluxes` through `setptoc`, and the P-limitation-off pin | 0.002 | kgP/m2 labile P | Parton, Stewart and Cole (1988) Fig. 3 p. 115, whose labile-P axis saturates at 2.0 gP/m2. DERIVED, and the driving pool is not the paper's | no, and the value is not the defect |
| `PCONC_SAT` | `somdynam.cpp` | `somfluxes` and `equilsom` through `setptoc` | 0.02 | phosphorus fraction of litter dry mass | UNDERIVED. Nitrogen's value exactly, on a ramp the cited source does not contain. Bounded above by 7.6e-4 | no |
| `USORB` | `somdynam.cpp` | `somfluxes` | 0.0067 / `VESPER_EARTH_YEAR_DAYS` | per absolute day | Wang et al. (2010) Appendix D | divisor, under the time-base contract |
| `USSORB` | `somdynam.cpp` | `somfluxes` | 0.0067 / `VESPER_EARTH_YEAR_DAYS` | per absolute day | the same, and equal to `USORB` in that source | divisor, under the time-base contract |
| `UOCC` | was `somdynam.cpp` | nowhere | removed | -- | no citation anywhere in the tree | deleted |
| `Soiltype::pwtr` | the soil input | `somdynam.cpp:soilpadd` | per BIO-5 | kgP/m2 per EARTH year | BIO-5 emits it against `biosphere/notes/time-base-unit-contract.md` | divisor, under that contract |

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
leaf P content falls 37.1 percent; and `ctop_root_avr` and `ctop_sap_avr`, which
are proportional to `ctop_leaf_max`, move with the latter.

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

## Sapwood: the form is refuted, so there is no scalar to derive

`PFRAC_LEAFTOSAP` is 6.9 from Friend et al. (1997), and a fixed proportion is a
defensible form for nitrogen. Heineman et al. (2016) regress wood on leaf
nutrient concentrations across 58 species and cannot distinguish the nitrogen
exponent from 1, which is what a fixed proportion asserts. For phosphorus the
same regression returns an exponent near 2, significantly above 1, and it is the
strongest of their wood-leaf relationships.

So wood C:P is not a fixed multiple of leaf C:P, and no scalar derived from that
measurement would mean what this constant claims to mean. Deriving one anyway
from the concentration means, 2557 ug/g wood N and 111 ug/g wood P against
tropical foliar N:P, gives about 8. It is recorded and not adopted: it rests on
one tropical gradient and on a form the same measurement rejects.

Choosing between a refuted scalar and a nonlinear wood-leaf phosphorus relation
is a modelling decision. It is what BIO-34 is still open on, and sapwood is a
large carbon pool, so it is not a small one.

## The saturation pair: one is its source's value, and neither ramp works

`setptoc` ramps a soil organic matter pool's C:P from its maximum down to its
minimum, linearly, as a driving quantity rises from `fmin` to `fmax`. `PMASS_SAT`
and `PCONC_SAT` are the two `fmax` values. Both ramps are inert in the same
direction the audit found them, and the reasons are different for the two
constants.

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

`fac` therefore exceeds `fmax` nearly everywhere and the slow, passive and soil
microbial pools sit at their MINIMUM C:P always, which is their most
phosphorus-rich end. That is not a wrong constant. It is a right constant
reading a pool its source did not define, and the two ways out -- converting the
threshold into this fork's labile-P currency, or driving `setptoc` with a
resin-equivalent fraction of that pool -- are modelling decisions rather than
arithmetic. Neither is taken here.

A second, independent reason the threshold cannot simply be rescaled: Parton's
is a 0 to 20 cm quantity, and LPJ-GUESS-CNP simulates soil organic matter as a
bulk pool with no explicit depth, which Dantas de Paula et al. (2025) state as a
known limitation. There is no depth on the receiving side to convert to.

### What the pinned ramp is worth, in the model's own reported stocks

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
What the arithmetic establishes is that the disabled ramp is worth up to a
factor of five on the largest phosphorus stock this model reports, so it is not
a tidiness question, and that no comparison against a measured soil organic P is
worth making until the driving pool is settled.

### The phosphorus-limitation-off pin is not a second defect

`PMASS_SAT` is also the value `pmass_labile` is pinned to whenever phosphorus
limitation is off, which is the configuration this project runs. That is the
right constant for the second job as well as the first: "saturated" here means
exactly "at the value where `setptoc` stops responding", the nitrogen side pins
`NH4_mass` to `NMASS_SAT` for the same reason, and the two uses move together if
the threshold ever moves.

What follows from it is a reporting hazard rather than a modelling one. Under
`ifplim 0` the labile P written to `PO4_mass` and `availp` is a CONSTANT meaning
"not limiting", identical in every simulated cell, and comparing it against a
measured labile P or against this fork's own phosphorus-limited run is
meaningless in either direction. It is eight times below the latter and fourteen
below the former, and neither number says anything. Read it as a flag.

The soil organic C:P ratios are at the same place in both configurations: under
`ifplim 0` because `fac` is pinned exactly to `fmax`, under `ifplim 1` because
the emergent labile P is far above it.

### `PCONC_SAT` has no source in that paper, or anywhere in the tree

`PCONC_SAT` is compared against `litter_pmass / (litter_cmass * 2)`, the
phosphorus fraction of litter dry mass, and sets the surface microbial pool's
C:P between 80 and 30.

Parton, Stewart and Cole (1988) contains no counterpart to that ramp. The model
has no surface microbial pool, and its only C:P ramps are the three soil pools
of Fig. 3, driven by labile P. Litter P in that model is not ramped at all: the
structural pool is fixed at C:P 500 and the metabolic pool receives the
remainder of the plant residue P (p. 115). So the ramp `PCONC_SAT` belongs to is
the nitrogen side's structure carried across, the pair 80 and 30 is Fig. 3's
ACTIVE SOIL line applied to a surface pool, and the value is `NCONC_SAT`
unchanged.

The bound stands, and the paper does not move it. Senesced-litter C:P is 660 to
1596 by mass across forest biomes (McGroddy et al. 2004, Table 1), which is a
litter phosphorus fraction of 7.6e-4 down to 3.1e-4. `PCONC_SAT` at 0.02 is 26
to 64 times above the richest litter the model can produce, so `fac` never
reaches `fmax` and the surface microbial pool sits at its MAXIMUM C:P of 80
always. Any replacement that lets the ramp span at all is at or below 7.6e-4.
Nothing in the cited source anchors it from below, so it stays a bound and not a
bracket, and `parameters.cpp` keeps refusing.

Two things would settle it, and choosing between them is the decision:

- A measurement of surface-litter microbial biomass C:P against litter
  phosphorus concentration, which is what the nitrogen constant has in Parton
  et al. (1993) Fig. 4 and what phosphorus does not.
- Deleting the constant, by driving the surface microbial pool's C:P from the
  same labile P the soil pools use, which is what Parton does for the active
  pool. The reason the nitrogen code does not is real -- surface litter is not
  in contact with the mineral soil's available nitrogen -- so this is a
  structural argument and not a simplification.

### The same-relative-position transfer is refuted for one and unnecessary for both

The transfer recorded earlier gave `PCONC_SAT` = 9.3e-4, bracketed 7.1e-4 to
1.7e-3, by placing the phosphorus threshold at the same fraction of observed
mean litter concentration that the nitrogen threshold sits at. For `PMASS_SAT`
it offered 0.03 to 0.05 kgP/m2, fifteen to twenty-five times the current value,
by placing the threshold at the upper end of the observed labile-P distribution.

That second construction is now refuted outright: the paper gives the value
directly, it is 0.002, and a threshold fifteen times higher would have been
wrong against its own source. The first is not refuted but is unnecessary in the
same way -- it would calibrate a phosphorus threshold against a nitrogen
threshold's position, and the phosphorus problem is not that the position is
unknown but that the ramp has no source at all. NEITHER IS ADOPTED.

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

### Two defects the same reading exposed

Neither is a constant, both are in the phosphorus path, and both bite only under
`ifplim 1`.

`SURFHUMUS` is in the nitrogen ramp and not the phosphorus one, and the
phosphorus immobilisation branch scales `sompool[SURFHUMUS].ptoc` down by
`ptoc_reduction` alongside `SLOWSOM` and `SOILMICRO`. Those two are re-derived
by `setptoc` at the top of every call to `som_dynamics_century`; `SURFHUMUS` is
not, because its `setptoc` line is commented out. So its P:C ratchets downward
without bound over a run, from the 1/150 it is initialised to in `soil.cpp`,
and the surface humus pool eventually receives carbon carrying no phosphorus.
Uncommenting the line is not the fix on its own: Parton, Stewart and Cole (1988)
has no humus pool to take a `(ctop_max, ctop_min)` pair from, so that would
import `SLOWSOM`'s pair without a source.

The nitrogen immobilisation branch short-circuits on `|| !ifnlim` and the
phosphorus branch has no `|| !ifplim`, so the phosphorus decay-rate reduction
and the `ptoc` ratchet can both fire in a run with phosphorus limitation off.
Reaching them needs daily immobilisation above the pinned 2 gP/m2, which is not
a realistic daily flux, so this is an asymmetry to close rather than an active
defect in the current configuration.

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
  is the same pool-definition question `PMASS_SAT` is caught in, on a flux
  nothing can undo.
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

- `PFRAC_LEAFTOSAP`, and the modelling decision behind it: a scalar on a form
  the measurement rejects, or a nonlinear wood-leaf phosphorus relation.
- `PMASS_SAT`'s driving pool, as WORLD-Z01O. The constant is settled and is not
  the defect; what is open is that `soil.pmass_labile` is a Hedley-labile pool
  and the threshold is a resin-extractable one, eight times smaller.
- `PCONC_SAT`, as WORLD-PIDX, which has no phosphorus source in the paper the
  ramp cites or anywhere else in the tree, and whose ramp that paper does not
  contain.
- The `SURFHUMUS` P:C ratchet and the missing `!ifplim` short-circuit beside it,
  as WORLD-16PB.
- Whether to represent terminal occlusion after all, now that the cited CENTURY
  submodel is held and does carry it, as WORLD-2LCW.
- Every derived value above is BRACKETED. A run that uses them has to say which
  end of each bracket it is on, and `parameters.cpp` refuses `ifplim 1`
  meanwhile.
- `MassBalance::check_patch_C`, `check_patch_N` and `check_patch_P` are declared
  and never called. A conservation check that does not run is why the strongly
  sorbed drain could stand.
