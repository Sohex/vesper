# The phosphorus cycle's parameters, and which of them are phosphorus

This is a worldbuilding project. Vesper is a simulated super-Earth around a
mid-K dwarf, and its biosphere is LPJ-GUESS-CNP, a vegetation model written for
Earth and ported to this world's calendar. Everything below is about that
model's phosphorus constants and the simulated soil and plant processes they
drive.

The C-N-P fork's phosphorus side was built by copying its nitrogen side. Each
copy is defensible on its own and together they mean the modelled phosphorus
cycle has no parameterisation of its own: uptake competition ranks the way
nitrogen competition ranks, leaf N:P is pinned, and the two soil saturation
thresholds are nitrogen's. This document says, for every one of those constants,
what it does to the modelled result, what the measured phosphorus value is or
what its bracket is, and which of them are still nitrogen's.

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
| `PMASS_SAT` | `somdynam.cpp` | `somfluxes` through `setptoc`, and the P-limitation-off pin | 0.002 | kgP/m2 labile P | UNDERIVED. Nitrogen's base value. Blocked on Parton, Stewart and Cole (1988) | no |
| `PCONC_SAT` | `somdynam.cpp` | `somfluxes` and `equilsom` through `setptoc` | 0.02 | half the litter P:C mass ratio | UNDERIVED. Nitrogen's value exactly. Blocked on the same paper | no |
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

## The saturation pair disables its own mechanism, twice, in opposite directions

`setptoc` ramps a soil organic matter pool's C:P from its maximum down to its
minimum as a driving quantity rises from `fmin` to `fmax`. `PMASS_SAT` and
`PCONC_SAT` are the two `fmax` values. Neither is derived for phosphorus, and
each does something worse than sit at the wrong magnitude.

`PCONC_SAT` is compared against `litter_pmass / (litter_cmass * 2)`, so 0.02 is
reached only at a litter C:P of 25 by mass. McGroddy et al. (2004) Table 1 gives
senesced-litter C:P of 3144 molar overall, which is 1219 by mass, and 1702 to
4116 molar across forest biomes, which is 660 to 1596 by mass. The threshold is
therefore 26 to 64 times below anything the model can produce, `fac` never
reaches `fmax`, and the surface microbial pool's C:P is pinned at its maximum of
80 always. The nitrogen constant is not in that position: the same 0.02 is a
litter C:N of 25 against an observed 57 by mass, a factor of 2.3, so the nitrogen
ramp does span its range at the rich end.

`PMASS_SAT` is compared against `soil.pmass_labile`, so 0.002 kgP/m2 is 2 gP/m2.
Yang et al. (2013) put global labile soil P at 3.6 plus or minus 3 PgP in the top
half metre, about 28 gP/m2 over ice-free land, and this fork's own published run
simulates 2.11 PgP, about 16 gP/m2. The threshold is 8 to 14 times BELOW the pool
it gates, so `fac` exceeds `fmax` nearly everywhere and the slow, passive and
microbial pools sit at their minimum C:P always.

`PMASS_SAT` has a second job: it is the value `pmass_labile` is pinned to
whenever phosphorus limitation is off, which is the configuration this project
runs. So the reported labile P of every simulated cell is currently an order of
magnitude below what this fork's own global run produces, for the same reason.

What would settle both is Parton, Stewart and Cole (1988), *Dynamics of C, N, P
and S in grassland soils: a model*, Biogeochemistry 5, 109-131,
`10.1007/BF02180320`, which `setptoc`'s own documentation cites and which this
project does not hold. It is not open access and was not obtainable.

The fork's published methods cite the wrong Parton for these values. Dantas de
Paula et al. (2025) section 2.2 attributes the C:P ratios of the slow, passive
and active pools and their variation with the labile P pool to "the CENTURY P
submodel (Parton et al., 2010)", which is the ForCent paper. That paper contains
no phosphorus at all. Its own reference list carries Parton, Stewart and Cole
(1988) separately, and that is where a CENTURY P submodel is.

A same-relative-position transfer off the nitrogen pair would give, for
`PCONC_SAT`: the nitrogen threshold sits at 25/57 = 0.44 of observed mean litter
C:N, and 0.44 of observed litter C:P is 536 by mass, which is
`PCONC_SAT` = 9.3e-4, bracketed 7.1e-4 to 1.7e-3 over the biome range. For
`PMASS_SAT` there is no equivalent, because the position of the nitrogen
threshold relative to the modelled mineral N pool is not stated anywhere; a
threshold placed at the upper end of the observed labile-P distribution would be
0.03 to 0.05 kgP/m2, fifteen to twenty-five times the current value. NEITHER IS
ADOPTED. Both would calibrate phosphorus against a nitrogen threshold whose own
position is not derived, and an undocumented value is a test case rather than an
anchor.

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
computation with little impact at decade to century scales. The vendored fork
carried a rate constant for it, `UOCC = 1.0e-5` per year, with no citation
anywhere in the tree, used nowhere.

The decision is that Vesper's simulated phosphorus cycle declares the absence
rather than closing with an occlusion loss, and `UOCC` is deleted rather than
wired up. Three reasons, and none of them is that omitting it improves an
agreement:

- The rate constant is undocumented. Occluded phosphorus is terminal, so an
  unbracketed rate would directly set this world's long-run soil phosphorus
  stock, which is the last place an invented number belongs.
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
- `PMASS_SAT` and `PCONC_SAT`, blocked on Parton, Stewart and Cole (1988),
  `10.1007/BF02180320`, which could not be obtained.
- Every derived value above is BRACKETED. A run that uses them has to say which
  end of each bracket it is on, and `parameters.cpp` refuses `ifplim 1`
  meanwhile.
- `MassBalance::check_patch_C`, `check_patch_N` and `check_patch_P` are declared
  and never called. A conservation check that does not run is why the strongly
  sorbed drain could stand.
