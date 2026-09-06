# What sets the simulated soil's C:N ratios, and where the citations run out

**Read:** 2026-08-30, from `references/pdf/parton1993-century-soil-organic-matter.pdf`,
`references/pdf/parton1988-century-c-n-p-s-grassland-model.pdf`,
`references/pdf/parton2010-forcent.pdf` and
`references/pdf/smith_2014_implications-of-incorporating-n-cycling-and-n-limitations-on-primary-p.pdf`,
against `vendor/lpj-guess/modules/somdynam.cpp` and `modules/soil.cpp`. No model
run: every number below is either read off a paper, measured off a published
figure against its own axis calibration, or evaluated from a declared line in
the source.

This is worldbuilding. Vesper is an invented super-Earth around a mid-K dwarf,
and the subject is the vegetation model that simulates its biosphere: the C:N
ratio at which each of that model's simulated soil organic matter pools receives
carbon, and what fixes it. A soil pool's C:N is a property of the material and
of the organisms that make it, so a measurement of one transfers to a planet
that was never measured; what does not transfer is a number with no measurement
behind it.

## The finding, in one line each

| Pool | The model runs | Its own paper says | Verdict |
| --- | --- | --- | --- |
| passive SOM | fixed 9, now a ramp 10 to 7 | Table C1: fixed 9, citing Parton et al. (2010) | **changed, SOURCED.** The cited paper has no C:N in it at all; Parton et al. (1993) Fig. 4(a) ramps this pool as one of its three lines, and the model called `setntoc` for the other two. WORLD-XFIF |
| soil microbial, which is CENTURY's active pool | ramp 15 to 6 | Table C1: 5 to 15, citing Parton et al. (2010) | **unsourced, and the code disagrees with its own paper by one unit at the floor.** Fig. 4(a) gives 15 to 2.14 and its text 3 to 15 |
| slow SOM | ramp 30 to 15 | Table C1: 15 to 30, citing Parton et al. (2010) | **unsourced.** Fig. 4(a) gives 20 to 12 in the figure and in the text alike, which is neither |
| surface humus | ramp 30 to 15 | Table C1: 15 to 30, citing Parton et al. (2010) | **unsourced.** CENTURY has no surface humus pool; the pair is the slow pool's |
| surface microbial | ramp 20 to 10 | Table C1: 10 to 20, citing Parton et al. (1993) | **sourced and correct.** Fig. 4(b) is exactly this line |
| `NMASS_SAT`, the fmax of all five | `0.002` kgN/m2 | the saturation point of Fig. 4, and Appendix C1's own cap | **changed, SOURCED.** The 0.002 is the figure's break point; the removed 0.05 had no derivation. WORLD-LNZN |
| the driver of the four soil ramps | `nmass_avail(NO)` | Fig. 4(a)'s abscissa is soil NO3 + NH4 | **changed, SOURCED.** `nmass_avail(NH4)` was the figure's quantity only under `ifntransform 0`, which this project does not run. WORLD-JUG1 |

## The citation that is empty, and it is the one the phosphorus side already found

Smith et al. (2014) is LPJ-GUESS's own documentation of this scheme. Its
Appendix C states that "the soil passive pool has a fixed C:N ratio of 9 (Parton
et al., 2010)", and that "for the soil microbial, surface humus and soil slow
pools, C:N ratio varies between upper and lower bounds depending on Navail
(Parton et al., 2010; Fig. C2)". Table C1 lists all four with the same
attribution.

`references/pdf/parton2010-forcent.pdf` is that paper, and this project holds it
and has read it. **The string "C:N" does not occur in it.** The word "nitrogen"
occurs twelve times, none of them giving a pool C:N ratio. The paper is about
the fraction of mineral soil organic matter in the slow and passive pools and
about their decay rates, fitted to radiocarbon; it carries no stoichiometry at
all.

That is the same defect SDEC-2 found on the phosphorus side, where Dantas de
Paula et al. (2025) attributes the CENTURY C:P ratios to the same 2010 paper,
which has no phosphorus in it either. One miscitation reaches both elements of
this model's soil organic matter, and in both cases the real source is an
earlier Parton.

## What the two Parton papers this project holds actually say

**Parton, Stewart and Cole (1988) fixes every soil pool's C:N and ramps only
C:P.** Its nitrogen submodel, p. 114: "We assume that the C:N ratio of
structural (150), active (8), slow (11), and passive (11) soil fractions remain
fixed." Its phosphorus submodel ramps the C:P of the active, slow and passive
pools against labile P, which is the Fig. 3 `setptoc` implements. So the paper
that settled this fork's C:P ramps deliberately does not ramp C:N at all, and if
it were the source of the passive pool's fixed value that value would be 11.

**Parton et al. (1993) Fig. 4 is the source `setntoc` and `NMASS_SAT` cite, and
it ramps three soil pools, not two.** Panel (a) plots the C:N ratio of the
active, slow and passive SOM pools against soil NO3 + NH4 in gN/m2 -- the whole
mineral pool, which is why the ramps read `nmass_avail(NO)`; panel (b)
plots newly formed surface microbial biomass against surface litter N content.
The passive pool is one of panel (a)'s three lines, and it is the one this model
never called `setntoc` for. The fork's own comment at the initialiser asked why:
"passive has a fixed value (why? passive SOM should also vary.)"

## Figure 4(a), measured

Rendered at 600 dpi and calibrated against its own axes. The vertical ticks sit
at 20, 15, 10, 5 and 0 at 76.8 px per unit; the horizontal ticks at 0.0, 0.5,
1.0, 1.5 and 2.0 gN/m2 at 888 px per gN/m2. The lines are sampled at 0.1 gN/m2
rather than at the origin, where the axis tick marks overlap them, and the
intercept is the linear segment extrapolated back:

| Line | at 0.1 gN/m2 | intercept | flat, from 2.0 gN/m2 to the axis end | p. 791's text |
| --- | --- | --- | --- | --- |
| slow SOM | 19.68 | 20.08 | 12.06 | 12-20 |
| active SOM | 14.41 | 15.06 | 2.14 | 3-15 |
| passive SOM | 9.68 | 10.02 | 3.17 | 7-10 |

All three lines break at 2.0 gN/m2 and are flat to the axis end at 2.5.

**The figure and the text of the same paper disagree, and it is not a reading
error.** The slow line agrees at both ends. Every upper end agrees exactly. What
is in dispute is only where two of the three lines are drawn to stop falling:
the active line's floor is 2.14 against a stated 3, and the passive line's is
3.17 against a stated 7.

**The text's arm is the one taken, on a check that could have gone the other
way.** This model resolves the same conflict the same way at the one pool where
it implements the ramp and both statements exist: Smith et al. (2014) Table C1
gives the soil microbial pool 5 to 15 where Fig. 4(a)'s active line reaches
2.14, so LPJ-GUESS reads the text's arm, and the code takes 6. Following the
figure for the passive pool instead would give it a C:N of 3.2, which is below
the C:N of soil microbial biomass itself and therefore below anything the
material the pool stands for is measured at; the text's 7 is not. The figure's
low ends are the part of it that is schematic.

## What the change is worth, and why it is not what the endpoints span

`sompool[PASSIVESOM].ntoc` is the N:C at which the passive pool RECEIVES carbon,
applied in `transferdecomp()` to the flows from the slow pool and from the soil
microbial pool. It therefore sets the nitrogen immobilised into the pool with
the longest residence time in the model: `K_MAX` 1.9e-6/day, of order 1400
years. What is locked there is out of circulation on that timescale, so it sets
the steady-state mineral nitrogen after `equilsom`.

The ramp spans, and what the change is worth is therefore a function of the
driver rather than a factor. Against the fixed 9 the passive pool locks
10/9 = 1.11 times LESS nitrogen per unit carbon at an empty mineral pool and
9/7 = 1.29 times more at saturation, crossing 9 at 0.667 gN/m2. Which side of
that a cell sits on is one distribution off the first run that reaches output.

It spans because `NMASS_SAT` is the figure's break point. The constant was
written `0.002 * 0.05`. The 0.002 is Fig. 4's break point in the model's own
units, and Smith et al. (2014) Appendix C1 states independently that the model's
mineral N pool "is capped at a saturation level of 2 g N m-2 following Parton et
al. (1993)" -- the same number for the same reason, and stated for the fixation
ceiling as well as for the ramps. The 0.05 had no derivation anywhere in the
tree; what stood beside it was one attributed opinion, "NMASS_SAT is too high
when considering BNF - Zaehle", with no argument and no number behind it. On
that fmax all five `setntoc` ramps were declared, evaluated, and effectively
inert.

**The objection the scalar encoded is answered by the document that supplies the
number.** Appendix C1 says BNF "is distributed equally throughout the year and
added directly to the soil-available mineral N pool, Navail, which is capped at
a saturation level of 2 g N m-2 following Parton et al. (1993). BNF in excess of
the saturation level is discarded (assumed not to have occurred)", and that
deposition above the same level goes to leaching. So the fixation ceiling the
comment worried about is documented AT the ramp threshold's value, and
`soilnadd()` is the code that implements it. A second mechanism reading the
symbol does not change what the symbol means.

## What this did not establish

- **Where 30/15 and 15/6 come from.** Table C1's 15-30 and 5-15 are attributed
  to a paper containing neither, and neither pair is Fig. 4(a)'s 20-12 and
  15-2.14. Three of this model's five nitrogen ramps have no reachable source
  for their endpoints, and the soil microbial floor in the code (6) is not even
  the one in the model's own paper (5).
- **The passive pool's base decay rate.** Smith et al. (2014) Table C1 gives
  3.9e-6/day; `somdynam.cpp`'s `K_MAX` array gives 1.9e-6. A factor of two on
  the turnover time of the slowest pool in the model, which is the residence
  time every number above is weighed by.
- **Where the simulated `nmin_mass` sits against 2 gN/m2.** That does not
  decide `NMASS_SAT`, which is stated by the figure and by the model's own
  documentation; it decides the SIGN and size of what these ramps are worth,
  pool by pool, and it is one distribution off any run that reaches output.
- **The nitrate share of the simulated mineral pool.** That is the size of the
  driver change, and it is one reported quantity off the same run.
- **The original CENTURY parameter values.** Parton, Schimel, Cole and Ojima
  (1987), `10.2136/sssaj1987.03615995005100050015x`, is what Parton et al.
  (1993) cites for the justification of its N submodel and is where the
  figure-against-text conflict would be settled from outside. It could not be
  fetched: no open-access route, and no reachable PDF at the publisher.

## Sources

- Parton, Stewart, Cole (1988). *Dynamics of C, N, P and S in grassland soils: a model.* Biogeochemistry 5(1), 109-131. `10.1007/BF02180320`. p. 114 for the fixed soil C:N ratios; Fig. 3 p. 115 for the C:P ramps
- Parton, Scurlock, Ojima, Gilmanov, Scholes, Schimel, Kirchner, Menaut, Seastedt, Garcia Moya, Kamnalrut, Kinyamario (1993). *Observations and modeling of biomass and soil organic matter dynamics for the grassland biome worldwide.* Global Biogeochem. Cycles 7(4), 785-809. `10.1029/93GB02042`. Fig. 4 and p. 791
- Parton, Hanson, Swanston, Torn, Trumbore, Riley, Kelly (2010). *ForCent model development and testing using the Enriched Background Isotope Study experiment.* J. Geophys. Res.-Biogeo. 115, G04001. `10.1029/2009JG001193`. Read, and it contains no C:N ratio
- Smith, Warlind, Arneth, Hickler, Leadley, Siltberg, Zaehle (2014). *Implications of incorporating N cycling and N limitations on primary production in an individual-based dynamic vegetation model.* Biogeosciences 11, 2027-2054. `10.5194/bg-11-2027-2014`. Appendix C and Table C1
