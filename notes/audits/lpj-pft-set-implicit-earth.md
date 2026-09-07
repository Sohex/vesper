# The LPJ-GUESS plant functional type set is implicit-Earth, and where that bites

Audited 2026-09-07 against `biosphere/generated/vesper_pfts.ins` as generated
for `lpj_1e6a2b9ca51a4eff9592992cad96677b`, the accepted run. This is a note
about a simulated world: every quantity is a modelled field of Vesper or a
parameter of the model that produces it.

**This ground was declared out of scope by both existing implicit-Earth
audits.** `notes/audits/inherited-earth-constants.md` says so in as many words
-- "I did not audit ... the LPJ-GUESS PFT parameter values themselves" -- and
`notes/audits/model-earth-centrism.md`, which took the same question into the
fork, does not mention `tcmin_surv`, phenology, rooting or establishment at
all. So the twelve types' parameter values have never been through the sweep
every other inherited constant in this tree has been through. This is that
sweep. The task is `world-orok`.

## What is NOT in this class, so the class is not overread

The audit found 78 distinct parameters across the group and PFT blocks. Three
groups of them carry no planetary content and are not findings:

- **Time base, and it is already correct.** `build_vesper_pfts.py` converts
  `gdd5min_est` and `greff_min` as ANNUAL_SUM, `longevity`, `leaflong`,
  `distinterval` and `freenyears` as YEAR_COUNT, and `turnover_leaf`,
  `turnover_root` and `turnover_sap` as ANNUAL_RATE by `1-(1-r)^f` rather than
  by multiplication. `phengdd5ramp` is deliberately untouched as a within-season
  accumulation. The registry is `biosphere/notes/time-base-unit-contract.md` and
  it is sound; nothing in this note is a time-base defect.
- **Model convention.** `include`, `landcover`, `intercrop`, `restrictpfts`,
  `lifeform`, `leafphysiognomy`, `phenology`, `pathway`, and the managed-land
  bookkeeping (`harv_eff`, `harvest_slow_frac`, `res_outtake`,
  `turnover_harv_prod`) select submodels or are unreachable in a natural run.
- **Inactive.** `eps_iso`, `eps_mon`, `storfrac_mon` and `seas_iso` are BVOC
  emission capacities and the run sets `ifbvoc 0`; `has_aerenchyma`,
  `inund_duration` and `wtp_max` are wetland and the run sets `run_peatland 0`.

Everything else -- the bioclimatic limits, the photosynthesis temperature
response, chilling and phenology, establishment, the water relations, allometry,
fire resistance, and the storage and uptake fractions -- is a value observed on
Earth's flora under Earth's forcing. That is the class.

## The finding: Earth's parameter set is ANTICORRELATED where this world is not

Cold tolerance and heat tolerance are not independent in Earth's flora, because
they are not independent in Earth's climate. Earth's cold places are LOW-ENERGY
places, so a plant that survives a hard winter has never been selected to
photosynthesise well in a hot summer. The shipped parameters carry that
correlation exactly:

| PFT | `pstemp_high` | `pstemp_max` | `tcmin_surv` |
| --- | --- | --- | --- |
| BNE, BINE | 25 | 38 | -31 |
| BNS | 25 | 38 | none |
| IBS | 25 | 38 | -30 |
| TeNE, TeBS, TeBE | 25 | 38 | -2, -14, -1 |
| TrBE, TrIBE, TrBR | 30 | 55 | +15.5 |
| **C3G** | **30** | **45** | **none** |
| C4G | 45 | 55 | +15.5 |

Every type with a 30 degC optimum needs a coldest month above +15.5 degC,
**except C3G**, which is the single row combining an unlimited cold survival
with a 30 degC optimum plateau.

**Vesper breaks the correlation.** Its polar cap runs a coldest month of
-68.75 degC and a warmest month of +31.11 degC, with 844 growing degree-days
above 5 degC. A plant there must survive a winter no Earth PFT but BNS, IBS and
C3G is parameterised for AND photosynthesise in a summer hotter than any of
those three is optimised for. Earth has no such place, so Earth's set contains
one row that spans it by accident.

That row is the one that wins. Poleward of 75 degrees, of 0.0300 total cover,
**C3G carries 0.0274 and BNS 0.0027**, and the other ten types are at exactly
zero.

## The cold threshold is Earth's STATISTIC as much as Earth's number

Larcher (2005), fetched for this row, gives hardened frost resistance BY TISSUE
rather than by climate index, and it reframes the filter above.

| hardened boreal and alpine conifers | threshold with injury |
| --- | --- |
| leaves | -40 to below -70 degC |
| shoot buds | -40 to below -70 degC |
| twigs and stems | -50 to below -70 degC |
| **roots** | **-20 to -30 degC** |

So a -69 degC winter is INSIDE Earth's own boreal envelope for the aerial
parts, not beyond it, and `tcmin_surv` of -31 is not a tissue threshold at all
-- it is a coldest-month MEAN, an index calibrated on Earth where a monthly mean
of -31 implies absolute minima far below it. Two Earth-derived things are
therefore imported at once: the number, and the statistic that relates it to
what tissue experiences.

On this world the statistic does not carry over, and it happens not to matter
here for the reason worth recording: at the cap the coldest-month mean is
-68.75 degC and the coldest-month MINIMUM is -69.59, a gap of **0.84 K**. The
polar night is so uniformly cold that mean and minimum coincide, where on Earth
they differ by tens of kelvin. The mean-based criterion is accidentally safe
here and would not be at a lower latitude.

**The root threshold is the one that bites, and the climate model does not
report the layer it bites in.** `dsoilz` gives the soil column five layers of
0.4, 0.8, 1.6, 3.2 and 6.4 m, so layer 1 spans 0 to 0.4 m and its midpoint is
Larcher's 20 cm. That layer is not written by the model at all, and this
project already records why: `exoplasim/scripts/close_state_energy.py` notes
that the model writes layers 2 to 5 "and not layer 1, whose 0.4 m is tied to the
surface temperature by the implicit top-layer solve. `ts` stands in for it".
`landmod.f90` bears that out -- the top-layer flux is solved as
`2*zdiff1/zsoilz1*(dts - dsoilt(:,1))`, so the two are coupled rather than
independent. **The surface temperature is therefore the stand-in for the root
zone**, and it is the project's own choice for the column's heat budget rather
than one made here:

| on the model's own land, poleward of 75 deg | coldest month |
| --- | --- |
| surface (`ts`), **the stand-in for layer 1** | **-79.53 degC** |
| layer 1, 0-0.4 m, the root zone | not written; tied to `ts` above |
| layer 2, 0.4-1.2 m (`tso2`) | -64.69 degC |
| layer 5, 6.0-12.4 m (`tsod`) | -31.96 degC, and its annual mean -31.85 matches the surface's -31.74, which is what says the column is equilibrated |

So the root zone runs at about **-79.5 degC** in the coldest month, colder than
either figure an earlier draft reached, against a boreal conifer root tolerance
of -20 to -30 degC. The conclusion holds and strengthens; the number does not,
and an
earlier draft of this note quoted -63.91 degC as "what roots see", which was
`tsod` -- the layer 6 to 12.4 m down -- and was additionally averaged over 31
cells the climate model calls ocean, where the soil array is fill. What protects
roots on Earth is snow, and the cap carries a deepest monthly cover of
**0.023 m**. Aerial tissue could survive this winter; roots could not, and
nothing in the model represents that separately, because `tcmin_surv` is one
number for a whole plant.

Larcher's section 4 bounds what is physically available, with qualifiers that
have to travel with it. The deepest hardening comes from a month or more of slow
cooling at -10 to -60 degC, and willow, birch, pine and black currant hardened
that way survive -196 degC -- but the organ is TWIGS, the protocol is
ARTIFICIAL, and he says the state "may not commonly be attained in nature". So
extreme freezing tolerance is available to woody TISSUE and the -69 degC air
temperature is not a physical barrier. It says nothing about roots.

Why roots are the exception is his Fig. 8, which tracks bud and root hardening
on separate courses in *Acer saccharum*: soil at 20 cm stayed between 0 and
-7 degC while the air fell below -20. **Earth's roots are buffered by soil and
snow and are therefore never selected for deep hardiness.** That buffer is what
this cap does not have -- a root zone near -79.5 degC under 0.023 m of snow,
against Earth's 0 to -7 -- so the root threshold binds
here for a reason no Earth analogue has been selected against. Deep supercooling alone caps at -30 to -50 degC before
homogeneous nucleation, and beyond that survival is by freezing TOLERANCE --
extracellular ice and cellular dehydration. Two mechanisms, not one scale, and
LPJ carries one number for both.

## The temperature response is a contributing suppression, not the binding one

BNS is the only tree that survives the polar winter, and its occupancy tracks
its Earth-calibrated optimum. Mean BNS foliar projective cover by the
warmest-month temperature of the cell:

| warmest month | cells | BNS | C3G |
| --- | --- | --- | --- |
| 0-15 degC | 442 | 0.0015 | 0.1008 |
| 15-20 | 235 | **0.0311** | 0.1005 |
| 20-25 | 347 | 0.0260 | 0.0508 |
| 25-30 | 227 | 0.0197 | 0.0530 |
| 30-35 | 143 | 0.0186 | 0.0449 |
| 35-40 | 83 | 0.0019 | 0.0287 |

BNS peaks at 15-20 degC and declines monotonically above its `pstemp_high` of
25; the median cell it occupies has a warmest month of 23.3 degC, just under
that plateau's top. The polar summer at 31.11 degC sits on the declining limb.

**But the size of that effect is smaller than the polar deficit.** BNS reaches
0.0186 in OTHER cells with a 30-35 degC warmest month, against 0.00265 at the
cap, a factor of seven apart. So the temperature response contributes and does
not bind; what binds at the cap is water capture, measured separately in
`biosphere/notes/polar-cover-cold-filter-and-capture.md`. Stating it the other
way round would be checking the instrument against the wrong effect.

## What this audit does NOT establish

- **That an adapted flora would do better, or by how much.** Nothing here
  measures that. The sizing arms named in `world-orok` -- capture-side
  parameters moved, and `tcmin_surv` relaxed separately -- are what would, and
  neither has been run. Until they are, "Earth's set is wrong here" is a
  classification and not a quantity.
- **That the parameters are wrong in the ordinary sense.** Implicit-Earth is a
  diagnosis, not a defect report: each value has a sound derivation for the
  planet it was derived on. What the diagnosis says is that the derivation does
  not transfer, and that Earth's figure is a distance to report rather than a
  target to solve onto. `docs/src/reference/vocabulary.md` argues that.
- **Anything about the parameters outside the polar cap.** The allometry, fire
  resistance, establishment and storage parameters are in the class by
  construction and were not tested against any measured effect here. They are
  enumerated above so a later sweep starts from a list rather than from the
  file.
