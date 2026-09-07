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
