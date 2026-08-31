# One mineral-reactivity contract, and why its second arm refuses

This is a worldbuilding project. Vesper is a simulated super-Earth around a
mid-K dwarf, and its biosphere is LPJ-GUESS-CNP. Everything below is about how
that model decides what fraction of decomposing carbon in the simulated soil
becomes stable organic matter, and what fraction of released phosphorus is
sorbed, and about the pedology artifact that feeds both.

The soil decomposition audit's finding 4 is that mineral protection is reduced
to Earth texture regressions while the pedology model already distinguishes
andic material and phosphate fixation. This document is the contract that puts
both arms in one place: what the texture-only arm is, what the mineral-aware arm
would need, and what crosses the interface today.

`biosphere/config/mineral_reactivity.yaml` is the declaration.
`biosphere/scripts/mineral_reactivity_gate.py` enforces it against
`modules/somdynam.cpp` and against the soil map itself, and it can fail.

Nothing here is verified by execution. LPJ-GUESS builds and runs on this tree,
and no run it has made has passed acceptance.

## The active arm is five linear functions of texture

Every mineral control on soil organic matter in the model is one of these, all
from Parton et al. (1993) and all read from `modules/somdynam.cpp`.

| what it controls | form | source |
| --- | --- | --- |
| the soil microbial pool's decay rate | `1 - 0.75*(clay + silt)` | eqn 5 |
| slow SOM to passive SOM | `max(0, 0.003 - 0.009*clay)` | Fig 1 |
| the microbial pool's respired share | `max(0, 0.85 - 0.68*(clay + silt))` | eqn 7 |
| the microbial pool to passive SOM | `0.003 + 0.032*clay` | eqn 9 |
| organic carbon leached out of the microbial pool | `0.03 + 0.12*sand`, scaled by percolation | eqn 8 |

Nothing else in the decomposition path is a function of the mineral phase. There
is no iron or aluminium oxide term, no allophane term, no aggregate capacity and
no polyvalent cation term, and the two columns pedology writes that bear on any
of it were read by nothing.

Parton et al. (1993) is read, and its calibration is what makes these Earth
central tendencies rather than laws: grassland sites, 0 to 20 cm soil stocks,
over Earth's range of textures.

## Two things the texture arm does that are worth knowing

**The clay control on passive SOM formation saturates, and it saturates inside
this world's range.** `max(0, 0.003 - 0.009*clay)` is exactly zero above a clay
fraction of 1/3, so on any soil past that the slow pool has no direct route to
passive SOM at all and reaches it only through the microbial pool. Measured on
`pedology/data/precarve-craton-10m/soilmap_T21.txt` on 2026-08-24: 248 of 1,019
land cells, 24.3 per cent, sit above that threshold. So a quarter of this world's
land is on the flat part of the one equation that is supposed to make clay-rich
soils protect carbon, and the andic contrast that would discriminate among those
cells is not read.

**The microbial partition can leave a negative remainder, and the range that
does it is not a soil.** What leaves the soil microbial pool is respiration plus
leaching plus the passive share plus the slow share, and the slow share is
computed as whatever is left. On pure sand at saturating percolation the first
three sum to 1.003, so the remainder is -0.003 and the transfer to slow SOM runs
backwards. Nothing in `somdynam.cpp` guards the sign.

The remainder form is not this fork's. Parton et al. (1993) eqn 10 on p. 788 is
`C_AS = (1 - C_AL - C_AP - F_t)`, a bare remainder with no sign guard either. But
the paper's own eqn 8 is `C_AL = (H2O_30/18)*(0.01 + 0.04*sand)`, and under those
coefficients the remainder is positive over the whole texture simplex, bottoming
out at pure sand. What admits a negative value is the CENTURY 5 leaching update
the code applies instead, `0.03 + 0.12*sand` saturating at 1.9 cm H2O per month
rather than 18: three times the coefficients over a ninth of the range. That
update arrived with the LPJ-GUESS subtree and is upstream, not a divergence this
project made.

Solving the remainder for where it reverses, at saturating percolation it is
`-0.003 + 0.8*(clay + silt) - 0.032*clay`, so a texture reverses only if its
clay-plus-silt fraction is below 0.0039 -- a sand fraction above 0.996 -- and
only if percolation is at least 0.98 of the leaching saturation point. Both have
to hold at once.

Neither is close to this world, and neither is close to anything the model can
be handed. This world's sandiest cell is 0.663 sand and the soil map's worst
remainder is 0.263, over its 1,019 land cells. The LPJ soil code table in
`soilinput.h`, which is what the model uses when no soil map is supplied, has
its coarse code at 0.90 sand and a worst remainder of 0.0754; the fixed-texture
fallback beside it gives 0.569. So the transfer runs forwards on every texture
any input path in this model can produce, by a margin of a factor of 90 in the
fine fraction, and `world-t67j` is not a defect on this world. The declaration
registers all four bounds, so a coefficient change that moves any of them fails
the gate.

One path does approach it. `somfluxes` reads the fine fractions through
`Soil::get_clayfrac` and `get_siltfrac`, which return the peat texture on a
high-latitude peatland stand, but reads the sand fraction for leaching straight
off `soiltype.sand_frac`, which is always the mineral one. The peat fractions
are never assigned and stay zero, so on such a stand the respired share is at
its 0.85 maximum while leaching still scales with mineral sand, and the coarse
soil code's remainder falls to 0.009. It reverses above 0.975 mineral sand,
which the soil code table's 0.90 does not reach. `run_peatland` is written as 0
whenever the wetland gate has not granted activation, so nothing runs this path
today.

## The mineral-aware arm refuses, and names four things

Cotrufo et al. (2013) and Lehmann and Kleber (2015), both read, put microbial
products, mineral association, aggregation and accessibility at the centre of
stable soil organic matter formation, and name the associating surfaces:
phyllosilicates, metal oxides, polyvalent cations, allophane, aggregates.
Neither supplies a parameterised transfer function, which is why replacing the
texture arm is not a coefficient swap. It needs state, and the state is not
there.

- **Iron and aluminium oxide content**, the standard measure of the reactive
  surface that sorbs phosphorus and stabilises carbon. Nothing produces it: the
  pedology model derives lithology fractions and a weathering intensity, and
  neither is an oxide content.
- **Allophane content.** This is the near miss. `build_soil.py`'s
  `andisol_properties` computes `andic`, the fraction of a gridcell whose
  volcanic glass has weathered to allophane under sufficient leaching. That is an
  AREAL fraction of andic material, not a concentration of allophane in the fine
  earth, and substituting one for the other needs a conversion nothing in this
  project supplies. Registering it as a partial producer rather than as the
  input is the whole point of the distinction.
- **Aggregate capacity.** Occlusion inside aggregates is physical protection,
  separate from mineral association, and nothing in the pipeline resolves soil
  structure.
- **Polyvalent cation saturation.** The soil map carries pH, from which base
  status can be argued but not derived. The non-N/P adequacy screen already
  refuses the trace set for want of a release table, and this is the same gap.

Each carries the `undeclared` sentinel and there is no default for any of them.
`--strict` refuses while any is undeclared.

## What crosses the interface today

`build_soil.py` writes `andic` and `pfixation` as the last two columns of the
soil map. `SoilInput::load_mineral_soils` now matches both by name, `get_mineral`
carries them into its `SoilProperties`, and both `get_soil` paths copy them onto
`Soiltype::andic_frac` and `Soiltype::p_fixation_frac`. An input path that has no
andic state -- every soil-code path, since an LPJ soil code is a texture class --
leaves `UNSET_SOIL_FRAC` rather than a zero, because "no andic material here" and
"nothing supplied one" are different claims and only the second should refuse.

**No equation reads either field**, and the gate checks that: it greps every
model translation unit outside the input module and fails on one that does. That
is the honest state of finding 4. The state now crosses the interface with its
sentinel discipline, so adding the mineral-aware arm is a change to
`somdynam.cpp` alone rather than a second soil-map interface, which is what
SDEC-4 asks for. What it is not yet is a mineral-aware model.

## Neither arm may be retuned

A coefficient in this declaration moves when its source moves, and for no other
reason. Fitting either arm so that simulated soil carbon or productivity comes
out where it is wanted would make every subsequent comparison circular, and the
prediction this component is scored against is registered in
`productivity-prediction.md` precisely so that it cannot be.

## Primary sources read

- Parton et al. (1993), *Observations and Modeling of Biomass and Soil Organic
  Matter Dynamics for the Grassland Biome Worldwide*, `10.1029/93GB02042`.
- Cotrufo et al. (2013), *The Microbial Efficiency-Matrix Stabilization (MEMS)
  framework integrates plant litter decomposition with soil organic matter
  stabilization*, `10.1111/gcb.12113`.
- Lehmann and Kleber (2015), *The contentious nature of soil organic matter*,
  `10.1038/nature16069`.
