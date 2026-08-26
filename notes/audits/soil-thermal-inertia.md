# The modelled soil's thermal inertia is set by water, not by rock

Worldbuilding. Vesper is an invented planet. Everything below is about the
simulation of its land column: `landmod.f90`'s soil heat solver, and a
published soil thermal model applied to the materials this world's land is
made of.

## The finding this tests

`landmod.f90:tands` builds the soil column's properties per layer as

    zcap(:,jlev)  = sicecap*dglac(:)  + soilcap*(1.-dglac(:))
    zdiff(:,jlev) = sicediff*dglac(:) + soildiff*(1.-dglac(:))

so outside glacier ice the whole simulated planet carries one heat capacity and
one thermal conductivity. Their thermal inertia is the quantity the diurnal and
seasonal surface temperature amplitude is inversely proportional to, and that
amplitude reaches evaporation, P minus E, the carve criterion and the dust
emission threshold.

The recorded finding was that the barren classes -- evaporite salt crust and
playa clastics, about a third of this world's land -- carry roughly a third of
the model's inertia and are damped several times too much, and that the fix is
two per-cell fields out of the LITHOLOGY map: two NHOR arrays, two surface
codes with their `surfcode` rows, two restart policy entries, a generator with
its pipeline row, and a literature table of conductivity, density and heat
capacity per rock class.

The magnitude is right. The attribution is not, and the deliverable it implies
would write a uniform field over exactly the cells it was specified for.

## What the standard soil thermal model says before it is run

Johansen's interpolation, read from Farouki (1981) CRREL Monograph 81-1,
Section 7.11 and Table 24 on pages 112 and 113, and the formulation
`lawrence_2007` blends organic material into:

    k_dry   = (0.135*rho_b + 64.7) / (2700 - 0.947*rho_b)     +/- 20 per cent
    k_sat   = k_solid^(1-phi) * k_water^phi
    k       = Ke*(k_sat - k_dry) + k_dry
    k_solid = k_quartz^q * k_other^(1-q)

Mineralogy enters in one place, `k_solid`, and `k_solid` enters in one place,
`k_sat`. The conductivity of dry soil is a function of bulk density and of
nothing else. So the parent material's effect on thermal inertia is largest
where the soil is wet and is identically zero where it is dry -- and the
classes the finding names are the dry ones.

**Two of Farouki's relations are texture-branched**, and the branch is what a
restatement loses. `Ke` is `log10(Sr) + 1` for a fine soil and
`0.7*log10(Sr) + 1` for a coarse one; `k_other` is 2.0 W/m/K in general and
3.0 W/m/K for a COARSE soil below 0.20 quartz, which is most of this world's
parent materials. Farouki's cut is his Figure 160's, a soil with more than 5 per
cent finer than 2 micrometres being fine, so the medium-textured column measured
below is fine and takes the general rows. The coarse branch is measured anyway,
because a coarse low-quartz soil is exactly the case where the mineralogical
term is largest.

**The dry end sits below the method's stated range.** Farouki gives the fine
relation for `Sr > 0.1` and the coarse one for `Sr > 0.05`, and says in 7.13.1
that below about 0.1 no method predicts well except Van Rooyen's. At and below
the floor the Kersten number is zero and `k` is `k_dry` exactly, so the answer
there is qualitatively right and quantitatively outside what Johansen claims.
Both points are reported.

## The measurement

`analysis/soil_thermal_inertia.py`, on a medium-textured mineral column at
1300 kg/m3, measured on 2026-08-25 and re-measured on 2026-08-26 against the
monograph itself. Thermal inertia in J/m2/K/s^0.5:

| Term swept | Range | Factor |
| --- | --- | --- |
| Saturation, air dry to saturated | 406 to 1992 | 4.91 |
| Saturation, from Farouki's own validity floor to saturated | 439 to 1992 | 4.54 |
| Bulk density at the dry end, 1100 to 1600 kg/m3 | 333 to 538 | 1.61 |
| Parent material at saturation, quartz 0 to 0.95 | 1807 to 2460 | 1.36 |
| Parent material at saturation, on a COARSE column | 1993 to 2460 | 1.23 |
| Parent material at the dry end, either column | 406 to 406 | 1.00 |

`landmod`'s single pair is 2078, which this column does not reach even
saturated: what the whole simulated planet carries is a wet, dense soil.

So the three terms rank moisture, then texture, then mineralogy, and the
mineralogical term is exactly nothing on the dry surfaces the finding is about.
The recorded 3.3x damping of the playa and evaporite classes is real and is a
MOISTURE difference wearing a lithology label: those surfaces are dry, and the
model gives them a saturated soil's inertia.

**Farouki's coarse branch makes the lithology term SMALLER, not larger.** It is
the one coefficient a restatement could have hidden that bears on this finding,
and it moves the wrong way for the deliverable it might have rescued: raising
`k_other` to 3.0 lifts the low-quartz classes by up to a tenth at saturation and
closes the gap to the quartz-rich ones, taking the saturated lithology spread
from 1.36 to 1.23. At the dry end it changes nothing at all, because it enters
`k_solid`, `k_solid` enters `k_sat`, and `k_sat` is multiplied by a Kersten
number that is zero there.

## What follows

**The lithology field should not be built.** Its own physics says it would
carry a factor of 1.36 where the soil is wet, 1.23 on a coarse column, and
nothing at all where it is dry, at the cost of two surface codes, two restart
policy entries and a per-class literature table this project does not hold. It
is the smallest of the three terms and it is absent on the cells it was
specified for. Reading the monograph rather than a restatement of it was the
test that could have overturned that, and it did not: the one branch the
restatement omitted moves the term down.

**The moisture dependence is the fix**, and it needs none of that. The model
already carries the soil water prognostically, as `dwatc` and as `dwatcl` under
the layered scheme, and the soil heat solver does not read it. No surface code,
no restart record, no generator and no new field: `zcap` and `zdiff` become
functions of state the column already has.

**It is blocked, and on a declaration rather than on effort.**
`pedology/config/land_column_properties.yaml` carries
`thermal.composition_dependent: undeclared`, owned by `lshy-5`, and the mapping
the change needs runs through it. `dwatc` is a bucket store in metres against
`dwmax`, a plant-available capacity; a degree of saturation is a fraction of
PORE volume. Turning one into the other needs a porosity and a retention
mapping that the land column property contract is the place for, and choosing
them here would be calibrating against numbers with no derivation.

**Texture, not lithology, is what a per-cell field would honestly carry** at
the dry end, and pedology's soil map is where texture lives. That is a smaller
term than moisture and it is second rather than third, so it is worth revisiting
after the moisture dependence exists and not before.

`world-jsfm` carries the moisture dependence; `world-yip` is closed on this
measurement, and `world-sy9`'s third acceptance clause is answered by it.
