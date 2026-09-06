# The three soil phosphorus inputs, and which of them this world can derive

This is a worldbuilding project. Vesper is an invented super-Earth around a
mid-K dwarf, and the subject is the vegetation model that simulates its
biosphere: the three per-cell soil phosphorus parameters LPJ-GUESS-CNP reads,
what each one is, and what this pipeline can put in it.

`vendor/lpj-guess/modules/soilinput.cpp` sets all three by hand, in both the
soil-code path and the mineral-texture path this world takes:

    soiltype.kplab = 0.010;
    soiltype.spmax = 0.145;
    soiltype.pwtr  = 0.000003;

**Those are not site defaults in the sense of being measured at a site. They are
Wang, Law and Pak (2010) Table 2's OXISOL row, verbatim**, converted from
gP/m2 to kgP/m2: `kplab` 10, `spmax` 145 and a P weathering rate of 0.003
gP/m2/yr are that row exactly, and no other row matches any of the three. The
comment above them reads `// FIXED FOR AMAZON FACE AT THE MOMENT`, and Amazonian
soils are Oxisols, so the fork shipped one soil order's parameters and applied
them planet-wide.

## The three are three different problems

| parameter | what it is | unit | standing |
| --- | --- | --- | --- |
| `pwtr` | phosphorus released from rock into the root zone per unit time | kgP/m2 per EARTH year | DERIVABLE, and derived below from the ANUT-2 absolute arm |
| `spmax` | the Langmuir maximum sorbed P, the sorption capacity of the column | kgP/m2 | BRACKETED, on Wang's own optimisation prior against a stock this project derives |
| `kplab` | the Langmuir half-saturation: the labile P at which sorbed P is half of `spmax` | kgP/m2 | BLOCKED. No reachable source constrains it independently |

**Wang's `kplab` and `spmax` are fitted, and their own table says so.** Table 2's
caption: "Parameters `kplab` and `spmax` are estimated in this study, P
weathering rate is prescribed in this study." Section 3.1 says how: "we
restricted the value of `spmax` to be between 50 and 100% of the total inorganic
soil P in the optimization. The value of `kplab` is consequently dependent on
the range of `spmax`". So both are the residual of a fit, keyed to Earth soil
orders, and adopting either is adopting a tuned value for a world that has no
soil order map. That is the reason this document does not convert Table 2.

Table 2 in full, for the comparison the derivations below are checked against:

| soil order | `kplab` gP/m2 | `spmax` gP/m2 | P weathering gP/m2/yr |
| --- | --- | --- | --- |
| Entisol | 64 | 50 | 0.05 |
| Inceptisol/Gelisol/Histosol | 65 | 77 | 0.05 |
| Aridisol/Andisol | 78 | 80 | 0.01 |
| Vertisol | 32 | 32 | 0.01 |
| Mollisol | 54 | 74 | 0.01 |
| Alfisol/Spodosol | 75 | 134 | 0.01 |
| Ultisol | 64 | 133 | 0.005 |
| Oxisol | 10 | 145 | 0.003 |

`kplab/spmax` falls from 1.28 on the least weathered order to 0.069 on the most,
which is the physical statement that a deeply weathered soil binds phosphorus
both harder and in greater quantity. That pattern is real and fitting a function
to eight fitted points would be fitting to a fit, so it is recorded and not used.

## `pwtr` is the same law this pipeline already runs

`somdynam.cpp:soilpadd` has three routes and the interesting two are the same
equation.

**The texture route**, which this world takes today, is
`daily_pwtr = soil.soiltype.pwtr / VESPER_EARTH_YEAR_DAYS`: a per-cell constant
rate, declared per EARTH year under `world-ith2` because rock weathering does not
know this world's orbit. It carries no runoff, no temperature and no moisture.

**The gridded route** is
`bi * (pcont / 100) * runoff * exp(-ea/R * (1/T - 1/284.15)) * shield / 1000`,
which is Hartmann et al. (2014)'s P release law with this model's own daily
runoff: `bi` the lithology's runoff-normalised major-element yield, `pcont` its
relative P content in percent, `shield` a soil-thickness factor.

**ANUT-2's absolute arm is that same law, offline.**
`pedology/scripts/phosphorus_budget.py` computes

    F_P = release_percent/100 * F_(Ca + Mg + Na + K + SiO2)

with the major-element flux as Meybeck (1987) Table 2C concentrations times the
climatology's runoff, and `release_percent` is
`pedology/config/pedogenesis.yaml`'s `phosphorus_release_relative`, which is
Hartmann's own Table A1-2 relative P content. So `pcont` and `bi` are two
declarations this project already holds, and the gridded route is where they
belong: it carries the runoff dependence prognostically where a single number
per cell cannot.

**The quantity is the right one for a root zone and the artifact says so.**
Meybeck's model is fluvial export from exorheic continents, and on land draining
to a closed basin the solute dissolves and then stays, so the concentration
describes what enters solution rather than what leaves. `phosphorus_budget.py`
reports the first deliberately, which is what `pwtr` wants.

### The field

    pwtr[cell] = F_P[cell] * (1 - pfixation[cell])

in kgP/m2 per Earth year, on the climatology's grid, matched to the model by
index and never by longitude.

`pfixation` is the soil map's fifteenth column, `andic_p_fixation` in
`build_soil.py`: the fraction of released P in that cell that andic material
takes out of circulation OVER AND ABOVE ordinary retention, from
`pedogenesis.yaml`'s `phosphate_retention_andic` 0.85 and
`phosphate_retention_vitric` 0.25, the first being Soil Taxonomy's own threshold
for an andic horizon. **The andic share and no more is exactly the right
correction here**, because the model carries ordinary sorption itself, downstream
of `pwtr`, in the Langmuir isotherm `spmax` and `kplab` parameterise. Applying a
whole retention model would double-count it; applying nothing would leave the one
part of retention the model has no term for.

### Two terms it does not carry, and neither is set to zero by default

**Soil shielding.** Hartmann applies a soil-thickness term that shields fresh
rock from water. Regolith depth exists here and the shielding FUNCTION does not:
no read source in this repository gives one as a response to a soil thickness, so
`shield` is 1.0 as a DECLARED ABSENCE and the field is an upper bound by whatever
the missing function is worth.

**Temperature.** `pwtr_ea` must be 0, and that is a consequence of where the
concentration comes from rather than a claim that dissolution is athermal.
Meybeck's Table 2C concentrations are observed means over rivers at the
temperatures those rivers had, so they already contain the temperature response
in the mean. Multiplying them by an Arrhenius factor normalised to 284.15 K
would count it twice. What would let the term be switched on is a
temperature-normalised concentration table, which this project does not have and
which would be a re-derivation of Meybeck rather than a parameter choice. The
absence is registered in `phosphorus_budget.py`'s `not_carried` block already,
where its wording is that a concentration-times-runoff flux has no temperature
term.

### What the derived field is worth against what it replaces

ANUT-2's land mean is 3.33 kgP/km2/yr, which is 0.00333 gP/m2/yr. The fork's
shipped 0.003 gP/m2/yr is Wang's Oxisol row, so the constant happens to sit
within a tenth of this world's land mean and carries none of its spatial
structure. What the derived field adds is that structure: the yield percentiles
`phosphorus_budget.py` reports run from the p10 to a maximum that ANUT-2's own
acceptance criterion holds below Hartmann and Moosdorf (2011)'s Japanese
390 kgP/km2/yr, with no cell above it. Against Wang's prescribed range across
Earth's soil orders, 0.003 to 0.05 gP/m2/yr, this world's land mean sits at the
bottom, where Earth's most weathered soils are.

## `spmax` is bracketed on the stock, not taken from the fit

Wang's own constraint is the one usable thing in the optimisation: `spmax`
between 50 and 100 per cent of the total inorganic soil P. That is a PRIOR and
not a measurement, and it is recorded here as one. What changes when it is
applied on this world is the quantity it is a fraction OF: ANUT-2 emits a
root-zone total-P stock per cell in kgP/m2, parent P content times bulk density
times the regolith depth cut at the land column property contract's rootable
base, so the bracket becomes per-cell and derived instead of per-soil-order and
fitted.

    spmax[cell] in [0.5, 1.0] * root_zone_stock[cell]

**Two things the bracket's upper end is loose about, stated because they are not
symmetric.** The stock is a TOTAL parent-derived P and Wang's fraction is of
total INORGANIC P, so the upper end is above what the prior intends by whatever
share has become organic; and the stock is itself an upper bound, being parent
content with no depletion history. Both push the same way, so the bracket's low
end is the defensible one and a run says which end it is on. `anut-10` owns the
initial mineral-P state and is where a depleted stock would come from.

## `kplab` is blocked, and the measurement that would unblock it is named

`kplab` is a half-saturation: the labile P at which sorbed P reaches half of
`spmax`. Nothing reachable constrains it. Wang fitted it jointly with `spmax` and
says so; his eight values are per Earth soil order; and this pipeline emits no
field that a published sorption isotherm is regressed on.

**The relation that would settle it exists, and it does not constrain the
half-saturation.** Gu et al. (2025) assemble 83 paired Langmuir maxima and
affinity constants from 16 publications against pH, soil organic matter, clay
and oxalate-extractable Fe and Al: Qmax runs 73 to 1865 mgP/kg with a mean of
548, K_L runs 0.006 to 3.7 L/mg with a mean of 0.41. On the independent test set
their linear model explains 52 per cent of Qmax and **12 per cent of K_L**, and
their gradient-boosted model 50 and 34 per cent. FeOX, AlOX and pH are the only
predictors the linear Qmax model retains; clay and pH carry what little K_L skill
there is. So the literature route is answered, and its answer is that a published
relation constrains the CAPACITY and not the half-saturation, and needs a field
this pipeline does not emit either way.

**And its Qmax is not Wang's `spmax`.** Gu's is a laboratory sorption maximum per
unit soil mass; Wang's is a land model's effective capacity per unit area, fitted
inside a model with a leaching flux and a strongly sorbed pool taking P out
beside it. Scaled at Gu's mean over 20 cm at a bulk density of 1500 kg/m3 it is
164 gP/m2 against Wang's fitted 32 to 145; over a metre it is 822. The 20 cm
figure landing inside Wang's range is a consistency check and not a depth
derivation: what it says is that the two quantities are the same order over a
plough layer and differ by a factor of six over a rooted column, so a lab
isotherm cannot be substituted for a model parameter without a depth argument
nobody here has made.

**What is left is the field, and it is one pedology emission.** The properties
that control phosphate sorption are iron and aluminium oxide content, allophane,
clay and pH, and of those the soil map emits clay, pH and the andic fraction. The
andic column is the closest thing here to a sorption measurement and it is the
wrong shape: Soil Taxonomy's phosphate retention is a single point on an isotherm
at one assay concentration, so it gives one equation in two unknowns and cannot
separate a capacity from a half-saturation. An oxalate-extractable Fe and Al
proxy from `pedology` would give Gu's Qmax model its predictors, and
`biosphere/notes/mineral-reactivity-contract.md` already records that same field
as missing for the other half of its argument, so one emission answers both. It
would still leave K_L where it is.

Nothing is gated on this today. `framework/parameters.cpp` refuses `ifplim 1` on
`PFRAC_LEAFTOSAP`'s level and shape, so no configuration of this model reads
`kplab` for anything.
