# A Vesperian polar plant functional type: the derivation, and what would refute it

This is a document about a simulated world. Vesper is a worldbuilding project,
the plant functional types below are parameter sets for LPJ-GUESS, and every
comparison to Earth is a distance to report rather than a target to solve onto.
It is a CANDIDATE and not a decision: `docs/src/reference/design-intent.md` is
where a standing decision lands, and this note is the evidence one would rest
on.

`world-orok` establishes that the twelve shipped types are Earth's and that the
polar cap is where that bites. `biosphere/notes/polar-cover-cold-filter-and-
capture.md` closes the sizing question negatively for every capture-side trait
proposed: the cap converts leaf area into captured water at percentile 72 of
this world's cells, so it is not inefficient, and `lai.out` carries the leaf
area at FULL display, so a faster phenology shows the same 0.067 sooner. What
is short is leaf carbon, in a loop that limits itself.

This note derives what a type shaped by THIS world's selective pressures would
carry, from those pressures as measured, and states in advance what result
would refute it.

## The pressures, each measured

Measured 2026-09-07 on `lpj_1e6a2b9ca51a4eff9592992cad96677b` by
`biosphere/scripts/polar_water_timing.py` and `polar_capture_efficiency.py`,
over the 234 cells poleward of 75 degrees.

| pressure | measurement |
| --- | --- |
| the water arrives in one pulse | 83.5% of annual bare-soil evaporation leaves in one interval of 15.08 days; transpiration is spread at 45% |
| the pulse is IN season | 99.6% of that evaporation is in months above freezing, 0.0% before the thaw month; thaw and the 5 degC threshold arrive 0.27 days apart |
| the season is short | 3.39 months above freezing and 2.87 above 5 degC, of twelve |
| the season is HOT | growing-season air temperature p5 8.46, median 22.76, p95 45.70, maximum 48.43 degC; duration-weighted mean 24.16 |
| the winter is lethal to tissue | coldest month -68.75 degC, below every finite `tcmin_surv` in the shipped set |
| capture per unit leaf is already good | the cap's median extinction coefficient ranks at percentile 72 of the planet's cells |
| allocation runs away from leaf | `growth.cpp:1379` sets `ltor = min(wscal_mean, npscal) * ltor_max`, and the cap holds 3.36 mm of liquid soil water against 64.09 mm elsewhere |

## The anticorrelation, in its sharpest form

The cap needs a thermal ceiling the shipped set puts only on types that cannot
survive its winter. Against C3G's declared band, the cap spends **30.5% of its
growing season above `pstemp_high` 30**, where the type's photosynthesis is
already declining, and **7.6% above `pstemp_max` 45, where it stops entirely**.
The south cap's growing season averages 32.40 degC and its warmest month
averages 43.87.

C4G carries the ceiling that regime wants, `pstemp_high` 45 and `pstemp_max`
55, and is excluded from the cap by a `tcmin_surv` of 15.5. C3G carries the
cold survival and a temperate ceiling. **No shipped type carries both, because
on Earth nothing selects for both**: Earth's cold places are low-energy places.
Vesper's cap is a cold place that receives 390.67 W/m2 in its warmest month
against the tropics' 211.84.

## The traits, and how each is fixed

Three dispositions are used, and each parameter says which: DERIVED from a
measurement above, BRACKETED for a sweep because it cannot be derived, or
INHERITED from C3G deliberately.

| parameter | value | disposition |
| --- | --- | --- |
| `tcmin_surv`, `tcmin_est`, `tcmax_est`, `twmin_est` | no limit | INHERITED from C3G: the type survives the winter below ground rather than as tissue, which is what an absent limit means here |
| `leaflong` | 0.2392 simulation years | DERIVED: the leaf lives the 2.87-month growing season and no longer |
| `sla` | 44.8 m2/kgC | DERIVED, and not free: `ifcalcsla 1` computes it from `leaflong` through the Reich regression in `guess.h:initsla()`. C3G's 0.998 simulation years gives 26.0, so the shorter leaf buys 1.72 times the leaf area per unit carbon |
| `turnover_leaf` | 1 | DERIVED: one growing season per orbit, one leaf cohort, shed complete |
| `pstemp_low` | 15 | DERIVED: the plateau opens at the growing season's p25 of 15.30 degC |
| `pstemp_high` | 46 | DERIVED: the plateau closes at its p95 of 45.70 degC |
| `pstemp_max` | 55 | DERIVED: above the season's measured maximum of 48.43 degC, so photosynthesis is not switched off inside the season |
| `pstemp_min` | 0 | DERIVED: below the season's p5 of 8.46 degC and above the frozen months, which carry no liquid water to work with |
| `ltor_max` | 1.0 to 2.0 | BRACKETED and swept. C3G's 0.5 encodes a long-season foraging strategy: at least two units of root per unit leaf, and water stress scales it further down. Where the extractable water passes in fifteen days, that investment cannot be recovered, but no measurement here fixes the replacement |
| `phengdd5ramp` | 100 | INHERITED from C3G, so the arm moves leaf economics and not display timing. It is not neutral once potential leaf area rises, and is the second thing to sweep |
| `rootdist`, `root_beta` | C3G's | INHERITED: 90% of roots already sit in the top five layers, and the note shows pulse depth reaches the root zone |

## What the model cannot express, and what that costs the test

**LPJ-GUESS has no plant carbon reserve.** `nstore_longterm` and
`nstore_labile` are nitrogen; the only `calculate_carbon_store` is in
`soilmethane.cpp` and is a soil pool for methane. Leaf carbon is built from
concurrent assimilation.

Earth's answer to a fifteen-day pulse at the head of a short season is the
perennial storage organ: a desert geophyte deploys leaf on the PREVIOUS
season's stored carbon, so the canopy is present when the water arrives rather
than being funded by it. That is the trait most likely to break the loop here,
and it is a fork change rather than a parameter block.

So the arm below tests the parameter-expressible half only. **A null result
refutes the leaf-economics and thermal-ceiling traits. It does not refute the
storage hypothesis**, which would need the reserve to exist before it could be
tried.

## What would refute this, fixed before the run

The type is added ALONGSIDE the twelve rather than replacing any, so the arm
answers whether a type shaped by this world's pressures takes ground the
shipped set leaves bare. Judged on summed foliar projective cover over the same
234 cells, against the accepted run's 0.030:

| outcome | reading |
| --- | --- |
| polar FPC above 0.10 | the leaf-economics and thermal traits matter, and the cap's bareness is an artefact of the shipped parameter set |
| polar FPC below 0.05 | REFUTED. The parameter-expressible traits do not break the loop, and the cap is bare for reasons no PFT block reaches |
| between | partial; report the number and the bracket that produced it |

The `ltor_max` bracket is swept across both ends, and a result that holds at one
end and not the other is reported as the bracket rather than as its better half.

Two further readings are recorded whatever the verdict, because they are what a
partial result would be diagnosed with: the share of the polar growing season
the new type spends above its own `pstemp_high`, which says whether the derived
ceiling was set high enough, and the capture fraction, which says whether the
extra leaf reached the water.
