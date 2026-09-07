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

## Result: the traits work and the types are not polar

Two arms, both on 1253 retained years over 8600 of spin-up, npatch 5, root seed
20260828 matching the accepted run, differing only in the `ltor_max` bracket.

| | `lpj_fb69ecbb8b5d468fbbf6815a03afcbf9` | `lpj_19140f59581f4c8ab9bf20b0fe3f80c1` |
| --- | --- | --- |
| `ltor_max` | 1.0, the low end | 2.0, the high end |
| acceptance | **PASS** | FAIL |

The high arm is refused by the per-cell trend test on `anpp.out`, which
`world-4hlw` records as the half of the contract whose limit is calibrated
rather than derived and which errs toward refusing. Its end-to-end drift is
-0.269% against a 5% limit, so nothing below rests on it; the accepted arm is
the one scored, and the high arm is quoted only where the bracket's direction
is the point.

### The declared criterion is met, by a wide margin

Summed foliar projective cover over the 234 cells poleward of 75 degrees,
against the accepted run's 0.030 and the criterion fixed above:

| | polar FPC | VPE | VPP | verdict |
| --- | --- | --- | --- | --- |
| accepted run, twelve Earth types | 0.030 | -- | -- | -- |
| `ltor_max` 1.0, ACCEPTED | **0.350** | 0.3245 | 0.0199 | above 0.10: the traits matter |
| `ltor_max` 2.0 | 0.473 | 0.3842 | 0.0847 | above 0.10 |

**The ephemeral wins and the perennial persists.** At the accepted arm the ratio
is sixteen to one; at the high end of the bracket four and a half to one, so a
leaf-biased allocation favours the perennial relatively even as both grow. The
two Earth survivors are displaced: C3G falls from 0.0274 to 0.0020 and BNS from
0.0027 to 0.0017. Cover is flat to four figures across the whole retained
record in both arms.

### The criterion did not ask the question that matters, and the answer is no

Native share of the summed cover, by band, on the accepted arm and the high arm:

| band | accepted run FPC | arm FPC | VPE | VPP | native share |
| --- | --- | --- | --- | --- | --- |
| poleward of 75 | 0.030 | 0.350 / 0.473 | 0.325 / 0.384 | 0.020 / 0.085 | 98.4% / 99.2% |
| 45 to 60 | 0.152 | 0.467 / 0.557 | 0.332 / 0.348 | 0.072 / 0.145 | 86.6% / 88.5% |
| equatorward of 15 | 0.835 | 0.899 / 0.948 | 0.188 / 0.260 | 0.013 / 0.042 | 22.4% / 31.9% |

**The two types took the planet.** They hold a fifth to a third of the
TROPICAL cover and raise total cover there from 0.835 to 0.899 and 0.948. A
type derived for the polar cap that outcompetes tropical broadleaf evergreens
in the tropics is not a polar type; it is a super-competitor, and what it says
about the cap is confounded by what it says about everywhere.

The bracket rules out allocation as the cause. At `ltor_max` 1.0 the takeover
is as complete at the cap and only somewhat weaker in the tropics, so the
thermal band and the absence of an establishment ceiling carry it, not the
leaf-root split.

WHAT WAS MISSING IS THE COST EVERY SHIPPED COLD TYPE PAYS. LPJ-GUESS confines
a cold-climate type to cold climates with `tcmax_est`, the warmest coldest
month it can establish under: BNE and BINE -1, BNS -2, TeBS 6, IBS 7, TeNE 10,
TeBE 18.8 degC. The tropical types and both grasses declare no ceiling, and
`vesper_polar` inherits C3G's. So the derivation gave the types a photosynthesis
plateau from 15 to 46 degC, unlimited cold survival, and no place they may not
establish. A plant with no weakness wins everywhere, and that is what the run
reports.

The refutation criterion above tested whether the traits move polar cover and
never tested whether the type is polar. That is the criterion's defect and is
recorded as one: a second criterion belongs beside the first, on the native
share of cover OUTSIDE the diagnosis's region, and it is stated below for the
next arm rather than applied to this one after the fact.

### What the next arm carries, and the criterion it is judged on

`tcmax_est` for `vesper_polar`, DERIVED from where the diagnosis holds. The
region `notes/underoccupied-niches.md` identifies is poleward of 60 degrees,
527 cells, and its coldest month runs from -75.34 to -30.92 degC. A ceiling at
**-30.92** bars establishment in 919 of 1617 cells, 57% of the land, and admits
the whole of the region. Every shipped cold type is confined the same way; this
is the derivation those types have and these did not.

Judged on both criteria, fixed now: polar summed FPC above 0.10 as before, AND
native share of cover equatorward of 45 degrees below 5%. A type meeting the
first and failing the second is again a super-competitor and the ceiling was
set too warm; a type failing the first with the second met is confined to a
niche it cannot fill, and the traits do not do what the cap needs.
