# The nitrogen closure tolerance was finer than the column it differenced

Found 2026-08-31, on run `lpj_de8a0c8bd4324639b4bd62b66867d23f`: the first
LPJ-GUESS run of the simulated Vesper biosphere to pass every equilibrium check,
which then refused on nitrogen conservation in four of 1617 gridcells.

Worldbuilding. Vesper is an invented planet; every nitrogen quantity below is a
modelled soil or vegetation quantity in a simulation of it, and the tolerances
are this project's own acceptance rule rather than anything measured.

## What was true

`biosphere/config/lpj_acceptance.yaml` closes cumulative nitrogen over ten
complete forcing cycles as

    residual = (npool.out Total[last] - npool.out Total[first]) * 1e4
             + sum(nflux.out NEE over the window)

and fails a gridcell whose residual exceeds `absolute_floor_kg_n_ha`, 2.0.

`vendor/lpj-guess/modules/commonoutput.cpp` wrote `npool.out` Total to FOUR
decimal places in kgN/m2. The rule converts that column by 1e4. **One
least-significant digit of the differenced column was therefore 1.0 kgN/ha.**
Each endpoint carries up to half a written quantum, so their difference carries
a whole one: a model conserving nitrogen exactly still reported up to
1.0 kgN/ha, against a tolerance of 2.0.

Measured on the run, over all three closure elements:

| element | quantum of the differenced column | residual an exact model reports | declared floor | margin |
| --- | --- | --- | --- | --- |
| nitrogen | 1.0e-4 kgN/m2, converted by 1e4 | 1.045 kgN/ha | 2.0 kgN/ha | 1.91 |
| water | 0.01 mm AET, 0.1 mm runoff | 0.55 mm | 2.0 mm | 3.64 |
| carbon | 1.0e-3 kgC/m2 | 0.001045 kgC/m2 | 0.01 kgC/m2 | 9.57 |

Nitrogen was the acute case and the only inverted one: its tolerance sat at
under twice its own quantisation, so the floor did all the discriminating and
what it discriminated on was rounding. Carbon and water were thin for the same
reason rather than wrong.

The release concedes the point itself. `commonoutput.cpp` carries a
`RUN_BENCHMARKS` switch that adds two decimals to these same columns, which is
mainline saying the shipped precision is too coarse for a budget check; two is
not enough for the nitrogen conversion.

## What the four refused gridcells were

Full-precision residuals over the contract's ten-cycle window, and the same
gridcells' residuals over longer windows ending at the same final cycle. A real
per-cycle leak accumulates linearly with the window; endpoint quantisation does
not depend on window length at all.

| lon | lat | L=10 | L=25 | L=100 | L=500 | L=1252 |
| --- | --- | --- | --- | --- | --- | --- |
| -174.38 | -13.84 | +2.43 | -0.21 | -0.39 | -0.09 | -0.43 |
| -135.00 | -19.38 | -4.95 | -4.61 | -5.98 | -5.50 | -5.75 |
| -28.12 | 8.31 | +2.03 | +0.42 | -0.54 | -0.05 | +0.32 |
| 0.00 | -13.84 | -2.43 | -3.29 | -3.42 | -3.57 | -3.62 |

Two of the four are pure endpoint quantisation: their residual collapses to
well under 1.0 kgN/ha at every window but the one the contract happened to
choose. The other two carry a flat offset of a few kgN/ha that does not grow
with the window, so it is not a leak either.

The whole population says the same thing. Per-gridcell-year residuals
`(pool[y] - pool[y-1]) * 1e4 + NEE[y]` over all 2,024,484 gridcell-years:
99.885 per cent fall within 1.0 kgN/ha, which is exactly the endpoint
quantisation bound, and the population |residual| over the ten-cycle window has
its shoulder at q0.99 = 0.96 with only 12 of 1617 gridcells above 1.0. The
worst gridcell's mean per-cycle residual is -0.0046 kgN/ha, which over the whole
1252-cycle record accumulates to -5.75 kgN/ha and is 0.035 per cent of that
gridcell's nitrogen throughput -- inside the contract's own
`relative_throughput_limit` of 0.001 by a factor of three.

The gridcell-years that do exceed the quantisation bound associate with fire and
not with the soil nitrogen operator: 71.7 per cent of them carry a fire gas
emission against a 44.4 per cent base rate, and their mean fire gas nitrogen is
2.51 kgN/ha against 0.87 over all gridcell-years.

## What the symmetry was

Two of the four carried residuals of exactly +2.43 and -2.43 at the same
latitude, one of them at longitude 0.00, which is the shape a coordinate
convention manufactures and CLAUDE.md rule 3 exists for. It was neither.

The two are distinct simulated gridcells: their `npool.out` and `nflux.out`
series differ in all 1253 retained cycles, `npool.out` holds 1617 distinct
(lon, lat) keys for 1617 simulated gridcells, and longitude 0.00 sits in the
regular T21 sequence beside -5.62 and 5.62. The tie is arithmetic. Residuals are
written to hundredths, so the 1617 gridcells take only 110 distinct magnitudes
and share 25,685 tied pairs; the tie also does not survive a change of window,
becoming -0.21 against -3.29 at twenty-five cycles.

## Whether the soil nitrogen operator's divergences are implicated

They are not, on two independent grounds.

**By construction.** `commonoutput.cpp:1464` accumulates `flux_ntot` from all
four soil gas terms the operator produces -- NH3_SOIL, NO_SOIL, N2O_SOIL,
N2_SOIL -- alongside the four fire terms, and `nflux.out` NEE is built from
`flux_ntot`. Every one of those terms is inside the budget the closure check
differences, so whatever the divergences do to their magnitude, they move
nitrogen from a pool to a flux the check already counts. They cannot open a
hole in it.

**By size.** The four soil gas terms together average 0.105 kgN/ha per cycle
over the last hundred cycles, against 0.874 for the fire terms and a total
nitrogen input of 2.199. Producing the -4.95 kgN/ha residual over ten cycles
from an unbudgeted soil gas would take about 0.5 kgN/ha per cycle, which is five
times the entire soil gas flux. Even a term omitted from the budget outright
could not reach it.

## What the check was replaced by

`assess_lpj_run.py:written_quantum` measures each closure column's decimal step
from the output table it read, `closure_resolution` forms the residual an
exactly conserving model would report from those quanta, and the gate refuses a
tolerance under `minimum_resolution_margin` times it before evaluating a single
gridcell. The quantum is measured from the artifact and not read from the model
source, so raising a column's precision and forgetting to record it still gives
the honest margin, and dropping one back refuses the next run.

`commonoutput.cpp` now writes the six columns the closure differences at a
precision that can resolve their tolerances. That is a change to the written
representation only; no simulated quantity moves.
