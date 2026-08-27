# The interval the carve verdict's open-water evaporation is evaluated over

Vesper's closed basins. Measured on 2026-08-27 against `canonical-10m-base`
(8,772 basins, `selectionSource: threshold`, so pre-carve and never carved)
forced by `bootstrap_regular_climatology.nc`, a T21 bootstrap.

Found while fixing the same defect one file over; the lake balance's version is
in `hydrography/notes/lake-balance-integration.md`.

## What was wrong

`carve_verdict.py:main()` read every climatology field through `annual_mean`
and called `reference_level_air` with no bin index, then evaluated Penman once
on that annual-mean air. `export_carve_list.py:climate_terms` did the same.
Penman is nonlinear in everything it reads, which is already the stated reason
it integrates the DIURNAL cycle rather than reading a daily mean; the seasonal
cycle is the same argument at a longer period, and
`config/land_water_ledger.yaml` holds `open_water_evaporation` at
`interval_floor: climatology_bin`. So the annual evaluation was against a
standing decision rather than a simplification anyone had chosen.

## The prediction, stated before measuring

A basin carves when

    (E - P) / runoff  <=  catchment / area_at_spill - 1

and the right-hand side is pure geometry that does not move. Raising `E` raises
the left-hand side, so the inequality is satisfied for strictly fewer basins.
**The corrected carve set had to be a SUBSET of the old one**: basins leaving,
none entering. A basin entering would have meant the fix went in backwards, and
a mixed set would have meant something other than the evaporation had moved.

Measured: 4,772 to 3,944. **828 left the carve set and 0 entered**, subset
confirmed on the basin ids rather than on the counts. That is 9.4% of the
catalogue and 17.4% of the old carve set.

## Why the per-bin evaluation is not simply the right answer

The bin mean of the per-bin evaluation is 1.569x the single annual evaluation
over land, 1.198x over ocean. That is large enough to check before believing,
and the check does not come out one-sided.

**The magnitude is not any single input.** Taking one input per bin at a time
and leaving the rest annual, over land: temperature 1.106, humidity 1.099,
radiation 1.008, wind 0.998, diurnal range 0.999, and all of temperature,
humidity and pressure together 1.127. None of them, and no product of them,
reaches 1.569. The effect is an interaction, and the interacting parts are
Penman's two hard clamps.

**Per bin is right because of rectification.** `max(net_radiation, 0)` and
`max(es_a - e_air, 0)` are rectifiers, because a lake does not evaporate a
negative amount in a dark month. Averaging the air over a year first cancels a
real summer against a winter that never physically happens. Over this
climatology the bin mean of `max(net_radiation, 0)` is 86.4 W/m2 over land
against 77.0 at the annual mean, and 684 of 1,019 land cells go negative in
some bin while staying positive in the annual mean.

**Annual is right because Penman carries no heat storage.** It closes the
surface energy budget instant by instant, so per bin it charges the water for
bright-season energy that actually went into the water column and comes back
out in the dark season. That term integrates to zero over a periodic year and
to nothing like zero over a bin.

**The ocean measures the second effect**, because it is the one surface whose
evaporation the model already computes with the true storage in it:

| | ocean mean, mm/day | against the model |
| --- | --- | --- |
| the model's own evaporation | 1.770 | -- |
| Penman on annual-mean air | 1.685 | 0.952 |
| Penman as the bin mean | 2.019 | 1.141 |

and bin by bin the per-bin estimate runs 1.43 to 1.46 times the model in the
three brightest bins and 0.89 to 0.96 in the two dimmest, ranging 0.887 to
1.463 over the twelve. That is the signature of a missing storage term and not
a constant bias: a constant bias would have shown as the same ratio in every
bin, and would have meant the annual agreement was luck.

## What is done about it

**The interval is bracketed rather than chosen**, because the two ends are
interpretable and the spread between them is a quantity this project does not
have:

- the **annual evaluation** is the limit for a lake deep enough to hold its
  temperature through the year;
- the **bin mean** is the limit for a lake with no heat capacity at all;
- what sits between them is the water body's own heat storage, which needs a
  lake depth and a mixed-layer model. Declared, not estimated.

An ocean stores far more heat than any lake, so the 14.1% is an UPPER bound on
the bias at the bin-mean end, and a shallow playa carries almost none of it.

Both ends are bounds alongside `wet`, and the carve list is the intersection --
which is what the file already did for the two evaporation estimators and what
`export_carve_list.py` already did for the two climate arms, by taking the
larger retain. Ties break toward not carving: a basin carved in error has lost
a depression from the terrain everything else is built on, and a basin left
uncarved is still there to carve next cycle.

| | carves |
| --- | --- |
| `carve_verdict.py`, bin mean | 3,944 |
| `carve_verdict.py`, annual evaluation | 4,772 |
| `carve_verdict.py`, wet bound | 6,260 |
| **the carve list, intersection of all three** | **3,944** |
| `export_carve_list.py`, bin mean | 3,967 |
| `export_carve_list.py`, annual evaluation | 5,337 |
| **the exported list, intersection** | **3,967** |

The two scripts differ by 23 basins because they are two formulations of one
criterion -- the ratio form and the discharge form -- which the export already
documents as not nesting. The interval bracket is 25.7% of the union in the
export's discharge form.

## Two other defects found in the same files

**The cold arm was the warm arm.** `climate_terms(clim_path, args, config,
basins)` read `args.climatology` in both places it opens a climate and never
touched `clim_path`, so `climate_terms(args.endmember_climatology, ...)`
evaluated the primary climatology. The two-climate bracket was one climate
compared with itself, and it did not fail: `cut_endmember` equalled
`cut_primary`, so the reported disagreement was 0 and the bracket width 0.0% of
the union, which reads as the warm and cold ends agreeing completely. `main`
now refuses two arms whose discharge fields come out bit-identical. It cannot be
exercised on real data yet: the tree holds one regular climatology.

**The rung guard was unreachable in practice.** `require_configured_grid` was
reached only through `climatology_path()`, and both scripts call that ONLY when
`--climatology` is absent -- so every re-take, every arm and every sensitivity
skipped it. The rung appears nowhere in a climatology's name and the carve list
header takes its resolution from `config/planet.yaml`, so a sidecar could have
recorded T21 beside a T85 climatology with nothing in the tree contradicting
itself. The guard now runs on whatever climatology is in force, the endmember
arm included.

## What this costs

Nothing on the active lineage. `canonical-10m-base` carries
`selectionSource: threshold`: it is pre-carve and no carve list has ever been
applied to it, so the instrument is corrected before its first use here.

The four builds that were carved from a verdict -- `carved-zoned`,
`carved-zoned-v2`, `carved-zoned-v4` and `carved-zoned-v5`, all
`selectionSource: preserve-list` -- were taken with the annual evaluation and
with the endmember arm silently equal to the primary. Each therefore over-carved
by the direction argued above, and each preserved fewer basins than the
threshold rule does. All four are archived, payload deleted, and all predate
this lineage. Under CLAUDE.md rule 7 nothing is owed to them: a build is
disposable until the canonical climatology lineage is declared, and it is not.
