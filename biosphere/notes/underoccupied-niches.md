# Where the shipped plant functional types leave this world's ground empty

This is a document about a simulated world. Vesper is a worldbuilding project,
the twelve plant functional types are parameter sets LPJ-GUESS reads, and every
comparison to Earth is a distance to report rather than a target.

`world-orok` found one niche the shipped types do not reach by looking at the
polar cap. This asks the question over the whole planet, by
`biosphere/scripts/pft_exclusion_sweep.py`, measured 2026-09-07 on
`lpj_1e6a2b9ca51a4eff9592992cad96677b`.

## The count of barred types is not the signature

A type absent from a cell is evidence of nothing on its own. Boreal types are
correctly absent from the tropics and tropical types from the caps, and in both
places the ground is full. Barred counts run 7.3 to 10.1 of twelve in every
band on this planet, including the band carrying the most cover.

What the polar diagnosis actually rested on is the second signature: **the cover
that remains runs outside the temperature band it declares.** That is what
generalises, and the barred count belongs beside it as context.

| band | cells | killed | barred | GDD5 | precip | FPC | above `pstemp_high` | above `pstemp_max` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| poleward of 75 | 222 | 10.00 | 10.02 | 892.5 | 79.2 | 0.032 | 31.8% | 7.5% |
| 60 to 75 | 305 | 9.99 | 10.10 | 734.1 | 111.5 | 0.095 | 20.7% | 1.2% |
| 45 to 60 | 226 | 7.90 | 9.88 | 268.7 | 294.9 | 0.208 | 0.0% | 0.0% |
| 30 to 45 | 276 | 4.90 | 7.81 | 1124.6 | 305.4 | 0.569 | 8.6% | 0.1% |
| 15 to 30 | 190 | 3.04 | 8.63 | 2722.9 | 242.4 | 0.498 | 16.8% | 0.2% |
| equatorward of 15 | 264 | 0.70 | 7.31 | 2915.6 | 631.9 | 0.839 | 1.7% | 0.0% |

The two shares are cover-weighted over the resident types' own declared limits:
the time the vegetation spends where its photosynthesis is already declining,
and the time it spends where photosynthesis stops entirely.

## The confirmed niche is polar, and it is larger than the cap

Poleward of 60, 10 of 12 types are killed outright by `tcmin_surv` and the two
survivors run above their own ceiling for a fifth to a third of the growing
season. Of the 77 cells losing more than 5% of the season to `pstemp_max`
shutdown, 71 lie poleward of 60; their mean warmest month is 45.0 degC and their
mean cover is 0.014.

**So the diagnosis in `biosphere/notes/vesperian-polar-type.md` covers 527 cells
and not the 234 it was derived on.** The 60-to-75 band carries the same
structure at lower intensity: the same ten types killed, 20.7% of the season
above `pstemp_high`, 1.2% above `pstemp_max`, and cover of 0.095 against the
cap's 0.032. Any arm testing a derived polar type should report both bands,
because the trait values were derived on the harsher one.

## The 45-to-60 band is a different case and is NOT the same finding

It is the obliquity's cold trough: 268.7 GDD5, a quarter of the cap's, against
294.9 mm of precipitation, which is nearly four times the cap's water. It reads
**0.0% on both thermal shares**, so its vegetation is inside its declared band
and the polar diagnosis does not apply.

What it has instead is water it does not use, and 7.90 types killed with 9.88
barred, the extra two on `gdd5min_est` rather than on cold survival. Whether an
adapted type would take that water is NOT established here, and one measurement
that looks like evidence is not: the band's capture per unit leaf area ranks at
percentile 18, the planet's worst, but that statistic is annual and does not
normalise for how long leaf is present. A band whose season is a fifth as long
will rank low on it whatever its leaf does per day. Settling this needs a
season-normalised capture measure, which does not exist yet.

## What holds the ground that IS occupied

Grass carries 27.9% of the planet's cover, which is an unremarkable share.
There are 325 grass-dominated cells against 714 tree-dominated, and grass
dominance appears in every band at 11% to 27% of its cells.

| band | grass-dominated cells | share of the band | grass FPC there | tree FPC there |
| --- | --- | --- | --- | --- |
| poleward of 75 | 39 | 16.7% | 0.143 | 0.009 |
| 60 to 75 | 92 | 27.1% | 0.144 | 0.014 |
| 45 to 60 | 66 | 21.2% | 0.299 | 0.017 |
| 30 to 45 | 59 | 21.2% | 0.512 | 0.086 |
| 15 to 30 | 39 | 20.5% | 0.332 | 0.071 |
| equatorward of 15 | 30 | 11.4% | 0.647 | 0.136 |

**This world's grassland inverts Earth's.** Earth puts its great grasslands in
warm semi-arid mid-latitudes, held open against trees by water and by fire.
Vesper's densest grass-dominated ground is TROPICAL, at 0.647 cover, and its
mid-latitude grass sits in the 45-to-60 cold trough at 0.299, thin because the
obliquity denies it degree-days rather than because anything denies it water:
that band receives 294.9 mm, more than any band outside the tropics.

The planet is sparsely vegetated in general, and that is the right reading of a
thin-looking map rather than any shortage of grass: total cover over all land
cells runs p25 0.005, median 0.238, p75 0.680, so a quarter of land cells carry
essentially nothing.

## What this rules out

- **Not a second thermal-ceiling niche.** The 6 cells outside the polar bands
  losing more than 5% of the season to shutdown are scattered across two bands
  and carry no common structure.
- **Not a grass deficit.** Grass's share of cover is ordinary and its dominance
  is planet-wide; what is unusual is where its densest stands are.
- **Not an establishment-limit artefact.** The 16 cells carrying cover under a
  violated establishment limit are types that recruited under conditions since
  departed, which is the model behaving correctly; only `tcmin_surv` kills.
