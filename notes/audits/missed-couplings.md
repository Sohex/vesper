# Audit: couplings nobody priced, and quantities read at the wrong level

*Audited 2026-08-17, against the baseline climatology
`baseline_regular_climatology.nc` on `precarve-craton`. Read-only audit; nothing
here was fixed. The tasks these findings became are named at the end.*

Commissioned after an outside reader pointed out that dust suppresses global
precipitation, which cuts runoff, which is the denominator of the carve
criterion. Dust's radiative forcing had been priced carefully and its effect on
the water cycle had never been asked about. Nobody was wrong; a seam between two
correct components went unexamined.

This audit hunts that shape. It found one thing larger than the dust coupling,
one already-flagged item that turns out to be quantifiable and big, one piece of
today's work that is silently absent from the repository, and one structural
one-way coupling. It also knocked down a fifth that looked worse than it is.

Every number below is measured on the current baseline unless labelled
otherwise. Findings are tagged **[numeric]** where I computed them,
**[inspection]** where I read them out of code or artifacts, and **[physics]**
where I am reasoning without checking.

---

## 1. Runoff amplifies everything by six, and nothing is budgeted in runoff

**[numeric]** The land water balance on the current baseline:

| | mm per Earth year |
| --- | ---: |
| land precipitation | 816.4 |
| land evaporation | 690.0 |
| **runoff, the residual** | **126.4** |

Runoff is **15.5% of precipitation**. It is the difference of two numbers five to
six times its own size, so the amplification factors are

    d(runoff)/runoff  =  6.5 x  dP/P        (from precipitation alone)
                      =  5.5 x  dE/E        (from evaporation alone)

**A 10% error in land evaporation is 55% of the runoff.** A 10% error in
precipitation is 65% of it. Both feed the carve criterion directly, as its
denominator.

That much is arithmetic. The finding is what sits beside it:
`analysis/error_budget.json` prices **every** item in kelvin or W/m2, and **not
one** in runoff, evaporation or basins. Its own stopping rule -- "refine an input
when its plausible range exceeds the effect of the thing you last refined" --
is therefore being applied to the wrong output. The temperature is a result this
project can revise. The carve list is not: it leaves the project, changes the
terrain, and cannot be undone.

Using the parent session's own sensitivity table, which maps evaporation changes
onto basin counts, the existing budget items convert roughly as:

| budget item | its kelvin | plausible runoff effect | basins |
| --- | ---: | ---: | ---: |
| biosphere, bare rock vs vegetated | 4.28 K | tens of percent | hundreds |
| lakes composited into albedo | 0.61 K | several percent | tens |
| playa_clastic albedo 0.25-0.33 | -0.52 K | several percent | tens |

The middle column is **[physics]**, not measured -- converting a kelvin to a
runoff change needs the model's own hydrological sensitivity, which nobody has
measured on this world. That measurement does not exist and is the gap. One
perturbation run reporting dP and dE against dT would supply a conversion factor
usable for every item in the budget forever, exactly as the flux slope of 150.2 K
per unit flux ratio already does for temperature.

**Why this is the largest finding here.** The dust coupling that prompted the
audit is one instance of it. So is the biosphere. So is every albedo item. The
budget has one currency and the irreversible decision reads another.

---

## 2. "Penman over a dry land column" is the budget's only unquantified item, and it is two errors of very different size

**[numeric]** The budget lists this as *"unquantified, one-signed"*. It is
quantifiable from artifacts already in the repository.

`penman_open_water` forms its vapour pressure deficit as `es(tas) - e_air`. Over
a subgrid lake the surrounding column is dry, so `e_air` is too low, the deficit
too high, and evaporation overstated -- which closes basins and makes the verdict
under-carve. Raising the column toward what a lake-influenced boundary layer
would hold:

| assumed column relative humidity | land-mean Penman | change |
| --- | ---: | ---: |
| the model's own `q_air` (current) | 3.107 mm/day | -- |
| at least 70% | 2.810 | -10% |
| at least 80% | 2.691 | -13% |
| at least 90% | 2.553 | -18% |

Against the parent session's measured basin sensitivity -- 144 basins per 10% of
lake evaporation, 235 per 20% -- this reads as worth 144 to 235 basins, in the
direction of under-carving. **Both halves of that turned out to be inflated and
the item is worth 1 to 35 basins**; `notes/audits/carve-criterion-terms.md`
finding 3 measures it. The level mixing below accounts for most of the humidity
change once it is fixed properly, and the sensitivity it is multiplied by was
itself measured with a floor at the land rate in place, which inflates any
evaporation sensitivity by about a factor of four. That floor is finding 1 of the
same document and is much the larger item.

**Part of it is a level mixing, which is the archetype this audit was asked to
hunt, and it is the part that was worth something.** `es_a` was saturation at
`tas`, the **2 m** temperature. `e_air` came from `q_air` at the **lowest model
level**, whose height `penman_open_water` computes itself as
`z_ref = (R T / g) ln(1/sigma)`, of order 300 m. A standard Penman takes both at
one reference height. Taking the saturation term low and warm and the actual term
high and dry inflates the deficit systematically.

The wind and the aerodynamic resistance were mutually consistent -- both at the
lowest model level -- so this was specifically the humidity/temperature pair, and
the fix is to read all four at the level the transfer coefficient is derived over.
The 2 m air is 1.85 K warmer than the lowest level over land, worth 13.7% on
`e_s`. The ocean validation is what says this is a fix rather than a preference:
the error against the model's own open-water evaporation goes from 8.45% high to
3.28% low. `notes/audits/carve-criterion-terms.md` finding 2.

---

## 3. Today's error-budget updates are not in the repository

**[inspection]** `analysis/error_budget.json` currently states, of dust:

> `"magnitude": "-0.5 to -2 K, cooling everywhere"`
> `"note": "Aerosols are off entirely (L_AERO = 0), and ExoPlaSim 3.4.2 cannot
> switch them on: radmod.f90 declares aero_nl but never reads it..."`

Both claims were superseded today. DUST-2 priced the forcing at +0.35 to +0.74
W/m2 global mean with the sign surface-dependent and positive over bright basin
fill, and the parent established that the aerosol **can** be switched on --
`aero_ini` at `plasim.f90:190` reads `aero_nl` with `l_aerorad` and `aerofile`
use-associated from `radmod`, before `radini`'s broadcast at line 439.

Commits `86baba7` and `2153dcf` both describe updating this file. Neither
contains it:

    $ git log --oneline -1 -- analysis/error_budget.json
    ce15a69 Dust sign resolved: cooling everywhere, ...     # 2026-08-16

The working tree is clean and the file's content matches `scripts/error_budget.py`
line 122 exactly, so the current file is the generator's output.

**The hazard underneath is that this is a generated artifact that was hand
edited.** Its own note says *"Generated by scripts/error_budget.py"*. An edit to
the JSON is not an edit to the project: the next run of the generator discards
it, and in the meantime `git` may see no change at all, so the loss is silent
both ways. The content belongs in the generator.

This is close to failure-modes class 5, "provenance written once and never
recomputed", but not the same: there the record went stale, here the record was
corrected and the correction evaporated. **I would call it a sixteenth class:
editing the output of a generator, where the edit is discarded and the discarding
is invisible.** `world_state.json` is protected against exactly this by an
explicit "never edited, re-run the generator" rule in CLAUDE.md; the error budget
carries no such rule and is otherwise identical in kind.

---

## 4. The biosphere reaches the climate through albedo, and not through water

**[inspection]** `WORKFLOW.md` section 4 states the biosphere-climate coupling
entirely in radiative terms: *"Bare rock and a vegetated surface differ in land
albedo by enough to be worth several kelvin."* That is the channel that is
modelled, and it is priced at 4.28 K, the largest single item in the budget.

What reaches ExoPlaSim from the vegetation state is albedo (174-176), forest
fraction (212) and roughness (173). Roughness does affect evaporation through the
transfer coefficient, so the coupling is not zero. **What does not exist is any
stomatal or LAI control on evaporation, and any rooting depth.** ExoPlaSim's land
surface is one bucket whose wetness factor saturates at 40% of `dwmax`
(`landmod.f90:52-53`); a vegetated cell and a bare cell with the same bucket
state evaporate identically apart from roughness.

Meanwhile `pedology/README.md` records, from Lapides et al. (2024), that adding a
bedrock vadose zone to LPJ-GUESS raises annual transpiration by **100-150 mm** --
which against a runoff of 126 mm/yr is the entire residual. The biosphere
component knows transpiration is worth about one runoff, and none of it is
represented in the model that produces the runoff.

**[physics]** Direction is not obvious and I am not going to guess: stomatal
closure under stress reduces evaporation relative to a wet bucket, while deep
roots accessing water the bucket cannot hold increase it. The magnitude is what
matters, and by finding 1 it is amplified 5.5x into the carve.

This is a structural one-way coupling rather than an error. It deserves recording
in WORKFLOW section 4 beside the albedo channel, because the section currently
reads as though albedo were the whole of it.

---

## 5. Knocked down: the clamped aerodynamic term

**[numeric]** `penman_open_water` clamps its deficit with
`np.maximum(es_a - e_air, 0.0)`. Given the level mixing in finding 2, I expected
that clamp to fire widely and silently delete the aerodynamic term -- which would
have made the whole justification for choosing Penman inoperative.

It fires on **29.3% of land by area**, which looked alarming. It is not, because
of where:

| | fraction of cells clamped |
| --- | ---: |
| land warmer than 280 K | 0.7% |
| land at or below 280 K | 87.1% |

The clamp is a cold-cell phenomenon: at low `tas`, saturation at 2 m falls below
the absolute humidity carried at 300 m. **Closed basins sit in the arid, warm
part of the land**, where the clamp essentially never fires. So this does not
touch the carve verdict, and finding 2's overstatement stands on its own.

Recorded because the 29.3% figure is alarming out of context and someone will
find it again.

---

## What I would chase first

**Finding 1**, and specifically the one measurement inside it: a perturbation run
reporting dP and dE against dT on this world. It is the conversion factor that
turns the existing error budget -- which is already careful, already ordered, and
already in one currency -- into a budget that speaks about the carve. Without it,
findings 2 and 4 and the dust coupling can each be argued about separately and
none of them can be ranked against the others.

It is also cheap relative to what it settles, and unlike most of this list it
does not need any new physics in the model.

---

## Tasks

Tracked in `TASKS.md`, not restated here. They sat in this document as full rows
for a day without ever being transferred, which is why they are a pointer now:
a findings document that carries its own to-do list is a second copy of the
tracker, and the copy that nobody works from is the one that drifts.

| finding | id |
| --- | --- |
| 1. runoff amplifies everything by six | `BUDG-1`, `BUDG-2` |
| 2. Penman over a dry land column | `HYD-11` |
| 3. the error budget is a hand-edited generator output | `BUDG-3` |
| 4. the biosphere reaches the climate only through albedo | `BIO-4` |
| 5. knocked down; no task |  |
