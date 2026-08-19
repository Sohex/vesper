# The prescribed-dust run: specification, gates and a prediction

DUST-11. Written 2026-08-18, before the run exists, so the prediction below is a
prediction and not a description. `notes/dust.md` carries why the catchment limb
matters and what the lake limb already measured; this file says what to run, what
would count as measuring it, and what would mean the patch is wrong rather than
the world being surprising.

## What the question is

The carve criterion is an aridity index `(E_lake - P) / R` compared against a
catchment's own geometry. `R` is catchment runoff, and runoff is a small residual
of two large numbers, so a change in precipitation arrives at the verdict
amplified. Radiatively dark mineral dust lowers the simulation's rainfall
through a fast adjustment to the modelled atmosphere's energy budget, which is
a response of the circulation. No
surface energy balance can reach it, which is why DUST-10 could settle the lake
limb offline and this one needs a climate run.

## What was built

| piece | where |
| --- | --- |
| the model change | `exoplasim/patches/exoplasim-3.4.2-prescribed-dust.patch`, resident and compiled into all five binaries since 2026-08-18 |
| the boundary field | `exoplasim/scripts/build_surface_dust.py` -> surface code 1811 |
| the run wiring | `model.dust_source` in `run_exoplasim.py`, absent by default |
| the prediction | `aeolian/scripts/dust_runoff_sensitivity.py` -> `aeolian/analysis/dust_runoff_sensitivity.json` |

The patch adds a prescribed column optical depth to the shortwave, a grey
absorber to the longwave, and the band-1 shortwave absorption term upstream
computed and never used. It does not touch `aerocore`: nothing is emitted,
transported or removed, so the bottom-level sink, the settling term and the
un-populated `apart` are all off the path.

**Measured on 2026-08-18**, from `aeolian/analysis/dust_baseline.nc` after the
chain was re-run on HYD-13's lake solution. Every row of it was superseded later
the same day when the chain was regenerated with the snapshot-fitted wind tail
instead of DUST-5's measured one and the burden fell by a factor of 70; see the
last section of `notes/dust.md`. Rebuild the field before running anything
against it:

| quantity | value |
| --- | ---: |
| land-mean band-1 column optical depth | 0.38028 |
| global-mean band-1 column optical depth | 0.24698 |
| maximum | 5.1019 |
| `DUSTQLW`, thermal-IR absorption per unit band-1 extinction | 0.113334 |
| Planck-weighted thermal-IR mass absorption, at 289.71 K | 88.33 m2/kg |
| `DUSTHSC`, dust scale height | 3000 m |

## The run

Two runs, not one, and the second is not optional.

**Control.** Branch from the baseline's final restart on the REBUILT binary with
`model.dust_source` absent. Same length, same segments, same output settings as
the dust run.

**Dust.** The same branch with `model.dust_source: prescribed`.

Everything else is held: flux ratio, CO2, surface fields, resolution, ranks and
the radiation configuration. One thing changes between the two runs and it is the
dust.

### Which radiation configuration this assumes, and why it has to be said

**The pair runs on whatever radiation the baseline ran on, and both members of
the pair run on the same one.** As of 2026-08-18 that is a 4965 K BLACKBODY, not
the k25v spectrum: `radmod.f90` takes the spectrum branch only when `NSTARFILE >
0`, no run's namelist carries it, and `solarini` has been building a blackbody
throughout. The spectrum fix is in flight, and a second defect in the same
subroutine means declaring the spectrum without patching the Rayleigh reference
would weaken Rayleigh scattering by 3.4x.

So there are now two things moving in the radiation, and WORKFLOW.md section 6
says to change one at a time. Concretely:

- If the spectrum lands BEFORE this pair, run both members after it. Fine.
- If it lands BETWEEN the control and the dust run, **the pair is void** and both
  members have to be re-run. Four hours, not a judgment call.
- If it lands after, the pair measures dust under a blackbody. That is still a
  valid measurement of the dust term, and it should be labelled as one.

The spectrum matters to this measurement in a specific way rather than a general
one: the band-1 flux share sets how much of the star's output meets the more
absorbing of the two dust bands, so the atmospheric absorption in the prediction
below moves with it. It does not change the sign of anything here.

    # once, before either run. The patch is already resident, so this is a check
    # rather than an application; --verify fails if a .venv reinstall dropped it.
    python exoplasim/scripts/build_surface_dust.py
    python exoplasim/scripts/rebuild_binaries.py --verify

    # then, per run
    python exoplasim/scripts/run_exoplasim.py --restart-from <baseline>/MOST_REST.<n>
    python exoplasim/scripts/build_climatology.py <run>

`model.dust_source` is deliberately NOT in `config/planet.yaml`. Adding a key
moves `config_sha256` and `continue_exoplasim.py` will not resume a run whose
config hash has changed, so the key goes in when the dust run is prepared and
comes out afterwards, exactly as `soil_water_source` was held back for the
bootstrap. `model.dust_scale` is the optional multiplier for the shelter bracket
and defaults to 1.

**Length: 30 orbits of adjustment, then 10 with `NLOWIO = 0` and seasonal output
for the climatology.** The baseline's own fitted relaxation is about 9.6 orbits,
so 30 is three e-foldings; the climatology block matches the baseline's, which is
also 10 orbits, so the two are averaged over the same window length. At the
baseline's measured rate -- 1.8 min per orbit on 16 ranks in spin-up, 6.4 with
`NLOWIO = 0` and snapshots -- that is about two hours per run and four for the
pair.

**Why a control run rather than the existing baseline climatology.** The baseline
is quasi-equilibrated with a mean TOA balance of -0.60 W/m2, so it is still
drifting slowly, and the dust run would be averaged over a later window than the
baseline climatology was. Differencing those two would measure dust plus drift
plus a rebuilt binary, and nothing would separate them. A control on the same
binary over the same orbits removes all three at the cost of two hours.

**Seeding from a restart is allowed here and the guard was taught why.**
`run_exoplasim.py` refuses `--restart-from` when the set of surface codes
changes, because landmod reads soil water, roughness and albedo from the restart
rather than from the `.sra` files. That does not apply to code 1811: `radini`
reads it with `mpsurfgp` outside any `nrestart` test, where landmod's block is
inside `if (nrestart == 0)`. The dust code is exempted by name in the guard, and
anything else added there has to be checked the same way, in the Fortran.

## Gates, before anything is measured

Each of these has a right answer and can fail.

1. **The null.** One orbit on the rebuilt binary with `dust_source` absent,
   against one orbit continued on the old binary from the same restart. The patch
   is gated on `ndustrad`, so the first output bin must agree to roundoff. If it
   does not, the patch is not the no-op it claims to be and nothing downstream
   compares anything.
2. **The column identity.** `dustprof` normalises its vertical weights, so the
   sum of the per-layer optical depths must equal the prescribed field exactly.
   The model reports `max |column - prescribed|` once, on the first radiation
   step. Anything above about 1e-5 means the vertical distribution is losing mass
   and every optical depth below it is wrong by that much.
3. **The band ratio.** `build_surface_dust.py` refuses to write unless the
   aerofile's `qex2/qex1` equals the optics' `MEE2/MEE1`. `apart` and `rhop`
   cancel out of that ratio, so it is an identity between two files that are
   generated separately and it is the only thing tying the boundary field to the
   optics the model reads.
4. **The surface shortwave.** `rss` must FALL, everywhere, including over the
   bright closed-basin fill where the top-of-atmosphere term is positive. That is
   the whole content of DUST-10's correction, and a positive surface shortwave
   perturbation over playa would mean the model is reporting a top-of-atmosphere
   quantity or the aerosol is not reaching the surface budget.
5. **The surface longwave.** `rls` must move UP: the dust emits downward. A
   negative change means the sign of the longwave term is wrong, which is the
   single most likely way for this patch to be broken, because it is the only
   piece that goes inside the longwave solver.
6. **The global mean temperature.** Between about -1.5 and +1.0 K. Shortwave-only
   dust is -3.8 to -3.3 W/m2 in the global mean against a true +0.55 to +0.67
   (`analysis/dust_forcing.json`, read on 2026-08-18, which is the authority for
   those and has moved once already), an error of 4.0 to 4.3 W/m2. On the
   canonical conversion in `lib/sensitivity.py`, 0.832 K per W/m2 absorbed, that
   is 3.3 to 3.6 K of spurious cooling, which is second in the whole error
   budget. **A cooling of more than about 2 K is the signature of the longwave
   term not working**, not of a surprising world. Check `DUSTQLW` reached the
   namelist and that `ztaudu` is multiplying into `ztau`.

## What counts as measuring the catchment term

On the two climatologies, differenced:

- land-mean precipitation, evaporation and `P - E`;
- catchment-mean runoff, from `carve_verdict.py`'s own aggregation, which is the
  quantity the criterion divides by and not `mrro`;
- the carve list, taken on the dust climatology with

      python hydrography/scripts/carve_verdict.py \
          --climatology <dust climatology> \
          --output hydrography/analysis/carve_verdict_dust.json

  and compared against the same command on the control climatology.

**Do not pass `--dust-forcing`.** That flag adds DUST-10's offline surface
perturbation to a climatology that does not contain dust. This one does. Adding
both double-counts the lake limb, and the two numbers must not be summed for the
same reason: the run contains both limbs at once, so DUST-10's +23 basins becomes
a cross-check on the run rather than a term to add to it.

## The prediction, made before the run

Recorded here so it can be wrong.

**Direction: the dust-free carve list OVER-carves.** Dust suppresses
precipitation, runoff falls, basins that were judged to overflow do not, and a
list taken without dust cuts outlets that should not be cut. This is the
irreversible direction, which is why it gates.

**Size.** The chain, with the number each step rests on:

| step | value | from |
| --- | ---: | --- |
| global-mean band-1 optical depth | 0.24698 | the field |
| flux-weighted absorption optical depth | 0.0069 | the aerofile's own single-scattering albedos, at a band-1 share near 0.38 |
| shortwave absorbed in the layer, global mean | 3 to 5 W/m2 | that, at an insolation-weighted slant path of 2 |
| plus longwave heating of the column | 1 to 3 W/m2 | smaller than the offline estimate, because the overlap is now explicit |
| against global latent heating | 77 W/m2 | the baseline |
| fast adjustment on precipitation | -5% to -10% | the ratio |
| minus the slow response to +0.2 to +0.5 K | +0.5% to +1.5% | the global-mean forcing |
| **net global precipitation** | **-4% to -9%** | |

Runoff amplifies that, because it is the residual. Measured on the baseline
climatology by `dust_runoff_sensitivity.py`, 2026-08-18:

| dP | evaporation follows | catchment runoff | basins carved, from 1,721 |
| ---: | ---: | ---: | ---: |
| -5% | 1.00 | -5.0% | -24 |
| -5% | 0.80 | -8.3% | -121 |
| -5% | 0.50 | -12.6% | -195 |
| -10% | 1.00 | -10.0% | -51 |
| -10% | 0.80 | -16.2% | -201 |
| -10% | 0.50 | -24.1% | -309 |
| -15% | 1.00 | -15.0% | -89 |
| -15% | 0.80 | -23.9% | -288 |
| -15% | 0.50 | -34.8% | -415 |

Land evaporation on this world is mostly supply-limited, which argues that it
follows precipitation closely; runoff is the overflow of a bucket sitting at a
median 15% of capacity, which is a threshold quantity and argues that it falls
faster than precipitation does. Those pull opposite ways and the second is the
stronger, so the expected following fraction is between 0.7 and 1.0.

**So: catchment runoff falls 6% to 20%, and the carve list shrinks by 50 to 250
basins**, against 1,721. The lake limb runs the other way at about 23 basins and
does not change the sign.

**What would falsify it.** The carve list growing; land precipitation moving by
less than 1% or more than 15%; or runoff changing by less than half or more than
three times the precipitation change. Any of those, with gates 1 to 6 passed, is
a result about this world. Any of them with a gate failed is a result about the
patch.

## The equivalence this replaces

`notes/dust.md` translated a runoff change into basins by asserting that -10% on
runoff is worth about what +10% on lake evaporation is worth, which is 144
basins. That was an assumption about a ratio and it does not survive measurement.
Scaling the criterion's denominator alone, with precipitation and evaporation
left where they are:

| runoff scaled by | basins carved, from 1,721 |
| ---: | ---: |
| -40% | -29 |
| -20% | -9 |
| -10% | -7 |
| +10% | +4 |

Twenty times smaller than the shorthand. The reason is that the aridity index is
a ratio whose SIGN varies across the population: the median basin on this world
has `P > E` and therefore a negative index, and scaling the denominator makes a
positive index larger and a negative one more negative, moving the two halves
toward opposite verdicts. Reducing precipitation raises the numerator for every
basin as well as cutting the denominator, so it does not cancel, and that is why
the two tables above differ by a factor of twenty at the same runoff change.

The shorthand's magnitude happened to land near the right answer; its mechanism
did not. Quote the measured table.

## What this run cannot do

- **No dust-climate feedback.** The field is a fixed annual mean, so the dust
  cannot respond to the winds it changes. That is DUST-2's objection and it is
  the reason DUST-3 exists; it is accepted here because the question is what a
  given burden does, one iteration deep.
- **No deposition.** Nothing settles, so loess and the phosphorus return leg
  still come from the offline chain. DUST-8.
- **One effective radius**, 0.98 um, optically right and mass-weighted wrong.
- **The shelter bracket is not spanned.** The central roughness is used. The low
  end is a land-mean optical depth of 5.2, which is a different world rather than
  an error bar on this one, and `model.dust_scale` exists for anyone who wants to
  test that claim cheaply.
- **One effective set of refractive indices, and it is the config's.** DUST-12,
  2026-08-18. The offline forcing used to compute everything from OPAC while the
  aerofile the model reads was built from the measured datasets, an absorption
  optical depth a factor of 2.5 apart for the same burden; both sides now resolve
  `aeolian/config/dust.yaml`'s `optics.indices` through
  `exoplasim/scripts/dust_indices.py`. The aerofile did not move, so the patch,
  the boundary field and the binaries are unaffected -- what moved is
  `analysis/dust_forcing.json` and `analysis/dust_surface_forcing.nc`, which are
  now the same dust the run will contain and can be compared against it. Gate 6's
  window still stands: it was set on the shortwave-only error, which is carried
  by the thermal term and is unchanged.
- **The thermal infrared is still OPAC's** on both sides, because no measured
  dataset in this repository spans 4-40 um. `aeolian/config/dust.yaml` says why
  and what the alternatives were.
