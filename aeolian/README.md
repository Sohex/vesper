# Aeolian

Mineral dust: emission, transport and deposition, computed offline from a
climatology rather than inside the climate model. `notes/dust.md` carries the
decision and why it is not in the GCM; this component is the implementation.

```bash
python aeolian/scripts/build_dust.py                          # the baseline source map
python aeolian/scripts/build_dust.py --variant arid_bare_ground   # the bracket on vegetation
python aeolian/scripts/build_dust_source_fields.py            # the same map as .sra boundary fields
python aeolian/scripts/build_dust_source_fields.py --self-test  # the identity, and that it can fail
```

Products are `analysis/dust_<variant>.json` and a matching `.nc` carrying annual
mean optical depth, deposition flux, emission and the erodible fraction.

`build_dust_source_fields.py` is the bridge to the in-model scheme: it writes the
same per-cell source map as ExoPlaSim surface codes 1801, 1802 and 1803, plus the
`aero_nl` group that goes with them. It is a separate script rather than a flag
on `build_dust.py` because it produces a boundary condition rather than an
answer, and because it must be regenerated whenever the terrain, the lake
solution or the soil moves, which is a different cadence from the offline run.

## STATUS: the chain is sound, the answer is undetermined by two parameters

The 1700x emission excess is found, and it was not in this component.

**Every orbit of the source run carries a corrupted first output bin.** In
`spd`, `ua` and `va` the lower troposphere is inflated, worsening downward from
1.02x at the model top to **7.5x at the bottom level**. It is systematic rather
than a restart shock: orbits 86, 87, 88 and 90 of `run_b014469b8091` give 7.50,
7.57, 7.54 and 7.58. The snapshot climatology built from the same run is clean,
so it is the regular 12-bin output's first interval specifically.

Emission goes as roughly u* cubed above a threshold, so that one bin was
**100.00% of the annual total**. Excluding it, and with the transport fixed
below, emission falls from 3.4e6 Tg per Earth year to a range that brackets
Earth's ~2000. Nothing in the physics or the unit conversions was wrong; the
driver was.

`flag_anomalous_bins` now detects it and the run excludes it loudly rather than
consuming it. That is a guard against a known defect, not defensive habit.

**Transport now converges, and the reason it did not is worth recording.** It
was neither a CFL violation nor a cycle nor a missing sink. The explicit step is
CFL-limited by the convergence of the longitude grid at the poles, which gave a
14 km cell and a 597 s step, so 4000 iterations covered 27.6 days -- against a
relaxation time of 1351 days for the 0.5 um bin where nothing rains. The
integration was being stopped at 2% of the way. A declared polar `cos(lat)`
floor of 0.2 relaxes the step to 6029 s and the iteration cap is now 40000; all
cases report converged.

### What it says now, with the wind tail measured and gravity in the threshold

Measured 2026-08-17 on the baseline climatology, at k = 2.012 fitted from 1,463
three-hourly samples (DUST-5) with the saltation threshold corrected for this world's
gravity (DUST-6), on a climatology built from ten clean NLOWIO = 0 orbits, and
re-run the same day once HYD-13 corrected the mapping the lake solution is built
through:

| | z0 = 3e-6 m | z0 = 1e-4 m | z0 = 1e-3 m |
| --- | ---: | ---: | ---: |
| emission, Tg per Earth year | 171158 | 14029 | 2.4 |
| land-mean optical depth | 5.235 | 0.376 | 0.00002 |
| deposition, g/m2 per Earth year | 370.5 | 27.4 | 0.002 |

Earth for scale: about 2000 Tg per year and a land-mean dust optical depth near
0.03. The central roughness is therefore about seven Earths of emission at
twelve times Earth's optical depth. This is a dusty world.

**Both of the large corrections went in opposite directions and neither cancelled
the other.** Measuring the wind tail raised emission 8.8x, because a Weibull
fitted to 32 snapshots 5.7 days apart is too narrow and biases the shape high.
Putting this world's gravity into the saltation threshold then cut it 30%, a
6.9% change in the threshold amplified by the same u* nonlinearity working the
other way. The threshold correction is the fourth root of the gravity ratio and
not the square root; `aeolian/config/dust.yaml` carries the derivation.

The roughness bracket is the only large uncertainty left and it no longer spans
the answer: the rough end is 2.4 Tg per year rather than the zero it once read.

### The reopening test, stated without tuning

`notes/dust.md` reopens the in-model question at a land-mean optical depth above
0.10. **It crosses by 3.8x at the central roughness** and by 52x at the smooth
end. Only the roughest end of the bracket is below it.

That was a pre-committed threshold, fixed before the answer was known, and what
it commits to is that a prescribed field is no longer defensible and the emission
scheme belongs inside the model. That is DUST-3's problem now, not a conclusion
this file should quietly draw.

### What would decide it, and what to ask for

The wind tail can be measured rather than fitted. ExoPlaSim writes high-cadence
output and every run directory already contains a `highcadence.nl`. A short
high-cadence segment off the settled baseline would give a real gust
distribution instead of a Weibull fitted to 32 samples.

**Specification: one orbit, sampled every four timesteps.** At the 45-minute
timestep that is 3-hourly, 1464 samples per cell over 183 days, which resolves
the diurnal cycle on a 30-hour day and gives enough independent samples to fit
the tail rather than the body. Bottom-level wind alone is sufficient. That is
roughly 190 MB per variable at T42, and about 25 minutes of model time at the
rate the baseline is running.

The baseline has settled, so this is runnable now. It is DUST-5, and it is the
one measurement that would turn the reopening test from crossed-with-unknown-
magnitude into a number.

## What the component does get right

**The source map, which is the thing the built-in GCM scheme gets wrong.**
`fcoeff * land_mask` emits as much from forest as from salt pan. This one starts
from `substrate_class`, which is consolidated lithology plus closed-basin fill,
so the only unconsolidated material is the fill; weights the two barren classes
separately, because a cemented salt crust is not a silicate soil and Kok's
fragmentation theory does not describe halite cement; removes standing water
from the solved lake extent; and removes snow. It gives **16.1% of land** as
bare erodible ground against a playa fraction of 23.9%, and the difference is
lakes and crust rather than an assumption.

**The physics is grounded rather than recited.** Every constant traces to a
fetched primary source, listed in `config/dust.yaml` with what it is worth:
Kok et al. (2014) for emission, Kok (2011) for the emitted size distribution,
Fecan et al. (1999) for the soil-moisture threshold, Marticorena and Bergametti
(1995) for the drag partition, Sportisse (2007) for below-cloud scavenging.

**Settling uses this world's gravity.** At 12.81 m/s2 a given particle falls
1.31x faster than terrestrial intuition, so dust lifetime here is shorter than
Earth analogues suggest. That is in `settling_velocity` rather than in a comment.

## What is declared rather than derived, and where the brackets are

`config/dust.yaml` carries all of it. The widest are, in order:

1. **The subgrid wind shape**, now measured but only to an upper bound.
2. **The aeolian roughness of the erodible surface**, bracketed 3e-6 to 1e-3 m,
   worth a factor of 40 in emission across that range. The grid-cell roughness
   field is deliberately NOT used: its median over source cells is 0.49 m, and
   feeding a 15 km orographic variance to a scheme built for centimetre-scale
   roughness elements returns zero emission everywhere. Sheltering of a patch by
   the terrain around it is therefore not represented, which biases emission up.
3. **The evaporite erodible weight**, 0.1 with a bracket of 0.0 to 0.3, for
   crust cementation. Declared suppression, not measured efficiency.
4. **Vegetation cover**, which is not modelled at all. Non-barren land is
   assumed to carry a canopy and not emit, following the project's existing
   declared position, which is generous on a world whose median land runoff is a
   few mm per Earth year. `--variant arid_bare_ground` brackets it.
5. **Which refractive indices this world's dust has**, `optics.indices`, chosen
   as the measured datasets in the shortwave and OPAC in the thermal infrared,
   with OPAC kept as the absorbing bracket. This is the one entry in the file
   that something OUTSIDE this component also reads: the aerofile ExoPlaSim's
   radiation loads and the offline forcing are both built from it, through
   `exoplasim/scripts/dust_indices.py`, so that the model and the pricing
   describe one particle. The config carries the argument; DUST-12.

## What it does not do

No dust-climate feedback: the climatology is an input and does not respond.
No vertical structure: a well-mixed column of declared scale height advected by
a single steering wind. No inter-bin microphysics, which is correct for mineral
dust because it neither coagulates nor grows appreciably.

## One-off tools

- `scripts/dust_runoff_sensitivity.py` -- one-off: converts a precipitation change into a basin count, which is what prices the catchment half of the dust question in WORKFLOW A4. Registered under `one_offs` in `config/pipeline.yaml`.
- `scripts/extract_high_cadence_wind.py` -- one-off: pulls instantaneous near-surface winds out of a high-cadence run segment. This is what DUST-5 measured and what `build_dust.py --gust-samples` must be given; the 32-sample snapshot climatology fits a Weibull shape about twice too high and costs a factor of 40 on emission. Registered under `one_offs` in `config/pipeline.yaml`.

## The in-model port, and why this component survives it

`aeolian/notes/in-model-dust.md` is the design for putting emission, deposition
and scavenging inside ExoPlaSim (DUST-3), with each piece's predicted effect and
what result would falsify it. Three of its patches are written and verified but
not yet applied, and they are listed in `PENDING_PATCHES` in
`exoplasim/scripts/rebuild_binaries.py`. They stack in this order:

- `exoplasim/patches/exoplasim-3.4.2-aerocore-defects.patch`, seven latent
  defects in `aerocore.f90` and `aeromod.f90`, unconditional.
- `exoplasim/patches/exoplasim-3.4.2-aerosol-deposition.patch`, a dry deposition
  velocity behind `ldepvel` and Sportisse below-cloud scavenging behind
  `lwetdep`, both defaulting to off and both reading their coefficients from
  `aeolian/config/dust.yaml` rather than carrying Fortran defaults.
- `exoplasim/patches/exoplasim-3.4.2-dust-emission.patch`, the three boundary
  fields, the gathers, the `aero_nl` calibration group and Kok (2014) equation
  18 in the source term, behind `ldustemit` and defaulting to off.

The run side is `model.dust_emission` in `config/planet.yaml`, which is absent
and therefore `none` by default. Setting it makes `run_exoplasim.py` stage the
three fields, write the whole `aero_nl` group from the provenance file beside
them, and set `L_AERO = 1` and `l_source = 2`. `model.dust_dry_deposition` and
`model.dust_wet_scavenging` switch on the two removal terms independently, which
is what lets each be its own A/B arm off one binary.

The emitted dust is RADIATIVELY INERT on that path and the driver sets
`l_aerorad = 0` to say so. That is not a preference: `radmod`'s own `apart` is
never populated from the namelist, which is upstream defect 1 and still open
because it lives in `radmod.f90`, and the longwave aerosol term does not exist
yet. Turning the radiation on before those land would price this world's dust at
a small fraction of its true optical depth and cool with it without warming.

**This component does not retire when that lands.** ExoPlaSim's aerosol is one
tracer with one radius and one density fixed at compile time, so the in-model
chain gets the burden and the optical depth and cannot get the size-resolved
DEPOSITION field that pedology and the phosphorus budget read. That is the
declared cost of DUST-8 and the note carries the numbers.

## Sea salt

The second aerosol, added 2026-08-19 for CLIM-27. `notes/audits/unpriced-terms.md`
finding 2 says why: this project priced mineral dust carefully and had never
asked what any other aerosol was worth, on a world whose ocean is more than half
the surface.

```bash
python aeolian/scripts/sea_salt_source.py     # the Earth check, and nothing else
python aeolian/scripts/sea_salt_optics.py     # step sea_salt_optics
python aeolian/scripts/build_sea_salt.py      # step sea_salt
python aeolian/scripts/build_sea_salt.py --bracket   # every declared end
```

`sea_salt_source.py` is the source function, Grythe et al. (2014) equation 7,
shared by the other two because one needs its mass and the other its shape.
`sea_salt_optics.py` writes `analysis/sea_salt_optics.json`, band-averaged
optics per dry size bin and per relative humidity. `build_sea_salt.py` writes
`aeolian/analysis/sea_salt_baseline.{nc,json}`: emission, burden, optical depth
and the top-of-atmosphere shortwave forcing, at both ends of every declared
bracket. Constants are in `aeolian/config/sea_salt.yaml`.

**Three things about it are not obvious from the dust component beside it.**

The particle is WET and the mass budget is DRY. A sea-salt particle at ambient
humidity is a solution droplet, twice the dry diameter at 80% relative humidity
and a quarter of the dry mass fraction, so its extinction is computed on the wet
particle and divided by the dry mass. The growth curve is OPAC's, in
`exoplasim/data/sea_salt/`, and the optics script checks that OPAC's growth
factor and its wet density say the same thing.

Wet removal is NOT dust's. Sea salt is the most cloud-condensation-active
aerosol there is and is removed in cloud rather than below it; using dust's
below-cloud coefficient gives the accumulation mode a three-week lifetime
against a measured day. The level is anchored to Jaegle et al. (2011).

The spume mode is reported apart. Grythe's third lognormal is centred at 30 um
and its lower tail dominates the emitted mass below the 10 um cut, so the
product carries an `all_modes` and a `no_spume` variant and any comparison with
an Earth number is against the second. Neither is a correction to the other.
