# Putting dust emission in the model: what it actually takes

Design for DUST-3, written 2026-08-17 from source inspection, revised 2026-08-18
when items 1 to 3 and then item 5 were written as patches. The decision to go
in-model is in
`notes/dust.md` and rests on DUST-2: the reopening test crosses by 3.8x, and a
prescribed field cannot respond to the winds the dust itself changes.

**The headline is that "replace one line" was wrong, and by a lot.** The hook at
`aerocore.f90:633` is one line, but around it sit a serial-on-NROOT transport
call, a resolution-dependent sink, a shortwave-only radiation scheme, no wet
deposition, one hard-wired particle size, and SEVEN upstream numerical defects.
The line is the smallest part.

## What is already done

`exoplasim/data/dust/vesper_dust_aerosol.dat`, from
`exoplasim/scripts/dust_aerofile.py`. Validated: `radmod` reconstructs the
single-scattering albedo and backscatter ratio of `analysis/dust_optics.json`,
and the optical depth round-trips against the offline chain's own value exactly.
The Mie plausibility check is on the UNSCALED Qext at the optical radius, 2.65
and 2.61 for a micron particle; the file itself carries those times the
burden-match rescale, which is not an efficiency. The script prints both.

Five aerosol-path patches, all RESIDENT in the `vendor/exoplasim` subtree
(DUST-3, closed 2026-08-18). The files under `exoplasim/patches/` are the
record of what each changed and why; applied state lives in the subtree's
history now -- the `PENDING_PATCHES` registry `rebuild_binaries.py` once
policed no longer exists. Each patch names its base shas, which is where the
stacking order is recorded:

| patch | items | touches |
| --- | --- | --- |
| `exoplasim-3.4.2-aerocore-defects.patch` | item 1, defects 2 to 7 | `aerocore.f90`, `aeromod.f90` |
| `exoplasim-3.4.2-aerosol-deposition.patch` | items 2 and 3 | `aerocore.f90`, `aeromod.f90` |
| `exoplasim-3.4.2-dust-emission.patch` | item 4 | `aerocore.f90`, `aeromod.f90`, `surfmod.f90`, `plasim.f90`, `make_plasim` |
| `exoplasim-3.4.2-aerosol-apart.patch` | item 1, defect 1 | `radmod.f90`, `aeromod.f90` |
| `exoplasim-3.4.2-aerosol-longwave.patch` | item 5 | `radmod.f90` |
`aeolian/scripts/build_dust_source_fields.py` is the generator for item 4's three
boundary fields, and `exoplasim/scripts/run_exoplasim.py` gained
`enable_dust_emission`, which writes the whole `aero_nl` group from the
provenance file beside them.

## Seven upstream defects, and none of them was reachable

Six of the seven are in `aerocore.f90` and `aeromod.f90`; defect 1 is in
`radmod.f90`. All are LATENT: every run this project has made sets `L_AERO = 0`
in `plasim_namelist`, `aero_main` runs only when
`nsela == 1 .and. nkits == 0 .and. l_aero > 0`, and so `aerocore` has never
executed. Defect 1's path needs `l_aerorad == 1` on top of that and it defaults
to 0. That is why they survived, and it is why fixing them needs no namelist
switch: on every configuration this project runs, a binary with the defect
patches must produce output BIT-IDENTICAL to one without them. That is the test
to run first, and it can fail.

1. **`radmod`'s `apart` is never populated from the namelist.** `apart` is
   declared twice, in `aeromod` and in `radmod`, and `aero_ini` use-associates
   only `l_aerorad` and `aerofile`, so the transport settled the particle the
   run asked for while the radiation kept its 50e-9 haze default. Optical depth
   goes as `apart**2` at fixed number density, so it came out 1/385 of intent at
   the optical effective radius and 1/1948 at the burden-matched one. **Fixed**
   in `exoplasim-3.4.2-aerosol-apart.patch`: `aero_ini` copies the value across
   under a rename and `radini` broadcasts it. It cannot go the other way round
   because `radmod` cannot `use aeromod` -- `aeromod` already uses `radmod` and
   `make_plasim` compiles it second. This is the gate DUST-8's retune was
   waiting on, and the retune is now applied.
2. `aerocore.f90:1156` and `:1183`, `mmr2n` and `n2mmr`: `(4/3)` is integer
   division and evaluates to 1, so the sphere volume is 4/3 too small and the
   number density exactly 4/3 = 33.3% high. **Fixed.**
3. `aerocore.f90:949`, `viscos`: `(4/25)` is integer division and evaluates to 0,
   so the collision integral's temperature dependence vanishes. This has a right
   answer rather than a plausible one, which is what makes it a test: with the
   exponent restored, `mu` for N2 at 288 K comes out 1.77e-5 Pa s against the
   measured 1.76e-5; without it, 1.48e-5. Stokes settling goes as 1/mu, so it was
   13 to 21 per cent fast over 200-320 K in N2. The design note previously said
   "a fifth to a quarter", which is the figure for `l_bulk = 2`, an H2
   atmosphere; this world is N2. **Fixed.**
4. **The settling flux never carried the timestep.** `gsettle` returns kg/m2/s;
   the FFSL fluxes fx, fy and fz carry dt inside their Courant numbers and gz did
   not, so one model step applied one SECOND of sedimentation. At a 2700 s step
   that is the whole factor. The corrected term is stable at this step -- the
   settling Courant number in the bottom layer is 0.004 for a 2 um particle and
   0.36 for a 20 um one -- which matters because "it would have gone unstable" is
   the only defence the omission could have had. **Fixed.**
5. **The bottom level was updated twice.** `do k=2,NL` already covered k=NL, and
   the "Now for k=nl" block below it applied fx, fy and fz to that level a second
   time. The two settling terms cancelled exactly, -gz(nl) from the loop and
   +gz(nl) from the block, so nothing ever left the atmosphere: dust settling out
   of the bottom layer was handed straight back to it. **That is what the
   99%-per-step sink exists to hide**, and it is why the sink cannot be removed
   until this is fixed. **Fixed**: the loop stops at NL-1 and the bottom block
   takes gz(nl-1) in and gz(nl) out, which is dry deposition by sedimentation.
6. **The polar caps never settled.** Both cap rows carried the fy and fz terms
   and no gz term at all, at every level. **Fixed.**
7. **The gathered boundary fields were in the wrong hemisphere.** `aero_main`
   gathers with `mpgagp`, which returns the MODEL's latitude order, and hands the
   result to `aerocore`, whose tracer array `plasim.f90` flips south-to-north on
   the way in. So `l_source = 2`'s land mask and `l_source = 1`'s solar zenith
   angle drove the source in the mirror hemisphere. This world's land is not
   hemispherically symmetric, so the test is direct: one step with `l_source = 2`
   and a nonzero `fcoeff`, and the bottom-level mixing ratio must correlate with
   the land mask rather than with its reflection. **Fixed.**

Defects 4 to 7 were found while fixing 2 and 3 and are not in the original
report. Expect more once it runs: this is code that was presumably believed to
work and has never been exercised.

### Predicted effects, stated before they run

Per `docs/src/pipeline/sequencing.md` A3, with what result would mean "wrong".

| change | prediction | falsified by |
| --- | --- | --- |
| all seven, on any configuration that exists today | output bit-identical, because `aerocore` is unreachable at `L_AERO = 0` and the radiation's aerosol block at `l_aerorad = 0` | any difference at all, which would mean the aerosol path is reachable and `notes/dust.md` is wrong about that |
| defect 1 | with the aerosol on, band-1 layer optical depth rises by exactly `(apart/50e-9)**2`: 384.6x at 0.98054 um, 1948.0x at the burden-matched 2.20682 um | any other ratio, and specifically 1.0, which means the copy or the broadcast is not reaching `swr` |
| defect 2 | `nrho` falls by exactly 0.75 everywhere, and shortwave aerosol optical depth with it | a ratio that is not 0.75 to round-off |
| defect 3 | `mu` at 288 K in N2 rises from 1.48e-5 to 1.77e-5 Pa s, within 1% of the measured 1.76e-5; terminal velocity falls 13 to 21% over 200-320 K | a ratio outside (T/95.5)**(-0.16) cell by cell |
| defect 4 | sedimentation mass flux rises by exactly `deltsec` | anything other than the timestep in seconds |
| defect 5 | with the bottom sink disabled and no wet term, total aerosol mass must DECREASE at the surface sedimentation flux; before the fix it cannot decrease at all | a non-decreasing total after the fix, or a decrease that does not match the integral of rhog*mmr*v_s over the surface |
| defect 7 | bottom-level mmr correlates with the land mask at 1.0 | correlation with the reflected mask instead |

## The sink at `aerocore.f90:869` was a stability hack, not a parameterisation

It is not undocumented -- there is a comment one line above, "put in a sink term
at the bottom level to avoid infinite build-up of haze particles" -- and that is
exactly what it is.

```fortran
mmr(:,:,nl,ic) = mmr(:,:,nl,ic)*10e-3
```

`10e-3` is 0.01, so it removes **99% of the bottom layer every timestep**, after
transport, inside the tracer loop. `aerocore` is called once per timestep from
`gridpointd`. Measured against this world's own bottom layer -- sigma 0.9524 to
1.0, 4647 Pa, 309 m thick at 2700 s:

| | |
| --- | --- |
| e-folding time | 9.8 minutes |
| implied deposition velocity | 52.7 cm/s |
| real dust dry deposition | 0.1 to 2 cm/s |

So it is 26 to 527 times too fast. **And it does not scale with the timestep**,
so the implied velocity is inversely proportional to it: 158 cm/s at a 15-minute
step, 40 at an hour. A physical process whose rate depends on the integration
step is not a physical process.

**Why it was necessary is defect 5, not the source term.** The earlier reading
here was that the existing dust source *sets* the bottom-level mixing ratio to
`fcoeff*land` every step rather than adding a flux, so source and sink degenerate
into "surface concentration is pinned". That is true but it is not the reason.
The reason is that the sedimentation flux out of the bottom layer was
algebraically cancelled, so aerosol reaching the ground was returned to the air
and the build-up the comment describes was unbounded by construction.

**Replaced, behind `ldepvel`.** `ldepvel = 0` keeps the legacy line and is the
default, so one executable runs both arms of the A/B. `ldepvel = 1` gives the
bottom layer `exp(-vdaero*dt/dz)` per step, with `dz` from that layer's own
pressure thickness and gas density.

`vdaero` is the NON-gravitational part of dry deposition, turbulent transfer and
impaction, and that is deliberate: sedimentation to the surface is now the
`gz(nl)` term, so a Stokes velocity here would deposit the same mass twice. It
has no Fortran default -- `aero_ini` aborts on `ldepvel = 1` without it, the way
`radini` aborts on `ndustrad = 1` with `dustqlw` at zero -- because the value
belongs in `aeolian/config/dust.yaml` with the rest of the removal budget and a
default here would be a fourth place for it to go stale.

The right long-run form is a resistance model, `1/(r_a + r_b)`, which needs the
friction velocity and the aerodynamic roughness. Neither is gathered today, and
item 4 has to gather them for the emission scheme anyway, so that is where it
belongs.

**Prediction.** At `vdaero = 2e-3` m/s the bottom layer's non-gravitational
e-folding time is 1.79 days against the legacy 9.8 minutes, a factor of 263, and
the implied total dry deposition velocity is 0.24 cm/s -- inside the measured
0.1 to 2 cm/s, against 52.7. The test with a right answer is not the burden but
the timestep independence: halve `deltsec` and the removal per unit time must not
move, to the exponential's second order. The legacy line fails that by
construction, which is the whole argument for replacing it.

## Wet scavenging

DUST-7: `aerocore`'s argument list contained no precipitation field of any kind.
Gravitational settling was the only removal.

`aeolian/README.md` says of the offline chain that without the wet term the fine
mode's lifetime is out by an order of magnitude, so this is a first-order sink
and not a refinement. **The LMD Generic PCM has no wet deposition for its
aerosols either** (`~/git/generic_pcm`, grepped 2026-08-18: one hit, in a
chemistry header). So this is a common omission in exoplanet GCMs rather than a
peculiarity of ExoPlaSim, and there is nothing to crib.

**Added, behind `lwetdep`**, default 0. `lwetdep = 1` applies below-cloud
scavenging at `Lambda = scava * p**scavb` with `p` the total precipitation rate
in mm/h, from `dprl + dprc` in m/s. That is Sportisse (2007) and it is the same
form and the same coefficients `aeolian/config/dust.yaml` already carries for the
offline chain, so the two chains can be compared rather than merely contrasted.
`scava` and `scavb` have no Fortran defaults for the same reason `vdaero` has
none.

Applied through the whole column rather than below a diagnosed cloud base. That
is the offline chain's own approximation and is kept deliberately; it overstates
the scavenging of aerosol sitting above the precipitating layer, and it is the
first thing to refine if the wet sink turns out to dominate more than the offline
chain says it should.

**Prediction.** At the single effective radius the dry settling lifetime over a
5 km dust scale height is 80 days, while the wet lifetime is 27 days at 0.1 mm/h,
6.7 days at 1 mm/h and 1.7 days at 10 mm/h. So at fixed emission, turning
`lwetdep` on must cut the global burden by a factor of roughly 10, and the stated
bracket is 4 to 15.

- **Below 1.5x means the precipitation is not arriving.** Check the gather, check
  the latitude flip, check that `dprl` and `dprc` hold what you think at the
  point `aero_main` runs.
- **Above 50x means the units are wrong.** Reading m/s as mm/h inflates Lambda by
  `(3.6e6)**0.61` = 9985, and the signature is a burden that collapses to zero
  rather than settling to a smaller steady state. That factor is specific enough
  to identify the error from the result alone.

## Size: DUST-8, decided

**Single mode. Not per-bin arrays.** `plasimmod.f90:99` fixes
`parameter(NAERO = 1)` at compile time, and the `do jc=1,NAERO` loops are already
written -- but `apart` and `rhop` are scalars shared by every bin, `mmr2n` takes
them as scalars, the settling chain (`chamfac`, `vterm`) takes them as scalars,
and the optics are single-valued in the aerofile and in `radmod`. Promoting them
threads arrays through `aerocore`'s whole argument list and through the whole of
`radmod`'s aerosol block in both solvers. The design note recommended single mode
first and nothing since argues otherwise.

The LMD Generic PCM is a precedent and it points the same way: dust there is ONE
aerosol kind with one effective radius (`aerosol_radius.F90`), and it carries no
per-bin size or density arrays for it.

**Does 0.98 um reproduce the mass extinction efficiency?** Yes, and trivially, so
this is the wrong question to have asked. `dust_aerofile.py` sets
`Qext = 4 * apart * rho_p * MEE / 3`, and `radmod` builds optical depth as
`nrho * PI * apart**2 * Qext * dz` with `nrho` from `mmr2n`, so the geometry
cancels and `aod = column_mass * 3 * Qext / (4 * apart * rho_p)`. The MEE
round-trips to the digit for ANY `apart`, and it does so for either refractive
index dataset DUST-12 is choosing between: the two differ by 0.32% in band 1 and
0.33% in band 2 in MEE, because MEE is set by the real part while the datasets
differ in absorption. **The optics do not constrain `apart` at all.**

**What `apart` does constrain is settling, and there the single mode is wrong by
a factor that is worth writing down.** Measured 2026-08-18, from
`aeolian/config/dust.yaml`'s Kok (2011) emitted distribution and
`build_dust.settling_velocity` at this world's gravity, 288 K, 1.15 kg/m3:

| single value | settling velocity | radius that gives it |
| --- | ---: | ---: |
| optical effective radius, 0.98 um | 0.043 cm/s | 0.98 um |
| burden-matched, `1/sum(f_b/v_b)` | 0.209 cm/s | 2.21 um |
| mass-weighted, `sum(f_b v_b)` | 0.706 cm/s | 4.09 um |

So a mode at the optical radius keeps the dust up **4.8x** too long for the
burden and settles **16x** too slowly for the emitted mass.

**The Generic PCM decouples the two radii and ExoPlaSim can too.** There,
`radius(iq)` comes from the tracer definition and drives sedimentation, while
`reffrad(:,:,iaer)` is set separately in `aerosol_radius.F90` and drives only the
optics; its dust `reffrad` is 10 um, a coarse settling-relevant size rather than a
fine optical one. ExoPlaSim couples them through the single `apart`, but the
identity above means the coupling is free to break: set `apart` from the settling
requirement and multiply the aerofile's four Q values per band by
`apart_new / apart_old`, and the optical depth, the single-scattering albedo, the
backscatter ratio and the band-2 ratio are all unchanged, because every one of
them is a ratio the rescale preserves.

**So the configuration this project is now in is `apart` = 2.21 um with the Q
values scaled by 2.251**, burden-matched because the in-model chain exists for
the AOD and the precipitation response (`docs/src/pipeline/sequencing.md` A4), not for the deposition
field. APPLIED 2026-08-18, once defect 1 landed and made it bite. Three things
carried with it:

- `dust_aerofile.py` derives the radius rather than carrying it: it imports the
  aeolian component's own `emitted_mass_fractions` and `settling_velocity`, so
  the in-model mode and the offline bins settle by one formula and the number
  moves if the size distribution does.
- The plausibility check survives but is applied to the UNSCALED number: `Qext`
  at the optical radius is a believable Mie efficiency for a micron particle,
  while the scaled value is not an efficiency at all and must not be read as
  one. Both are printed and both are in the provenance sidecar.
- The invariance is CHECKED rather than asserted. The script recomputes the mass
  extinction efficiency at both radii and refuses to write the file if they
  differ, because "the algebra says the rescale is free" is exactly the kind of
  claim a sign error survives. Against the file actually written, the ratios
  move only in the seventh significant figure, which is the rounding of the
  six-decimal format on a larger mantissa and not a change in the optics.

**What the single mode costs deposition, stated rather than implied.** In steady
state the total deposition equals the total emission whatever the settling
velocity, so the cost is not in the total but in WHERE it lands. One lifetime
means one travel distance: a burden-matched mode carries the coarse 73% of the
emitted mass about 3.4x too far and the sub-micron fraction about 30x too short.
The near-source playa rim and the far-field upland both come out wrong, in
opposite directions, and those are exactly the two ends of the phosphorus
argument in `notes/dust.md`. **So the offline chain keeps producing the
deposition field for pedology and the phosphorus budget, size-resolved, and the
in-model chain is not asked for it.** Bins are a later iteration, and only if
that stops being acceptable.

## Boundary fields: three, and no fewer

The mechanism is `surfcode(code,'name')` to register and
`mpsurfgp('name',array,NHOR,1)` to read. `code_surf_file` builds the filename as
`N%03d_surf_%04d.sra` from NLAT and the code, so at T42 that is
`N064_surf_1801.sra`. Registered codes run to 1741 and 1811, so **1801-1803 is
clear**.

**The read cannot live in `aero_ini`, and that is not obvious from anywhere.**
`plasim.f90` calls `aero_ini` inside an `if (mypid == NROOT)` block, while
`mpsurfgp` is collective -- it broadcasts the read flag and scatters the field --
so a rank that never enters it hangs. The three reads therefore go in a new
`aero_surf`, called from `plasim.f90` by every rank after the `l_aero` broadcast.
They are gathered ONCE there rather than per step, because `aerocore` runs serial
on NROOT over the global grid and these fields never change.

The split is decided by where each quantity enters the flux, not by tidiness. Two
of the offline chain's terms are linear prefactors and can be multiplied together
outside the model; one enters inside the nonlinearity and cannot.

| code | name | quantity | why separate |
| ---: | --- | --- | --- |
| 1801 | `dsrcw` | erodible fraction x clipped clay fraction | both are linear prefactors on the flux, so their product is sufficient |
| 1802 | `ddrage` | drag efficiency from the aeolian roughness | scales `u*` itself, so it is inside the threshold and the exponent |
| 1803 | `dwpr` | Fecan residual moisture `w'`, percent | a pure function of clay; the gate also needs the model's own soil water, which is already in memory |

All three are dimensionless-or-percent pure functions of terrain, lithology, the
lake solution and the soil. Everything time-varying -- friction velocity, soil
moisture, snow -- comes from the model. That is the boundary this project already
draws for albedo, roughness and soil-water capacity.

A composite `E_eff` folding in the roughness would be wrong: the drag partition
multiplies `u*` before it is cubed and compared against a threshold, so it cannot
be commuted outside.

**Anything gathered here has to be flipped in latitude**, per defect 7. That is
not optional and it is not visible from the code around it.

**Written, and the generator with it.** `aeolian/scripts/build_dust_source_fields.py`
produces all three from the same `source_fractions` map `build_dust.py` uses, so
the in-model and offline chains cannot describe different ground. Its
`--self-test` is the check that the split above is right: it reconstructs the
offline flux from the three fields and the namelist alone, using a Python
transcription of the Fortran, and requires the two to agree to 1e-10. Measured
2026-08-18 at 4.9e-16, which is round-off. It then re-runs the comparison with
the drag partition removed, the gravity scaling removed, the Weibull collapsed
and the Fecan residual zeroed, and requires every one of them to BREAK the
identity, because a check that cannot be made to fail is not evidence.

## The emission law

Port `build_dust.py:emission_over_weibull`, which is Kok (2014) equation 18:

    F = Cd * f_clay * rho_a * (u*^2 - u*t^2) / u*st * (u* / u*t)^alpha

**Not the White/Kawamura saltation flux.** That form -- `(rho/g) u*^3 (1 -
u*t^2/u*^2)(1 + u*t/u*)` -- has units of kg/m/s and is the HORIZONTAL saltation
flux, which needs a sandblasting efficiency to become a vertical dust emission.
Kok 18 gives the vertical flux directly and its coefficients are already fitted
and sourced in `aeolian/config/dust.yaml`.

Carry across with it:

- **The gravity term in the threshold.** DUST-6: Kok's standardized thresholds
  are Earth-fitted, and the correction is `(g/g_earth)^(1/4)` = 1.0691, the
  fourth root and not the square root. Both `u*st0` and `u*st_typical` take the
  same factor.
- **The subgrid wind distribution.** Emission is evaluated at quadrature points
  of a Weibull and mass-weighted, never at the gridbox mean; at `u*` cubed above
  a threshold those differ by orders of magnitude. Shape k = 2.012 from DUST-5.
  In-model this could instead use the model's own instantaneous `u*`, which is
  the entire point of going in-model -- but the gridbox is still 15 km and still
  needs a subgrid distribution, so the Weibull stays.
- **Snow suppression** from the model's own `dsnow`, against a depth threshold
  in the namelist. **Lake suppression is already inside field 1801**: the
  offline `source_fractions` zeroes the erodible weight wherever the solved lake
  extent stands, so the model does not need the lake fraction and could not
  reconstruct it from its own land mask anyway.

**The source injects a mixing ratio, not a flux.** Converting is
`d(mmr) = F * g * dt / dp` for the bottom layer, which is arithmetic, but it is a
change of kind from what `case(2)` does now and it is what made the old sink
intolerable.

**The namelist is the other half of this item.** `aero_namelist` is written by
ExoPlaSim's own Python API and this project had never touched it, so `ldepvel`,
`vdaero`, `lwetdep`, `scava` and `scavb` were unreachable from the drivers and
`aero_ini` would have aborted rather than run. `run_exoplasim.py:enable_dust_emission`
writes the whole group from `aeolian/config/dust.yaml` through the provenance
file the field generator leaves beside the `.sra` files, so the run cannot be
given a threshold that disagrees with the map it is applied to.

**And the shipped `aero_namelist` cannot be read at all.** It carries
`aerofile = 0`, an unquoted integer for a `character(len=80)`, and gfortran
rejects it with iostat 5010 while `aero_ini` reads without an iostat. Nothing had
noticed because `aero_ini` runs only at `L_AERO > 0` and every run this project
has made set it to zero -- the same reason the seven `aerocore` defects survived.
The first run to enable the aerosol would have died in the namelist read before
reaching any of this. `enable_dust_emission` rewrites the key.

**The emitted dust defaults to `l_aerorad = 0`.** The apart fix (defect 1) and
the longwave term are RESIDENT, and `AEROQLW` is now derived beside `DUSTQLW` in
the prescribed field's provenance, from the same optics and for the reason
recorded there, so the thermal term has a writer as well as a home. Enabling the
radiation is therefore `model.dust_emission_radiative`, taken deliberately
rather than defaulted, and the aerofile is staged and named on either setting.

### Predicted effects, stated before they run

The comparator is `aeolian/analysis/dust_baseline.json`, whose central-roughness
total emission was computed offline from the SAME climatology, the SAME map and
the SAME `dust.yaml`. Take the number from the file, not from here.

| change | prediction | falsified by |
| --- | --- | --- |
| `ldustemit = 0`, any existing configuration | output bit-identical: the legacy source line is untouched, nothing extra is gathered and no surface field is read | any difference at all |
| `ldustemit = 1` at the measured Weibull shape | total emission EXCEEDS the offline central figure, by 3x to 10x. The sign is Jensen's inequality on a convex flux and is not a guess: the offline chain's Weibull carries all the sub-daily variance while the model resolves that variance AND applies the same Weibull on top | a total BELOW the offline one, which no correct implementation can produce; check the latitude flip first |
| `ldustemit = 1` with `dustwk = 50`, the near-delta arm | within a factor of 5 of the offline central figure, either side. Collapsing the quadrature substitutes the model's own resolved variability for the fitted Weibull, which is what that Weibull was standing in for -- it was fitted to this model's own winds | a factor beyond 10, which is structural rather than distributional: units, the mixing-ratio conversion, or the threshold |
| doubling boundary field 1801 | total emission doubles, to round-off. An exact identity, and the entire argument for packing erodible fraction and clay fraction into one field | any ratio other than 2 |
| halving `deltsec` | emission per unit time does not move. The Kok flux carries no timestep and the `d(mmr) = F g dt / dp` conversion's `dt` cancels against twice as many steps | a rate that scales with the step, which means `deltsec` applied twice or a flux read as a concentration |
| the hemisphere | annual-mean bottom-level dust correlates with field 1801 more strongly than with its north-south mirror | the reverse, which is defect 7 reintroduced |
| the UNSCALED Earth thresholds | emission RISES, by about 30%, which is what the offline chain measured when the gravity term went in | the opposite sign, or more than a factor of two |

Two known small disagreements, recorded so they are not rediscovered as defects:
the offline chain uses dry-air `R = 287.05` for air density while the model uses
its declared `gascon`, and the offline reference height is built from surface air
temperature while the model uses bottom-level temperature. Both are a few percent
on `u*` and neither approaches the brackets above.

## The longwave term is mandatory, and it was the largest single piece

DUST-2: `radmod.f90` put the aerosol in bands 1 and 2 only and the longwave
solver had no aerosol term at all. Shortwave-only would apply -4.5 to -5.3 W/m2
of global-mean cooling against a true +0.35 to +0.74 -- about 4 to 5 K of
spurious cooling on the canonical flux-to-kelvin conversion in
`lib/sensitivity.py`.

**AUTHORED as `exoplasim-3.4.2-aerosol-longwave.patch`.** The physics was already
settled by `exoplasim-3.4.2-prescribed-dust.patch` item 5 -- a grey absorber, no
longwave scattering, at the same 1.66 diffusivity the cloud term uses, multiplied
into the total layer transmissivity so the overlap with water vapour and CO2 is
handled by construction -- and this patch adds one branch to that same term. It
changes no coefficient and no form.

What it took, and the reason it was priced as the dominant risk is that none of
it is confined to one statement:

- `aeroprof`, beside `dustprof` and called from the same place in `radstep`,
  builds `daerod` from `nrho`, `apart`, `qex1` and the layer thickness. That is
  the arithmetic `swr` used to do inline.
- `swr` now READS `daerod` instead of rebuilding it, so the shortwave and the
  longwave cannot drift apart. The LMD Generic PCM has the same structure: one
  per-layer optical depth at a reference wavelength, and every band including
  the infrared is that field times a per-band ratio held with the optics
  (`aerosol_opacity.F90`, `rad_correlatedk_ini_aerosol.F90`).
- `aeroqlw` in `radmod_nl`, separate from the prescribed path's `dustqlw`
  because the two paths carry different particles, with no default and a
  `radini` abort when the aerosol is on without it. And a second abort when both
  paths are enabled at once, because a prescribed column and a transported one
  are two aerosols whose optical depths would add.

The whole thing is a no-op at `l_aerorad = 0`, which is the default, so one
binary runs both arms of the A/B. The test with a right answer is the identity:
set `aeroqlw = dustqlw` and give the interactive path an `nrho` whose `daerod`
equals `ddustod` layer for layer, and the two arms' longwave fluxes must agree to
round-off, because they evaluate the same expression.

## The structural surprise: `aerocore` is serial

`aero_main` calls it inside `if (mypid == NROOT)`, having gathered the fields it
needs with `mpgagp`. So transport runs on one rank with the whole global grid.

Nothing the emission scheme needs was gathered. This was called the item most
likely to be underestimated, because it is invisible from the hook. **It came in
under that estimate, and the reason is worth recording: two of the three fields
it was budgeted for turned out not to be needed.**

- **The wind is already there.** `prepare_uvps` gathers `du` and `dv` and flips
  them in latitude on the way, so `zu(:,:,NLEV)` and `zv(:,:,NLEV)` inside
  `aerocore` are the bottom-level physical wind in the right hemisphere with no
  work at all.
- **The friction velocity must NOT be gathered.** `dtaux` and `dtauy` are the
  grid-cell surface stress, over the grid-cell roughness, and
  `aeolian/config/dust.yaml` argues at length that this scheme must not be fed
  that: MB95 partitions stress between a bed and centimetre-scale roughness
  ELEMENTS, and 15 km of orographic variance is not a roughness element. The
  scheme builds `u*` from the bottom-level wind through the aeolian roughness of
  the patch instead, which is what the offline chain does.
- **Air density and temperature** are `rhog` and `temp`, already computed inside
  `aerocore` for the settling term.

So the per-step cost is two `mpgagp` calls, `dsnow` and `dwatc`, both made only
when `ldustemit = 1`. The three boundary fields are read and gathered once at
initialisation.

What did NOT come in under estimate is the serialisation itself: emission is now
computed on NROOT over the whole global grid every timestep, 96 quadrature points
per source cell. That is bounded by the source-cell count rather than by the grid
(`dustsrc` cycles on `gsrcw <= 0`), which is about a sixth of land, but it is
still a per-timestep serial calculation inside a 16-rank run and it is the first
thing to look at if the model slows.

One consequence worth having: because `aero_ini` and `aerocore` both run on NROOT
only, the removal switches and the emission calibration need no broadcast.
`ldustemit` itself is the exception and does get one, because `aero_surf` runs on
every rank and they all have to agree about whether the collective read happens.

## Order, and scale

1. **Upstream defects.** Everything downstream is uncalibratable without
   them. RESIDENT: defects 2 to 7 in `exoplasim-3.4.2-aerocore-defects.patch`,
   defect 1 in `exoplasim-3.4.2-aerosol-apart.patch`.
2. **Replace the bottom-level sink** with a deposition velocity. RESIDENT,
   behind `ldepvel`.
3. **Wet scavenging.** RESIDENT, behind `lwetdep`. Both 2 and 3 are in
   `exoplasim-3.4.2-aerosol-deposition.patch`.
4. **Boundary fields and the emission law.** RESIDENT, behind `ldustemit`, as
   `exoplasim-3.4.2-dust-emission.patch` with
   `aeolian/scripts/build_dust_source_fields.py` and
   `run_exoplasim.py:enable_dust_emission`. The largest mechanical piece, and
   the one that had to be done for the switches added in 2 and 3 to be
   reachable at all: `aero_namelist` is written by ExoPlaSim's own Python API
   and nothing in this project had ever touched it, so `ldepvel`, `vdaero`,
   `lwetdep`, `scava` and `scavb` were unreachable from the drivers and the
   model aborted rather than run. The aerofile's provenance sidecar carries
   `namelist_values` with `APART` and `RHOP` so that the run gets its radius
   from the file the optics were built for rather than from a config that can
   drift.
5. **The longwave aerosol term.** The largest risk, inside the radiation
   solver. RESIDENT as `exoplasim-3.4.2-aerosol-longwave.patch`, written
   together with defect 1 on one branch because both live in `radmod.f90`.
6. **Size bins**, only if the single mode proves insufficient. DUST-8 says it
   does not, and says what that costs. CLOSED.

Honest scale: this is not an afternoon. It is five or six separate patches to a
compiled model, each needing a rebuild and verification per CLAUDE.md rule 4,
against a model whose aerosol path has demonstrably never been exercised -- three
defects found by inspection before touching it, and four more found while fixing
those three.

**Sequencing.** What `docs/src/pipeline/sequencing.md` A3 asks for applies
cleanly here. Each of the six pieces states its
predicted effect before it runs, and each is tested by a short A/B off a common
restart against that prediction rather than by an iteration of its own. The
no-op-until-enabled convention is what makes that possible, so every piece here
MUST default off, and the whole set can then ride one converged run. The defect
patch is the one exception and it is argued above: it changes nothing that any
existing configuration can reach, so its A/B is a bit-identity check rather than
a forcing measurement.

Two constraints on the ORDER are real and are unaffected. The upstream defects
come first because everything downstream is uncalibratable without them. And the
longwave term sits inside the radiation solver, so it cannot be developed in
parallel with any other patch to `radmod.f90` without the two being tested
together rather than independently, which `docs/src/reference/environment.md`
warns about. That is why defect 1 and item 5 waited for PHYS-6 and were then
written together, on one branch, against one base.

**A patch's place in the STACK is not its place in the ordering above**, and
conflating the two costs regeneration for nothing. Defect 1 is item 1 and its
patch applies FOURTH, after the deposition patch, because both edit the same two
lines of `aero_ini` and the deposition patch was written first. "Defects first"
is an argument about calibration -- nothing downstream can be measured against a
model that is wrong upstream -- and every one of these patches lands in the same
rebuild, so within that rebuild the order is mechanical. The cheap order is the
one that regenerates nothing already written.
