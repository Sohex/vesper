# Putting dust emission in the model: what it actually takes

Design for DUST-3, written 2026-08-17 from source inspection, revised 2026-08-18
when items 1 to 3 were written as patches. The decision to go in-model is in
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
single-scattering albedo and backscatter ratio of `analysis/dust_optics.json` to
four decimals, and the optical depth round-trips against the offline chain's own
0.5122 exactly. Qext comes out 2.65 and 2.61, which is a plausible Mie efficiency
for a micron particle and an independent check that the mass-to-efficiency
conversion is right.

`exoplasim/patches/exoplasim-3.4.2-aerocore-defects.patch`, which is item 1 of
the ordering below, and
`exoplasim/patches/exoplasim-3.4.2-aerosol-deposition.patch`, which is items 2
and 3. Both are AUTHORED and VERIFIED but not applied: they sit in
`PENDING_PATCHES` in `exoplasim/scripts/rebuild_binaries.py` and move to
`RESIDENT_PATCHES` in the commit that applies them and rebuilds. Neither touches
`radmod.f90`, and neither file they do touch is touched by any other patch in the
stack, so they are independent of the resident set.

## Seven upstream defects, and none of them was reachable

All seven are in `aerocore.f90` and `aeromod.f90`. All are LATENT: every run this
project has made sets `L_AERO = 0` in `plasim_namelist`, `aero_main` runs only
when `nsela == 1 .and. nkits == 0 .and. l_aero > 0`, and so `aerocore` has never
executed. That is why they survived, and it is why fixing them needs no namelist
switch: on every configuration this project runs, a binary with the defect patch
must produce output BIT-IDENTICAL to one without it. That is the test to run
first, and it can fail.

Defect 1 is the exception and is still open, because it lives in `radmod.f90`.

1. **OPEN.** `radmod`'s `apart` is never populated from the namelist -- `aero_ini`
   use-associates only `l_aerorad` and `aerofile` -- so it stays at its 50e-9
   haze default while `aerocore` uses the namelist value. Optical depth comes out
   1/385 of intent at our effective radius. Deferred: `radmod.f90` is being
   edited on another branch for PHYS-6, and two patches to that file cannot be
   developed in parallel without being tested together.
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

Per `WORKFLOW.md` A3, with what result would mean "wrong".

| change | prediction | falsified by |
| --- | --- | --- |
| all seven, on any configuration that exists today | output bit-identical, because `aerocore` is unreachable at `L_AERO = 0` | any difference at all, which would mean the aerosol path is reachable and `notes/dust.md` is wrong about that |
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
threads arrays through `aerocore`'s whole argument list and through `radmod.f90`,
which is off limits while PHYS-6 is in flight. The design note recommended single
mode first and nothing since argues otherwise.

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

**So the configuration this project should end up in is `apart` = 2.21 um with
the Q values scaled by 2.251**, burden-matched because the in-model chain exists
for the AOD and the precipitation response (`WORKFLOW.md` A4), not for the
deposition field. Two things to carry with it:

- `dust_aerofile.py`'s plausibility check survives, but it has to be applied to
  the unscaled number: `Qext` at the optical radius is 2.65 and 2.61, which is a
  believable Mie efficiency for a micron particle, while the scaled value of 5.96
  is not an efficiency at all and must not be read as one.
- **Not applied yet, on purpose.** It only bites once `apart` reaches `radmod`,
  which is upstream defect 1, which is deferred. Applying it before then would
  change a tracked aerofile for no effect and collide with DUST-12.

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

The mechanism is `surfcode(code,'name')` to register, a module-level array, and
`mpsurfgp('name',array,NHOR,1)` alongside the existing calls at
`landmod.f90:312` onward. `code_surf_file` builds the filename as
`N%03d_surf_%04d.sra` from NLAT and the code, so at T42 that is
`N064_surf_1801.sra`. Registered codes run to 1741, so **1801-1803 is clear**.

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
- **Snow and lake suppression**, from the model's own `dsnow` and the lake
  fraction already in the land mask.

**The source injects a mixing ratio, not a flux.** Converting is
`d(mmr) = F * g * dt / dp` for the bottom layer, which is arithmetic, but it is a
change of kind from what `case(2)` does now and it is what made the old sink
intolerable.

**The namelist is the other half of this item.** `aero_namelist` is written by
ExoPlaSim's own Python API and this project has never touched it, so
`ldepvel`, `vdaero`, `lwetdep`, `scava` and `scavb` are unreachable from the
drivers today. Writing them from `aeolian/config/dust.yaml` belongs here, with
the rest of the configuration the emission scheme needs.

## The longwave term is mandatory, and it is the largest single piece

DUST-2: `radmod.f90` puts the aerosol in bands 1 and 2 only and the longwave
solver has no aerosol term at all. Shortwave-only would apply -4.5 to -5.3 W/m2
of global-mean cooling against a true +0.35 to +0.74 -- of order nine kelvin of
spurious cooling, against 21 W/m2 for the entire stellar sweep that produced a
33 K range.

What it needs: a longwave absorption optical depth per layer, built the same way
`aod1`/`aod2` are at `radmod.f90:1844`, and an emissivity term in the longwave
flux solver. The optics exist -- `dust_forcing.py` already computes a
Planck-weighted thermal-infrared absorption from the OPAC indices, which run to
40 um -- so this is plumbing rather than new physics, but it is plumbing inside
the radiation solver, which is the part of the model it is least comfortable to
touch. **Budget this as the dominant risk.**

The prescribed-dust patch already has a longwave dust term for its own path
(`exoplasim-3.4.2-prescribed-dust.patch`, item 5), so the physics is settled and
what is left is wiring the interactive `nrho` path into it.

## The structural surprise: `aerocore` is serial

`aero_main` calls it inside `if (mypid == NROOT)`, having gathered the fields it
needs with `mpgagp`. So transport runs on one rank with the whole global grid.

Nothing the emission scheme needs is currently gathered -- not friction velocity,
not soil wetness, not snow depth. Three more `mpgagp` calls, and the emission
computed on NROOT with everything else, or computed distributed and gathered.
Either way it is more plumbing than the hook suggests, and it serialises a
per-timestep calculation across a 16-rank run.

This is the item most likely to be underestimated, because it is invisible from
the hook.

One consequence worth having: because `aero_ini` and `aerocore` both run on NROOT
only, the removal switches need no broadcast. Anything the emission scheme adds
that is read outside `aerocore` will.

## Order, and scale

1. **Upstream defects.** Everything downstream is uncalibratable without them.
   AUTHORED as `exoplasim-3.4.2-aerocore-defects.patch`; defect 1 remains open
   and is blocked on the `radmod.f90` branch.
2. **Replace the bottom-level sink** with a deposition velocity. AUTHORED, behind
   `ldepvel`.
3. **Wet scavenging.** AUTHORED, behind `lwetdep`. Both 2 and 3 are in
   `exoplasim-3.4.2-aerosol-deposition.patch`.
4. **Boundary fields and the emission law.** The three `.sra` codes, the
   generator on the `aeolian/` side, `surfcode`/`mpsurfgp` registration, the
   gathers, the namelist writing, and Kok 18 in the source case. The largest
   mechanical piece, and the one that has to be done for the switches added in
   2 and 3 to be reachable at all.
5. **The longwave aerosol term.** The largest risk, inside the radiation solver.
   Blocked with defect 1 on the same branch.
6. **Size bins**, only if the single mode proves insufficient. DUST-8 says it
   does not, and says what that costs.

Honest scale: this is not an afternoon. It is five or six separate patches to a
compiled model, each needing a rebuild and verification per CLAUDE.md rule 4,
against a model whose aerosol path has demonstrably never been exercised -- three
defects found by inspection before touching it, and four more found while fixing
those three.

**Sequencing**, corrected 2026-08-18. This used to read that `WORKFLOW.md` A3
forbade landing in the same iteration as the carve, and A3 no longer says that.
The constraint was never real: attribution cannot gate a correct term, so it
bought nothing that a converged run's cost could be justified against.

What A3 asks for instead applies cleanly here. Each of the six pieces states its
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
together rather than independently, which CLAUDE.md's Environment section warns
about directly.
