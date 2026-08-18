# Putting dust emission in the model: what it actually takes

Design for DUST-3, written 2026-08-17 from source inspection. Nothing here is
applied. The decision to go in-model is in `notes/dust.md` and rests on DUST-2:
the reopening test crosses by 3.8x, and a prescribed field cannot respond to the
winds the dust itself changes.

**The headline is that "replace one line" was wrong, and by a lot.** The hook at
`aerocore.f90:633` is one line, but around it sit a serial-on-NROOT transport
call, a resolution-dependent sink, a shortwave-only radiation scheme, no wet
deposition, one hard-wired particle size, and three upstream numerical defects.
The line is the smallest part.

## What is already done

`exoplasim/data/dust/vesper_dust_aerosol.dat`, from
`exoplasim/scripts/dust_aerofile.py`. Validated: `radmod` reconstructs the
single-scattering albedo and backscatter ratio of `analysis/dust_optics.json` to
four decimals, and the optical depth round-trips against the offline chain's own
0.5122 exactly. Qext comes out 2.65 and 2.61, which is a plausible Mie efficiency
for a micron particle and an independent check that the mass-to-efficiency
conversion is right.

## The sink at `aerocore.f90:869` is a stability hack, not a parameterisation

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

It survives in the shipped code because the existing dust source *sets* the
bottom-level mixing ratio to `fcoeff*land` every step rather than adding a flux,
so the two degenerate together into "surface concentration is pinned, everything
above is what leaks up". That degeneracy is exactly what a real flux-based source
breaks.

**It must be removed and replaced with a deposition velocity** before any
emission is calibrated, because it multiplies whatever is injected by 0.01 every
step. Replace with `mmr = mmr * exp(-v_d * dt / dz)` and a `v_d` from the size
bin, which is the same Stokes route `settling_velocity` already uses offline.

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
change of kind from what `case(2)` does now and it is what makes the sink above
intolerable.

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

## Wet scavenging

DUST-7: `aerocore`'s argument list contains no precipitation field of any kind.
Gravitational settling is the only removal.

`aeolian/README.md` says of the offline chain that without the wet term the fine
mode's lifetime is out by an order of magnitude, so this is a first-order sink
and not a refinement. It needs `dprl` and `dprc` passed in and a loss rate
proportional to them, which is the Sportisse (2007) route the offline chain
already implements. Mechanically the easiest of the physics additions: the fields
exist, they just are not handed over.

## Size

DUST-8: `plasimmod.f90:99` fixes `parameter(NAERO = 1)` at compile time, and the
`do jc=1,NAERO` loops are already written -- but `apart` and `rhop` are scalars
shared by every bin, `mmr2n` takes them as scalars, and the optics are
single-valued.

Two routes:

**Single mode.** Ship what `dust_aerofile.py` already produces: an effective
radius of 0.98 um reproducing the optical depth exactly. Optically right, and
mass-weighted wrong -- the emitted volume-median diameter is 3.4 um, so the
coarse fraction settles far too slowly and the burden stays up too long. The
deposition field loess and the phosphorus leg want would be wrong, which is why
the offline chain has to survive regardless.

**Bins.** Raise NAERO, promote `apart` and `rhop` to arrays, extend `aerofile`
to per-bin optics, and give `mmr2n` and the settling term per-bin values. More
invasive than it sounds because the arrays thread through `aerocore`'s whole
argument list.

Recommend single mode first, with the deposition field left to the offline chain
and that stated, and bins as a later iteration if the deposition matters enough.

## Three upstream defects are prerequisites

Reported upstream; all three bite as soon as `l_aerorad = 1`.

1. `radmod`'s `apart` is never populated from the namelist -- `aero_ini`
   use-associates only `l_aerorad` and `aerofile` -- so it stays at its 50e-9
   haze default while `aerocore` uses the namelist value. Optical depth comes out
   **1/385** of intent at our effective radius.
2. `aerocore.f90:1156` and `:1183`: `(4/3)` is integer division, so sphere volume
   is 4/3 too small and number density 33% high.
3. `aerocore.f90:949`: `(4/25)` is integer division, so the collision-integral
   temperature dependence vanishes and viscosity is a fifth to a quarter low,
   making Stokes settling correspondingly fast.

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

## Order, and scale

1. **Upstream defects 1-3.** Everything downstream is uncalibratable without
   them. Small, mechanical, and they are one-character changes plus a namelist
   read.
2. **Replace the bottom-level sink** with a deposition velocity. Small, and it
   has to precede any calibration.
3. **Wet scavenging.** Fields exist; pass and apply. Moderate.
4. **Boundary fields and the emission law.** The three `.sra` codes, the
   generator on the `aeolian/` side, `surfcode`/`mpsurfgp` registration, the
   gathers, and Kok 18 in the source case. The largest mechanical piece.
5. **The longwave aerosol term.** The largest risk, inside the radiation solver.
6. **Size bins**, only if the single mode proves insufficient.

Honest scale: this is not an afternoon. It is five or six separate patches to a
compiled model, each needing a rebuild and verification per CLAUDE.md rule 4,
against a model whose aerosol path has demonstrably never been exercised -- three
defects found by inspection in code that was presumably believed to work. Expect
to find more once it runs.

**Sequencing** stands as `WORKFLOW.md` A3: this moves land-surface forcing, so it
must not land in the same iteration as the carve. The natural slot is its own
iteration after the baseline re-run.
