# The surface longwave: what emissivity is worth, and where the reflection went

Worldbuilding. Vesper is an invented planet. Everything below is about the
simulation of its surface and the radiation code that emits from it; every
number is a property of a laboratory spectrum, of `radmod.f90`, or of a run of
this model.

Two findings, both in `radmod.f90:lwr`, both about the same three lines.
Finding 1 is a defect in the solver and is fixed. Finding 2 is a measurement
that settles what shape the emissivity should have, and the answer is not the
shape the backlog assumed.

## Finding 1. The reflected surface longwave was never debited from the surface

`lwr` builds the surface emissivity and the surface emission

    zeps(:) = elwland*dls(:) + elwsea*(1.-dls(:))
    zbu(:,NLEP) = zeps(:)*zst4(:,NLEP)

and, at the foot of the routine, corrects for the surface not being black:

    zeps(:) = (1.-zeps(:))*dftd(:,NLEP)
    do jlev = 1,NLEV
      dftu(:,jlev)   = dftu(:,jlev)   - ztausf(:,jlev)*zeps(:)
      dftue2(:,jlev) = dftue2(:,jlev) - ztausf(:,jlev)*zeps(:)
    enddo

`ztausf` is dimensioned `(NHOR,NLEV)`, so that loop cannot reach `NLEP` and the
surface level never received the term. The reflected part of the downward
longwave was propagated up through every atmospheric level and was never taken
off the surface.

What reads the surface level is not a diagnostic. `dlwfl(:,NLEP)` is the net
surface longwave that `landmod` and `seamod` settle the surface energy budget
with and that output code 177 reports, so the modelled surface was absorbing
`(1-eps)*LWdown` that it had just reflected, while the lowest atmospheric layer
paid for it. The column conserved and the partition did not.

At `eps = 1` the term is identically zero, which is why nothing showed while
the land emissivity was a blackbody. It was never zero over water. `ELWSEA` is
0.98, so the modelled ocean and sea-ice surface has been carrying a spurious
0.02 of the downward longwave as absorbed energy for every run in the tree.

Fixed by applying the term at `NLEP` with a transmissivity of one, which is
what the transmissivity from the surface to itself is:

    dftu(:,NLEP)   = dftu(:,NLEP)   - zeps(:)
    dftue2(:,NLEP) = dftue2(:,NLEP) - zeps(:)

The fix is a no-op wherever the surface is black and changes the ocean surface
budget everywhere else. It is upstream behaviour, not something this fork
introduced.

## Finding 2. The land emissivity wants a better number, not a field

`config/planet.yaml` declared the land emissivity as a blackbody, at the value
the compiled literal carried, and recorded that a per-cell field out of the
lithology map was the honest answer. The field was built, measured and removed.
This is the measurement.

### The derivation

`analysis/rock_emissivity.py` gives one broadband thermal emissivity per Orogen
rock class: one minus the directional-hemispherical reflectance of the
ECOSTRESS spectra (Meerdink et al. 2019, which absorbs the ASTER library of
Baldridge et al. 2009), Planck-weighted over 8 to 14 micrometres. Only
hemispherical measurements are read, because Kirchhoff needs the whole
hemisphere and a bidirectional reflectance gives an upper bound on emissivity
rather than an estimate of it.

Three things about it are declared rather than assumed.

**Preparation.** Every class carries a bracket from its solid samples to its
particulate ones, because the reststrahlen bands that pull silicate emissivity
down are surface-scattering features of a coherent interface and weaken when
the material is broken up. A grid cell is a mixture of outcrop and regolith
that this project cannot measure, so the value is the midpoint and the ends are
the arm.

**Halite is held out.** NaCl has no absorption band anywhere in the thermal
window, so an optically thick pure halite powder measures as a poor emitter and
the library says so. A salt crust on a playa is neither optically pure nor
optically thick there. The library carries no crust sample, so the pure mineral
is the wrong preparation rather than a low answer, and the `evaporite` class is
derived from the evaporite minerals that are opaque in the band.

**The spectra stop.** They reach about 14 micrometres and a land surface at
these temperatures radiates most of its energy beyond that. The in-band value
is carried to the whole thermal spectrum; the opposing bound is a blackbody
outside the band, which multiplies every between-class contrast by the in-band
Planck fraction. That is the bound that compresses the signal, so it is the one
a verdict has to survive.

### The criterion, fixed before any spectrum was read

The surface net longwave is `-eps*(sigma*Ts^4 - LWdown)`, exactly linear in the
emissivity and in nothing else the surface carries, so a contrast of `d` in
emissivity is worth `d` times the model's own surface longwave loss, per cell,
with no run needed to find out.

A per-cell field must therefore buy, over a scalar already set at the field's
own land mean, an area-weighted standard deviation of more than **1.4 W/m2** in
that flux. 1.4 is the top of the model's own dry adiabatic energy sink, which
the filter hides and the surface fluxes pay for: a surface energy term smaller
than that sits inside the model's own non-conservation and cannot be claimed as
a resolved effect. The criterion is applied at both out-of-band bounds, and the
answer is the weaker of the two.

### The measurement

`analysis/emissivity_contrast.py`, run on the active build at T21 against the
bootstrap climatology (a bootstrap is not a baseline; what is taken from it is
a flux scale, not a state), measured on 2026-08-25:

| Quantity | Value |
| --- | --- |
| Surface longwave loss, land mean | 51.9 W/m2 |
| Lithology land-mean emissivity | 0.9357, arms 0.9191 and 0.9524 |
| Worth of the declared blackbody against that mean | 3.36 W/m2, one-signed |
| Field contrast over a scalar at the same mean, sd | 0.47 W/m2 |
| The same at the compressing out-of-band bound | 0.17 W/m2 |
| Largest single-cell contrast | 2.23 W/m2 |
| Land area above the 1.4 W/m2 criterion | 2.9 per cent |

The field fails the criterion by a factor of three at the optimistic bound and
by a factor of eight at the compressing one. The scalar's own error is seven
times the whole contrast the field would have bought, and it is one-signed
across every land cell rather than a redistribution.

So the finding the backlog recorded was right about the number and wrong about
the shape. A blackbody land surface is a real several-W/m2 error; the
variability around it is not resolvable here.

### What would reopen it

Two things, and both are measurable rather than arguable.

**Resolution.** The contrast is a spatial one and cell averaging removes most
of it. On the native mesh the area-weighted spread is 0.0236 in emissivity,
against 0.0077 once reduced to T21 -- a factor of three, worth about 1.2 W/m2
if a fine enough grid recovered it, which is at the criterion rather than under
it. `analysis/emissivity_contrast.py` takes a grid and re-runs the same test,
so a finer rung answers this rather than reopening the argument.

**The band.** A spectral library reaching past 20 micrometres would settle the
out-of-band bound and could raise the compressing answer to the carried one.
The library held here does not, for the hemispherical samples.

The contrast is mildly concentrated where it would matter most: cells that are
at least half closed-basin carry a mean absolute contrast of 0.65 W/m2 against
0.26 over the rest of the land, and the cells above the criterion are enriched
in endorheic area by about half again. That is the argument
`build_surface_roughness.py` makes for its own existence, and here it is not
enough on its own -- the absolute size on those cells is still at or under the
criterion.

`world-vhhs` carries the re-run; `world-38y` is closed on this measurement.
