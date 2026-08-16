# Ozone: what the model does, and why it is wrong here in a known direction

`config/planet.yaml` sets `atmosphere.ozone: true` and has said nothing else
about it. This records what that actually selects.

## What ExoPlaSim does

`no3 = 1`, the default, computes a **synthetic Earth ozone distribution** from six
fitted parameters in `radmod.f90`:

    column = a0o3 + a1o3*|sin(lat)| + aco3*sin(lat)*cos(2*pi*(day - toffo3))
           = 0.25 + 0.11*|sin(lat)| + 0.08*seasonal

That is Earth's ozone column, roughly 250 to 360 Dobson units, distributed
vertically through `bo3 = 20000` and `co3 = 5000`. It is prescribed, not
computed: nothing about the star enters it. `no3 = 2` would read a climatology
file instead, which is also Earth's.

## What it gets right

The shortwave absorption is not blind to the spectrum. Ozone transmissivity
divides by `zsolar1`, the fraction of stellar flux below 0.75 microns, which
`k25v` correctly sets to 0.384 against the solar 0.517. So the scheme accounts
for how much of the flux arrives in the band ozone absorbs in.

## What it gets wrong, twice, in the same direction

**The column is Earth's.** Ozone is produced by O2 photolysis in the far
ultraviolet and destroyed by photolysis at longer ultraviolet wavelengths, so the
equilibrium column depends on the host star's UV-to-visible ratio. A K2 dwarf has
a *lower* UV/visible ratio than the Sun, and modelled ozone abundance rises with
that ratio, so an Earth-like O2 atmosphere here should carry **less** ozone than
Earth, not the same.

**The absorption is shaped for the Sun.** Dividing by `zsolar1` corrects for the
band-1 total but not for the distribution *within* band 1. Ozone absorbs in the
Hartley and Huggins bands, 200 to 350 nm, and a 4965 K star puts proportionally
far less of its band-1 flux there than a 5772 K one does. The absorption
coefficients were fitted against a solar-shaped band 1, so absorption per unit
ozone is overestimated.

Both errors overstate stratospheric heating. They do not cancel.

## What to do about it

`o3scale` is a direct multiplier on the ozone field, in `radmod_nl`. Like
`nenergy`, it is not exposed by ExoPlaSim's Python API, so setting it needs the
same direct namelist edit `run_exoplasim.py` already performs for `STARFILE` and
`NENERGY`.

Not yet changed, because the magnitude is unknown to this project: the direction
is established but no number here is defensible without a photochemical model or
a literature value for a K2V host specifically. What is cheap is a sensitivity
test -- one T21 run at `o3scale` 0.5 against the baseline -- which would say
whether this is worth caring about for surface climate at all. In a ten-layer
model there is not much stratosphere to heat, so the surface effect may well be
small even though the bias is real.

## Surface ultraviolet: measured, not chosen

Settled by Segura et al. (2003), "Ozone Concentrations and Ultraviolet Fluxes on
Earth-Like Planets Around Other Stars", Astrobiology 3, 689-708, read directly.
It models a K2V host explicitly, which is very nearly this star, with an
Earth-like O2 atmosphere.

    ozone column, 1 PAL O2    6.64e18 cm-2   against the Sun's 8.36e18   = 0.794x
    surface UV-B              "about 0.4 times Earth's flux"

So this world is **substantially better protected from ultraviolet than Earth**,
not comparably. The star's lower ultraviolet output beats the thinner ozone
column it produces.

`config/planet.yaml` now carries `model.ozone_scale: 0.794` and
`star.surface_uv_relative_to_earth: 0.4`, and `run_exoplasim.py` sets `O3SCALE`
through the namelist directly, as it already does for `NENERGY` and `STARFILE`.

### What this replaced, and why it was wrong

An earlier version of this file closed the question by *choosing*: declaring
surface UV Earth-like so that neither the stellar emission nor the ozone column
had to be pinned. The arithmetic behind it took photospheric UV-B at 0.45x
Earth's and a halved ozone column transmitting 2.3x, multiplying to 1.05x.

Both inputs were wrong and the answer was wrong by a factor of 2.6. The column is
0.794x rather than 0.5x, and the surface flux is 0.4x rather than 1.05x. The
cancellation that motivated the choice was two rough estimates landing near one,
which was flagged at the time as not being evidence -- correctly, and the flag
should have been treated as a work item rather than a disclaimer.

What survives is the mechanism: ozone is UV-produced, so the column does track
incident ultraviolet and does buffer the surface. It simply does not buffer it
all the way to parity.

### Consequences

- **The biosphere faces less ultraviolet stress than Earth's**, so a UV damage
  term in LPJ-GUESS remains unnecessary -- now for a measured reason rather than
  an assumed one.
- **The stellar activity level is no longer pinned by this.** The earlier
  reasoning made activity a free parameter fixed by requiring Earth-like surface
  UV. With the column and the surface flux both measured for a quiescent K2V,
  activity returns to being undetermined, and the stellar-cycle block in
  `planet.yaml` stays PROVISIONAL on its own merits.
- **The radiative half is still only half corrected.** `ozone_scale` fixes the
  column. The absorption coefficients are still fitted to a solar-shaped band 1,
  and a 4965 K star puts proportionally less of its band-1 flux in the 200-350 nm
  Hartley and Huggins bands, so absorption per unit ozone remains overestimated.
  The sensitivity test now measures a known correction rather than a guess.

Nothing goes into LPJ-GUESS. It models no ultraviolet, and the surface
environment is measured as gentler than Earth's.


## The photosphere/chromosphere gap, and how to close it

Worth stating precisely, because it is an input gap rather than a model defect
and it has three compounding layers.

A star's ultraviolet comes from two places. The **photosphere** radiates roughly
as a blackbody at the effective temperature, and its output collapses toward
short wavelengths on the Wien tail. The **chromosphere** is hot, magnetically
heated plasma above it, at 10^4 K and more, radiating in emission lines rather
than a continuum. For the Sun at 5772 K the photosphere still supplies useful
flux at 250-300 nm. At 4965 K it does not, so the chromosphere supplies
proportionally far more of a K dwarf's ultraviolet than of the Sun's.

Our spectrum has none of it:

1. `k25v` is built from **BT-Settl**, a photospheric model in radiative-convective
   equilibrium. It has no mechanism that produces a chromosphere, so it has no
   chromospheric emission, by construction rather than by omission.
2. The file starts at 0.34 microns regardless.
3. `radmod.f90:226` zeroes flux below 0.316 microns on the star-file path.

ExoPlaSim would consume ultraviolet if it were given any. Nothing is broken; the
information was never supplied.

### Where it actually bites

Only in one place that matters: the Hartley-Huggins weight in the ozone
absorptance, `model.ozone_uv_weight`. Everything else is either insensitive --
ultraviolet is a few percent of total flux even for the Sun, so the energy budget
barely notices -- or already answered from a proper photochemical model, as the
column and the surface flux now are.

### The correction, and it is a known construction

Use an **observed** ultraviolet spectrum instead of a modelled photosphere, and
splice it onto the photospheric model where the two overlap. That is exactly what
Segura et al. did: coadded IUE observations of epsilon Eridani from 115 to 335 nm,
merged onto a Kurucz photosphere at 320 nm. `build_stellar_spectrum.py` already
performs the analogous merge for the BT-Settl grid, so the machinery exists.

Two sources would serve: the IUE archive for epsilon Eridani directly, or the
MUSCLES survey, which publishes panchromatic spectra of low-mass stars including
the ultraviolet.

### What this tells us that is not a caveat

Epsilon Eridani is young and chromospherically active, and Segura notes a quieter
K2V would receive about 2.5 times less ultraviolet. So the column of 0.794x and
the surface flux of 0.4x describe an **active** K2V, while our 0.469 UV weight
describes a star with no chromosphere at all. The two bracket the answer, and
what sits between them is set by the stellar activity level -- which this project
has left PROVISIONAL and undecided.

So the ozone numbers, the surface ultraviolet, and the stellar cycle amplitude
are one decision, not three. Deciding the activity level resolves all of them,
and until it is decided the honest form is a bracket rather than a value.
