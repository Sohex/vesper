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

## Surface ultraviolet: closed by decision, not by calculation

The biosphere half of this question is now settled by choosing rather than
computing. `star.surface_uv: earth_like` in `config/planet.yaml`.

The reasoning, and its honest status. A photospheric estimate puts UV-B at the
planet at about 0.45 times Earth's -- a cooler star emits far less at 300 nm,
partly offset by a geometric factor of 1.73 from a smaller star at a closer
orbit -- and a halved ozone column would transmit about 2.3 times as much, which
multiplies to roughly Earth's. **That near-cancellation is not a result.** It is
two rough estimates multiplying to about one, and the 0.5 column was chosen with
nothing behind it. It should not be quoted as a finding.

What is defensible is the choice it motivated. Surface UV is a product of stellar
emission and ozone absorption, and neither factor is known here: BT-Settl carries
no chromosphere, which is where a K dwarf's 200-320 nm flux mostly originates,
and the ozone column is prescribed as Earth's. Fixing the product means neither
factor has to be pinned separately, and it removes a free dimension rather than
adding one.

It is also physically self-consistent, because ozone is UV-produced: a more
active star makes more ozone, which absorbs more of the extra ultraviolet. The
column tracks the incident flux and buffers the surface. That is what allows this
world to carry a dramatic bolometric activity cycle and a terrestrial surface UV
environment at once, where those two would otherwise pull against each other.

**This does not settle ozone's radiative effect**, which is a different quantity
from its shielding. The column that produces Earth-like surface UV under a
weaker, redder star is probably not Earth's column, and its effect on
stratospheric heating is still unmeasured. The `o3scale` sensitivity test
stands.

Nothing goes into LPJ-GUESS. It models no ultraviolet, adding a damage term would
need a calibration this project cannot supply, and the answer above is that the
surface environment is terrestrial anyway.
