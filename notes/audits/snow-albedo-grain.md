# A grain-size axis for the modelled snow albedo, on this star

**Derived:** 2026-08-25, by `analysis/snow_albedo_grain.py`, from the three
grain-resolved snow reflectance arrays `vendor/exoplasim/exoplasim/plasim/src/specblock.f90`
compiles, the configured stellar spectrum, and `lib/stellar.py`'s band integral.
No World Orogen generation and no ExoPlaSim run: every number here is a spectral
integral over data already in the tree.

This is worldbuilding. Vesper is an invented super-Earth around a K2.5V host;
every quantity below is a property of that planet's simulated snow surface or of
the two-band shortwave scheme that reflects light off it. The reflectance
library is an Earth laboratory's measurement of a MATERIAL, and re-weighting a
material's reflectance for a different star is the whole of what this note does.

WORLD-DD48 asked for a two-band snow albedo on a grain-radius axis for this
star, derived rather than imported, with a check that can fail. This note is
that. It is an EXPLORATION product: nothing in the model was changed, and no
model term rests on it yet.

## What the model ships and what it uses

`specblock.f90` carries `fsnowalb`, `msnowalb` and `csnowalb` over the same 965
wavelengths every other endmember uses. They are the JHU becknic library's fine,
medium and coarse granular snow, and each library entry states an effective
grain size in its own header. **No code path selects among them.** The one line
that would have combined the three is commented out, and `snowalbedos`, the
array it would have written, is initialised to zero and read nowhere.

What `radmod.f90` weights instead is the `iceblend` family, and that family is
not a grain axis even though it looks like one.
`plasim/src/specs/combinedspec-snow_ice.ipynb` built each blend as a weighted
mixture of five component spectra -- clear ice, frost, and the same three snows
-- with the snow components held at fixed proportions and the CLEAR ICE fraction
solved by bisection against a declared broadband albedo target under a 5772 K
blackbody. So `iceblendmin` and `iceblendmax`, the endpoints `landmod`'s snow
albedo ramp interpolates between on surface temperature, differ in how much
clear ice is mixed in, tuned to hit a solar broadband number. They do not differ
in grain size, and the ramp between them is not an ageing law.

## The result: this star's snow is darker, and by how much

The three spectra weighted against the configured spectrum, and against the 5772
K blackbody the blends were fitted under, are in `analysis/snow_albedo_grain.json`.
The shape of it:

- In BAND 1 the two stars agree to within a few thousandths at every grain size.
  Snow reflectance is near unity below 0.75 um whatever the illumination, so
  moving flux around inside that band changes almost nothing.
- In BAND 2 this star's weighting is materially lower at every grain size,
  because a K dwarf's band-2 flux sits further into the infrared, where snow
  goes on darkening.
- In BROADBAND the two effects compound with the band-1 flux share, which is
  much smaller here than under the Sun, and the result is one-signed and large:
  the same snow, at the same grain size, is roughly a tenth of an albedo darker
  under this star. It grows with grain size.

**That is not a property of this world's snow. It is a property of this world's
light**, and it is the reason an Earth-tuned broadband snow albedo is wrong here
in a direction that is calculable. The material never changed; the flux moved
into the half of the spectrum where the material absorbs.

## The form, and why the coefficients are fitted rather than imported

`references/climber-x/src/smb/smb_surface_par.f90:snow_albedo_dang` implements
Dang, Brandt and Warren (2015): a quadratic in `rn = log10(r/r0)` with `r0` at
100 um, one coefficient triple per band and per illumination, and a dust
darkening whose black-carbon-equivalent loading scales as `(r/r0)**s`. That
exponent is why a dust coefficient cannot be grafted onto a grain-free ramp, and
it is why this row unblocks `dust-14`.

**The exponent's provenance, which the paper gives.** It is not a fitted
constant of the darkening law but the exponent that COLLAPSES the family: Dang's
equation (7) combines the impurity mixing ratio and the grain radius into one
predictor, and `s` is varied until the albedo-reduction curves for different
radii merge onto a single curve. The paper reports two values, one for the
all-wave and visible bands and a much smaller one for the near infrared, and
CLIMBER-X uses the visible one because it applies the darkening in the visible
only. That choice is the paper's own: it states that mineral dust does not
significantly change the near-infrared band albedo, because dust grains and snow
grains absorb similarly when averaged over 0.7 to 4.0 um. So `dust-14`'s term is
a VISIBLE-band term, and the visible is the band where this note's comparison
against the published coefficients came out small at both ends.

The paper also states a validity range for the collapsed form -- a range in the
combined predictor, and a floor on the radius below which the single-predictor
version should not be used -- and a ceiling on impurity loading. Those are
declarations `dust-14` inherits with the exponent, and they belong in the
implementation rather than here.

The published coefficients are a fit to Earth's solar spectrum, so importing
them would carry the Sun's band weighting into a K dwarf's radiation, which is
exactly what this row exists to avoid. The FORM is taken and the COEFFICIENTS
are fitted to the three star-weighted points. Three points and three
coefficients: the fit is an exact interpolation, not a regression, and what it
buys is a smooth differentiable statement of the same three numbers in the shape
a dust term can scale against.

## What the published coefficients are worth as a comparison, split rather than asserted

The obvious check -- does the solar arm reproduce the published values -- is not
a check as it stands, because the two do not average over the same wavelengths.
Dang's Table 1 bands are 0.3 to 0.7 um and 0.7 to 4.0; the model's are 0.34 to
the row boundary near 0.75 and that boundary out to 100. So each grain is
re-integrated over Dang's OWN band edges under the same blackbody, which splits
the difference into a band-definition part and a data part.
`band_edge_mismatch.per_grain` in the JSON carries both, per grain and per band.
The shape:

- In the visible both parts are small, a few thousandths each. The library's own
  spectra begin at 0.34 um, so the narrow 0.30 to 0.34 slice of Dang's visible
  band is held at the first measured value; snow is flat there and the slice is
  narrow, so this is stated rather than corrected.
- In the near infrared the band definition is worth about five hundredths, and
  almost all of that is the 0.70 to 0.75 um slice, which Dang counts as near
  infrared and the model counts as band 1, and over which snow is still bright.
  The 4 to 100 um tail the model's band 2 includes and Dang's excludes runs the
  other way and is much smaller, because under either star well under one per
  cent of the shortwave flux lands there.
- What remains after the band definition is removed is a difference between the
  JHU library's snow and Dang's model of pure snow. It is about a twentieth of
  an albedo, comparable to the band-definition term, and it FALLS SLOWLY WITH
  GRAIN SIZE: across the factor of seven in radius the three library entries
  span, it moves by about a fifth of itself.

**The constant term does not transfer and the grain-size slope does**: the
fitted near-infrared linear coefficient on both arms sits within a few
thousandths of Dang's published one, while the constant terms are several
hundredths apart. A dust term scales against the slope, so the part this row
needs is the part that survives the comparison, and the part that does not
survive it is the part the row was right to fit locally.

## Why the constant does not transfer, which the paper settles

Dang's pure-snow albedo is a DISORT computation on a semi-infinite
plane-parallel snowpack of Mie spheres. Its three ingredients, each named in the
paper's section 3:

- **The optical constants** are Warren and Brandt (2008)'s revised compilation of
  the complex refractive index of ice. The paper states what the revision was
  worth: its pure-snow albedos are slightly higher than the same group's earlier
  work because the revised absorption coefficient of ice is smaller between 300
  and 600 nm.
- **The grain shape is a sphere**, with Mie theory for radii 5 to 2500 um. A
  nonspherical crystal is represented by a collection of spheres of the same
  volume-to-area ratio, and a size distribution by its area-weighted effective
  radius. The paper cites the evaluations of that representation and states
  their result: good accuracy for the extinction efficiency and the
  single-scattering albedo, DISCREPANCIES FOR THE ASYMMETRY PARAMETER.
- **The geometry** is a deep opaque snowpack under a direct beam at a standard
  zenith cosine of 0.65, or under diffuse illumination for the overcast set,
  with incident spectra measured at the Arctic sea surface and extended to 4 um
  with a radiative transfer model.

Two of those three are AXIS CONVENTIONS rather than physics -- what counts as
the radius, and at what illumination angle the albedo is stated -- and the paper
supplies the transformation that makes their effect calculable. Its equation (5)
says a change of illumination angle can be MIMICKED by a change of grain radius,
which is to say that a geometry difference acts as a rigid shift of
`rn = log10(r/r0)`. Under a quadratic in `rn` a shift moves the constant term,
moves the linear coefficient by twice the quadratic coefficient times the shift,
and does not move the quadratic coefficient at all. Dang's own near-infrared
coefficients put the ratio of those two effects at about four: **any axis
offset, whatever caused it, lands four times harder on the constant than on the
slope.** That is the mechanism, and it is algebra rather than an observation.
`band_edge_mismatch.axis_leverage` in the JSON carries it for the two axis
questions this comparison actually has -- the library's stated near-normal
illumination against Dang's 49.5 degrees, and the library's unqualified
"effective size" if it is a diameter rather than a radius.

Applying equation (5) between the two illumination geometries removes a little
under half of what the near-infrared comparison had left unexplained, and it is
the largest identified term in it. What remains after both the band definition
and the geometry is the part with no axis representation: the asymmetry
parameter is the one single-scattering quantity a sphere representation is
documented to get wrong, and it is not degenerate with the radius, so a
difference in grain shape shows up as an offset that no rescaling of the axis
can absorb. It cannot be pushed further from this side, because the library
entries state that their 0.3 to 2.08 um portion was MODELLED and do not say with
what optical constants, what shape or what method; the introductory text that
would say is not held here.

## The zenith correction is the paper's, not CLIMBER-X's

`snow_albedo_dang` writes the effective-radius correction as
`r * (1 + a*(coszm - 0.65)**2)`, with the square on the angle difference alone.
Dang's equation (5), after Marshall (1989), is `r * (1 + a*(mu - 0.65))**2`,
with the square on the whole bracket. The two agree only at the reference angle.
Away from it CLIMBER-X's is flat to first order where the paper's is linear, and
it raises the effective radius for a low sun as well as for a high one, which
reverses the sign of the effect on the low side: a low sun makes snow brighter,
not darker, and only the paper's form says so. At the library's near-normal
incidence the two differ by about a factor of one and a half in effective
radius.

The coefficients here were read from CLIMBER-X's implementation before this
paper was held, and this is what holding it changed. It is failure-mode class 9
in the small: a form taken from an implementation of a citation rather than from
the citation. Nothing this project runs consumed the wrong form -- the
comparison in this note is its only reader -- and `analysis/snow_albedo_grain.py`
now carries the paper's.

## The checks, all fixed before any was run

1. **Reproduction.** `lib/stellar.py`'s band integral must be the arithmetic
   `radmod` runs, or the axis would be a second opinion rather than an
   extension. It is checked against all eight band pairs `radini` printed under
   "Finalized Albedos", which `analysis/ice_albedo.py` reproduces from
   `radmod`'s own loops to 1e-9. The bar was set at one per cent of a
   reflectance and the worst miss is four orders of magnitude inside it, which
   is the grid-resolution difference between integrating on the 965-point blend
   grid and interpolating onto the 2048-point stellar grid.
2. **Agreement of the three copies.** Each reflectance is read from
   `specblock.f90`, which is what the model compiles, and checked against
   `basespecs/*.npy`, which is what `surfacespecs.py` builds from, and against
   the library file both came from. All three agree to the library's last
   printed digit. The grain radii are PARSED from the library headers in the
   same pass, so a spectrum cannot be paired with another one's radius.
3. **Monotonicity.** Albedo must fall with grain radius in both bands and on
   both arms, because a larger grain is a longer path in ice per scattering
   event. It does, everywhere.

## What is read as a radius, and what that costs if it is wrong

The library states an "effective size" in micrometres and does not say radius or
diameter. It is read as a RADIUS, and the paper supports that reading twice
over: Dang's `r0` is a radius, defined as the area-weighted effective radius of
a distribution of spheres, and the paper states that the effective radii of
surface snow on Earth are rarely smaller than 30 um. Read as radii the three
library entries straddle that lower end and sit in an ordinary seasonal-snow
range; read as diameters they would all fall below it.

If the reading is wrong the three albedos are untouched and only the AXIS moves,
by `log10(2)`. What that costs is now arithmetic rather than an assertion, and
it is the same arithmetic as the geometry offset above: the constant term moves
by roughly the linear coefficient times the shift, the linear coefficient moves
by twice the QUADRATIC coefficient times the shift, and the quadratic does not
move. On the star's band 2 the quadratic coefficient is small, so a
radius-for-diameter error would move the near-infrared constant by about six
hundredths and the slope by under a thousandth. It is a statement about where
the axis is pinned, not about the shape of the curve on it.

## What this does not establish

The grain radius is an axis and not a prediction. Nothing in this model evolves
a snow grain, so a run is still told where on the axis its snow sits; how the
radius evolves is a surface mass balance's question and arrives with one, as
`exoplasim/notes/semi-scoping.md` sections 1 and 5 set out. The three points are
a laboratory's fine, medium and coarse and bracket seasonal snow, not firn and
not whatever metamorphism this world's own snow would reach.

## Sources

| Source | Status |
| --- | --- |
| `references/ecospeclib-all/water.snow.{fine,medium,coarse}granular.*.jhu.becknic.spectrum.txt` | **read** -- the three spectra with their stated effective grain sizes, their measurement geometry (directional-hemispherical at 10 degrees), and the statement that 0.3 to 2.08 um was modelled and 2.08 to 14 um measured |
| `references/climber-x/src/smb/smb_surface_par.f90` | **read** -- the implementation of Dang's form, with all four published coefficient triples and the zenith-angle effective-radius correction. This is the source the coefficients are quoted from |
| `references/dang2015-parameterizations-for-narrowband-and-broadband-albedo-of-snow.pdf` -- Dang, Brandt, Warren (2015). J. Geophys. Res. Atmos. 120, 5446-5468. `10.1002/2014JD022646` | **read** 2026-08-25, sections 1, 3, 3.1, 6 and 10 and Tables 1 and 3. SUPPLIED BY THE USER after open-access routes and Sci-Hub returned 403. It settles all three parts of what lay behind the near-infrared offset -- Warren and Brandt (2008) optical constants, Mie spheres under the volume-to-area equivalence with the asymmetry parameter called out as the quantity that representation gets wrong, and a semi-infinite DISORT geometry at a standard zenith cosine of 0.65 -- and its equation (5) supplies the transformation that makes an axis difference calculable, which is what turns the offset from fitted into understood. It also corrects two things read out of CLIMBER-X: the band edges of Table 1, and the FORM of the zenith correction |
