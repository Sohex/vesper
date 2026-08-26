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
darkening whose black-carbon-equivalent loading scales as `(r/r0)**0.73`. That
exponent is why a dust coefficient cannot be grafted onto a grain-free ramp, and
it is why this row unblocks `dust-14`.

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
Dang's bands are 0.2 to 0.7 um and 0.7 to 5.0; the model's are 0.34 to the row
boundary near 0.75 and that boundary out to 100. So each grain is re-integrated
over Dang's OWN band edges under the same blackbody, which splits the difference
into a band-definition part and a data part. `band_edge_mismatch.per_grain` in
the JSON carries both, per grain and per band. The shape:

- In the visible both parts are small, a few thousandths each.
- In the near infrared the band definition is worth several hundredths, and
  almost all of that is the 0.70 to 0.75 um slice, which Dang counts as near
  infrared and the model counts as band 1, and over which snow is still bright.
  The 5 to 100 um tail the model's band 2 includes and Dang's excludes is worth
  under a hundredth of the band, because under either star well under one per
  cent of the shortwave flux lands there.
- What remains after the band definition is removed is a difference between the
  JHU library's snow and Dang's two-stream model of pure snow. It is roughly
  twice the band-definition term in the near infrared and it is NEARLY CONSTANT
  IN GRAIN SIZE.

That last point is the useful one. **The constant term does not transfer and the
grain-size slope does**: the fitted near-infrared linear coefficient on both
arms sits within a few thousandths of Dang's published one, while the constant
terms are a tenth apart. A dust term scales against the slope, so the part this
row needs is the part that survives the comparison, and the part that does not
survive it is the part the row was right to fit locally.

The residual is not resolved here and this note does not claim to. The library
entries state that their 0.3 to 2.08 um portion was MODELLED and only the longer
wavelengths measured, so the two sides are two models of pure snow using
different optical constants, grain shapes and geometries, and separating those
would need Dang's own paper. It could not be fetched; see the reference row
below.

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
diameter. It is read as a RADIUS: that is what Dang's `r0` is, and it puts the
three points in an ordinary seasonal-snow range, whereas read as diameters they
would be finer than fresh snow. If that reading is wrong, every fitted LINEAR
coefficient shifts by `log10(2)` times itself and the three albedos are
untouched -- it is a statement about the axis, not about the values on it.

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
| Dang, Brandt, Warren (2015). *Parameterizations for narrowband and broadband albedo of pure snow and snow containing mineral dust and black carbon.* J. Geophys. Res. Atmos. 120, 5446-5468. `10.1002/2014JD022646` | **NOT OBTAINED.** Open-access routes and Sci-Hub all returned 403. It would have settled what optical constants, grain shape and geometry lie behind the near-constant near-infrared offset between the published coefficients and the JHU library, which is the one thing this note leaves open. Nothing here depends on it: the coefficients used for comparison are read from the CLIMBER-X implementation, and the deliverable is fitted locally rather than imported |
