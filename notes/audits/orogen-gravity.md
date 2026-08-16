# Audit: how gravity enters Orogen

*Audited 2026-08-16, against `vendor/orogen` at `cf-fork` 3907c58. Tasks in
`TASKS.md` under `GRAV`.*

Two questions were put: whether gravity being a post-hoc scaling factor warrants
fixing, and whether the unscaled bathymetric side is a bug. The answers differ.

## What the code actually does

`planet-params.js` sets `reliefScale = REFERENCE_GRAVITY_MS2 / gravityMS2`,
exactly 1 at Earth gravity. It is applied at the model-unit-to-kilometres
conversion and nowhere else. The whole pipeline -- tectonics, erosion, basin
detection -- runs in dimensionless model units. Two builds differing only in
gravity are bit-identical in every hash.

The rule is: **positive heights scale, negative ones do not.**

## Bathymetry: not a bug, and now stated once

Ocean depth is an isostatic balance between the water-plus-oceanic-crust column
and the continental column. Write the mass balance out and g multiplies every
term, so it cancels. Ocean depth is set by density contrasts and crustal
thicknesses; a planet at 2 g has the same ocean depth as one at 1 g. **Scaling it
would be wrong**, and the manifest already said so.

The real defect was not the physics but its distribution: the same rule was
written out at three call sites -- `js/basins.js`, `js/data-export.js`,
`tools/export-maps.mjs` -- and only one carried the reasoning. The site that
computes the exported `elevation_km` had no comment on the sign asymmetry at all,
where a bare `h > 0 ? h * reliefScale : h` reads like a sign bug to anyone who
does not already know why.

Fixed: `scaledHeightKm()` in `js/color-map.js` is now the single definition,
carries the argument, and all three sites call it. All 106 fork tests pass.

## Land: the leading-order amplitude is right, the texture is not

This is the part that does not fully hold up, and the intuition that gravity
should affect slopes and drainage is correct.

Orogen scales height and leaves the planform alone. Since `S = |grad h|`, that
scales every **slope** by 1/gamma too. At this world's gamma = 1.3058, slope
angles come out **23.4% gentler** than Earth's everywhere. Whether that is right
depends on which process is setting the slope, and the two disagree:

| control | prediction | at gamma = 1.3058 | Orogen gives |
| --- | --- | --- | --- |
| stream power at steady state, `E = K A^m S^n = U`, `K ~ rho g`, n = 1 | `S ~ g^(-1/n)` | x0.766 | x0.766 |
| same, n = 2 | `S ~ g^(-1/n)` | x0.875 | x0.766 |
| threshold hillslope, `tan(theta_c) = tan(phi)` | **gravity-independent** | x1.000 | x0.766 |

The hillslope row is the problem. Driving stress is `rho g h sin(theta)` and
resisting stress is `mu rho g h cos(theta)`; g cancels, so the angle of repose
does not depend on gravity -- which is why repose angles on Mars and the Moon are
much like Earth's. What does scale as 1/g is the critical slope *height*,
`H_c ~ c/(rho g)`.

So the physically expected high-gravity landscape has **similar slope angles at a
finer horizontal texture**: shorter hillslopes, more closely spaced ridges,
higher drainage density. Orogen instead produces the same landscape vertically
squashed. Both routes agree that relief amplitude falls roughly as 1/g, which is
why the leading-order result is defensible; they disagree about the shape.

A second, smaller point: the 1/g argument is about **strength-supported** relief.
Isostatically compensated land -- a plateau floating on thick crust -- is
gravity-independent for exactly the reason the ocean is. Separating the two needs
a crustal-thickness field the model does not carry, so the uniform 1/g on land is
an upper bound and over-suppresses plateaus.

## Does it warrant fixing?

**Not cheaply.** Gravity cannot enter erosion in a principled way while Orogen's
erosion constants are dimensionless -- there is no rho, no g and no real K to
attach it to. Giving them physical units is a different project from a bug fix,
and it would invalidate every build. `planet-params.js` also records that scaling
the model parameter *before* the nonlinear height curve was tried and damps the
effect badly: a 4% parameter change came out as 1.3% in height.

What this costs us, in the direction of the error:

- **Basin statistics.** A real high-g world should have smaller and more numerous
  closed basins. Ours are the same basins made shallower. Since the carve verdict
  is a depth-versus-water-balance test, this is not neutral. `GRAV-2`.
- **Drainage density**, which hydrography consumes, is gravity-invariant here and
  should not be.
- **Slope-dependent downstream work**: the pedology catena, desert pavement in
  the derived surface classes, and the glacier rough pass, which keys on
  `q99 - mean` relief and inherits the 23.4% directly. `GRAV-3`.

The recommendation is to carry this as a stated limitation with a known direction
rather than to change the generator, and to re-examine it if basin geometry
starts driving conclusions. That is `GRAV-1`.

## What is NOT affected

- Land/sea area split: scaling is monotonic and sea level sits at elevation 0, so
  the mask does not move.
- Ocean volume and sea level: bathymetry is unscaled and the basin is unchanged,
  so the water budget stays consistent.
- Basin ids and carve verdicts across a gravity change: the catalogue is
  bit-identical, so an existing verdict replays. This is already in `CLAUDE.md`.
