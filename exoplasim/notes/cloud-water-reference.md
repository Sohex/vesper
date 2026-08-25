# The cloud water reference, read from CCM3 and bracketed for arms

**Context first, because this document can be read standalone: this is
WORLDBUILDING.** Vesper is a fictional super-Earth. Everything below is about
one constant in a toy GCM's diagnostic cloud scheme and what a simulated
atmosphere does when that constant moves. Nothing here is a statement about
the real world.

`rainmod`'s `mkclouds` gives the modelled in-cloud liquid water density an
exponential profile, `rho_l(z) = clwref * exp(-z/hl)`, and converts it to a
mixing ratio by the local air density before `radmod` consumes it. `clwref` is
`rho_l0`, the reference density at the surface, and `world-ofn` exposed it as a
`rainmod_nl` key at CCM3's own value while deriving its neighbour `clwhsc` from
the model's own `gascon` and `ga`. The direction was settled there: `clwref` is
a condensate concentration set by microphysics, carries no gravity dependence,
and so is not scaled to this world. What was missing was a bracket, and
`world-8h6` was blocked on the source.

## What the sources settle

Both primary CCM3 sources are now in `references/` and read.

**Kiehl et al. (1998), J. Climate 11, 1131-1149, p. 1132-1133.** Eq. 3 is the
exponential profile and states "the reference value `rho_l0` is equal to 0.21 g
m-3". Eq. 4 is the scale height, `hl = a ln(1 + (b/g) integral q dp)`, with "the
parameters have been empirically determined to be a = 700 m and b = 1 m2 kg-1".
Cloud optical properties come from Slingo (1989) for liquid and Ebert and Curry
(1992) for ice, with a droplet effective radius that differs over ocean and
land. Longwave cloud emissivity is Eq. 12-13, `A_c' = eps * A_c` with
`eps = 1 - exp(-D k_abs CWP)` and `D = 1.66`.

**Kiehl et al. (1996), NCAR/TN-420+STR, p. 49-50.** The same profile as
Eq. 4.a.11 and the same scale height as Eq. 4.a.12, plus Eq. 4.a.14, the
ANALYTIC layer integral CCM3 uses to turn the profile into a layer cloud water
path. It also gives the value's provenance, and that is the part the journal
paper does not: the profile is an inheritance from CCM2, where in-cloud liquid
water paths "were evaluated from a prescribed, meridionally and height varying,
but time independent, cloud liquid water density profile ... which was
analytically determined on the basis of a meridionally specified liquid water
scale height". CCM3 kept the reference value and replaced only the meridional
scale-height specification with the precipitable-water diagnosis.

**Neither source gives a range, an uncertainty, or a sensitivity for
`rho_l0`.** It is a single number in both, and the technical note says it
anchors a prescribed analytic profile rather than a measurement. There is no
bracket to inherit. That is the paper's answer to `world-8h6`, and it is a real
answer rather than a failure to find one.

The one nearby quantity either source does put against observations is the
global annual mean cloud water path: Kiehl et al. (1998) Table 1 gives CCM3 at
0.0465 mm against 0.0813 mm observed (Greenwald et al. 1995). CCM3's scheme, at
`rho_l0` = 0.21 g m-3, therefore delivered about 0.57 of the observed global
mean. That is a comparison, not a bracket, and it is not used here to choose a
value; it is recorded because it says which end of the bracket the value's own
authors' comparison sits toward, and it is written down before any arm is run.

## Does this fork's formulation match the one the value came from

Checked against the code rather than against the documentation, because a
sibling batch found this fork's published methods citing the wrong paper for a
different constant.

**The vertical distribution matches.** `mkclouds` builds `dqvi` as
`sum(q dsigma) * dp / ga`, which is the precipitable water Eq. 4 integrates,
and `zzh = clwhsc * ALOG(1 + dqvi)`, which is Eq. 4 with `b` = 1 m2 kg-1
implicit and `a` = `clwhsc`. The default `clwhsc` derives from `gascon/ga` and
returns 700 m exactly at Earth's values, so an Earth configuration reproduces
CCM3. `zzf` is a hypsometric height above the surface, which is the `z` in
Eq. 3.

**The layer integral differs, and the difference does not matter.** CCM3
integrates the exponential across the layer analytically (Eq. 4.a.14);
`mkclouds` evaluates `rho_l` at the mid-layer height and multiplies by
`dp/ga * dsigma`. Midpoint evaluation of a convex exponential is low. Computed
on the `NEQSIG = 4` sigma set at ten layers with 25 kg/m2 of precipitable water
and a 6.5 K/km profile from 288 K, the ratio of the two is 0.999, 0.994, 0.987,
0.977, 0.966 and 0.951 from the surface layer upward, and only leaves 1.3 per
cent in layers carrying more than 3 g/m2. Not a defect.

**The longwave consumer matches.** `radmod`'s
`ztaucc0 = 1 - dcc*(1 - exp(-1.66*acllwr*1000*dql*dp*dsigma/ga))` is Eq. 12-13
with `D` = 1.66 written as a literal, `acllwr` in the place of CCM3's `k_abs`,
and the in-cloud water path in g m-2. CCM3 weights `k_abs` between liquid and
ice phase by an ice fraction; this model has no phase split, so `acllwr` is
used unweighted.

**The shortwave consumer does not match, and this is the finding.** CCM3 puts
the cloud water path into Slingo (1989) delta-Eddington, in which extinction
optical depth is LINEAR in the path and inverse in a droplet effective radius
(technical note Eq. 4.b.3). `radmod`'s `swr` says of itself that its cloud
transmissivities are "Stephens (1978) + Stephens et al. (1984)", and what it
computes is

    ztau = 2.0 * ALOG10(zlwp + 1.5)**3.9

with no effective radius anywhere in the model. The two schemes have different
curvature in the path. The elasticity `d ln ztau / d ln CWP` of the Stephens fit
is 1.28 at 0.4 g m-2, peaks near 1.8 around 2 g m-2, passes 1.0 near 45 g m-2,
and falls to 0.84 at 94 and 0.73 at 200; Slingo's is exactly 1 everywhere. So
the fit is SUPER-linear in thin cloud and sub-linear in thick, and
`opaque-constants.md` finding 8 and `world-8h6` were both half right in calling
the response sublinear: that holds only for the thick low-cloud layers.

The consequence for this issue is direct. `rho_l0` was fixed inside a scheme
where optical depth is linear in the water path; it is consumed here by one
where it is not. Whatever cloud water path CCM3 chose 0.21 g m-3 to deliver,
the optical depth that path produces in this model is not the optical depth it
produced in CCM3, and the value carries no calibration across the boundary.
That is an argument for measuring the lever rather than trusting the inherited
number, which is what the arms below do.

CCM3 further splits the condensate into liquid and ice by temperature and gives
each its own optics, and modifies the layer optical depth by a power of the
cloud fraction for overlap; this model does neither. Those are separate
structural gaps, recorded here so they are not re-derived, and they are not
`world-8h6`'s to close.

## The bracket

**Bracket: `clwref` in [0.000105, 0.00042] kg m-3, a factor of two either side
of the CCM3 value. Fixed before any arm is run.**

**The bracket's source is not the value's source, and the two are different
kinds of quantity.** The value is the reference of a prescribed analytic
profile. The bracket is taken from what in-cloud liquid water content is
observed to be, because that is the only observable the constant is
commensurate with.

The profile puts its reference density at the surface and decays it with a
scale height that is about 1.7 km at 25 kg/m2 of precipitable water on this
world, so nearly all of the water sits in the lowest few kilometres and the
reference is a low-stratiform in-cloud density. The width of the bracket is
what the observations set, and two independent facts fix it.

**In-cloud liquid water content varies by more than a factor of two within a
single cloud field, so a globally uniform reference cannot be pinned tighter
than that by observation.** Lloyd et al. (2018) measure peak in-cloud liquid
water content across profiles in the same stratocumulus deck at 0.8 to 1 g m-3
in the overcast region (p. 17194), around 0.8 g m-3 in a second case's
stratiform region, 0.7 g m-3 in a third, and 0.2 to 0.3 g m-3 in profiles near
the deck's break-up. Those are in-profile peaks near cloud top rather than
profile means, so the corresponding profile-mean densities are smaller by
roughly the adiabatic factor, which puts 0.21 g m-3 in the middle of the range
rather than at an end. Covert et al. (2022) find liquid water path in one
simulated stratocumulus domain spanning below 25 to over 300 g m-2 (p. 1165),
a factor of twelve inside one deck.

**CCM3's own comparison points toward the upper end, and by less than the
bracket allows.** Table 1 of Kiehl et al. (1998) puts CCM3's global mean cloud
water path at 0.57 of the observed value beside it, so 1.75x would be needed to
close that gap. That figure is inside the 2.0x arm rather than outside it,
which is the check that the bracket is wide enough to contain the answer.

**It is not used to set the value.** `clwref` stays at CCM3's value in the
baseline configuration and only the arms move it. A value fitted to improve an
agreement would be a knob, and this constant is a physical property.

## The arms, and what each predicts

Three points, exactly as PHYS-11 is run: `clwref` at 0.000105, at the compiled
0.00021 as control, and at 0.00042, same binary and same restart. The 1.0x arm
must be bit-identical to the control, because the key defaults to the value the
binary already carries.

The table below is the offline prediction, computed on the `NEQSIG = 4` ten
layer sigma set with 25 kg/m2 precipitable water and a 6.5 K/km profile from
288 K, using this fork's own `clwhsc` derivation, `tswr1`, `zmu00` = 0.5 and
`acllwr`. `R` is the band-1 cloud reflectance `1 - 1/(1 + zb3*ztau)` the
shortwave scheme builds from `ztau`; `eps` is the longwave cloud emissivity.

| layer, height | CWP g m-2, 0.5x / 1x / 2x | R, 0.5x / 1x / 2x | eps, 0.5x / 1x / 2x |
| --- | --- | --- | --- |
| 9, 0.5 km | 47 / 94 / 188 | 0.625 / 0.757 / 0.844 | 1.000 / 1.000 / 1.000 |
| 8, 1.3 km | 46 / 92 / 184 | 0.621 / 0.754 / 0.842 | 1.000 / 1.000 / 1.000 |
| 5, 5.3 km | 8.8 / 17.6 / 35.2 | 0.187 / 0.364 / 0.556 | 0.768 / 0.946 / 0.997 |
| 4, 7.3 km | 3.5 / 6.9 / 13.9 | 0.050 / 0.139 / 0.298 | 0.437 / 0.683 / 0.900 |
| 3, 9.7 km | 1.1 / 2.2 / 4.3 | 0.007 / 0.023 / 0.071 | 0.164 / 0.301 / 0.512 |

**The sign IS assignable here, and the reason is in the table.** The modelled
longwave cloud emissivity of the low layers is already saturated at 1.000 in
all three arms, so low cloud's longwave effect does not respond to `clwref` at
all; only its shortwave reflectance does, and that moves by +0.09 on the 2.0x
arm and -0.13 on the 0.5x arm. The low-cloud response is therefore purely
shortwave and one-signed. The longwave gain opposes it but is confined to the
mid and high layers, where the water path is one to two orders of magnitude
smaller. **Predict: the 2.0x arm cools the simulated global mean and the 0.5x
arm warms it.** This is what separates this term from PHYS-11, where two
shortwave effects opposed and the sign could not be assigned.

**Predicted magnitude: both arms move the simulated global mean by more than
2 K, and the plausible range is 5 to 10 K.** The estimate: modelled low cloud
band-1 reflectance moves 0.09 absolute on the 2.0x arm; at a low cloud fraction
of order a third and this world's mean insolation, that alone is several W m-2
at the top of the atmosphere, and PHYS-11 measured 3.2 W m-2 as worth 2.6 K in
this model. The bound is deliberately wide because the offline estimate that
priced PHYS-11 was falsified by a factor of five, and it was falsified in the
direction of the model responding MORE than the estimate allowed.

**What would mean wrong, in the A/B:**

- The 1.0x arm differing from the control at all. The key defaults to the
  compiled value.
- Either arm moving the simulated global mean by less than 1 K. The reflectance
  table does not allow it, and it would mean `CLWREF` is not reaching
  `rainmod`. Check the run's namelist and `rainini`'s echo before believing it.
- The 2.0x arm warming, or the 0.5x arm cooling. That contradicts the low-cloud
  saturation argument above and points at the implementation, not at the
  physics.

**If an arm fails to equilibrate, or crosses the sea-ice transition, run the
0.75x and 1.33x pair on the same argument and read the slope from those.** A
response that large is not a failure of the arm; it is the result, and it says
the cloud water reference is a larger lever in this model's radiation than any
term in the forcing bundle. The project's `lib/sensitivity.py` slope is LOCAL
and bracketed over about 7 K, so an arm this size is read from the run's own
mean and never converted through that slope.

## What is not settled, and what would settle it

The arms have not been run. The route exists: `config/planet.yaml` declares
`model.cloud_water_reference_kg_m3`, `run_exoplasim.py` writes
`CLWREF@rainmod_namelist` unconditionally the way it writes
`TSWR3@radmod_namelist` for PHYS-11, `continue_exoplasim.py` reapplies it per
segment, and `expected_namelist_keys` covers it, so a run that dropped the key
fails its own staging check rather than reverting in silence. The declared
value is read out of `rainmod.f90` rather than copied, so the control arm
reproduces the compiled value by construction and the bit-identity the arms are
read against is a property of the route rather than a coincidence to verify
each time.

What remains is the three points themselves, on one binary and one restart,
against the predictions and the three falsifying conditions above.

The Stephens fit itself remains uncited beyond `radmod`'s own header naming
Stephens (1978) and Stephens et al. (1984); neither is in `references/`. Its
validity outside the water paths it was fitted over is untested, and that is
the second half of what `opaque-constants.md` finding 8 named. That is
`world-jimc`.
