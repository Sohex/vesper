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

**The longwave consumer matches, and `acllwr` now carries CCM3's own value.**
`radmod`'s `ztaucc0 = 1 - dcc*(1 - exp(-1.66*acllwr*1000*dql*dp*dsigma/ga))` is
Eq. 12-13 with `D` = 1.66 written as a literal, `acllwr` in the place of
`k_abs`, and the in-cloud water path in g m-2. Eq. 14 makes `k_abs` a
phase-weighted mean and gives `k_l` = 0.090361 m2/g for liquid cloud water,
chosen so that `D k_l` = 0.15 is "in the range of observations and theory";
Eq. 15's `k_i` = 0.005 + 1/r_ei is smaller than `k_l` at every effective radius,
so no weighted `k_abs` can exceed `k_l`. This model has no phase split, so the
liquid-only limit is both the coefficient consistent with the one condensate
`mkclouds` diagnoses and the ceiling on any weighted value.

Upstream's 0.100 sat ABOVE that ceiling, so it could not be a rounded weighted
`k_abs`, and it carried no comment, no unit and no source. It is now 0.090361,
and the difference is confined to thin cloud: the modelled longwave cloud
emissivity moves by at most 0.0373, at a cloud water path of 6.35 g m-2, and by
nothing at all in the two lowest layers of the table below, which are saturated
at either coefficient.

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

## The Stephens fit, read, and what it changed

Both papers are in `references/` and both have been read against the code rather
than against `radmod`'s header. They settle three questions the section above
could only name, and one of them changed the model.

**The band-1 optical depth was band 2's, and now is not (world-jgen).**
Stephens (1978) p. 2124 gives TWO broadband optical depths, one per spectral
region, each fitted by least squares to Mie calculations on the eight standard
cloud models of Part 1:

    log10(tau_N1) = 0.2633 + 1.7095 ln(log10 W)   0.3 to 0.75 um, Eq. (10a)
    log10(tau_N2) = 0.3492 + 1.6518 ln(log10 W)   0.75 to 4.0 um, Eq. (10b)

with W the cloud liquid water path in g m-2. Both are a power law in log10 W:
10^(a + b ln x) is 10^a x^(b ln 10), so (10a) is a prefactor 1.8336 with an
exponent 3.9363 and (10b) is 2.2346 with 3.8034. `radmod` carried one value for
both bands, `2*ALOG10(1.5 + zlwp)**3.9`. That value tracks (10b) to within 3 per
cent at every path from 20 to 10,000 g m-2 and runs 4 to 19 per cent ABOVE (10a)
over the same range. The split it was applied across is Stephens's own: `radmod`
divides its shortwave at 0.75 um, exactly where (10a) gives way to (10b), and
`zrcl1` is Eq. (1), the conservative-scattering solution the paper reserves for
the region below 0.75 um. So band 2 had its own equation and band 1 had band
2's, one-signed, and the visible band's cloud was too bright. Each band now
takes its own fit.

**The `1.5 +` offset is gone, and what replaced it is a physical statement
rather than a guard.** Stephens's form is singular at W = 1 g m-2 and undefined
below it; the offset kept the logarithm's argument above 1 and `radmod`'s own
comment said so. It also returned 0.055 of optical depth in a layer holding no
cloud water at all. Removing it is not enough by itself, because the fitted
range is Figs. 1a and 1b, 10 to 10,000 g m-2, and this model's top three layers
sit one to three decades below that. What happens under 10 g m-2 is an
extrapolation, and the choice of extrapolation is the question the offset was
hiding.

Extending the fitted form is the wrong answer to it. The fit's elasticity
d ln tau / d ln W is b ln10 / (ln10 log10 W), which is 1.7 at 10 g m-2 and grows
without bound as W falls towards 1, so the extrapolated optical depth collapses
far faster than linearly. Read back through Stephens Eq. (7),
tau_N = 1.5 W / r_e, the extrapolated fit at 2.15 g m-2 implies a droplet
effective radius of 134 um and at 1.5 g m-2 over a millimetre. `radmod` now
continues each band linearly in the water path below 10 g m-2, matching that
band's own fitted value there, which is Eq. (7) at the effective radius the fit
itself implies at the bottom of its range: 8.2 um in band 1 and 6.7 um in
band 2. It is continuous at the join, it returns zero optical depth at zero
cloud water, and it introduces no constant the fit does not already carry. It
does not pretend to be the thin-cloud parameterization Stephens, Ackerman and
Smith (1984) p. 690 ask for below an optical depth of about 2, quoted in full in
the next section; it is the weakest continuation that does not assert something
false about the droplets.

Because the offset dominated the inherited optical depth at the thin end, the
correction does not have one sign everywhere. Over the paths the model's cloud
water actually reaches -- the bottom seven layers, 7 to 90 g m-2 -- band-1 cloud
gets less bright, which is a warming. In the top three layers, under 2.5 g m-2,
it gets brighter, because the linear continuation returns more optical depth
than the offset did. The second is much the smaller of the two, and
`exoplasim/analysis/cloud_optical_depth_bracket.json` carries the layer table.

**`tswr1`, `tswr2` and `tswr3` were not the formulas the 1984 paper withdrew,
and they sat much further from the tables than those formulas did. world-f9ig
deleted all three and `swr` now interpolates the 1984 tables.** Stephens et al.
(1984) p. 687 say of the 1978 analytic
representations: "Unfortunately, these formulas contained errors and better
results are obtained when the parameter values are taken directly from the
tables provided." The formulas they mean are Eqs. (11a) and (11b), surface
polynomials in ln(tau_N) and mu0 whose coefficients are Table 1 of the 1978
paper. `radmod` evaluated nothing of that shape. It carried three fits of its
own, one coefficient each:

    beta1  = tswr1 sqrt(mu0)                    tswr1 = 0.077
    beta2  = tswr2 sqrt(mu0) / ln(3 + 0.1 tau)  tswr2 = 0.065
    1-om0  = tswr3 mu0^2 ln(1000/tau)           tswr3 = 0.0055

`radmod`'s own comments called all three tunings of cloud albedo and scattering,
which is what they were. The 1984 withdrawal therefore did not reach them, and
what retired them was the size of their disagreement with the tables rather than
the withdrawal.

What the 1984 paper does change is the reference they should be measured
against. Its Table 1(a) replaces the 1978 single-scattering albedo table
outright, at 99 of its 108 entries, and is monotone in mu0 where the 1978 table
was not: 1 - omega at tau_N = 1 and mu0 = 0.1 goes from 0.0173 to 0.0017, a
factor of ten, and at tau_N = 10 from 0.0028 to 0.0003. Its Tables 1(b) and 1(c)
are the 1978 backscatter tables with two points altered in each, which is what
"this smoothing has altered only a couple of points in each table" (p. 689)
means. `exoplasim/scripts/stephens_tables_vs_fits.py` asserts those counts
rather than trusting the transcription, and asserts the monotonicity too, so a
slipped digit in either table fails the script rather than moving a result.

Measured against the tables at the optical depths this model's cloud layers
reach, the three fits are not approximations to them. `beta2` from the fit is
about half the tabulated value from tau_N 10 to 200, and the fitted co-albedo is
about twice the tabulated one at high sun; both push the same way, so the
tabulated cloud reflects far more in the near infrared than the fitted cloud
does. `beta1` is the closest of the three, within about 20 per cent over the
same range. Interpolating the tables in place of the fits is therefore not a
correction to a port; it is a different shortwave cloud scheme, and its price is
in the section below. It is the scheme the model now runs: `swr` interpolates
Tables 1(a), 1(b) and 1(c) bilinearly, as p. 689 prescribes, with both axes
clamped to the tabulated range. `smoke_test` reads the Fortran tables back
against the transcription in `stephens_tables_vs_fits.py` on every invocation.

### Do this model's cloud water paths sit inside the fitted range

The fitted range is Fig. 1's axes, 10 to 10,000 g m-2, over eight terrestrial
cloud models (St I, St II, Sc I, Sc II, As, Ns, Cu, Cb); the tuned tables span
tau_N 1 to 500, which inverts through (10a) to about 7 to 14,000 g m-2; and
p. 2130 calls 10 to 1000 g m-2 "typical of what may be expected in stratiform
clouds observed in the real atmosphere". `radmod` caps `zlwp` at 1000, which
lands on the top of that range at tau 145, well inside the tuned table. The
CEILING is safe.

The floor is not, and this is answerable without a climatology because the
surface pressure cancels out of it. `mkclouds` builds `dql` as a specific
quantity, `clwref*exp(-z/hl)*gascon*T/(sigma*dp)`, and `swr` turns it into a
path as `1000*dql*dp/ga*dsigma`, so

    zlwp(j) = 1000 clwref gascon T(j) (dsigma(j)/sigma(j)) exp(-z(j)/hl) / ga

carries the sigma grid, the temperature profile and the column water and
nothing else. On the ten-level quartic grid at a 50 hPa top, with the derived
`clwhsc`, the bottom three layers carry 43 to 105 g m-2 across column water
from 10 to 50 kg m-2, which is inside the fitted range; the layers above sigma
0.55 carry under 30 g m-2 and the top three carry under 2.2 g m-2 at every
column water a moist column reaches. So the model's high cloud layers run one
to three decades BELOW the range the fit was made over, in exactly the regime
Stephens, Ackerman and Smith (1984) p. 690 excludes: "the reflection errors
become large when the liquid water path is small (i.e., when cloud optical
depth is less than about 2). Special parameterization of reflection (and for
that matter, absorption) are required for optically thin cloud such as cirrus
cloud."

Two further stated limits apply to cells this model has. Stephens (1978)
p. 2127 on the cloud-over-reflecting-surface correction: "care must be observed
when (12) is applied for larger surface albedos", with the error growing sharply
above a = 0.75, which is every snow and sea-ice cell. And p. 2132: "The
parameterization is less effective at lower solar elevations", which at 32
degrees of obliquity is a larger share of the year than on Earth.

### The absent effective radius is a scheme boundary, not a gap

Stephens does not eliminate the droplet effective radius; he buries a particular
one. His Eq. (7) is tau_N = 1.5 W / r_e, and p. 2125 says of the fits that "the
dependence of tau_N on r_e as described in (7) has been inherently parameterized
in (10a) and (10b)", with p. 2132 adding that a formal parameterization of
optical depth against water path "is academic at this point due to our current
lack of knowledge of the behavior of effective radius". The calculations also
"assume that the cloud is vertically uniform with respect to drop-size
distribution" (p. 2124).

So this model does not lack an input its own shortwave scheme wants. It carries,
silently, the r_e(W) relation of eight terrestrial stratiform cloud models,
which through (10a) implies an effective radius near 5 um at 100 g m-2 and near
11 um at 1000. THAT is what to declare: this world's cloud droplet effective
radius is Earth's stratiform relation, because nothing in this project
determines a cloud condensation nucleus population and nothing measures one.
Filling the gap instead would mean replacing Stephens with a scheme that carries
r_e explicitly, which is Slingo, and that is a scheme decision rather than a
constant.

### `clwref` does not move under Stephens, and that is falsifiable

The question the mismatch invites is what `clwref` should be under the optics
that actually consume it. It has no answer, and the reason is arithmetic rather
than preference.

`clwref` is a condensate density fixed by microphysics. Stephens and Slingo are
two mappings from a water path to an optical depth, and changing which one reads
the water does not change how much water is in the cloud. Re-fitting `clwref` so
that Stephens returns Slingo's optical depth would be a physical constant used
to absorb a difference between two optics schemes.

It also cannot be done. Slingo's optical depth is linear in the path; the
Stephens fit's elasticity runs from about 1.9 at 5 g m-2 to 0.73 at 200. Ask for
the scale factor s on `clwref` that makes `2*log10(1.5 + sW)**3.9` equal
`1.5 W / r_e`, and s is a different number at every path: over 20 to 300 g m-2
it spreads 36 per cent at r_e = 14 um, 63 per cent at 10 um and 86 per cent at
8 um. A single value matches at one path and misses at every other, so there is
no `clwref` that carries CCM3's calibration across the boundary. The bracket
below stands as the answer instead, and it is a bracket on the WATER rather than
on the optics.
### What this is worth, and what would settle it

**world-jgen, the band-1 optical depth. Bracket: +0.9 to +2.1 W m-2 at the top
of the atmosphere, +0.7 to +1.7 K, warming.** Measured on the bootstrap
climatology of run_2b20e3324bb0, 2026-08-25, by
`exoplasim/scripts/cloud_optical_depth_bracket.py`. It is a bracket and not a
measurement, and the reasons are in the script.

The instrument carries the model's own per-layer cloud water path, built from
its temperature, humidity and sigma grid; the model's own total cloud cover with
its own random overlap; the model's zenith-angle geometry integrated over the
day; and the band-1 surface albedo through Stephens Eq. (12). Its own check is
that the band-1 incident flux it reconstructs agrees with the climatology's
`rst + rsut` to better than a part in ten thousand. What it cannot carry is the
per-layer cloud fraction, which is not written out at any cadence and cannot be
recovered from a time mean of `hur` -- the diagnosis is the square of a
threshold difference, so the mean of the diagnosis is not the diagnosis of the
mean, and on this field it leaves 72 per cent of columns cloudless against a
written cover of 0.61. The bracket's width is therefore the spread over where in
the column the model's own cover is placed. It still leans HIGH, because band-1
gas absorption and Rayleigh scattering are not in the chain and neither is
multiple reflection between layers.

That replaces an earlier hand chain of 0.6 to 3.4 W m-2 and 0.5 to 3.1 K, which
had no cloud fraction, no overlap and no surface in it.

**world-f9ig, the tuned coefficients against the tables. Bracket: -6.5 to
-25.4 W m-2, -5.3 to -20.6 K, cooling**, on the same cover and insolation and
over a black surface, so it is a magnitude upper bound. Measured 2026-08-25 by
`exoplasim/scripts/stephens_tables_vs_fits.py`, with `cloud_absorption_scale`
applied to both sides, which is what the model does. It is an order of magnitude
larger than anything else in this document, and that is why the tables were
adopted rather than the coefficients kept: three numbers whose provenance is a
tuning were standing between the model and a published quantity, and the size of
the gap is information about the tuning rather than a reason to preserve it. The
whole shortwave balance moves with them, so the design flux is re-derived on the
next commissioning cycle. `cloud_absorption_scale` was NOT stranded: it is a
star-over-Sun ratio of range-2 flux-weighted co-albedo, the tables carry the
same solar weighting in the denominator, and it now multiplies the interpolated
co-albedo through the namelist key `cloudabs`.

WHAT SETTLES world-jgen IS A RUN, and the run is nameable: a T21 pair on a
settled baseline, control against an arm carrying the correction as it now
stands -- band 1 on Eq. (10a), band 2 on Eq. (10b), each continued linearly in
the water path below 10 g m-2. `tswr1` is HELD at its inherited value in that
arm, on purpose: it was tuned against the optical depth being corrected, and an
arm that moves both reports nothing about either. If the corrected optical depth
at the inherited tuning makes an agreement worse, that is information about
`tswr1`. world-f9ig settled that by deleting it: the band-1 backscatter is now
interpolated from Table 1(b) at the layer's own range-1 optical depth, so the
arm's control and its arm no longer share a coefficient that was tuned against
the quantity being corrected.

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

## Measured: the three arms, 2026-08-26

T21, `most_plasim_t21_l10_p8.x` at model source `a041a1e9`, all three branched
from `run_14906cb7b914`'s `MOST_REST.00034`, 25 orbits each, each pinned to its
own eight cores, host load 25 to 50 over 32 cores. The energy fixer and the
energy diagnostics are OFF in all three, which is a departure from
`config/planet.yaml` and is forced: `nenergy > 0` segfaults at every orbit
boundary on this source, per `exoplasim/notes/energy-diagnostic-restart-defect.md`.
It is common to all three arms and the fixer is worth about 0.26 K, against
separations of four kelvin.

| arm | run | `d` mean surface T, K | TOA forcing, first 3 orbits | `d` sea ice |
| --- | --- | ---: | ---: | ---: |
| 0.5x, `clwref` 0.000105 | `run_7536d6b7288a` | **+4.34** | +4.04 W/m2 | -0.040 |
| 1.0x, `clwref` 0.00021 | `run_57a43e1fc3f4` | control | control | control |
| 2.0x, `clwref` 0.00042 | `run_ccb8ee8384e6` | **-4.64** | -3.42 W/m2 | +0.058 |

Both separations plateau by about orbit 15 and oscillate after it, the 0.5x arm
around +4.2 and the 2.0x arm around -4.75.

**The three falsifying conditions, and all three hold.**

- **The 1.0x arm does not differ from a run that declares nothing.** A run with
  `cloud_water_reference_kg_m3` deleted from its config, `run_0730a12ecfbd`,
  stages `CLWREF = 0.00021` from `rainmod_default("clwref")` and writes a
  `MOST_REST.00000` bit-identical to the control's, sha `77c3a916`. This is the
  check that could have failed and it is the route's property, not a
  coincidence.
- **Neither arm moves the simulated global mean by less than 1 K.** Both move
  more than four.
- **The 2.0x arm cools and the 0.5x arm warms.** The sign this note assigned
  from the low cloud's saturated longwave emissivity is the sign measured.

**The predicted magnitude half holds and half does not.** "Both arms move the
simulated global mean by more than 2 K" is confirmed at better than twice the
threshold. "The plausible range is 5 to 10 K" is missed low: the measured
separations are 4.34 and 4.64, so the range's lower end is 0.4 to 0.7 K above
what the arms did. The bound was widened upward because the offline estimate
that priced PHYS-11 had been falsified in the direction of the model responding
MORE; this time the model responded slightly LESS, so the widening was in the
wrong direction. That is the finding and the range is left as it was written.

**The response is asymmetric, and the reflectance table does not predict that
half.** Doubling `clwref` cools by more than halving it warms, 4.75 against 4.2,
where the band-1 reflectance moves +0.09 on the 2.0x arm and -0.13 on the 0.5x
arm and would give the opposite asymmetry. What is available to explain it is the
sea ice, which the 2.0x arm grows by +0.058 and the 0.5x arm melts by -0.040, and
the ice-albedo feedback is stronger on the cooling side. That is an argument and
not a measurement; separating it needs the arms repeated with the sea ice held.

**What the arms do NOT establish.** They are not converged. Every arm, control
included, is inside a common relaxation of at least 12.5 K from the donor's state,
because that state was integrated by the superseded source and config;
`exoplasim/notes/forcing-bundle-predictions.md` measures that drift. The
separations above are paired and share it, so they are attributable to `clwref`;
the absolute temperatures are not this world's.

**The conclusion the bracket was for.** `clwref` moves the modelled global mean
by about nine kelvin across a factor of four, which is larger than any other term
in the forcing bundle and larger than the bundle's own registered sum. It is a
first-order lever on this model's radiation and a value with one source, no
range and no sensitivity behind it.

The Stephens fit is read and the section above is what it settled, so the
second half of `opaque-constants.md` finding 8 is answered: this model's low
cloud layers sit inside the range the fit was made over and its high ones sit
one to three decades below it, in the regime the 1984 revision excludes by name.

One thing the papers cannot settle remains, and it needs a run rather than a
reading: what the two changes together are worth to the modelled balance. Both
are in the model now, so the pair that settles it is one T21 pair on the current
source, not two. The scheme question is closed. Stephens's tuned tables replaced
`tswr1`, `tswr2` and `tswr3` because those three had no derivation anywhere and
sat a factor of two from published values in the near-infrared backscatter at
the optical depths this model reaches; a sourced quantity that disagrees with a
tuned one is the reason to take the sourced one.
