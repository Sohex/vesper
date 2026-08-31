# The radiation used a different particle radius than the transport

Read 2026-08-18, against `vendor/exoplasim` and Cohen et al. (2024),
`references/pdf/2307.10931v2.pdf`.

**Fixed in this fork.** `aero_ini` use-associates radmod's copy under an alias,
`use radmod, only: l_aerorad, aerofile, rad_apart => apart`
(`aeromod.f90:192`), and assigns `rad_apart = apart` after the namelist read
(`:277`), with the defect recorded at `:266-276`. The two copies cannot diverge.
Everything below is what the divergence was and what it would have cost, which
is still the argument for offering the fix upstream.

## The mechanism

`apart`, the aerosol particle radius, is declared twice: in `aeromod` at 50e-9
and in `radmod` at 50e-09. Only aeromod's is in a namelist. `aero_ini` read
`aero_nl` into aeromod's copy and use-associated `l_aerorad` and `aerofile` from
radmod, but NOT `apart`, so radmod's copy kept its compiled default however the
run was configured.

Both halves of the optical depth then disagree:

- `aerocore` calls `mmr2n(mmr, apart, rhop, ...)` with AEROMOD's value, so the
  number density `nrho` is right for the particle the run asked for.
- `radmod`'s `aeroprof` computes
  `daerod = nrho * PI * apart**2 * qex1 * zdz` with RADMOD's value.

For fixed mass, `nrho` goes as `1/r^3` and the optical depth as `nrho * r^2`, so
a correct calculation goes as `1/r`. With the two copies disagreeing the model
instead computes `nrho(r_true) * r_radmod^2`, and

    tau_computed / tau_intended = (r_radmod / r_true)^2 = (50e-9 / r_true)^2

Both defaults are 50 nm, so the defect was dormant until a run set `apart`.
This project would have set it: `world-906` found that `enable_dust_emission`
wrote the fourteen `DUST*` keys and neither `APART` nor `RHOP`, and made it read
both out of the aerofile's own sidecar, where `dust_aerofile.py` had derived
`APART = 2.2068e-06` for this planet's gravity. Through the divergence that
would have been `(50e-9/2.2068e-6)**2`, and through Stokes settling with the
grain left at 50 nm it was worth about 1790 times too slow a fall speed.
`notes/audits/model-earth-centrism.md` has the settling half.

## Why this is not a theoretical concern

Cohen et al. (2024) is, per ExoPlaSim's maintainer, the only published science
using this module. Its haze is radiatively active, not a passive tracer, and its
Appendix eq. A10 states the haze optical depth as

    tau_N = N * Qext * pi * r^2 * dz,  "r is the particle radius in meters"

which is `aeroprof` line for line. The paper uses ONE particle size for all 64
simulations: **500 nm**, with `Qext` computed offline for that size from He et
al. (2023) refractive indices and reported per band in its Table 1.

So on stock 3.4.2 there are two possibilities and both are large:

| what the run did | consequence |
| --- | --- |
| set `apart = 500e-9` in `aero_nl` | `nrho` right, cross-section uses 50 nm: **tau is 1% of intended**, and `Qext` is for a particle the geometry does not describe |
| never set `apart` | everything is 50 nm, which is not the 500 nm the paper documents, and tau comes out ~10x the intended value |

The first is far more likely, because the gravitational settling scheme needs
the real radius and the paper is explicit that particle density is "a necessary
input" to it.

## What is NOT established

Which of the two the published runs did, and whether they ran stock 3.4.2 at
all. Maureen Cohen authored the module and may carry a local fix that never
reached the release. Nothing here is a claim about the paper's results; it is a
claim about what stock ExoPlaSim computes, and a question about which code the
runs used. That question is answerable in one line of their namelist.

## What the cited settling papers do and do not supply

Cohen et al. describe the settling as "based on Steinrueck et al. (2021) and
Parmentier et al. (2013)", and both are now in `references/`. Read against the
source they narrow the question rather than widening it.

**They supply the settling velocity and nothing else.** Parmentier's Appendix A
equates gravity and drag, `Vf^2 CD = (8a/3)(rho_p - rho)/rho`, with tabulated
drag coefficients running from the Stokes limit `CD = 24` to the asymptotic
`CD = 0.45`. Neither paper carries a surface-deposition treatment, because both
model gas giants that have no surface.

**Steinrueck's sink is a rate; ExoPlaSim's is not.** Its lower boundary
condition is a deep sink, `L = -chi/tau_loss` for `p > p_deep`, with
`tau_loss = 1e3 s` and `p_deep = 100 mbar`, standing in for cloud nucleation and
thermal ablation. The timescale is explicit and independent of the timestep, and
it acts over a deep region. ExoPlaSim's `mmr = mmr*0.01` acts in ONE layer, per
STEP, so its implied timescale is `dt/ln(100)`. At the 15-minute timestep of
Cohen et al. that is about 195 s.

So the bottom-layer sink is ExoPlaSim's own and is not inherited from either
citation. That matters for how the fix is offered: it is not contradicting a
published design, it is replacing an undocumented implementation detail with the
treatment a rocky planet's lower boundary actually calls for.

## The settling is separately affected, by defects already in PR #58

Stokes settling goes as `1/mu`. `aerocore`'s viscosity used `(4/25)` as INTEGER
division, which evaluates to 0 and deleted the collision integral's temperature
dependence, leaving `mu` 11 to 18 per cent low over 200-320 K in N2 -- so the
settling velocity is 12 to 22 per cent too fast. `mmr2n` used `(4/3)` the same
way, evaluating to 1, so the sphere volume was 4/3 too small and the number
density 33 per cent too high, which propagates straight into optical depth.

Both are in the aerocore-defects branch already offered as PR #58. Combined with
the particle radius above, a configuration using 500 nm particles carries an
optical depth of roughly `(50/500)^2 * 1.33`, about 1.3 per cent of intent, and
PR #58 alone corrects only the 1.33.

## The paper's declared limitation cuts both ways

Cohen et al. state their haze scheme "does not account for sinks (chemical
sinks, wet deposition) aside from gravitational settling onto the surface".

That REMOVES wet deposition from the erratum question. Its absence is declared,
acknowledged and scoped, so the wet-scavenging half of the deposition work is a
capability offer and not a defect report against that paper.

It SHARPENS the other half. The sentence asserts that gravitational settling
onto the surface is operative, and in stock ExoPlaSim it is not: the bottom
layer's update ADDS the settling flux `gz(nl)` where it must subtract it, so
nothing ever leaves the atmosphere by sedimentation and the bottom layer builds
up without bound. The 99%-per-step line is what holds it down. So the only sink
the paper names did not run, and the sink that did run is not mentioned, because
it is an implementation detail rather than a modelling choice.

## The two fixes are coupled and must not be split

The sign fix is in `aerocore-defects`, offered as PR #58. The 99%-per-step line
is in the deposition branch, not yet offered. Their interaction is the thing to
say out loud:

| state | what removes haze at the surface |
| --- | --- |
| stock | the 99%-per-step scrub only; sedimentation is sign-inverted and removes nothing |
| PR #58 alone | sedimentation AND the scrub -- two sinks, more removal than before |
| #58 plus deposition | sedimentation, with the scrub replaceable by a timestep-independent velocity |

Merging #58 on its own therefore does not simply restore the intended scheme; it
leaves a band-aid in place over a wound that has been closed. Whether the
`ldepvel = 0` default should keep that line at all, once the sign is fixed, is a
judgement for the maintainer rather than something to decide here -- keeping it
preserves existing results bit for bit, dropping it makes the model match what
the paper describes.

## The paper's own appendix confirms two of the fixes

Cohen et al. lay their scheme out in Appendix A, and it can be read against the
source line for line. Two of the defects are contradictions between that
appendix and the code that implements it.

**Eq. A10, the optical depth.** `tau_N = N * Qext * pi * r^2 * dz`, "r is the
particle radius in meters". That is `aeroprof` exactly, except that the code's
`r` is radmod's `apart`, which the namelist never reaches. Covered above.

**Eq. A5, the viscosity.** The Rosner (1986) parameterisation, with the exponent
written explicitly:

    eta = 5*sqrt(pi*m*kb*T) * (kb*T/eps)**0.16 / (16 * 1.22 * pi * d^2)

for N2 with m = 4.652e-26 kg, d = 3.64e-10 m, eps/kb = 95.5 K. The code writes
that exponent as `(4/25)`, which in Fortran is INTEGER division and evaluates to
0, so `(temp/eps)**0 = 1` and the temperature dependence vanishes. Every other
factor matches A5 -- the 5/16, the 1/1.22, the 1/(pi d^2), the sqrt(pi m kb T)
-- so it is only the exponent that is lost, which is why it survives a reading.

Using their own eps, and settling going as 1/eta:

| T (K) | (T/95.5)^0.16 | viscosity low by | settling fast by |
| ---: | ---: | ---: | ---: |
| 200 | 1.126 | 11.2% | 12.6% |
| 250 | 1.167 | 14.3% | 16.6% |
| 300 | 1.201 | 16.7% | 20.1% |

So for temperate rocky planets the particles fall about 17% too fast. It is a
one-signed bias rather than scatter, so it shifts the equilibrium burden
coherently. The fix is in PR #58.

## The paper names the limitation the shortwave weights remove

On its own model limitations, Cohen et al. write that "the absorptivity of water
vapour in shortwave band 2 (lambda > 0.75 um) relies on a parameterisation that
is tuned to the solar spectrum and is likely less accurate for M-class stellar
spectra".

That is precisely what `h2osww` exists for, and it means the shortwave weights
offered as PR #62 are not a speculative capability: a published user has already
identified the defect and carried it as a known inaccuracy. The band-2
absorptance is a fit expressed as a fraction of total incident SOLAR flux, so it
carries the Sun's share of flux in that band; an M dwarf puts far more of its
output beyond 0.75 um and the share is different.

The same construction underlies the Lacis & Hansen ozone terms, which is why
`o3visw` and `o3uvw` exist alongside it.

## Does any of this port to an aquaplanet

Cohen et al. model aquaplanets; this project models a world with land, deserts
and lithology-dependent dust. The split is clean.

**Ports unchanged**, because it is mechanism rather than surface physics: the
particle radius handover, the settling sign, the viscosity and number-density
integer divisions, and the longwave term.

**Ports as a capability, not as a number**: the deposition velocity and the wet
scavenging coefficients. A deposition velocity to open ocean is not the one to
vegetated land or to desert pavement, and the scavenging coefficients are
aerosol- and rain-specific. This is exactly why the patch refuses to run without
explicit values rather than carrying a default, and that refusal is what makes
it safe to offer across the boundary.

**Does not port at all**: the dust emission scheme, which is a land-surface
parameterisation and is meaningless over ocean. It is being held anyway.

## The cheap diagnostic

Whether the bottom sink matters for a given result depends on whether it is
rate-limiting for the column burden, which is not something this repository can
determine. It is cheaply testable by whoever owns the runs: **re-run at a
different timestep.** The sink's implied timescale scales with `dt` and nothing
else physical does, so if the haze burden moves when only the timestep changes,
the sink is both rate-limiting and timestep-dependent. If it does not move, this
defect does not reach that result.

## Why it matters now

ExoPlaSim's maintainer has said that the aerocore defects reported in upstream
PR #58 may affect published results, and that Cohen will re-run her published
grid with the fixes on a university cluster in autumn 2026, with an erratum if
the differences are large enough to matter.

If that re-run carries #58's fixes but not this one, it re-runs with the
radiation still seeing the wrong particle size, and a conclusion of "the
differences are small" would be drawn from an incomplete fix set. This one is
thirty-one lines and is the largest single error in the aerosol path.

The same argument applies, less sharply, to the two other aerosol changes not
yet offered upstream: the bottom-level sink that removes 99% of the layer per
timestep, which is a rate that depends on the integration step and therefore is
not a rate, and the absence of any longwave aerosol term, which gives a model
that can cool with haze and cannot warm with it.

**The longwave half is no longer absent in this fork.** `aeroqlw`
(`radmod.f90:324`) is the thermal-IR absorption optical depth per unit band-1
extinction optical depth for the INTERACTIVE aerosol, beside `dustqlw` for the
prescribed field; the two are kept apart because the two paths carry different
particles, and there is no defensible default for either, so `radini` ABORTS on
an interactive aerosol with `aeroqlw` at zero (`radmod.f90:1382`). world-24v gave
it a writer and turned the `l_aerorad` pin into a switch. The exclusion that made
the two paths mutually exclusive is also gone: they no longer share one optical
depth slot, so a prescribed column and a transported one can both act.
