# The radiation used a different particle radius than the transport

Read 2026-08-18, against `vendor/exoplasim` and Cohen et al. (2024),
`references/2307.10931v2.pdf`.

## The mechanism

`apart`, the aerosol particle radius, is declared twice: in `aeromod` at 50e-9
and in `radmod` at 50e-09. Only aeromod's is in a namelist. `aero_ini` reads
`aero_nl` into aeromod's copy and use-associates `l_aerorad` and `aerofile` from
radmod, but NOT `apart`, so radmod's copy keeps its compiled default however the
run is configured.

Both halves of the optical depth then disagree:

- `aerocore` calls `mmr2n(mmr, apart, rhop, ...)` with AEROMOD's value, so the
  number density `nrho` is right for the particle the run asked for.
- `radmod`'s `aeroprof` computes
  `daerod = nrho * PI * apart**2 * qex1 * zdz` with RADMOD's value.

For fixed mass, `nrho` goes as `1/r^3` and the optical depth as `nrho * r^2`, so
a correct calculation goes as `1/r`. With the two copies disagreeing the model
instead computes `nrho(r_true) * r_radmod^2`, and

    tau_computed / tau_intended = (r_radmod / r_true)^2 = (50e-9 / r_true)^2

Both defaults are 50 nm, so the defect is dormant until a run sets `apart`.

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
