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
