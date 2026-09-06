#!/usr/bin/env python
"""What a droplet number would be worth, through the model's own cloud optics.

    python analysis/aerosol_indirect_sensitivity.py

WORLDBUILDING FRAME. Vesper is a simulated super-Earth. Everything below is a
property of this project's climate model or of its offline aerosol products.
Nothing here is about the real world.

CLIM-64 asks what the aerosol indirect effect costs. Half of that answer is what
the term is WORTH, and this computes it. The other half -- what it costs to
build -- is in `notes/audits/aerosol-indirect-effect-cost.md`, which also carries
the criterion this number is judged against and the date it was fixed.

## What is computed, and why it is a sensitivity rather than a forcing

A forcing needs a droplet number, and this project has none: the aerosol products
are offline and mass-based, and the model carries one monodisperse radius per
species slot. What CAN be computed without inventing one is the derivative --
what the model's own cloud albedo does per factor of two in droplet number -- and
then a bracket on the factor the aerosol side could supply.

The chain is the model's, not a textbook's, at every step:

  1. `rainmod`'s CCM3 cloud water, `dql = clwref exp(-z/zzh) R T / (sigma ps)`
     with `zzh = clwhsc ln(1 + prw)`, reconstructed from the climatology's own
     `prw`, `ta` and `ps`.
  2. `rainmod`'s stratiform cloud fraction, `((rh - rcrit)/(1 - rcrit))^2`, on
     the climatology's own `hur`, with `rcrit` scaled off the T21 anchor exactly
     as `rainini` scales it.
  3. `radmod`'s liquid water path and its two Stephens (1978) optical depth fits.
  4. `radmod`'s cloud reflectance: the Stephens, Ackerman and Smith (1984)
     backscatter and co-albedo TABLES, bilinearly interpolated, read from
     `exoplasim/scripts/stephens_tables_vs_fits.py`, which already checks that
     transcription against the model source entry for entry.

## The one identity that makes this cheap, and it is the model's own

At fixed liquid water path a cloud's optical depth goes as the cube root of
droplet number, because Stephens Eq. (7) is `tau = 1.5 W / (rho_w r_e)` and
`r_e` goes as `(LWC/N)^(1/3)`. So a doubling of N multiplies every optical depth
by `2^(1/3)`, and nothing else in the chain moves. The reflectance is then
re-evaluated on the SAME tables at the new optical depth, so no linearisation is
taken and the answer carries the tables' own curvature.

`radmod`'s conservative branch is `A = 1 - 1/(1 + beta tau / mu0)`, whose
logarithmic derivative is `A(1 - A)`. That is Twomey's sensitivity exactly, and
its appearance here is a check rather than a coincidence: the two-stream form the
model happens to use and the form Twomey wrote it in are the same function. The
report carries both the table answer and that closed form so the two can be
compared, and a disagreement between them is a coding error rather than a
physical claim.

## Cloud fraction is the model's own total, distributed by an identity

The model emits the column total `clt` and no per-level cover, and the per-level
cover cannot be reconstructed from the climatology: `mkclouds` adds a convective
branch keyed on `icclev`, which is not an output, to the stratiform branch that
`hur` does support. Reconstructing the stratiform half alone was tried and the
check refused it -- the median ratio of reconstructed total to emitted `clt` came
out at 0.0 over 24,575 cells with cloud, so the missing branch is most of the
cloud rather than a correction to it.

What is used instead needs no reconstruction. `clt` is distributed over the
levels as

    dcc_k = 1 - (1 - clt)^s_k,    sum_k s_k = 1

whose random-overlap total is `1 - (1-clt)^(sum s_k) = clt` EXACTLY, for any
shape `s`. So the emitted total is reproduced by construction and the vertical
shape is the only assumption left. Two shapes are run and reported as a bracket:
the reconstructed liquid water profile, which puts cloud where the model's own
CCM3 diagnostic puts water, and a uniform split, which puts it everywhere. The
answer is quoted across both.

## What this deliberately does not do

It does not compute a droplet number, from this project's aerosol or from
anything else. The factor of two is a unit of the derivative and not a
prediction, and the bracket on how many factors of two the aerosol side could
supply is argued in the note rather than computed here, because the number
concentration those products would have to carry does not exist yet.

It also stops at the SHORTWAVE. The longwave cloud absorption is grey in liquid
water path with no effective radius at all (`acllwr`), so a droplet number does
not reach it without a second derivation, and pricing that is part of the note's
cost side rather than of this measurement.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import xarray as xr
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[0] / ".." / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "exoplasim" / "scripts"))

ROOT = Path(__file__).resolve().parents[1]

import paths  # noqa: E402
import sensitivity  # noqa: E402
import gridding  # noqa: E402
import stellar  # noqa: E402

# radmod.f90:2667-2670, the two Stephens (1978) Eq. (10a)/(10b) fits.
ZTAUA1, ZTAUP1 = 1.8336, 3.9363
ZTAUA2, ZTAUP2 = 2.2346, 3.8034
ZWFIT = 10.0            # g/m2, the bottom of the fitted range
# rainmod.f90:212-214
CLWREF = 0.00021        # kg/m3
CLWHSC_EARTH = 700.0
RCRITREF = 0.85         # rainmod.f90, the T21 floor
RCNLATREF = 32
EPSC = math.sqrt(np.finfo(np.float32).eps)   # radmod's conservative floor


def stephens():
    """The 1984 tables and radmod's two reflectance branches, from the checker."""
    import stephens_tables_vs_fits as s
    return s


def rcrit_profile(sigma: np.ndarray, nlat: int) -> np.ndarray:
    """`rainini`'s critical relative humidity, level by level.

    `rcrit(:) = MAX(RCRITREF, MAX(sigma, 1-sigma))`, then narrowed by the
    subgrid humidity width factor `(RCNLATREF/NLAT)^(1/3)`. The unit factor is
    skipped rather than applied, as the model skips it, so a rung at the anchor
    is bit-identical.
    """
    limb = np.maximum(RCRITREF, np.maximum(sigma, 1.0 - sigma))
    width = (RCNLATREF / float(nlat)) ** (1.0 / 3.0)
    if width == 1.0:
        return limb
    return 1.0 - (1.0 - limb) * width


def reconstruct(ds: xr.Dataset, cfg: dict) -> dict:
    model = cfg["model"]
    ga = float(cfg["planet"]["gravity_m_s2"])
    gascon = float(model.get("gas_constant_j_kg_k", 287.0))
    nlat = int(model["latitudes"])
    clwhsc = CLWHSC_EARTH * (gascon / ga) / (287.0 / 9.80665)

    sigma = ds["lev"].values.astype(float)
    nlev = sigma.size
    # Half levels either side of each full level, so dsigma is the model's.
    half = np.empty(nlev + 1)
    half[0] = 0.0
    half[1:-1] = 0.5 * (sigma[:-1] + sigma[1:])
    half[-1] = 1.0
    dsigma = np.diff(half)

    ps = ds["ps"].values * 100.0            # hPa -> Pa
    ta = ds["ta"].values                    # (t, lev, lat, lon)
    hur = ds["hur"].values / 100.0
    prw = ds["prw"].values                  # kg/m2
    clt = ds["clt"].values
    czen = np.clip(ds["czen"].values, 0.0, 1.0)

    # Level height by the hypsometric integral on the model's own sigma levels,
    # which is what `zzf` is in the model.
    z = np.zeros_like(ta)
    p_half = half[None, :, None, None] * ps[:, None, :, :]
    p_full = sigma[None, :, None, None] * ps[:, None, :, :]
    zz = np.zeros_like(p_half)
    for k in range(nlev - 1, -1, -1):
        lower = np.where(p_half[:, k + 1] > 0.0, p_half[:, k + 1], p_full[:, k])
        upper = np.maximum(p_half[:, k], 1.0)
        zz[:, k] = zz[:, k + 1] + gascon * ta[:, k] / ga * np.log(lower / upper)
        z[:, k] = 0.5 * (zz[:, k] + zz[:, k + 1])

    zzh = clwhsc * np.log(1.0 + prw)                      # (t, lat, lon)
    dql = (CLWREF * np.exp(np.clip(-z / np.maximum(zzh[:, None], 1e-30), -80, 80))
           * gascon * ta / (sigma[None, :, None, None] * ps[:, None, :, :]))
    dql = np.maximum(dql, 1e-9)

    # The stratiform-only reconstruction, kept ONLY as the refused control that
    # established why the identity below is used instead.
    rcrit = rcrit_profile(sigma, nlat)
    strat = np.minimum(1.0, np.maximum(
        0.0, (hur - rcrit[None, :, None, None])
        / (1.0 - rcrit[None, :, None, None])) ** 2)
    strat_total = 1.0 - np.prod(1.0 - strat, axis=1)

    return dict(sigma=sigma, dsigma=dsigma, ps=ps, dql=dql,
                strat_total=strat_total, clt=clt, czen=czen, ga=ga, nlev=nlev)


def cloud_reflectance(s, tau1, tau2, mu0, cloud_absorption_scale):
    """`swr`'s two cloud reflectances and transmittances, on its own tables.

    Band 1 is conservative so its transmittance is `1 - r`. Band 2 absorbs, and
    `reflect_absorbing` returns the absorptance beside the reflectance, so its
    transmittance is `1 - r - a`.
    """
    b1 = s.interp_table(s.BETA1_1984, tau1, mu0)
    b2 = s.interp_table(s.BETA2_1984, tau2, mu0)
    coalb = np.maximum(EPSC, cloud_absorption_scale
                       * s.interp_table(s.COALB_1984, tau2, mu0))
    # The tables clip mu0 to their own last column at 0.1, so the reflectance
    # is floored at the same place rather than dividing by a zero the table
    # never saw. Unlit columns are masked out of the answer by the caller.
    m = np.maximum(mu0, 0.1)
    r1 = s.reflect_conservative(b1, tau1, m)
    r2, a2 = s.reflect_absorbing(b2, coalb, tau2, m)
    return r1, r2, 1.0 - r1, np.maximum(0.0, 1.0 - r2 - a2)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--climatology", type=Path, default=None,
                    help="default: lib/paths.py's baseline climatology, which "
                         "config/planet.yaml names and which has no fallback")
    ap.add_argument("--factors", type=float, nargs="+", default=[2.0, 4.0],
                    help="droplet-number multipliers to evaluate")
    ap.add_argument("--out", type=Path,
                    default=ROOT / "analysis/aerosol_indirect_sensitivity.json")
    args = ap.parse_args()

    cfg = yaml.safe_load((ROOT / "config/planet.yaml").read_text())
    clim = args.climatology or paths.climatology_path()
    ds = xr.open_dataset(clim)
    s = stephens()
    r = reconstruct(ds, cfg)

    # The refused control, reported as the reason the identity is used.
    ok = r["clt"] > 0.01
    ratio = np.where(ok, r["strat_total"] / np.maximum(r["clt"], 1e-6), np.nan)
    check = dict(
        what="stratiform-only cloud fraction reconstruction, REFUSED",
        reconstructed_over_emitted_median=float(np.nanmedian(ratio)),
        reconstructed_over_emitted_p90=float(np.nanpercentile(ratio, 90)),
        cells_compared=int(ok.sum()),
        why_refused="mkclouds adds a convective branch keyed on icclev, which "
                    "is not an output; the stratiform half alone is not the "
                    "model's cloud",
    )

    lwp = np.minimum(1000.0, 1000.0 * r["dql"] * r["ps"][:, None]
                     / r["ga"] * r["dsigma"][None, :, None, None])
    mu0 = np.broadcast_to(r["czen"][:, None], lwp.shape)
    lit = r["czen"] > 1e-6

    band1 = stellar.band1_fraction()
    cas = float(cfg["model"].get("cloud_absorption_scale", 1.0))
    alb1 = ds["alb1"].values
    alb2 = ds["alb2"].values

    # The area weight through the one door. `cos(lat)` is the metric factor of
    # an equally spaced band and a Gaussian row is not one, so the axis is
    # handed over and checked rather than reconstructed here.
    weights = gridding.gaussian_area_weights(ds["lat"].values, ds.sizes["lon"])

    shapes = {}
    tot = lwp.sum(axis=1, keepdims=True)
    shapes["liquid water profile"] = np.where(tot > 0, lwp / np.maximum(tot, 1e-30),
                                              1.0 / r["nlev"])
    shapes["uniform over levels"] = np.full_like(lwp, 1.0 / r["nlev"])

    def taus(factor):
        f = factor ** (1.0 / 3.0)
        zwl = np.log10(np.maximum(ZWFIT, np.maximum(0.0, lwp)))
        ramp = np.minimum(1.0, np.maximum(0.0, lwp) / ZWFIT)
        return (f * ZTAUA1 * zwl ** ZTAUP1 * ramp,
                f * ZTAUA2 * zwl ** ZTAUP2 * ramp)

    def column_albedo(factor, cc, with_surface):
        """The cloud column's albedo, per band, optionally over the surface.

        Random overlap in reflectance: one minus the product of what each layer
        fails to reflect, weighted by its own cover. That is the same overlap
        rule `mkclouds` uses for the cover itself.

        WITH THE SURFACE UNDERNEATH is the arm that matters for the size of the
        answer rather than for its existence. A cloud over a bright surface adds
        less at the top of the atmosphere than the same cloud over a dark one,
        because the surface was already reflecting what the cloud now reflects.
        The adding is the standard one, `R + T^2 A_s / (1 - R A_s)`, taken on the
        model's own two-band surface albedos.
        """
        t1, t2 = taus(factor)
        r1, r2, tr1, tr2 = cloud_reflectance(s, t1, t2, mu0, cas)
        out = []
        for rb, tb, asurf in ((r1, tr1, alb1), (r2, tr2, alb2)):
            rcol = 1.0 - np.prod(1.0 - cc * rb, axis=1)
            if not with_surface:
                out.append(rcol)
                continue
            tcol = np.prod(1.0 - cc * (1.0 - tb), axis=1)
            out.append(rcol + tcol * tcol * asurf
                       / np.maximum(1e-6, 1.0 - rcol * asurf))
        return band1 * out[0] + (1.0 - band1) * out[1]

    alpha, _ = sensitivity.planetary_albedo()
    results = []
    for shape_name, shp in shapes.items():
      for with_surface in (False, True):
        cc = 1.0 - (1.0 - r["clt"][:, None]) ** shp
        # The identity: the random-overlap total of cc is clt for any shape.
        overlap_err = float(np.nanmax(np.abs(
            (1.0 - np.prod(1.0 - cc, axis=1)) - r["clt"])))
        base = column_albedo(1.0, cc, with_surface)
        for factor in args.factors:
            d_alb = column_albedo(factor, cc, with_surface) - base
            # Insolation weighting: an albedo change where the star is not up is
            # worth nothing, and mu0 is the model's own.
            num = (np.where(lit, d_alb * r["czen"], 0.0)).mean(axis=0)
            den = (np.where(lit, r["czen"], 0.0)).mean(axis=0)
            d_toa = float((weights * num).sum() / max((weights * den).sum(), 1e-12))
            # The one conversion: an albedo change IS a forcing, through the same
            # denominator, and lib/sensitivity.py owns both directions.
            d_k = sensitivity.planetary_albedo_to_kelvin(d_toa, alpha)
            w_m2 = -d_k / sensitivity.kelvin_per_w_m2(alpha)
            results.append(dict(
                vertical_shape=shape_name,
                surface_underneath=bool(with_surface),
                overlap_identity_max_error=overlap_err,
                droplet_number_factor=factor,
                optical_depth_factor=factor ** (1.0 / 3.0),
                insolation_weighted_planetary_albedo_change=d_toa,
                global_mean_shortwave_forcing_w_m2=float(w_m2),
                kelvin=float(d_k),
            ))

    # The closed form, elementwise, as a check on the table path rather than as
    # an answer. radmod's band-1 branch is A = x/(1+x), whose derivative in
    # ln(tau) is A(1-A); a doubling of N is ln(2)/3 in ln(tau). The tables carry
    # curvature the closed form does not, so the two agree to that curvature and
    # a larger gap is a coding error.
    t1a, t2a = taus(1.0)
    t1b, _ = taus(2.0)
    r1a = cloud_reflectance(s, t1a, t2a, mu0, cas)[0]
    r1b = cloud_reflectance(s, t1b, t2a, mu0, cas)[0]
    sel = (t1a > 1.0) & (mu0 > 0.05)
    closed = (r1a * (1.0 - r1a) * math.log(2.0) / 3.0)[sel]
    tabled = (r1b - r1a)[sel]
    rel = np.abs(tabled - closed) / np.maximum(np.abs(closed), 1e-12)

    # The inversion, which is what a reader actually needs: how big a droplet
    # number change the threshold corresponds to. The forcing is close to linear
    # in ln(N) over this range, so the crossing is read off the doubling.
    per_doubling = [x["global_mean_shortwave_forcing_w_m2"]
                    for x in results if x["droplet_number_factor"] == 2.0]
    crossing = sorted(2.0 ** (1.5 / abs(v)) for v in per_doubling if v)

    report = dict(
        measured_on="2026-09-05",
        climatology=str(clim),
        run_id=str(ds.attrs.get("vesper_run_id", "")),
        source_build=str(ds.attrs.get("vesper_source_build", "")),
        band_1_flux_fraction=band1,
        cloud_absorption_scale=cas,
        reconstruction_check=check,
        results=results,
        closed_form_cross_check=dict(
            note="band-1 layer reflectance, tabled difference against Twomey "
                 "A(1-A)ln2/3, elementwise on layers with tau > 1 and mu0 > 0.05",
            layers=int(sel.sum()),
            median_relative_gap=float(np.nanmedian(rel)),
            p90_relative_gap=float(np.nanpercentile(rel, 90)),
        ),
        droplet_number_factor_reaching_threshold=dict(
            note="the factor in droplet number whose shortwave forcing reaches "
                 "the 1.5 W/m2 reopening threshold, read off the doubling on "
                 "the assumption of linearity in ln(N) over this range",
            low=crossing[0], high=crossing[-1],
        ),
        threshold_w_m2=1.5,
        threshold_source="notes/dust.md reopening test, "
                         "aeolian/config/dust.yaml",
    )
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
