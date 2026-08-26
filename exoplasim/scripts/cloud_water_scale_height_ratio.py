#!/usr/bin/env python
"""What an Earth cloud-water scale height does to the modelled liquid water.

WORLDBUILDING FRAME. Vesper is a simulated super-Earth. Every number below is a
property of a toy GCM's diagnostic cloud scheme run on that simulated planet, or
of the Earth configuration of the same scheme that the constant came from.
Nothing here is a measurement of anything outside the simulation.

THE QUESTION. `rainmod`'s `mkclouds` distributes the modelled in-cloud liquid
water as `rho_l(z) = clwref * exp(-z/hl)` with `hl = clwhsc * ln(1 + PW)`. CCM3's
`clwhsc` is 700 m, a length fitted in an Earth configuration, while the heights
`z` it is measured against are built hypsometrically and so scale as `gascon/ga`.
`rainini` now DERIVES `clwhsc` as `700*(gascon/ga)/(287.0/9.80665)`, which
reproduces 700 m exactly at Earth's constants. This script measures the two
ratios that difference produces, per layer, on the sigma grid every run on
record integrates.

THE TWO RATIOS, kept apart because they answer different questions:

  inherited_over_derived   the modelled layer liquid water with `clwhsc` held at
                           700 m, over the same layer with `clwhsc` derived.
                           Same planet, same sigma level, same temperature
                           profile, same column water. This is what the defect
                           cost, and it is >= 1 in every layer by construction:
                           the inherited length is longer, so it lifts water up
                           the column and takes none away.

  vesper_over_earth        the modelled layer liquid water on Vesper's gravity
                           and gas constant, over the SAME SCHEME evaluated at
                           Earth's `gascon` = 287.0 and `ga` = 9.80665 -- the
                           configuration in which the derivation returns 700 m
                           -- at the same sigma level, the same temperature
                           profile and the same precipitable water. It carries
                           the gravity ratio that the derived configuration also
                           carries, so it is NOT the cost of anything; it is what
                           a shallower atmosphere does to a layer's water.

THE IDENTITY THIS RESTS ON, and the check that can fail. `zzf` scales as
`gascon/ga` and the derived `hl` scales as `gascon/ga`, so the exponent
`-zzf/hl` is INVARIANT between the two planets at the same sigma level. Under
the derived coefficient `vesper_over_earth` is therefore exactly
`(gascon/ga)/(287.0/9.80665)` in EVERY layer, for any temperature profile and
any column water. The script computes both sides independently and raises if
they differ. That constant is also why the two ratios above are the same
quantity up to a fixed factor -- and why a table that reports one while its
prose reads the other is wrong by that factor.

Run it as `python exoplasim/scripts/cloud_water_scale_height_ratio.py`.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset

import _paths  # noqa: F401  anchors every path on this file and adds lib/

import cloud_optical_depth_bracket as codb  # noqa: E402
import lapse  # noqa: E402
from paths import rel  # noqa: E402

ROOT = _paths.PROJECT_ROOT
DEFAULT_CLIM = _paths.ANALYSIS / "climatology" / "bootstrap_regular_climatology.nc"
DEFAULT_CROSS = _paths.ANALYSIS / "cloud_optical_depth_bracket.json"

# CCM3's own constants, the ones the derivation reproduces. Kiehl et al. (1998)
# Eq. 4: a = 700 m. The gas constant and gravity beside them are the Earth
# configuration of the same model, which is what `rainini` divides by.
CLWHSC_CCM3_M = 700.0
GASCON_EARTH = 287.0
GA_EARTH = 9.80665

# The column the finding's table is stated at. Fixed here rather than chosen
# from a result: it is the condition the superseded table carried, so the redo
# is comparable to it.
PRECIPITABLE_WATER_KG_M2 = 25.0


def neqsig4_sigma(nlev: int, ptop_pa: float, psurf_pa: float):
    """(sigma, sigmah) from plasim.f90's `neqsig == 4` construction.

    Transcribed from `plasim.f90` rather than read off a run, so that the run's
    own axis can be checked against it. The quartic has no quadratic term and
    two of its three coefficients are fixed by `sigmah(NLEV) = 1` and a zero
    slope at the surface; the three lines after it are what `neqsig == 4` adds
    over the `neqsig == 0` fallback: shift the top half-level to zero,
    normalise, and map affinely onto [ptop/psurf, 1].
    """
    zsk = np.arange(1, nlev + 1, dtype=float) / nlev
    sigmah = 0.75 * zsk + 1.75 * zsk ** 3 - 1.5 * zsk ** 4
    top = ptop_pa / psurf_pa
    sigmah = sigmah - sigmah[0]
    sigmah = sigmah / sigmah[-1]
    sigmah = sigmah * (1.0 - top) + top
    sigma = np.empty(nlev)
    prev = 0.0
    for j in range(nlev):
        sigma[j] = 0.5 * (prev + sigmah[j])
        prev = sigmah[j]
    return sigma, sigmah


def mid_layer_heights(ta, sigma, sigmah, gascon, ga):
    """`mkclouds`'s `zzf` for a single column. ta is (nlev,), metres out."""
    return codb.layer_heights(np.asarray(ta, dtype=float)[None, :, None, None],
                              sigma, sigmah, gascon, ga)[0, :, 0, 0]


def layer_water_path(ta, sigma, sigmah, dsigma, gascon, ga, ps_pa,
                     clwref, clwhsc, prw):
    """Layer cloud liquid water path, g m-2, for a single column.

    `mkclouds` builds the mixing ratio and `swr` turns it into a path, so the
    surface pressure cancels exactly as it does in the model. No floor and no
    cap: both are applied in the model and neither binds on a ratio taken at
    the same layer of two configurations that differ only in `clwhsc`.
    """
    zzf = mid_layer_heights(ta, sigma, sigmah, gascon, ga)
    hl = clwhsc * math.log(1.0 + prw)
    dql = clwref * np.exp(-zzf / hl) * gascon * np.asarray(ta, dtype=float) \
        / (sigma * ps_pa)
    return 1000.0 * dql * ps_pa / ga * dsigma


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--climatology", type=Path, default=DEFAULT_CLIM)
    ap.add_argument("--run", type=Path, default=None,
                    help="run directory whose namelists supply the constants")
    ap.add_argument("--cross-check", type=Path, default=DEFAULT_CROSS,
                    help="layer water paths this recomputation must reproduce")
    args = ap.parse_args()

    cfg = yaml.safe_load(_paths.CONFIG.read_text(encoding="utf-8"))
    ga = float(cfg["planet"]["gravity_m_s2"])
    gascon, _cp = lapse.gas_properties(cfg)
    nlev = int(cfg["model"]["layers"])
    ptop_pa = float(cfg["model"]["model_top_hpa"]) * 100.0
    vertical_grid = int(cfg["model"]["vertical_grid"])
    psurf_pa = 1.0e5 * sum(v for k, v in cfg["atmosphere"].items()
                           if k.startswith("p") and k.endswith("_bar"))

    # ---------------------------------------------------------------- checks
    # 1. The configured planet is the one the runs integrated. A mismatch here
    #    means the table would be computed for a planet no run has.
    run_dir = args.run
    if run_dir is None:
        runs = sorted((ROOT / "exoplasim" / "runs").glob("run_*"))
        runs = [r for r in runs if (r / "planet_namelist").exists()]
        if not runs:
            raise SystemExit("no run directory carries namelists; pass --run")
        run_dir = runs[-1]
    nl = codb.read_namelists(run_dir)
    for key, want in (("GA", ga), ("GASCON", gascon)):
        got = float(nl[key])
        if abs(got - want) > 1e-6 * abs(want):
            raise SystemExit(
                f"{key} in {rel(run_dir)} is {got!r}, config says {want!r}")
    if int(float(nl["NEQSIG"])) != vertical_grid:
        raise SystemExit(f"{rel(run_dir)} ran NEQSIG = {nl['NEQSIG']}, "
                         f"config declares vertical_grid {vertical_grid}")
    if abs(float(nl["PTOP"]) - ptop_pa) > 1e-6 * ptop_pa:
        raise SystemExit(f"{rel(run_dir)} ran PTOP = {nl['PTOP']}, "
                         f"config declares {ptop_pa} Pa")
    clwref = float(nl.get("CLWREF", 0.00021))
    clwhsc_derived = CLWHSC_CCM3_M * (gascon / ga) / (GASCON_EARTH / GA_EARTH)
    if "CLWHSC" in nl and float(nl["CLWHSC"]) > 0.0:
        got = float(nl["CLWHSC"])
        if abs(got - clwhsc_derived) > 1e-6 * clwhsc_derived:
            raise SystemExit(f"{rel(run_dir)} ran CLWHSC = {got!r}, the "
                             f"derivation gives {clwhsc_derived!r}")

    # 2. The transcribed sigma construction against the axis the model wrote.
    sigma, sigmah = neqsig4_sigma(nlev, ptop_pa, psurf_pa)
    dsigma = np.diff(np.concatenate([[0.0], sigmah]))
    with Dataset(args.climatology) as ds:
        lev = np.asarray(ds.variables["lev"][:], dtype=float)
        centres = np.asarray(ds.variables["time"][:], dtype=float)
        ta_field = np.asarray(ds.variables["ta"][:], dtype=float)
        prw_field = np.asarray(ds.variables["prw"][:], dtype=float)
        ps_field = np.asarray(ds.variables["ps"][:], dtype=float) * 100.0
        run_id = getattr(ds, "vesper_run_id", None)
    if lev.size != nlev:
        raise SystemExit(f"the climatology has {lev.size} layers, config {nlev}")
    axis_err = float(np.max(np.abs(sigma - lev)))
    if axis_err > 1e-5:
        raise SystemExit(
            "the NEQSIG = 4 construction does not reproduce the model's own "
            f"layer midpoints: worst difference {axis_err}")

    # 3. THE IDENTITY. Under the derived coefficient the Vesper-over-Earth layer
    #    water path is exactly (gascon/ga)/(287.0/9.80665) in every layer, for
    #    any temperature profile and any column water, because zzf and hl carry
    #    the same gascon/ga. Computed both ways and compared.
    ntime, _, nlat, nlon = ta_field.shape
    gw = codb.gaussian_weights(nlat, nlon)
    edges = np.concatenate([[0.0], 0.5 * (centres[1:] + centres[:-1]),
                            [centres[-1] + 0.5 * (centres[-1] - centres[-2])]])
    tw = np.diff(edges)
    w = tw[:, None, None] * gw[None, :, :]
    w = w / w.sum()
    ta_mean = np.einsum("tlyx,tyx->l", ta_field, w)
    ps_mean = float(np.einsum("tyx,tyx->", ps_field, w))
    scale_ratio = (gascon / ga) / (GASCON_EARTH / GA_EARTH)
    derived_v = layer_water_path(ta_mean, sigma, sigmah, dsigma, gascon, ga,
                                 ps_mean, clwref, clwhsc_derived,
                                 PRECIPITABLE_WATER_KG_M2)
    derived_e = layer_water_path(ta_mean, sigma, sigmah, dsigma,
                                 GASCON_EARTH, GA_EARTH, ps_mean, clwref,
                                 CLWHSC_CCM3_M, PRECIPITABLE_WATER_KG_M2)
    identity_err = float(np.max(np.abs(derived_v / derived_e - scale_ratio)))
    if identity_err > 1e-9:
        raise SystemExit(
            "the derived coefficient does not make the Vesper-over-Earth layer "
            f"water path constant: worst departure {identity_err}")

    # 4. CROSS-CHECK against a quantity the other side already knows. The layer
    #    water paths in cloud_optical_depth_bracket.json are this same scheme at
    #    the derived coefficient over the same climatology, computed by a
    #    different script from the model's own per-cell fields.
    #    It is run on the grid THAT script inverts out of the climatology's
    #    `lev` axis rather than on the transcription above, so that the
    #    comparison is of the scheme and not of the axis's float32 storage; the
    #    two grids agree to `axis_err`, which is checked separately.
    #    It also takes `gascon` from the namelist rather than from the config,
    #    which is what the other side read: the namelist carries the derived gas
    #    constant to fewer digits than the composition does, and a comparison
    #    across that difference measures a namelist's decimal precision instead
    #    of the scheme. The two agree to `gascon_precision`, which is recorded.
    cross = json.loads(args.cross_check.read_text(encoding="utf-8"))
    gascon_nl = float(nl["GASCON"])
    gascon_precision = abs(gascon_nl / gascon - 1.0)
    sigma_x, sigmah_x, dsigma_x = codb.sigma_grid(lev)
    clwhsc_x = CLWHSC_CCM3_M * (gascon_nl / ga) / (GASCON_EARTH / GA_EARTH)
    lwp_field = codb.water_paths(ta_field, ps_field, prw_field, sigma_x,
                                 sigmah_x, dsigma_x, gascon_nl, ga, clwref,
                                 clwhsc_x)
    # The reduction is the other script's, an unweighted mean over time and
    # cell, because what is being reproduced is ITS number and not a global
    # mean: a different reduction would leave the two disagreeing for a reason
    # that says nothing about the scheme.
    lwp_mean = lwp_field.mean(axis=(0, 2, 3))
    cross_mean = np.array([row["water_path_g_m2_mean"] for row in cross["layers"]])
    cross_err = float(np.max(np.abs(lwp_mean / cross_mean - 1.0)))
    if cross_err > 1e-12:
        raise SystemExit(
            "the recomputed layer water paths disagree with "
            f"{rel(args.cross_check)}: worst relative difference {cross_err}")

    # ----------------------------------------------------------------- table
    inherited_v = layer_water_path(ta_mean, sigma, sigmah, dsigma, gascon, ga,
                                   ps_mean, clwref, CLWHSC_CCM3_M,
                                   PRECIPITABLE_WATER_KG_M2)
    r_cost = inherited_v / derived_v
    r_ve = inherited_v / derived_e

    # The same ratio over the climatology's own columns, so the table's single
    # condition can be placed against the spread the modelled atmosphere has.
    lwp_inherited = codb.water_paths(ta_field, ps_field, prw_field, sigma,
                                     sigmah, dsigma, gascon, ga, clwref,
                                     CLWHSC_CCM3_M)
    wet = prw_field > 1.0
    cell_ratio = np.where(wet[:, None], lwp_inherited / lwp_field, np.nan)
    lo = np.nanpercentile(cell_ratio, 10, axis=(0, 2, 3))
    hi = np.nanpercentile(cell_ratio, 90, axis=(0, 2, 3))

    zzf_mean = mid_layer_heights(ta_mean, sigma, sigmah, gascon, ga)

    # NOTHING IS WRITTEN. This is a reproducer for one table in a dated finding,
    # not a generator: its numbers live in notes/audits/model-earth-centrism.md,
    # which is where a finding's numbers belong, and no step consumes them. An
    # artifact no step in config/pipeline.yaml generates does not exist, so the
    # honest form of a calculation with no consumer is one that prints.
    print(f"climatology     {rel(args.climatology)}  ({run_id})")
    print(f"constants from  {rel(run_dir)}")
    print(f"condition       {PRECIPITABLE_WATER_KG_M2:g} kg/m2 precipitable "
          f"water, global mean ta, ps {ps_mean:.0f} Pa")
    print(f"clwhsc          inherited {CLWHSC_CCM3_M:.1f} m, "
          f"derived {clwhsc_derived:.2f} m, ratio {scale_ratio:.6f}")
    print(f"hl at this PW   inherited "
          f"{CLWHSC_CCM3_M * math.log(1.0 + PRECIPITABLE_WATER_KG_M2):.1f} m, "
          f"derived "
          f"{clwhsc_derived * math.log(1.0 + PRECIPITABLE_WATER_KG_M2):.1f} m")
    print("checks")
    print(f"  NEQSIG = 4 construction against the model's lev axis  {axis_err:.3e}")
    print(f"  the identity, worst departure from the constant       {identity_err:.3e}")
    print(f"  layer water paths against {rel(args.cross_check)}"
          f"  {cross_err:.3e}")
    print(f"  config gascon over the namelist's, minus one          {gascon_precision:.3e}")
    print("\n lev    sigma    z (m)   derived  inherited   inh/der    V/E"
          "    p10    p90")
    for j in range(nlev):
        print(f" {j + 1:3d}  {sigma[j]:.5f}  {zzf_mean[j]:7.0f}  "
              f"{derived_v[j]:8.3f}  {inherited_v[j]:9.3f}  "
              f"{r_cost[j]:7.3f}  {r_ve[j]:5.3f}  {lo[j]:5.3f}  {hi[j]:5.3f}")
    print(f" col                      {derived_v.sum():8.1f}  "
          f"{inherited_v.sum():9.1f}  "
          f"{inherited_v.sum() / derived_v.sum():7.3f}  "
          f"{inherited_v.sum() / derived_e.sum():5.3f}")


if __name__ == "__main__":
    main()
