#!/usr/bin/env python
"""What the band-1 cloud optical depth correction is worth, over a real field.

WORLDBUILDING FRAME. This is a simulated planet, Vesper, and every number below
is a property of the model or of the paper the model's shortwave cloud optics
are ported from. Nothing here is a measurement of anything outside the
simulation.

THE QUESTION. `radmod`'s shortwave computed one cloud optical depth,
`2*ALOG10(1.5 + W)**3.9`, and used it for both shortwave bands. Stephens (1978)
p. 2124 fits two, Eq. (10a) below 0.75 um and Eq. (10b) above it, and `radmod`
splits its own bands at the same 0.75 um. world-jgen replaced the single value
with the paper's two and continued each below the fitted range linearly in the
water path. This script asks what that swap is worth at the top of the
atmosphere, and it is a BRACKET rather than a measurement: only a run settles
it, and the run is named in exoplasim/notes/cloud-water-reference.md.

WHY IT IS BETTER THAN THE HAND CHAIN IT REPLACES. The estimate in the issue that
opened this had no cloud fraction field, no overlap and no surface in it, and
leaned high for want of all three. This one carries the model's own per-layer
cloud water path, its own total cloud cover with its own random overlap, its own
zenith-angle geometry integrated over the day, and the band-1 surface albedo
through Stephens Eq. (12).

WHAT THE BRACKET IS OVER. The per-layer cloud fraction is not written out at any
cadence, so the AMOUNT of cloud is taken from `clt` and the arms differ only in
WHERE in the column it sits. The spread between them is the bracket's width.

WHAT ELSE MAKES IT A BRACKET AND NOT A MEASUREMENT:

  1. The field is a BOOTSTRAP climatology, on terrain-only surface fields, and
     the cloud in it was made by the optical depth being corrected. Under the
     correction the cloud and the humidity move too, which is exactly what only
     a run can say.
  2. Band-1 gas absorption and Rayleigh scattering above and below the cloud
     are not in this chain. Both attenuate the change, so this leans HIGH.
  3. The column is composed by random overlap without multiple reflection
     between layers, which also leans high.
  4. Stephens (1978) p. 2127 warns that the cloud-over-reflecting-surface
     correction is unreliable above a surface albedo of about 0.75, which is
     every snow and sea-ice cell.

Run it as `python exoplasim/scripts/cloud_optical_depth_bracket.py`.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

import _paths  # noqa: F401  anchors every path on this file and adds lib/

import sensitivity  # noqa: E402
import stellar  # noqa: E402
from paths import rel  # noqa: E402

ROOT = _paths.PROJECT_ROOT
DEFAULT_CLIM = _paths.ANALYSIS / "climatology" / "bootstrap_regular_climatology.nc"
DEFAULT_OUT = _paths.ANALYSIS / "cloud_optical_depth_bracket.json"

# --------------------------------------------------------------------------
# The two optical depths, and the one they replace.
#
# Stephens (1978) p. 2124: log10(tau_N) = a + b ln(log10 W), W in g m-2, with
# (a, b) = (0.2633, 1.7095) for 0.3 to 0.75 um and (0.3492, 1.6518) for 0.75 to
# 4.0 um. 10**(a + b ln x) is 10**a x**(b ln 10), so each fit is a power law in
# log10 W. The fitted range is Figs. 1a and 1b: 10 to 10000 g m-2.
# --------------------------------------------------------------------------
A1, P1 = 10.0 ** 0.2633, 1.7095 * math.log(10.0)   # Eq. (10a), band 1
A2, P2 = 10.0 ** 0.3492, 1.6518 * math.log(10.0)   # Eq. (10b), band 2
WFIT = 10.0                                        # g m-2, bottom of the fit


def tau_inherited(w):
    """The single optical depth radmod used for both bands before world-jgen."""
    return 2.0 * np.log10(1.5 + np.maximum(0.0, w)) ** 3.9


def tau_band(w, prefactor, exponent):
    """One band's optical depth, fit above WFIT and linear in W below it.

    Below WFIT the fit is extrapolation, and extending its own form makes the
    optical depth fall faster than linearly, which through Stephens Eq. (7),
    tau_N = 1.5 W / r_e, is a droplet effective radius running past 100 um as
    the cloud thins. The linear branch is Eq. (7) at the effective radius the
    fit itself implies at WFIT, and it is continuous there.
    """
    w = np.maximum(0.0, w)
    x = np.log10(np.maximum(WFIT, w))
    return prefactor * x ** exponent * np.minimum(1.0, w / WFIT)


# --------------------------------------------------------------------------
# The model's own constants, read from the run rather than restated.
# --------------------------------------------------------------------------

def read_namelists(run_dir: Path) -> dict:
    """Every KEY = VALUE across a run's namelists, upper-cased keys."""
    out: dict[str, str] = {}
    for path in sorted(run_dir.glob("*namelist")):
        for line in path.read_text(encoding="utf-8").splitlines():
            if "=" not in line or line.strip().startswith("&"):
                continue
            key, _, value = line.partition("=")
            out[key.strip().upper()] = value.strip().rstrip(",")
    return out


def sigma_grid(lev: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(sigma, sigmah, dsigma) recovered from the layer midpoints.

    plasim.f90 builds sigma from sigmah as sigma(1) = 0.5*sigmah(1) and
    sigma(j) = 0.5*(sigmah(j-1) + sigmah(j)), which inverts exactly. The
    inversion alternates in sign, so its error shows up entirely in the closure
    sigmah(NLEV) = 1: that is checked rather than assumed, because a silently
    wrong sigma grid would put the cloud water in the wrong layers.
    """
    sigmah = np.empty_like(lev)
    prev = 0.0
    for j, s in enumerate(lev):
        sigmah[j] = 2.0 * s - prev
        prev = sigmah[j]
    closure = abs(sigmah[-1] - 1.0)
    if closure > 1e-4:
        raise SystemExit(
            f"the sigma grid does not close: sigmah(NLEV) = {sigmah[-1]!r}, "
            "so `lev` is not this model's layer midpoints")
    dsigma = np.diff(np.concatenate([[0.0], sigmah]))
    return lev.copy(), sigmah, dsigma


def layer_heights(ta, sigma, sigmah, gascon, ga):
    """Mid-layer height above the surface, as rainmod's mkclouds builds it.

    ta is (time, lev, lat, lon). Returns the same shape.
    """
    nlev = ta.shape[1]
    zzf = np.zeros_like(ta)
    running = np.zeros(ta.shape[:1] + ta.shape[2:], dtype=ta.dtype)
    for j in range(nlev - 1, 0, -1):
        zdh = -ta[:, j] * gascon / ga * math.log(sigmah[j - 1] / sigmah[j])
        zzf[:, j] = running + 0.5 * zdh
        running = running + zdh
    zdh = -ta[:, 0] * gascon / ga * math.log(sigma[0] / sigmah[0]) * 0.5
    zzf[:, 0] = running + zdh
    return zzf


def water_paths(ta, ps_pa, prw, sigma, sigmah, dsigma, gascon, ga,
                clwref, clwhsc):
    """Per-layer cloud liquid water path in g m-2, as mkclouds and swr make it.

    The surface pressure cancels between dql and the path, exactly as it does
    in the model: dql carries 1/dp and swr multiplies by dp.
    """
    zzf = layer_heights(ta, sigma, sigmah, gascon, ga)
    zzh = clwhsc * np.log(1.0 + prw)                       # (time, lat, lon)
    arg = np.clip(-zzf / np.maximum(1e-30, zzh)[:, None], -80.0, 80.0)
    dql = clwref * np.exp(arg) * gascon * ta \
        / (sigma[None, :, None, None] * ps_pa[:, None])
    dql = np.maximum(dql, 1e-9)
    lwp = 1000.0 * dql * ps_pa[:, None] / ga * dsigma[None, :, None, None]
    return np.minimum(1000.0, lwp)


def cloud_fraction(clt, shape):
    """Per-layer cloud fraction carrying the model's own total cover.

    THE PER-LAYER COVER IS NOT WRITTEN OUT and it cannot be recovered from what
    is. rainmod diagnoses non-convective cover as ((rh - rcrit)/(1 - rcrit))**2
    with rcrit = max(0.85, max(sigma, 1 - sigma)), and `hur` in a climatology is
    a mean over time bins and orbits. That square of a threshold difference is
    strongly convex, so the mean of the diagnosis is not the diagnosis of the
    mean: run on this field it leaves 72 per cent of columns with no cloud at
    all against a written total cover of 0.61. The convective branch is worse
    off still, since the flags it needs are not written at any cadence.

    So the AMOUNT of cloud is taken from `clt`, which the model does write, and
    the unknown is reduced to WHERE in the column it sits. `shape` is a relative
    vertical profile, scaled here so that the column's random-overlap total
    is `clt` exactly. The arms that call this differ only in that profile, and
    the spread between them is the bracket's width.

    clt is (time, lat, lon); shape is (time, lev, lat, lon) or (lev,).
    """
    shape = np.asarray(shape, dtype=float)
    if shape.ndim == 1:
        shape = shape[None, :, None, None] * np.ones((1, 1) + clt.shape[1:])
    peak = np.maximum(1e-30, shape.max(axis=1, keepdims=True))
    unit = np.clip(shape / peak, 0.0, 1.0)
    # One scalar k per column, from 1 - prod(1 - k*unit) = clt. The left side
    # rises monotonically from 0 to 1 - prod(1 - unit) as k goes 0 to 1, so a
    # bisection on k is exact to the tolerance it is run to. Where the target is
    # above what k = 1 can reach, k saturates and the shortfall is reported.
    target = np.asarray(clt, dtype=float)
    lo = np.zeros_like(target)
    hi = np.ones_like(target)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        got = 1.0 - np.prod(1.0 - mid[:, None] * unit, axis=1)
        lo = np.where(got < target, mid, lo)
        hi = np.where(got < target, hi, mid)
    dcc = np.clip(0.5 * (lo + hi)[:, None] * unit, 0.0, 1.0)
    return dcc, 1.0 - np.prod(1.0 - dcc, axis=1)


def band1_column_reflectance(tau, dcc, mu0, tswr1, zmu00=0.5):
    """Band-1 cloud reflectance of the column, composed the way swr composes it.

    Per layer, Stephens Eq. (1) in radmod's form: Re = 1 - 1/(1 + beta tau/mu0),
    with beta = tswr1 sqrt(mu0) for the direct stream and the diffuse stream
    evaluated at zmu00. The two are blended by the clear-sky fraction above the
    layer, which is radmod's `zcs`, and the column is the random-overlap
    product of what each layer leaves unreflected.

    tau and dcc are (lev, lat, lon); mu0 is (lat,) and may hold zeros.
    """
    mu = np.maximum(1e-30, mu0)[None, :, None]
    b1 = tswr1 * np.sqrt(mu)
    b3 = tswr1 * math.sqrt(zmu00) / zmu00
    direct = 1.0 - 1.0 / (1.0 + b1 * tau / mu)
    diffuse = 1.0 - 1.0 / (1.0 + b3 * tau)
    clear_above = np.concatenate(
        [np.ones_like(dcc[:1]), np.cumprod(1.0 - dcc, axis=0)[:-1]], axis=0)
    layer = clear_above * direct + (1.0 - clear_above) * diffuse
    return 1.0 - np.prod(1.0 - dcc * layer, axis=0)


def over_surface(reflectance, albedo):
    """Stephens (1978) Eq. (12), cloud over a reflecting surface, band 1.

    Band 1 is the conservative-scattering band, so Tr = 1 - Re, and the system
    albedo is Re + a Tr**2 / (1 - a Re). This is what turns a change in cloud
    reflectance into a change at the top of the atmosphere: the surface returns
    part of whatever the cloud now lets through, so it always reduces the
    magnitude. Stephens p. 2127 warns that (12) is least reliable above a
    surface albedo of about 0.75, which is every snow and sea-ice cell here.
    """
    tr = 1.0 - reflectance
    return reflectance + albedo * tr * tr / np.maximum(1e-6, 1.0 - albedo * reflectance)


def gaussian_weights(nlat, nlon):
    from numpy.polynomial.legendre import leggauss
    return leggauss(nlat)[1][::-1][:, None] * np.ones((1, nlon))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--climatology", type=Path, default=DEFAULT_CLIM)
    ap.add_argument("--run", type=Path, default=None,
                    help="run directory whose namelists supply the constants")
    ap.add_argument("--hours", type=int, default=48,
                    help="hour angles per day in the zenith integration")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    run_dir = args.run
    if run_dir is None:
        runs = sorted((ROOT / "exoplasim" / "runs").glob("run_*"))
        runs = [r for r in runs if (r / "planet_namelist").exists()]
        if not runs:
            raise SystemExit("no run directory carries namelists; pass --run")
        run_dir = runs[-1]
    nl = read_namelists(run_dir)
    ga = float(nl["GA"])
    gascon = float(nl["GASCON"])
    gsol0 = float(nl["GSOL0"])
    clwref = float(nl.get("CLWREF", 0.00021))
    clwhsc = float(nl["CLWHSC"]) if "CLWHSC" in nl \
        else 700.0 * (gascon / ga) / (287.0 / 9.80665)
    tswr1 = float(nl.get("TSWR1", 0.077))

    zsolar1 = float(stellar.band1_fraction())

    with Dataset(args.climatology) as ds:
        lev = np.asarray(ds.variables["lev"][:], dtype=float)
        lat = np.asarray(ds.variables["lat"][:], dtype=float)
        centres = np.asarray(ds.variables["time"][:], dtype=float)
        zdec = np.asarray(ds.variables["zdec"][:], dtype=float)
        rdist = np.asarray(ds.variables["rdist"][:], dtype=float)
        ta = np.asarray(ds.variables["ta"][:], dtype=float)
        prw = np.asarray(ds.variables["prw"][:], dtype=float)
        ps = np.asarray(ds.variables["ps"][:], dtype=float) * 100.0
        clt = np.asarray(ds.variables["clt"][:], dtype=float)
        alb1 = np.asarray(ds.variables["alb1"][:], dtype=float)
        rst = np.asarray(ds.variables["rst"][:], dtype=float)
        rsut = np.abs(np.asarray(ds.variables["rsut"][:], dtype=float))
        run_id = getattr(ds, "vesper_run_id", None)
        flux_ratio = float(getattr(ds, "vesper_flux_ratio", float("nan")))

    ntime, nlev, nlat, nlon = ta.shape
    sigma, sigmah, dsigma = sigma_grid(lev)
    lwp = water_paths(ta, ps, prw, sigma, sigmah, dsigma, gascon, ga,
                      clwref, clwhsc)

    gw = gaussian_weights(nlat, nlon)
    # The climatology's bins hold unequal numbers of raw records; weight the
    # annual mean by the bin widths the time centres imply.
    edges = np.concatenate([[0.0], 0.5 * (centres[1:] + centres[:-1]),
                            [centres[-1] + 0.5 * (centres[-1] - centres[-2])]])
    tw = np.diff(edges)
    tw = tw / tw.sum()

    tau_old = tau_inherited(lwp)
    tau_new1 = tau_band(lwp, A1, P1)

    # THE ARMS DIFFER ONLY IN WHERE THE COVER SITS. Every one carries the
    # model's own `clt` and the model's own water paths; the spread between
    # them is what not writing the per-layer cover costs this estimate.
    #
    #   uniform         the cover spread evenly over the ten layers
    #   water_weighted  the cover in proportion to each layer's cloud water,
    #                   which is where the model's own condensate is
    #   thickest_layer  all of it in the layer holding the most water, the
    #                   optically thickest cloud and so the least sensitive
    #   largest_move    all of it in the layer whose band-1 cloud reflectance
    #                   moves most, evaluated at the diffuse-stream cosine
    b3 = tswr1 * math.sqrt(0.5) / 0.5
    per_layer_move = np.abs(
        (1.0 - 1.0 / (1.0 + b3 * tau_new1))
        - (1.0 - 1.0 / (1.0 + b3 * tau_old))).mean(axis=(0, 2, 3))
    shapes = {
        "uniform": np.ones(nlev),
        "water_weighted": lwp,
        "thickest_layer": np.eye(nlev)[int(np.argmax(lwp.mean(axis=(0, 2, 3))))],
        "largest_move": np.eye(nlev)[int(np.argmax(per_layer_move))],
    }

    results = {}
    for label, shape in shapes.items():
        dcc, total = cloud_fraction(clt, shape)
        d_reflected = np.zeros((ntime, nlat, nlon))
        d_reflected_black = np.zeros((ntime, nlat, nlon))
        incident = np.zeros((ntime, nlat, nlon))
        r_old_mean = np.zeros((ntime, nlat, nlon))
        r_new_mean = np.zeros((ntime, nlat, nlon))
        hours = (np.arange(args.hours) + 0.5) / args.hours * 2.0 * math.pi \
            - math.pi
        for t in range(ntime):
            dec = math.radians(zdec[t])
            f1 = zsolar1 * gsol0 / rdist[t] ** 2
            for h in hours:
                mu0 = np.sin(np.radians(lat)) * math.sin(dec) \
                    + np.cos(np.radians(lat)) * math.cos(dec) * math.cos(h)
                mu0 = np.where(mu0 > 0.0, mu0, 0.0)
                if not np.any(mu0 > 0.0):
                    continue
                # The weight IS the incident flux: a night lane carries mu0 = 0
                # and contributes nothing, so no separate day mask is needed.
                w = f1 * (mu0[:, None] * np.ones((1, nlon))) / args.hours
                r_old = band1_column_reflectance(tau_old[t], dcc[t], mu0, tswr1)
                r_new = band1_column_reflectance(tau_new1[t], dcc[t], mu0, tswr1)
                d_toa = over_surface(r_new, alb1[t]) - over_surface(r_old, alb1[t])
                d_reflected[t] += w * d_toa
                d_reflected_black[t] += w * (r_new - r_old)
                incident[t] += w
                r_old_mean[t] += w * r_old
                r_new_mean[t] += w * r_new

        def gmean(field):
            return float((np.tensordot(tw, field, axes=(0, 0)) * gw).sum()
                         / gw.sum())

        # Flux-weighted means of the column's band-1 cloud reflectance, so the
        # forcing can be read beside the reflectance change that produced it.
        # Both were accumulated against the same incident-flux weight.
        results[label] = {
            "cloud_fraction_total_mean": gmean(total),
            "band1_incident_w_m2": gmean(incident),
            "band1_cloud_reflectance_old": gmean(r_old_mean) / gmean(incident),
            "band1_cloud_reflectance_new": gmean(r_new_mean) / gmean(incident),
            "d_reflected_w_m2": gmean(d_reflected),
            "d_reflected_black_surface_w_m2": gmean(d_reflected_black),
            "forcing_w_m2": -gmean(d_reflected),
            "forcing_black_surface_w_m2": -gmean(d_reflected_black),
        }

    # The instrument's own check: the reconstructed band-1 incident flux against
    # the climatology's own top-of-atmosphere shortwave. rst is net down and
    # rsut the upward part, so their sum is the incident flux, and zsolar1 of it
    # is band 1. A reconstruction that missed the geometry fails here.
    incident_model = float(((np.tensordot(tw, rst + rsut, axes=(0, 0))) * gw).sum()
                           / gw.sum()) * zsolar1
    incident_recon = results["uniform"]["band1_incident_w_m2"]
    geometry_error = incident_recon / incident_model - 1.0

    alpha = sensitivity.planetary_albedo_from_fluxes(
        float((np.tensordot(tw, rst, axes=(0, 0)) * gw).sum() / gw.sum()),
        float((np.tensordot(tw, rsut, axes=(0, 0)) * gw).sum() / gw.sum()))

    # The bracket is over WHERE the cover sits, with the surface in every arm.
    # The black-surface figures stay in the JSON beside them because they say
    # what the surface is worth, but a bracket end without a surface would be a
    # claim about a planet this is not.
    forcings = [r["forcing_w_m2"] for r in results.values()]
    lo, hi = min(forcings), max(forcings)
    kelvin = [sensitivity.forcing_to_kelvin(f, alpha) for f in (lo, hi)]

    # The per-layer table is the finding, and the bracket is what it sums to.
    # A layer thinner than about 4 g m-2 moves the OTHER way: its inherited
    # optical depth came almost entirely from the `1.5 +` offset, and the
    # continuation that replaces the offset returns more than the offset did.
    layers = []
    for j in range(nlev):
        w_mean = float(lwp[:, j].mean())
        layers.append({
            "level": j + 1,
            "sigma": float(sigma[j]),
            "water_path_g_m2_mean": w_mean,
            "water_path_g_m2_p10": float(np.percentile(lwp[:, j], 10)),
            "water_path_g_m2_p90": float(np.percentile(lwp[:, j], 90)),
            "tau_band1_inherited": float(tau_inherited(w_mean)),
            "tau_band1_eq10a": float(tau_band(w_mean, A1, P1)),
            "tau_band2_eq10b": float(tau_band(w_mean, A2, P2)),
            "d_reflectance_diffuse": float(
                (1.0 - 1.0 / (1.0 + b3 * tau_band(w_mean, A1, P1)))
                - (1.0 - 1.0 / (1.0 + b3 * tau_inherited(w_mean)))),
        })

    out = {
        "what": "band-1 cloud optical depth, Stephens Eq. (10b) replaced by "
                "Eq. (10a), at the top of the atmosphere",
        "issue": "world-jgen",
        "climatology": rel(args.climatology),
        "climatology_is_bootstrap": True,
        "run_id": run_id,
        "flux_ratio": flux_ratio,
        "constants_from": rel(run_dir),
        "constants": {
            "ga_m_s2": ga, "gascon_j_kg_k": gascon, "gsol0_w_m2": gsol0,
            "clwref_kg_m3": clwref, "clwhsc_m": clwhsc, "tswr1": tswr1,
            "zsolar1": zsolar1,
        },
        "planetary_albedo": alpha,
        "geometry_check": {
            "band1_incident_reconstructed_w_m2": incident_recon,
            "band1_incident_from_climatology_w_m2": incident_model,
            "relative_error": geometry_error,
        },
        "layers": layers,
        "arms": results,
        "bracket_w_m2": [lo, hi],
        "bracket_kelvin": kelvin,
        "sign": "positive forcing warms; over the paths this model's cloud "
                "water actually reaches the correction makes band-1 cloud "
                "less bright, so it is a warming",
        "leans": "high: no band-1 gas absorption or Rayleigh above or below "
                 "the cloud, and no multiple reflection between layers",
        "settles_it": "nothing here. A T21 pair, named in "
                      "exoplasim/notes/cloud-water-reference.md",
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")

    print(f"climatology   {rel(args.climatology)} (bootstrap)")
    print(f"constants     {rel(run_dir)}")
    print(f"zsolar1       {zsolar1:.4f}")
    print(f"geometry      reconstructed band-1 incident "
          f"{incident_recon:.2f} W/m2 against the climatology's "
          f"{incident_model:.2f} ({geometry_error * 100:+.2f} per cent)")
    print(f"albedo        {alpha:.4f} planetary, from rst and rsut")
    for label, r in results.items():
        print(f"{label:<13} cloud total {r['cloud_fraction_total_mean']:.3f}  "
              f"forcing {r['forcing_w_m2']:+.3f} W/m2  "
              f"(black surface {r['forcing_black_surface_w_m2']:+.3f})")
    print(f"bracket       {lo:+.2f} to {hi:+.2f} W/m2, "
          f"{kelvin[0]:+.2f} to {kelvin[1]:+.2f} K")
    print(f"written       {rel(args.out)}")


if __name__ == "__main__":
    main()
