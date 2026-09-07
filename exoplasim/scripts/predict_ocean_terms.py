#!/usr/bin/env python3
"""Predict what the two ocean namelist terms will do, before any run measures it.

WORLDBUILDING CONTEXT, stated first because this file's vocabulary invites
misreading: Vesper is a FICTIONAL planet, and this script is planetary climate
modelling for that fiction. It prices two switches of a toy GCM's slab ocean --
horizontal heat diffusion, and the freezing point implied by a declared ocean
salinity -- offline, against the current baseline climatology. Nothing here is
Earth science and nothing here is biology.

WHY IT EXISTS. `docs/src/pipeline/sequencing.md` A3: every forcing change lands with a quantitative
prediction of its own effect, stated before it is run, with what result would
mean "wrong". The predictions and their wrongness bounds are written up in
`exoplasim/notes/forcing-bundle-predictions.md`; this script is how the numbers
there were computed, so they can be recomputed when the climatology moves.

CLIM-16, `nhdiff`/`hdiffk`. Reproduces `oceanmod.f90:hdiffo`'s own
discretisation -- anomaly from the ocean mean, flux form, no-flux across
coastlines -- and applies it to the baseline annual-mean SST. The operator is
VALIDATED before it is used: on P2(sin phi), an eigenfunction of the spherical
Laplacian, it must return -6/a^2 times the field, which is a check with a right
answer and it already caught one real error (the latitude axis runs north to
south, which sign-flips the meridional term if forgotten).

CLIM-17, `TFREEZE` from salinity. Counts the ocean whose coldest output bin sits
between the current freezing point and the bracket value, which is the area
newly able to freeze seasonally, and converts through an ice albedo contrast
carried as a BRACKET (0.25 to 0.40 of planetary albedo per unit area: surface
contrast ~0.5 times the budget's measured attenuation, up to weakly attenuated).
That attenuation is 0.38 with a measured span of 0.29 to 0.48, so the low end of
this bracket now sits slightly above what the surface-contrast route gives and
the bracket is kept as it is: it was fixed before the measurement and its width
is about how attenuated sea ice is, not about the constant.

This prints. It writes nothing.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset
from numpy.polynomial.legendre import leggauss

from _paths import CONFIG, PROJECT_ROOT  # noqa: F401  (puts lib/ on sys.path)

import climatology  # noqa: E402  from lib/
import sea_water  # noqa: E402  from lib/
import sensitivity  # noqa: E402  from lib/
from paths import climatology_path  # noqa: E402
from run_exoplasim import freezing_point_k  # noqa: E402
import nc_geometry

ALBEDO_PER_ICE_AREA = (0.25, 0.40) # planetary albedo per unit new ice area; see module docstring

# THE SLAB'S HEAT CAPACITY COMES FROM THE MODEL. Density and specific heat are
# icemod_nl keys that icemod passes to oceanini, so oceanmod does not own them
# and a run can set them. The pair stood here as literals attributed to
# oceanmod.f90 with CPS at 4180, fresh water at about 25 C, while the model had
# moved to sea water's value at S = 34.7 and its freezing point; the prediction
# below is a heat flux per unit temperature change, so it carried their ratio.
_SEA_WATER = sea_water.constants()
CRHOS, CPS = _SEA_WATER["CRHOS"], _SEA_WATER["CPS"]


def hdiffo_tendency(T: np.ndarray, sea: np.ndarray, k: float, a: float,
                    phi: np.ndarray, gw: np.ndarray) -> np.ndarray:
    """dT/dt from oceanmod.f90:hdiffo's discretisation. phi must run SOUTH TO NORTH."""
    nlat, nlon = T.shape
    dlam = 2.0 * np.pi / nlon
    cphi = np.cos(phi)
    both = sea & np.roll(sea, -1, axis=1)
    fx = np.where(both, (np.roll(T, -1, axis=1) - T) / dlam, 0.0)
    phih = 0.5 * (phi[:-1] + phi[1:])
    dphi = phi[1:] - phi[:-1]
    bothy = sea[:-1, :] & sea[1:, :]
    fy = np.zeros((nlat + 1, nlon))
    fy[1:-1, :] = np.where(bothy, (T[1:, :] - T[:-1, :])
                           * np.cos(phih)[:, None] / dphi[:, None], 0.0)
    div = ((fx - np.roll(fx, 1, axis=1)) / dlam / cphi[:, None] ** 2
           + (fy[1:, :] - fy[:-1, :]) / gw[:, None])
    return np.where(sea, (k / a ** 2) * div, 0.0)


def validate_operator(a: float, phi: np.ndarray, gw: np.ndarray, nlon: int) -> float:
    """P2(sin phi) is an eigenfunction: the operator must return -6/a^2 * T."""
    mu = np.sin(phi)
    p2 = np.broadcast_to((0.5 * (3 * mu ** 2 - 1))[:, None], (len(phi), nlon)).copy()
    allsea = np.ones_like(p2, dtype=bool)
    got = hdiffo_tendency(p2, allsea, 1.0, a, phi, gw)
    want = -6.0 / a ** 2 * p2
    err = float(np.abs(got[2:-2] / want[2:-2] - 1.0).max())
    if err > 0.02:
        raise SystemExit(f"operator FAILS its eigenfunction check: {err:.3%} from "
                         "-6/a^2 on P2. Do not quote anything below.")
    return err


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--climatology", type=Path, default=None)
    ap.add_argument("--hdiffk", type=float, nargs="+", default=[300.0, 1000.0, 3000.0],
                    help="m2/s; the middle value is oceanmod.f90's default")
    ap.add_argument("--salinity", type=float, nargs="+", default=[30.0, 25.0, 20.0],
                    help="psu bracket; the declared value is in config/planet.yaml")
    args = ap.parse_args()

    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    a = nc_geometry.planet_radius_m(cfg)
    mld = float(cfg["surface"]["mixed_layer_depth_m"])
    heatcap = CRHOS * CPS * mld
    incident = (float(cfg["orbit"]["earth_solar_constant_w_m2"])
                * float(cfg["orbit"]["baseline_flux_earth"]) / 4.0)
    s_declared = float(cfg["ocean"]["salinity_psu"])
    tf_declared = freezing_point_k(s_declared)
    alpha, alpha_src = sensitivity.planetary_albedo(cfg)

    path = args.climatology or climatology_path()
    with Dataset(path) as ds:
        lat = np.asarray(ds["lat"][:], dtype=float)
        nlat, nlon = len(lat), len(ds["lon"][:])
        order = np.argsort(lat)                     # hdiffo wants south to north
        gw = leggauss(nlat)[1][::-1][order]
        phi = np.deg2rad(lat[order])
        sea = climatology.annual_mean_of(ds, "lsm")[order] < 0.5
        sst = climatology.annual_mean_of(ds, "ts")[order]
        tmin = np.asarray(ds["ts"][:], dtype=float).min(axis=0)[order]
        run_id = getattr(ds, "vesper_run_id", "?")

    area_w = np.broadcast_to(gw[:, None], (nlat, nlon)) / (2.0 * nlon)
    w_sea = area_w * sea
    latx = lat[order]

    err = validate_operator(a, phi, gw, nlon)
    print(f"operator check: P2 eigenfunction reproduced to {err:.2%}  "
          f"(climatology {path.name}, run {run_id})")
    print(f"planetary albedo {alpha:.4f} from {alpha_src}; incident mean {incident:.1f} W/m2\n")

    print(f"CLIM-16  nhdiff heating implied by the baseline SST, {mld:g} m slab:")
    t_anom = np.where(sea, sst - float((sst * w_sea).sum() / w_sea.sum()), 0.0)
    for k in args.hdiffk:
        h = heatcap * hdiffo_tendency(t_anom, sea, k, a, phi, gw)
        rms = float(np.sqrt((h ** 2 * w_sea).sum() / w_sea.sum()))
        gmean = float((h * w_sea).sum() / w_sea.sum())
        print(f"  hdiffk {k:6.0f} m2/s: rms {rms:6.3f} W/m2 of ocean, "
              f"|max| {np.abs(h[sea]).max():7.2f}, global mean {gmean:+.1e} (conserves)")
    h = heatcap * hdiffo_tendency(t_anom, sea, 1000.0, a, phi, gw)
    for lo, hi in ((0, 20), (20, 40), (40, 60), (60, 90)):
        m = sea & (np.abs(latx)[:, None] >= lo) & (np.abs(latx)[:, None] < hi)
        print(f"    |lat| {lo:2d}-{hi:2d}: {float((h * area_w * m).sum() / (area_w * m).sum()):+.3f} "
              f"W/m2 at the default 1000 m2/s")

    print(f"\nCLIM-17  ocean newly able to freeze seasonally, against "
          f"TFREEZE({s_declared:g} psu) = {tf_declared:.4f} K:")
    for s in args.salinity:
        tf = freezing_point_k(s)
        newly = float((area_w * (sea & (tmin > tf_declared) & (tmin <= tf))).sum())
        lo_k, hi_k = (sensitivity.planetary_albedo_to_kelvin(newly * c, alpha, cfg)
                      for c in ALBEDO_PER_ICE_AREA)
        lo_w, hi_w = (newly * c * incident for c in ALBEDO_PER_ICE_AREA)
        print(f"  S {s:5.1f} psu, TFREEZE {tf:.4f} (+{tf - tf_declared:.3f} K): "
              f"newly freezing {newly:.5f} of planet, "
              f"-{lo_w:.2f} to -{hi_w:.2f} W/m2 absorbed, {lo_k:+.2f} to {hi_k:+.2f} K")
    print("\nThis prints. The predictions and their wrongness bounds are in "
          "exoplasim/notes/forcing-bundle-predictions.md.")


if __name__ == "__main__":
    main()
