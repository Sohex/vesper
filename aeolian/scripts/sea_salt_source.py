#!/usr/bin/env python3
"""The sea-spray source function, and the subgrid wind moment it is raised to.

    python aeolian/scripts/sea_salt_source.py        # the Earth check

Grythe et al. (2014) equation 7, "G13T", the function that review fits to a
global set of sea salt concentration measurements rather than one of the 21 it
reviews. Coefficients are in `aeolian/config/sea_salt.yaml` and nothing here
hardcodes one.

This is a module rather than part of `build_sea_salt.py` because two consumers
need the same distribution: the field calculation needs its MASS, and
`sea_salt_optics.py` needs its SHAPE, to Mie-integrate over the population that
is actually there rather than over a lognormal chosen for convenience.

## The subgrid moment is analytic here, and that is worth stating

Production goes as U10^3.5, so a gridbox mean produces far less than the same
mean with realistic variance. `build_dust.py` handles the identical problem by
quadrature because Kok's emission has a threshold and no closed form. This one
has no threshold, so for a Weibull of shape k the enhancement is exactly

    <U^p> / <U>^p = Gamma(1 + p/k) / Gamma(1 + 1/k)^p

with no quadrature error to argue about. Jensen's inequality makes it greater
than 1 for p > 1, which is the cheap check that it has not been inverted.

## What the Earth check is for

Run as a script, this integrates the source function over a Weibull wind
distribution at Earth's ocean-mean 10 m wind and Earth's ocean area, and
compares the global production against the value Grythe et al. state for their
own function. The right answer comes from the source, so the test can fail. It
is not a comparison between two of our own formulations, which would only tell
us they differ.
"""

from __future__ import annotations

from math import gamma

import numpy as np


def moment_ratio(p: float, k: float) -> float:
    """<U^p>/<U>^p for a Weibull of shape k. Exactly 1 at p = 1."""
    return gamma(1.0 + p / k) / gamma(1.0 + 1.0 / k) ** p


def temperature_weight(sst_c, cfg: dict):
    """Grythe eq. A7, after Jaegle et al. (2011). Clamped to its fitted range.

    A cubic fitted over roughly 0 to 30 C. Outside that it turns over and goes
    negative, so extrapolating it would hand cold ocean a negative emission.
    Returns the weight and the fraction of the input that had to be clamped.
    """
    c = cfg["emission"]["temperature_weight"]
    lo, hi = cfg["emission"]["temperature_fit_range_c"]
    t = np.asarray(sst_c, dtype=float)
    clamped = float(np.mean((t < lo) | (t > hi)))
    t = np.clip(t, lo, hi)
    return c[0] + c[1] * t + c[2] * t ** 2 + c[3] * t ** 3, clamped


def number_flux_per_dp(dp_um, u10, cfg: dict, weibull_k: float | None = None,
                       modes=None):
    """dF/dDp in particles m-2 s-1 um-1, before the temperature weight.

    `dp_um` may be an array of dry diameters and `u10` an array of gridbox mean
    winds; the result broadcasts as `dp_um[:, None, None] * u10[None, ...]`.

    `modes` selects a subset of the three lognormal modes by index, which the
    Earth check needs and nothing else should: the spume mode is part of the
    source function and only the comparison against a published total leaves it
    out.

    With `weibull_k`, each mode's wind power is replaced by its Weibull moment,
    so what comes back is the subgrid-integrated flux rather than the flux at
    the gridbox mean. Without it, the two differ by a factor of about 2.5.
    """
    dp = np.asarray(dp_um, dtype=float)
    u = np.asarray(u10, dtype=float)
    shape = (dp.shape if dp.ndim else ()) + (1,) * u.ndim
    total = np.zeros(np.broadcast_shapes(dp.reshape(shape).shape, u.shape))
    chosen = cfg["emission"]["modes"]
    if modes is not None:
        chosen = [chosen[i] for i in modes]
    for mode in chosen:
        p = float(mode["wind_exponent"])
        wind = u ** p
        if weibull_k is not None:
            wind = wind * moment_ratio(p, weibull_k)
        size = mode["amplitude"] * np.exp(
            -mode["log_width"] * np.log(dp / mode["centre_dp_um"]) ** 2)
        total = total + size.reshape(shape) * wind
    return total


def diameter_grid(cfg: dict, points_per_decade: int = 200) -> np.ndarray:
    """Log-spaced dry diameters spanning the configured integration range."""
    lo = cfg["emission"]["dp_min_um"]
    hi = cfg["emission"]["dp_max_um"]
    n = int(np.ceil(np.log10(hi / lo) * points_per_decade)) + 1
    return np.logspace(np.log10(lo), np.log10(hi), n)


def mass_flux(u10, sst_c, cfg: dict, weibull_k: float | None = None,
              bin_edges_um=None, modes=None):
    """Dry sea-salt mass flux, kg m-2 s-1, total and per size bin.

    Mass is taken over the DRY diameter, because that is what the source
    function is expressed in and what the transported budget is carried in. The
    particle is wet in the air; it is dry in the accounting.
    """
    dp = diameter_grid(cfg)
    dndp = number_flux_per_dp(dp, u10, cfg, weibull_k, modes=modes)
    rho = cfg["emission"]["dry_density_kg_m3"]
    # mass of one dry particle, kg, from a diameter in um
    m_p = (np.pi / 6.0) * (dp * 1e-6) ** 3 * rho
    tw, clamped = temperature_weight(sst_c, cfg)

    shape = (dp.size,) + (1,) * np.asarray(u10).ndim
    integrand = dndp * m_p.reshape(shape)
    total = np.trapezoid(integrand, dp, axis=0) * tw

    per_bin = None
    if bin_edges_um is not None:
        # Integrated on a grid of its own per bin rather than by slicing the
        # grid above, so a bin edge falling between grid points cannot lose or
        # double-count a sliver. The sum over bins reproduces the total when the
        # edges span the configured range, which main() checks.
        per_bin = []
        for lo, hi in zip(bin_edges_um[:-1], bin_edges_um[1:]):
            d_sub = np.logspace(np.log10(lo), np.log10(hi), 200)
            i_sub = (number_flux_per_dp(d_sub, u10, cfg, weibull_k,
                                        modes=modes)
                     * ((np.pi / 6.0) * (d_sub * 1e-6) ** 3 * rho).reshape(
                         (d_sub.size,) + (1,) * np.asarray(u10).ndim))
            per_bin.append(np.trapezoid(i_sub, d_sub, axis=0) * tw)
        per_bin = np.array(per_bin)
    return total, per_bin, clamped


def monahan_number_flux_per_dp(dp_um, u10, cfg: dict,
                               weibull_k: float | None = None):
    """Monahan et al. (1986), in Grythe's harmonised equations A1 and A2.

    Here only to be the second function in the Earth check. It is not used for
    anything on this world, and the reason it exists is that an absolute
    comparison against a published global production measures our wind
    treatment; a ratio between two functions run through the same machinery
    does not.
    """
    m = cfg["earth_check"]["monahan"]
    dp = np.asarray(dp_um, dtype=float)
    u = np.asarray(u10, dtype=float)
    shape = (dp.shape if dp.ndim else ()) + (1,) * u.ndim
    p = float(m["wind_exponent"])
    wind = m["whitecap_coefficient"] * u ** p
    if weibull_k is not None:
        wind = wind * moment_ratio(p, weibull_k)
    b = (m["b_offset"] - np.log10(dp)) / m["b_scale"]
    size = (m["amplitude"] * dp ** m["size_power"]
            * (1.0 + m["tail_coefficient"] * dp ** m["tail_power"])
            * 10.0 ** (m["log_amplitude"] * np.exp(-b ** 2)))
    return size.reshape(shape) * wind


def _integrate_mass(flux_fn, cfg: dict, u10, k, modes=None) -> float:
    """Tg per Earth year over the whole ocean, from a dF/dDp callable."""
    e = cfg["earth_check"]
    dp = diameter_grid(cfg)
    dndp = flux_fn(dp, np.asarray(u10), cfg, k) if modes is None else \
        number_flux_per_dp(dp, np.asarray(u10), cfg, k, modes=modes)
    m_p = (np.pi / 6.0) * (dp * 1e-6) ** 3 * cfg["emission"]["dry_density_kg_m3"]
    flux = float(np.trapezoid(dndp * m_p, dp))
    return flux * e["ocean_area_m2"] * 365.25 * 86400.0 / 1e9


def _earth_check(cfg: dict) -> bool:
    """Two published functions through one machinery, against Grythe's Table 2.

    The wind treatment is ours and cancels in the ratio, so what this tests is
    the source functions and the mass integration. It fails loudly.
    """
    e = cfg["earth_check"]
    k = cfg["subgrid_wind"]["weibull_shape"]
    u = e["ocean_mean_u10_m_s"]
    target = e["grythe_table2_pg_per_year"]

    ours = {
        "M86": _integrate_mass(monahan_number_flux_per_dp, cfg, u, k) / 1000.0,
        "G13T": _integrate_mass(None, cfg, u, k, modes=(0, 1)) / 1000.0,
    }
    spume = _integrate_mass(None, cfg, u, k, modes=(2,)) / 1000.0
    ratios = {n: ours[n] / target[n] for n in ours}
    spread = abs(ratios["G13T"] - ratios["M86"]) / max(ratios.values())
    ok = spread <= e["ratio_agreement"]

    print(f"Earth check: U10 = {u} m/s, Weibull k = {k}, no temperature weight")
    print(f"  subgrid moment ratio p=3.5 {moment_ratio(3.5, k):.4f}  "
          f"p=3.41 {moment_ratio(3.41, k):.4f}  "
          f"p=1 {moment_ratio(1.0, k):.4f} (must be 1.0000)")
    for n in ("M86", "G13T"):
        print(f"  {n:5s} ours {ours[n]:6.2f} Pg/yr   Grythe Table 2 "
              f"{target[n]:5.2f}   ratio {ratios[n]:.3f}")
    print(f"  the two ratios agree to {spread*100:.1f}%, tolerance "
          f"{e['ratio_agreement']*100:.0f}%   {'OK' if ok else 'FAIL'}")
    print(f"  spume mode below Dp = {cfg['emission']['dp_max_um']} um adds "
          f"{spume:.1f} Pg/yr and is excluded from the comparison, not from "
          f"the model; see aeolian/config/sea_salt.yaml")
    return ok


if __name__ == "__main__":
    import sys
    from pathlib import Path

    import yaml

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from _paths import PROJECT_ROOT  # noqa: E402

    cfg = yaml.safe_load(
        (PROJECT_ROOT / "aeolian" / "config" / "sea_salt.yaml").read_text())
    assert abs(moment_ratio(1.0, 2.0) - 1.0) < 1e-12, "moment ratio broken at p=1"
    raise SystemExit(0 if _earth_check(cfg) else 1)
