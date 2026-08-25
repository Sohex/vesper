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

## What the Earth check is for, and which function it is about

Run as a script, this integrates over a Weibull wind distribution at an ocean
mean 10 m wind and Earth's ocean area and compares against Grythe et al.'s
Table 2. It is about MONAHAN, not about equation 7, because Table 2 prices each
function over that row's own validity range and applies the temperature weight
only to the rows whose acronym ends in T. Equation 7's row carries a weight
applied cell by cell over a reanalysis, which one mean temperature cannot
reproduce, and Grythe table no weight-free total for it. Monahan's rows carry
none, and Monahan appears in Table 2 twice over two different size ranges.

That gives three gates, each with a right answer the project did not choose: the
absolute production against the range published implementations of Monahan span,
which is the range Grythe attribute to the wind treatment; the ratio of the two
Monahan rows, which carries no wind at all because that function factorises into
a whitecap term and a size term, so it tests the mass integration and the size
grid; and the sum over the configured size bins against the total from the same
call, which is an identity. `aeolian/config/sea_salt.yaml` carries every number
and the argument for each.

Both arms run through `mass_flux`, the entry point `build_sea_salt.py` uses to
make the field, so the integrator under test is the one that ships.
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
              bin_edges_um=None, modes=None, flux_fn=None):
    """Dry sea-salt mass flux, kg m-2 s-1, total and per size bin.

    Mass is taken over the DRY diameter, because that is what the source
    function is expressed in and what the transported budget is carried in. The
    particle is wet in the air; it is dry in the accounting.

    `flux_fn` substitutes another dF/dDp for this world's, and the ONLY caller
    that passes one is the Earth check: it is what lets the Monahan arm run
    through the integrator that ships rather than through a second trapezoid
    written beside it. It takes `(dp_um, u10, cfg, weibull_k)` and has no mode
    structure, so `flux_fn` and `modes` are mutually exclusive rather than the
    second being silently ignored.
    """
    if flux_fn is not None and modes is not None:
        raise ValueError(
            "a supplied flux function has no mode structure; pass `modes` only "
            "with the default source function")

    def dndp_of(d):
        if flux_fn is not None:
            return flux_fn(d, u10, cfg, weibull_k)
        return number_flux_per_dp(d, u10, cfg, weibull_k, modes=modes)

    dp = diameter_grid(cfg)
    dndp = dndp_of(dp)
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
        # double-count a sliver. The sum over bins is then the same integral as
        # `total` computed a second way whenever the edges span the configured
        # range, which is an IDENTITY and is what `_earth_check` gates on.
        per_bin = []
        for lo, hi in zip(bin_edges_um[:-1], bin_edges_um[1:]):
            d_sub = np.logspace(np.log10(lo), np.log10(hi), 200)
            i_sub = (dndp_of(d_sub)
                     * ((np.pi / 6.0) * (d_sub * 1e-6) ** 3 * rho).reshape(
                         (d_sub.size,) + (1,) * np.asarray(u10).ndim))
            per_bin.append(np.trapezoid(i_sub, d_sub, axis=0) * tw)
        per_bin = np.array(per_bin)
    return total, per_bin, clamped


def monahan_number_flux_per_dp(dp_um, u10, cfg: dict,
                               weibull_k: float | None = None):
    """Monahan et al. (1986), in Grythe's harmonised equations A1 and A2.

    Here only for the Earth check. It is not used for anything on this world,
    and the reason it exists is that Grythe's Table 2 prices this function
    twice, over two different size ranges and with no temperature weight on
    either row. That makes it the one function in the table this project can
    compare against absolutely and in shape at the same time.
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


def _pg_per_year(flux_kg_m2_s: float, cfg: dict) -> float:
    """Ocean-wide production in Pg per Earth year, from a mean areal flux."""
    return (flux_kg_m2_s * cfg["earth_check"]["ocean_area_m2"]
            * 365.25 * 86400.0 / 1e12)


def _arm(cfg: dict, u10, k, dp_um, modes=None, flux_fn=None):
    """One Table 2 row, THROUGH THE SHIPPED INTEGRATOR. Pg/yr, and the identity.

    `mass_flux` is the entry point `build_sea_salt.py` uses to make the field,
    so this exercises the code that runs rather than a trapezoid written beside
    it: the bin loop, the per-bin grid and the sliver handling are all in the
    path. The row's own size range is passed as a single bin, because that is
    how the integrator takes an interval that is not the configured one.

    The temperature weight is divided back out because the Monahan rows of
    Table 2 carry none. Dividing it out afterwards rather than switching it off
    inside `mass_flux` keeps the integrator whole, and at one temperature the
    weight is a scalar so the division is exact.

    Returns the row's production, and separately the sum over the CONFIGURED
    size bins of the same call, which is the same integral computed a second
    way and is the identity `_earth_check` gates on.
    """
    sst = cfg["earth_check"]["sst_mean_c"]
    tw, _ = temperature_weight(sst, cfg)
    # One bin spanning the row's own range: `total` is always the configured
    # 0.01-10 um integral, so the row's interval has to come through the bin
    # path, which is also the path the field calculation uses.
    _, row_bins, _ = mass_flux(u10, sst, cfg, k, bin_edges_um=list(dp_um),
                               modes=modes, flux_fn=flux_fn)
    row = row_bins[0]
    _, per_bin, _ = mass_flux(u10, sst, cfg, k,
                              bin_edges_um=cfg["size"]["bin_edges_um"],
                              modes=modes, flux_fn=flux_fn)
    total, _, _ = mass_flux(u10, sst, cfg, k, modes=modes, flux_fn=flux_fn)
    identity = abs(float(np.sum(per_bin)) / float(total) - 1.0)
    return _pg_per_year(float(row) / float(tw), cfg), identity


def _earth_check(cfg: dict) -> bool:
    """Monahan through the shipped integrator, against Grythe's Table 2.

    Three gates, and each has a right answer supplied by the source rather than
    by this project. The ABSOLUTE production of both Monahan rows must land
    inside the range published implementations of that function span, which is
    the range Grythe attribute to the wind treatment and is the only absolute
    statement this check's wind can support. The RATIO of the two rows must
    reproduce theirs, and that one carries no wind at all because Monahan
    factorises into a whitecap term and a size term, so it is the integrator and
    the size grid that are under test. And the sum over the configured size bins
    must reproduce the total from the same call, which is an identity.

    G13T is reported and not gated: its Table 2 row carries the equation A7
    temperature weight, applied cell by cell over a reanalysis, and this check
    has one temperature. `aeolian/config/sea_salt.yaml` argues it.
    """
    e = cfg["earth_check"]
    # EARTH's pooled ocean wind-speed shape, not this world's.
    # `subgrid_wind.weibull_shape` is Vesper's and is refitted at build time
    # from a high-cadence extract, so reading it here made an Earth verdict move
    # for a Vesper reason. The check is sensitive to it: at the sourced wind the
    # gated arms pass at 1.6 and 2.0 and fail at 2.6.
    # `aeolian/config/sea_salt.yaml` carries the source.
    k = e["weibull_shape"]
    u = e["ocean_mean_u10_m_s"]
    table = e["grythe_table2"]

    ours, identity = {}, {}
    for name, row in table.items():
        fn = monahan_number_flux_per_dp if name.startswith("M86") else None
        # Equation 7 is compared without its spume mode; the Monahan form has
        # no mode structure and takes none.
        modes = None if fn else (0, 1)
        ours[name], identity[name] = _arm(cfg, u, k, row["dp_um"],
                                          modes=modes, flux_fn=fn)

    band_lo, band_hi = e["absolute_pg_per_year_band"]
    gated = [n for n, row in table.items() if row["gated"]]
    absolute_ok = all(band_lo <= ours[n] <= band_hi for n in gated)

    theirs_ratio = table["M86E"]["pg_per_year"] / table["M86"]["pg_per_year"]
    our_ratio = ours["M86E"] / ours["M86"]
    ratio_miss = abs(our_ratio / theirs_ratio - 1.0)
    ratio_ok = ratio_miss <= e["size_range_ratio_tolerance"]

    worst_identity = max(identity.values())
    identity_ok = worst_identity <= e["bin_sum_identity_tolerance"]
    ok = absolute_ok and ratio_ok and identity_ok

    print(f"Earth check: U10 = {u} m/s, Weibull k = {k}, "
          f"temperature weight divided out")
    print(f"  subgrid moment ratio p=3.5 {moment_ratio(3.5, k):.4f}  "
          f"p=3.41 {moment_ratio(3.41, k):.4f}  "
          f"p=1 {moment_ratio(1.0, k):.4f} (must be 1.0000)")
    for name, row in table.items():
        lo, hi = row["dp_um"]
        tag = "gated" if row["gated"] else "REPORTED, not gated: its Table 2 " \
                                          "row carries the A7 weight"
        print(f"  {name:5s} Dp {lo:5.2f}-{hi:5.1f} um   ours {ours[name]:7.2f} "
              f"Pg/yr   Grythe {row['pg_per_year']:5.2f}   "
              f"ratio {ours[name] / row['pg_per_year']:6.3f}   {tag}")
    print(f"  ABSOLUTE, the gated arms inside {band_lo}-{band_hi} Pg/yr, the "
          f"range published implementations of Monahan span: "
          f"{'OK' if absolute_ok else 'FAIL'}")
    print(f"  SIZE RANGE, M86E/M86 ours {our_ratio:.4f} against Grythe's "
          f"{theirs_ratio:.4f}, off by {ratio_miss*100:.2f}%, tolerance "
          f"{e['size_range_ratio_tolerance']*100:.0f}%   "
          f"{'OK' if ratio_ok else 'FAIL'}")
    print(f"  IDENTITY, the configured bins summing to the total, worst "
          f"{worst_identity:.2e}, tolerance "
          f"{e['bin_sum_identity_tolerance']:.0e}   "
          f"{'OK' if identity_ok else 'FAIL'}")
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
