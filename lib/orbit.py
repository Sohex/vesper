"""Orbital period, derived the same way ExoPlaSim's run scripts derive it.

Duplicated literals for the year length drift the moment the flux changes, which
is exactly what happened: 189.6145 d is the 0.90-flux year and was still in the
hydrography scripts when the baseline moved to 0.96 and 180.655 d.
"""

from __future__ import annotations

import math

EARTH_SIDEREAL_YEAR_DAYS = 365.2568983


def orbital_year_days(config: dict, flux_ratio: float | None = None) -> float:
    """Orbital period in Earth days for a given stellar flux ratio."""
    star, orbit = config["star"], config["orbit"]
    if flux_ratio is None:
        flux_ratio = float(orbit["baseline_flux_earth"])
    semimajor_au = math.sqrt(float(star["luminosity_solar"]) / flux_ratio)
    return EARTH_SIDEREAL_YEAR_DAYS * math.sqrt(
        semimajor_au ** 3 / float(star["mass_solar"]))
