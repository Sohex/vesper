"""Orbital period, derived the same way ExoPlaSim's run scripts derive it.

Duplicated literals for the year length drift the moment the flux changes, which
is exactly what happened: 189.6145 d is the 0.90-flux year and was still in the
hydrography scripts when the baseline moved to 0.96 and 180.655 d.
"""

from __future__ import annotations

import math

EARTH_SIDEREAL_YEAR_DAYS = 365.2568983
"""Earth's orbital period. Use for ratios between orbits, never for calendars."""

EARTH_CALENDAR_YEAR_DAYS = 365.2425
"""The Gregorian mean year. Use for annualising rates to "per Earth year".

Two constants because two different questions. They differ by 0.014 days, which
never matters numerically, but keeping one name for each stops a reader having to
work out which was meant. Both used to be scattered as literals across six
scripts in two components, which is exactly how the 189.6145-day year survived a
flux change."""


def orbital_year_days(config: dict, flux_ratio: float | None = None) -> float:
    """Orbital period in Earth days for a given stellar flux ratio."""
    star, orbit = config["star"], config["orbit"]
    if flux_ratio is None:
        flux_ratio = float(orbit["baseline_flux_earth"])
    semimajor_au = math.sqrt(float(star["luminosity_solar"]) / flux_ratio)
    return EARTH_SIDEREAL_YEAR_DAYS * math.sqrt(
        semimajor_au ** 3 / float(star["mass_solar"]))


def earth_years_per_orbit(config: dict, flux_ratio: float | None = None) -> float:
    """How many Earth years one orbit of this world lasts.

    The conversion every downstream consumer needs, in one place. LPJ-GUESS
    reports per simulation year, which is one orbit; multiply by this to compare
    against anything quoted per Earth year.
    """
    return orbital_year_days(config, flux_ratio) / EARTH_SIDEREAL_YEAR_DAYS
