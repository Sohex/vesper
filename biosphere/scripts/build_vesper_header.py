"""Generate `vesper.h`, the planetary constants LPJ-GUESS compiles against.

Everything LPJ-GUESS assumes about Vesper's orbit, spin and star resolves to one
of these constants, and every one of them is derived here from
`config/planet.yaml` rather than written down. That matters because the year
length is a function of `orbit.baseline_flux_earth`: move the flux and the orbit
moves, and a hardcoded year silently desynchronises the model from its forcing.

This project has already been bitten by exactly that. `lib/orbit.py` exists
because 189.6145 d, the 0.90-flux year, was still sitting in the hydrography
scripts after the baseline moved to 0.96 and 180.655 d.

    python biosphere/scripts/build_vesper_header.py

Writes `biosphere/generated/vesper.h` plus a provenance JSON, and copies the
header into the LPJ-GUESS source tree. Rebuild afterwards:
VESPER_YEAR_LENGTH_DAYS sizes arrays, so it is a compile-time constant and a
stale binary is a silently wrong one, not a failing one.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from _paths import CONFIG, GENERATED, GUESS_SOURCE, PROJECT_ROOT, climatology_path

import orbit  # lib/orbit.py, the single source of truth for the year length

EARTH_SOLAR_CONSTANT_DEFAULT = 1361.0
SOLAR_EFFECTIVE_TEMPERATURE_K = 5772.0

# LPJ-GUESS's FRADPAR is 0.5 for Earth. It is a surface quantity, while what we
# can measure is top-of-atmosphere, so the star is applied as a ratio against the
# Sun rather than as an absolute.
EARTH_FRADPAR = 0.5

# Earth's photosystem window. Lehmer et al. 2021 put the optimum for a K2V at
# 675, 711 and 746 nm, so a biosphere evolved here would plausibly run to about
# 0.75 um. Widening this is a worldbuilding decision, not a correction: pass
# --par-window to make it, and re-register the productivity prediction.
PAR_WINDOW_UM = (0.40, 0.70)


def planck(wavelength_um: np.ndarray, temperature_k: float) -> np.ndarray:
    """Spectral radiance, in whatever units; only ratios are ever taken."""
    h, c, k_b = 6.62607015e-34, 2.99792458e8, 1.380649e-23
    wl = wavelength_um * 1e-6
    return (1.0 / wl**5) / (np.exp(h * c / (wl * k_b * temperature_k)) - 1.0)


def par_fraction(wavelength_um: np.ndarray, flux: np.ndarray,
                 window: tuple[float, float]) -> float:
    inside = (wavelength_um >= window[0]) & (wavelength_um <= window[1])
    return float(np.trapezoid(flux[inside], wavelength_um[inside])
                 / np.trapezoid(flux, wavelength_um))


def resolve_spectrum(config: dict) -> Path | None:
    """Find the configured spectrum the same way run_exoplasim.py does."""
    name = config.get("radiation", {}).get("stellar_spectrum")
    if not name:
        return None
    roots = [PROJECT_ROOT / "exoplasim" / "inputs" / "stellarspectra"]
    try:
        import exoplasim as exo

        roots.append(Path(exo.__file__).resolve().parent / "stellarspectra")
    except ImportError:
        pass
    for base in roots:
        candidate = base / f"{name}.dat"
        if candidate.is_file():
            return candidate
    raise SystemExit(f"stellar spectrum {name!r} not found under {roots}")


def derive_fradpar(config: dict, window: tuple[float, float]) -> tuple[float, dict]:
    """Scale Earth's FRADPAR by this star's PAR fraction against the Sun's."""
    spectrum = resolve_spectrum(config)
    if spectrum is None:
        raise SystemExit(
            "radiation.stellar_spectrum is unset, so the PAR fraction cannot be "
            "measured. Set it, or pass --fradpar explicitly and say why."
        )
    wavelength, flux = np.loadtxt(spectrum, skiprows=1, unpack=True)
    star = par_fraction(wavelength, flux, window)
    sun = par_fraction(wavelength, planck(wavelength, SOLAR_EFFECTIVE_TEMPERATURE_K),
                       window)
    return EARTH_FRADPAR * star / sun, {
        "spectrum": str(spectrum.relative_to(PROJECT_ROOT)),
        "spectrum_sha256": hashlib.sha256(spectrum.read_bytes()).hexdigest(),
        "par_window_um": list(window),
        "star_par_fraction": star,
        "solar_par_fraction_same_method": sun,
        "ratio": star / sun,
        "earth_fradpar": EARTH_FRADPAR,
    }


def month_lengths(year_length: int) -> list[int]:
    """Twelve months as even as the year divides, remainder on the last."""
    base = year_length // 12
    lengths = [base] * 12
    lengths[-1] += year_length - base * 12
    if min(lengths) < 3:
        raise SystemExit(f"year of {year_length} days gives unusable months")
    return lengths


def fit_solstice_offset(year_length: int, obliquity_deg: float,
                        climatology: Path) -> tuple[float, dict]:
    """Fit the declination phase against what the climate model actually did.

    `driver.cpp` models declination as -obliquity * cos(2pi(day + offset)/year).
    Rather than deriving the offset from orbital elements and hoping the
    conventions line up, it is fitted to the `zdec` the ExoPlaSim climatology
    reports. The residual is returned so a bad fit is visible rather than silent.

    This pins day 0 of the LPJ-GUESS year to the first climatology bin. The input
    module must emit bins in that same order.
    """
    import netCDF4 as nc

    with nc.Dataset(climatology) as data:
        zdec = np.asarray(data["zdec"][:], dtype=float)
    bins = len(zdec)
    day = (np.arange(bins) + 0.5) * year_length / bins
    offsets = np.arange(0.0, year_length, 0.05)
    predicted = -obliquity_deg * np.cos(
        2.0 * np.pi * (day[None, :] + offsets[:, None]) / year_length)
    residuals = np.sqrt(np.mean((predicted - zdec[None, :]) ** 2, axis=1))
    best = int(np.argmin(residuals))
    fitted = predicted[best]
    return float(offsets[best]), {
        "climatology": str(climatology.relative_to(PROJECT_ROOT)),
        "bins": bins,
        "rms_residual_deg": float(residuals[best]),
        "max_abs_residual_deg": float(np.max(np.abs(fitted - zdec))),
        "note": "day 0 of the simulation year is the first climatology bin",
    }


def render(constants: dict) -> str:
    months = ", ".join(str(n) for n in constants["month_lengths"])
    return f"""// GENERATED by biosphere/scripts/build_vesper_header.py. Do not edit.
//
// Every assumption LPJ-GUESS makes about Vesper's orbit, spin and star. Derived
// from config/planet.yaml, because the year length is a function of the stellar
// flux and moves whenever the flux does.
//
// generated {constants['generated']}
// config   sha256 {constants['config_sha256'][:16]}
// flux     {constants['flux_earth']} S-Earth
// orbit    {constants['orbital_year_earth_days']:.4f} Earth days

#ifndef LPJ_GUESS_VESPER_H
#define LPJ_GUESS_VESPER_H

/// Days in a simulation year.
/** The orbit is {constants['orbital_year_earth_days']:.4f} Earth days; {constants['year_length_days']} steps of 24 h is
 *  {constants['modelled_year_hours']:.1f} h against the true {constants['true_year_hours']:.1f}, an error of {constants['year_length_error_percent']:+.2f}%.
 *
 *  The timestep stays 24 h deliberately. Every per-day rate constant in
 *  LPJ-GUESS was calibrated against 24 h of absolute time, so stepping in real
 *  {constants['rotation_hours']:.1f} h Vesper days would put {constants['rotation_error_percent']:.1f}% into respiration,
 *  decomposition, snowmelt and phenology alike. Stepping in Earth days confines
 *  the error to daylength, where it has a known sign.
 *
 *  This sizes arrays, so it is a compile-time constant. Regenerate and rebuild
 *  after any change to orbit.baseline_flux_earth.
 */
const int VESPER_YEAR_LENGTH_DAYS = {constants['year_length_days']};

/// Twelve month lengths summing to VESPER_YEAR_LENGTH_DAYS.
#define VESPER_MONTH_LENGTHS {{{months}}}

/// Stellar constant at Vesper's orbit, W/m2. {constants['flux_earth']} x {constants['earth_solar_constant']}.
const double VESPER_STELLAR_CONSTANT = {constants['stellar_constant_w_m2']:.4f};

/// Orbital eccentricity. Earth 0.01675.
const double VESPER_ECCENTRICITY = {constants['eccentricity']};

/// Obliquity in degrees. Earth 23.4.
const double VESPER_OBLIQUITY_DEG = {constants['obliquity_deg']};

/// Phase of the declination cycle, in days before day 0 of the simulation year.
/** Earth's 10.5 puts the December solstice ten days before 1 January. This is
 *  fitted to the solar declination the ExoPlaSim climatology actually reports:
 *  rms residual {constants['solstice_fit']['rms_residual_deg']:.2f} deg over {constants['solstice_fit']['bins']} bins, max {constants['solstice_fit']['max_abs_residual_deg']:.2f} deg.
 *
 *  IT ASSUMES day 0 of the LPJ-GUESS year is the first climatology bin. The
 *  input module must emit bins in that order or the seasons run out of phase
 *  with the forcing.
 */
const double VESPER_SOLSTICE_OFFSET_DAYS = {constants['solstice_offset_days']};

/// Fraction of surface shortwave that is photosynthetically active.
/** Earth's {constants['fradpar_detail']['earth_fradpar']} scaled by this star's {constants['fradpar_detail']['par_window_um'][0]}-{constants['fradpar_detail']['par_window_um'][1]} um fraction against the
 *  Sun's measured the same way: {constants['fradpar_detail']['star_par_fraction']:.4f} / {constants['fradpar_detail']['solar_par_fraction_same_method']:.4f} = {constants['fradpar_detail']['ratio']:.4f}.
 *
 *  The window is Earth's photosystem transplanted unchanged, which is a
 *  deliberate conservative choice rather than a physical claim. See
 *  biosphere/notes/productivity-prediction.md.
 */
const double VESPER_FRADPAR = {constants['fradpar']:.6f};

#endif // LPJ_GUESS_VESPER_H
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--climatology", type=Path, default=None,
                        help="climatology NetCDF to fit the declination phase against")
    parser.add_argument("--par-window", type=float, nargs=2, default=PAR_WINDOW_UM,
                        metavar=("LO_UM", "HI_UM"),
                        help="photosystem window in microns (default Earth's 0.40 0.70)")
    parser.add_argument("--no-install", action="store_true",
                        help="write the header but do not copy it into the LPJ-GUESS tree")
    args = parser.parse_args()

    config = yaml.safe_load(CONFIG.read_text())
    flux = float(config["orbit"]["baseline_flux_earth"])
    earth_solar = float(config["orbit"].get("earth_solar_constant_w_m2",
                                            EARTH_SOLAR_CONSTANT_DEFAULT))
    orbital_days = orbit.orbital_year_days(config)
    year_length = int(round(orbital_days))
    rotation_hours = float(config["planet"]["rotation_hours"])

    climatology = args.climatology or climatology_path()
    if not climatology.is_file():
        raise SystemExit(
            f"{climatology} does not exist, so the declination phase cannot be "
            f"fitted. Point --climatology at a built climatology."
        )
    obliquity = float(config["planet"]["obliquity_degrees"])
    offset, solstice_fit = fit_solstice_offset(year_length, obliquity, climatology)
    fradpar, fradpar_detail = derive_fradpar(config, tuple(args.par_window))

    constants = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "flux_earth": flux,
        "earth_solar_constant": earth_solar,
        "orbital_year_earth_days": orbital_days,
        "year_length_days": year_length,
        "month_lengths": month_lengths(year_length),
        "modelled_year_hours": year_length * 24.0,
        "true_year_hours": orbital_days * 24.0,
        "year_length_error_percent": 100.0 * (year_length - orbital_days) / orbital_days,
        "rotation_hours": rotation_hours,
        "rotation_error_percent": 100.0 * (rotation_hours - 24.0) / 24.0,
        "stellar_constant_w_m2": flux * earth_solar,
        "eccentricity": float(config["planet"]["eccentricity"]),
        "obliquity_deg": obliquity,
        "solstice_offset_days": offset,
        "solstice_fit": solstice_fit,
        "fradpar": fradpar,
        "fradpar_detail": fradpar_detail,
        "source_build": config.get("source_build"),
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }

    GENERATED.mkdir(parents=True, exist_ok=True)
    header = GENERATED / "vesper.h"
    header.write_text(render(constants))
    constants["header_sha256"] = hashlib.sha256(header.read_bytes()).hexdigest()
    (GENERATED / "vesper_provenance.json").write_text(
        json.dumps(constants, indent=2) + "\n")

    installed = ""
    if not args.no_install:
        target = GUESS_SOURCE / "framework" / "vesper.h"
        if not target.parent.is_dir():
            raise SystemExit(f"{target.parent} does not exist; is LPJ-GUESS unpacked?")
        shutil.copyfile(header, target)
        installed = f"\ninstalled to {target}"

    print(f"orbit          {orbital_days:.4f} Earth days at {flux} S-Earth")
    print(f"year length    {year_length} days of 24 h "
          f"({constants['year_length_error_percent']:+.2f}%)")
    print(f"months         {constants['month_lengths']}")
    print(f"solstice phase {offset:.2f} d "
          f"(rms {solstice_fit['rms_residual_deg']:.2f} deg over "
          f"{solstice_fit['bins']} bins)")
    print(f"FRADPAR        {fradpar:.4f} "
          f"({args.par_window[0]}-{args.par_window[1]} um)")
    print(f"\nwrote {header.relative_to(PROJECT_ROOT)}{installed}")
    print("REBUILD LPJ-GUESS: the year length sizes arrays at compile time.")


if __name__ == "__main__":
    main()
