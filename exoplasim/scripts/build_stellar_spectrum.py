"""Build an ExoPlaSim stellar spectrum for this world's star from BT-Settl.

ExoPlaSim's shipped `k2.dat` is the star K2-18, an M2.5V at about 3450 K, not
the spectral type K2. It ships no K dwarf spectrum at all. See
`../notes/stellar-spectrum-audit.md` for how that was established and what it
cost. This script builds the spectrum that should have been used.

The source is the BT-Settl (CIFIST2011) grid, served by the SVO Theoretical
Spectra Server because phoenix.ens-lyon.fr is unreliable. The grid steps 100 K
in this range, so the declared 4965 K is interpolated between the 4900 K and
5000 K models rather than rounded to either.

Conversion to ExoPlaSim's two-file format is done by ExoPlaSim's own
`makestellarspec.convert`, so the output is format-compatible by construction
rather than by reimplementation. That function needs
`../patches/exoplasim-3.4.2-makestellarspec.patch` applied first, for two
reasons: upstream adds a `Path` to a `str` and calls `np.trapz`, which NumPy
2.0 removed, and upstream resamples with `np.interp`, which point-samples a
line-blanketed spectrum and reads high in the blue. `ensure_patched` refuses to
build without both, because an unpatched module produces a plausible file
rather than an error.

## Which parameters live where

`config/planet.yaml` is the source of truth for the star, so the spectral type,
the effective temperature and the METALLICITY are read from it. Metallicity
belongs there because nothing else determines it: it is a free choice worth
0.0114 in the band-1 share per 0.5 dex, which is the largest single lever on
this file. `GRID_METALLICITY` below is not a second declaration of it -- it is
what the pinned SVO record ids describe, and the build refuses if the config
asks for a metallicity those ids do not serve.

Surface gravity is NOT declared and must not be. It FOLLOWS from the mass,
luminosity and effective temperature the config already carries, so declaring
it would create a fourth copy of a derived quantity that could silently
disagree with the three it comes from. `derive_log_g` computes it and snaps to
the BT-Settl grid, and the build refuses if the snap lands anywhere but the
grid point the fids are pinned to.

## The check this build leaves behind

The blend is integrated at SOURCE resolution, under `radmod.f90:solarini`'s own
band definition, and the band-1 share and `zcross/z1` go into the provenance
record. `check_consistency.py` compares `lib/stellar.py`'s reading of the
shipped file against those, so the file is checked against what it is supposed
to represent rather than against other copies of itself. Nothing outside the
repository is needed at check time. See `notes/audits/stellar-spectrum-oracle.md`.

    python exoplasim/scripts/build_stellar_spectrum.py

Writes `<name>.dat`, `<name>_hr.dat` and `<name>_provenance.json` into
`inputs/stellarspectra/`. Downloads are cached; pass --refresh to refetch.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from _paths import CONFIG, INPUTS, PATCHES  # noqa: F401  (also puts lib/ on sys.path)

import stellar

SSAP = "http://svo2.cab.inta-csic.es/theory/newov2/ssap.php"
MODEL = "bt-settl"

# BT-Settl CIFIST2011, log g = 4.5, [M/H] = 0, no alpha enhancement. The fids
# are SVO's record identifiers for those grid points; they are pinned so a
# silent server-side renumbering becomes a checksum failure rather than a
# different star. Verified against the header each file carries.
ENDPOINTS = {4900: 3697, 5000: 3850}

# What the pinned fids above describe. These are NOT declarations of the star:
# `config/planet.yaml` declares the metallicity and the mass, luminosity and
# temperature that fix the gravity, and the two checks below refuse a config
# that has moved away from what these record ids serve.
GRID_LOGG = 4.5
GRID_METALLICITY = 0.0

# BT-Settl steps 0.5 dex in surface gravity in this range.
LOGG_GRID_STEP = 0.5

# IAU 2015 nominal solar values, for the gravity derivation only:
# GM_sun = 1.3271244e20 m3/s2 and R_sun = 6.957e8 m give g_sun = 274.20 m/s2,
# and 5772 K is the nominal effective temperature `radmod.f90:207` also uses.
SOLAR_LOG_G_CGS = 4.437968
SOLAR_EFFECTIVE_TEMPERATURE_K = 5772.0

# The name deliberately does not collide with ExoPlaSim's misleading `k2`.
OUTPUT_NAME = "k25v"

CACHE = Path(tempfile.gettempdir()) / "vesper-btsettl-cache"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def fetch(teff: int, fid: int, refresh: bool) -> Path:
    """Download one BT-Settl grid point, or reuse the cached copy."""
    CACHE.mkdir(parents=True, exist_ok=True)
    target = CACHE / f"btsettl_{teff}_logg{GRID_LOGG}_m{GRID_METALLICITY}.txt"
    if target.is_file() and not refresh:
        return target
    url = f"{SSAP}?model={MODEL}&fid={fid}&format=ascii"
    print(f"fetching {teff} K from {url}")
    with urllib.request.urlopen(url, timeout=600) as response:
        target.write_bytes(response.read())
    return target


def read_btsettl(path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    """Read SVO's ASCII BT-Settl export: Angstrom against erg/cm2/s/A."""
    header: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            break
        if "=" in line:
            # "# teff = 4900 K (value for ...)" -> {"teff": "4900"}
            key, _, rest = line[1:].partition("=")
            value = rest.split("(")[0].strip().split()
            header[key.strip()] = value[0] if value else ""
    data = np.loadtxt(path, comments="#")
    return data[:, 0], data[:, 1], header


def interpolate(teff: float, refresh: bool) -> tuple[np.ndarray, np.ndarray, dict]:
    """Log-linear blend of the two bracketing grid points onto the upper grid.

    Flux is interpolated in the log because a spectrum varies exponentially with
    temperature at fixed wavelength; a linear blend would bias the blue end low.
    """
    lo_t, hi_t = min(ENDPOINTS), max(ENDPOINTS)
    if not lo_t <= teff <= hi_t:
        raise SystemExit(f"{teff} K is outside the pinned bracket {lo_t}-{hi_t} K")

    sources = {}
    for t in (lo_t, hi_t):
        path = fetch(t, ENDPOINTS[t], refresh)
        wave, flux, header = read_btsettl(path)
        if int(float(header.get("teff", -1))) != t:
            raise SystemExit(f"{path} declares teff={header.get('teff')}, expected {t}")
        sources[t] = (wave, flux, path, header)

    wave_hi, flux_hi, path_hi, _ = sources[hi_t]
    wave_lo, flux_lo, path_lo, _ = sources[lo_t]

    # The two grid points do not share a wavelength grid exactly; put the cooler
    # one onto the hotter one's grid before blending.
    flux_lo_on_hi = np.interp(wave_hi, wave_lo, flux_lo)

    fraction = (teff - lo_t) / (hi_t - lo_t)
    floor = 1e-300
    log_blend = (1.0 - fraction) * np.log(np.maximum(flux_lo_on_hi, floor)) + (
        fraction
    ) * np.log(np.maximum(flux_hi, floor))
    flux = np.exp(log_blend)
    flux[(flux_lo_on_hi <= floor) | (flux_hi <= floor)] = 0.0

    meta = {
        "interpolation_fraction": fraction,
        "endpoints": {
            str(lo_t): {"fid": ENDPOINTS[lo_t], "sha256": sha256(path_lo)},
            str(hi_t): {"fid": ENDPOINTS[hi_t], "sha256": sha256(path_hi)},
        },
    }
    return wave_hi, flux, meta


def write_btsettl(path: Path, wave: np.ndarray, flux: np.ndarray, teff: float) -> None:
    """Re-emit in SVO's own layout so `makestellarspec.readspec` parses it.

    That reader keys off a leading '#' block and counts its lines, then drops
    everything below 0.1 Angstrom. It also discards the final two split entries,
    so a trailing newline costs one row near 999 um; that is 1e-9 of the flux.
    """
    lines = [
        "# BT-Settl",
        f"# teff = {teff} K (interpolated between grid points by "
        "build_stellar_spectrum.py)",
        f"# logg = {GRID_LOGG} log(cm/s2)",
        f"# meta = {GRID_METALLICITY} ",
        "# alpha = 0 ",
        "#",
        "# column 1: WAVELENGTH (ANGSTROM), Wavelength in Angstrom",
        "# column 2: FLUX (ERG/CM2/S/A), Flux in erg/cm2/s/A",
    ]
    body = "\n".join(f"{w:15.3f}  {f:.6e} " for w, f in zip(wave, flux))
    path.write_text("\n".join(lines) + "\n" + body + "\n")


def derive_log_g(config: dict) -> float:
    """Surface gravity in cgs dex, from what `config/planet.yaml` already says.

    R/Rsun = sqrt(L/Lsun) (Tsun/Teff)^2 from the Stefan-Boltzmann law, then
    log g = log g_sun + log(M/Msun) - 2 log(R/Rsun). Nothing here is a new
    parameter: mass, luminosity and effective temperature are all declared, and
    a declared gravity would be a fourth value free to disagree with them.
    """
    star = config["star"]
    teff = float(star["effective_temperature_k"])
    radius_solar = (float(star["luminosity_solar"]) ** 0.5
                    * (SOLAR_EFFECTIVE_TEMPERATURE_K / teff) ** 2)
    return (SOLAR_LOG_G_CGS + np.log10(float(star["mass_solar"]))
            - 2.0 * np.log10(radius_solar))


def star_parameters(config: dict) -> dict:
    """What the config asks for, checked against what the pinned fids serve.

    The record ids in `ENDPOINTS` name one point in a four-dimensional grid.
    Reading the metallicity from the config without checking it against them
    would let a config edit silently produce a spectrum for a different star,
    which is the failure the fids were pinned to prevent in the first place.
    """
    star = config["star"]
    metallicity = star.get("metallicity")
    if metallicity is None:
        raise SystemExit(
            "config/planet.yaml declares no star.metallicity. It is a free "
            "parameter worth 0.0114 in the band-1 share per 0.5 dex, and "
            "rule 2 puts the star in the config; add it there.")
    metallicity = float(metallicity)
    if metallicity != GRID_METALLICITY:
        raise SystemExit(
            f"config asks for [M/H] = {metallicity}, but the pinned SVO record "
            f"ids {sorted(ENDPOINTS.values())} serve [M/H] = {GRID_METALLICITY}. "
            "Re-pin ENDPOINTS to the fids for the metallicity you want; do not "
            "relax this check.")

    log_g = derive_log_g(config)
    snapped = round(log_g / LOGG_GRID_STEP) * LOGG_GRID_STEP
    if snapped != GRID_LOGG:
        raise SystemExit(
            f"the config's mass, luminosity and temperature give log g = "
            f"{log_g:.4f}, whose nearest BT-Settl grid point is {snapped}, not "
            f"the {GRID_LOGG} the pinned fids serve. Re-pin ENDPOINTS.")
    return {
        "spectral_type": star["spectral_type"],
        "effective_temperature_k": float(star["effective_temperature_k"]),
        "log_g": GRID_LOGG,
        "log_g_derived": float(log_g),
        "log_g_grid_step": LOGG_GRID_STEP,
        "metallicity": metallicity,
    }


def _segment_integral(wave_um: np.ndarray, flux: np.ndarray,
                      lo_um: float, hi_um: float, exponent: float = 0.0) -> float:
    """Trapezoidal integral of `flux / wave_um**exponent` from `lo_um` to `hi_um`.

    The limits are inserted as their own samples rather than snapped to the
    nearest source point, so the support is the band definition's and not the
    grid's.
    """
    inside = (wave_um > lo_um) & (wave_um < hi_um)
    ww = np.concatenate(([lo_um], wave_um[inside], [hi_um]))
    ff = np.concatenate(([float(np.interp(lo_um, wave_um, flux))], flux[inside],
                         [float(np.interp(hi_um, wave_um, flux))]))
    if exponent:
        ff = ff / ww ** exponent
    return float(np.trapezoid(ff, ww))


def source_resolution_integrals(wave_angstrom: np.ndarray,
                                flux: np.ndarray) -> dict:
    """`solarini`'s band quantities for the blend, at ITS OWN resolution.

    This is the number the shipped file has to reproduce, and the reason it can
    fail: the file is a resampling of this blend, so the blend integrated
    directly is a right answer rather than a second opinion. It is computed
    from the arrays this script blended, never from the converter's output, so
    a converter defect cannot cancel out of both sides.

    The band definition is `lib/stellar.py`'s constants, which quote
    `radmod.f90`. Everything below `minwavel` is discarded, the split is at
    0.75 um, and the support ends where the model's grid does. Absolute units
    cancel in both quantities, so the erg/cm2/s/A the source carries is fine.
    """
    wave_um = wave_angstrom * 1.0e-4
    lo = stellar.MIN_WAVELENGTH_NM * 1.0e-3
    split = stellar.BAND_SPLIT_UM
    top = stellar.MAX_WAVELENGTH_UM
    band1 = _segment_integral(wave_um, flux, lo, split)
    band2 = _segment_integral(wave_um, flux, split, top)
    zcross = (_segment_integral(wave_um, flux, lo, split, exponent=4.0)
              + _segment_integral(wave_um, flux, split, top, exponent=4.0))
    return {
        "note": (
            "The blend integrated at BT-Settl's own resolution under "
            "radmod.f90:solarini's band definition. lib/stellar.py reads the "
            "same quantities off the shipped 2048-point file, and "
            "scripts/check_consistency.py compares them against these within "
            "stellar.SOURCE_BAND1_TOLERANCE and "
            "stellar.SOURCE_CROSS_SECTION_TOLERANCE. Disagreement means the "
            "file does not represent the source it was resampled from."),
        "rows": int(wave_angstrom.size),
        "min_wavelength_nm": stellar.MIN_WAVELENGTH_NM,
        "band_split_um": split,
        "max_wavelength_um": top,
        "band1_fraction": band1 / (band1 + band2),
        "cross_section_ratio_um4": zcross / band1,
    }


def ensure_patched() -> Path:
    """Confirm the makestellarspec fixes are applied to the vendored copy.

    Two separate failures, and only the first announces itself. Without the
    NumPy-2 fixes `convert` raises. Without `_rebin_conserving` it runs and
    writes a file that is biased in the blue by a couple of parts in a
    thousand, which nothing downstream would notice; that is why the marker is
    checked here rather than left to the patch stack, and why a `.venv`
    reinstall must be followed by `exoplasim/scripts/rebuild_binaries.py`.
    """
    import exoplasim

    module = Path(exoplasim.__file__).resolve().parent / "makestellarspec.py"
    text = module.read_text()
    unpatched = 'Path(__file__).parent.resolve()+"/wvref.txt"'
    missing = []
    if unpatched in text or "np.trapz(" in text:
        missing.append("the NumPy 2 fixes, without which convert() raises")
    if "_rebin_conserving" not in text:
        missing.append("the flux-conserving rebin, without which convert() "
                       "point-samples and the band-1 share reads high")
    if missing:
        patch = PATCHES / "exoplasim-3.4.2-makestellarspec.patch"
        raise SystemExit(
            f"{module} is missing " + " and ".join(missing) + ".\n"
            f"Apply it with:\n"
            f"  patch --forward --strip=1 --directory={module.parent} < {patch}"
        )
    return module


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true", help="refetch the grid points")
    parser.add_argument("--name", default=OUTPUT_NAME)
    args = parser.parse_args()

    module = ensure_patched()
    config = yaml.safe_load(CONFIG.read_text())
    star = star_parameters(config)
    teff = star["effective_temperature_k"]
    spectral_type = star["spectral_type"]

    wave, flux, meta = interpolate(teff, args.refresh)
    source = source_resolution_integrals(wave, flux)

    out_dir = INPUTS / "stellarspectra"
    out_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        blended = Path(tmp) / "blended.txt"
        write_btsettl(blended, wave, flux, teff)

        # convert() writes <name>.dat and <name>_hr.dat into the cwd.
        from exoplasim.makestellarspec import convert

        cwd = Path.cwd()
        try:
            import os

            os.chdir(tmp)
            convert(str(blended), args.name, plot=False, numwavelengths=2048,
                    normalize=False)
        finally:
            os.chdir(cwd)

        # Before anything is copied out: does the file the converter just
        # wrote reproduce the blend it was resampled from? A failure here
        # leaves the previous spectrum in place rather than replacing it with
        # one nothing downstream would notice was wrong.
        hires = Path(tmp) / f"{args.name}_hr.dat"
        shipped_band1 = stellar.band_fractions(path=hires)[0]
        shipped_ratio = stellar.cross_section_ratio(path=hires)
        d_band1 = shipped_band1 - source["band1_fraction"]
        d_ratio = shipped_ratio / source["cross_section_ratio_um4"] - 1.0
        print(f"\nband 1    : {shipped_band1:.6f} against "
              f"{source['band1_fraction']:.6f} at source resolution "
              f"({d_band1:+.2e}, tolerance "
              f"{stellar.SOURCE_BAND1_TOLERANCE:.0e})")
        print(f"zcross/z1 : {shipped_ratio:.6f} against "
              f"{source['cross_section_ratio_um4']:.6f} "
              f"({d_ratio:+.2e} relative, tolerance "
              f"{stellar.SOURCE_CROSS_SECTION_TOLERANCE:.0e})")
        if (abs(d_band1) > stellar.SOURCE_BAND1_TOLERANCE
                or abs(d_ratio) > stellar.SOURCE_CROSS_SECTION_TOLERANCE):
            raise SystemExit(
                "the converted file does not reproduce the blend it was "
                "resampled from, so nothing was written. The resampler is the "
                "first thing to look at: see "
                "notes/audits/stellar-spectrum-oracle.md.")

        products = {}
        for suffix in ("", "_hr"):
            src = Path(tmp) / f"{args.name}{suffix}.dat"
            dst = out_dir / src.name
            dst.write_bytes(src.read_bytes())
            products[dst.name] = {
                "sha256": sha256(dst),
                "rows": sum(1 for _ in dst.read_text().splitlines()) - 1,
            }

    provenance = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "star": star,
        "source_resolution": source,
        "grid": {
            "name": "BT-Settl (CIFIST2011)",
            "reference": "Allard, Homeier, 2012",
            "server": "SVO Theoretical Spectra Server",
            "endpoint": SSAP,
            **meta,
        },
        "converter": {
            "tool": "exoplasim.makestellarspec.convert",
            "module_sha256": sha256(module),
            "patch": "patches/exoplasim-3.4.2-makestellarspec.patch",
            "numwavelengths": 2048,
            "normalize": False,
        },
        "products": products,
        "replaces": (
            "ExoPlaSim's k2.dat, which is the star K2-18 (M2.5V, ~3450 K), not "
            "the spectral type K2. See notes/stellar-spectrum-audit.md."
        ),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=CONFIG.parent,
        ).stdout.strip()
        or None,
    }
    report = out_dir / f"{args.name}_provenance.json"
    report.write_text(json.dumps(provenance, indent=2) + "\n")

    print(f"\nwrote {out_dir}/{args.name}.dat and {args.name}_hr.dat")
    print(f"      {report.name}")
    print(f"\n{spectral_type}, {teff} K, interpolated at fraction "
          f"{meta['interpolation_fraction']:.2f} between "
          f"{min(ENDPOINTS)} and {max(ENDPOINTS)} K")
    print(f"log g {star['log_g_derived']:.4f} derived, built on the "
          f"{star['log_g']} grid point; [M/H] {star['metallicity']}")


if __name__ == "__main__":
    main()
