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
`../patches/exoplasim-3.4.2-makestellarspec.patch` applied first: upstream adds
a `Path` to a `str` and calls `np.trapz`, which NumPy 2.0 removed.

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

from _paths import CONFIG, INPUTS, PATCHES

SSAP = "http://svo2.cab.inta-csic.es/theory/newov2/ssap.php"
MODEL = "bt-settl"

# BT-Settl CIFIST2011, log g = 4.5, [M/H] = 0, no alpha enhancement. The fids
# are SVO's record identifiers for those grid points; they are pinned so a
# silent server-side renumbering becomes a checksum failure rather than a
# different star. Verified against the header each file carries.
ENDPOINTS = {4900: 3697, 5000: 3850}

LOGG = 4.5
METALLICITY = 0.0

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
    target = CACHE / f"btsettl_{teff}_logg{LOGG}_m{METALLICITY}.txt"
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
        f"# logg = {LOGG} log(cm/s2)",
        f"# meta = {METALLICITY} ",
        "# alpha = 0 ",
        "#",
        "# column 1: WAVELENGTH (ANGSTROM), Wavelength in Angstrom",
        "# column 2: FLUX (ERG/CM2/S/A), Flux in erg/cm2/s/A",
    ]
    body = "\n".join(f"{w:15.3f}  {f:.6e} " for w, f in zip(wave, flux))
    path.write_text("\n".join(lines) + "\n" + body + "\n")


def ensure_patched() -> Path:
    """Confirm the makestellarspec fixes are applied to the vendored copy."""
    import exoplasim

    module = Path(exoplasim.__file__).resolve().parent / "makestellarspec.py"
    text = module.read_text()
    unpatched = 'Path(__file__).parent.resolve()+"/wvref.txt"'
    if unpatched in text or "np.trapz(" in text:
        patch = PATCHES / "exoplasim-3.4.2-makestellarspec.patch"
        raise SystemExit(
            f"{module} is unpatched and will fail on this NumPy.\n"
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
    teff = float(config["star"]["effective_temperature_k"])
    spectral_type = config["star"]["spectral_type"]

    wave, flux, meta = interpolate(teff, args.refresh)

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
        "star": {
            "spectral_type": spectral_type,
            "effective_temperature_k": teff,
            "log_g": LOGG,
            "metallicity": METALLICITY,
        },
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


if __name__ == "__main__":
    main()
