"""The one place a refractive-index dataset NAME becomes n(lambda) and k(lambda).

`aeolian/config/dust.yaml` declares which optical dataset this world's dust is,
per spectral region, under `optics.indices`. Every consumer of dust optics
resolves that declaration through this module rather than loading a table of its
own, which is the whole content of DUST-12: the offline forcing computed
everything from OPAC while the aerofile the model reads was built from the
measured datasets, so the two priced the same dust at absorption optical depths
differing by a factor of 2.5.

The names are the `indices` labels in `analysis/dust_optics.json`, so a name
that appears in the config, in that file, and in a provenance record is the same
name throughout.

    from dust_indices import selection, indices
    sel = selection()                       # {'band1': ..., 'band2': ..., 'longwave': ...}
    n_of, k_of = indices(sel["band1"])      # callables of wavelength in um

## What each dataset covers, and why there are three names rather than one

  Di Biagio et al. (2019)   0.37-0.95 um   measured, 19 natural soils
  Rocha-Lima et al. (2018)  0.95-2.45 um   measured, but FIGURE ONLY, digitised
  OPAC (Hess et al. 1998)   0.25-40.0 um   continuous, and the only dataset here
                                           that reaches the thermal infrared

The measured pair is one curve in two pieces, so `Di Biagio 2019 measured` and
`Rocha-Lima <panel> fine` both return the SAME composite k -- Di Biagio below
0.95 um, the named Rocha-Lima panel above it, held flat past 2.45 um -- and
differ only in which panel is used above 0.95 um. The names are per band because
that is how the config selects them and how `dust_optics.json` reports them; the
band boundary at 0.75 um is ExoPlaSim's, not the datasets'.

`n` is 1.52 for the measured case, which is Di Biagio's own wavelength-
independent value; Rocha-Lima assumed 1.56 and did not retrieve it.

Above 40 um and below 0.25 um nothing here is valid and `indices` says so rather
than extrapolating: `np.interp` would hold the end value flat and silently
return a number.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Callable

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "exoplasim" / "data" / "dust"
OPAC_PATH = DATA / "opac_mineral_refractive_index.dat"
ROCHALIMA_PATH = DATA / "rochalima2018_fig10_digitized.csv"
DUST_CONFIG = ROOT / "aeolian" / "config" / "dust.yaml"

# Di Biagio et al. (2019) Table 4, population mean over 19 natural soils. The
# real part is wavelength-independent at 1.52 +/- 0.04 by the paper's own
# statement, so only k is spectral.
DIBIAGIO_N = 1.52
DIBIAGIO_LAM = np.array([0.370, 0.470, 0.520, 0.590, 0.660, 0.880, 0.950])
DIBIAGIO_K = np.array([0.0033, 0.0024, 0.0018, 0.0012, 0.0010, 0.0009, 0.0009])

# Where the measured composite stops being a measurement. Above this k is held
# at its last digitised value, which carries about 5% of the stellar flux.
MEASURED_FLAT_ABOVE_UM = 2.45

# name -> which Rocha-Lima panel supplies k above 0.95 um. `Di Biagio 2019
# measured` is a band-1 label and band 1 ends at 0.75 um, so the panel named for
# it is never reached in the use it exists for; it is named rather than left
# undefined so the dataset is a total function of wavelength like the others.
_MEASURED_PANEL = {
    "Di Biagio 2019 measured": "algeria",
    "Rocha-Lima Algeria fine": "algeria",
    "Rocha-Lima Mauritania fine": "mauritania",
}
OPAC_NAME = "OPAC"

NAMES = tuple(_MEASURED_PANEL) + (OPAC_NAME,)

# What each name is valid over, in micrometres. Used to refuse rather than
# extrapolate; the measured composite is flat-extended to 4 um deliberately,
# which is ExoPlaSim's own shortwave ceiling.
COVERAGE_UM = {name: (0.34, 4.0) for name in _MEASURED_PANEL}
COVERAGE_UM[OPAC_NAME] = (0.25, 40.0)


@lru_cache(maxsize=1)
def opac_table() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """OPAC mineral wavelength, n and k. k ships negative; the sign is flipped."""
    d = np.loadtxt(OPAC_PATH)
    lam, n, k = d[:, 0], d[:, 1], np.abs(d[:, 2])
    order = np.argsort(lam)
    return lam[order], n[order], k[order]


@lru_cache(maxsize=4)
def rochalima_table(panel: str, series: str = "fine_mie"):
    """Binned median k for one Rocha-Lima panel, on a coarse wavelength grid."""
    lam, k = [], []
    for line in ROCHALIMA_PATH.read_text(encoding="ascii").splitlines():
        if line.startswith("#") or line.startswith("panel"):
            continue
        p, s, w, kk = line.split(",")
        if p == panel and s == series:
            lam.append(float(w) / 1000.0)
            k.append(float(kk))
    lam, k = np.asarray(lam), np.asarray(k)
    grid = np.array([0.95, 1.05, 1.25, 1.45, 1.65, 1.85, 2.05, 2.25, 2.45])
    out = np.array([np.median(k[np.abs(lam - g) < 0.05])
                    if (np.abs(lam - g) < 0.05).any() else np.nan for g in grid])
    ok = np.isfinite(out)
    return grid[ok], out[ok]


def _measured_k(lam_um: float, panel: str) -> float:
    """Di Biagio below 0.95 um, the named Rocha-Lima panel above, flat past 2.45."""
    if lam_um < 0.95:
        return float(np.interp(lam_um, DIBIAGIO_LAM, DIBIAGIO_K))
    grid, k = rochalima_table(panel)
    return float(np.interp(lam_um, grid, k)) if lam_um <= grid[-1] else float(k[-1])


def indices(name: str) -> tuple[Callable[[float], float], Callable[[float], float]]:
    """(n_of_lambda, k_of_lambda) for a dataset name, wavelength in micrometres.

    Raises on an unknown name rather than falling back to a default: a typo in
    the config is then a stop, not a silent switch of optical dataset, which is
    the failure this module exists to make impossible.
    """
    if name == OPAC_NAME:
        lam_t, n_t, k_t = opac_table()
        return (lambda lam: float(np.interp(lam, lam_t, n_t)),
                lambda lam: float(np.interp(lam, lam_t, k_t)))
    if name in _MEASURED_PANEL:
        panel = _MEASURED_PANEL[name]
        return (lambda lam: DIBIAGIO_N,
                lambda lam: _measured_k(lam, panel))
    raise SystemExit(
        f"unknown refractive-index dataset {name!r}. Known: {', '.join(NAMES)}. "
        f"The name has to match an `indices` label in analysis/dust_optics.json.")


def check_coverage(name: str, lo_um: float, hi_um: float) -> None:
    """Raise if a dataset is being asked for wavelengths it does not have."""
    a, b = COVERAGE_UM[name]
    if lo_um < a - 1e-9 or hi_um > b + 1e-9:
        raise SystemExit(
            f"{name} covers {a}-{b} um and is being integrated over "
            f"{lo_um}-{hi_um} um. Interpolation would hold the end value flat "
            f"and return a number that looks like a measurement.")


def load_config(path: Path | None = None) -> dict:
    return yaml.safe_load((path or DUST_CONFIG).read_text(encoding="utf-8"))


def selection(dust_cfg: dict | None = None, end: str = "central") -> dict:
    """The declared dataset per spectral region, from `aeolian/config/dust.yaml`.

    `end` is `central` for the declared choice and `bracket` for the absorbing
    end that DUST-12 keeps as a bound. Both are declared in the config; nothing
    here has a default of its own.
    """
    cfg = dust_cfg if dust_cfg is not None else load_config()
    key = {"central": "indices", "bracket": "indices_bracket"}[end]
    try:
        sel = cfg["optics"][key]
    except KeyError:
        raise SystemExit(
            f"aeolian/config/dust.yaml has no optics.{key}. It is the file that "
            f"declares which refractive indices this world's dust has; see "
            f"DUST-12 and notes/dust.md.") from None
    missing = {"band1", "band2", "longwave"} - set(sel)
    if missing:
        raise SystemExit(f"optics.{key} is missing {sorted(missing)}")
    for region, name in sel.items():
        if name not in NAMES:
            raise SystemExit(
                f"optics.{key}.{region} names {name!r}, which is not a dataset "
                f"this repository has. Known: {', '.join(NAMES)}")
    return dict(sel)
