"""The flux-to-kelvin conversion, in one place, because it was in three.

Anything that turns a radiative quantity into a temperature reads this module.
Three sensitivities were in simultaneous use and they are not interchangeable:
one measured at T21 on a superseded terrain across a bracket that no longer
contains the baseline, one measured at T42 on this build, and one taken from a
stellar sweep that crosses the ice transition. They differ by a factor of 2.2,
and the error budget's ranking is computed with whichever one the caller picked.

Two things are settled here, and they are separate errors.

## 1. The slope

`SLOPE_K_PER_FLUX_RATIO` is dT/df at the baseline flux, where f is the stellar
constant in Earth units. It is a LOCAL slope, bracketed between two converged
points that span the baseline, and it must not be used outside the regime it was
measured in. `verify()` recomputes it from the run index and complains if the
runs it names have moved.

## 2. The denominator

A top-of-atmosphere forcing is an ABSORBED-flux change, so converting one to a
flux ratio divides by absorbed flux per unit flux ratio, not by incident:

    absorbed per unit f  =  S_earth / 4 * (1 - planetary_albedo)

Dividing by incident insolation instead understates every answer by
1/(1 - planetary_albedo). Dividing by the incident flux AT the baseline rather
than per unit f overstates it by 1/f, in the other direction and much smaller.
Both were being done at once, and the two do not cancel.

The same denominator makes an albedo change and a forcing the same quantity.
Absorbed shortwave is S/4 * (1 - alpha), so a change of `d` in PLANETARY albedo
is worth `-d / (1 - alpha)` in flux ratio. That identity is what
`planetary_albedo_to_kelvin` is; `test_identity()` checks the two routes against
each other, and they are the same arithmetic written twice, so a disagreement is
a coding error rather than a physical claim.

## What this module deliberately does not do

It does not convert a SURFACE albedo change into a planetary one. That
attenuation is a modelled quantity with a factor-of-two uncertainty and it
belongs with the budget that owns the assumption, not here.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG = PROJECT_ROOT / "config" / "planet.yaml"
RUN_INDEX = PROJECT_ROOT / "exoplasim" / "runs" / "INDEX.json"
ARCHIVED_RUNS = PROJECT_ROOT / "archive" / "runs"

# --------------------------------------------------------------------------
# The slope. Measured 2026-08-18, on the active build.
# --------------------------------------------------------------------------
#
# Two T42 runs on `precarve-craton`, same executable, same config hash, same
# geography hash, bracketing downward from the baseline flux:
#
#     f = 0.910   run_bfa3f5269660   fitted asymptote 281.968 K
#     f = 0.945   run_524fbed77a9a   fitted asymptote 289.030 K
#
# giving 201.8 K per unit flux ratio. Run means rather than asymptotes give
# 202.1, so the fit contributes nothing to the uncertainty.
#
# Corroborated on a DIFFERENT terrain, which is the check that could have
# failed: three T42 points on the superseded `precarve-zoned-g1281` at 0.9125,
# 0.945 and 0.968 give a chord of 201.3 across the pair that spans the baseline,
# and 203.4 from the derivative at 0.945 of a quadratic through all three. Five
# estimates, two terrains, span 201 to 204. The declared value is the middle of
# that and the spread is 1.5%, which is far inside the factor of two the error
# budget claims for itself.
#
# THIS REPLACES 150.2, which was measured at T21 on `precarve-zoned-g1281`
# between f = 0.95 and f = 1.00. That bracket does not contain the 0.945
# baseline, so it violated the rule its own comment stated. It is 26% low.
#
# THIS IS NOT the stellar sweep's sensitivity. The 0.85-to-0.95 sweep spans 21
# W/m2 absorbed and 33 K, which is 330 K per unit flux ratio: 1.6x this, because
# it crosses the ice transition and this does not. Do not use one for the other.
SLOPE_K_PER_FLUX_RATIO = 202.0
SLOPE_SPREAD_K_PER_FLUX_RATIO = (201.0, 204.0)

SLOPE_BRACKET_RUNS = {
    "cold": {"run_id": "run_bfa3f5269660", "flux_ratio": 0.910,
             "asymptote_k": 281.968, "status": "quasi_equilibrated"},
    "warm": {"run_id": "run_524fbed77a9a", "flux_ratio": 0.945,
             "asymptote_k": 289.030, "status": "equilibrated_for_worldbuilding"},
}
SLOPE_TOLERANCE_K_PER_FLUX_RATIO = 3.0

# The bracket sits on geography 3a17c498, one surface revision behind the
# geography the baseline climatology was run on. The surface change is worth
# +0.80 K at fixed flux and does not cross the ice transition -- sea ice moves
# 1.8% to 1.5% -- so it moves the intercept and not the slope. Stated because
# pairing the 0.910 run against the CURRENT baseline run instead would give 224
# K per unit flux ratio, and that number is a surface change wearing a slope's
# units.
SLOPE_GEOGRAPHY = "3a17c498"


def config(cfg: dict | None = None) -> dict:
    return cfg if cfg is not None else yaml.safe_load(
        CONFIG.read_text(encoding="utf-8"))


def incident_w_m2_per_flux_ratio(cfg: dict | None = None) -> float:
    """Global-mean insolation per unit of stellar flux ratio, S_earth / 4."""
    return float(config(cfg)["orbit"]["earth_solar_constant_w_m2"]) / 4.0


def _gaussian_weights(nlat: int, nlon: int) -> np.ndarray:
    from numpy.polynomial.legendre import leggauss
    return leggauss(nlat)[1][::-1][:, None] * np.ones((1, nlon))


def planetary_albedo(cfg: dict | None = None,
                     climatology: Path | None = None) -> tuple[float, dict]:
    """Measured, not assumed: reflected over incident at the top of atmosphere.

    `rst` is net downward shortwave at the top and `rsut` is the upward part, so
    their difference is the incident flux and no separate insolation figure is
    needed. Taking it from the climatology rather than declaring it means it
    moves when the surface does, which is the whole reason the budget's 0.266
    went stale.
    """
    from netCDF4 import Dataset

    cfg = config(cfg)
    path = climatology or (PROJECT_ROOT / str(cfg["baseline_climatology"]))
    with Dataset(path) as ds:
        w = _gaussian_weights(len(ds.dimensions["lat"]),
                              len(ds.dimensions["lon"]))
        def mean(name: str) -> float:
            v = np.asarray(ds.variables[name][:], dtype=float)
            if v.ndim == 3:
                v = v.mean(axis=0)
            return float((v * w).sum() / w.sum())
        net_down = mean("rst")
        upward = abs(mean("rsut"))
        run_id = getattr(ds, "vesper_run_id", None)
        flux = float(getattr(ds, "vesper_flux_ratio", float("nan")))
    incident = net_down + upward
    alpha = upward / incident
    return alpha, {
        "source": str(Path(path).relative_to(PROJECT_ROOT)),
        "run_id": run_id,
        "flux_ratio": flux,
        "incident_w_m2": round(incident, 3),
        "reflected_w_m2": round(upward, 3),
        "absorbed_w_m2": round(net_down, 3),
    }


def absorbed_w_m2_per_flux_ratio(alpha: float,
                                 cfg: dict | None = None) -> float:
    """The denominator. Absorbed shortwave per unit stellar flux ratio."""
    return incident_w_m2_per_flux_ratio(cfg) * (1.0 - alpha)


def kelvin_per_w_m2(alpha: float, cfg: dict | None = None,
                    slope: float | None = None) -> float:
    """Kelvin per W/m2 of top-of-atmosphere ABSORBED-flux forcing."""
    return (SLOPE_K_PER_FLUX_RATIO if slope is None else slope) \
        / absorbed_w_m2_per_flux_ratio(alpha, cfg)


def forcing_to_kelvin(w_m2: float, alpha: float, cfg: dict | None = None,
                      slope: float | None = None) -> float:
    """A top-of-atmosphere forcing in W/m2 to kelvin. Positive warms."""
    return w_m2 * kelvin_per_w_m2(alpha, cfg, slope)


def flux_ratio_to_kelvin(d_flux: float, slope: float | None = None) -> float:
    return d_flux * (SLOPE_K_PER_FLUX_RATIO if slope is None else slope)


def planetary_albedo_to_kelvin(d_alpha: float, alpha: float,
                               cfg: dict | None = None,
                               slope: float | None = None) -> float:
    """A change in PLANETARY albedo to kelvin. Brighter cools, so the sign flips."""
    return flux_ratio_to_kelvin(-d_alpha / (1.0 - alpha), slope)


# --------------------------------------------------------------------------
# Checks that can fail
# --------------------------------------------------------------------------

def test_identity(alpha: float = 0.3, cfg: dict | None = None) -> None:
    """The albedo route and the forcing route are one conversion, so they agree.

    Not a physical claim and not evidence about the slope: it is a definition,
    which is the only kind of thing worth asserting about arithmetic. It catches
    a sign flip or a dropped `1 - alpha` in either path, which is what it is for.
    The albedo is an argument with an arbitrary default because the identity holds
    for every value of it; using the measured one would suggest otherwise.
    """
    d_alpha = 0.01
    # The forcing that a planetary albedo change of d_alpha exerts on one unit of
    # flux ratio, which is the normalisation `slope` is stated in.
    forcing = -d_alpha * incident_w_m2_per_flux_ratio(cfg)
    a = planetary_albedo_to_kelvin(d_alpha, alpha, cfg)
    b = forcing_to_kelvin(forcing, alpha, cfg)
    if not math.isclose(a, b, rel_tol=1e-12):
        raise AssertionError(f"albedo route {a} != forcing route {b}")


def _index_entries() -> list[dict]:
    entries = []
    if RUN_INDEX.is_file():
        entries += json.loads(RUN_INDEX.read_text(encoding="utf-8"))["runs"]
    if ARCHIVED_RUNS.is_dir():
        for p in sorted(ARCHIVED_RUNS.glob("*/INDEX_ENTRY.json")):
            entries.append(json.loads(p.read_text(encoding="utf-8")))
    return entries


def verify() -> list[str]:
    """Recompute the declared slope from the runs it names. Empty means agreed.

    The point is that `SLOPE_K_PER_FLUX_RATIO` stops being a literal somebody has
    to remember to update. The runs are named, their temperatures are read from
    the index rather than from this file, and a run that is re-run, extended or
    re-assessed moves the recomputed value and trips this. The predecessor
    constant was called FALLBACK_SLOPE and nothing ever fell back to it or
    checked it against anything for two terrains.
    """
    by_id = {e.get("run_id"): e for e in _index_entries()}
    problems, points = [], {}
    for role, want in SLOPE_BRACKET_RUNS.items():
        entry = by_id.get(want["run_id"])
        if entry is None:
            problems.append(f"{role} run {want['run_id']} is not in any run index")
            continue
        flux = float(entry.get("physical", {}).get("flux_ratio", float("nan")))
        if not math.isclose(flux, want["flux_ratio"], abs_tol=1e-9):
            problems.append(f"{want['run_id']} is at flux {flux}, "
                            f"not the {want['flux_ratio']} recorded here")
        geo = entry.get("physical", {}).get("geography")
        if geo != SLOPE_GEOGRAPHY:
            problems.append(f"{want['run_id']} is on geography {geo}, "
                            f"not the {SLOPE_GEOGRAPHY} the bracket was measured on")
        metrics = entry.get("convergence_metrics") or {}
        asymptote = metrics.get("temperature_asymptote_k")
        if asymptote is None:
            problems.append(f"{want['run_id']} has no fitted asymptote")
            continue
        if abs(asymptote - want["asymptote_k"]) > 0.05:
            problems.append(f"{want['run_id']} asymptote is {asymptote:.3f} K, "
                            f"not the {want['asymptote_k']:.3f} recorded here")
        points[role] = (flux, asymptote)
    if len(points) == 2:
        (fc, tc), (fw, tw) = points["cold"], points["warm"]
        measured = (tw - tc) / (fw - fc)
        if abs(measured - SLOPE_K_PER_FLUX_RATIO) > SLOPE_TOLERANCE_K_PER_FLUX_RATIO:
            problems.append(
                f"the bracket now measures {measured:.1f} K per unit flux ratio "
                f"against the declared {SLOPE_K_PER_FLUX_RATIO}")
    return problems


def provenance(cfg: dict | None = None) -> dict:
    """The conversion and where every number in it came from, for an artifact."""
    cfg = config(cfg)
    alpha, alpha_src = planetary_albedo(cfg)
    absorbed = absorbed_w_m2_per_flux_ratio(alpha, cfg)
    return {
        "module": "lib/sensitivity.py",
        "slope_k_per_flux_ratio": SLOPE_K_PER_FLUX_RATIO,
        "slope_spread_k_per_flux_ratio": list(SLOPE_SPREAD_K_PER_FLUX_RATIO),
        "slope_bracket_runs": SLOPE_BRACKET_RUNS,
        "slope_geography": SLOPE_GEOGRAPHY,
        "slope_regime": "local to the baseline flux, on the low-ice branch. The "
                        "0.85-to-0.95 stellar sweep gives 330 K per unit flux "
                        "ratio because it crosses the ice transition; that value "
                        "is not interchangeable with this one.",
        "planetary_albedo": round(alpha, 5),
        "planetary_albedo_source": alpha_src,
        "incident_w_m2_per_flux_ratio": round(incident_w_m2_per_flux_ratio(cfg), 3),
        "absorbed_w_m2_per_flux_ratio": round(absorbed, 3),
        "kelvin_per_w_m2_absorbed": round(SLOPE_K_PER_FLUX_RATIO / absorbed, 4),
        "denominator_note": "A top-of-atmosphere forcing is an absorbed-flux "
                            "change, so it divides by absorbed flux per unit "
                            "flux ratio. Dividing by incident insolation "
                            "understates by 1/(1 - planetary_albedo).",
    }


def main() -> None:
    cfg = config()
    test_identity()
    p = provenance(cfg)
    alpha = p["planetary_albedo"]
    print(f"slope          {SLOPE_K_PER_FLUX_RATIO} K per unit flux ratio "
          f"({SLOPE_SPREAD_K_PER_FLUX_RATIO[0]}-{SLOPE_SPREAD_K_PER_FLUX_RATIO[1]})")
    print(f"planet albedo  {alpha:.5f}  from {p['planetary_albedo_source']['source']}")
    print(f"incident       {p['incident_w_m2_per_flux_ratio']} W/m2 per unit flux ratio")
    print(f"absorbed       {p['absorbed_w_m2_per_flux_ratio']} W/m2 per unit flux ratio")
    print(f"conversion     {p['kelvin_per_w_m2_absorbed']} K per W/m2 absorbed")
    problems = verify()
    print("\nverify: " + ("agrees with the run index"
                          if not problems else "; ".join(problems)))


if __name__ == "__main__":
    main()
