#!/usr/bin/env python3
"""Compare the albedo-endmember runs and decide whether the bracket resolved.

The experiment asks whether this world has one stable climate or more than one.
Two land-surface endmembers, bare rock at 0.315 land-mean albedo and vegetated at
0.223, are integrated separately at each flux. If they converge, the vegetation
feedback is a correction and either branch can be taken forward. If they stay
apart, the land-surface assumption is a first-order term and the coupled loop
has to resolve it.

This measures *sensitivity*, not bistability. Albedo is a fixed input here, since
vegetation is not interactive, so the two cases are separately forced problems
with separately correct answers. Bistability would mean one set of boundary
conditions admitting two stable states, which needs the vegetation feedback
actually closed and cannot be read off a wide spread.

The thresholds below are fixed here rather than chosen after seeing the numbers.

  CONVERGED_K      the two endmembers at one flux are the same climate within
                   this much mean surface temperature. 2 K, comfortably above
                   the residual drift a 50-orbit spin-up leaves behind and well
                   below the 15 K their forcing difference could produce.

  TARGET_K         the design range this project has been aiming at since the
                   first sweep, unchanged: 290 to 293 K.

  A flux is MARGINAL when its warmer endmember falls short of the target range,
  because that means the habitable case lies at a higher flux and has not been
  bracketed yet. Evaluate it on equilibrated values: a run still carrying a
  negative TOA balance is reporting an upper bound on its own temperature.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

from netCDF4 import Dataset
import numpy as np
from numpy.polynomial.legendre import leggauss

from _paths import ANALYSIS, RUNS

CONVERGED_K = 2.0
TARGET_K = (290.0, 293.0)
MM_PER_DAY = 86400.0 * 1000.0


def gauss_weights(nlat: int) -> np.ndarray:
    return leggauss(nlat)[1][::-1]


def orbit_means(run_dir: Path, last: int = 5) -> dict:
    """Global means over the final `last` completed orbits."""
    files = sorted(run_dir.glob("MOST.*.nc"))
    if not files:
        raise RuntimeError(f"{run_dir} has no output")
    chosen = files[-last:]
    acc: dict[str, list[float]] = {}
    for f in chosen:
        with Dataset(f) as ds:
            w = gauss_weights(ds.dimensions["lat"].size)
            for name in ("ts", "pr", "evap", "ntr", "hfns", "sic", "snd", "glac", "alb"):
                if name not in ds.variables:
                    continue
                a = np.asarray(ds[name][:])
                while a.ndim > 2:
                    a = a.mean(axis=0)
                acc.setdefault(name, []).append(
                    float((a.mean(axis=1) * w).sum() / w.sum()))
    out = {k: float(np.mean(v)) for k, v in acc.items()}
    out["orbits_used"] = len(chosen)
    out["orbits_total"] = len(files)
    for k, scale in (("pr", MM_PER_DAY), ("evap", MM_PER_DAY)):
        if k in out:
            out[k] *= scale
    if "evap" in out:
        out["evap"] = -out["evap"]      # code 182 is negative upward
    return out


def albedo_mode(run_dir: Path, nlat: int, nlon: int) -> tuple[str, float]:
    """Identify the endmember from the albedo field the run actually staged.

    Each run directory holds the SRA it was given, so the mode is recoverable
    from the run itself rather than from `inputs/`, which later cases overwrite.
    """
    from run_exoplasim import read_sra
    # Named exactly, not globbed. ExoPlaSim writes these as N<nlat>_surf_<code>,
    # so the filename is fully determined by the grid the run already declares --
    # there is nothing to choose between, and sorting to pick [0] would silently
    # take another resolution's file if one were ever present.
    apath = run_dir / f"N{nlat:03d}_surf_0174.sra"
    if not apath.is_file():
        return "uniform", float("nan")
    field = read_sra(apath, 174, nlat, nlon)
    mpath = run_dir / f"N{nlat:03d}_surf_0172.sra"
    land = (read_sra(mpath, 172, nlat, nlon) > 0.5 if mpath.is_file()
            else np.ones_like(field, bool))
    w = gauss_weights(nlat)[:, None] * np.ones_like(field)
    mean = float((field[land] * w[land]).sum() / w[land].sum())
    # The two endmembers are 0.12 apart, so a midpoint split is unambiguous.
    name = "lithology" if mean > 0.256 else "vegetated"
    return name, mean


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="*", default=None,
                    help="explicit run directories; otherwise all t21 glacier runs")
    ap.add_argument("--pattern", default="t21*glac*")
    ap.add_argument("--last", type=int, default=5)
    args = ap.parse_args()

    dirs = [Path(r) for r in args.runs] if args.runs else sorted(
        p for p in RUNS.glob(args.pattern) if (p / "run_manifest.json").is_file())
    if not dirs:
        raise SystemExit(f"no runs matched {args.pattern}")

    rows = []
    for d in dirs:
        man = json.loads((d / "run_manifest.json").read_text(encoding="utf-8"))
        mcfg = man["source_config"]["model"]
        mode, land_albedo = albedo_mode(d, int(mcfg["latitudes"]), int(mcfg["longitudes"]))
        rows.append({
            "run": d.name,
            "flux": float(man["derived_parameters"]["stellar_flux_ratio_earth"]),
            "status": man.get("status"),
            "albedo_mode": mode,
            "land_mean_albedo": land_albedo,
            **orbit_means(d, args.last),
        })
    rows.sort(key=lambda r: (r["flux"], r["albedo_mode"]))

    print(f"{'flux':>5} {'mode':>10} {'orbits':>7} {'ts K':>8} {'pr':>7} "
          f"{'ntr':>7} {'sic':>7} {'glac':>7} {'alb':>6}")
    for r in rows:
        print(f"{r['flux']:5.2f} {r['albedo_mode']:>10} {r['orbits_total']:7d} "
              f"{r.get('ts', float('nan')):8.2f} {r.get('pr', float('nan')):7.3f} "
              f"{r.get('ntr', float('nan')):7.2f} {r.get('sic', float('nan')):7.4f} "
              f"{r.get('glac', float('nan')):7.4f} {r.get('alb', float('nan')):6.3f}")

    verdict = {"converged_threshold_k": CONVERGED_K, "target_k": list(TARGET_K),
               "by_flux": {}}
    print()
    for flux in sorted({r["flux"] for r in rows}):
        at = [r for r in rows if r["flux"] == flux]
        temps = {r["albedo_mode"]: r.get("ts") for r in at}
        if len(temps) < 2:
            print(f"  {flux:.2f}: only {list(temps)} present, cannot bracket")
            continue
        spread = max(temps.values()) - min(temps.values())
        warmest = max(temps.values())
        converged = spread <= CONVERGED_K
        marginal = warmest < TARGET_K[0]
        verdict["by_flux"][f"{flux:.2f}"] = {
            "temperatures_k": temps, "spread_k": spread,
            "converged": bool(converged), "warmest_k": warmest,
            "marginal": bool(marginal),
        }
        print(f"  {flux:.2f}: spread {spread:5.2f} K -> "
              f"{'converged' if converged else 'SENSITIVE'}; "
              f"warmest {warmest:.2f} K -> "
              f"{'MARGINAL, needs a warmer flux' if marginal else 'reaches the target'}")

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    out = ANALYSIS / "albedo_bracket.json"
    out.write_text(json.dumps({"runs": rows, "verdict": verdict}, indent=2) + "\n",
                   encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
