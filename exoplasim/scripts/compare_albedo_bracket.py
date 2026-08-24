#!/usr/bin/env python3
"""Compare the albedo-endmember runs and decide whether the bracket resolved.

The experiment asks whether this world has one stable climate or more than one.
Two land-surface endmembers, bare rock and vegetated, are integrated separately
at each flux; which run is which is read off `model.land_albedo_source` on its
own manifest, and the land-mean albedo each was handed is measured from its own
staged field and reported beside it. If they converge, the vegetation
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

  TARGET_K         290 to 293 K, the design range in force when this experiment
                   was registered and run. It was retired on 2026-08-16 for
                   having no derivation beyond its own persistence, and the mean
                   is now chosen by habitability per latitude band; see docs/src/pipeline/state.md,
                   section 5b. The constant stays as it was, because a threshold
                   is not rewritten after the results it judged have been seen.
                   Read the verdicts this script produced against that band, not
                   against the current one.

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

import climatology  # noqa: E402  from lib/, put on sys.path by _paths

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
                # The LEADING axis is time and its bins hold unequal numbers of
                # raw records, so it is weighted; any axis under it is a level
                # and averages plainly. CLIM-13.
                if a.ndim > 2:
                    a = climatology.annual_mean(a, np.asarray(ds["time"][:]))
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


def albedo_mode(run_dir: Path, manifest: dict, nlat: int,
                nlon: int) -> tuple[str, float]:
    """The endmember the run was CONFIGURED as, and the land-mean albedo it got.

    THE MODE IS THE DECLARATION, NOT A THRESHOLD OVER THE RESULT. This used to
    classify on `mean > 0.256`, a number placed between endmembers that were
    0.315 and 0.223 when it was written. The endmembers have moved: measured
    2026-08-24, `inputs/t21/albedo_report.json` gives 0.2483 for bare rock and
    0.1659 for vegetated. Both are BELOW 0.256, so the classifier returned
    "vegetated" for every run, `temps` is a dict keyed on the name, the two arms
    of the bracket collapsed into one entry, and the script reported "cannot
    bracket" on the experiment that decides whether the vegetation feedback is
    first-order. The separation between the endmembers survived -- 0.082 against
    the 0.092 the threshold was placed in -- so what failed was not resolution
    but PLACEMENT: a fixed cut point cannot track two means that both move.

    `model.land_albedo_source` is on every run manifest under `source_config`,
    it is the thing the run was asked for, and it moves with nothing. There is
    no fallback: a manifest without it is a manifest this script cannot read,
    and guessing from the field is what produced the collapse.

    The land-mean albedo is still measured from the run's own staged SRA and
    reported beside the mode, because it is what makes the two arms a bracket
    and `main` now checks that they actually differ.
    """
    from run_exoplasim import read_sra

    mode = (manifest.get("source_config", {}).get("model", {})
            .get("land_albedo_source"))
    if not mode:
        raise SystemExit(
            f"{run_dir.name} has no source_config.model.land_albedo_source on "
            "its manifest, so what endmember it is cannot be read off the run. "
            "Classifying it from the albedo field is what this script used to "
            "do and it is how the bracket collapsed.")

    # Named exactly, not globbed. ExoPlaSim writes these as N<nlat>_surf_<code>,
    # so the filename is fully determined by the grid the run already declares --
    # there is nothing to choose between, and sorting to pick [0] would silently
    # take another resolution's file if one were ever present.
    apath = run_dir / f"N{nlat:03d}_surf_0174.sra"
    if not apath.is_file():
        return str(mode), float("nan")
    field = read_sra(apath, 174, nlat, nlon)
    mpath = run_dir / f"N{nlat:03d}_surf_0172.sra"
    land = (read_sra(mpath, 172, nlat, nlon) > 0.5 if mpath.is_file()
            else np.ones_like(field, bool))
    w = gauss_weights(nlat)[:, None] * np.ones_like(field)
    mean = float((field[land] * w[land]).sum() / w[land].sum())
    return str(mode), mean


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
        mode, land_albedo = albedo_mode(d, man, int(mcfg["latitudes"]),
                                        int(mcfg["longitudes"]))
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
        # THE ARMS MUST DIFFER IN THE THING BEING BRACKETED. Two runs labelled
        # with different modes that were handed the same land albedo are one
        # experiment run twice, and their spread measures drift rather than
        # sensitivity. Checked on the field each run actually staged, so a
        # mislabelled manifest cannot pass as a bracket.
        albedos = {r["albedo_mode"]: r.get("land_mean_albedo") for r in at}
        finite = [v for v in albedos.values() if v == v]
        if len(finite) > 1 and max(finite) - min(finite) <= 0.0:
            print(f"  {flux:.2f}: the arms were handed the same land-mean "
                  f"albedo ({albedos}), so this is not a bracket")
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
