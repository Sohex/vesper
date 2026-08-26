#!/usr/bin/env python3
"""Measure the eddy wind the hyperdiffusion rule divides by.

    python exoplasim/scripts/measure_eddy_wind.py <run_dir> [<run_dir> ...]
    python exoplasim/scripts/measure_eddy_wind.py <run_dir> --first 10 --last 24

Worldbuilding frame: a numerical diagnostic on the Vesper project's climate
model. Every quantity below is a property of that model's simulated atmosphere.

WHY THIS EXISTS. `config/planet.yaml`'s `model.hyperdiffusion.eddy_wind_m_s` is
the single free quantity in

    tau_vorticity = pi * radius / (NTRU * eddy_wind)

and that timescale divides into every entry of the `timescales_days` table, five
rungs by four fields, which reach the model as TDISSD, TDISSZ, TDISST and
TDISSQ. The value is therefore linear in every damping timescale the model runs
with, at every rung of the resolution ladder. It had been recorded as a
measurement on a run that no longer exists, so the derivation could not be
inspected, checked or re-run. This script is the instrument, so that the next
reader can re-run it rather than trust a sentence.

THE DEFINITION, stated exactly, because the note it supports records a 3.3x
spread across defensible readings of "the eddy wind":

  * INSTANTANEOUS fields, from the snapshot stream. The quantity wanted is an
    advective time, and what advects a feature at the truncation scale is the
    wind that is there at that moment. A time-mean field carries the STATIONARY
    eddies only and drops the transients, which is a different quantity and a
    smaller one; this script reports it separately so the two are never
    confused.
  * EDDY means departure from the ZONAL MEAN at the same instant, level and
    row. A uniform zonal flow Doppler-shifts a wave rather than straining it,
    so the zonal-mean jet is excluded by construction. The jet-inclusive total
    is reported alongside as the upper bracket on the definition.
  * The magnitude is the VECTOR one, `u'^2 + v'^2`, so the result is the RMS
    speed of the eddy wind and not of one component.
  * The average is AREA-weighted horizontally, by the Gaussian row weights the
    model's own quadrature uses, and MASS-weighted vertically by `dsigma`.
  * The run's value is the root of the MEAN VARIANCE over samples, not the mean
    of per-sample roots. Variance is what averages; a mean of roots is smaller
    and is not the RMS of anything.

THE VERTICAL WEIGHT IS THE MODEL'S, NOT THE FILE'S. The postprocessor writes
`levp` as midpoints interpolated from `lev`, which is not the model's `sigmah`.
`plasim.f90:2160-2166` sets `sigma(1) = 0.5*sigmah(1)` and
`sigma(k) = 0.5*(sigmah(k-1) + sigmah(k))`, so the half levels are recovered
from `lev` by that recursion and the recovery is checked against
`sigmah(NLEV) = 1`.

WHETHER THE INSTRUMENT CAN SEE THE EFFECT, which is CLAUDE.md's rule on
checking the instrument against the size of the effect. Two scatters are
reported and both are against the same value:

  * SAMPLING. Per-orbit means, their integrated autocorrelation time, and the
    standard error of the run mean through `lib/autocorrelation.py`. The
    criterion fixed before the run: the instrument resolves the value if that
    standard error is under 5 per cent of it, which is what makes the bracket a
    statement about the definition rather than about the sample.
  * TRUNCATION. The eddy variance is split by zonal wavenumber and the share
    sitting above 0.7 of the truncation is reported. One `eddy_wind` is applied
    at every rung, so a measurement at one rung is transferable only in so far
    as the variance is not piled against that rung's truncation. The criterion
    fixed before the run: under 5 per cent is a measurement within 5 per cent
    of a truncation-converged one.

WHAT IT WRITES. One JSON per run under `exoplasim/analysis/eddy-wind/`, holding
the numbers, the criteria, the sample, and the run and executable identity so
the measurement survives the run directory.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path

import netCDF4 as nc
import numpy as np

from _paths import ANALYSIS, PROJECT_ROOT, RUNS      # noqa: F401
import autocorrelation as ac
import gridding

# Fixed before any result was seen. Both are shares of the measured value.
SAMPLING_CRITERION = 0.05        # standard error of the run mean, fraction
TRUNCATION_CRITERION = 0.05      # eddy variance above 0.7 N, fraction
TRUNCATION_BAND = 0.7            # "near the truncation" starts here

_SNAP = re.compile(r"MOST_SNAP\.(\d+)\.nc$")
_MEAN = re.compile(r"MOST\.(\d+)\.nc$")


def half_levels(sigma: np.ndarray) -> np.ndarray:
    """The model's `sigmah`, recovered from the `sigma` the file writes.

    `plasim.f90:2165-2166`. The recursion is exact and its own check is that
    the last half level is the surface.
    """
    sigma = np.asarray(sigma, dtype=np.float64)
    sigmah = np.empty_like(sigma)
    sigmah[0] = 2.0 * sigma[0]
    for k in range(1, sigma.size):
        sigmah[k] = 2.0 * sigma[k] - sigmah[k - 1]
    # The stream is single precision and the recursion accumulates its rounding
    # with alternating sign, about 1e-6 over ten levels. The bar is loose against
    # that and tight against a `lev` belonging to a different vertical grid.
    if abs(sigmah[-1] - 1.0) > 1e-4:
        raise SystemExit(
            f"the half levels recovered from `lev` end at {sigmah[-1]:.6f} and not at "
            "the surface, so `lev` is not this model's sigma and the column weight "
            "would be a weight for a different vertical coordinate")
    return sigmah


def column_weight(sigma: np.ndarray) -> np.ndarray:
    """`dsigma`, normalised. `plasim.f90:2160-2161`."""
    sigmah = half_levels(sigma)
    dsigma = np.diff(np.concatenate([[0.0], sigmah]))
    return dsigma / dsigma.sum()


def row_weight(lat: np.ndarray) -> np.ndarray:
    """The Gaussian row weights, from `lib/gridding.py` and checked against the axis.

    `lib/gridding.py:require_gaussian_rows` is the shared guard and its bar is
    1e-6 degrees, which is the agreement of two float64 evaluations of the same
    expression. The model's output stream stores its latitude axis in SINGLE
    precision, so the same comparison in float64 is 3e-6 out on rounding alone
    and the shared guard cannot be used against this stream. The check here is
    the same test taken through the same rounding: the axis on disk must be the
    single-precision image of the constructed rows. That is tighter than the
    shared bar in relative terms and is 2e-6 of a T21 row, so it still refuses
    an axis that belongs to a different grid, which is what CLAUDE.md rule 3 is
    protecting.
    """
    lat = np.asarray(lat, dtype=np.float64)
    spec = gridding.gaussian_grid(int(lat.size))
    want = gridding.gaussian_latitudes(spec.nlat).astype(np.float32).astype(np.float64)
    off = float(np.abs(lat - want).max())
    if off > 1e-5:
        raise SystemExit(
            f"the snapshot stream's latitude axis differs from the constructed "
            f"Gaussian rows by up to {off:.3g} degrees. Do NOT reconcile them by "
            "matching nearest centres; construct both from one source. "
            "CLAUDE.md rule 3.")
    weights = np.abs(np.diff(spec.sin_edges))
    return weights / weights.sum()


def _weighted_variance(field_sq, w_lev, w_lat) -> np.ndarray:
    """Reduce `(time, lev, lat, lon)` squares to one number per time."""
    per_row = field_sq.mean(axis=3)                       # longitude is uniform
    per_level = np.tensordot(per_row, w_lat, axes=([2], [0]))
    return np.tensordot(per_level, w_lev, axes=([1], [0]))


def file_variances(path: Path) -> dict:
    """Eddy and total wind variance per record in one file, plus the zonal split."""
    with nc.Dataset(path) as ds:
        lat = np.asarray(ds.variables["lat"][:], dtype=np.float64)
        sigma = np.asarray(ds.variables["lev"][:], dtype=np.float64)
        u = np.asarray(ds.variables["ua"][:], dtype=np.float64)
        v = np.asarray(ds.variables["va"][:], dtype=np.float64)
    w_lev = column_weight(sigma)
    w_lat = row_weight(lat)

    up = u - u.mean(axis=3, keepdims=True)
    vp = v - v.mean(axis=3, keepdims=True)
    eddy = _weighted_variance(up ** 2 + vp ** 2, w_lev, w_lat)
    total = _weighted_variance(u ** 2 + v ** 2, w_lev, w_lat)
    # The top model level on its own, as the note's second definition.
    top = (((up[:, 0] ** 2 + vp[:, 0] ** 2).mean(axis=2)) * w_lat).sum(axis=1)

    # Zonal-wavenumber split of the same eddy variance. Parseval on a real
    # transform of `nlon` points: variance = sum over m >= 1 of the one-sided
    # power, so the shares below sum to one against `eddy` above.
    nlon = u.shape[3]
    power = np.abs(np.fft.rfft(up, axis=3)) ** 2 + np.abs(np.fft.rfft(vp, axis=3)) ** 2
    power = power * (2.0 / nlon ** 2)
    if nlon % 2 == 0:
        power[..., -1] *= 0.5                            # Nyquist is not doubled
    power[..., 0] = 0.0                                  # the zonal mean is not an eddy
    per_row = np.moveaxis(power, 3, 0)                   # (m, time, lev, lat)
    by_m = np.tensordot(np.tensordot(per_row, w_lat, axes=([3], [0])),
                        w_lev, axes=([2], [0]))          # (m, time)
    return {"eddy": eddy, "total": total, "top": top, "by_zonal_wavenumber": by_m}


def measure(run_dir: Path, first: int | None, last: int | None) -> dict:
    manifest_path = run_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}

    snaps = sorted((p for p in (run_dir / "snapshots").glob("MOST_SNAP.*.nc")),
                   key=lambda p: int(_SNAP.search(p.name).group(1)))
    if not snaps:
        raise SystemExit(f"{run_dir}: no snapshot stream, so no instantaneous field to "
                         "measure. The time-mean stream carries stationary eddies only.")
    indices = [int(_SNAP.search(p.name).group(1)) for p in snaps]
    keep = [(i, p) for i, p in zip(indices, snaps)
            if (first is None or i >= first) and (last is None or i <= last)]
    if not keep:
        raise SystemExit(f"{run_dir}: no orbits in [{first}, {last}]")

    per_orbit, totals, tops, spectra, samples = [], [], [], [], 0
    for _, path in keep:
        out = file_variances(path)
        per_orbit.append(float(out["eddy"].mean()))
        totals.append(float(out["total"].mean()))
        tops.append(float(out["top"].mean()))
        spectra.append(out["by_zonal_wavenumber"].mean(axis=1))
        samples += int(out["eddy"].size)

    per_orbit = np.asarray(per_orbit)
    spectrum = np.mean(np.stack(spectra), axis=0)

    # The run's value: the root of the mean variance.
    eddy_wind = math.sqrt(per_orbit.mean())

    # Sampling scatter, on the per-orbit VARIANCE series because that is what
    # averages; converted to the wind through the derivative of the root.
    tau = ac.integrated_time(per_orbit)
    se_var = ac.mean_standard_error(per_orbit, tau=tau["tau"])
    se_wind = se_var / (2.0 * eddy_wind)

    ntru = int(str(manifest.get("hyperdiffusion", {}).get("rung", "T21")).lstrip("Tt"))
    edge = int(math.ceil(TRUNCATION_BAND * ntru))
    total_var = float(spectrum.sum())
    near_truncation = float(spectrum[edge:].sum()) / total_var if total_var else float("nan")

    return {
        "run_id": manifest.get("run_id", run_dir.name),
        "run_directory": run_dir.name,
        "executable_sha256": (manifest.get("executable") or {}).get("sha256"),
        "config_sha256": manifest.get("config_sha256"),
        "source_build": manifest.get("source_build"),
        "rung": manifest.get("hyperdiffusion", {}).get("rung"),
        "eddy_wind_m_s_in_the_run": manifest.get("hyperdiffusion", {}).get("eddy_wind_m_s"),
        "orbits": [i for i, _ in keep],
        "snapshots": samples,
        "measured": {
            "eddy_wind_m_s": eddy_wind,
            "top_level_eddy_wind_m_s": math.sqrt(float(np.mean(tops))),
            "total_wind_including_jet_m_s": math.sqrt(float(np.mean(totals))),
        },
        "sampling": {
            "per_orbit_variance_m2_s2": per_orbit.tolist(),
            "integrated_autocorrelation_time_orbits": tau["tau"],
            "effective_sample_size_orbits": tau["effective_sample_size"],
            "tau_reliable": tau["reliable"],
            "standard_error_m_s": se_wind,
            "standard_error_fraction": se_wind / eddy_wind,
            "criterion_fraction": SAMPLING_CRITERION,
            "resolves": bool(se_wind / eddy_wind < SAMPLING_CRITERION),
        },
        "truncation": {
            "zonal_wavenumber_variance_m2_s2": spectrum.tolist(),
            "band_starts_at_zonal_wavenumber": edge,
            "share_above_0.7N": near_truncation,
            "criterion_fraction": TRUNCATION_CRITERION,
            "transferable": bool(near_truncation < TRUNCATION_CRITERION),
        },
        "measured_on": datetime.now(timezone.utc).isoformat(),
        "instrument": "exoplasim/scripts/measure_eddy_wind.py",
    }


def stationary_only(run_dir: Path, first: int | None, last: int | None) -> dict | None:
    """The same reduction on the TIME-MEAN stream: stationary eddies alone.

    Reported so that a value measured on monthly means is recognisable as the
    smaller, different quantity it is rather than mistaken for this one.
    """
    means = sorted(run_dir.glob("MOST.*.nc"),
                   key=lambda p: int(_MEAN.search(p.name).group(1)))
    keep = [p for p in means
            if (first is None or int(_MEAN.search(p.name).group(1)) >= first)
            and (last is None or int(_MEAN.search(p.name).group(1)) <= last)]
    if not keep:
        return None
    out = [file_variances(p) for p in keep]
    def rms(key):
        return math.sqrt(float(np.mean([float(o[key].mean()) for o in out])))
    return {"eddy_wind_m_s": rms("eddy"),
            "top_level_eddy_wind_m_s": rms("top"),
            "total_wind_including_jet_m_s": rms("total"),
            "orbits": len(keep)}


def report(result: dict, stationary: dict | None) -> None:
    m, s, t = result["measured"], result["sampling"], result["truncation"]
    print(f"\n{result['run_id']}  {result['rung']}  "
          f"{len(result['orbits'])} orbits, {result['snapshots']} snapshots")
    print(f"  eddy wind, instantaneous, column mass-weighted   {m['eddy_wind_m_s']:.3f} m/s")
    print(f"    +/- {s['standard_error_m_s']:.3f} "
          f"({s['standard_error_fraction']*100:.2f}% of it), "
          f"tau {s['integrated_autocorrelation_time_orbits']:.2f} orbits, "
          f"ESS {s['effective_sample_size_orbits']:.1f}")
    print(f"    resolves the value at {SAMPLING_CRITERION*100:.0f}%: "
          f"{'yes' if s['resolves'] else 'NO'}")
    print(f"  top model level only                             {m['top_level_eddy_wind_m_s']:.3f} m/s")
    print(f"  total including the zonal-mean jet               {m['total_wind_including_jet_m_s']:.3f} m/s")
    if stationary:
        print(f"  stationary eddies only, time-mean stream          "
              f"{stationary['eddy_wind_m_s']:.3f} m/s "
              f"({stationary['eddy_wind_m_s']/m['eddy_wind_m_s']:.2f} of the instantaneous)")
        print(f"    its top level                                  "
              f"{stationary['top_level_eddy_wind_m_s']:.3f} m/s")
        print(f"    its total including the jet                    "
              f"{stationary['total_wind_including_jet_m_s']:.3f} m/s")
    print(f"  eddy variance above zonal wavenumber {t['band_starts_at_zonal_wavenumber']}: "
          f"{t['share_above_0.7N']*100:.2f}% "
          f"(criterion {TRUNCATION_CRITERION*100:.0f}%: "
          f"{'transferable' if t['transferable'] else 'NOT transferable'})")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("run_dirs", nargs="+", type=Path)
    ap.add_argument("--first", type=int, default=None, help="first orbit index to use")
    ap.add_argument("--last", type=int, default=None, help="last orbit index to use")
    ap.add_argument("--out", type=Path, default=ANALYSIS / "eddy-wind",
                    help="where the JSON verdicts land")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    for run_dir in args.run_dirs:
        run_dir = run_dir.resolve()
        result = measure(run_dir, args.first, args.last)
        stationary = stationary_only(run_dir, args.first, args.last)
        result["stationary_eddies_only"] = stationary
        report(result, stationary)
        target = args.out / f"{result['run_directory']}.json"
        target.write_text(json.dumps(result, indent=2) + "\n")
        print(f"  written: {target}")


if __name__ == "__main__":
    main()
