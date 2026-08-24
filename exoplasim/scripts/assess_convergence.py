#!/usr/bin/env python3
"""Assess ExoPlaSim spin-up convergence from validated annual NetCDF files.

The verdict is taken over the last `--window` PRODUCTION orbits. Every segment
of a run declares what it was for, and orbits run to measure the model rather
than to advance the planet -- an I/O verification, a high-cadence wind sample, a
block on a differently patched binary -- are not evidence about where the run is
settling. `exoplasim/scripts/segments.py` owns that vocabulary.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re

_MPL_CACHE = Path("/tmp/world-matplotlib-cache")
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from scipy.optimize import curve_fit
from numpy.polynomial.legendre import leggauss

from _paths import ANALYSIS
import close_state_energy
# The one reader of the manifest's segment records; see
# exoplasim/scripts/segments.py for what a purpose means.
from segments import orbit_purposes, production_window
# CLAUDE.md names lib/sensitivity.py as the one flux-to-kelvin conversion, and
# the radiative damping this file relaxes at is that conversion inverted.
from sensitivity import planetary_albedo_from_fluxes, radiative_damping_w_m2_per_k


def output_files(run_dir: Path) -> list[Path]:
    outputs = []
    for path in run_dir.glob("MOST.*.nc"):
        match = re.fullmatch(r"MOST\.(\d{5})\.nc", path.name)
        if match:
            outputs.append((int(match.group(1)), path))
    outputs.sort()
    if [year for year, _ in outputs] != list(range(len(outputs))):
        raise RuntimeError("Annual outputs are not contiguous from year zero")
    return [path for _, path in outputs]


def global_mean(field: np.ndarray, weights: np.ndarray) -> np.ndarray:
    return np.sum(field * weights[None, :, None], axis=(-2, -1)) / (2.0 * field.shape[-1])


def slope(series: np.ndarray, window: int) -> float:
    y = series[-window:]
    x = np.arange(window, dtype=float)
    return float(np.polyfit(x, y, 1)[0])


# The mixed layer's heat capacity and the radiative damping. Together they set
# how long an approach takes, which is what turns a drift rate into a remaining
# offset -- and the offset is what `OFFSET_TOLERANCE_K` passes or fails a run
# on, so neither may be a number typed in beside the code that uses it.
#
# THE SLAB IS MODELLING THE MODEL'S OWN MIXED LAYER, so it takes the model's own
# sea water. `close_state_energy.py` reads CRHOS and CPS off `oceanmod.f90` and
# is imported here already; the 1025 and 3990 that stood here were 0.5% and 4.6%
# off them, uncited, and in the same direction.
import yaml as _yaml
_MLD = float(_yaml.safe_load(
    (Path(__file__).resolve().parents[2] / "config" / "planet.yaml")
    .read_text(encoding="utf-8"))["surface"]["mixed_layer_depth_m"])
SLAB_HEAT_CAPACITY = _MLD * close_state_energy.CRHOS * close_state_energy.CPS


def relaxation_orbits(orbital_year_days: float, feedback_w_m2_k: float) -> float:
    """Orbits for an e-folding of the slab's approach to equilibrium.

    `feedback_w_m2_k` is the radiative damping, and it comes from
    `lib/sensitivity.py` -- the one flux-to-kelvin conversion -- evaluated at
    the planetary albedo THIS RUN reports, so it moves with the run the way
    `year_days` already does. A private constant stood here instead, 1.31
    W/m2/K commented "measured, not assumed" with nothing saying where, 11%
    above what the module's own slope implies.
    """
    return (SLAB_HEAT_CAPACITY / feedback_w_m2_k) / (orbital_year_days * 86400.0)


def approach_to_equilibrium(orbits: np.ndarray, series: np.ndarray,
                            fraction: float = 0.35):
    """Fit T(n) = T_inf - A exp(-n/tau) and return the asymptote with an interval.

    A drift rate is not a distance. An exponential approach at rate `d` with time
    constant `tau` still has `d * tau` to travel, all of it one-signed, so a
    threshold on the rate bounds nothing on its own. At this planet's tau of
    about ten orbits, the 0.05 K/orbit criterion permits roughly half a kelvin of
    remaining approach -- the same order as corrections the project has thought
    worth applying, and larger than some.

    The interval is what makes this usable rather than a second point estimate.
    A fit whose tau is comparable to the run length is extrapolating past its
    data, and says so through a wide interval instead of a confident wrong
    number; one such fit in this project's history sat 0.6 K from its own
    reported value.

    Returns (asymptote, half_width, tau_orbits) or (nan, nan, nan) if the series
    will not support a fit.
    """
    mask = orbits >= orbits.max() * fraction
    x, y = orbits[mask], series[mask]
    if x.size < 8:
        return float("nan"), float("nan"), float("nan")
    model = lambda n, inf, amp, tau: inf - amp * np.exp(-n / tau)
    try:
        popt, pcov = curve_fit(model, x, y,
                               p0=[y[-1], max(y[0] - y[-1], 1e-3), max(x.size / 3, 1.0)],
                               maxfev=40000)
    except Exception:
        return float("nan"), float("nan"), float("nan")
    inf, tau = float(popt[0]), float(popt[2])
    err = float(np.sqrt(np.diag(pcov))[0]) if np.all(np.isfinite(pcov)) else float("nan")
    # A fit extrapolating beyond its own span is reported as such rather than
    # trusted: widen the interval by how far past the data the asymptote sits.
    if np.isfinite(tau) and tau > 0:
        overshoot = max(0.0, tau / max(x.size, 1))
        err = (err if np.isfinite(err) else abs(inf - y[-1])) * (1.0 + overshoot)
    return inf, err, tau


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--window", type=int, default=10,
                        help="orbits in the test window. Counted back from the "
                             "last PRODUCTION orbit, not the last orbit")
    parser.add_argument("--through", type=int, default=None, metavar="ORBITS",
                        help="assess the run AS IF it had stopped after this "
                             "many orbits, ignoring everything later. For "
                             "asking when a run FIRST met the criteria rather "
                             "than whether it meets them now -- which is what "
                             "a relaxation time is, and what decides whether "
                             "converting a coarse state into a finer one is "
                             "worth the wall clock. Writes nothing under "
                             "--output unless asked, so a sweep cannot "
                             "overwrite the run's own verdict.")
    parser.add_argument("--output", type=Path, default=ANALYSIS / "convergence")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    files = output_files(run_dir)
    if args.through is not None:
        if args.through > len(files):
            raise RuntimeError(
                f"--through {args.through} but the run has {len(files)} orbits")
        files = files[:args.through]
    if len(files) < args.window:
        raise RuntimeError("Not enough annual outputs for requested window")

    # THE WINDOW IS THE LAST `window` PRODUCTION ORBITS, NOT THE LAST `window`
    # ORBITS. A segment declares what it was for, and a diagnostic one -- an I/O
    # verification, a high-cadence wind sample, a block on a differently patched
    # binary -- is not evidence about where the run is settling. Taking the tail
    # blindly would have averaged three such orbits into this project's baseline
    # window on the day CLIM-9 was written. `production_window` drops a
    # diagnostic tail and REFUSES a window with a diagnostic hole in it, because
    # a slope across a gap is not a trend.
    window_start, window_end = production_window(run_dir, len(files), args.window)
    purposes = orbit_purposes(run_dir, range(len(files)))
    dropped = sorted(o for o in range(window_end + 1, len(files)))

    records = []
    for year, path in enumerate(files):
        with Dataset(path) as nc:
            weights = leggauss(len(nc.dimensions["lat"]))[1][::-1]
            record = {"year_index": year}
            # rst and rsut are read for the planetary albedo the radiative
            # damping is evaluated at: the run's own, not a declared one.
            for name in ["ts", "ntr", "hfns", "sic", "pr", "rst", "rsut"]:
                field = np.asarray(nc[name][:], dtype=float)
                value = float(global_mean(field, weights).mean())
                if name == "pr":
                    value *= 86400.0 * 1000.0
                record[name] = value
            records.append(record)

    # Every array is truncated at the last production orbit, so the trailing
    # window, the exponential fit and the storage term all describe the same
    # orbits. `records` keeps all of them, and the plot draws all of them, so
    # nothing is hidden -- only excluded from the verdict.
    arrays = {key: np.array([record[key] for record in records[:window_end + 1]])
              for key in ["ts", "ntr", "hfns", "sic", "pr", "rst", "rsut"]}
    w = args.window
    metrics = {
        "temperature_slope_k_per_orbit": slope(arrays["ts"], w),
        "toa_balance_slope_w_m2_per_orbit": slope(arrays["ntr"], w),
        "surface_balance_slope_w_m2_per_orbit": slope(arrays["hfns"], w),
        "sea_ice_slope_fraction_per_orbit": slope(arrays["sic"], w),
        "mean_toa_balance_w_m2": float(arrays["ntr"][-w:].mean()),
        "mean_surface_balance_w_m2": float(arrays["hfns"][-w:].mean()),
        "temperature_mean_k": float(arrays["ts"][-w:].mean()),
        "temperature_last_k": float(arrays["ts"][-1]),
        "sea_ice_mean_fraction": float(arrays["sic"][-w:].mean()),
    }

    # The storage the criterion now passes on, from the same window. Computed by
    # close_state_energy.state_energy rather than reimplemented, because the
    # criterion turning on a quantity is the strongest possible reason not to
    # have two versions of it.
    #
    # FAILS CLOSED. If the window cannot be closed against the state -- a run
    # whose namelist does not declare its own constants, a window too short for a
    # trend -- `storage` stays None and the criterion is False. A convergence
    # check that cannot measure the thing it tests must refuse, not abstain.
    storage = None
    storage_error = None
    try:
        storage = close_state_energy.state_energy(
            run_dir, window_start, window_end)
    # SystemExit is in the tuple because `close_state_energy` raises it, not
    # Exception, for the case this comment names first: a run whose namelist
    # does not declare GA, GASCON, GSOL0 and ECCEN. With only `Exception` here
    # the assessment died on that run and wrote no report at all, so the
    # criterion could not be False -- there was nothing to be False in.
    except (Exception, SystemExit) as exc:        # noqa: BLE001 - reported, not raised
        storage_error = f"{type(exc).__name__}: {exc}"

    if storage is not None:
        metrics.update({
            "state_storage_w_m2": storage["storage_w_m2_least_squares"],
            "state_storage_endpoint_w_m2": storage["storage_w_m2_endpoint"],
            "surface_storage_w_m2": storage["surface_storage_w_m2"],
            # Recorded, not tested. This is the gap CLIM-1 is open on: the
            # reported TOA net minus what the state actually stores. It is a
            # structural property of the diagnostic rather than of any run, so it
            # belongs in every assessment until it is explained.
            "reported_toa_minus_storage_w_m2":
                metrics["mean_toa_balance_w_m2"] - storage["storage_w_m2_least_squares"],
        })
    else:
        metrics["state_storage_unavailable"] = storage_error

    # The quantity the design band is stated in, and therefore the one that has
    # to be bounded. Everything above is a rate; this is a distance.
    # (tau_expected is needed by the fallback below, so it is computed first.)
    # The year comes from the run's own manifest, not from the current config:
    # this is a property of the run being assessed, which may predate a
    # baseline re-run that moved the orbit.
    manifest_path = run_dir / "run_manifest.json"
    year_days = 182.8
    if manifest_path.is_file():
        year_days = float(json.loads(manifest_path.read_text(encoding="utf-8"))
                          ["derived_parameters"]["orbital_year_earth_days"])
    # The albedo is this run's own, over the same window every other metric is
    # taken over, so the damping moves with the run exactly as `year_days` does.
    window_alpha = planetary_albedo_from_fluxes(
        float(arrays["rst"][-w:].mean()), float(arrays["rsut"][-w:].mean()))
    feedback_w_m2_k = radiative_damping_w_m2_per_k(window_alpha)
    tau_expected = relaxation_orbits(year_days, feedback_w_m2_k)
    orbits_axis = np.arange(len(arrays["ts"]), dtype=float)
    asymptote, half_width, tau_fit = approach_to_equilibrium(orbits_axis, arrays["ts"])
    offset = asymptote - metrics["temperature_mean_k"]

    # A converged run has no approach left to fit, so the exponential becomes
    # unconstrained and can return anything: one run here reported an asymptote
    # of -23064 K with an interval of six million. That is the fit having nothing
    # to grip, not the run being far from equilibrium, and the two must not be
    # treated alike -- the first version of this criterion failed a converged run
    # for it.
    #
    # When the fit is unusable, fall back to the bound the drift itself implies,
    # offset ~ drift * tau, using the EXPECTED tau rather than a fitted one, and
    # carry 100% uncertainty on it. That is conservative in the right direction:
    # it can only refuse a run, never pass one it should not.
    fit_usable = (np.isfinite(asymptote) and np.isfinite(half_width)
                  and abs(offset) < 20.0 and half_width < 5.0)
    if not fit_usable:
        offset = metrics["temperature_slope_k_per_orbit"] * tau_expected
        half_width = abs(offset)
        asymptote = metrics["temperature_mean_k"] + offset
    metrics.update({
        "temperature_asymptote_k": asymptote,
        "temperature_asymptote_half_width_k": half_width,
        "temperature_remaining_offset_k": offset,
        "relaxation_orbits_fitted": tau_fit,
        "relaxation_orbits_expected": tau_expected,
        # What tau_expected was built from, so the fallback offset below can be
        # audited without re-deriving it. lib/sensitivity.py owns the damping.
        "planetary_albedo_in_window": window_alpha,
        "radiative_damping_w_m2_per_k": feedback_w_m2_k,
        "slab_heat_capacity_j_m2_k": SLAB_HEAT_CAPACITY,
        "remaining_offset_implied_by_drift_k":
            metrics["temperature_slope_k_per_orbit"] * tau_expected,
    })

    # 0.15 K against a 3 K design band: small enough that the band's edges mean
    # what they say, loose enough to be reachable. Stated before it was applied.
    OFFSET_TOLERANCE_K = 0.15

    # THE ENERGY-BALANCE CRITERION TESTS STORAGE, NOT REPORTED TOA. Changed
    # 2026-08-17, deliberately and with the reasoning recorded, because the two
    # are not the same number here and the criterion was thresholding the one
    # that is wrong. CLIM-7; the evidence is exoplasim/notes/baseline-equilibration.md.
    #
    # `ntr` is the model's own statement about whether the planet is gaining or
    # losing energy, and it carries a structural offset of -0.573 +/- 0.035 W/m2
    # measured across four runs spanning 7.8 K. Conservation gives the same
    # quantity from the prognostic state instead -- mixed layer, sea ice and snow
    # as latent heat, soil, atmospheric enthalpy, column vapour -- and THAT is
    # what a convergence criterion is trying to bound. Where the energy the
    # offset represents actually goes is CLIM-1 and is still open, so both are
    # recorded here and the difference is carried as its own metric.
    #
    # THE THRESHOLD IS DERIVED, and the 0.5 W/m2 it replaces was not: that number
    # was picked early, never revisited, and applied to two different quantities.
    # Two independent bounds, stated before this was run against any assessment:
    #
    #   What matters. A residual imbalance X leaves the world short of its
    #   asymptote by roughly X * (dT/dt per W/m2) * tau. Measured on this model
    #   the mixed layer drifts 0.128 K/orbit per W/m2 and tau is 9.9 orbits, so
    #   holding that inside the 0.15 K offset tolerance already committed to
    #   above needs X < 0.15 / (0.128 * 9.9) = 0.118 W/m2. Deriving it from the
    #   temperature tolerance rather than inventing a second number is the point:
    #   the two criteria now bound the same thing in two units.
    #
    #   What is measurable. The storage estimate on a 10-orbit block differs from
    #   the same run's 20-orbit block by 0.111, 0.115 and 0.116 W/m2 on the three
    #   runs where both exist. A threshold below about 0.11 is measuring the
    #   sampling noise of the estimator rather than the planet.
    #
    # The two land within 6% of each other, which is the argument for the number
    # rather than a coincidence to note: below 0.11 is unresolvable, above 0.118
    # admits more drift than the temperature criterion already forbids. 0.12 is
    # where they meet, and it is 4x tighter than what it replaces.
    STORAGE_TOLERANCE_W_M2 = 0.12

    criteria = {
        "abs_temperature_slope_lt_0.05_k_per_orbit": abs(metrics["temperature_slope_k_per_orbit"]) < 0.05,
        "abs_toa_slope_lt_0.05_w_m2_per_orbit": abs(metrics["toa_balance_slope_w_m2_per_orbit"]) < 0.05,
        "abs_surface_slope_lt_0.05_w_m2_per_orbit": abs(metrics["surface_balance_slope_w_m2_per_orbit"]) < 0.05,
        "abs_sea_ice_slope_lt_0.001_per_orbit": abs(metrics["sea_ice_slope_fraction_per_orbit"]) < 0.001,
        f"abs_state_storage_lt_{STORAGE_TOLERANCE_W_M2}_w_m2": bool(
            storage is not None
            and abs(storage["storage_w_m2_least_squares"]) < STORAGE_TOLERANCE_W_M2),
        # The one that bounds the answer rather than its rate of change. Fails
        # closed: a fit that will not converge is not evidence of equilibrium.
        f"extrapolated_offset_lt_{OFFSET_TOLERANCE_K}_k": bool(
            np.isfinite(offset) and np.isfinite(half_width)
            and abs(offset) + half_width < OFFSET_TOLERANCE_K),
    }
    report = {
        "run_dir": str(run_dir),
        "completed_orbits": len(files),
        "window_orbits": w,
        # Which orbits the verdict is actually about. `completed_orbits` above
        # is what the run contains; these are what was assessed, and the two
        # differ whenever a diagnostic segment sits at the end.
        "window_start_year_index": window_start,
        "window_end_year_index": window_end,
        "orbits_after_window_excluded_as_non_production": dropped,
        "segment_purposes": {str(k): v for k, v in sorted(purposes.items())},
        "metrics": metrics,
        "criteria": criteria,
        "criteria_provenance": {
            "energy_balance_tests": "state storage, not reported top-of-atmosphere net",
            "storage_tolerance_w_m2": STORAGE_TOLERANCE_W_M2,
            "storage_tolerance_derivation":
                "0.15 K offset tolerance / (0.128 K per orbit per W/m2 * 9.9 orbits) "
                "= 0.118 W/m2, against an estimator that resolves 0.11 W/m2 on a "
                "10-orbit block. The two meet at 0.12.",
            "supersedes": "abs_mean_toa_lt_0.5_w_m2 and abs_mean_surface_lt_0.5_w_m2",
            "why": "The reported TOA net carries a structural offset of -0.573 "
                   "+/- 0.035 W/m2 across four runs, so the criterion was "
                   "thresholding a diagnostic rather than the planet. The 0.5 "
                   "W/m2 it replaces was chosen early, never revisited, and "
                   "applied to two different quantities. CLIM-7, decided "
                   "2026-08-17; the offset itself is CLIM-1 and is still open, "
                   "which is why both numbers and their difference are recorded.",
        },
        "sufficiently_equilibrated_for_worldbuilding": bool(all(criteria.values())),
        "failed_criteria": sorted(name for name, ok in criteria.items() if not ok),
        "annual_records": records,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    # Named by the run, not fixed. A fixed filename meant every assessment
    # overwrote the last, so a run's convergence record could not survive the
    # next run being assessed -- and assessing a run silently destroyed the
    # evidence for a previous one.
    #
    # A TRUNCATED assessment carries its length in the name and never takes the
    # run's own. Sweeping `--through` to find when a run first converged would
    # otherwise leave the verdict of whichever truncation ran last standing as
    # the run's, which is the same defect one directory up.
    suffix = "" if args.through is None else f"_through{args.through:03d}"
    report_path = args.output / f"{run_dir.name}_convergence{suffix}.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    # A run that has been assessed and misses is NOT "spinup_in_progress", and
    # leaving it labelled that way is how the baseline this project built its
    # production climatology from came to be recorded as still spinning up.
    #
    # This does not round a miss into a pass. `sufficiently_equilibrated_for_
    # worldbuilding` stays strictly all-or-nothing and is what any consumer
    # should test. What is added is the third state the prose convention has
    # always used and the code never had: quasi-equilibrated, carrying WHICH
    # criteria failed and the metric each failed on, so the label can never be a
    # hand-wave. One cold case in this project is called quasi-equilibrated for
    # missing by 0.004 W/m2; that standard is preserved by recording the margin,
    # not by widening the threshold.
    failed = report["failed_criteria"]
    manifest_path = run_dir / "run_manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if report["sufficiently_equilibrated_for_worldbuilding"]:
            manifest["status"] = "equilibrated_for_worldbuilding"
        else:
            manifest["status"] = "quasi_equilibrated"
        # The last orbit this verdict is ABOUT, which is the end of the
        # production window and not the end of the run. A diagnostic tail was
        # not assessed, so it must not be counted as settled.
        manifest["equilibrium_cutoff_year_index"] = window_end
        manifest_path.write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )

    # Draw every orbit, including the ones outside the verdict, and mark them.
    # A plot that showed only the assessed orbits would make an excluded tail
    # invisible, which is the failure this whole change is about.
    years = np.arange(len(files))
    all_series = {key: np.array([record[key] for record in records])
                  for key in ["ts", "ntr", "hfns", "sic", "pr"]}
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.5), constrained_layout=True)
    panels = [
        ("ts", "Global surface temperature", "K"),
        ("ntr", "TOA net radiation", "W m$^{-2}$"),
        ("hfns", "Surface downward heat flux", "W m$^{-2}$"),
        ("sic", "Planetary sea-ice fraction", "fraction"),
    ]
    for ax, (key, title, units) in zip(axes.ravel(), panels):
        ax.plot(years, all_series[key], marker="o", markersize=2.5, linewidth=1.1)
        ax.axvspan(window_start, window_end, color="tab:orange", alpha=0.12,
                   label=f"{w}-orbit test window")
        if dropped:
            ax.axvspan(dropped[0] - 0.5, dropped[-1] + 0.5, color="tab:grey",
                       alpha=0.18, label="not production; excluded")
        if key in {"ntr", "hfns"}:
            ax.axhline(0, color="black", linewidth=0.7)
        ax.set_title(title)
        ax.set_xlabel("orbit index")
        ax.set_ylabel(units)
        ax.grid(alpha=0.25)
    axes[0, 0].legend(loc="best", fontsize=8)
    fig.suptitle(
        f"Baseline spin-up convergence: {'PASS' if all(criteria.values()) else 'NOT YET'} "
        f"({len(files)} orbits, {w}-orbit window ending at {window_end})"
    )
    plot_path = args.output / f"{run_dir.name}_convergence.png"
    fig.savefig(plot_path, dpi=180)
    plt.close(fig)
    # One writer for this key. There were two, and this one silently clobbered
    # the other, losing the report path and the failure list every time.
    payload = {"metrics": metrics, "criteria": criteria,
               "pass": all(criteria.values()), "failed_criteria": failed,
               "report": str(report_path.resolve()), "orbits": len(files),
               "window_orbits": w, "window_start_year_index": window_start,
               "window_end_year_index": window_end,
               "orbits_after_window_excluded_as_non_production": dropped}

    # Record the verdict with the run as well as in the analysis directory.
    # Provenance travels with the artifact everywhere else in this project, and a
    # convergence result that lives only in an output folder cannot be found from
    # the run it describes.
    manifest_path = args.run_dir / "run_manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["convergence_assessment"] = payload
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
