#!/usr/bin/env python
"""Difference a paired forcing arm against its control, orbit by orbit. A3.

WORLDBUILDING CONTEXT: Vesper is a fictional super-Earth. Every quantity below
is a global mean of the climate model that simulates it, and every orbit is a
modelled orbit.

WHAT A3 BUYS AND WHAT IT DOES NOT. A short A/B measures a FORCING in W/m2; a
converged run measures a RESPONSE in kelvin and is expensive. These arms branch
from one restart and differ by one namelist key, so the difference between them
is the term's, but a settling block is not a window: `lib/run_lengths.py` buys
the orbits it takes for a perturbation to decay, and a mean whose interval
covers needs a production span on top of that. Both numbers are reported and
they are not interchangeable.

THE FORCING IS READ EARLY AND THE RESPONSE LATE. The top-of-atmosphere
imbalance in the first orbits after branching is the term acting on a state that
has not yet moved, which is the forcing; the surface temperature separation at
the end is the response so far, and it is labelled by how much of the settling
block it has used rather than called converged.

WHAT WOULD MEAN THE ANSWER IS NOT A NUMBER. Every arm carries the paired
difference's own scatter across orbits. Where the separation does not exceed it,
the result is "below what this configuration resolves" and is reported that way
rather than as a small number, which is `docs/src/practice/failure-modes.md`
class 34.
"""
import argparse
import json
import pathlib
import re
import subprocess
import sys

import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[3]))
import lib.autocorrelation as ac  # noqa: E402
import lib.run_lengths as rl  # noqa: E402


def orbit_files(run_dir: pathlib.Path) -> list[pathlib.Path]:
    found = []
    for path in run_dir.glob("MOST.*.nc"):
        m = re.fullmatch(r"MOST\.(\d{5})\.nc", path.name)
        if m:
            found.append((int(m.group(1)), path))
    found.sort()
    return [p for _, p in found]


def gaussian_weights(nlat: int) -> np.ndarray:
    """Area weights for the model's own Gaussian latitudes."""
    nodes, weights = np.polynomial.legendre.leggauss(nlat)
    return weights[np.argsort(-nodes)]


def series(run_dir: pathlib.Path) -> dict:
    """Per-orbit global means of the metrics a forcing arm is read on."""
    files = orbit_files(run_dir)
    out = {"orbits": len(files), "tas": [], "ts": [], "toa_net": [],
           "sea_ice": [], "precip": [], "cloud": [], "toa_sw_out": [],
           "water_vapour": []}
    weights = None
    for path in files:
        with Dataset(path) as nc:
            if weights is None:
                weights = gaussian_weights(nc.variables["ts"].shape[-2])
                weights = weights / weights.sum()

            def gmean(name: str) -> float:
                # Time mean over the orbit's bins, then an area mean on the
                # model's own Gaussian weights. Longitude is uniform, so it
                # averages plainly; latitude does not and must not.
                field = np.asarray(nc.variables[name][:], dtype=float)
                zonal = field.mean(axis=-1)
                return float((zonal.mean(axis=0) * weights).sum())

            out["tas"].append(gmean("tas"))
            out["ts"].append(gmean("ts"))
            out["toa_net"].append(gmean("rst") + gmean("rlut"))
            out["sea_ice"].append(gmean("sic"))
            # IN MM/DAY, NOT THE MODEL'S OWN UNITS. `pr` is a rate in m/s and a
            # global mean of it is about 2.5e-8, which every fixed-width format
            # in this file rendered as 0.0000: a real ten per cent difference
            # between two arms read as an exact zero. The conversion is here
            # rather than at the print because the JSON carries these numbers too.
            out["precip"].append(gmean("pr") * 1000.0 * 86400.0)
            out["cloud"].append(gmean("clt"))
            # Negative-up in this model, so a MORE POSITIVE value is LESS
            # reflection. Checked against the file rather than assumed.
            out["toa_sw_out"].append(gmean("rsut"))
            out["water_vapour"].append(gmean("prw"))
    return out


def paired(arm: dict, control: dict, forcing_orbits: int, window_start: int = 0) -> dict:
    """Difference the two, and report the window mean the caller declared.

    `window_start` is the first orbit of the MEASUREMENT window: everything
    before it is the settling block and is excluded, because a pair branched
    from one restart starts at zero separation and averaging the approach into
    the answer understates it. The two numbers come from different places and
    must: the settling block from `lib/run_lengths.py`, which prices how long a
    perturbation takes to decay, and the window from `assess_convergence.py`'s
    own derivation of the shortest span its criteria can be tested on.
    """
    n = min(arm["orbits"], control["orbits"])
    result = {"orbits_compared": n, "forcing_window_orbits": min(forcing_orbits, n),
              "measurement_window_starts_at_orbit": window_start}
    for metric in ("tas", "ts", "toa_net", "sea_ice", "precip", "cloud",
                   "toa_sw_out", "water_vapour"):
        a = np.asarray(arm[metric][:n], dtype=float)
        c = np.asarray(control[metric][:n], dtype=float)
        diff = a - c
        k = result["forcing_window_orbits"]
        # The paired difference's own scatter is the instrument. Its standard
        # error carries the memory in the series, via lib/autocorrelation.py,
        # because consecutive orbits of this model are not independent samples.
        stats = ac.integrated_time(diff) if diff.size >= 4 else {"tau": float("nan")}
        tau = float(stats.get("tau", float("nan")))
        scatter = float(diff.std(ddof=1)) if diff.size > 1 else float("nan")
        # lib/autocorrelation.py owns this arithmetic and is the only place in
        # the tree that may: sigma * sqrt(tau / n), not sigma / sqrt(n), because
        # consecutive orbits carry memory and the raw count overstates how many
        # independent samples a window holds.
        sem = (ac.mean_standard_error(diff, tau)
               if diff.size > 1 and np.isfinite(tau) and tau > 0 else float("nan"))
        window = diff[window_start:] if window_start < diff.size else diff[:0]
        if window.size > 3:
            wtau = float(ac.integrated_time(window).get("tau", float("nan")))
            wsem = (ac.mean_standard_error(window, wtau)
                    if np.isfinite(wtau) and wtau > 0 else float("nan"))
            wmean = float(window.mean())
        else:
            wtau, wsem, wmean = float("nan"), float("nan"), float("nan")
        result[metric] = {
            "forcing_mean_first_%d_orbits" % k: float(diff[:k].mean()),
            "final_orbit_difference": float(diff[-1]),
            "mean_difference": float(diff.mean()),
            "window_mean": wmean,
            "window_standard_error": wsem,
            "window_tau_orbits": wtau,
            "window_orbits": int(window.size),
            "paired_scatter": scatter,
            "tau_orbits": tau,
            "standard_error": sem,
            "resolved": bool(np.isfinite(sem) and abs(diff.mean()) > 2.0 * sem),
            "trajectory": [float(x) for x in diff],
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control", type=pathlib.Path, required=True)
    parser.add_argument("--arm", type=pathlib.Path, action="append", required=True)
    parser.add_argument("--label", action="append", required=True)
    parser.add_argument("--forcing-orbits", type=int, default=3,
                        help="orbits at the start read as the forcing, before the "
                             "state has moved far from the shared restart")
    parser.add_argument("--window-start", type=int, default=0,
                        help="first orbit of the measurement window; everything before "
                             "it is the settling block and is excluded from the mean")
    parser.add_argument("--perturbation-k", type=float, required=True,
                        help="the predicted response, so the settling block this "
                             "arm would need is reported beside what it bought")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    if len(args.arm) != len(args.label):
        raise SystemExit("--arm and --label must come in pairs")

    control = series(args.control.resolve())
    settling = rl.settling_bracket(args.perturbation_k)
    arms = {}
    for arm_dir, label in zip(args.arm, args.label):
        s = series(arm_dir.resolve())
        arms[label] = {
            "run_directory": arm_dir.resolve().name,
            "orbits_on_disk": s["orbits"],
            "absolute_final_tas": s["tas"][-1] if s["tas"] else None,
            **paired(s, control, args.forcing_orbits, args.window_start),
        }

    payload = {
        "provenance": {
            "generator": "exoplasim/analysis/arms/compare_forcing_arms.py",
            "git_commit": subprocess.run(["git", "rev-parse", "HEAD"],
                                         capture_output=True, text=True).stdout.strip(),
            "python": sys.version.split()[0],
        },
        "instrument": {
            "control_run": args.control.resolve().name,
            "control_orbits": control["orbits"],
            "control_final_tas": control["tas"][-1] if control["tas"] else None,
            "settling_block_orbits_required": list(settling),
            "settling_block_bought": control["orbits"],
            "production_span_orbits_required":
                [rl.production_span_orbits(t) for t in rl.TAU_MEMORY_ORBITS_BRACKET],
            "note": ("The settling block is what these arms bought. A window whose "
                     "mean has an interval that covers needs the production span on "
                     "top of it, and none of these arms has one, so a separation is "
                     "reported as a trajectory and never as a converged response."),
        },
        "arms": arms,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    for label, a in arms.items():
        print(f"\n=== {label}  ({a['run_directory']}, {a['orbits_on_disk']} orbits)")
        for metric in ("tas", "toa_net", "toa_sw_out", "cloud", "water_vapour",
                       "sea_ice", "precip"):
            m = a[metric]
            key = [k for k in m if k.startswith("forcing_mean_first_")][0]
            print(f"  {metric:9s} forcing {m[key]:+.4f}  all {m['mean_difference']:+.4f}"
                  f"  WINDOW[{m['window_orbits']}] {m['window_mean']:+.4f}"
                  f" +/- {m['window_standard_error']:.4f}")


if __name__ == "__main__":
    main()
