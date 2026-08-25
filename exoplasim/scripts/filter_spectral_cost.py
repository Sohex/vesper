#!/usr/bin/env python3
"""What does `filter_kappa` cost the resolved spectrum, and can this be seen?

    python exoplasim/scripts/filter_spectral_cost.py <run_dir> --first 0 --last 39
    python exoplasim/scripts/filter_spectral_cost.py --compare <run_a> <run_b>

Worldbuilding frame: a numerical diagnostic on the Vesper project's climate
model. Nothing here is about the simulated planet.

WHY THIS EXISTS. `model.filter_kappa` is 8.0 and was chosen on the trap
boundary of `exoplasim/notes/physics-filter-stability.md`: the longest stable
step at every rung measured. The physics filter is this model's small-scale
damping -- two to three orders of magnitude stronger than the derived
hyperdiffusion at every scale -- so kappa sets the amplitude of the whole
damping, and nothing had priced it against the spectrum it damps.
`spectral_tail.py` answers a different question: it PASSES OR FAILS one run
against a confinement floor. This one asks how much the reported number would
MOVE under a change in kappa, and compares that against the instrument's own
window-to-window scatter, which is CLAUDE.md's rule on checking the instrument
against the size of the effect.

THE RELAXATION MODEL, and the one check that makes it more than an assumption.
The filter multiplies each coefficient by `f(n) = exp(-kappa x^gamma)`, `x =
n/NTRU`, at BOTH transform directions every timestep, so its rate is
`r(x) = 2 kappa x^gamma / dt`. The flow's own cascade rate at the same scale is
`c(x) = x / tau_vorticity`. A spectrum held between a source that fills it at
`c` and a sink that empties it at `r` sits at

    spec(x) / law(x) = c / (c + r) = 1 / (1 + R x^(gamma-1)),  R = 2 kappa tau / dt

where `law` is the run's own inertial range extrapolated. THE CHECK: this makes
the ratio one half exactly where `R x^(gamma-1) = 1`, i.e. at
`x = (dt / (2 kappa tau))^(1/(gamma-1))` -- which is both `spectral_tail.py`'s
bite point (its EXCESS_FACTOR is 2) and the crossover
`scripts/check_consistency.py` already gates on. So the model reproduces two
independent things that were derived without it, and is used here only to
differentiate: what a change in kappa is worth in each reported unit.

WHAT IT SAYS ABOUT KAPPA, and it is the point of the script. The bite point
goes as `kappa^(-1/(gamma-1))`. At gamma 16 that is `kappa^(-1/15)`: the
confinement criterion is very nearly BLIND to kappa. The depth of the spectrum
AT the truncation goes as `1/(1+R)`, which for `R >> 1` is linear in `1/kappa`.
The two reported quantities therefore differ by a factor of `gamma-1` in their
leverage on kappa, and only one of them can carry a kappa decision.

THE CRITERIA, FIXED HERE BEFORE ANY ARM RAN.

  1. RESOLVABLE. A quantity resolves a change in kappa by a factor
     KAPPA_FACTOR only if the predicted change exceeds RESOLVE_SIGMA times the
     instrument's own scatter, measured as the standard deviation across
     disjoint windows of one settled run.
  2. QUANTISATION. The bite point is an integer wavenumber, so a predicted
     change below `1/NTRU` is below what the instrument can report at all,
     whatever its scatter.
  3. WINDOWS. At least MIN_WINDOWS disjoint windows, or the scatter is not
     measured and no verdict is issued.

WHAT IT DOES NOT ANSWER. The weakest kappa that still RUNS. That is the trap
boundary and it needs the stability grid, not a spectrum.
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _paths  # noqa: F401
from paths import rel  # noqa: E402  from lib/, put on sys.path by _paths
from spectral_tail import ke_spectrum, truncation  # noqa: E402  one reader, one convention

ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "exoplasim" / "analysis"

# -- criteria, fixed before results were seen -------------------------------
KAPPA_FACTOR = 2.0      # the change in kappa a sweep would plausibly take
RESOLVE_SIGMA = 3.0     # an effect must clear this many window sigmas
MIN_WINDOWS = 4         # fewer than this and the scatter is not measured
WINDOW_ORBITS = 2       # the window the gamma trade was taken on
EXCESS_FACTOR = 2.0     # spectral_tail.py's bite criterion, restated for the model
# WHERE THE DEPTH IS READ, and why not at the truncation. The first version of
# this script read the depth at `m = NTRU`, which is where kappa acts most
# directly -- and at every rung on disk that wavenumber sits AT floating-point
# roundoff, so the number it returned was the roundoff floor and not the filter.
# Two things put it there and neither is the filter's amplitude: a triangular
# truncation gives `m = NTRU` exactly one meridional mode against `NTRU` at
# `m = 1`, and the last few wavenumbers of every rung on disk are dead. So the
# depth is read at a fraction of the truncation that is still resolved, and the
# reading is REFUSED unless the spectrum there stands clear of the measured
# roundoff floor by FLOOR_MARGIN. The fraction is 0.8 because the confinement
# floor is 0.6 and the fit band ends at a third: 0.8 is inside the damped band
# at every rung and outside the mode-count collapse at all of them.
DEPTH_FRACTION = 0.8
FLOOR_MARGIN = 100.0    # the spectrum must stand this far above roundoff to be read


def namelist_keys(run_dir: Path) -> dict[str, str]:
    path = run_dir / "plasim_namelist"
    if not path.is_file():
        raise SystemExit(f"no plasim_namelist in {run_dir}: the filter this run "
                         f"used is not recoverable and must not be assumed")
    out = {}
    for line in path.read_text(encoding="latin-1").splitlines():
        m = re.match(r"\s*([A-Za-z_]\w*)\s*=\s*(.*?)\s*$", line)
        if m:
            out[m.group(1).upper()] = m.group(2)
    return out


def rung_of(run_dir: Path) -> str:
    manifest = run_dir / "run_manifest.json"
    if not manifest.is_file():
        raise SystemExit(f"no run_manifest.json in {run_dir}")
    return str(json.loads(manifest.read_text())["physical"]["resolution"]).upper()


def cascade_time_s(rung: str) -> float:
    """tau_vorticity for this rung, from config/planet.yaml and nowhere else."""
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    hd = cfg["model"]["hyperdiffusion"]["timescales_days"]
    if rung not in hd:
        raise SystemExit(f"config/planet.yaml names no cascade time for {rung}")
    return float(hd[rung]["vorticity"]) * 86400.0


def measure(full: np.ndarray, n: int) -> dict:
    """The reported quantities, on one spectrum, in spectral_tail's convention.

    `full` runs out to the FFT's Nyquist, which the model represents nothing of
    above `n`. That dead band is not waste here: it MEASURES the roundoff floor,
    which is what tells the depth reading whether it is looking at the filter or
    at nothing.
    """
    floor = float(np.median(full[n + 1:])) if len(full) > n + 1 else 0.0
    spec = full[:n + 1]
    m = np.arange(len(spec))
    lo, hi = max(2, n // 8), max(3, n // 3)
    band = (m >= lo) & (m <= hi) & (spec > 0)
    fit = np.polyfit(np.log(m[band]), np.log(spec[band]), 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        law = np.exp(np.polyval(fit, np.log(np.where(m > 0, m, 1))))
        ratio = np.where(law > 0, spec / law, np.nan)
    upper = (m > hi) & (m <= n)
    bitten = [int(k) for k in m[upper] if ratio[k] < 1.0 / EXCESS_FACTOR]
    bite = bitten[0] if bitten else n
    # The depth below the run's own inertial range, read where the spectrum is
    # still above roundoff. At the truncation itself it is not: reported too, as
    # a diagnostic, with the flag that says so.
    probe = int(round(DEPTH_FRACTION * n))
    live = [int(k) for k in m[1:] if spec[k] > FLOOR_MARGIN * floor]
    readable = bool(floor > 0 and spec[probe] > FLOOR_MARGIN * floor)
    depth = float(np.log10(ratio[probe])) if readable and ratio[probe] > 0 else float("nan")
    at_trunc = float(np.log10(ratio[n])) if ratio[n] > 0 else float("nan")
    # The energy the filter has taken out of the resolved spectrum, as a
    # fraction of what the inertial range would have carried over the same
    # wavenumbers. Only the deficit counts: a pile-up is a different failure and
    # spectral_tail.py is what reports it.
    over = (m >= lo) & (m <= n)
    deficit = float(np.clip(law[over] - spec[over], 0.0, None).sum())
    removed = deficit / float(law[over].sum())
    # And as a fraction of the whole eddy kinetic energy the run carries, which
    # is the number "how much of this rung survives" should be read in.
    removed_of_total = deficit / float(spec[m >= 1].sum() + deficit)
    return {"inertial_slope": float(fit[0]),
            "bite_wavenumber": int(bite), "bite_fraction": bite / n,
            "depth_dex": depth, "depth_wavenumber": probe,
            "depth_readable": readable,
            "roundoff_floor": floor,
            "last_live_wavenumber": max(live) if live else 0,
            "last_live_fraction": (max(live) if live else 0) / n,
            "depth_at_truncation_dex": at_trunc,
            "ke_removed_fraction_of_law": float(removed),
            "ke_removed_fraction_of_total": float(removed_of_total)}


def windows(first: int, last: int, length: int) -> list[tuple[int, int]]:
    out, a = [], first
    while a + length - 1 <= last:
        out.append((a, a + length - 1))
        a += length
    return out


def scatter(run_dir: Path, first: int, last: int, n: int, length: int,
            level: int | None) -> tuple[dict, dict, int]:
    """Mean and standard deviation of each quantity across disjoint windows."""
    wins = windows(first, last, length)
    rows = [measure(ke_spectrum(run_dir, a, b, level), n) for a, b in wins]
    keys = [k for k in rows[0] if isinstance(rows[0][k], float) or k.endswith("fraction")]
    mean = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}
    sd = {k: float(np.std([r[k] for r in rows], ddof=1)) if len(rows) > 1 else float("nan")
          for k in rows[0]}
    return mean, sd, len(rows)


def sensitivity(kappa: float, gamma: int, dt_s: float, tau_s: float,
                n: int, factor: float) -> dict:
    """What a change in kappa by `factor` is worth, in the reported units."""
    def bite(k: float) -> float:
        return (dt_s / (2.0 * k * tau_s)) ** (1.0 / (gamma - 1))

    def depth_dex(k: float) -> float:
        # At DEPTH_FRACTION of the truncation, which is where the measurement is
        # read; at the truncation itself the measurement is roundoff.
        return -math.log10(1.0 + (2.0 * k * tau_s / dt_s)
                           * DEPTH_FRACTION ** (gamma - 1))

    weaker, stronger = kappa / factor, kappa * factor
    return {"kappa": kappa, "gamma": gamma, "dt_seconds": dt_s,
            "tau_vorticity_seconds": tau_s,
            "R_at_truncation": 2.0 * kappa * tau_s / dt_s,
            "predicted_bite_fraction": bite(kappa),
            "depth_fraction": DEPTH_FRACTION,
            "predicted_depth_dex": depth_dex(kappa),
            "factor": factor,
            "d_bite_fraction": abs(bite(weaker) - bite(kappa)),
            "d_bite_fraction_stronger": abs(bite(stronger) - bite(kappa)),
            "d_depth_dex": abs(depth_dex(weaker) - depth_dex(kappa)),
            "d_depth_dex_stronger": abs(depth_dex(stronger) - depth_dex(kappa)),
            "bite_quantisation": 1.0 / n}


def verdict(sens: dict, sd: dict, n_windows: int) -> dict:
    if n_windows < MIN_WINDOWS:
        return {"resolvable_bite": None, "resolvable_depth": None,
                "reason": f"only {n_windows} windows; {MIN_WINDOWS} needed"}
    bar_bite = RESOLVE_SIGMA * sd["bite_fraction"]
    bar_depth = RESOLVE_SIGMA * sd["depth_dex"]
    ok_bite = (sens["d_bite_fraction"] >= bar_bite
               and sens["d_bite_fraction"] >= sens["bite_quantisation"])
    ok_depth = sens["d_depth_dex"] >= bar_depth
    return {"resolvable_bite": bool(ok_bite), "resolvable_depth": bool(ok_depth),
            "bar_bite_fraction": bar_bite, "bar_depth_dex": bar_depth}


def report(run_dir: Path, first: int, last: int, length: int,
           level: int | None, factor: float) -> dict:
    nl = namelist_keys(run_dir)
    rung = rung_of(run_dir)
    n = truncation(run_dir, None)
    kappa = float(nl.get("FILTERKAPPA", "nan"))
    gamma = int(float(nl.get("NFILTEREXP", "nan")))
    dt_s = float(nl["MPSTEP"]) * 60.0
    tau_s = cascade_time_s(rung)
    mean, sd, k = scatter(run_dir, first, last, n, length, level)
    sens = sensitivity(kappa, gamma, dt_s, tau_s, n, factor)
    ver = verdict(sens, sd, k)
    row = {"run": run_dir.name, "rung": rung, "truncation": n,
           "orbits": [first, last], "window_orbits": length, "windows": k,
           "namelist": {"FILTERKAPPA": nl.get("FILTERKAPPA"),
                        "NFILTEREXP": nl.get("NFILTEREXP"),
                        "MPSTEP": nl.get("MPSTEP"), "NDEL": nl.get("NDEL"),
                        "NHDIFF": nl.get("NHDIFF"), "TDISSZ": nl.get("TDISSZ")},
           "mean": mean, "sd": sd, "sensitivity": sens, "verdict": ver}
    print(f"{run_dir.name}  {rung}  kappa {kappa:g}  gamma {gamma}  "
          f"dt {dt_s/60:g} min  orbits {first}-{last} in {k} windows of {length}")
    print(f"  inertial slope            {mean['inertial_slope']:+.3f} "
          f"+/- {sd['inertial_slope']:.3f}")
    print(f"  bite fraction             {mean['bite_fraction']:.3f} "
          f"+/- {sd['bite_fraction']:.3f}   (predicted {sens['predicted_bite_fraction']:.3f}, "
          f"quantisation {1.0/n:.3f})")
    print(f"  last live wavenumber      m={mean['last_live_wavenumber']:.1f} = "
          f"{mean['last_live_fraction']:.3f} of the truncation "
          f"(above roundoff by {FLOOR_MARGIN:g}x)")
    print(f"  depth at {DEPTH_FRACTION:g}N, dex        {mean['depth_dex']:+.3f} "
          f"+/- {sd['depth_dex']:.3f}   "
          f"(predicted {sens['predicted_depth_dex']:+.3f}; "
          f"at the truncation {mean['depth_at_truncation_dex']:+.2f}, which is roundoff)")
    print(f"  eddy KE removed           {mean['ke_removed_fraction_of_total']*100:.3f}% "
          f"+/- {sd['ke_removed_fraction_of_total']*100:.3f}% of the run's own eddy KE")
    print(f"  a factor {factor:g} in kappa is worth: "
          f"bite {sens['d_bite_fraction']:.4f}, depth {sens['d_depth_dex']:.3f} dex")
    if ver.get("reason"):
        print(f"  verdict: no scatter -- {ver['reason']}")
    else:
        print(f"  bar ({RESOLVE_SIGMA:g} sigma): bite {ver['bar_bite_fraction']:.4f}, "
              f"depth {ver['bar_depth_dex']:.3f} dex")
        print(f"  verdict: confinement resolves a factor {factor:g} in kappa: "
              f"{'YES' if ver['resolvable_bite'] else 'NO'};  "
              f"depth at {DEPTH_FRACTION:g}N resolves it: "
              f"{'YES' if ver['resolvable_depth'] else 'NO'}")
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", type=Path, nargs="?")
    ap.add_argument("--compare", type=Path, nargs="+", default=None,
                    help="two or more run directories to report side by side")
    ap.add_argument("--first", type=int, default=0)
    ap.add_argument("--last", type=int, default=None,
                    help="default is the last orbit the run wrote")
    ap.add_argument("--window", type=int, default=WINDOW_ORBITS,
                    help="orbits per disjoint window; the scatter is across these")
    ap.add_argument("--level", type=int, default=None)
    ap.add_argument("--factor", type=float, default=KAPPA_FACTOR)
    ap.add_argument("--out", type=Path, default=ANALYSIS / "filter_spectral_cost.json")
    args = ap.parse_args()

    runs = list(args.compare) if args.compare else ([args.run_dir] if args.run_dir else [])
    if not runs:
        raise SystemExit("give a run directory or --compare")

    rows = []
    for run in runs:
        last = args.last
        if last is None:
            outs = sorted(run.glob("MOST.*.nc"))
            if not outs:
                raise SystemExit(f"no MOST.*.nc in {run}")
            last = int(outs[-1].stem.split(".")[-1])
        rows.append(report(run, args.first, last, args.window, args.level, args.factor))
        print()

    if len(rows) > 1:
        print("pairwise, against the pooled window scatter:")
        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                a, b = rows[i], rows[j]
                for key, unit in (("bite_fraction", ""),
                                  ("depth_dex", " dex"),
                                  ("last_live_fraction", ""),
                                  ("ke_removed_fraction_of_total", "")):
                    d = b["mean"][key] - a["mean"][key]
                    pooled = math.hypot(a["sd"][key], b["sd"][key])
                    sigmas = abs(d) / pooled if pooled > 0 else float("inf")
                    print(f"  {a['run']} -> {b['run']}  {key:32s} "
                          f"{d:+.4f}{unit}  = {sigmas:.1f} sigma "
                          f"({'SEPARATED' if sigmas >= RESOLVE_SIGMA else 'not separated'})")

    payload = {"note": "Spectral price of the physics filter, and whether the "
                       "instrument can resolve a change in filter_kappa. "
                       "Criteria fixed in the script before any arm ran. "
                       "exoplasim/scripts/filter_spectral_cost.py",
               "generated": datetime.now(timezone.utc).isoformat(),
               "criteria": {"kappa_factor": KAPPA_FACTOR,
                            "resolve_sigma": RESOLVE_SIGMA,
                            "min_windows": MIN_WINDOWS,
                            "window_orbits": WINDOW_ORBITS,
                            "excess_factor": EXCESS_FACTOR},
               "runs": rows}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2))
    print(f"wrote {rel(args.out)}")


if __name__ == "__main__":
    main()
