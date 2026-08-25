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

THE LAW, AND WHY IT IS THIS ONE. `legmod.f90` folds `f(n) = exp(-kappa x^gamma)`,
`x = n/NTRU`, into every per-mode transform weight in both directions, and what
passes through those weights is the nonlinear TENDENCY: no filtered value is
ever stored back into the prognostic spectral state. So the filter attenuates
the SOURCE that fills a wavenumber by a fixed factor rather than damping the
state at a rate, and a source scaled by `f^2` against an unchanged sink
equilibrates at

    spec(x) / law(x) = f(x)^2 = exp(-2 kappa x^gamma)

where `law` is the run's own inertial range extrapolated. It carries no
timestep, and the tendency's own `dt` factor cancels exactly because a tendency
applied for a shorter step is applied more often.

THE LAW IT REPLACES, because this script argued for it and was wrong. The
relaxation model `spec/law = 1 / (1 + R x^(gamma-1))`, `R = 2 kappa tau / dt`,
treated the filter as a rate `2 kappa x^gamma / dt` against the flow's cascade
rate `x / tau_vorticity`. Two timestep pairs differing in `MPSTEP` alone, same
binary by sha and same cold start, move the measured spectrum by 0.1 sigma at
T21 and 0.2 sigma at T42 where that law predicts several times the instrument's
own three-sigma bar. The gamma pair separates the two laws where they disagree
by nearly a factor of two -- the bite point at gamma 8, `(ln2/2 kappa)^(1/8)` =
0.675 against the rate law's 0.381 -- and it measures 0.688.
`exoplasim/notes/filter-spectral-price.md` carries both.

WHAT IT SAYS ABOUT KAPPA, and it is the point of the script. The bite point
`x_f = (ln2 / (2 kappa))^(1/gamma)` goes as `kappa^(-1/gamma)`. At gamma 16 that
is `kappa^(-1/16)`: the confinement criterion is very nearly BLIND to kappa, and
a factor of two in it moves the bite by less than one wavenumber at either rung
on disk. The depth at `DEPTH_FRACTION` of the truncation is `2 kappa x^gamma /
ln 10` dex and is LINEAR in kappa. The two reported quantities therefore differ
by a factor of gamma in their leverage, and only one of them can carry a kappa
decision.

WHAT THE LAW DOES NOT COVER, stated because the residual is not small. The
measured response to a change in gamma is 0.81 of what `f^2` alone predicts, and
the shortfall is the same at both gammas: what is left after the filter is the
hyperdiffusion and the spectrum's own steepening away from a fit taken over the
middle band. `f^2` is the filter's share and not the whole depletion, so a
predicted depth here is a LOWER BOUND on what a run will measure.

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
import rungs  # noqa: E402  from lib/, the one rung-to-dimension mapping
from spectral_tail import FLOOR_MARGIN, ke_spectrum, truncation  # noqa: E402  one reader, one convention

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
# FLOOR_MARGIN -- how far above the measured roundoff floor a wavenumber must
# stand to be read -- is imported from `spectral_tail.py` rather than restated.
# It bounds the depth reading here and the tail fit there, and a second copy of
# it is a number that can drift away from the one the note cites.


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


def carried(n: int, kappa: float, gamma: int) -> dict:
    """How many zonal wavenumbers a truncation carries, and what the rest cost.

    NO RUN IS READ. Both halves are arithmetic on the truncation and on the
    filter the config declares, which is the only way to answer the question at
    a rung nothing has been integrated at.

    CARRIED, and the criterion is the one already declared above rather than a
    new one: a wavenumber is carried when the filter leaves at least
    `1/EXCESS_FACTOR` of what the flow's own inertial range would put there,
    which is the same factor two `spectral_tail.py` calls a bite. Under
    `spec/law = exp(-2 kappa x^gamma)` that is `m <= (ln2/(2 kappa))^(1/gamma) *
    NTRU`, a fixed FRACTION of every rung.

    IT IS AN UPPER BOUND, and the reason is the spectrum's coordinate. The
    filter is a function of TOTAL wavenumber `n` and this count is in ZONAL
    wavenumber `m`, which sums over `n >= m`. `f` decreases in `n`, so `f(m)^2`
    is the largest factor anything in the zonal wavenumber `m` survives by:
    a wavenumber this says is not carried is certainly not carried, and one it
    says is carried may still be depleted by the modes above the diagonal.

    WHAT THE DEAD TOP COSTS IS NOT ITS SHARE OF THE WAVENUMBERS. A triangular
    truncation holds `n = m .. NTRU`, so zonal wavenumber `m` carries
    `NTRU - m + 1` meridional modes and the Legendre work at the top of a rung
    is the small end of a triangle. The uncarried band is 18 per cent of the
    wavenumbers at every rung and a twentieth of that in modes, and because the
    filter is scale-free in `n/NTRU` the fraction is the SAME at every rung --
    so it cancels out of every rung-to-rung cost ratio rather than distorting
    one.
    """
    x_f = (math.log(EXCESS_FACTOR) / (2.0 * kappa)) ** (1.0 / gamma)
    last = int(math.floor(x_f * n))
    modes_total = (n + 1) * (n + 2) // 2
    above = n - last                      # wavenumbers m = last+1 .. n
    modes_above = above * (above + 1) // 2
    return {"truncation": n, "kappa": kappa, "gamma": gamma,
            "carried_fraction": x_f,
            "last_carried_wavenumber": last,
            "carried_wavenumbers": last + 1,
            "wavenumbers_total": n + 1,
            "uncarried_wavenumbers": above,
            "uncarried_share_of_wavenumbers": above / (n + 1),
            "modes_total": modes_total,
            "modes_uncarried": modes_above,
            "uncarried_share_of_legendre_work": modes_above / modes_total}


def rung_budget(kappa: float, gamma: int) -> list[dict]:
    """`carried` for every ladder rung, with the Legendre work each one buys."""
    rows = []
    base = None
    for rung in rungs.RUNGS:
        nlat, nlon, n = rungs.geometry(rung)
        row = carried(n, kappa, gamma)
        # Legendre work is modes x latitudes: every mode integrates over every
        # Gaussian latitude. That is the accounting SPAT-11's recorded ratios
        # already use -- 7.9x at T85 and 26.2x at T127 against T42 -- and it is
        # restated here so the wavenumber count and the work are not confused.
        row |= {"rung": rung, "latitudes": nlat, "longitudes": nlon,
                "legendre_work": row["modes_total"] * nlat}
        base = base or row["legendre_work"]
        row["legendre_work_vs_coarsest"] = row["legendre_work"] / base
        rows.append(row)
    return rows


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


def sensitivity(kappa: float, gamma: int, n: int, factor: float,
                dt_s: float) -> dict:
    """What a change in kappa by `factor` is worth, in the reported units.

    Under `spec/law = exp(-2 kappa x^gamma)`. `dt_s` is carried into the result
    for the record and is NOT used: the filter is folded into the transform
    weights and acts on the tendency, so its spectral effect has no timestep in
    it. An earlier version of this function solved a rate law that did, and the
    number it produced disagreed with the measured bite in both directions.
    """
    def bite(k: float) -> float:
        return (math.log(2.0) / (2.0 * k)) ** (1.0 / gamma)

    def depth_dex(k: float) -> float:
        # At DEPTH_FRACTION of the truncation, which is where the measurement is
        # read; at the truncation itself the measurement is roundoff.
        return -2.0 * k * DEPTH_FRACTION ** gamma / math.log(10.0)

    weaker, stronger = kappa / factor, kappa * factor
    # THE BITE MOVE IS AN UPPER BOUND, and the gamma pair says by how much. What
    # a run measures as its bite point is where the TOTAL departure from the
    # inertial range reaches EXCESS_FACTOR, and the hyperdiffusion contributes
    # to that departure as well -- so the measured point moves by less than the
    # filter's own does. The gamma pair is the one arm that moves it: the
    # filter's bite goes 0.675 to 0.822 between gamma 8 and gamma 16 and the
    # measured bite goes 0.688 to 0.714, a transfer of about a fifth. A
    # `resolvable_bite` verdict below is therefore the optimistic half of the
    # answer and is labelled as one.
    return {"law": "spec/law = exp(-2*kappa*x**gamma); carries no timestep",
            "d_bite_is_an_upper_bound": True,
            "kappa": kappa, "gamma": gamma, "dt_seconds": dt_s,
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
    return {"resolvable_bite": bool(ok_bite),
            "resolvable_bite_is_upper_bound": True,
            "resolvable_depth": bool(ok_depth),
            "bar_bite_fraction": bar_bite, "bar_depth_dex": bar_depth}


def report(run_dir: Path, first: int, last: int, length: int,
           level: int | None, factor: float) -> dict:
    nl = namelist_keys(run_dir)
    rung = rung_of(run_dir)
    n = truncation(run_dir, None)
    kappa = float(nl.get("FILTERKAPPA", "nan"))
    gamma = int(float(nl.get("NFILTEREXP", "nan")))
    dt_s = float(nl["MPSTEP"]) * 60.0
    mean, sd, k = scatter(run_dir, first, last, n, length, level)
    sens = sensitivity(kappa, gamma, n, factor, dt_s)
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
              f"{'YES at most' if ver['resolvable_bite'] else 'NO'};  "
              f"depth at {DEPTH_FRACTION:g}N resolves it: "
              f"{'YES' if ver['resolvable_depth'] else 'NO'}")
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", type=Path, nargs="?")
    ap.add_argument("--rungs", action="store_true",
                    help="how many zonal wavenumbers each ladder rung carries "
                         "under the declared filter, and what the rest cost in "
                         "Legendre work. Reads no run: both halves are "
                         "arithmetic on the truncation and on config/planet.yaml.")
    ap.add_argument("--restate", type=Path, default=None,
                    help="recompute the `sensitivity` and `verdict` blocks of an "
                         "existing result file from its own recorded namelist "
                         "and truncation, under the current law, and rewrite it. "
                         "The measured `mean` and `sd` are untouched. This is "
                         "how a refuted law is corrected out of a result whose "
                         "run directories are gone.")
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

    if args.rungs:
        cfg = yaml.safe_load(
            (ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))["model"]
        kappa = float(cfg["filter_kappa"])
        gamma = int(cfg["filter_power"])
        rows = rung_budget(kappa, gamma)
        print(f"kappa {kappa:g}, gamma {gamma}: the filter leaves at least "
              f"1/{EXCESS_FACTOR:g} of the flow's own amplitude up to "
              f"{rows[0]['carried_fraction']:.3f} of the truncation, at every rung.\n")
        print(f"{'rung':>6} {'NTRU':>5} {'carried m':>10} {'of':>4} "
              f"{'modes':>7} {'uncarried':>10} {'work':>9}")
        for r in rows:
            print(f"{r['rung']:>6} {r['truncation']:>5} "
                  f"{r['carried_wavenumbers']:>10} {r['wavenumbers_total']:>4} "
                  f"{r['modes_total']:>7} "
                  f"{r['uncarried_share_of_legendre_work']*100:>9.2f}% "
                  f"{r["legendre_work_vs_coarsest"]:>8.1f}x")
        payload = {"note": "How many zonal wavenumbers each rung carries under "
                           "the declared filter, and what the uncarried top "
                           "costs in Legendre work. Arithmetic on the "
                           "truncation and config/planet.yaml; no run is read. "
                           "exoplasim/scripts/filter_spectral_cost.py --rungs",
                   "generated": datetime.now(timezone.utc).isoformat(),
                   "criteria": {"excess_factor": EXCESS_FACTOR},
                   "rungs": rows}
        out = args.out if args.out != ANALYSIS / "filter_spectral_cost.json" \
            else ANALYSIS / "rung_spectral_budget.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=2))
        print(f"\nwrote {rel(out)}")
        return

    if args.restate is not None:
        payload = json.loads(args.restate.read_text())
        for row in payload.get("runs", []):
            nl = row["namelist"]
            n = int(row["truncation"])
            k = float(nl["FILTERKAPPA"])
            g = int(float(nl["NFILTEREXP"]))
            dt_s = float(nl["MPSTEP"]) * 60.0
            row["sensitivity"] = sensitivity(k, g, n, args.factor, dt_s)
            row["verdict"] = verdict(row["sensitivity"], row["sd"], row["windows"])
            print(f"{row['run']}  {row['rung']}  bite "
                  f"{row['sensitivity']['predicted_bite_fraction']:.3f}, depth "
                  f"{row['sensitivity']['predicted_depth_dex']:+.3f} dex, "
                  f"confinement resolves a factor {args.factor:g}: "
                  f"{row['verdict'].get('resolvable_bite')}, depth: "
                  f"{row['verdict'].get('resolvable_depth')}")
        payload["restated"] = {
            "at": datetime.now(timezone.utc).isoformat(),
            "why": "The `sensitivity` and `verdict` blocks were computed from a "
                   "rate law two timestep pairs refute. They are recomputed "
                   "from this file's own recorded namelist under "
                   "spec/law = exp(-2*kappa*x**gamma). The measured `mean` and "
                   "`sd` are as taken and are untouched."}
        args.restate.write_text(json.dumps(payload, indent=2))
        print(f"rewrote {rel(args.restate)}")
        return

    runs = list(args.compare) if args.compare else ([args.run_dir] if args.run_dir else [])
    if not runs:
        raise SystemExit("give a run directory, --compare, --rungs or --restate")

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
