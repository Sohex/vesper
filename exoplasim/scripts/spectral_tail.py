#!/usr/bin/env python3
"""Is the damping absorbing the cascade, or letting it pile up at the grid scale?

    python exoplasim/scripts/spectral_tail.py <run_dir> --first 0 --last 1

Worldbuilding frame: a numerical diagnostic on the Vesper project's climate
model. Nothing here is about the simulated planet.

WHY. `model.hyperdiffusion` is derived rather than inherited
(`exoplasim/notes/resolution-tuned-parameters.md`), and a derivation needs
something that can refute it. "It ran" cannot: the ladder's own arms show a
configuration surviving one orbit and failing the next. What the damping is FOR
is absorbing the enstrophy cascade at the truncation without reaching the scales
the model was bought to resolve, and both halves of that are visible in the
kinetic-energy spectrum.

  TOO LITTLE and energy accumulates at the smallest resolved scales: the
  spectrum flattens or turns up near the truncation, which is the classic
  spectral blocking that precedes a blow-up.

  TOO MUCH and the spectrum bends down well below the truncation: resolution
  that was paid for in wall clock is being damped away.

IN ZONAL WAVENUMBER, not total. The postprocessed output carries gridpoint winds
rather than spectral coefficients, so the transform available without a full
spherical harmonic analysis is an FFT along longitude. That is enough for this
question -- a pile-up at the grid scale shows in `m` as clearly as in `n` -- and
it is stated because the two are not the same spectrum and the slopes are not
directly comparable to a quoted `n^-3`.

**THE TRUNCATION IS THE MODEL'S, NOT THE FFT'S**, and getting this wrong is how
the first run of this diagnostic produced three wrong conclusions. A T42 run on
128 longitudes gives an FFT out to m=64, but the model represents nothing above
m=42: the spectrum there sits at roundoff and not physics. Normalising by 64
made a bite at m=28 read as 0.44 of the truncation when it is 0.67, and turned a
filter that was inside its confinement requirement into one that appeared to be
damping away a third of the resolved spectrum.

**BUT THE BAND ABOVE THE TRUNCATION IS KEPT, because it is the instrument's own
noise floor.** Everything above `m = NTRU` is what this arithmetic does to zero,
so its median MEASURES the roundoff the whole spectrum sits on. Nothing is fitted
there and no criterion is read there; it is what says whether a wavenumber inside
the truncation is carrying a number or carrying nothing. Cutting the spectrum at
`m = NTRU` first, which this script used to do, throws away the only calibration
available and leaves no way to tell a steep tail from a fit through the floor.

THE CRITERIA, AND WHY THEY ARE THESE. The question is not the tail's slope in
isolation -- a steep tail is what damping is supposed to produce -- but WHERE the
spectrum departs from the flow's own inertial range. So fit the power law the
model itself produces over the middle band, extrapolate it, and find the
wavenumber at which the measured spectrum leaves it. Call that the bite point.

  1. NO PILE-UP. The spectrum must nowhere exceed its own extrapolated power law
     by more than a factor of two. Energy above the inertial range is the
     blocking that precedes a blow-up.
  2. CONFINEMENT. The bite point must sit at or above 0.6 of the truncation.
     This follows the confinement requirement the damping was designed to --
     damping no more than about a tenth of the local cascade rate at 0.7 of the
     truncation -- with margin for measurement slop, and it is NOT read off any
     arm's result.

An earlier version of this script tested the tail slope against a floor and the
middle slope against the tail. That pair was ill-posed in two ways at once. The
comparison could not fire in the direction it was written, and the tail slope it
compared was not a property of the damping: fitted to `m = NTRU` on a spectrum
whose top sits at roundoff, it measures how far the run's inertial range is above
the floating-point floor. The replacement asks where the damping starts rather
than how hard it finishes, and the tail slope is now REPORTED rather than tested,
fitted only over wavenumbers that stand clear of the measured floor and labelled
with the band it was fitted over.
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _paths  # noqa: F401
from gridding import gaussian_row_weights  # noqa: E402  from lib/, via _paths
from paths import rel  # noqa: E402  from lib/, put on sys.path by _paths

ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "exoplasim" / "analysis"

EXCESS_FACTOR = 2.0        # above its own power law, this is a pile-up
MIN_BITE_FRACTION = 0.6    # the damping may not start below this fraction of N
# A wavenumber is READABLE when it stands this far above the measured roundoff
# floor. It bounds the tail fit here and the depth reading in
# `filter_spectral_cost.py`, which imports it: one bar, defined once, because a
# second copy of it is a number that can drift away from the note both scripts
# cite. It bounds NEITHER criterion above -- the bite point is the first
# departure BELOW the inertial range and the pile-up test looks for excess
# above it, so roundoff at the top cannot move either one.
FLOOR_MARGIN = 100.0


def ke_spectrum(run_dir: Path, first: int, last: int, level: int | None):
    """Kinetic energy by zonal wavenumber, area-weighted, time-averaged."""
    import netCDF4 as nc
    files = [f for f in sorted(glob.glob(str(run_dir / "MOST.*.nc")))
             if first <= int(Path(f).stem.split(".")[-1]) <= last]
    if not files:
        raise SystemExit(f"no MOST.*.nc in {run_dir} for orbits {first}-{last}")
    total = None
    for path in files:
        d = nc.Dataset(path)
        ua, va = d.variables["ua"][:], d.variables["va"][:]
        lat = d.variables["lat"][:]
        sel = slice(None) if level is None else slice(level, level + 1)
        ua, va = ua[:, sel], va[:, sel]
        # EDDY field only: the zonal mean is wavenumber zero and carries the jet,
        # which is not what cascades and would dominate the first bin.
        ua = ua - ua.mean(axis=3, keepdims=True)
        va = va - va.mean(axis=3, keepdims=True)
        fu = np.fft.rfft(ua, axis=3) / ua.shape[3]
        fv = np.fft.rfft(va, axis=3) / va.shape[3]
        power = 0.5 * (np.abs(fu) ** 2 + np.abs(fv) ** 2)
        w = gaussian_row_weights(np.asarray(lat, dtype=float),
                                 what=str(path))[None, None, :, None]
        spec = (power * w).sum(axis=2) / w.sum()      # over latitude
        spec = spec.mean(axis=(0, 1))                 # over time and level
        total = spec if total is None else total + spec
        d.close()
    return total / len(files)


def truncation(run_dir: Path, override: int | None) -> int:
    """The model's spectral truncation, which bounds the meaningful spectrum.

    Not the FFT's Nyquist: a T42 run on 128 longitudes transforms out to m=64
    and represents nothing above m=42. Taken from the run's own manifest so it
    cannot disagree with what ran.
    """
    if override:
        return int(override)
    manifest = run_dir / "run_manifest.json"
    if manifest.is_file():
        res = json.loads(manifest.read_text())["physical"]["resolution"]
        return int(str(res).lstrip("Tt"))
    raise SystemExit(
        f"no run_manifest.json in {run_dir} and no --truncation given; the "
        f"spectrum cannot be normalised by a truncation nobody named")


def slope(spec: np.ndarray, lo: int, hi: int) -> float:
    m = np.arange(len(spec))
    band = (m >= lo) & (m <= hi) & (spec > 0)
    if band.sum() < 3:
        return float("nan")
    return float(np.polyfit(np.log(m[band]), np.log(spec[band]), 1)[0])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir", type=Path)
    ap.add_argument("--first", type=int, default=0)
    ap.add_argument("--last", type=int, default=0)
    ap.add_argument("--truncation", type=int, default=None,
                    help="spectral truncation; default is the run manifest's")
    ap.add_argument("--level", type=int, default=None,
                    help="model level index; default is the column mean")
    ap.add_argument("--out", type=Path, default=ANALYSIS / "spectral_tail.json")
    args = ap.parse_args()

    full = ke_spectrum(args.run_dir, args.first, args.last, args.level)
    n = truncation(args.run_dir, args.truncation)
    # The dead band above the model's truncation, kept and used before the
    # spectrum is cut down to what the model represents.
    floor = float(np.median(full[n + 1:])) if len(full) > n + 1 else float("nan")
    spec = full[:n + 1]
    m = np.arange(len(spec))
    # The fit band ends at a THIRD of the truncation, not a half, so the search
    # above it can resolve a bite point well below 0.6N. With the band ending at
    # N/2 the smallest detectable bite fraction is 0.52 by construction, which
    # would have made the confinement criterion untestable against itself.
    lo, hi = max(2, n // 8), max(3, n // 3)
    middle = slope(spec, lo, hi)

    # The flow's own inertial range, extrapolated across the whole spectrum.
    band = (m >= lo) & (m <= hi) & (spec > 0)
    fit = np.polyfit(np.log(m[band]), np.log(spec[band]), 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        law = np.exp(np.polyval(fit, np.log(np.where(m > 0, m, 1))))
        ratio = np.where(law > 0, spec / law, np.nan)

    upper = (m > hi) & (m <= n)
    excess = float(np.nanmax(ratio[upper])) if upper.any() else float("nan")
    bitten = [int(k) for k in m[upper] if ratio[k] < 1.0 / EXCESS_FACTOR]
    bite = bitten[0] if bitten else n
    bite_fraction = bite / n

    pile_up = excess > EXCESS_FACTOR
    too_broad = bite_fraction < MIN_BITE_FRACTION
    verdict = ("pile-up above the inertial range" if pile_up else
               "damping starts too far down" if too_broad else "clean")

    # The tail slope, over the wavenumbers that carry a number. Fitting to the
    # truncation measures the distance from the inertial range to the roundoff
    # floor and reports it as damping.
    if floor > 0:
        live = [int(k) for k in m[1:] if spec[k] > FLOOR_MARGIN * floor]
        last_live = max(live) if live else 0
        tail_hi = min(n, last_live)
        tail = slope(spec, hi, tail_hi)
        refusal = ("" if tail == tail else
                   f"fewer than three wavenumbers between m={hi} and the last "
                   f"one standing {FLOOR_MARGIN:g}x clear of roundoff")
    else:
        # No band above the truncation, so the floor is unmeasured and a tail
        # slope cannot be told from a fit through it.
        last_live, tail_hi, tail = 0, 0, float("nan")
        refusal = (f"the FFT reaches only m={len(full) - 1} and the model "
                   f"truncates at m={n}, so there is no dead band to measure "
                   f"the roundoff floor in")

    print(f"{args.run_dir.name}, orbits {args.first}-{args.last}, "
          f"{n} zonal wavenumbers")
    print(f"  inertial-range slope, m {lo}-{hi} : {middle:+.2f}")
    if refusal:
        print(f"  tail slope                      : refused, {refusal}")
    else:
        print(f"  tail slope, m {hi}-{tail_hi}             : {tail:+.2f}"
              + (f"   ({n - tail_hi} wavenumber(s) at roundoff dropped)"
                 if tail_hi < n else ""))
    print(f"  roundoff floor, median above m={n}: {floor:.3g}")
    print(f"  peak excess over the power law  : {excess:.2f}x  "
          f"(pile-up above {EXCESS_FACTOR:g}x)")
    print(f"  bite point                      : m={bite} = {bite_fraction:.2f} "
          f"of the truncation (must be >= {MIN_BITE_FRACTION:g})")
    print(f"  verdict: {verdict}")

    payload = {"note": "KE by zonal wavenumber; criteria fixed in the script "
                       "before any arm ran. exoplasim/scripts/spectral_tail.py",
               "generated": datetime.now(timezone.utc).isoformat(),
               "run": args.run_dir.name, "first": args.first, "last": args.last,
               "level": args.level, "wavenumbers": n,
               "tail_slope": tail if tail == tail else None,
               "tail_band": None if refusal else [int(hi), int(tail_hi)],
               "tail_refused": refusal or None,
               "roundoff_floor": floor if floor == floor else None,
               "last_live_wavenumber": last_live,
               "last_live_fraction": last_live / n,
               "inertial_slope": middle,
               "peak_excess": excess, "bite_wavenumber": bite,
               "bite_fraction": bite_fraction,
               "criteria": {"excess_factor": EXCESS_FACTOR,
                            "min_bite_fraction": MIN_BITE_FRACTION,
                            "floor_margin": FLOOR_MARGIN},
               "verdict": verdict,
               "spectrum": [float(x) for x in spec]}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2))
    print(f"wrote {rel(args.out)}")
    raise SystemExit(0 if verdict == "clean" else 1)


if __name__ == "__main__":
    main()
