#!/usr/bin/env python3
"""Which ORDER in the flow amplitude does the dry adiabatic energy sink carry?

    python exoplasim/scripts/dry_energy_order.py 1.0=run_a 0.5=run_b 0.25=run_c

Worldbuilding frame: a diagnostic on the Vesper project's climate model. The
"atmosphere" is a prognostic array and the sink is a numerical one.

WHAT IT ANSWERS. `world-bxr` measures an energy sink in the dynamical core with
every diabatic source removed, and several candidate terms could carry it. They
are not distinguishable by size, but they ARE distinguishable by how the sink
responds when the whole flow is multiplied by a factor, because each term enters
the energy tendency at a different power of that factor:

    eps^1   a term linear in the state against the FIXED boundary -- here the
            orographic geopotential, which `scale_restart.py` deliberately does
            not scale
    eps^2   the linear scheme: the semi-implicit split, the leapfrog, the
            vertical discretisation; or a product of the state with the
            unscaled reference profile t0
    eps^3   the quadratic nonlinear terms, whose spectral truncation is the
            standing hypothesis
    eps^4   the cubic terms. This grid does not dealias them: NLON = 3*NTRU+1
            dealiases a product of two fields, and the sigma-coordinate
            tendencies in `calcgp` carry products of three.

THE ATTRIBUTION BANDS ARE FIXED HERE, BEFORE ANY ARM RAN, and an exponent that
lands between them is reported as unattributed rather than rounded into the
nearest one.

THE ENERGY. For a hydrostatic column the internal and potential energies combine
into the enthalpy, so the conserved total is

    E = column enthalpy + column kinetic energy + surface geopotential * mass

with the enthalpy taken from the model's own `denergy01` rather than rebuilt
from `ta`, so the specific heat, the humidity correction and the sigma weights
are the model's and cannot disagree with it. The kinetic term needs `dsigma`,
which the output does not carry; it is recovered from the full levels by the
model's own midpoint rule, sigma(j) = (sigmah(j-1) + sigmah(j)) / 2 with
sigmah(0) = 0, and the recursion is CHECKED against sigmah(NLEV) = 1.

Global means use Gaussian weights, not cos(lat): the output is on the model's
Gaussian grid and the quadrature that is exact there is the one the model
integrates with.

HOW THE ARMS ARE COMPARED. A dry adiabatic atmosphere has no source: with the
radiation emptied there is nothing to maintain the temperature gradient the
eddies feed on, so the flow DECAYS -- by a factor of four over one orbit at full
amplitude. Nothing regenerates and no arm returns to a common state, which means
a single trend per arm would compare arms at different stages of their own
spin-down.

So the exponent is fitted ACROSS ARMS AT MATCHED MODEL TIME, sample by sample.
At the first sample the arms are exact scaled copies of one another by
construction, and they drift apart from there, so the exponent is reported as a
series and attributed from its median over the FIRST QUARTER of the record.
That rule is fixed here and looks only at the model time and the scale factor,
never at the sink, so it cannot be tuned toward an answer.

The instantaneous sink is a centred difference of the total energy, one-sided at
the ends.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _paths  # noqa: F401

ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "exoplasim" / "analysis"

# Attribution bands, fixed before any arm ran. See the module docstring.
BANDS = [
    (0.70, 1.40, "linear against the fixed orography"),
    (1.70, 2.40, "the linear scheme, or a product with the reference profile"),
    (2.70, 3.40, "the quadratic nonlinear terms"),
    (3.70, 4.40, "the cubic terms, which this grid does not dealias"),
]
ATTRIBUTION_FRACTION = 0.25   # the leading quarter of the record, where the
                              # arms are still near-scaled copies of each other


def gaussian_weights(nlat: int) -> np.ndarray:
    _, w = np.polynomial.legendre.leggauss(nlat)
    return w[::-1] / w.sum()          # model latitudes run north to south


def dsigma_from_levels(lev: np.ndarray) -> np.ndarray:
    """The model's half levels, recovered from its full ones and then checked."""
    sigmah = np.empty(len(lev) + 1)
    sigmah[0] = 0.0
    for j, s in enumerate(lev):
        sigmah[j + 1] = 2.0 * s - sigmah[j]
    if abs(sigmah[-1] - 1.0) > 1e-4:
        raise SystemExit(
            f"the midpoint recursion gives sigmah(NLEV) = {sigmah[-1]:.6f}, not "
            f"1: these are not the levels this rule describes and dsigma cannot "
            f"be trusted")
    return np.diff(sigmah)


def arm_series(run_dir: Path, gravity: float):
    """Total energy per unit area and eddy kinetic energy, per output time."""
    import netCDF4 as nc
    files = sorted(run_dir.glob("MOST.*.nc"))
    if not files:
        raise SystemExit(f"no MOST.*.nc in {run_dir}")
    times, etot, enth, kin, oro, eddy, mass = [], [], [], [], [], [], []
    for path in files:
        d = nc.Dataset(path)
        lev = np.asarray(d.variables["lev"][:], dtype=float)
        ds = dsigma_from_levels(lev)
        w = gaussian_weights(d.dimensions["lat"].size)[None, :, None]
        ua = np.asarray(d.variables["ua"][:], dtype=float)
        va = np.asarray(d.variables["va"][:], dtype=float)
        ps = np.asarray(d.variables["ps"][:], dtype=float) * 100.0   # hPa -> Pa
        h = np.asarray(d.variables["denergy01"][:], dtype=float)     # J/m2
        gz = np.asarray(d.variables["grnz"][:], dtype=float)         # m2/s2
        # column kinetic energy, J/m2
        ke_col = (ps / gravity) * np.einsum("tzyx,z->tyx",
                                            0.5 * (ua * ua + va * va), ds)
        # the mass column's potential energy on the topography it stands on
        oro_col = gz * ps / gravity
        # the eddy part, as the window criterion and the amplitude axis
        ue = ua - ua.mean(axis=3, keepdims=True)
        ve = va - va.mean(axis=3, keepdims=True)
        eddy_col = (ps / gravity) * np.einsum("tzyx,z->tyx",
                                              0.5 * (ue * ue + ve * ve), ds)

        def gm(x):
            return (x * w).sum(axis=(1, 2)) / (w.sum() * x.shape[2])

        times.append(np.asarray(d.variables["time"][:], dtype=float))
        enth.append(gm(h)); kin.append(gm(ke_col)); oro.append(gm(oro_col))
        eddy.append(gm(eddy_col)); mass.append(gm(ps) / gravity)
        d.close()
    cat = lambda xs: np.concatenate(xs)
    enth, kin, oro = cat(enth), cat(kin), cat(oro)
    return {"time": cat(times), "enthalpy": enth, "kinetic": kin,
            "orographic": oro, "total": enth + kin + oro,
            "eddy_ke": cat(eddy), "mass": cat(mass)}


def instantaneous_sink(s: dict, seconds_per_step: float) -> np.ndarray:
    """W/m2 at each sample: a centred difference of the total energy."""
    t = s["time"] * seconds_per_step
    return np.gradient(s["total"], t)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("arms", nargs="+",
                    help="EPS=run_dir, e.g. 0.5=exoplasim/runs/run_abc")
    ap.add_argument("--out", type=Path, default=ANALYSIS / "dry_energy_order.json")
    args = ap.parse_args()

    import yaml
    planet = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text())
    gravity = float(planet["planet"]["gravity_m_s2"])

    rows = []
    for spec in args.arms:
        eps, _, path = spec.partition("=")
        run = Path(path)
        if not run.is_absolute():
            run = ROOT / run
        manifest = json.loads((run / "run_manifest.json").read_text())
        dt = float(manifest["source_config"]["model"]["timestep_minutes"]) * 60.0
        s = arm_series(run, gravity)
        sink = instantaneous_sink(s, dt)
        total = s["total"]
        span = (s["time"][-1] - s["time"][0]) * dt
        rows.append({
            "eps": float(eps), "run": run.name, "samples": int(len(s["time"])),
            "mean_sink_w_m2": float((total[-1] - total[0]) / span),
            "eddy_ke_initial_j_m2": float(s["eddy_ke"][0]),
            "eddy_ke_final_over_initial": float(s["eddy_ke"][-1] / s["eddy_ke"][0]),
            "mass_drift_ppm": float((s["mass"][-1] / s["mass"][0] - 1) * 1e6),
            "enthalpy_trend_j_m2": float(s["enthalpy"][-1] - s["enthalpy"][0]),
            "kinetic_trend_j_m2": float(s["kinetic"][-1] - s["kinetic"][0]),
            "orographic_trend_j_m2": float(s["orographic"][-1] - s["orographic"][0]),
            "eddy_ke_series": [float(x) for x in s["eddy_ke"]],
            "total_series": [float(x) for x in total],
            "sink_series": [float(x) for x in sink],
        })

    rows.sort(key=lambda r: -r["eps"])
    print(f"{'eps':>6} {'run':>18} {'mean sink':>10} {'eddy KE J/m2':>13} "
          f"{'decayed to':>11} {'mass ppm':>9}")
    for r in rows:
        print(f"{r['eps']:>6.2f} {r['run']:>18} {r['mean_sink_w_m2']:>+10.4f} "
              f"{r['eddy_ke_initial_j_m2']:>13.4g} "
              f"{r['eddy_ke_final_over_initial']:>11.3f} "
              f"{r['mass_drift_ppm']:>+9.1f}")

    exponent = None
    attribution = "not attributable: fewer than two arms"
    series = []
    if len(rows) >= 2:
        n = min(len(r["sink_series"]) for r in rows)
        x = np.log([r["eps"] for r in rows])
        for k in range(n):
            y = np.array([r["sink_series"][k] for r in rows])
            if np.any(y >= 0):
                series.append(None)
                continue
            series.append(float(np.polyfit(x, np.log(-y), 1)[0]))
        lead = [v for v in series[:max(1, int(n * ATTRIBUTION_FRACTION))]
                if v is not None]
        if lead:
            exponent = float(np.median(lead))
            attribution = ("no single band: the sink is a mixture of orders, or "
                           "the arms have already diverged")
            for lo, hi, label in BANDS:
                if lo <= exponent <= hi:
                    attribution = label
                    break
        whole = [v for v in series if v is not None]
        print(f"\nexponent in eps, per sample (leading quarter): "
              + " ".join(f"{v:+.2f}" if v is not None else "  n/a"
                         for v in series[:max(1, int(n * ATTRIBUTION_FRACTION))]))
        print(f"median over the leading quarter : {exponent:+.2f}"
              if exponent is not None else "median: unavailable")
        if whole:
            print(f"median over the whole record   : {np.median(whole):+.2f}")
        print(f"attribution: {attribution}")
        print("\nfor reference, the exponent from the whole-record mean sinks:")
        y = np.log([-r["mean_sink_w_m2"] for r in rows])
        print(f"  {float(np.polyfit(x, y, 1)[0]):+.2f}")

    payload = {"note": "dry adiabatic energy sink against flow amplitude, for "
                       "world-bxr; bands fixed in the script before any arm ran. "
                       "exoplasim/scripts/dry_energy_order.py",
               "generated": datetime.now(timezone.utc).isoformat(),
               "gravity_m_s2": gravity,
               "attribution_fraction": ATTRIBUTION_FRACTION,
               "exponent_per_sample": series,
               "bands": [{"low": lo, "high": hi, "term": t} for lo, hi, t in BANDS],
               "arms": rows, "exponent": exponent, "attribution": attribution}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
