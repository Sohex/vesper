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

THE WINDOW, ALSO FIXED IN ADVANCE. A scaled state is out of equilibrium and the
flow regenerates toward its natural amplitude, which would contaminate a trend
fitted across the whole orbit. So the sink is fitted only while the global eddy
kinetic energy stays within +/-25 percent of its value at the start of the arm,
and an arm whose window holds fewer than five samples is reported as
unmeasurable rather than fitted.
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
WINDOW_TOLERANCE = 0.25   # eddy KE may drift this far before the arm is regenerating
MIN_WINDOW_SAMPLES = 5


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


def fit_sink(s: dict, seconds_per_step: float):
    """W/m2, over the window where the arm is not yet regenerating."""
    e0 = s["eddy_ke"][0]
    inside = np.abs(s["eddy_ke"] / e0 - 1.0) <= WINDOW_TOLERANCE
    stop = len(inside)
    for i, ok in enumerate(inside):
        if not ok:
            stop = i
            break
    if stop < MIN_WINDOW_SAMPLES:
        return None, stop
    t = s["time"][:stop] * seconds_per_step
    slope = float(np.polyfit(t, s["total"][:stop], 1)[0])
    return slope, stop


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
        sink, stop = fit_sink(s, dt)
        rows.append({
            "eps": float(eps), "run": run.name, "samples": int(len(s["time"])),
            "window_samples": int(stop),
            "sink_w_m2": sink,
            "eddy_ke_initial_j_m2": float(s["eddy_ke"][0]),
            "mass_drift_ppm": float((s["mass"][-1] / s["mass"][0] - 1) * 1e6),
            "enthalpy_trend_j_m2": float(s["enthalpy"][-1] - s["enthalpy"][0]),
            "kinetic_trend_j_m2": float(s["kinetic"][-1] - s["kinetic"][0]),
            "orographic_trend_j_m2": float(s["orographic"][-1] - s["orographic"][0]),
            "eddy_ke_series": [float(x) for x in s["eddy_ke"]],
            "total_series": [float(x) for x in s["total"]],
        })

    rows.sort(key=lambda r: -r["eps"])
    print(f"{'eps':>6} {'run':>18} {'window':>7} {'sink W/m2':>11} "
          f"{'eddy KE J/m2':>13} {'mass ppm':>9}")
    for r in rows:
        sink = "unmeasurable" if r["sink_w_m2"] is None else f"{r['sink_w_m2']:+.4f}"
        print(f"{r['eps']:>6.2f} {r['run']:>18} {r['window_samples']:>7d} "
              f"{sink:>11} {r['eddy_ke_initial_j_m2']:>13.4g} "
              f"{r['mass_drift_ppm']:>+9.1f}")

    usable = [r for r in rows if r["sink_w_m2"] is not None
              and r["sink_w_m2"] < 0 and r["eps"] > 0]
    exponent = None
    attribution = "not attributable: fewer than two arms with a measurable sink"
    if len(usable) >= 2:
        x = np.log([r["eps"] for r in usable])
        y = np.log([abs(r["sink_w_m2"]) for r in usable])
        exponent = float(np.polyfit(x, y, 1)[0])
        resid = float(np.max(np.abs(y - np.polyval(np.polyfit(x, y, 1), x))))
        attribution = ("no single band: the sink is a mixture of orders, or the "
                       "window is contaminated")
        for lo, hi, label in BANDS:
            if lo <= exponent <= hi:
                attribution = label
                break
        print(f"\nexponent in eps: {exponent:+.2f}  "
              f"(worst log residual {resid:.3f})")
        print(f"attribution: {attribution}")
        print("pairwise exponents:")
        for a, b in zip(usable, usable[1:]):
            p = np.log(abs(a["sink_w_m2"]) / abs(b["sink_w_m2"])) / np.log(a["eps"] / b["eps"])
            print(f"  {a['eps']:.2f} -> {b['eps']:.2f} : {p:+.2f}")

    payload = {"note": "dry adiabatic energy sink against flow amplitude, for "
                       "world-bxr; bands fixed in the script before any arm ran. "
                       "exoplasim/scripts/dry_energy_order.py",
               "generated": datetime.now(timezone.utc).isoformat(),
               "gravity_m_s2": gravity,
               "window_tolerance": WINDOW_TOLERANCE,
               "min_window_samples": MIN_WINDOW_SAMPLES,
               "bands": [{"low": lo, "high": hi, "term": t} for lo, hi, t in BANDS],
               "arms": rows, "exponent": exponent, "attribution": attribution}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
