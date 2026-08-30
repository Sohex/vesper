#!/usr/bin/env python3
"""Price OCN-1's dormant multilayer slab with oceanmod's exact time step."""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config/planet.yaml"
DECLARATION = ROOT / "ocean/config/column_response.yaml"
OUTPUT = ROOT / "analysis/ocean_column_response.json"
SOURCE = ROOT / "vendor/exoplasim/exoplasim/plasim/src/oceanmod.f90"
sys.path.insert(0, str(ROOT / "exoplasim/scripts"))

from run_exoplasim import derive


def diffusion_step(temperature: np.ndarray, depths: np.ndarray,
                   diffusivity: float, dt: float) -> np.ndarray:
    """One call to oceanmod:vdiffo, reduced from its elimination verbatim."""
    t = np.asarray(temperature, dtype=float)
    h = np.asarray(depths, dtype=float)
    if t.shape != h.shape or t.ndim != 1 or np.any(h <= 0.0):
        raise ValueError("temperature and positive layer depths must be 1-D peers")
    n = len(h)
    if n == 1:
        return t.copy()
    # vdiffkl is uniform in this experiment, so oceanini's thickness-weighted
    # interface mean is the same coefficient at every interface.
    zk = dt * float(diffusivity) * 2.0 / (h[:-1] + h[1:])
    zebs = np.zeros(n)
    ztn = np.zeros(n)
    out = t.copy()
    zebs[0] = zk[0] / (h[0] + zk[0])
    ztn[0] = h[0] * t[0] / (h[0] + zk[0])
    for j in range(1, n - 1):
        den = h[j] + zk[j] + zk[j - 1] * (1.0 - zebs[j - 1])
        zebs[j] = zk[j] / den
        ztn[j] = (t[j] * h[j] + zk[j - 1] * ztn[j - 1]) / den
    out[-1] = ((t[-1] * h[-1] + zk[-1] * ztn[-2])
               / (h[-1] + zk[-1] * (1.0 - zebs[-2])))
    for j in range(n - 2, -1, -1):
        out[j] = ztn[j] + zebs[j] * out[j + 1]
    return out


def response(depths: list[float], diffusivity: float, dt: float,
             period: float, rho: float, cp: float) -> dict:
    """Exact discrete steady harmonic response to one W m-2 at the surface."""
    h = np.asarray(depths, dtype=float)
    n = len(h)
    operator = np.column_stack([
        diffusion_step(np.eye(n)[:, j], h, diffusivity, dt)
        for j in range(n)
    ])
    input_vector = np.zeros(n)
    input_vector[0] = dt / (rho * cp * h[0])
    phase_step = np.exp(1j * 2.0 * np.pi * dt / period)
    complex_response = np.linalg.solve(
        phase_step * np.eye(n) - operator, operator @ input_vector)
    # A closed no-flux column must conserve the thickness-weighted heat state.
    closure = max(float(abs(h @ operator[:, j] - h[j])) for j in range(n))
    return {
        "layers_m": [float(x) for x in h],
        "total_depth_m": float(h.sum()),
        "surface_amplitude_k_per_w_m2": float(abs(complex_response[0])),
        "surface_phase_degrees": float(np.angle(complex_response[0], deg=True)),
        "bottom_amplitude_k_per_w_m2": float(abs(complex_response[-1])),
        "diffusion_skin_depth_m": float(math.sqrt(
            2.0 * diffusivity / (2.0 * np.pi / period))),
        "thickness_weighted_closure_residual_m": closure,
    }


def main() -> None:
    cfg = yaml.safe_load(CONFIG.read_text())
    dec = yaml.safe_load(DECLARATION.read_text())
    source = SOURCE.read_text()
    model = derive(cfg, float(cfg["orbit"]["baseline_flux_earth"]))
    period = float(model["orbital_year_seconds"])
    solar_day = float(cfg["planet"]["rotation_hours"]) * 3600.0
    ntspd_match = re.search(r"integer\s*::\s*ntspd\s*=\s*(\d+)", source)
    if not ntspd_match:
        raise RuntimeError("cannot derive ocean steps per solar day from oceanmod")
    ntspd = int(ntspd_match.group(1))
    dt = solar_day / ntspd
    wc = dec["water_column"]
    rho, cp = float(wc["density_kg_m3"]), float(wc["heat_capacity_j_kg_k"])
    rows = []
    for name, profile in wc["candidates"].items():
        for k in wc["vertical_diffusivity_m2_s"]:
            intended = response(profile, float(k), dt, period, rho, cp)
            intended.update({"candidate": name,
                             "profile_interpretation": "declared_total_50m",
                             "vertical_diffusivity_m2_s": float(k)})
            rows.append(intended)
            # This is what the unmodified source actually executes after a
            # user supplies the same dlayer list: it replaces the final entry
            # with mldepth. Keep both answers; silently pricing the intended
            # subdivision would conceal the structural reason for rejection.
            if len(profile) > 1:
                overwritten = [*profile[:-1], float(wc["total_depth_m"])]
                actual = response(overwritten, float(k), dt, period, rho, cp)
                actual.update({"candidate": name,
                               "profile_interpretation": "unmodified_source_last_layer_overwrite",
                               "vertical_diffusivity_m2_s": float(k)})
                rows.append(actual)
    slab = next(x for x in rows if x["candidate"] == "slab"
                and x["vertical_diffusivity_m2_s"]
                == float(wc["source_default_vertical_diffusivity_m2_s"]))
    for row in rows:
        row["surface_amplitude_relative_to_slab"] = (
            row["surface_amplitude_k_per_w_m2"]
            / slab["surface_amplitude_k_per_w_m2"])

    # N=1 has an analytic integrator response. This catches a timestep, heat
    # capacity or harmonic-convention error independently of vdiffo.
    analytic = 1.0 / (rho * cp * float(wc["total_depth_m"])
                      * (2.0 * np.pi / period))
    analytic_residual = abs(slab["surface_amplitude_k_per_w_m2"] / analytic - 1.0)
    last_overwrite = bool(re.search(
        r"dlayer\s*\(\s*NLEV_OCE\s*\)\s*=\s*mldepth", source, re.IGNORECASE))
    top_only = bool(re.search(
        r"zsst\s*\(:,\s*1\s*\)\s*=\s*zsst\s*\(:,\s*1\s*\).*yheat",
        source, re.IGNORECASE))
    max_closure = max(x["thickness_weighted_closure_residual_m"] for x in rows)
    default_rows = [x for x in rows if x["vertical_diffusivity_m2_s"]
                    == float(wc["source_default_vertical_diffusivity_m2_s"])]
    report = {
        "contract_version": dec["contract_version"],
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "forcing": {
            "orbital_period_days": period / 86400.0,
            "solar_day_hours": solar_day / 3600.0,
            "ocean_steps_per_solar_day": ntspd,
            "timestep_seconds": dt,
            "amplitude_w_m2": 1.0,
        },
        "checks": {
            "single_slab_against_analytic_relative_residual": analytic_residual,
            "maximum_thickness_weighted_closure_residual_m": max_closure,
            "source_overwrites_last_layer_with_mldepth": last_overwrite,
            "source_deposits_all_surface_heat_in_level_1": top_only,
        },
        "responses": rows,
        "default_diffusivity_summary": default_rows,
        "cost": {
            "persistent_real_fields_per_horizontal_cell_per_layer": 2,
            "peak_additional_scratch_bytes_per_horizontal_cell_per_layer": 40,
            "restart_reals_per_horizontal_cell_per_layer": 1,
            "vertical_solver_complexity": "O(horizontal_cells * layers)",
            "interpretation": ("computationally viable at the declared profiles; "
                               "the rejection is physical/interface scope, not cost"),
        },
        "verdict": dec["verdict"],
    }
    if analytic_residual > 1e-3 or max_closure > 1e-12 or not last_overwrite or not top_only:
        raise SystemExit("OCN-1 source/response gate failed")
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n")
    for row in default_rows:
        print(f"{row['candidate']} ({row['profile_interpretation']}): "
              f"{row['surface_amplitude_k_per_w_m2']:.6f} "
              f"K/(W m-2), {row['surface_amplitude_relative_to_slab']:.2f}x slab")
    print(f"verdict: {dec['verdict']['status']}; wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
