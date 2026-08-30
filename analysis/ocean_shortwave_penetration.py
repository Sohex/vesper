#!/usr/bin/env python3
"""OCN-6 clear-water shortwave deposition cost and bound."""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))
import stellar

DECLARATION = ROOT / "ocean/config/shortwave_penetration.yaml"
OUTPUT = ROOT / "analysis/ocean_shortwave_penetration.json"


def deposition(wavelength_um: np.ndarray, extinction_k: np.ndarray,
               layer_depths_m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-wavelength Beer-Lambert absorption fractions and transmission."""
    lam = np.asarray(wavelength_um, dtype=float)
    k = np.asarray(extinction_k, dtype=float)
    depths = np.asarray(layer_depths_m, dtype=float)
    if lam.shape != k.shape or np.any(lam <= 0) or np.any(k < 0):
        raise ValueError("positive wavelengths and nonnegative peer k are required")
    if depths.ndim != 1 or np.any(depths <= 0):
        raise ValueError("layer depths must be a positive vector")
    absorption_m1 = 4.0 * np.pi * k / (lam * 1.0e-6)
    edges = np.concatenate([[0.0], np.cumsum(depths)])
    transmission = np.exp(-absorption_m1[:, None] * edges[None, :])
    absorbed = transmission[:, :-1] - transmission[:, 1:]
    return absorbed, transmission[:, -1]


def integrate_spectrum(wavelength_um: np.ndarray, spectral_flux: np.ndarray,
                       absorbed: np.ndarray, transmitted: np.ndarray) -> dict:
    flux = np.asarray(spectral_flux, dtype=float)
    norm = float(np.trapezoid(flux, wavelength_um))
    layers = np.array([
        np.trapezoid(flux * absorbed[:, j], wavelength_um) / norm
        for j in range(absorbed.shape[1])])
    below = float(np.trapezoid(flux * transmitted, wavelength_um) / norm)
    return {"absorbed_fraction_by_layer": layers.tolist(),
            "transmitted_below_column": below,
            "closure_residual": float(abs(layers.sum() + below - 1.0)),
            "absorbed_top_1m": None}


def main() -> None:
    dec = yaml.safe_load(DECLARATION.read_text())
    optics_path = ROOT / dec["optical_source"]
    optics = np.loadtxt(optics_path)
    _, high = stellar.spectrum_paths()
    wavelength_m, k_flux = stellar.read_hires(high)
    stellar.assert_model_grid(wavelength_m)
    wavelength_um = wavelength_m * 1.0e6
    keep = wavelength_um >= stellar.MIN_WAVELENGTH_NM * 1e-3
    wavelength_um, k_flux = wavelength_um[keep], k_flux[keep]
    # Hale and Querry ends at 4 um in the held extract. At longer wavelengths
    # water is more opaque, so assigning that small spectral tail to level 1 is
    # the conservative exact limit rather than extrapolating an absorption line.
    k_interp = np.interp(np.minimum(wavelength_um, 4.0), optics[:, 0], optics[:, 1])
    layers = np.asarray(dec["hypothetical_layers_m"], dtype=float)
    absorbed, transmitted = deposition(wavelength_um, k_interp, layers)
    beyond = wavelength_um > 4.0
    absorbed[beyond] = 0.0
    absorbed[beyond, 0] = 1.0
    transmitted[beyond] = 0.0

    kstar = integrate_spectrum(wavelength_um, k_flux, absorbed, transmitted)
    # Solar comparison on the IDENTICAL grid and truncation, so only spectral
    # shape changes. Arbitrary Planck units cancel in the normalisation.
    solar_flux = stellar._planck(wavelength_m[keep], 5772.0)
    solar = integrate_spectrum(wavelength_um, solar_flux, absorbed, transmitted)

    # Separate top-metre diagnostic, not a special model layer.
    a1, t1 = deposition(wavelength_um, k_interp, np.array([1.0]))
    a1[beyond, 0], t1[beyond] = 1.0, 0.0
    for row, flux in ((kstar, k_flux), (solar, solar_flux)):
        norm = float(np.trapezoid(flux, wavelength_um))
        row["absorbed_top_1m"] = float(
            np.trapezoid(flux * a1[:, 0], wavelength_um) / norm)

    # Reduction identity: infinite opacity deposits exactly the entire beam in
    # the surface layer, which is oceanmod's existing behaviour.
    reduction_absorbed, reduction_transmitted = deposition(
        np.array([0.5, 1.0]), np.array([1e6, 1e6]), layers)
    reduction_residual = max(
        float(np.max(np.abs(reduction_absorbed[:, 0] - 1.0))),
        float(np.max(np.abs(reduction_absorbed[:, 1:]))),
        float(np.max(np.abs(reduction_transmitted))))

    source = (ROOT / "vendor/exoplasim/exoplasim/plasim/src/icemod.f90").read_text()
    report = {
        "contract_version": dec["contract_version"],
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sources": {
            "optics": str(optics_path.relative_to(ROOT)),
            "optics_sha256": hashlib.sha256(optics_path.read_bytes()).hexdigest(),
            "stellar_spectrum": str(high.relative_to(ROOT)),
            "stellar_sha256": hashlib.sha256(high.read_bytes()).hexdigest(),
        },
        "layers_m": layers.tolist(),
        "k25v_clear_water": kstar,
        "solar_5772k_clear_water": solar,
        "k25v_over_solar_top_1m_absorption": (
            kstar["absorbed_top_1m"] / solar["absorbed_top_1m"]),
        "checks": {
            "surface_deposition_reduction_residual": reduction_residual,
            "maximum_energy_closure_residual": max(
                kstar["closure_residual"], solar["closure_residual"]),
            "icemod_already_receives_separate_shortwave": "xswfl(:)=pswfl(:)" in source,
        },
        "implementation_cost": {
            "new_spectral_radiation_call": False,
            "new_horizontal_prognostic_fields_per_layer": 0,
            "required_interface_change": (
                "accumulate the already available xswfl separately from total xcflux, "
                "pass net shortwave into oceanstep, and distribute it by layer"),
            "operator_complexity": "O(horizontal_ocean_cells * layers) per ocean coupling",
        },
        "excluded": {
            "pigment_and_particle_absorption": "pending OCN-14",
            "surface_reflection": "owned by OCN-7",
        },
        "verdict": dec["verdict"],
    }
    if reduction_residual > 1e-14 or report["checks"]["maximum_energy_closure_residual"] > 1e-12:
        raise SystemExit("OCN-6 deposition gate failed")
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"k25v top 1 m absorption: {100*kstar['absorbed_top_1m']:.2f}%")
    print(f"solar top 1 m absorption: {100*solar['absorbed_top_1m']:.2f}%")
    print(f"k25v below 50 m: {100*kstar['transmitted_below_column']:.2f}%")
    print(f"wrote {OUTPUT.relative_to(ROOT)}; verdict {dec['verdict']['adopted']=}")


if __name__ == "__main__":
    main()
