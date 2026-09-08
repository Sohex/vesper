#!/usr/bin/env python3
"""Close the model's INTERNAL energy budget against PlaSim's own 28 terms.

`close_state_energy.py` asks whether the planet stores what the top of the
atmosphere reports. When it does not -- and on this planet it does not, by about
half a watt -- that script cannot say where the missing watt is made, because
every flux diagnostic it can reach is derived from the same radiative profile.
This one asks the next question: does the model conserve energy internally?

It uses three identities, each with a right answer and each able to fail.

1. GRIDPOINT PHYSICS CLOSURE. `denergy04` is the column enthalpy tendency built
   from `dtdt`, the total gridpoint temperature tendency, at the end of
   `gridpointd` (plasim.f90:3462). The parameterisations each book their own
   heating: 6 miscmod, 7 sensible flux, 8 vertical diffusion of heat, 9 longwave,
   10 shortwave, 11 large-scale condensation, 12 convection, 13 and 14 mkrain,
   21 surface-stress dissipation, 22 vertical-diffusion dissipation. Their sum
   must equal term 4. A gap is a heating applied to the atmosphere that no term
   books.

2. ADIABATIC ENERGY CONSERVATION. `denergy26` is the column enthalpy change
   across `spectrala`, the adiabatic spectral step, and `denergy27` is MINUS the
   kinetic energy change across the same step (plasim.f90:3125-3145). The
   adiabatic dynamics moves energy between enthalpy and kinetic energy and
   creates none, so `26 - 27` must be zero. It is not, and that is this script's
   result: the spectral core is a net energy source.

3. KINETIC ENERGY STEADY STATE. In a settled run the kinetic energy is flat, so
   the adiabatic generation `-denergy27` must equal the frictional dissipation
   returned as heat: 21 (surface stress) + 22 (vertical momentum diffusion) +
   23 (Rayleigh friction) + 25 (biharmonic diffusion). This is what makes
   identity 2 interpretable rather than a mis-pairing: if term 27 were not the
   adiabatic generation, this would not close.

WHY THE WINDOW MATTERS. `denergy` is written UNACCUMULATED (outmod.f90:1162)
while every flux it is compared against is a time mean, so under `NLOWIO = 1` the
28 terms are instantaneous snapshots and disagree with the fluxes by several
percent -- the sensible-heat identity, which is exact algebra, was out by 1.6
W/m2 that way. This script therefore REFUSES a window whose segments did not run
with low I/O off, rather than returning numbers that look like measurements.

Usage:

    python exoplasim/scripts/close_term_energy.py <run_dir> --first 67 --last 76

Evidence and interpretation: exoplasim/notes/water-and-energy-closure.md.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from _paths import ANALYSIS
import gridding  # noqa: E402  from lib/, via _paths: the one Gaussian quadrature

# The terms that carry a heating applied inside `gridpointd`. 15 and 16 are the
# moisture counterparts of 11 and 12 rather than second heatings, 17 to 20 are
# nested sub-splits of the longwave, 1 and 28 are a column enthalpy and a column
# mean temperature rather than fluxes, and 2, 3, 5, 24, 26, 27 belong to the
# spectral half. Summing any of those with the rest double-counts.
GRIDPOINT_TERMS = [6, 7, 8, 9, 10, 11, 12, 13, 14, 21, 22]
DISSIPATION_TERMS = [21, 22, 23, 25]
TERMS = [f"denergy{i:02d}" for i in range(1, 29)]
FLUXES = ["ntr", "rst", "rsut", "rlut", "rss", "rls", "hfss", "hfls", "hfns"]


def global_mean(field: np.ndarray, weights: np.ndarray) -> float:
    """Gaussian-weighted mean over latitude, longitude and the output bins.

    The bins are equal-length time averages, so their plain mean is the orbit
    mean and no time weighting is needed.
    """
    return float(np.sum(field * weights[None, :, None])
                 / (2.0 * field.shape[-1] * field.shape[0]))


def low_io_off(manifest: dict, first: int, last: int) -> tuple[bool, list[str]]:
    """Did every segment covering the window run with NLOWIO = 0?

    Fails closed: a segment that does not say is treated as low I/O ON, because
    every run made before 2026-08-17 was, and those records carry no key.
    """
    problems = []
    covered = set()
    for segment in manifest.get("segments", []):
        start = int(segment.get("start_year_index", 0))
        end = int(segment.get("end_year_index", -1))
        orbits = set(range(start, end + 1)) & set(range(first, last + 1))
        if not orbits:
            continue
        covered |= orbits
        if segment.get("low_io", True):
            problems.append(f"orbits {start}-{end} ran with low I/O on "
                            f"({segment.get('purpose', 'no purpose recorded')})")
    missing = sorted(set(range(first, last + 1)) - covered)
    if missing:
        problems.append(f"no segment record covers orbits {missing}")
    return not problems, problems


def term_energy(run_dir: Path, first: int, last: int) -> dict:
    manifest_path = run_dir / "run_manifest.json"
    manifest = (json.loads(manifest_path.read_text(encoding="utf-8"))
                if manifest_path.is_file() else {})

    per_orbit = []
    for index in range(first, last + 1):
        path = run_dir / f"MOST.{index:05d}.nc"
        if not path.is_file():
            raise SystemExit(f"missing annual output {path}")
        with Dataset(path) as nc:
            missing = [n for n in TERMS + FLUXES if n not in nc.variables]
            if missing:
                raise SystemExit(
                    f"{path} does not carry {missing[0]}; this window was run "
                    "without model.energy_diagnostics")
            weights = gridding.gaussian_row_weights(
                np.asarray(nc["lat"][:], dtype=float),
                what=f"{path}'s latitude axis")
            record = {"orbit": index}
            for name in TERMS + FLUXES:
                record[name] = global_mean(
                    np.asarray(nc[name][:], dtype=float), weights)
        per_orbit.append(record)

    s = {name: np.array([r[name] for r in per_orbit])
         for name in per_orbit[0] if name != "orbit"}

    def stat(series):
        return {"mean_w_m2": float(series.mean()),
                "spread_w_m2": float(series.std(ddof=1))}

    gridpoint = sum(s[f"denergy{i:02d}"] for i in GRIDPOINT_TERMS)
    generation = -s["denergy27"]
    dissipation = sum(s[f"denergy{i:02d}"] for i in DISSIPATION_TERMS)
    adiabatic = s["denergy26"] - s["denergy27"]

    # The atmosphere's total energy budget, closed against the fluxes. The
    # latent asymmetry is real rather than an error: convection books the latent
    # heat of sublimation below freezing, so the atmosphere receives the heat of
    # fusion that the surface pays back when the snow melts, which is the melt
    # term inside `hfns`.
    latent_asymmetry = (s["denergy11"] + s["denergy12"] + s["denergy14"]
                        + s["hfls"])
    radiative = s["ntr"] - s["rss"] - s["rls"]
    atmosphere_flux = radiative - s["hfss"] - s["hfls"]

    return {
        "manifest": manifest,
        "per_orbit": per_orbit,
        "identities": {
            "gridpoint_physics_closure": {
                "statement": "sum(denergy 6,7,8,9,10,11,12,13,14,21,22) = denergy04",
                "right_answer_w_m2": 0.0,
                "sum_of_terms": stat(gridpoint),
                "denergy04": stat(s["denergy04"]),
                "residual": stat(s["denergy04"] - gridpoint),
            },
            "adiabatic_energy_conservation": {
                "statement": "denergy26 - denergy27 = 0 across the adiabatic "
                             "spectral step (enthalpy plus kinetic energy)",
                "right_answer_w_m2": 0.0,
                "enthalpy_tendency": stat(s["denergy26"]),
                "minus_kinetic_tendency": stat(s["denergy27"]),
                "residual": stat(adiabatic),
            },
            "kinetic_energy_steady_state": {
                "statement": "-denergy27 = denergy 21+22+23+25 in a settled run",
                "right_answer_w_m2": 0.0,
                "adiabatic_generation": stat(generation),
                "frictional_heat_returned": stat(dissipation),
                "residual": stat(generation - dissipation),
            },
        },
        "atmosphere_budget_w_m2": {
            "radiative_convergence": stat(radiative),
            "sensible_from_surface": stat(-s["hfss"]),
            "latent_from_surface": stat(-s["hfls"]),
            "latent_asymmetry": stat(latent_asymmetry),
            "adiabatic_non_conservation": stat(adiabatic),
            "spectral_diffusion_of_heat": stat(s["denergy24"]),
            "flux_route_total": stat(atmosphere_flux),
            "closed_total": stat(atmosphere_flux + latent_asymmetry + adiabatic
                                 + s["denergy24"]),
            "column_enthalpy_j_m2": float(s["denergy01"].mean()),
        },
        "terms_w_m2": {name: stat(s[name]) for name in TERMS},
        "fluxes_w_m2": {name: stat(s[name]) for name in FLUXES},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--first", type=int, required=True)
    parser.add_argument("--last", type=int, required=True)
    parser.add_argument("--output", type=Path, default=ANALYSIS / "energy_terms")
    parser.add_argument("--allow-low-io", action="store_true",
                        help="compute anyway on a low-I/O window. The 28 terms "
                             "are then instantaneous snapshots and every number "
                             "below is unquotable; see the module docstring")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()

    manifest_path = run_dir / "run_manifest.json"
    manifest = (json.loads(manifest_path.read_text(encoding="utf-8"))
                if manifest_path.is_file() else {})
    clean, problems = low_io_off(manifest, args.first, args.last)
    if not clean and not args.allow_low_io:
        raise SystemExit(
            "the 28 energy terms are written unaccumulated, so they are only "
            "comparable with the fluxes on a window run with NLOWIO = 0:\n  "
            + "\n  ".join(problems)
            + "\npass --allow-low-io to compute anyway and label the result "
              "as a snapshot reading.")

    result = term_energy(run_dir, args.first, args.last)
    report = {
        "schema_version": 1,
        "generator": "exoplasim/scripts/close_term_energy.py",
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "run_dir": str(run_dir),
        "run_id": manifest.get("run_id"),
        "source_build": manifest.get("source_build"),
        "config_sha256": manifest.get("config_sha256"),
        "geography": manifest.get("physical", {}).get("geography"),
        "window": {"first_orbit": args.first, "last_orbit": args.last,
                   "orbits": len(result["per_orbit"]),
                   "low_io_off": clean, "low_io_problems": problems},
        "identities": result["identities"],
        "atmosphere_budget_w_m2": result["atmosphere_budget_w_m2"],
        "terms_w_m2": result["terms_w_m2"],
        "fluxes_w_m2": result["fluxes_w_m2"],
        "executable_sha256": manifest.get("executable", {}).get("sha256"),
        "software": manifest.get("software"),
        "per_orbit": result["per_orbit"],
    }

    args.output.mkdir(parents=True, exist_ok=True)
    # Named by the run AND the window, for the same reason close_state_energy.py
    # is: the window is a parameter of the measurement, not a detail of it.
    path = args.output / f"{run_dir.name}_term_energy_{args.first}-{args.last}.json"
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({"identities": report["identities"],
                      "atmosphere_budget_w_m2": report["atmosphere_budget_w_m2"],
                      "report": str(path)}, indent=2))


if __name__ == "__main__":
    main()
