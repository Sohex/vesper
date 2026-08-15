#!/usr/bin/env python3
"""Compare equilibrated climatology reports across the stellar-flux sweep."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

_MPL_CACHE = Path("/tmp/world-matplotlib-cache")
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from _paths import ANALYSIS


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case", action="append", nargs=2, metavar=("FLUX", "REPORT"), required=True,
        help="Repeat for each case: flux ratio and baseline_climate_report.json",
    )
    parser.add_argument("--output", type=Path, default=ANALYSIS / "sweep")
    args = parser.parse_args()

    cases = []
    for flux_text, report_text in args.case:
        report_path = Path(report_text).resolve()
        report = json.loads(report_path.read_text(encoding="utf-8"))
        cases.append(
            {
                "flux_earth": float(flux_text),
                "report": str(report_path),
                "label": report.get("label"),
                "year_indices": report["year_indices"],
                "global_metrics": report["global_metrics"],
                "koppen_land_area_fractions": report["koppen_land_area_fractions"],
                "biome_land_area_fractions": report["biome_land_area_fractions"],
            }
        )
    cases.sort(key=lambda item: item["flux_earth"])
    flux = np.array([case["flux_earth"] for case in cases])
    metrics = [case["global_metrics"] for case in cases]

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    combined = {"cases": cases}
    (output / "flux_sweep_report.json").write_text(
        json.dumps(combined, indent=2) + "\n", encoding="utf-8"
    )

    panels = [
        ("surface_temperature_k", "Global mean surface temperature", "K"),
        ("precipitation_mm_day", "Global mean precipitation", "mm day⁻¹"),
        ("planetary_sea_ice_fraction", "Planetary sea-ice fraction", "fraction"),
        ("toa_net_radiation_w_m2", "TOA net radiation", "W m⁻²"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    for ax, (key, title, units) in zip(axes.ravel(), panels):
        values = np.array([metric[key] for metric in metrics])
        ax.plot(flux, values, marker="o", linewidth=1.5)
        for x, y in zip(flux, values):
            ax.annotate(f"{y:.3g}", (x, y), xytext=(0, 7), textcoords="offset points", ha="center")
        if key == "toa_net_radiation_w_m2":
            ax.axhline(0, color="black", linewidth=0.7)
        ax.set_title(title)
        ax.set_xlabel("incident stellar flux / modern Earth")
        ax.set_ylabel(units)
        ax.grid(alpha=0.25)
    fig.suptitle("T42 stellar-flux sensitivity climatologies")
    fig.savefig(output / "flux_sweep_summary.png", dpi=180)
    plt.close(fig)
    print(json.dumps(combined, indent=2))


if __name__ == "__main__":
    main()
