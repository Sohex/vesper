#!/usr/bin/env python3
"""What controls the modelled soil's thermal inertia: moisture, or lithology?

    python analysis/soil_thermal_inertia.py

Worldbuilding. Vesper is an invented planet; every number below is a property
of the simulated soil column, of `landmod.f90`, or of a published soil thermal
model applied to them.

`landmod.f90:tands` builds the soil column's heat capacity and conductivity per
layer as

    zcap(:,jlev)  = sicecap*dglac(:)  + soilcap*(1.-dglac(:))
    zdiff(:,jlev) = sicediff*dglac(:) + soildiff*(1.-dglac(:))

so outside glacier ice the whole simulated planet has one thermal conductivity
and one heat capacity. Their thermal inertia sqrt(k*rho*c) is the number the
diurnal and seasonal surface temperature amplitude is inversely proportional
to, and that amplitude reaches evaporation, P minus E, the carve criterion and
the dust emission threshold.

The recorded finding is that the barren classes -- salt crust and playa mud,
about a third of this world's land -- carry roughly a third of the model's
inertia and are therefore damped several times too much, and that the fix is a
per-cell field out of the lithology map. The first half is right. This script
is the test of the second half, and the answer is that lithology is not what
the difference is made of.

## The model, and why it is this one

Johansen's interpolation as Farouki tabulates it, the same formulation
`lawrence_2007` describes blending organic material into. The monograph itself
is not held -- `references/INDEX.md` records that DTIC does not resolve from
this host -- so the relations here come from that secondary statement, which is
enough for the form and not enough to check the coefficients against the
original. Three relations:

    k_dry = (0.135*rho_b + 64.7) / (2700 - 0.947*rho_b)
    k_sat = k_solid^(1-phi) * k_water^phi
    k     = Ke*(k_sat - k_dry) + k_dry

with the Kersten number `Ke` a function of saturation alone, and

    k_solid = k_quartz^q * k_other^(1-q)

carrying the only mineralogical dependence there is. The heat capacity is the
volume-weighted sum over solid and water.

**Read what that says before running it.** The conductivity of DRY soil is a
function of bulk density and nothing else: mineralogy enters only through
`k_sat`, so it is largest where the soil is wet and vanishes where it is dry.
The classes the finding names are the dry ones. That is the shape of the answer
and this script is the size of it.

## What is measured

Two spreads, in the same units, on the same column:

**The moisture spread.** One texture, one mineralogy, saturation swept from air
dry to saturated.

**The lithology spread.** One texture, one saturation, the quartz fraction
swept across the parent materials this world's land is made of -- and reported
at both ends of the moisture range, because the whole question is whether the
two are separable.

The comparison decides which term a per-cell field should carry, and it is
falsifiable: if the lithology spread at the dry end were the larger, the
recorded finding's shape would be right and this would say so.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "analysis" / "soil_thermal_inertia.json"

# Johansen / Farouki, the formulation lawrence_2007 blends organic material
# into. W/m/K unless stated.
K_QUARTZ = 7.7
K_OTHER = 2.0
K_WATER = 0.57
RHO_SOLID = 2700.0              # kg/m3, the mineral grain density the dry
                                # relation is written against
C_SOLID = 2.0e6                 # J/m3/K, mineral grains
C_WATER = 4.18e6                # J/m3/K

# The column this is evaluated on. A medium-textured mineral soil, declared
# rather than fitted: what is being compared is two SPREADS, and both move
# together with the column, so the comparison is insensitive to it. The
# sensitivity is reported by re-running at the bracket ends.
BULK_DENSITY_KG_M3 = 1300.0
BULK_DENSITY_BRACKET = (1100.0, 1600.0)

# Volume fraction of quartz in the parent material, per the rock classes this
# world's land carries. Quartz is the only mineral the model distinguishes,
# because it is the only one whose conductivity is far from the rest.
QUARTZ_FRACTION = {
    "basalt and carbonate": 0.00,
    "evaporite salt crust": 0.00,
    "andesite and granodiorite": 0.15,
    "playa mud and fan fill": 0.30,
    "granite and gneiss": 0.30,
    "schist and melange": 0.30,
    "continental clastics": 0.65,
    "quartzite": 0.95,
}

# Saturation, as a fraction of pore volume. The low end is air dry rather than
# oven dry: a soil in a climate model never reaches zero.
SATURATION_RANGE = (0.02, 1.00)


def kersten(sr: np.ndarray) -> np.ndarray:
    """Johansen's interpolation between the dry and saturated conductivities."""
    return np.clip(np.log10(np.maximum(sr, 1e-6)) + 1.0, 0.0, 1.0)


def conductivity(sr: np.ndarray, quartz: float, rho_b: float) -> np.ndarray:
    phi = 1.0 - rho_b / RHO_SOLID
    k_solid = K_QUARTZ ** quartz * K_OTHER ** (1.0 - quartz)
    k_dry = (0.135 * rho_b + 64.7) / (2700.0 - 0.947 * rho_b)
    k_sat = k_solid ** (1.0 - phi) * K_WATER ** phi
    return kersten(sr) * (k_sat - k_dry) + k_dry


def heat_capacity(sr: np.ndarray, rho_b: float) -> np.ndarray:
    phi = 1.0 - rho_b / RHO_SOLID
    return (1.0 - phi) * C_SOLID + sr * phi * C_WATER


def inertia(sr, quartz: float, rho_b: float):
    return np.sqrt(conductivity(sr, quartz, rho_b) * heat_capacity(sr, rho_b))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bulk-density", type=float, default=BULK_DENSITY_KG_M3)
    ap.add_argument("--json", type=Path, default=REPORT)
    args = ap.parse_args()
    rho_b = args.bulk_density

    # The model's own pair, read from the source rather than copied.
    src = (ROOT / "vendor" / "exoplasim" / "exoplasim" / "plasim" / "src"
           / "landmod.f90").read_text(encoding="utf-8", errors="replace")
    import re
    def declared(name):
        m = re.search(rf"^\s*real\s*::\s*{name}\s*=\s*([-+0-9.eEdD]+)", src,
                      re.IGNORECASE | re.MULTILINE)
        if not m:
            raise SystemExit(f"landmod.f90 no longer declares {name}")
        return float(m.group(1).replace("D", "e").replace("d", "e"))
    soildiff, soilcap = declared("soildiff"), declared("soilcap")
    model_inertia = float(np.sqrt(soildiff * soilcap))

    lo, hi = SATURATION_RANGE
    mid_quartz = QUARTZ_FRACTION["playa mud and fan fill"]

    moisture = {
        "air_dry": float(inertia(lo, mid_quartz, rho_b)),
        "saturated": float(inertia(hi, mid_quartz, rho_b)),
    }
    moisture["ratio"] = moisture["saturated"] / moisture["air_dry"]

    lithology = {}
    for end, sr in (("air_dry", lo), ("saturated", hi)):
        vals = {name: float(inertia(sr, q, rho_b))
                for name, q in QUARTZ_FRACTION.items()}
        lithology[end] = {
            "by_parent_material": {k: round(v, 1) for k, v in sorted(vals.items())},
            "ratio": max(vals.values()) / min(vals.values()),
        }

    # What the model's single pair actually corresponds to, on this column.
    sr_grid = np.linspace(1e-3, 1.0, 20001)
    curve = inertia(sr_grid, mid_quartz, rho_b)
    model_saturation = (float(sr_grid[int(np.argmin(np.abs(curve - model_inertia)))])
                        if curve.max() >= model_inertia else None)

    sens = {}
    for r in BULK_DENSITY_BRACKET:
        sens[str(int(r))] = {
            "air_dry": round(float(inertia(lo, mid_quartz, r)), 1),
            "saturated": round(float(inertia(hi, mid_quartz, r)), 1),
            "moisture_ratio": round(float(inertia(hi, mid_quartz, r)
                                          / inertia(lo, mid_quartz, r)), 2),
        }
    # The THIRD term, and the one a per-cell field could honestly carry at the
    # dry end: bulk density. It is a texture property, which pedology's soil
    # map holds and the lithology map does not, and unlike mineralogy it does
    # not vanish where the soil is dry.
    ends = [inertia(lo, mid_quartz, r) for r in BULK_DENSITY_BRACKET]
    bulk_density_spread = {
        "air_dry_ratio": round(float(max(ends) / min(ends)), 2),
        "note": "bulk density is the one term besides moisture that acts at "
                "the dry end. It is texture, not parent material, and the "
                "component that carries it is pedology",
    }

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "question": "is the modelled soil's thermal inertia controlled by "
                    "moisture or by parent material",
        "model": "Johansen's interpolation as Farouki tabulates it; the "
                 "mineralogy enters only through k_sat, so the dry "
                 "conductivity is a function of bulk density alone",
        "column": {"bulk_density_kg_m3": rho_b,
                   "porosity": round(1.0 - rho_b / RHO_SOLID, 3)},
        "landmod_scalars": {"soildiff_w_m_k": soildiff,
                            "soilcap_j_m3_k": soilcap,
                            "thermal_inertia_j_m2_k_s05": round(model_inertia, 1),
                            "saturation_that_reproduces_it": (
                                round(model_saturation, 3)
                                if model_saturation is not None else None)},
        "moisture_spread": {k: (round(v, 1) if k != "ratio" else round(v, 2))
                            for k, v in moisture.items()},
        "lithology_spread": {
            end: {"by_parent_material": v["by_parent_material"],
                  "ratio": round(v["ratio"], 2)}
            for end, v in lithology.items()},
        "bulk_density_sensitivity": sens,
        "bulk_density_spread": bulk_density_spread,
        "ranking": [
            "moisture, on one parent material and one texture",
            "bulk density at the dry end, which is texture and not lithology",
            "parent-material mineralogy, which is zero at the dry end because "
            "Johansen's dry conductivity carries no mineralogy at all",
        ],
        "units": "thermal inertia in J/m2/K/s^0.5",
    }
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    if model_saturation is not None:
        print(f"landmod's single pair is a thermal inertia of "
              f"{model_inertia:.0f}, which this column reaches at a saturation "
              f"of {model_saturation:.2f}")
    else:
        print(f"landmod's single pair is a thermal inertia of "
              f"{model_inertia:.0f}, which this column does not reach even "
              f"saturated: it is a wetter or denser soil than this one")
    print(f"moisture, one parent material:  {moisture['air_dry']:.0f} air dry "
          f"to {moisture['saturated']:.0f} saturated, a factor "
          f"{moisture['ratio']:.2f}")
    for end in ("air_dry", "saturated"):
        v = lithology[end]["by_parent_material"]
        print(f"parent material at {end:10s}: {min(v.values()):.0f} to "
              f"{max(v.values()):.0f}, a factor {lithology[end]['ratio']:.2f}")
    print(f"bulk density at air dry:        a factor "
          f"{bulk_density_spread['air_dry_ratio']:.2f}, over "
          f"{BULK_DENSITY_BRACKET[0]:.0f} to {BULK_DENSITY_BRACKET[1]:.0f} kg/m3")
    print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
