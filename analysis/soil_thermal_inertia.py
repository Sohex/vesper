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

Johansen's interpolation as Farouki (1981) tabulates it, the same formulation
`lawrence_2007` describes blending organic material into. **The coefficients are
checked against the monograph itself**, Section 7.11 and Table 24 on pages 112
and 113, rather than through that secondary statement: a restatement is enough
for the FORM of a relation and not for its coefficients, and the whole finding
below turns on where mineralogy enters, which is a coefficient question.

    k_dry = (0.135*rho_b + 64.7) / (2700 - 0.947*rho_b)     +/- 20 per cent
    k_sat = k_solid^(1-phi) * k_water^phi
    k     = Ke*(k_sat - k_dry) + k_dry

with the Kersten number `Ke` a function of saturation and texture, and

    k_solid = k_quartz^q * k_other^(1-q)

carrying the only mineralogical dependence there is. The heat capacity is the
volume-weighted sum over solid and water.

**Two of Farouki's relations are TEXTURE-BRANCHED and the branch is taken here
by the declared column rather than left implicit.** `Ke` is `log10(Sr) + 1` for
a fine soil and `0.7*log10(Sr) + 1` for a coarse one, and `k_other` is 2.0 W/m/K
in general but 3.0 W/m/K for a COARSE soil whose quartz fraction is below 0.20.
Farouki's cut is his Figure 160's: a soil with more than 5 per cent of material
finer than 2 micrometres is fine. The column below is a medium-textured mineral
soil and is fine by that rule, so the fine `Ke` and `k_other = 2.0` are the rows
that apply to it. Both alternatives are evaluated anyway, because a coarse
low-quartz soil is the one case where the mineralogical term is largest and this
finding is about how large that term is.

**Two things read out of the monograph that a later reader will otherwise
re-derive wrongly.** Table 24 misprints two coefficients the facing body text
gives correctly -- `0.137` for `0.135` in the dry relation's numerator, and
`0.39` for `0.039` in the crushed-rock relation -- so the table is not the row to
copy from. And Johansen's `k_water = 0.57` is water at about 5 C on Farouki's own
Table 12, which is part of the method as specified and is a few per cent low for
a warm column; it enters only through `k_sat`.

**The dry end of the sweep is outside the method's stated range and is reported
twice for that reason.** Farouki gives the fine relation for `Sr > 0.1` and the
coarse one for `Sr > 0.05`, and says that below about 0.1 no method predicts well
except Van Rooyen's. Air dry is below that floor, so the answer there is set by
the clip rather than by the relation. Both the air-dry point and the floor are
reported, and the finding is the same at either.

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

# Johansen's method as Farouki (1981) tabulates it, Section 7.11 and Table 24,
# read from the monograph. W/m/K unless stated.
K_QUARTZ = 7.7                  # Horai's mean at about 25 C. Farouki notes
                                # later Norwegian measurements suggesting nearer
                                # 10 and wants independent confirmation first.
K_OTHER = 2.0                   # feldspar or mica, on Horai's measurements
K_OTHER_COARSE_LOW_QUARTZ = 3.0 # Johansen's own branch, for a COARSE soil at
                                # q < 0.20, standing in for the probable mineral
                                # composition of such soils
LOW_QUARTZ_BRANCH = 0.20        # the quartz fraction the branch is taken below
K_WATER = 0.57                  # water at about 5 C on Farouki's Table 12; part
                                # of the method as specified
RHO_SOLID = 2700.0              # kg/m3, the mineral grain density the dry
                                # relation is written against
C_SOLID = 2.0e6                 # J/m3/K, mineral grains
C_WATER = 4.18e6                # J/m3/K

# Farouki's stated validity floors on the degree of saturation, below which he
# says no method predicts well except Van Rooyen's. They are the reason the dry
# end is reported at two points.
SR_FLOOR_FINE = 0.10
SR_FLOOR_COARSE = 0.05

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


def kersten(sr: np.ndarray, texture: str = "fine") -> np.ndarray:
    """Johansen's interpolation between the dry and saturated conductivities.

    Farouki Table 24 rows (d) and (e): `0.7*log10(Sr) + 1` for a coarse unfrozen
    soil, `log10(Sr) + 1` for a fine one. The clip has no counterpart in Farouki
    -- the relations simply are not defined below their floors -- so below
    `SR_FLOOR_*` it is the clip that sets the answer and `k` collapses to
    `k_dry`. That is qualitatively the right behaviour and outside the stated
    range, which is why the dry end is reported at the floor as well.
    """
    slope = 0.7 if texture == "coarse" else 1.0
    return np.clip(slope * np.log10(np.maximum(sr, 1e-6)) + 1.0, 0.0, 1.0)


def solid_conductivity(quartz: float, texture: str = "fine") -> float:
    """Farouki Table 24 rows (i) and (j), and the branch is the load-bearing one.

    `k_other` is 2.0 W/m/K in general and 3.0 W/m/K for a COARSE soil below
    `LOW_QUARTZ_BRANCH` quartz. Farouki's 7.13.5 says why it matters: where the
    quartz content is high `k_other` barely moves `k_solid`, and where it is low
    it dominates. The parent materials this world's land carries are mostly low
    in quartz, so the branch would be the difference between two answers for
    them -- if they were coarse. It changes nothing at the dry end, because
    `k_dry` carries no mineralogy at all.
    """
    k_other = (K_OTHER_COARSE_LOW_QUARTZ
               if texture == "coarse" and quartz < LOW_QUARTZ_BRANCH
               else K_OTHER)
    return float(K_QUARTZ ** quartz * k_other ** (1.0 - quartz))


def conductivity(sr: np.ndarray, quartz: float, rho_b: float,
                 texture: str = "fine") -> np.ndarray:
    phi = 1.0 - rho_b / RHO_SOLID
    k_solid = solid_conductivity(quartz, texture)
    k_dry = (0.135 * rho_b + 64.7) / (2700.0 - 0.947 * rho_b)
    k_sat = k_solid ** (1.0 - phi) * K_WATER ** phi
    return kersten(sr, texture) * (k_sat - k_dry) + k_dry


def heat_capacity(sr: np.ndarray, rho_b: float) -> np.ndarray:
    phi = 1.0 - rho_b / RHO_SOLID
    return (1.0 - phi) * C_SOLID + sr * phi * C_WATER


def inertia(sr, quartz: float, rho_b: float, texture: str = "fine"):
    return np.sqrt(conductivity(sr, quartz, rho_b, texture)
                   * heat_capacity(sr, rho_b))


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
        # Farouki's stated floor for the fine relation this column takes. Below
        # it the clip and not the relation sets the conductivity, so the sweep
        # is reported from here as well as from air dry, and the ratio that
        # stays inside the method's own range is the second one.
        "validity_floor": float(inertia(SR_FLOOR_FINE, mid_quartz, rho_b)),
    }
    moisture["ratio"] = moisture["saturated"] / moisture["air_dry"]
    moisture["ratio_within_stated_validity"] = (
        moisture["saturated"] / moisture["validity_floor"])

    lithology = {}
    for end, sr in (("air_dry", lo), ("validity_floor", SR_FLOOR_FINE),
                    ("saturated", hi)):
        vals = {name: float(inertia(sr, q, rho_b))
                for name, q in QUARTZ_FRACTION.items()}
        lithology[end] = {
            "by_parent_material": {k: round(v, 1) for k, v in sorted(vals.items())},
            "ratio": max(vals.values()) / min(vals.values()),
        }

    # THE ONE BRANCH THAT WOULD MOVE THE LITHOLOGY TERM, evaluated rather than
    # argued. Farouki's k_other rises to 3.0 W/m/K for a COARSE soil below 0.20
    # quartz, and most of this world's parent materials are below it. The
    # declared column is fine by his own 5-per-cent-below-2-micrometre rule, so
    # this is what the finding would look like on a coarse one instead.
    coarse = {}
    for end, sr in (("air_dry", lo), ("validity_floor", SR_FLOOR_COARSE),
                    ("saturated", hi)):
        vals = {name: float(inertia(sr, q, rho_b, "coarse"))
                for name, q in QUARTZ_FRACTION.items()}
        coarse[end] = {
            "by_parent_material": {k: round(v, 1) for k, v in sorted(vals.items())},
            "ratio": round(max(vals.values()) / min(vals.values()), 2),
        }
    coarse_low_quartz_lift = {
        name: round(float(inertia(hi, q, rho_b, "coarse")
                          / inertia(hi, q, rho_b, "fine")), 3)
        for name, q in QUARTZ_FRACTION.items() if q < LOW_QUARTZ_BRANCH
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
        "farouki_texture_branches": {
            "declared_column_is": "fine, by Farouki's Figure 160 rule that a "
                                  "soil with more than 5 per cent finer than 2 "
                                  "micrometres is fine. So the rows that apply "
                                  "are Ke = log10(Sr) + 1 and k_other = 2.0.",
            "kersten_coarse_alternative": "0.7*log10(Sr) + 1, which reaches "
                                          "k_dry at a HIGHER saturation and "
                                          "therefore widens the moisture "
                                          "spread rather than narrowing it",
            "k_other_coarse_low_quartz_w_m_k": K_OTHER_COARSE_LOW_QUARTZ,
            "low_quartz_branch_at": LOW_QUARTZ_BRANCH,
            "on_a_coarse_column": coarse,
            "saturated_inertia_lift_on_the_low_quartz_classes":
                coarse_low_quartz_lift,
            "why_it_does_not_change_the_finding": "the branch enters k_solid, "
                                                  "k_solid enters k_sat, and "
                                                  "k_sat is multiplied by a "
                                                  "Kersten number that is zero "
                                                  "at and below the dry floor. "
                                                  "So it lifts the WET end of "
                                                  "the lithology spread for "
                                                  "the low-quartz classes and "
                                                  "leaves the dry end at "
                                                  "exactly 1.00, which is the "
                                                  "row the finding rests on.",
        },
        "validity": {
            "kersten_floor_fine": SR_FLOOR_FINE,
            "kersten_floor_coarse": SR_FLOOR_COARSE,
            "dry_relation_accuracy": "+/- 20 per cent, stated with the "
                                     "equation on Farouki p. 112",
            "below_the_floor": "Farouki 7.13.1: below a saturation of about "
                               "0.1 no method predicts well except Van "
                               "Rooyen's. The air-dry point is below that "
                               "floor, so the clip and not the relation sets "
                               "the conductivity there; the floor row is the "
                               "one inside the method's stated range and the "
                               "ranking is the same at either.",
            "coefficients_checked_against": "Farouki (1981) CRREL Monograph "
                                            "81-1, Section 7.11 and Table 24, "
                                            "pp. 112-113. Table 24 misprints "
                                            "0.137 for 0.135 in the dry "
                                            "numerator and 0.39 for 0.039 in "
                                            "the crushed-rock relation; the "
                                            "body text is the row to copy.",
            "crushed_rock_not_used": "Farouki gives a separate dry relation "
                                     "for crushed rock, 0.039*n**-2.2 +/- 25 "
                                     "per cent, a function of porosity alone, "
                                     "and says Johansen's method does NOT "
                                     "apply well to dry crushed rocks. This "
                                     "column is a natural soil and takes the "
                                     "natural-soil relation.",
        },
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
    print(f"  inside Farouki's stated range: "
          f"{moisture['validity_floor']:.0f} at Sr = {SR_FLOOR_FINE} to "
          f"{moisture['saturated']:.0f}, a factor "
          f"{moisture['ratio_within_stated_validity']:.2f}")
    for end in ("air_dry", "validity_floor", "saturated"):
        v = lithology[end]["by_parent_material"]
        print(f"parent material at {end:10s}: {min(v.values()):.0f} to "
              f"{max(v.values()):.0f}, a factor {lithology[end]['ratio']:.2f}")
    print(f"bulk density at air dry:        a factor "
          f"{bulk_density_spread['air_dry_ratio']:.2f}, over "
          f"{BULK_DENSITY_BRACKET[0]:.0f} to {BULK_DENSITY_BRACKET[1]:.0f} kg/m3")
    print(f"on a COARSE column instead, where Farouki's k_other = "
          f"{K_OTHER_COARSE_LOW_QUARTZ} branch applies below "
          f"{LOW_QUARTZ_BRANCH} quartz:")
    print(f"  parent material at saturated : a factor "
          f"{coarse['saturated']['ratio']:.2f}")
    print(f"  parent material at air_dry   : a factor "
          f"{coarse['air_dry']['ratio']:.2f}")
    if coarse_low_quartz_lift:
        worst = max(coarse_low_quartz_lift.values())
        print(f"  the low-quartz classes gain up to {worst:.3f} of their "
              f"saturated inertia, and nothing at the dry end")
    print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
