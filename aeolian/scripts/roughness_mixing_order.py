#!/usr/bin/env python3
"""Does the ORDER of the aeolian roughness mixture change the emission? SPAT-7.

    python aeolian/scripts/roughness_mixing_order.py

Vesper is a simulated super-Earth. Every surface below is a modelled patch of
that planet's erodible ground, and the field measurements the endmembers rest on
are terrestrial and are what those modelled surfaces are grounded on.

## The question

`dust.yaml` collapses a patchwork of erodible surfaces to ONE roughness per
class and then one per cell, and it does it as a geometric mean -- the mean of
`ln z0`. Two different laws then read that single number, and they are not
nonlinear in the same way, so one reduction cannot be right for both by
assumption. This measures which it is right for.

**MB95's drag partition is AFFINE in `ln z0`.** `feff = 1 - ln(z0/z0s) / D`,
where `D = ln(a (X/z0s)^b)` contains no `z0` at all. So the geometric mean
reproduces the area-mean drag efficiency EXACTLY wherever the clip does not
bind, and the first arm here is that identity rather than a comparison. It can
fail: a clip that binds on any endmember breaks the affinity, and the arm
reports whether it does.

**The friction velocity is not.** `u* = k U / ln(z_ref/z0)` goes as `1/L` with
`L = ln(z_ref/z0)`, and emission goes as roughly `u*` cubed above a threshold.
`1/L^p` is convex in `L` for every positive `p`, so the area mean of `u*^p` is
at or above `u*^p` at the mean `L`, and the geometric mean UNDERSTATES it. That
sign is pre-registrable and one-signed, for the same reason the regolith depth
law's is: a single convex function of one variable does not change curvature.

`dust.yaml` gives the geometric mean the reason "the drag goes as `1/ln(z/z0)`,
so it is `ln(z0)` that averages". That is not the reason and the distinction is
not pedantic: the same sentence would justify a geometric mean for the exchange
coefficient, where `ce = k^2 / ln(z_ref/z0)^2` and
`notes/audits/nonlinear-spatial-reductions.md` section 3 measures a land-mean
error of 8 to 10 per cent from doing exactly that.

## The bar

The instrument each class already reports itself against: its own declared
`bracket` in `aeolian_z0_by_class_m`, which is the level uncertainty the
measurement sets carry. Converted into the arm's own units -- a factor on
`u*^p` -- that is the ignorance the step declares before any question of
operation order arises, and a gap inside it cannot move a decision the bracket
does not already move.

## What is measured, and what bounds what

The mixture inside `playa_clastic` is the WORST CASE available, because its
three landform bands span a factor of twenty in `z0` where the two erodible
CLASSES span under three. The class-level arm runs over a declared share sweep
rather than per-cell shares, which needs no export: if the worst case is inside
the instrument, the narrower mixture is too, and the arm reports the spread each
mixture carries so the bound is checkable rather than asserted.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[2]
sys.path.insert(0, str(ROOT / "lib"))

import numpy as np  # noqa: E402
import yaml  # noqa: E402

import lapse  # noqa: E402
import paths as paths_lib  # noqa: E402
import provenance  # noqa: E402

DUST_CONFIG = ROOT / "aeolian" / "config" / "dust.yaml"
MIX = ROOT / "aeolian" / "analysis" / "playa_roughness_mix.json"
OUT = ROOT / "aeolian" / "analysis" / "roughness_mixing_order.json"

# Powers of the friction velocity the arm reports. Emission is not a pure cube
# -- `build_dust.py:emission_over_weibull` carries `(u^2 - u_t^2)` times a
# power-law ratio -- so the three bracket it rather than pretending to be it.
POWERS = (1, 2, 3)

# The air temperature the reference height is taken at. `lapse.reference_height_m`
# is linear in it and declines to pick one, so this is the liquid-water span the
# roughness builder already brackets over rather than a value chosen here.
REFERENCE_AIR_BRACKET_K = (273.15, 313.15)

# Shares for the two-class arm. The endpoints are the pure cases, which must
# give a gap of exactly one, and the middle is the mixture that maximises the
# spread of `ln z0` over a two-point distribution.
CLASS_SHARE_SWEEP = (0.0, 0.25, 0.5, 0.75, 1.0)

# The exactness bar for the drag-partition control. Float64 round-off on a
# dimensionless fraction, which is what an affine identity leaves behind.
AFFINE_TOLERANCE = 1.0e-12


def geometric_mean(values: np.ndarray, weights: np.ndarray) -> float:
    """The reduction `dust.yaml` performs: the weighted mean of `ln z0`."""
    w = np.asarray(weights, dtype=float)
    return float(np.exp((np.log(np.asarray(values, dtype=float)) * w).sum() / w.sum()))


def drag_efficiency(z0_m, dp: dict) -> np.ndarray:
    """MB95's drag partition. The same expression `build_dust.py` evaluates."""
    z0_cm = np.maximum(np.asarray(z0_m, dtype=float) * 100.0, dp["z0s_cm"] * 1.0001)
    denom = np.log(dp["a"] * (dp["x_cm"] / dp["z0s_cm"]) ** dp["b"])
    return np.clip(1.0 - np.log(z0_cm / dp["z0s_cm"]) / denom, 0.0, 1.0)


def mixture_arm(name: str, z0: np.ndarray, shares: np.ndarray,
                dp: dict, z_ref_m: tuple[float, float]) -> dict:
    """Both orders of one mixture, for the drag partition and for `u*^p`."""
    z0 = np.asarray(z0, dtype=float)
    w = np.asarray(shares, dtype=float)
    mixed = geometric_mean(z0, w)

    feff = drag_efficiency(z0, dp)
    clipped = bool(((feff <= 0.0) | (feff >= 1.0)).any())
    pta = float((w * feff).sum() / w.sum())
    ata = float(drag_efficiency(mixed, dp))
    control = {
        "process_then_aggregate": pta,
        "aggregate_then_process": ata,
        "absolute_difference": abs(pta - ata),
        "tolerance": AFFINE_TOLERANCE,
        "clip_binds_on_an_endmember": clipped,
        # The identity holds only where the clip does not bind, so a bound clip
        # is a FAIL of the arm's premise and not of its arithmetic.
        "passes": (not clipped) and abs(pta - ata) <= AFFINE_TOLERANCE,
    }

    friction = {}
    for z_ref in z_ref_m:
        length = np.log(z_ref / z0)
        length_mixed = np.log(z_ref / mixed)
        friction[f"{z_ref:.6g}"] = {
            f"u_star_power_{p}_ratio": float(
                (w * length ** -p).sum() / w.sum() / (length_mixed ** -p))
            for p in POWERS
        }
    return {
        "mixture": name,
        "z0_m": [float(v) for v in z0],
        "shares": [float(v) for v in w],
        "geometric_mean_z0_m": mixed,
        "ln_z0_spread": float(np.log(z0.max() / z0.min())),
        "drag_partition_control": control,
        "friction_velocity_by_reference_height_m": friction,
    }


def instrument(bracket: list[float], z_ref_m: tuple[float, float]) -> dict:
    """What the class's own declared level bracket is worth on `u*^p`.

    This is the ignorance the step declares before any question of operation
    order arises, and it is what a gap has to exceed to be material. It is read
    from `dust.yaml` rather than restated.
    """
    lo, hi = float(min(bracket)), float(max(bracket))
    out = {}
    for z_ref in z_ref_m:
        out[f"{z_ref:.6g}"] = {
            f"u_star_power_{p}_factor": float(
                (np.log(z_ref / lo) ** -p) / (np.log(z_ref / hi) ** -p))
            for p in POWERS
        }
    return {"z0_bracket_m": [lo, hi], "factor_by_reference_height_m": out}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()

    config = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    dust = yaml.safe_load(DUST_CONFIG.read_text(encoding="utf-8"))
    dp = dust["drag_partition"]
    table = dp["aeolian_z0_by_class_m"]
    mix = json.loads(MIX.read_text(encoding="utf-8"))

    z_ref = tuple(lapse.reference_height_m(t, config)
                  for t in REFERENCE_AIR_BRACKET_K)

    build = mix["configured_build"]
    bands = mix["endmembers"]
    band_z0 = np.array([b["z0_geometric_mean_m"] for b in bands])
    band_shares = np.array(mix["by_build"][build]["band_shares"])
    arms = [mixture_arm("playa_clastic landform bands", band_z0, band_shares,
                        dp, z_ref)]

    names = sorted(table)
    if len(names) == 2:
        a, b = (float(table[n]["z0"]) for n in names)
        for share in CLASS_SHARE_SWEEP:
            arms.append(mixture_arm(
                f"{names[0]} at {share:g} against {names[1]}",
                np.array([a, b]), np.array([share, 1.0 - share]), dp, z_ref))

    report = {
        "issue": "SPAT-7",
        "what_this_is": (
            "whether the order of the aeolian roughness mixture changes the "
            "drag partition or the friction velocity, and by how much against "
            "the level bracket each class already declares"),
        "reference_height_m": {
            "air_temperature_bracket_k": list(REFERENCE_AIR_BRACKET_K),
            "height_m": [float(v) for v in z_ref],
            "source": "lib/lapse.py:reference_height_m, this planet's gravity",
        },
        "powers": list(POWERS),
        "band_share_build": build,
        "arms": arms,
        "instrument_by_class": {
            name: instrument(table[name]["bracket"], z_ref) for name in names},
        "affine_tolerance": AFFINE_TOLERANCE,
        "provenance": provenance.config_stamp(
            config, "aeolian/scripts/roughness_mixing_order.py"),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    failed = 0
    for arm in arms:
        c = arm["drag_partition_control"]
        if c["clip_binds_on_an_endmember"] and abs(
                c["process_then_aggregate"] - c["aggregate_then_process"]) <= AFFINE_TOLERANCE:
            mark = " skip "        # the premise does not hold; not a result
        elif c["passes"]:
            mark = "  ok  "
        else:
            mark = " FAIL "
            failed += 1
        print(f"[{mark}] {arm['mixture']}: drag partition differs by "
              f"{c['absolute_difference']:.3g}, must be <= {AFFINE_TOLERANCE:.0e}")
    worst = max(
        (v[f"u_star_power_{max(POWERS)}_ratio"]
         for arm in arms
         for v in arm["friction_velocity_by_reference_height_m"].values()),
        default=1.0)
    print(f"\nworst u*^{max(POWERS)} ratio across every mixture and reference "
          f"height: {worst:.4f}")
    print(paths_lib.rel(args.out))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
