"""What CHOOSING a weathering scheme costs this world's thermostat.

`pedogenesis.yaml` computes weathering from one parameterisation and attaches no
structural uncertainty to it. Three published schemes sit on disk with their
constants, and they disagree about the FORM of the law rather than about a
constant inside one form. That disagreement is a thermostat strength, and this
measures it.

    python pedology/scripts/weathering_scheme_bracket.py

THE COMPARISON IS ANALYTIC AND NEEDS NO CLIMATOLOGY, which is deliberate. All
three schemes are multiplicatively separable into a power law in runoff and an
exponential in temperature, so

    thermostat strength  S = d ln W / dT = n * gamma + 1 / T_e

with `n` the runoff exponent, `T_e` the temperature e-folding, and `gamma` the
runoff sensitivity d ln(runoff) / dT. Only `n` and `T_e` belong to the scheme.
`gamma` belongs to the climate, this project cannot yet measure it, and it is
therefore swept rather than assumed: every number reported names its gamma.

WHAT THE MAPPING DOES. GKWM's exponent is per lithology, so the effective
exponent is the SILICATE-FLUX-weighted mean over this world's own rock classes,
which requires mapping seventeen land classes onto six. That mapping is declared
in `weathering_schemes.yaml` and is part of the bracket, not a detail under it.
`evaporite` maps to nothing on purpose and is excluded from the weighting rather
than being given a rock's exponent.

    Orogen's rock classes are the SIMULATED world's; every fraction here is a
    pre-carve LIMIT, because basin fill fires on the endorheic flag and opening a
    basin removes it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml

from _paths import ANALYSIS, CONFIG, PROJECT_ROOT, WEATHERING_SCHEMES

import builds  # noqa: E402  from lib/, via _paths.
from orogen import Export  # noqa: E402
from paths import rel  # noqa: E402

GAS_CONSTANT = 8.314472          # J / mol / K, the value rokgem itself uses
KELVIN = 273.15
UNDECLARED = "undeclared"
LAND = 1


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def e_folding_k(scheme: dict) -> float:
    """Kelvin per e-fold of silicate weathering.

    Declared directly where the scheme states one, and derived from the
    activation energy where it states that instead. The linearised Arrhenius
    d ln W / dT = E_a / (R T0^2) is the same functional form in both cases, which
    is what makes the two comparable at all.
    """
    if "temperature_e_folding_k" in scheme:
        return float(scheme["temperature_e_folding_k"])
    e_a = float(scheme["temperature_activation_energy_kj_mol"]) * 1000.0
    t0 = float(scheme["temperature_reference_c"]) + KELVIN
    return GAS_CONSTANT * t0 * t0 / e_a


def land_fractions(build: str | None) -> tuple[dict[str, float], str, str]:
    """Area share of each rock class over SUBAERIAL land, from the raw mesh.

    surface_class, never land_mask: the dry closed-basin floors below sea level
    carry the basin fill, so an elevation-sign denominator drops exactly the
    classes this bracket is most sensitive to.
    """
    root = builds.mesh_export_of(builds.build_root({"source_build": build}))
    export = Export(root)
    rock = export.substrate_class
    area = export.cell_area.astype(np.float64)
    land = export.surface_class == LAND
    total = area[land].sum()
    codes = {c["id"]: c["code"] for c in export.manifest["lithology"]["rockClasses"]}
    out = {}
    for rock_id, code in codes.items():
        selected = land & (rock == rock_id)
        if selected.any():
            out[code] = float(area[selected].sum() / total)
    return out, export.terrain_hash, str(rel(root))


def effective_exponent(scheme: dict, mapping: dict, fractions: dict) -> dict:
    """Silicate-flux-weighted mean runoff exponent over this world's lithology.

    Weighting by area alone would be wrong: what the thermostat responds to is
    the SILICATE flux, and a class contributes to it in proportion to its base
    rate times its silicate fraction as well as its area. Base rates are not
    comparable between schemes, but this weighting only ever compares within one.
    """
    if not scheme.get("runoff_exponent_is_per_lithology"):
        return {"effective_runoff_exponent": float(scheme["runoff_exponent"]),
                "weighting": "none; the exponent is uniform over lithology",
                "unmapped_land_fraction": 0.0}
    exps = scheme["runoff_exponent_by_class"]
    rates = scheme["base_rate_by_class"]
    fsi = scheme["silicate_fraction_by_class"]
    num = den = 0.0
    unmapped = 0.0
    per_target: dict[str, float] = {}
    for code, share in fractions.items():
        target = mapping.get(code)
        if target is None:
            raise SystemExit(
                f"rock class {code!r} has no entry in the class mapping. Add it "
                f"to {rel(WEATHERING_SCHEMES)} and record the reasoning; a class "
                f"silently dropped is a silently smaller bracket.")
        if target == UNDECLARED:
            unmapped += share
            continue
        per_target[target] = per_target.get(target, 0.0) + share
        w = share * float(rates[target]) * float(fsi[target])
        num += w * float(exps[target])
        den += w
    if den <= 0:
        raise SystemExit("no mapped land carries silicate weathering; the bracket is empty")
    return {"effective_runoff_exponent": num / den,
            "weighting": "area x base rate x silicate fraction",
            "mapped_land_fraction_by_target": per_target,
            "unmapped_land_fraction": unmapped}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-build")
    ap.add_argument("--declaration", type=Path, default=WEATHERING_SCHEMES)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    if not args.declaration.is_file():
        raise SystemExit(f"WEATHERING-SCHEMES-MISSING: {rel(args.declaration)} does "
                         f"not exist. There is no default set of schemes.")
    decl = yaml.safe_load(args.declaration.read_text(encoding="utf-8"))
    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    build = args.source_build or config.get("source_build")

    fractions, terrain_hash, export_root = land_fractions(build)
    gammas = [float(g) for g in decl["runoff_sensitivity_per_k_sweep"]]

    results = {}
    for key, scheme in decl["schemes"].items():
        # A scheme with a uniform exponent needs no mapping and gets none;
        # effective_exponent short-circuits on that and never reads one.
        mapping = decl["class_mapping"].get(key, {})
        if scheme.get("runoff_exponent_is_per_lithology") and not mapping:
            raise SystemExit(
                f"scheme {key!r} has a per-lithology exponent and no class "
                f"mapping. A bracket cannot be built without one.")
        eff = effective_exponent(scheme, mapping, fractions)
        te = e_folding_k(scheme)
        n = eff["effective_runoff_exponent"]
        results[key] = {
            "name": scheme["name"],
            "source": scheme["source"],
            "lithology_resolved": bool(scheme.get("lithology_resolved")),
            "temperature_e_folding_k": te,
            "temperature_term_per_k": 1.0 / te,
            **eff,
            "thermostat_strength_per_k": {f"{g:g}": n * g + 1.0 / te for g in gammas},
        }

    base = "whak_geocarb"
    spread = {}
    for g in gammas:
        vals = {k: v["thermostat_strength_per_k"][f"{g:g}"] for k, v in results.items()}
        lo, hi = min(vals.values()), max(vals.values())
        spread[f"{g:g}"] = {
            "min_per_k": lo, "max_per_k": hi, "ratio_max_over_min": hi / lo,
            "relative_to_current": {k: v / vals[base] for k, v in vals.items()},
        }

    # Which axis the spread is on. Held at the current scheme's other term so the
    # two are separated rather than confounded.
    n_lo = min(v["effective_runoff_exponent"] for v in results.values())
    n_hi = max(v["effective_runoff_exponent"] for v in results.values())
    t_lo = min(v["temperature_term_per_k"] for v in results.values())
    t_hi = max(v["temperature_term_per_k"] for v in results.values())
    axes = {f"{g:g}": {
        "runoff_axis_alone": (n_hi * g + t_lo) / (n_lo * g + t_lo),
        "temperature_axis_alone": (n_lo * g + t_hi) / (n_lo * g + t_lo),
    } for g in gammas}

    report = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "generator": "pedology/scripts/weathering_scheme_bracket.py",
        "declaration": str(rel(args.declaration)),
        "declaration_sha256": sha256(args.declaration),
        "source_build": build,
        "terrain_hash": terrain_hash,
        "mesh_export": export_root,
        "land_fractions_are_pre_carve_limits": True,
        "lithology_land_fractions": fractions,
        "runoff_sensitivity_per_k_sweep": gammas,
        "schemes": results,
        "thermostat_strength_spread": spread,
        "which_axis_carries_the_spread": axes,
        "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT,
                                     capture_output=True, text=True).stdout.strip() or None,
    }
    out = args.output or ANALYSIS / "weathering_scheme_bracket.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"build {build}   terrain {terrain_hash[:12]}")
    for key, v in results.items():
        print(f"  {key:14s} n_eff {v['effective_runoff_exponent']:.4f}   "
              f"T_e {v['temperature_e_folding_k']:.3f} K   "
              f"unmapped land {v['unmapped_land_fraction']:.4f}")
    print("\nthermostat strength d lnW/dT, per kelvin, by runoff sensitivity gamma")
    hdr = "  gamma  " + "".join(f"{k:>16s}" for k in results)
    print(hdr)
    for g in gammas:
        row = "".join(f"{results[k]['thermostat_strength_per_k'][f'{g:g}']:16.5f}" for k in results)
        print(f"  {g:<7g}" + row)
    print("\nspread, max over min, and which axis it is on")
    for g in gammas:
        s = spread[f"{g:g}"]; a = axes[f"{g:g}"]
        print(f"  gamma {g:<5g} total {s['ratio_max_over_min']:.3f}x   "
              f"runoff axis alone {a['runoff_axis_alone']:.3f}x   "
              f"temperature axis alone {a['temperature_axis_alone']:.3f}x")
    print(f"\nwrote {rel(out)}")


if __name__ == "__main__":
    main()
