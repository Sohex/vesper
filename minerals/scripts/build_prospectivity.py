#!/usr/bin/env python3
"""Tectonic and magmatic ore prospectivity, as a 0-1 field per deposit type.

    python minerals/scripts/build_prospectivity.py

Writes `minerals/data/<source_build>/prospectivity.nc` on the native mesh, plus a
provenance report beside it.

## What this is, and what it deliberately is not

It is a FIELD, not deposits. A porphyry system is one to two kilometres against a
mesh cell several times that, and a vein is orders below one, so an individual
body is invisible at every resolution this pipeline runs at. The mesh spacing
follows the build's region count; `analysis/orogen_resolution.json` carries it
per build. Placing one would invent
detail the grid cannot hold. Discrete deposits belong with the downscaling pass,
which is also where glacial overdeepening goes, for the same reason and at almost
the same scale ratio. See `docs/src/reference/economic-minerals.md`.

It is NOT a lithology and NOT an erodibility modifier. `substrate_class` sets
erodibility and therefore terrain, albedo and therefore climate, texture and
therefore the biosphere, and solutes and therefore the weathering fluxes.
Anything entering it enters the climate path, and a mesh cell would be painted
with ore properties on the strength of a deposit occupying a fraction of a
percent of it.

## Why this is downstream rather than inside Orogen

`docs/src/reference/economic-minerals.md` assigns tectonic and magmatic genesis to Orogen,
meaning the PROCESSES that concentrate these deposits are the ones Orogen models.
It does not require the arithmetic to happen there, and every input but one is
already exported: craton weight, fold-belt weight, stress,
substrate and basement class, cover thickness and erosion delta. The missing one,
`lipV`, has its expression in the `flood_basalt` class.

Computing it here buys two things. Iterating on a rule costs no terrain rebuild,
so the layer stays freely rebuildable as the design says it should. And a value
that does not exist during generation cannot feed back into erosion or albedo, so
the constraint above becomes structural rather than a thing to remember.

## What the rules do

They are first order, they live in `config/prospectivity.yaml` rather than here,
and each names the control it encodes. The axis worth understanding is
exhumation, because it is what a rock map alone cannot tell you: a porphyry forms
1-5 km down and is destroyed by deep erosion, while an orogenic gold system forms
5-15 km down and is revealed by it. The same erosion that removes one exposes the
other.

## What 1.0 means

Each rule's own ATTAINABLE MAXIMUM, computed in closed form from the config by
`ceiling` below, and not this build's best cell. `prospectivity_scale.py` argues
why and `config/prospectivity.yaml`'s `input_ranges` block declares the tops of
the ranges it multiplies. The field therefore reaches 1.0 only where a cell
attains every modifier at once, and `max_over_land` in the report says how close
this build comes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import netCDF4 as nc
import numpy as np
import yaml

from _paths import CONFIG, DATA, PROJECT_ROOT, PROSPECTIVITY
from paths import rel  # noqa: E402
from prospectivity_scale import host_ceiling, host_weight, normalise

import builds
from orogen import LAND, Export


def ceiling(spec: dict, ranges: dict) -> float:
    """The largest score this rule can award, in closed form over the config.

    Mirrors the score expression in `main` term for term and in the same order,
    which is what makes it checkable: `normalise` raises if a score exceeds it.

    `require_cover` is absent because it multiplies by a boolean and so cannot
    raise the maximum, and the craton branch REPLACES the score rather than
    scaling it, which is why it returns rather than multiplying.
    """
    if spec.get("craton_weight"):
        return float(spec["craton_weight"]) * ranges["craton_weight"]
    top = host_ceiling(spec.get("hosts"), spec.get("name", ""))
    if spec.get("fold_belt_weight"):
        top *= 1.0 + spec["fold_belt_weight"] * ranges["fold_belt_weight"]
    if spec.get("stress_weight"):
        top *= 1.0 + spec["stress_weight"] * ranges["stress_norm"]
    if spec.get("favour_deep_exhumation"):
        top *= 1.0 + ranges["deep_exhumation"]
    return top


# `vendor/orogen/js/elevation.js` computes `r_t_foldBelt` as
# `min(1, stressNorm * FOLD_BELT_MULT)`, with FOLD_BELT_MULT declared in
# `terrain-config.js`. Restated here because the identity below IS that line,
# and a check that read the multiplier from the config it is checking would be
# checking nothing.
FOLD_BELT_MULT = 3.0


def check_frozen_references(rules: dict, mesh, land, stress_raw, fold,
                            erosion, substrate, basement, codes) -> None:
    """Hold the two DELIBERATELY frozen references to what they are frozen FROM.

    `stress_reference` and `deep_exhumation_delta` are quantiles of exported
    fields, frozen on purpose so that a score means the same thing on every
    build. That decision is sound and it is not self-enforcing: a frozen value
    still has to be told when the thing it was measured on has moved out from
    under it, and neither carried anything that could say so.

    Each is checked against what its own freeze is a decision ABOUT, and the two
    are different questions:

    `stress_reference` is Orogen's own divisor, so the check is the identity the
    config states. It does not ask whether the reference is still this build's
    0.97 quantile -- it is frozen so that it need not be -- it asks whether it is
    still the number Orogen normalised `r_t_foldBelt` by. A generation that
    changed the plate model or the propagation constants breaks that and nothing
    else does.

    `deep_exhumation_delta` has no identity, so the check is in QUANTILE space:
    the frozen threshold must still sit inside a declared band of the granite
    exhumation distribution. Comparing it against a re-taken p75 would refuse on
    any distribution shift, which is the thing a frozen threshold exists to
    absorb; leaving it unchecked lets "deep" quietly become "most of the map" or
    "almost nothing".
    """
    reference = float(rules["stress_reference"])
    tol = float(rules["stress_reference_identity_tolerance"])
    predicted = np.minimum(1.0, np.clip(stress_raw / reference, 0.0, 1.0)
                           * FOLD_BELT_MULT)
    resid = np.abs(predicted[land] - fold[land])
    worst = float(resid.max()) if resid.size else 0.0
    if worst > tol:
        raise SystemExit(
            f"stress_reference {reference} no longer reproduces "
            f"r_t_foldBelt on {rules['stress_reference_measured_on']!r}'s "
            f"successor: the identity "
            f"{rules['stress_reference_identity']} is out by {worst:.3e} at "
            f"worst over land, against a declared tolerance of {tol:g}.\n"
            "The reference is frozen deliberately, so this is not a licence to "
            "recompute it per build. It says the plate model or the stress "
            "propagation constants moved: re-take the 0.97 quantile of "
            "propagated stress over cells above STRESS_PROPAGATE_MIN, the way "
            "elevation.js does, and re-declare it with the build it came from.")

    ex = rules["exhumation"]
    delta = float(ex["deep_exhumation_delta"])
    lo, hi = (float(v) for v in ex["deep_exhumation_delta_quantile_band"])
    granite = codes.get("granite")
    if granite is None:
        raise SystemExit(
            "this export carries no `granite` rock class, so "
            "deep_exhumation_delta cannot be checked against the population "
            f"{ex['deep_exhumation_delta_population']!r} its declaration "
            "names. A build without granite is not one this rule was derived "
            "on.")
    pop = land & ((substrate == granite) | (basement == granite))
    if not pop.any():
        raise SystemExit(
            "no land cell on this build has granite as substrate or basement, "
            "so the population deep_exhumation_delta was measured over is "
            "empty here.")
    quantile = float((erosion[pop] <= delta).mean())
    if not lo <= quantile <= hi:
        raise SystemExit(
            f"deep_exhumation_delta {delta} sits at quantile {quantile:.3f} of "
            f"granite exhumation on this build, outside the declared band "
            f"[{lo}, {hi}]. It was taken as the p75 on "
            f"{ex['deep_exhumation_delta_measured_on']!r}.\n"
            "It is frozen on purpose, so the repair is not to re-take it "
            "silently: at this quantile the threshold no longer selects deep "
            "exhumation on the terrain being scored, and whether to re-take it "
            "or to keep it and accept what it now selects is a decision, taken "
            "before the scores it produces are read.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-build", default=None,
                    help="override config's source_build")
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()

    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    if args.source_build:
        config["source_build"] = args.source_build
    rules = yaml.safe_load(PROSPECTIVITY.read_text(encoding="utf-8"))

    mesh = Export(builds.mesh_export(config))
    land = mesh.field("surface_class") == LAND
    substrate = mesh.field("substrate_class").astype(int)
    basement = mesh.field("basement_rock").astype(int)
    thickness = mesh.field("cover_thickness").astype(float)
    erosion = mesh.field("erosionDelta").astype(float)
    craton = mesh.field("r_t_craton").astype(float)
    fold = mesh.field("r_t_foldBelt").astype(float)
    # Stress against the DECLARED reference, not against this build's own
    # maximum. The maximum is a single cell and it is a seafloor one, because
    # ocean stress runs higher than land stress and `r_stress` spans both; the
    # config says where the reference comes from and how to re-derive it.
    stress_raw = mesh.field("r_stress").astype(float)
    stress = np.clip(stress_raw / rules["stress_reference"], 0.0, 1.0)
    area = mesh.field("cell_area").astype(float)

    ranges = rules["input_ranges"]
    codes = {c["code"]: c["id"] for c in mesh.manifest["lithology"]["rockClasses"]}
    ex = rules["exhumation"]

    check_frozen_references(rules, mesh, land, stress_raw, fold,
                            erosion, substrate, basement, codes)

    cover_ok = thickness >= ex["cover_preserved_km"]
    deep = erosion >= ex["deep_exhumation_delta"]

    fields, summary = {}, {}
    for key, spec in rules["deposits"].items():
        # Host rock. Substrate is what is at the top of the stack; basement is
        # what a deposit sat in. Both count, because a deposit hosted in the
        # basement is still there when cover survives above it.
        score = host_weight(spec, substrate, basement, codes)

        if spec.get("fold_belt_weight"):
            score = score * (1.0 + spec["fold_belt_weight"] * fold)
        if spec.get("stress_weight"):
            score = score * (1.0 + spec["stress_weight"] * stress)

        # Exhumation, the axis that separates the hydrothermal families.
        if spec.get("require_cover"):
            score = score * cover_ok
        if spec.get("favour_deep_exhumation"):
            score = score * (1.0 + deep)

        if spec.get("craton_weight"):
            gate = craton >= spec.get("craton_threshold", 0.0)
            score = spec["craton_weight"] * craton * gate

        score[~land] = 0.0
        top = ceiling(spec, ranges)
        fields[key] = normalise(score, land, top, key)
        lit = land & (fields[key] > 0)
        summary[key] = {
            "name": spec["name"],
            # The divisor, so a reader can recover the raw score and see which
            # modifiers a rule had available to it.
            "attainable_maximum": top,
            "land_fraction_nonzero": float(area[lit].sum() / area[land].sum()),
            "land_fraction_above_half": float(
                area[land & (fields[key] > 0.5)].sum() / area[land].sum()),
            # How close the best cell on this build comes to what the rule can
            # award. Below 1 means no cell attains every modifier at once, which
            # the old land-maximum divisor reported as 1.0 whatever was true.
            "max_over_land": float(fields[key][land].max()),
            "mean_over_land": float(np.average(fields[key][land],
                                               weights=area[land])),
        }

    out = args.output or (DATA / str(config["source_build"]) / "prospectivity.nc")
    out.parent.mkdir(parents=True, exist_ok=True)
    with nc.Dataset(out, "w", format="NETCDF4") as data:
        data.createDimension("region", land.size)
        for key, field in fields.items():
            var = data.createVariable(key, "f4", ("region",), zlib=True)
            var[:] = field
            var.long_name = rules["deposits"][key]["name"]
            var.units = "1"
            var.description = ("relative prospectivity, 0-1 over land; NOT a "
                               "deposit and NOT a lithology")
            var.attainable_maximum = summary[key]["attainable_maximum"]
            var.scaling = ("fraction of this rule's attainable maximum, the "
                           "host weight times every modifier at the top of its "
                           "declared input range. A property of the rule, so "
                           "the same number means the same thing on any build; "
                           "the field reaches 1.0 only where some cell attains "
                           "every modifier at once")
        data.setncattr("vesper_source_build", str(config["source_build"]))
        data.setncattr("vesper_terrain_hash", mesh.terrain_hash)
        data.setncattr("note", "Tectonic and magmatic ore prospectivity. Read "
                               "only downstream of climate; never an input to "
                               "erodibility, albedo or any climate path.")

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_build": config.get("source_build"),
        "terrain_hash": mesh.terrain_hash,
        "config_sha256": hashlib.sha256(PROSPECTIVITY.read_bytes()).hexdigest(),
        "generator": "minerals/scripts/build_prospectivity.py",
        "git_commit": subprocess.run(["git", "rev-parse", "HEAD"],
                                     capture_output=True, text=True,
                                     cwd=PROJECT_ROOT).stdout.strip() or None,
        "deposits": summary,
        "note": ("Each deposit type is scaled by its own ATTAINABLE MAXIMUM, "
                 "computed in closed form from the config: the host weight "
                 "times every modifier at the top of its declared input range. "
                 "A 1.0 is a cell attaining everything the rule can award, not "
                 "an absolute grade or tonnage and not merely the best cell on "
                 "this build, so the same number means the same thing on any "
                 "build. Cross-type comparison of the numbers is still "
                 "meaningless: the maxima are different quantities."),
        "input_ranges": rules["input_ranges"],
        "stress_reference": rules["stress_reference"],
    }
    (out.parent / "prospectivity_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"{'deposit':28}{'% land':>9}{'% >0.5':>9}{'max':>8}{'mean':>8}"
          f"{'ceiling':>9}")
    for key, s in summary.items():
        print(f"  {s['name'][:26]:26}{100*s['land_fraction_nonzero']:9.2f}"
              f"{100*s['land_fraction_above_half']:9.2f}"
              f"{s['max_over_land']:8.3f}{s['mean_over_land']:8.3f}"
              f"{s['attainable_maximum']:9.3f}")
    print(f"\nwrote {rel(out)}")


if __name__ == "__main__":
    main()
