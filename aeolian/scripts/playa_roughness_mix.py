#!/usr/bin/env python3
"""What aerodynamic roughness `playa_clastic` gets, from measured surfaces.

    python aeolian/scripts/playa_roughness_mix.py

Vesper is a simulated super-Earth. `playa_clastic` is one of this project's
derived lithology classes and everything here is a property of that class and of
the terrestrial field measurements it is being grounded on, not of any place on
this world.

WHY THIS EXISTS
---------------
`aeolian/config/dust.yaml` gave `playa_clastic` a roughness bracket of 1.0e-5 to
1.0e-4 m, and `notes/audits/dust-intensity-levers.md` section 4 measures that
this one class carries 10.415 of the 10.501 factor the whole roughness bracket
puts on dust emission. Neither end of that bracket was a measurement. The smooth
end was a modelling convention quoted by Prigent et al. (2005), "a low surface
roughness on the order of 0.001 cm is used to describe active dust sources"; the
rough end was that paper's sand-desert boundary, "z0 below ~0.02 cm", stepped
down one order on the argument that a clay plain is smoother than a sand sheet.

THE CLASS IS NOT A CLAY PLAIN, WHICH IS WHAT MAKES A NARROWING POSSIBLE.
`notes/audits/orogen-lithology.md` derives it end to end: it is the residual of
two geometric conditions, inside a preserved closed basin and not within the
lowest quarter of its relief, and Orogen's own name for it is "Playa mud /
alluvial fan fill". Measured on the mesh, half of it sits above 83 per cent of
its basin's relief, and on steepest descent it splits into three landform bands.
So the question the literature has to answer is not "what is a clastic playa's
z0" but "which measured surfaces is this class a mixture of, and in what
proportion". The mesh cannot supply a roughness -- roughness elements are
centimetric and the mesh edge is kilometric -- but it does supply the mixture.

THE RULES, FIXED BEFORE THE MIX WAS COMPUTED
--------------------------------------------
`notes/audits/dust-intensity-levers.md` carries the six criteria this was
declared under. The three that determine a number:

1. ENDMEMBER MEMBERSHIP is by the measuring paper's OWN printed surface type,
   and no site inside a type is dropped. The playa band takes every row printed
   "playa"; the sand-flat and distal-fan band takes every row printed "distal
   alluvial fan", "alluvial fan" or "interdune flats"; the fan and bajada band
   takes every row printed "proximal alluvial fan". Selecting sites inside a
   type is how a bracket gets chosen by the answer it is wanted to give.

2. THE COMBINATION IS THE GEOMETRIC MEAN, weighted by the band shares, which is
   the operation `dust.yaml` already uses between classes and for the same
   reason: the drag goes as 1/ln(z/z0), so it is ln(z0) that averages. Within a
   class and between classes are the same operation.

3. THE BRACKET IS THE MIX EVALUATED WITH EVERY ENDMEMBER AT ITS MEASUREMENT
   SET'S MINIMUM AND AT ITS MAXIMUM. That is the convention `evaporite` already
   carries in the same file, where MacKinnon et al. (2004)'s 0.0053 to 0.039 cm
   over one measurement set is taken as the bracket. The support dependence, the
   2.5M band shares against the 10M ones, is reported beside it and widens the
   bracket if it falls outside.

WHAT THE BRACKET IS AND IS NOT. It is a LEVEL uncertainty on the class's central
roughness. The WITHIN-class spatial spread is a separate quantity and
`dust.yaml` carries it separately as `drag_partition.within_class_z0`, so that
`dust_intensity_levers.py` can price the collapse without counting the same
spread twice.

THE CLASS EXTENT IS NOT USED ANYWHERE HERE, and that is deliberate. The band
shares are a statement about the class's COMPOSITION and survive the carve; the
class's area share is a pre-carve LIMIT and does not. Failure-modes class 8.

Writes `aeolian/analysis/playa_roughness_mix.json`.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from math import exp, log
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import PROJECT_ROOT  # noqa: E402

# --- the measurements --------------------------------------------------------
#
# Greeley, Blumberg, McHone, Dobrovolskis, Iversen, Lancaster, Rasmussen, Wall,
# White (1997). Applications of spaceborne radar laboratory data to the study of
# aeolian processes. JGR Planets 102(E5), 10971-10983. doi 10.1029/97JE00518.
# Table 2, "Mean z0 Measured, m". Every row is an AERODYNAMIC determination from
# a vertical anemometer array solved iteratively from the log law, which is the
# same quantity `dust.yaml` takes for `evaporite` from MacKinnon et al. (2004),
# so the two classes are grounded on one kind of measurement rather than two.
#
# TAKEN FROM GREELEY'S OWN TABLE AND NOT FROM PRIGENT ET AL. (2005) TABLE 1,
# which reproduces it. The reproduction has a transcription error: Greeley's
# Golden Canyon NE (April 1989) is 0.001100 m, which is 1.10e-1 cm, and Prigent
# prints 1.00e-0 cm, a factor of about 9 high. Failure-modes class 9 is exactly
# this, and this project has already been caught by it once.
#
# `surface` is the type as PRINTED in Greeley Table 2 and section 4. Membership
# is decided by that string and by nothing else.
GREELEY_1997_TABLE_2 = [
    # Death Valley, California
    dict(site="Stovepipe Wells Main (S), April 1994", z0_m=0.00347, surface="distal alluvial fan"),
    dict(site="Stovepipe Wells N1", z0_m=0.002080, surface="distal alluvial fan"),
    dict(site="Stovepipe Wells S1", z0_m=0.002880, surface="distal alluvial fan"),
    dict(site="Stovepipe Wells S2", z0_m=0.003100, surface="distal alluvial fan"),
    dict(site="Stovepipe Wells Main, March-April 1991", z0_m=0.00076, surface="distal alluvial fan"),
    dict(site="Kit Fox Fan T2", z0_m=0.000850, surface="alluvial fan"),
    dict(site="Kit Fox Fan T1", z0_m=0.000696, surface="alluvial fan"),
    dict(site="Golden Canyon NE, April 1994", z0_m=0.0106, surface="proximal alluvial fan"),
    dict(site="Golden Canyon SE, April 1989", z0_m=0.002450, surface="proximal alluvial fan"),
    dict(site="Golden Canyon NE, April 1989", z0_m=0.001100, surface="proximal alluvial fan"),
    # Lunar Lake, Nye County, Nevada. "The playa covers 14 km2 and has a
    # silty-clay surface." Dry and fractured when measured, desiccation cracks
    # 15 mm deep at 20 cm spacing near the centre passing to 3 mm at 2-9 cm
    # spacing southward. THE ONLY AERODYNAMICALLY MEASURED CLASTIC PLAYA IN THE
    # SOURCES THIS PROJECT HOLDS.
    dict(site="Lunar Lake LL-1, near centre, Aug 1994", z0_m=0.000182, surface="playa"),
    dict(site="Lunar Lake LL-2, near SW edge, Jul-Aug 1990", z0_m=0.000126, surface="playa"),
    # Gobabeb, Namibia. "The northern station is on a gravel lag, and the
    # southern station is on a thin sand sheet."
    dict(site="Gobabeb NB, gravel lag", z0_m=0.000420, surface="interdune flats"),
    dict(site="Gobabeb NA, thin sand sheet", z0_m=0.000040, surface="interdune flats"),
]

# The steepest-descent bands, and which printed surface types stand for each.
# The band definitions and the shares are `notes/audits/orogen-lithology.md`,
# section "The class spans a basin from its sump to its rim", measured on the
# raw mesh of both builds and area-weighted over `surface_class`.
BANDS = [
    dict(name="playa_and_mud_flat", gradient="below 1e-3",
         surfaces=("playa",)),
    dict(name="sand_flat_and_distal_fan", gradient="1e-3 to 1e-2",
         surfaces=("distal alluvial fan", "alluvial fan", "interdune flats")),
    dict(name="alluvial_fan_and_bajada", gradient="above 1e-2",
         surfaces=("proximal alluvial fan",)),
]

# Area-weighted share of `playa_clastic` in each band, per build. The CONFIGURED
# build is `precarve-craton-10m` and its column is the central weighting; the
# 2.5M column is the support dependence. Both are COMPOSITION and survive the
# carve, unlike the class's area share, which is not used here.
BAND_SHARES = {
    "precarve-craton": [0.6379, 0.3058, 0.0563],
    "precarve-craton-10m": [0.5328, 0.3731, 0.0940],
}
CONFIGURED_BUILD = "precarve-craton-10m"

# --- the within-class spread, from the same table -----------------------------
#
# `dust.yaml` also carries a geometric spread of z0 WITHIN one arid surface
# class, derived from the three groups of repeat in-situ values Prigent et al.
# (2005) Table 1 reproduces from this same Greeley table. It is derived HERE
# rather than beside it because it rests on the same rows and therefore on the
# same transcription: taking it from Prigent carries the factor-9 error on
# Golden Canyon NE 1989 into the Death Valley group, which is the group with
# nine of the eleven degrees of freedom.
#
# The estimator is the one `dust.yaml` states: pool the within-group variance of
# ln z0 over all three groups for the central value, and take the smallest and
# largest single-group estimates as the bracket.
REPEAT_GROUPS = {
    "death_valley": ["Stovepipe Wells Main (S), April 1994", "Stovepipe Wells N1",
                     "Stovepipe Wells S1", "Stovepipe Wells S2",
                     "Stovepipe Wells Main, March-April 1991",
                     "Kit Fox Fan T2", "Kit Fox Fan T1",
                     "Golden Canyon NE, April 1994", "Golden Canyon SE, April 1989",
                     "Golden Canyon NE, April 1989"],
    "namib_gobabeb": ["Gobabeb NB, gravel lag", "Gobabeb NA, thin sand sheet"],
    "central_nevada_lunar_lake": ["Lunar Lake LL-1, near centre, Aug 1994",
                                  "Lunar Lake LL-2, near SW edge, Jul-Aug 1990"],
}

# What is being replaced, so the report can say what moved and by how much.
DECLARED_BEFORE = dict(z0=3.0e-5, bracket=[1.0e-5, 1.0e-4],
                       sigma_g=2.736, sigma_g_bracket=[1.297, 5.095])

# The prediction fixed before any measured value was read, in
# `notes/audits/dust-intensity-levers.md`: the mix-derived HIGH end comes out
# above the declared 1.0e-4 m, because a gravel apron is rougher than a sand
# sheet rather than smoother. At or below it the prediction has failed.
PREDICTED_HIGH_END_EXCEEDS_M = 1.0e-4


def geometric_mean(values) -> float:
    return exp(sum(log(v) for v in values) / len(values))


def endmembers() -> list[dict]:
    """One entry per band: the measurements in it and their geometric mean."""
    out = []
    for band in BANDS:
        rows = [r for r in GREELEY_1997_TABLE_2 if r["surface"] in band["surfaces"]]
        if not rows:
            raise SystemExit(f"no measured surface stands for band {band['name']}")
        z0 = [r["z0_m"] for r in rows]
        out.append({
            "band": band["name"],
            "gradient": band["gradient"],
            "surface_types": list(band["surfaces"]),
            "sites": [{"site": r["site"], "surface": r["surface"], "z0_m": r["z0_m"]}
                      for r in rows],
            "n": len(z0),
            "z0_geometric_mean_m": geometric_mean(z0),
            "z0_min_m": min(z0),
            "z0_max_m": max(z0),
        })
    return out


def mix(shares, values) -> float:
    """Share-weighted GEOMETRIC mean, the operation dust.yaml uses on z0.

    The shares are transcribed from a table that prints them to four decimal
    places, so they sum to 1 only to that precision and are renormalised here.
    A sum further from 1 than the printing can explain is a transcription
    error rather than rounding, and stops the run.
    """
    total = sum(shares)
    if abs(total - 1.0) > 1.0e-3:
        raise SystemExit(f"band shares sum to {total}, which rounding at four "
                         "decimal places cannot explain")
    return exp(sum(f * log(v) for f, v in zip(shares, values)) / total)


def within_class_sigma_g() -> dict:
    """Geometric spread of z0 within one arid area, pooled over repeat groups."""
    by_site = {r["site"]: r["z0_m"] for r in GREELEY_1997_TABLE_2}
    groups, ss, df = {}, 0.0, 0
    for name, sites in REPEAT_GROUPS.items():
        ln = [log(by_site[s]) for s in sites]
        mu = sum(ln) / len(ln)
        s2 = sum((x - mu) ** 2 for x in ln) / (len(ln) - 1)
        groups[name] = {"n": len(ln), "sigma_g": exp(s2 ** 0.5)}
        ss += s2 * (len(ln) - 1)
        df += len(ln) - 1
    pooled = exp((ss / df) ** 0.5)
    singles = [g["sigma_g"] for g in groups.values()]
    return {"by_group": groups, "degrees_of_freedom": df,
            "sigma_g": pooled, "sigma_g_bracket": [min(singles), max(singles)]}


def main() -> None:
    members = endmembers()
    central_z0 = [m["z0_geometric_mean_m"] for m in members]
    low_z0 = [m["z0_min_m"] for m in members]
    high_z0 = [m["z0_max_m"] for m in members]

    per_build = {}
    for build, shares in BAND_SHARES.items():
        per_build[build] = {
            "band_shares": shares,
            "central_m": mix(shares, central_z0),
            "endmembers_at_their_minimum_m": mix(shares, low_z0),
            "endmembers_at_their_maximum_m": mix(shares, high_z0),
        }

    central = per_build[CONFIGURED_BUILD]["central_m"]
    # The bracket is the endmember-range mix on the configured build, widened by
    # the support dependence if the other build's central falls outside it.
    lo = min(per_build[CONFIGURED_BUILD]["endmembers_at_their_minimum_m"],
             *(b["central_m"] for b in per_build.values()))
    hi = max(per_build[CONFIGURED_BUILD]["endmembers_at_their_maximum_m"],
             *(b["central_m"] for b in per_build.values()))

    prediction_held = hi > PREDICTED_HIGH_END_EXCEEDS_M

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "what_this_is": (
            "the aerodynamic roughness of the derived lithology class "
            "playa_clastic, as a geometric mixture of measured arid surfaces "
            "weighted by the landform bands the class is made of"),
        "source": {
            "measurements": (
                "Greeley et al. (1997) Table 2, aerodynamic z0 from wind "
                "profiles, doi 10.1029/97JE00518; read from that table and not "
                "from Prigent et al. (2005) Table 1, which reproduces it with a "
                "factor-9 transcription error on Golden Canyon NE 1989"),
            "band_shares": (
                "notes/audits/orogen-lithology.md, steepest descent on "
                "elevation_km, area-weighted over playa_clastic, both builds"),
            "criteria": "notes/audits/dust-intensity-levers.md section 4",
        },
        "endmembers": members,
        "within_class_z0": within_class_sigma_g(),
        "by_build": per_build,
        "configured_build": CONFIGURED_BUILD,
        "result": {
            "z0_m": central,
            "bracket_m": [lo, hi],
            "bracket_factor": hi / lo,
        },
        "what_it_replaces": {
            "z0_m": DECLARED_BEFORE["z0"],
            "bracket_m": DECLARED_BEFORE["bracket"],
            "sigma_g": DECLARED_BEFORE["sigma_g"],
            "sigma_g_bracket": DECLARED_BEFORE["sigma_g_bracket"],
            "central_moves_by_factor": central / DECLARED_BEFORE["z0"],
        },
        "declared_prediction": {
            "statement": (
                "the mix-derived high end exceeds 1.0e-4 m, because a third to "
                "a half of the class lies on sand-flat and fan gradients and a "
                "gravel apron is rougher than a sand sheet rather than smoother"),
            "threshold_m": PREDICTED_HIGH_END_EXCEEDS_M,
            "high_end_m": hi,
            "held": bool(prediction_held),
        },
    }

    out = PROJECT_ROOT / "aeolian" / "analysis" / "playa_roughness_mix.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")

    print("endmembers, geometric mean of the measured sites in each band:")
    for m in members:
        print(f"  {m['band']:26s} n={m['n']:2d}  "
              f"{m['z0_geometric_mean_m']:.3e} m   "
              f"({m['z0_min_m']:.3e} to {m['z0_max_m']:.3e})")
    print("\nthe mix, per build:")
    for build, b in per_build.items():
        mark = "  <- configured" if build == CONFIGURED_BUILD else ""
        print(f"  {build:22s} shares {b['band_shares']}  "
              f"central {b['central_m']:.3e} m{mark}")
    print(f"\nplaya_clastic z0 {central:.3e} m, bracket "
          f"{lo:.3e} to {hi:.3e} m, a factor of {hi / lo:.1f}")
    print(f"  replaces {DECLARED_BEFORE['z0']:.1e} m over "
          f"{DECLARED_BEFORE['bracket'][0]:.1e} to "
          f"{DECLARED_BEFORE['bracket'][1]:.1e} m; the central moves by a "
          f"factor of {central / DECLARED_BEFORE['z0']:.1f}")
    print(f"\nthe prediction fixed before the measurements were read: the high "
          f"end exceeds {PREDICTED_HIGH_END_EXCEEDS_M:.1e} m")
    print(f"  high end {hi:.3e} m -- "
          f"{'HELD' if prediction_held else 'FAILED'}")
    wc = report["within_class_z0"]
    print(f"\nwithin-class geometric spread, pooled over the three repeat "
          f"groups on {wc['degrees_of_freedom']} degrees of freedom:")
    for name, g in wc["by_group"].items():
        print(f"  {name:26s} n={g['n']:2d}  sigma_g {g['sigma_g']:.3f}")
    print(f"  pooled sigma_g {wc['sigma_g']:.3f}, bracket "
          f"{wc['sigma_g_bracket'][0]:.3f} to {wc['sigma_g_bracket'][1]:.3f}")
    print(f"  replaces {DECLARED_BEFORE['sigma_g']}, bracket "
          f"{DECLARED_BEFORE['sigma_g_bracket'][0]} to "
          f"{DECLARED_BEFORE['sigma_g_bracket'][1]}, which was pooled over "
          f"Prigent's reproduction and carried its factor-9 error")

    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
