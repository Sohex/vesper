#!/usr/bin/env python3
"""Put every uncertainty in one currency, so they can be ranked against each other.

    python scripts/error_budget.py

The project accumulates uncertainties in whatever unit each was measured in --
albedo, W/m2, a percentage of land, a bracket on a rock class -- and then has no
way to answer "is refining this worth more than refining that". This converts
them at zero compute, from artifacts already on disk.

## The conversion

`lib/sensitivity.py` owns it. Nothing here declares a slope, a planetary albedo
or an insolation, because three different flux-to-kelvin conversions were in
simultaneous use across this repository and the budget's ranking depended on
which one the caller happened to pick. Read that module before changing anything
about the kelvin column.

A land-mean albedo change reaches temperature through three steps:

    surface, planet-wide   d_surface = d_land * land_fraction
    top of atmosphere      d_toa     = d_surface * attenuation
    temperature            d_T       = sensitivity.planetary_albedo_to_kelvin(d_toa)

**The attenuation is the step it is tempting to skip, and skipping it doubles
every answer.** A surface albedo change does not arrive at the top of the modelled
atmosphere intact; everything above the surface scatters and absorbs. It was a
declared 0.5 with a declared factor of two around it; it is now MEASURED on two
terms across three builds, and `ATTENUATION_PAIRS` below names the arms.
`notes/audits/albedo-attenuation.md` carries the measurement with its evidence,
its null and what it does not cover.

**The denominator is the STAGED land-mean albedo and never the diagnosed one.**
Every item in `albedo_items` is a delta of a boundary condition read from
`albedo_report.json`. A run's diagnosed `alb` is that boundary condition as
realised through the modelled snow; it is a different quantity, nothing can set
it, and back-solving against it instead moves the answer by a factor of 1.45 and
would refuse a single constant outright. Both columns are still printed, because
a bracket of 1.62 is still a bracket.

## The second currency, and why kelvin is not enough

The temperature is a result this project can revise. The carve list is not: it
leaves the project, changes the terrain, and cannot be undone. It is decided by a
water balance, so an item priced only in kelvin cannot be ranked against it.

Runoff over land is the residual of two numbers five to six times its own size,
so it amplifies: a fractional change in precipitation moves runoff by `P/R` times
as much, and one in land evaporation by `E/R`. Both are measured here from the
baseline climatology rather than declared.

`basin_response` then measures how many basins change verdict under a uniform
fractional perturbation of precipitation, land evaporation and lake evaporation.
It READS the three water terms the carve list records against each basin and
re-evaluates the criterion on them. It used to invert the recorded aridity index
and evaporation margin to recover them instead, which put the verdict's formula
in a second place and broke the first time the criterion moved; BUDG-5.

**The carve half has two unavailabilities and they are reported separately.**
Its outputs read different things, so each is reported exactly when the thing it
reads is present. The runoff percent needs the configured build's water balance
and `HYDROLOGICAL_RESPONSE_PER_KELVIN`. The basin count needs a carve list on
the configured build on top of both. The `CARVE_ITEMS` table and the basin
verdict table are denominated in the criterion's own channels, so no kelvin
reaches them and they need the carve list alone. The albedo and forcing halves
read no hydrography artifact whatsoever and neither unavailability touches them.
A build whose carve list has not been regenerated -- which is the ordinary state
of the build a carve pass PRODUCES -- therefore gets everything but the basin
column, and the reason where the basin tables would be. It used to be a
traceback that denied a reader every column that was fine.

**The kelvin reaches the runoff and basin columns through
`HYDROLOGICAL_RESPONSE_PER_KELVIN`**, which is the one input and is declared
once, at the top of this file. Two of its three channels are measured and the
third is bracketed, so the basin column is a RANGE and the runoff column is not:
runoff is `P - E_land` and does not contain lake evaporation at all.

**The trap this conversion exists to defuse.** Runoff amplifies a precipitation
change by `P/R`, which is over six on this world, and a naive reading multiplies
that by the precipitation response and calls the answer the temperature
sensitivity. It is not. A kelvin moves land evaporation slightly FASTER than
precipitation, so the two amplified terms nearly cancel and runoff moves by
about a fifth of a percent of what the precipitation-only amplification
suggests. The P/R amplification below is right for a precipitation-only
perturbation -- a change in the hydrological cycle that leaves the temperature
alone, which is what a dust or a circulation item can be -- and an order of
magnitude wrong for a temperature one. Every item in this budget is priced in
kelvin, so every item takes the second path.

## What it is for

Not precision. Ordering. A factor-of-two budget is enough to answer whether the
subgrid roughness work buys anything next to the slab depth, and the answer is
visible immediately because the surface items land under a kelvin and the
structural ones do not.

## The stopping rule that falls out

**Refine an input when its plausible range exceeds the effect of the thing you
last refined.** That is a per-component stopping rule rather than a global
optimisation: it says which slot to upgrade next, not that the wrong thing was
built. It also says when to stop, which this project has not previously had a way
to say.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
CONVERGENCE = ROOT / "exoplasim" / "analysis" / "convergence"
sys.path.insert(0, str(ROOT / "lib"))

import numpy as np  # noqa: E402
import yaml  # noqa: E402

import builds  # noqa: E402
import rungs  # noqa: E402
import orbit as orbit_lib  # noqa: E402
import climatology
import paths  # noqa: E402
import sensitivity  # noqa: E402

# The surface-to-planetary attenuation. MEASURED 2026-08-30, and the value is
# the MIDPOINT of the measured span rather than an average over the rows below:
# a mean would weight the term that happens to have three builds behind it.
# Median 0.411, mean 0.396 and 0.362 weighting the two terms equally all sit
# inside 13 per cent of it, against a bracket that spans 62 per cent, so the
# choice of rule is below what the measurement resolves.
#
# `notes/audits/albedo-attenuation.md` is the finding. What it replaced was a
# declared 0.5 carrying a declared factor of two, which is to say 0.25 to 1.00;
# the measured bracket is a factor of 1.62 and contains the value that was
# declared. Every albedo item therefore falls to 0.76 of the kelvin it used to
# report, and the ranking between them does not move, because the constant
# multiplies all of them.
#
# The span is the measured endpoints 0.2917 and 0.4725 rounded OUTWARD, so that
# it contains every row rather than clipping the one it was taken from.
# `verify_attenuation` checks exactly that, and it fired on a span rounded the
# other way.
DEFAULT_ATTENUATION = 0.38
ATTENUATION_SPREAD = (0.29, 0.48)

# WHERE IT CAME FROM, named so that something can check it, on the pattern
# `HYDROLOGICAL_RESPONSE_SECANT` below sets. Each row is a paired pair of
# equilibrated runs differing in one thing, with the temperature difference
# taken between the two convergence reports' FITTED ASYMPTOTES -- the estimator
# `lib/sensitivity.py` measures its own slope with.
#
# `d_land_staged` is NOT re-derivable by this file: it needs the staged surface
# fields the runs were handed, and PHYS-15's two arms have been deleted. So the
# albedo half of each row is declared and the KELVIN half is re-read, which is
# exactly the split the hydrological secant lives with and for the same reason.
# `verify_attenuation` re-reads the asymptotes, re-derives each attenuation from
# them, and refuses when a row no longer reproduces its declared value.
#
# THE BUILDS ARE NOT THE CONFIGURED ONE and that is deliberate rather than
# overlooked. There is no bare-rock endmember run on `canonical-10m-carve2`, so
# no pair on it exists to prefer. The cross-build read is defensible here in a
# way it is not for the hydrological response: that response composes with an
# amplification taken from the configured build's own water balance, so two
# builds are two worlds, while the attenuation is a property of the modelled
# atmosphere and the builds it is measured across differ only in carved
# closed-basin floors. The agreement BETWEEN the two builds is what makes that
# evidence instead of an assumption.
ATTENUATION_PAIRS = {
    "soil wetting": {
        "cold": "run_598eb57c5a34", "warm": "run_a1c35075747c",
        "source_build": "canonical-10m-base",
        "d_land_staged": -0.006745, "delta_t_k": 0.1958, "attenuation": 0.292,
        "note": "PHYS-15's NWETSOIL pair. Both run directories are deleted; "
                "their convergence diagnostics survive and are what is re-read. "
                "The staged delta is the staged pair mixed at the donor's own "
                "measured skin fill fraction, which is the figure the arm's "
                "prediction was registered on and which the wet arm's own "
                "restarts reproduced to half a percent.",
    },
    "biosphere endmember, uniform soil water": {
        "cold": "run_4b7281ee038a", "warm": "run_76e441e0a761",
        "source_build": "canonical-10m-carve1",
        "d_land_staged": -0.082856, "delta_t_k": 3.1381, "attenuation": 0.381,
        "note": "The CLEAN pair: neither arm stages a soil-water field capacity, "
                "so the arms differ in the six albedo codes and in the forest "
                "fraction, and NVEG = 0 makes the forest fraction an albedo "
                "effect alone through landmod's snowcanopymask.",
    },
    "biosphere endmember, pedology soil water": {
        "cold": "run_bd7de9ba67a4", "warm": "run_8ff97d5e189a",
        "source_build": "canonical-10m-carve1",
        "d_land_staged": -0.080385, "delta_t_k": 3.5256, "attenuation": 0.441,
        "note": "Carries a soil-water field capacity difference between its "
                "arms, because the pedology field depends on the vegetation. "
                "That is a non-albedo forcing and it is why the clean pair "
                "above sits lower.",
    },
    "biosphere endmember, base": {
        "cold": "run_58f467b0872d", "warm": "run_893e276ee029",
        "source_build": "canonical-10m-base",
        "d_land_staged": -0.069569, "delta_t_k": 3.2714, "attenuation": 0.473,
        "note": "Same soil-water field capacity difference as the row above.",
    },
}
# The land area fraction the staged deltas were meaned over: the T21 land-sea
# mask every one of these runs shares. It is NOT `land_fraction()` below, which
# is the mesh's surface_class figure the budget composes with, and the two
# differ by 1.2 per cent. Declared here so the back-solve is reproducible in the
# units it was taken in rather than in the units it is consumed in.
ATTENUATION_LAND_FRACTION = 0.427692
# The asymptotes are declared to four decimals, so a unit in the last place is
# what agreement can mean. Fixed on that precision, not on today's residual.
ATTENUATION_TOLERANCE = 0.005

# BUDG-1. Fractional change in each land-mean quantity per kelvin of global-mean
# surface temperature. This is the whole of the kelvin-to-carve conversion; the
# derivation and its evidence are in `notes/audits/hydrological-sensitivity.md`.
#
# `precipitation` and `land_evaporation` are MEASURED, as a secant between the
# two converged fluxes that share a surface, over the last ten orbits of each.
# No run was bought for them. They are a secant across about seven kelvin rather
# than a local slope, and design intent forbids reusing a sensitivity measured
# in one regime in another, so the converged points loop A's flux re-bracket
# produces are what localises them. That is a step, not a task.
#
# They must be measured with the SAME land mean `land_water_balance` uses, since
# the two are multiplied together: Gaussian quadrature weights, land from
# `lsm > 0.5`, annualised to Earth years. A response measured under a different
# land mean does not compose with an amplification measured under this one, and
# that is what `verify_hydrological_response` below tests: it is a statement
# about which BUILD each side sits on, not a style note.
#
# `lake_evaporation` is a BRACKET, because the Penman open-water response has
# not been measured on this world and the project's convention is to bracket
# rather than guess. The ends are declared:
#
#   low   the land rate. A lake responding no faster than the moisture-limited
#         ground beside it is the weakest response that is physically arguable.
#   high  6.7 %/K, the convexity of saturation vapour pressure the project
#         already measures (closed PHYS-5). At fixed relative
#         humidity the deficit driving Penman's aerodynamic term grows at that
#         rate and the radiative term does not grow at all, so the Penman
#         response at fixed radiation and wind is strictly below it.
#
# Runoff does not contain lake evaporation, so the runoff column is a single
# number and only the basin column is a range. Closing the bracket needs Penman
# evaluated on two converged climatologies, which loop A produces anyway.
HYDROLOGICAL_RESPONSE_PER_KELVIN = {
    "precipitation": 0.0219,
    "land_evaporation": 0.0242,
    "lake_evaporation": (0.0242, 0.067),
}

# WHERE THE TWO MEASURED CHANNELS CAME FROM, named so that something can check
# them. They were three bare floats: `lib/sensitivity.py`'s slope sits beside
# them, names its runs and re-reads their asymptotes, and this -- the conversion
# that decides the CARVE LIST, the one output that leaves the project and cannot
# be revised by re-running anything -- named nothing at all.
#
# The temperatures below are the ten-orbit window means each run's own index
# entry carries, so the kelvin the secant is divided by is re-derived rather
# than trusted. The fractional responses are NOT re-derivable: they need
# per-orbit `pr` and `evap` from the raw output, and both runs have been
# archived, so the numerator of this secant has no surviving input anywhere in
# this tree. The recipe in `notes/audits/hydrological-sensitivity.md` section 6
# is complete and has nothing left to read.
#
# THAT IS THE HONEST STATE OF THIS QUANTITY. What re-measures it is named below
# and is not a task: loop A's flux re-bracket puts two converged points on
# whatever build is configured, and the same recipe runs on those. Until it
# does, `verify_hydrological_response` refuses the pairing and every column
# priced through it reports as unavailable rather than reporting a number whose
# provenance is a deleted directory.
HYDROLOGICAL_RESPONSE_SECANT = {
    "cold": {"run_id": "run_bfa3f5269660", "flux_ratio": 0.910,
             "window_mean_k": 281.9216},
    "warm": {"run_id": "run_524fbed77a9a", "flux_ratio": 0.945,
             "window_mean_k": 288.9964},
    "window_orbits": 10,
    "delta_t_k": 7.0748,
    "source_build": "precarve-craton",
    "geography": "3a17c498",
    "resolution": "T42",
    "measured": "2026-08-18",
    "land_mean": "Gaussian quadrature weights, land from lsm > 0.5, annualised "
                 "to Earth years; the same mean land_water_balance takes",
    "re_measured_by": "the secant across loop A's flux re-bracket on the "
                      "configured build, by the recipe in "
                      "notes/audits/hydrological-sensitivity.md section 6",
}
# The window means are declared to four decimals, so a unit in the last place is
# what agreement can mean. Fixed on that precision, not on today's residual.
SECANT_TEMPERATURE_TOLERANCE_K = 0.0001


def verify_hydrological_response(config) -> list[str]:
    """Re-derive what survives of the carve conversion. Empty means it holds.

    Three questions, in the order that decides whether the answer is usable.

    1. Do the two runs still exist, at the flux and geography recorded, with the
       window means the secant was divided by? Re-read rather than trusted.
    2. Is the secant's kelvin still the declared one? An extended, re-run or
       re-assessed run moves it and trips this.
    3. **Does the response sit on the build the amplification sits on?** This is
       the one that bites today. `land_water_balance` takes P, E and runoff from
       the CONFIGURED build's baseline climatology, and the response is
       multiplied through it. The rule this module's own docstring states is
       that a response measured under one land mean does not compose with an
       amplification measured under another, and two builds are two land means:
       one world's amplification scaled by another world's sensitivity, with
       nothing about the product looking wrong.
    """
    # `sensitivity` owns run lookup, and it keeps the LIVE index and the
    # ARCHIVED identities APART on purpose: merging them is what let a bracket
    # measured on two deleted runs go on reproducing its own constant. Both are
    # read here and the archived one is reported as a problem in its own right
    # rather than merged in silently, so every reason this secant is unusable
    # comes out at once instead of the first one. Reading
    # `exoplasim/runs/INDEX.json` directly here would be the same class of
    # defect this function exists to close, so it is not done.
    live, archived = sensitivity._live_entries(), sensitivity._archived_entries()
    problems, means = [], {}
    configured = config.get("source_build")
    for role in ("cold", "warm"):
        want = HYDROLOGICAL_RESPONSE_SECANT[role]
        entry = live.get(want["run_id"])
        if entry is None:
            entry = archived.get(want["run_id"])
            if entry is None:
                problems.append(f"the {role} secant run {want['run_id']} is in "
                                "no run index, live or archived")
                continue
            problems.append(
                f"{want['run_id']} survives only as archived identity, on build "
                f"{entry.get('source_build')!r}. Its raw output is gone, so the "
                "two fractional responses this secant carries cannot be "
                "recomputed from it. Re-measure by "
                + HYDROLOGICAL_RESPONSE_SECANT["re_measured_by"])
        physical = entry.get("physical", {})
        flux = physical.get("flux_ratio")
        if flux is None or abs(float(flux) - want["flux_ratio"]) > 1e-9:
            problems.append(f"{want['run_id']} is at flux {flux}, not the "
                            f"{want['flux_ratio']} the secant was taken across")
        geography = physical.get("geography")
        if geography != HYDROLOGICAL_RESPONSE_SECANT["geography"]:
            problems.append(
                f"{want['run_id']} is on geography {geography}, not the "
                f"{HYDROLOGICAL_RESPONSE_SECANT['geography']} the secant was "
                "measured on, so the two ends no longer share a surface")
        build = entry.get("source_build")
        if build != configured:
            problems.append(
                f"{want['run_id']} is on build {build!r} and the water balance "
                f"it multiplies is on {configured!r}. The response and the "
                "amplification are then two different worlds and the product is "
                "a sensitivity of neither. Re-measure by "
                + HYDROLOGICAL_RESPONSE_SECANT["re_measured_by"])
        recorded = (entry.get("convergence_metrics") or {}).get("temperature_mean_k")
        if recorded is None:
            problems.append(f"{want['run_id']} carries no window-mean temperature")
            continue
        if abs(recorded - want["window_mean_k"]) > SECANT_TEMPERATURE_TOLERANCE_K:
            problems.append(
                f"{want['run_id']} now reports a window mean of "
                f"{recorded:.4f} K against the {want['window_mean_k']:.4f} the "
                "secant was divided by")
        means[role] = recorded
    if len(means) == 2:
        delta = means["warm"] - means["cold"]
        declared = HYDROLOGICAL_RESPONSE_SECANT["delta_t_k"]
        if abs(delta - declared) > 2 * SECANT_TEMPERATURE_TOLERANCE_K:
            problems.append(
                f"the secant now spans {delta:.4f} K against the declared "
                f"{declared:.4f}, so every fractional response divided by it "
                "has moved")
    return problems

def _asymptote(run_id: str) -> float | None:
    """The fitted asymptote from a run's convergence report, or None.

    Two filenames, because an A/B arm's orbits are assessed as an EXPERIMENT and
    written to `_convergence_diagnostic.json` while a spin-up's verdict goes to
    `_convergence.json`. Two of the four pairs below are arms, so a resolver that
    knew only one name would find half the rows and report the other half as
    missing runs, which is a different failure from the one it would be
    describing.
    """
    for suffix in ("_convergence.json", "_convergence_diagnostic.json"):
        path = CONVERGENCE / f"{run_id}{suffix}"
        if path.is_file():
            metrics = json.loads(path.read_text(encoding="utf-8"))["metrics"]
            return float(metrics["temperature_asymptote_k"])
    return None


def verify_attenuation(planetary_albedo: float) -> list[str]:
    """Re-derive every attenuation row from the reports it names. Empty agrees.

    The albedo half of each row cannot be re-read here -- it needs staged fields,
    and PHYS-15's arms are deleted -- so this checks the half that CAN move. An
    extended, re-run or re-assessed arm moves its asymptote, which moves the
    separation, which moves the attenuation that separation implies. Three of the
    four rows would go on agreeing with a declaration nobody had re-derived
    otherwise, which is what `lib/sensitivity.py` learned the expensive way.

    It also checks the declared span still contains every row and the value, so
    the bracket cannot be left describing a set it no longer covers.
    """
    problems = []
    lo, hi = ATTENUATION_SPREAD
    if not lo <= DEFAULT_ATTENUATION <= hi:
        problems.append(f"DEFAULT_ATTENUATION {DEFAULT_ATTENUATION} is outside "
                        f"its own declared span {ATTENUATION_SPREAD}")
    for label, row in ATTENUATION_PAIRS.items():
        cold, warm = _asymptote(row["cold"]), _asymptote(row["warm"])
        missing = [r for r, a in ((row["cold"], cold), (row["warm"], warm))
                   if a is None]
        if missing:
            problems.append(f"{label}: no convergence report for "
                            f"{', '.join(missing)}, so the separation the "
                            "attenuation was back-solved from cannot be re-read")
            continue
        delta = warm - cold
        if abs(delta - row["delta_t_k"]) > ATTENUATION_TOLERANCE:
            problems.append(
                f"{label} now separates by {delta:.4f} K against the declared "
                f"{row['delta_t_k']:.4f}, so its attenuation has moved")
        measured = albedo_to_attenuation(row["d_land_staged"], delta,
                                         planetary_albedo)
        if abs(measured - row["attenuation"]) > ATTENUATION_TOLERANCE:
            problems.append(
                f"{label} re-derives to {measured:.4f} against the declared "
                f"{row['attenuation']:.4f}")
        if not lo <= row["attenuation"] <= hi:
            problems.append(
                f"{label}'s {row['attenuation']} is outside the declared span "
                f"{ATTENUATION_SPREAD}, which is meant to be that set's own")
    return problems


def albedo_to_attenuation(d_land, delta_t_k, planetary_albedo) -> float:
    """Invert `albedo_to_kelvin` for the attenuation. The one back-solve.

    Written once because the measurement and the check are the same arithmetic,
    and a second copy is how the two would drift apart. The land fraction is
    `ATTENUATION_LAND_FRACTION`, the mean the deltas were taken over, and not
    `land_fraction()`, the mesh figure the budget consumes: the back-solve has to
    run in the units it was measured in.
    """
    unit = albedo_to_kelvin(d_land, ATTENUATION_LAND_FRACTION, planetary_albedo,
                            1.0)
    return delta_t_k / unit


# Perturbation sizes the basin response is tabulated at. Both signs, because the
# criterion is a threshold and the basin population is not symmetric about it.
PERTURBATIONS = [-0.10, -0.05, -0.02, -0.01, 0.01, 0.02, 0.05, 0.10]


# Items already denominated in a fractional change of one of the criterion's own
# water terms. They reach basins EXACTLY, with no kelvin anywhere in the path,
# and no temperature can be claimed for them. Keeping them apart is the point:
# an item that perturbs the water cycle without perturbing the temperature takes
# the full `P/R` amplification, and an item priced in kelvin never does.
# (label, channel, (low, high) fractional change, note)
CARVE_ITEMS = [
    ("Penman over a dry land column", "lake_evaporation", (-0.18, -0.10),
     "Validated over ocean, where the air is equilibrated with the surface. Over "
     "a subgrid lake the column is dry, so the vapour pressure deficit is too "
     "high, E is overstated and the verdict under-carves. Quantified by "
     "notes/audits/missed-couplings.md finding 2 by raising the assumed column "
     "humidity: -10% at 70%, -13% at 80%, -18% at 90%. It is the largest "
     "carve-side item in this budget and it is opposite in sign to the "
     "albedo-driven over-carve. It does not touch runoff, because runoff is "
     "P - E_land."),
]

# The criterion channel each CARVE_ITEMS entry perturbs, as the argument
# position `basin_response`'s `overflowing` takes.
CARVE_CHANNELS = {"land_precipitation": (1, 0, 0),
                  "land_evaporation": (0, 1, 0),
                  "lake_evaporation": (0, 0, 1)}


def albedo_report_path(config) -> Path:
    """The albedo report for the CONFIGURED rung.

    It was a module-level constant naming `inputs/t42/`, which is a directory
    this tree no longer has: every surface input family is rooted by
    resolution. SPAT-2.
    """
    return (ROOT / "exoplasim" / "inputs"
            / rungs.model_grid(config)[0].lower() / "albedo_report.json")


# A ROCK CLASS'S SHARE OF THE STAGED LAND-MEAN ALBEDO: how far the land mean
# `build_surface_albedo.py` writes moves per unit of that class's albedo. It is
# NOT the class's area fraction of land, and the difference is what a literal
# 0.147 here was: an area-shaped number standing where a staged-mean sensitivity
# belongs, overstating the playa item by a factor of 1.727. That is the same
# error one level down from the one `notes/audits/albedo-attenuation.md` settles
# -- the denominator of an albedo item is the STAGED land mean, so its numerator
# has to be a derivative of that same staged mean.
#
# IT IS READ FROM THE REPORT, and that is the whole of WORLD-SU9O. It was a
# declared constant measured by differencing three generator runs, and a measured
# constant with no re-derivation is a number that drifts from the terrain it
# describes with nothing objecting: the same declaration read 0.0851423, taken on
# a staging that is not the one this budget reads, while the generator on the
# report it does read gives another figure. `build_surface_albedo.py` now carries
# a per-class indicator through the identical arithmetic the albedo takes and
# emits the exact derivative under each `lithology_albedo_overrides` entry, so
# the consumer reads the emitted value and there is no second copy to disagree
# with.
#
# A REPORT WITHOUT ONE IS REFUSED rather than fallen back on, for the reason the
# endmembers above are: a share from before the key existed is a share for
# another staging, and it would price this item as though it had been measured.


def land_mean_share(report: dict, rock_class: str) -> float:
    """The staged land mean's derivative with respect to one class's albedo.

    Read from the albedo report, never declared. `--mode scaled` clips the
    reduced field, which makes the response piecewise and a finite step able to
    cross a rail, so the generator writes a refusal there instead of a number and
    this passes that refusal on rather than pricing the item at a value that is
    exact only for an infinitesimal step.
    """
    entry = ((report.get("lithology_albedo_overrides") or {}).get(rock_class)
             or {})
    if "land_mean_share" not in entry:
        raise SystemExit(
            f"the staged albedo report records no land_mean_share for "
            f"{rock_class}: it predates the key. Rebuild it with "
            "`build_surface_albedo.py --lakes <surface_water.nc>` under "
            "model.land_albedo_source: vegetated. There is deliberately no "
            "declared value to fall back on -- one measured on another staging "
            "would price this item as though it had been measured on this one."
        )
    share = entry["land_mean_share"]
    if share is None:
        raise SystemExit(
            f"the staged albedo report refuses a land_mean_share for "
            f"{rock_class}: "
            f"{entry.get('land_mean_share_refused', 'no reason recorded')}"
        )
    return float(share)


def albedo_items(config) -> list[tuple[str, float, str]]:
    """(label, land-mean albedo delta, note) for the albedo half of the budget.

    Three of these are MEASURED by `build_surface_albedo.py` and are read from
    the report it writes rather than transcribed. They were literals until 2026-08-18
    and both had gone stale against the generator: the lakes item was -0.0139
    against a measured -0.0229, understated by 64%, and the biosphere item paired
    a bare-rock mean of 0.276 that no longer existed with a vegetated one that
    did. A number the pipeline computes has no business being typed in here.

    The biosphere pair is the subtle one and the report now settles it rather
    than leaving it to be inferred. Both endmembers are taken on the SAME rock
    table, overrides applied, and both PRE-LAKE. `land_mean_bare_rock` in the
    same report is a different quantity -- the raw export mean, still carrying
    the exported playa albedo the override replaces -- and pairing that with a
    vegetated mean compares two rock tables. Pairing a lake-composited vegetated
    mean with an uncomposited bare one double-counts the lakes, which are the
    line below.

    Deltas are the plausible RANGE of the item, not its current value: the budget
    ranks what refining each would buy.
    """
    ALBEDO_REPORT = albedo_report_path(config)
    if not ALBEDO_REPORT.exists():
        raise SystemExit(f"{ALBEDO_REPORT} is missing; run the surface_albedo step")
    report = json.loads(ALBEDO_REPORT.read_text(encoding="utf-8"))
    ends = report.get("endmembers") or {}
    missing = [k for k in ("bare_rock", "vegetated") if k not in ends]
    if missing:
        # No fallback. A guessed comparand here is a wrong number that reads as
        # a measured one, which is the whole reason this stopped being a literal.
        raise SystemExit(
            f"{ALBEDO_REPORT} carries no endmembers{tuple(missing)}: it predates "
            "the key, or was written in a mode that does not define them. "
            "Rebuild it with `build_surface_albedo.py --lakes <surface_water.nc>` "
            "under model.land_albedo_source: vegetated."
        )
    if report.get("lakes") is None:
        # The KEY is written whether or not the block is, so a membership test
        # passed on a report built without --lakes and the item died on a None
        # a few lines later. A refusal that names the rebuild is the point of
        # having one at all.
        raise SystemExit(
            f"{ALBEDO_REPORT} has no lakes block, so it was built without "
            "--lakes and the lake item cannot be read from it."
        )
    # The wetting item is CONDITIONAL on the saturated pair having been staged.
    # `--flat-bands` writes no pair for the mixing to run between, and a report
    # from before PHYS-15 armed carries no `at_full_surface_layer`; in both
    # cases the term is genuinely not in the field this budget prices, so it is
    # omitted rather than entered at zero. A zero row would read as "measured
    # and negligible", which is a different claim.
    wetting = []
    full = (report.get("wetting") or {}).get("at_full_surface_layer")
    if full:
        wetting = [
            ("soil wetting, dry to a full surface layer",
             float(full["broadband_fall"]),
             "PHYS-15, and a BOUND rather than a range: the fall at a surface "
             "layer permanently at its own capacity, which the tipping-bucket "
             "cascade cannot hold. Read from the same report; what refining it "
             "buys is the measurement of how wet the modelled skin actually "
             "is, since the term is this bound times that fraction.")]
    return wetting + [
        ("biosphere: bare rock vs vegetated",
         -(float(ends["bare_rock"]) - float(ends["vegetated"])),
         "the assumption LPJ-GUESS exists to replace. Measured: "
         f"{ends['bare_rock']:.5f} bare against {ends['vegetated']:.5f} "
         "vegetated, both pre-lake and on the overridden rock table."),
        ("playa_clastic albedo, 0.25 to 0.33",
         0.08 * land_mean_share(report, "playa_clastic"),
         "the largest single rock-class lever: a mixture of clay playa and "
         "varnished fan. The 0.08 step is the declared range of the class "
         "albedo and the share it multiplies is READ FROM THE REPORT, where "
         "the generator emits it as the exact derivative of the land mean it "
         "just wrote. It was 0.147, which is an area-shaped number where a "
         "staged-mean sensitivity belongs and overstated this item by 1.727. "
         "The staged value is 0.23, which sits BELOW this range, so the range "
         "is the class's plausible spread and not a bracket around what is "
         "staged."),
        ("lakes composited into albedo", float(report["lakes"]["delta"]),
         "solved lakes reaching the climate at all. Measured, from the same "
         "report: the cells carrying water are the bright playa and salt crust."),
        ("carve iteration, fill 16.5% to 12.35%", 0.0077,
         "one carve pass, measured across v4 to v5"),
        ("salt crust albedo, 0.40 to 0.50", 0.10 * 0.019,
         "degenerate with crust extent from the climate's side"),
    ]

# Items already in global-mean top-of-atmosphere W/m2. These convert exactly,
# because that is the quantity `lib/sensitivity.py` is denominated in.
# (label, (low, high) W/m2, note)
FORCING_ITEMS = [
    ("dust, radiative", (0.34, 0.61),
     "PRICED 2026-08-17, DUST-2, at the burden DUST-1 produces: land-mean "
     "optical depth 0.376. See analysis/dust_forcing.json and "
     "exoplasim/scripts/dust_forcing.py. Shortwave -8.7 W/m2 over ocean, -5.6 "
     "over vegetated land, -0.1 over playa and +2.0 over salt crust; longwave "
     "+6.7 W/m2 everywhere as a clear-sky-window UPPER bound. The sign is "
     "surface-dependent and POSITIVE over the bright closed-basin fill that "
     "makes this world unusual. The global mean is below the 1.5 W/m2 reopening "
     "threshold only because ocean cooling cancels land warming: the net spans "
     "11 W/m2 by surface, and that spread is what drives circulation. The "
     "global-mean kelvin here is therefore the least informative thing about "
     "this item."),
    ("sea salt, radiative", (-0.89, -0.16),
     "PRICED 2026-08-19, CLIM-27, from aeolian/analysis/sea_salt_baseline.json. "
     "The counterpart the budget had carried dust without: the same magnitude "
     "and the OPPOSITE sign. Sea salt's single-scattering albedo is 1 to within "
     "1e-5 in both bands, so the two-stream expression has no absorbing term "
     "and the sign cannot come out positive over any surface. The range spans "
     "two things that are not the same kind of uncertainty and are deliberately "
     "not separated here: whether Grythe's 30 um spume mode counts inside the "
     "10 um cut, worth a factor of 2.3, and the wet-removal lifetime bracket, "
     "worth a factor of 2.6. Central estimates are -0.30 excluding the spume "
     "mode and -0.68 including it. Longwave NOT computed: the layer sits in a "
     "boundary layer near the surface temperature, so the term is expected "
     "small, and expected is not measured. notes/audits/unpriced-terms.md "
     "finding 2. "
     "PAIRED WITH THE DUST ROW, AND NOT AN ALTERNATIVE TO IT. Both aerosols "
     "are off in every run on this build, so the two are simultaneous "
     "OMISSIONS rather than competing estimates of one thing, and the honest "
     "combination is their sum. Centrally that sum is about -0.07 W/m2, which "
     "is to say they very nearly cancel; the spread of the sum is -0.58 to "
     "+0.44. Do not read the cancellation as licence to drop either. It is a "
     "GLOBAL-MEAN cancellation between an aerosol that acts over land and one "
     "that acts over ocean, so it does not hold by surface, it does not hold "
     "in runoff -- where dust is several times sea salt because sea salt "
     "reaches land only through temperature -- and it is exactly the reading "
     "the dust row already warns against for dust alone."),
    ("volcanic sulfate, radiative", (-0.032, -0.013),
     "PRICED 2026-08-19, CLIM-28, from aeolian/analysis/volcanic_sulfate.json. "
     "Carried because it is MEASURED SMALL rather than assumed so, and because "
     "the same measurement bounds two sulfur sources this project cannot "
     "compute at all. The source is a coupling rather than a scaling: it is "
     "Carn et al.'s satellite-measured passive volcanic SO2 flux scaled by "
     "weathering_fluxes.py's own implied outgassing requirement, distributed "
     "over the arc classes. The sulfate PATHWAY is weak here -- a four-day "
     "aerosol from a source two orders of magnitude below the wind-driven "
     "ones -- so reaching sea salt's forcing would need a sulfur flux more "
     "than ten times Earth's, which is the bound that also covers explosive "
     "eruptions and the marine biogenic sulfur this project has no biosphere "
     "for. notes/audits/unpriced-terms.md finding 2."),
    ("CH4 and N2O, two estimates of one term that differ by 2.5x", (0.76, 1.18),
     "REVISED 2026-08-20, CLIM-42, when the band stopped being absent. "
     "`radmod.f90` now carries three bands on Donner and Ramanathan (1980), so "
     "this prices a DISAGREEMENT rather than an omission. The band lowers "
     "global mean outgoing longwave by 0.799 W/m2 at this world's 1.600 ppmv "
     "of CH4 and 0.300 of N2O against zero; the offline pricing for the same "
     "change is a little above 2. The range is that difference. "
     "NEITHER ESTIMATE IS REFUTED AND THEY ARE NOT THE SAME CALCULATION. Byrne "
     "and Goldblatt is line-by-line for an EARTH atmosphere -- its temperature "
     "profile, its water vapour, its clouds -- and 0.799 is this world's "
     "atmosphere in this scheme. Every step between them has been checked "
     "against a line list and none is where the difference lives: the CH4 band "
     "absorptance is within 7% of HITRAN 2020 at every model level, the "
     "narrow-band conversion to a fractional absorptivity reproduces the true "
     "Planck-weighted value to 1.004, and Sasamori's own CO2 absorptivity "
     "comes out 0.175 against a line-list 0.165. "
     "SO DO NOT READ THIS AS THE SCHEME BEING WEAK. An earlier revision of "
     "this row said the scheme runs at 40 to 45 percent of line-by-line, on an "
     "Earth CO2 greenhouse of 25 to 30 W/m2 that was quoted from general "
     "knowledge, never sourced, and is an ATTRIBUTION with overlap shared out "
     "rather than the REMOVAL measured here. Against a verified absorptivity "
     "the model's 12.28 W/m2 is what this atmosphere gives. "
     "NOT the thinner column either: 1 bar at 12.81 m/s2 is 0.756 of Earth's, "
     "and running at the 595 ppm that restores Earth's column moves the CO2 "
     "greenhouse by 0.642 W/m2. Gravity does that and radius does not -- the "
     "column over a square metre is P/g. "
     "THE MODEL HALF IS THE WEAK ONE: one 56-timestep segment on one bed, not "
     "a converged climate. Re-measure on the first baseline that carries the "
     "band. NO RUN ON DISK CARRIES IT, and that does not get a second row "
     "pricing what those runs lack -- a re-commission makes them worthless "
     "rather than stale, CLAUDE.md rule 7, and this budget prices the world "
     "going forward. "
     "THE SIGN IS KNOWN if the offline half is right: the band would be "
     "under-delivering, the modelled world short of greenhouse, and "
     "derive_design_flux.py still meeting its thresholds with too much "
     "starlight. The mixing ratios themselves are measured rather than "
     "assumed; CLIM-43 and analysis/trace_gas_forcing.json."),
    ("dust, shortwave-only if switched on as shipped", (-4.0, -3.4),
     "DUST-3's constraint, not an uncertainty in the world. Aerosols ARE "
     "switchable: aero_ini is called at plasim.f90:190 and reads aero_nl with "
     "l_aerorad and aerofile use-associated from radmod, and Model.configure "
     "writes l_aerorad, aerofile, L_AERO, l_source, apart, rhop and fcoeff. But "
     "the aerosol acts in the two SHORTWAVE bands only and the longwave solver "
     "has no aerosol term, so enabling it applies -3.4 to -4.0 W/m2 of "
     "global-mean cooling against a true +0.34 to +0.61. The fork needs a "
     "longwave aerosol term, not only an emission scheme."),
    ("stellar spectrum, k2 vs k25v", (0.0, 0.04),
     "Measured null on a warm nearly ice-free state, NOT on the cold branch, "
     "where it is the whole point. Re-measure before the cycle work."),
    ("energy closure residual", (-0.489, -0.400),
     "Consumes most of the 0.5 W/m2 convergence tolerance. Was recorded as a "
     "fixed offset; it is not. The kelvin is the equilibrium equivalent of the "
     "imbalance, which is warming still owed rather than a bias in the state."),
]

# Items in their own units. No conversion invented.
OTHER_ITEMS = [
    ("slab depth, 25 m and 100 m arms", "seasonal amplitude 1.51x and 0.66x",
     "STRUCTURAL and MEASURED, archive CLIM-33. Amplitude scales as 1/(omega*C) "
     "and this world's year is half Earth's, so 50 m damps seasonality about "
     "twice as hard as Earth's ocean does; coldest-month mean is what PFT "
     "survival gates on. The arms bracket the declared 50 m rather than "
     "searching for a better value. The annual global mean is 0.000 K BY "
     "CONSTRUCTION and measured so -- a heat capacity cannot move an "
     "equilibrium -- with the 100 m arm under a millikelvin. The one route to "
     "the mean is SEA ICE: the 25 m arm grew seasonal ice by 0.53% of the "
     "planet, worth -0.428 W/m2 against the 80 to 129 W/m2 per unit ice "
     "fraction that archive CLIM-17 derived for a different task and did not "
     "fit to this one. Amplitude reached 1.51x at 25 m and 0.66x at 100 m "
     "against a 2.00x/0.50x asymptote, compressed because two orbits from a "
     "50 m-equilibrated state is a transient. What remains for this term is "
     "the COLUMN rather than the slab: OCN-1."),
    ("no q-flux", "gradients too strong, ice too extensive",
     "STRUCTURAL in the sense that no ocean circulation is solved, but the "
     "claim that no cheap version exists was FALSE and is corrected here: "
     "oceanmod.f90 declares `nhdiff` and `hdiffk` in oceanmod_namelist, "
     "broadcasts both, and `mksst` acts on `nhdiff > 0`, so horizontal "
     "heat diffusion is one namelist key on the binary already built. That "
     "satisfies every one of `docs/src/pipeline/sequencing.md` A3's four A/B conditions with no further "
     "work, which makes this the CHEAPEST structural item to price rather than "
     "the impossible one. Read it as physics-is-not-a-knob: ocean heat "
     "transport exists, and the argument for switching it on is not that it "
     "improves an agreement. Bracket `hdiffk` and report the spread rather "
     "than tuning it -- a constant diffusivity is a BOUND on the missing "
     "transport, not the transport. CLIM-16. "
     "TWO THINGS ABOUT THIS TERM MOVED. `hdiffo` divided the requested "
     "coefficient by Earth's radius squared on an angular grid, so a bracket "
     "arm realised (a_earth/a)**2 times what it declared; it now runs on this "
     "planet's radius, and hdiffk means what it says. Any arm measured before "
     "that carries the wrong label on its coefficient, though not on its "
     "conservation. And the term is now STRUCTURALLY PERMANENT rather than "
     "merely unpriced: `nfluko` is ExoPlaSim's actual q-flux, both of its "
     "branches RELAX the modelled state toward a sea surface temperature and "
     "sea-ice climatology, and a world that has none has nothing to relax "
     "toward. icemod and oceanmod now refuse `nfluko` outright on that test "
     "rather than relaxing toward a field constructed from a sentinel. "
     "Diffusion is the bound this term gets; a q-flux is not available here "
     "at any price. world-6fh and world-mll."),
    ("ocean salinity, 20 to 40 psu", "freezing point 272.07 to 270.94 K",
     "DECLARED, not derived: config/planet.yaml carries Earth's 34.7 because no "
     "salt budget exists here, and the model's TFREEZE, CRHOS and CPS are now "
     "set from it rather than compiled in -- salinity reaches this model "
     "through four constants and not one, so a bracket arm moves the mixed-"
     "layer heat capacity and the snow-ice flooding threshold as well as the "
     "freezing point. world-9hb. OCN-9 corrected the old one-signed inference: "
     "exorheic delivery per unit ocean area ranges from 0.51x Earth's at zero "
     "carving through about 0.91x at the median basin threshold to 2.14x in the "
     "all-carve spill graph, before ocean-age sources and sinks. Fresher water "
     "freezes warmer and saltier water colder, so the sea-ice sign follows the "
     "accepted salinity rather than endorheism alone. The span is about 1.1 K "
     "on one threshold. Archive CLIM-35 measured the S=20 arm growing seasonal "
     "ice by 0.164% of the planet; its -0.105 W/m2 flux was below the two-orbit "
     "resolution, while the S=30 arm moved -0.474 W/m2 without an ice change. "
     "The mechanism is measured and the two-sided planetary prior remains "
     "declared rather than inferred (setting the key was archive CLIM-17)."),
    ("sea ice does not move", "ice too extensive, sign taken from the peer",
     "STRUCTURAL, and DECLARED here rather than priced: icemod.f90 runs "
     "thermodynamics with no advection, transport, drift or dynamics -- ice "
     "forms and melts in place. This is the configuration BIG-MITgcm runs and "
     "names as a stated limitation, reporting excessive ice as its consequence "
     "at a comparable resolution, so the direction comes from a peer report "
     "rather than a guess and the term is not negligible on that evidence. No "
     "kelvin is offered because none can be honestly computed here: the "
     "correction requires an ocean surface velocity this slab does not have. "
     "The eventual fix is small -- cGENIE's is a second-order explicit "
     "transport on two prognostic fields using upper-ocean velocity, with a "
     "diffusion term beside it, and no ice momentum equation at all -- so "
     "declaring the term is not deferring an expensive decision. It is "
     "one-signed against the ice-albedo feedback in the same direction as "
     "PHYS-14's spectral error, and the two are independent. OCN-21."),
    ("roughness distribution", "land median 0.502 m under a 2.0 m mean",
     "Anchored to ExoPlaSim's tuned land mean, which the distribution says is "
     "carried by a rough tail. Anchoring inflates mid-range cells; direction "
     "declared, magnitude not measured."),
    ("biosphere reaches the climate only through albedo",
     "transpiration worth about one runoff",
     "STRUCTURAL one-way coupling, notes/audits/missed-couplings.md finding 4. "
     "ExoPlaSim's land surface is one bucket: a vegetated cell and a bare cell "
     "with the same bucket state evaporate identically apart from roughness. "
     "There is no stomatal or LAI control and no rooting depth. pedology/ "
     "records from Lapides et al. (2024) that a bedrock vadose zone raises "
     "annual transpiration by 100-150 mm, against a runoff of the same size. "
     "The direction is not obvious -- stomatal closure cuts E, deep roots raise "
     "it -- and the magnitude is amplified by E/R into the carve."),
]


def albedo_to_kelvin(d_land, land_fraction, planetary_albedo, attenuation):
    d_toa = d_land * land_fraction * attenuation
    return sensitivity.planetary_albedo_to_kelvin(d_toa, planetary_albedo)


def land_fraction(config) -> float:
    """From the build manifest, by surface_class, never from a literal."""
    manifest = builds.mesh_export(config) / "manifest.json"
    m = json.loads(manifest.read_text(encoding="utf-8"))
    return float(m["landSeaMask"]["landFractionBySurfaceClass"])


def land_water_balance(config) -> dict:
    """P, E and the runoff residual over land, and what each amplifies into.

    `mrro` is deliberately not used. It is river-routed net divergence rather
    than local generation, and `carve_verdict.py` says why at length; the
    criterion divides by P - E, so that is what is budgeted here.
    """
    from netCDF4 import Dataset
    from numpy.polynomial.legendre import leggauss

    # Through the one resolver, which raises when no baseline is named; the
    # str() here turned a null into the literal path "None".
    #
    # IT ALSO REFUSES A CLIMATOLOGY ON THE WRONG RUNG, which is what keeps this
    # half and the albedo half describing one world. The albedo half resolves
    # its report through `rungs.model_grid(config)` and this one takes its grid
    # from the file, so before `require_configured_grid` existed a budget could
    # combine a T21 hydrology with a T42 albedo and say nothing.
    path = paths.climatology_path(root=ROOT)
    with Dataset(path) as ds:
        nlat, nlon = len(ds.dimensions["lat"]), len(ds.dimensions["lon"])
        w = leggauss(nlat)[1][::-1][:, None] * np.ones((1, nlon))
        land = climatology.annual_mean_of(ds, "lsm") > 0.5
        wl = w * land
        def land_mean(name):
            v = climatology.annual_mean_of(ds, name)
            return float((v * wl).sum() / wl.sum())
        # m/s to mm per Earth year. Annualised to Earth years, not orbits,
        # because every other rate in this budget is.
        scale = orbit_lib.EARTH_CALENDAR_YEAR_DAYS * 86400.0 * 1000.0
        precipitation = land_mean("pr") * scale
        evaporation = abs(land_mean("evap")) * scale
        run_id = getattr(ds, "vesper_run_id", None)
    runoff = precipitation - evaporation
    return {
        "source": str(Path(path).relative_to(ROOT)),
        "run_id": run_id,
        # Recorded so the product says which rung it was computed on, rather
        # than leaving that to be inferred from a filename that does not carry
        # it. It is the configured one because the resolver refuses any other.
        "grid": f"{nlat}x{nlon}",
        "land_precipitation_mm_per_earth_year": round(precipitation, 1),
        "land_evaporation_mm_per_earth_year": round(evaporation, 1),
        "runoff_residual_mm_per_earth_year": round(runoff, 1),
        "runoff_fraction_of_precipitation": round(runoff / precipitation, 4),
        "d_runoff_per_d_precipitation": round(precipitation / runoff, 2),
        "d_runoff_per_d_land_evaporation": round(-evaporation / runoff, 2),
        "note": "Fractional. A 1% change in land precipitation moves runoff by "
                "d_runoff_per_d_precipitation percent, and one in land "
                "evaporation by d_runoff_per_d_land_evaporation percent. The "
                "carve criterion divides by this runoff.",
    }


def basin_response(config) -> tuple[dict, int | None, "callable | None"]:
    """How many basins change verdict under a uniform fractional perturbation.

    The criterion is `export_carve_list.py`'s, and its three water terms are
    READ from the carve list rather than solved for:

        runoff = max(P - E_land, 0)
        Q      = runoff * (catchment - area_at_spill) - (E_lake - P) * area_at_spill

    **Read, because the previous version inverted.** It recovered P and E_lake
    from the recorded aridity index, evaporation margin and critical index,
    which is exact algebra and still the wrong design: it wrote the verdict's
    formula down a second time, so it went wrong every time the criterion moved,
    which is every iteration. It did -- the Penman land-rate floor came out, the
    evaporation reference level was corrected -- and the reconstruction then
    missed six basins while continuing to print a table. The generator now
    records the terms it means, and this reads them.

    **Two checks that can fail**, both enforced rather than described:

      * the recorded runoff must equal `max(P - E_land, 0)`. That is the
        verdict's own definition of runoff evaluated against fields it recorded
        separately, so it fails if the runoff source is not `p_minus_e` or if
        the catchment aggregation is not the linear one assumed here.
      * the criterion re-evaluated at zero perturbation must reproduce the
        recorded `overflow_km3_per_year` per basin, and hence the overflow
        count exactly.

    Both are identities of the recorded artifact rather than agreements between
    two formulations, which is the standard `docs/src/practice/failure-modes.md` class 17
    asks for.

    Basins carried forward from a previous carve pass have no water balance to
    read, because they were decided on a terrain that no longer exists. They are
    excluded and counted, not silently folded in.

    **When the carve list cannot be read this REFUSES rather than raising**, and
    the summary it returns carries `unavailable` instead of the tables. The
    build the budget is configured on is often a build the carve list has not
    been regenerated for -- a carve pass writes the list on the build it was
    decided from, and the build it produces gets one only when hydrography is
    taken round again -- and the albedo half of this budget reads no hydrography
    artifact at all. Dying here denied a reader every column that was fine.
    """
    path = builds.component_data("hydrography", config, strict=True) / "carve_list.json"

    def unavailable(reasons):
        return {"source": str(path.relative_to(ROOT)),
                "unavailable": reasons}, None, None

    if not path.exists():
        return unavailable([
            f"{path.relative_to(ROOT)} does not exist, so there is no basin "
            "population to evaluate the criterion against. Re-run "
            "hydrography/scripts/export_carve_list.py on this build."])
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = [r for r in data["basins"] if not r.get("carried_from_previous_pass")]
    carried = len(data["basins"]) - len(rows)
    if not rows:
        return unavailable([
            f"every one of the {carried} basins in "
            f"{path.relative_to(ROOT)} is carried from a previous carve pass, "
            "so none of them has a water balance on this terrain to "
            "re-evaluate the criterion on."])

    required = ("precipitation_km_per_year", "lake_evaporation_km_per_year",
                "land_evaporation_km_per_year")
    absent = [k for k in required if k not in rows[0]]
    if absent:
        return unavailable([
            f"{path.relative_to(ROOT)} does not record {', '.join(absent)}, so "
            "the carve criterion cannot be re-evaluated from it. This budget no "
            "longer inverts the aridity index to recover them; see BUDG-5. "
            "Re-run hydrography/scripts/export_carve_list.py, which writes them."])

    def column(key):
        return np.array([np.nan if r.get(key) is None else r[key] for r in rows],
                        dtype=float)

    catchment = column("catchment_km2")
    spill = column("area_at_spill_km2")
    runoff = column("runoff_km_per_year")
    overflow = column("overflow_km3_per_year")
    precip = column("precipitation_km_per_year")
    lake_e = column("lake_evaporation_km_per_year")
    land_e = column("land_evaporation_km_per_year")
    # `max(catchment - area_at_spill, 0)`, the same clamp the generator applies.
    # Taking the difference raw is a different criterion on any basin whose lake
    # at spill covers its whole catchment.
    dry = np.maximum(catchment - spill, 0.0)

    def discharge(d_precip=0.0, d_land_e=0.0, d_lake_e=0.0):
        p = (1.0 + d_precip) * precip
        r = np.maximum(p - (1.0 + d_land_e) * land_e, 0.0)
        e = (1.0 + d_lake_e) * lake_e
        return r * dry - (e - p) * spill

    def overflowing(d_precip=0.0, d_land_e=0.0, d_lake_e=0.0) -> int:
        return int((discharge(d_precip, d_land_e, d_lake_e) > 0.0).sum())

    # The two checks, and their tolerances are the carve list's own rounding
    # rather than a slack fitted until they passed. A depth is written to 8
    # decimal places of km/year and an area to 2 of km2, so the largest
    # disagreement the rounding alone can produce is a per-basin quantity and it
    # is computed as one: a basin with a 2.7e6 km2 dry catchment carries 0.013
    # km3/year of it, and a small basin carries nothing. A single global epsilon
    # would have to be set by the largest basin and would then wave through any
    # real error on every smaller one.
    # Half the last recorded digit of a depth, an area and a discharge. Each
    # difference of two recorded values carries twice one of them, and `max` is
    # 1-Lipschitz so the clamp does not widen the runoff term.
    #
    #   |dq| <= 2 h_depth (dry + spill) + 2 h_area runoff + h_area |E - P| + h_q
    HALF_DEPTH = 0.5e-8      # km
    HALF_AREA = 0.5e-2       # km2
    HALF_DISCHARGE = 0.5e-6  # km3/year
    runoff_residual = float(np.nanmax(np.abs(
        runoff - np.maximum(precip - land_e, 0.0))))
    q0 = discharge()
    q_tolerance = (2.0 * HALF_DEPTH * (dry + spill)
                   + 2.0 * HALF_AREA * runoff
                   + HALF_AREA * np.abs(lake_e - precip)
                   + HALF_DISCHARGE)
    q_excess = float(np.nanmax(np.abs(q0 - overflow) - q_tolerance))
    discharge_residual = float(np.nanmax(np.abs(q0 - overflow)))
    baseline = int((overflow > 0).sum())
    reconstructed = overflowing()

    channels = {}
    for label, keys in (("land_precipitation", (1, 0, 0)),
                        ("land_evaporation", (0, 1, 0)),
                        ("lake_evaporation", (0, 0, 1)),
                        ("both_evaporations", (0, 1, 1))):
        channels[label] = {
            f"{x:+.0%}": overflowing(*[x * k for k in keys]) - baseline
            for x in PERTURBATIONS
        }
    summary = {
        "source": str(path.relative_to(ROOT)),
        "terrain_hash": data["terrain_hash"],
        "overflowing_basins": baseline,
        "basins_total": len(data["basins"]),
        "basins_evaluated": len(rows),
        "basins_carried_from_a_previous_pass": carried,
        "basins_without_catchment_runoff": int((runoff <= 0).sum()),
        "reconstruction_check": "reproduces the recorded overflow count"
                                if reconstructed == baseline else
                                f"MISMATCH: {reconstructed} against {baseline}",
        "runoff_identity_check":
            "runoff == max(P - E_land, 0)" if runoff_residual <= 2e-8 else
            f"VIOLATED by up to {runoff_residual:.3e} km/year; the carve list "
            "was not written with runoff_source p_minus_e",
        "discharge_check":
            "the criterion re-evaluated from the recorded terms reproduces "
            "overflow_km3_per_year within the recorded rounding" if q_excess <= 0
            else f"EXCEEDS the rounding bound by {q_excess:.3e} km3/year, so "
                 "this file's criterion and export_carve_list.py's have "
                 "diverged",
        "max_discharge_residual_km3_per_year": round(discharge_residual, 8),
        "change_in_overflowing_basins": channels,
        "note": "Uniform fractional perturbations of the catchment-mean fields "
                "the carve list records. A real climate change is patterned and "
                "these are not a substitute for one; they are the criterion's "
                "own sensitivity, which is what a budget item has to be "
                "multiplied through. The response is strongly asymmetric "
                "because the basin population piles up just above the "
                "threshold.",
        "reading_note": "Prefer both_evaporations to land_evaporation alone. "
                        "Land evaporation enters only through the runoff "
                        "denominator and lake evaporation only through the "
                        "numerator, and a kelvin moves both; the single-channel "
                        "rows are the criterion's partial derivatives and not "
                        "climates.",
    }
    return summary, baseline, overflowing


def per_item_carve_currency(kelvin, water, baseline, overflowing):
    """A kelvin into a runoff percent and a range of basins.

    Every channel moves at once, because one kelvin moves all of them, and the
    basin count is the criterion re-evaluated at the JOINT perturbation rather
    than read off the per-channel table. The channels partly cancel: a
    precipitation change is in the criterion's numerator as well as its
    denominator, and a warmer world evaporates more from both the ground and the
    lake.

    The runoff figure is a single number and the basin figure is a range. Runoff
    is `P - E_land` and contains no lake evaporation, so the one bracketed
    channel cannot reach it; it reaches the basin count through the numerator
    only. Returning a single basin number would hide exactly the term that is
    not measured.

    **The two columns fail separately, because they read different things.** The
    runoff percent needs the configured build's water balance and the
    hydrological response; the basin count needs a carve list on top of both. So
    a build whose carve list has not been regenerated still gets a runoff column
    and reports `--` for basins alone.
    """
    if kelvin is None:
        return None, None
    h = HYDROLOGICAL_RESPONSE_PER_KELVIN
    d_precip = h["precipitation"] * kelvin
    d_land_e = h["land_evaporation"] * kelvin
    runoff_percent = 100.0 * (water["d_runoff_per_d_precipitation"] * d_precip
                              + water["d_runoff_per_d_land_evaporation"] * d_land_e)
    if overflowing is None:
        return round(runoff_percent, 2), None
    basins = sorted(int(overflowing(d_precip, d_land_e, rate * kelvin) - baseline)
                    for rate in h["lake_evaporation"])
    return round(runoff_percent, 2), basins


def runoff_per_kelvin(water) -> float:
    """Percent change in land runoff per kelvin, from the same one input.

    Written out because it is the number the trap turns on. It is NOT
    `d_runoff_per_d_precipitation` times the precipitation response: land
    evaporation rises slightly faster than precipitation on this world, and the
    two amplified terms very nearly cancel.
    """
    h = HYDROLOGICAL_RESPONSE_PER_KELVIN
    return 100.0 * (water["d_runoff_per_d_precipitation"] * h["precipitation"]
                    + water["d_runoff_per_d_land_evaporation"] * h["land_evaporation"])


def main() -> None:
    # Parsed first, so `--help` needs no artifacts. See the note in
    # lib/sensitivity.py:main -- both were answering `--help` by running to
    # completion, which only worked while a climatology existed.
    argparse.ArgumentParser(
        description="The error budget: what each unpriced term is worth in "
                    "kelvin and in carve currency. Takes no arguments. Writes "
                    "analysis/error_budget.json; needs the baseline "
                    "climatology named in config/planet.yaml.").parse_args()
    config = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    sensitivity.test_identity()
    problems = sensitivity.verify()

    conversion = sensitivity.provenance(config)
    planetary_albedo = conversion["planetary_albedo"]
    fraction = land_fraction(config)
    water = land_water_balance(config)
    basins, basin_baseline, overflowing = basin_response(config)

    # THE CARVE HALF HAS TWO UNAVAILABILITIES AND THEY ARE NOT THE SAME ONE.
    # Each output on that side is reported exactly when the thing it reads is
    # present, rather than all of them being gated on one flag:
    #
    #   runoff percent   the configured build's water balance and the
    #                    hydrological response. Absent iff `carve_problems`.
    #   basin count      those, AND a carve list on the configured build.
    #                    Absent iff either fails.
    #   the CARVE_ITEMS table and the basin verdict table
    #                    the carve list alone. There is no kelvin anywhere in
    #                    their path, so the response cannot reach them; absent
    #                    iff `basin_problems`.
    #
    # Folding the two together would have printed `--` for a runoff column that
    # is fully determined, and folding them together in the other direction is
    # what made a build with no carve list a traceback. The albedo and forcing
    # halves read no hydrography artifact at all and are never affected by
    # either.
    #
    # A kelvin priced into basins through a response measured on another world
    # is a number, and a number is what makes it dangerous; a basin count taken
    # from a carve list decided on another terrain is the same defect one step
    # earlier. Both report the reason on the artifact instead.
    carve_problems = verify_hydrological_response(config)
    basin_problems = basins.get("unavailable") or []

    def carve(kelvin):
        if carve_problems:
            return None, None
        return per_item_carve_currency(kelvin, water, basin_baseline, overflowing)

    def carve_columns(runoff_pct, basin_range):
        """Both carve columns as strings, so a null prints as a null."""
        runoff_col = "--" if runoff_pct is None else f"{runoff_pct:+.1f}%"
        if basin_range is None:
            return runoff_col, "--"
        low, high = basin_range
        return (runoff_col,
                f"{low:+d}" if low == high else f"{low:+d} to {high:+d}")

    # The shares the items were actually priced through, taken from the same
    # report they came out of rather than restated. WORLD-SU9O.
    land_mean_shares_recorded = {
        name: entry["land_mean_share"]
        for name, entry in (json.loads(
            albedo_report_path(config).read_text(encoding="utf-8")
        ).get("lithology_albedo_overrides") or {}).items()
        if "land_mean_share" in entry}

    rows = []
    for label, delta, note in albedo_items(config):
        naive = albedo_to_kelvin(delta, fraction, planetary_albedo, 1.0)
        atten = albedo_to_kelvin(delta, fraction, planetary_albedo,
                                 DEFAULT_ATTENUATION)
        rows.append((label, delta, naive, atten, note))
    rows.sort(key=lambda r: -abs(r[3]))

    forcings = []
    for label, (low, high), note in FORCING_ITEMS:
        k_low = sensitivity.forcing_to_kelvin(low, planetary_albedo, config)
        k_high = sensitivity.forcing_to_kelvin(high, planetary_albedo, config)
        forcings.append((label, (low, high), (k_low, k_high), note))
    forcings.sort(key=lambda r: -max(abs(r[2][0]), abs(r[2][1])))

    print(f"slope {conversion['slope_k_per_flux_ratio']} K per unit flux ratio, "
          f"{conversion['kelvin_per_w_m2_absorbed']} K per W/m2 absorbed, land "
          f"{fraction:.4f}, planetary albedo {planetary_albedo:.4f}, "
          f"attenuation {DEFAULT_ATTENUATION} {ATTENUATION_SPREAD}, measured on "
          f"{len(ATTENUATION_PAIRS)} pairs")
    if problems:
        print("SENSITIVITY DISAGREES WITH THE RUN INDEX: " + "; ".join(problems))
    attenuation_problems = verify_attenuation(planetary_albedo)
    if attenuation_problems:
        print("THE ATTENUATION NO LONGER MATCHES THE ARMS IT WAS MEASURED ON, "
              "so every kelvin in the albedo table below is stale:")
        for problem in attenuation_problems:
            print(f"  {problem}")
    if carve_problems:
        print("CARVE CURRENCY UNAVAILABLE, so every runoff and basin column "
              "below reads --:")
        for problem in carve_problems:
            print(f"  {problem}")
    if basin_problems:
        print("BASIN CRITERION UNAVAILABLE, so every basin column below reads "
              "--, and the two basin tables are not printed. The albedo, "
              "forcing and runoff columns read no carve list and are "
              "unaffected:")
        for problem in basin_problems:
            print(f"  {problem}")
    print()
    print(f"{'item':38} {'d(alb)':>8} {'K naive':>8} {'K':>7} "
          f"{'runoff':>8} {'basins':>16}")
    print("-" * 92)
    for label, delta, naive, atten, _ in rows:
        r_col, b_col = carve_columns(*carve(atten))
        print(f"{label:38} {delta:+8.4f} {naive:+8.2f} {atten:+7.2f} "
              f"{r_col:>8} {b_col:>16}")

    print(f"\n{'item':46} {'W/m2':>18} {'K':>18} {'runoff':>8} {'basins':>16}")
    print("-" * 110)
    for label, (low, high), (k_low, k_high), _ in forcings:
        r_col, b_col = carve_columns(*carve(k_high))
        print(f"{label:46} {low:+7.2f} to {high:+6.2f} "
              f"{k_low:+8.2f} to {k_high:+6.2f} {r_col:>8} {b_col:>16}")

    print(f"\n{'item, priced in the criterion\'s own channel':46} "
          f"{'perturbation':>22} {'basins':>16}")
    print("-" * 86)
    carve_rows = []
    for label, channel, (low, high), note in CARVE_ITEMS:
        if basin_problems:
            counts = None
        else:
            keys = CARVE_CHANNELS[channel]
            counts = sorted(overflowing(*[x * k for k in keys]) - basin_baseline
                            for x in (low, high))
        carve_rows.append((label, channel, (low, high), counts, note))
        rendered = ("--" if counts is None
                    else f"{counts[0]:+d} to {counts[1]:+d}")
        print(f"{label:46} {channel} {low:+.0%} to {high:+.0%} "
              f"{rendered:>16}")

    print("\nnot in one currency, no conversion invented:")
    for label, magnitude, _ in OTHER_ITEMS:
        print(f"  {label:44} {magnitude}")

    print(f"\nland water balance, mm per Earth year: "
          f"P {water['land_precipitation_mm_per_earth_year']}, "
          f"E {water['land_evaporation_mm_per_earth_year']}, "
          f"runoff {water['runoff_residual_mm_per_earth_year']} "
          f"({water['runoff_fraction_of_precipitation']:.1%} of P). "
          f"Runoff amplifies dP by {water['d_runoff_per_d_precipitation']}x "
          f"and dE by {abs(water['d_runoff_per_d_land_evaporation'])}x.")
    if basin_problems:
        print("\nbasins that change verdict: NOT REPORTED. "
              + " ".join(basin_problems))
    else:
        print(f"\nbasins that change verdict, from "
              f"{basins['overflowing_basins']} overflowing "
              f"({basins['reconstruction_check']}):")
        header = "  " + f"{'channel':20}" + "".join(
            f"{k:>8}" for k in
            basins["change_in_overflowing_basins"]["land_precipitation"])
        print(header)
        for label, table in basins["change_in_overflowing_basins"].items():
            print("  " + f"{label:20}"
                  + "".join(f"{v:+8d}" for v in table.values()))
    h = HYDROLOGICAL_RESPONSE_PER_KELVIN
    lo, hi = h["lake_evaporation"]
    if carve_problems:
        print(f"\nper kelvin: NOT REPORTED. The measured channels come from "
              f"{HYDROLOGICAL_RESPONSE_SECANT['cold']['run_id']} and "
              f"{HYDROLOGICAL_RESPONSE_SECANT['warm']['run_id']} on build "
              f"{HYDROLOGICAL_RESPONSE_SECANT['source_build']}, and the "
              "conversion does not hold here. Re-measure by "
              + HYDROLOGICAL_RESPONSE_SECANT["re_measured_by"] + ".")
    else:
        print(f"\nper kelvin: land P {h['precipitation']:+.2%}, land E "
              f"{h['land_evaporation']:+.2%}, lake E {lo:+.2%} to {hi:+.2%} "
              f"(bracketed), so runoff {runoff_per_kelvin(water):+.2f}%.")
        print(f"  NOT {water['d_runoff_per_d_precipitation'] * h['precipitation'] * 100:+.2f}%: "
              "the amplification applies to the precipitation channel alone, and "
              "a kelvin also raises land evaporation, slightly faster. The two "
              "amplified terms nearly cancel.")

    print("\nRefine an input when its plausible range exceeds the effect of the")
    print("thing you last refined. Everything under a kelvin here is below the")
    print("structural items and below the biosphere question.")

    payload = {
        "note": "Generated by scripts/error_budget.py; DO NOT EDIT, re-run the "
                "generator. A hand edit here is silently reverted the next time "
                "anything regenerates it, which has happened. Ordering, not precision; "
                "a factor of two is expected and is enough to rank refinements.",
        "conversion": dict(conversion, land_fraction=fraction,
                           attenuation=DEFAULT_ATTENUATION,
                           attenuation_spread=list(ATTENUATION_SPREAD),
                           attenuation_note="A surface albedo change does not "
                                            "reach the top of the modelled "
                                            "atmosphere intact. Skipping this "
                                            "step doubles every answer, which "
                                            "this project did once and paid a "
                                            "kelvin for. MEASURED on the pairs "
                                            "below, against STAGED land-mean "
                                            "deltas, which is the quantity "
                                            "every albedo item here is. "
                                            "notes/audits/albedo-attenuation.md.",
                           attenuation_pairs={
                               k: dict(v) for k, v in ATTENUATION_PAIRS.items()},
                           attenuation_verify=(attenuation_problems
                                               or "re-derives from its arms"),
                           land_mean_shares=land_mean_shares_recorded,
                           land_mean_shares_note=
                               "How far the STAGED land-mean albedo moves per "
                               "unit of a rock class's albedo, which is what an "
                               "item denominated in that class has to be "
                               "multiplied by. Not the class's area fraction of "
                               "land: that is an area-shaped number where a "
                               "staged-mean sensitivity belongs, and it "
                               "overstated the playa item by 1.727. READ FROM "
                               "THE ALBEDO REPORT the items are priced against, "
                               "where build_surface_albedo.py emits it as the "
                               "exact derivative of the land mean in the same "
                               "record; there is no declared copy here to go "
                               "stale against the terrain.",
                           verify=problems or "agrees with the run index"),
        "albedo_items": [],
        "forcing_items": [],
        "other_items": [{"item": l, "magnitude": m, "note": n}
                        for l, m, n in OTHER_ITEMS],
        "carve_items": [
            {"item": l, "channel": c, "fractional_change": [lo, hi],
             "basins": counts, "runoff_percent": None,
             "kelvin": None, "note": n}
            for l, c, (lo, hi), counts, n in carve_rows],
        "carve_items_note":
            "Items already denominated in a fractional change of one of the "
            "criterion's own water terms. They convert to basins exactly, with "
            "no kelvin in the path, and none should be given one. A "
            "lake_evaporation item has a null runoff_percent because runoff is "
            "P - E_land and does not contain it. A null `basins` means the "
            "carve list itself was unreadable; carve_currency.basin_response "
            "carries the reason.",
        "carve_currency": {
            "unavailability_is_per_currency":
                "The carve half has two unavailabilities and they are reported "
                "separately, because its outputs read different things. The "
                "runoff percent needs the configured build's water balance and "
                "the hydrological response, and is null when "
                "hydrological_response_verify lists problems. The basin count "
                "needs those AND a carve list on the configured build, and is "
                "null when either fails. The carve_items and the basin_response "
                "tables have no kelvin in their path at all and need the carve "
                "list alone. The albedo and forcing halves read no hydrography "
                "artifact and are never affected by either.",
            "why": "The temperature is a result this project can revise. The "
                   "carve list is not: it leaves the project, changes the "
                   "terrain, and cannot be undone. An item priced only in "
                   "kelvin cannot be ranked against it.",
            "hydrological_response_per_kelvin":
                {k: (list(v) if isinstance(v, tuple) else v)
                 for k, v in HYDROLOGICAL_RESPONSE_PER_KELVIN.items()},
            "hydrological_response_secant":
                {k: (dict(v) if isinstance(v, dict) else v)
                 for k, v in HYDROLOGICAL_RESPONSE_SECANT.items()},
            "hydrological_response_verify":
                carve_problems or "the named runs still carry the temperatures "
                                  "this secant was divided by, and they sit on "
                                  "the build the water balance sits on",
            "hydrological_response_is_re_derivable":
                "the KELVIN is, from the window means in the run index. The two "
                "fractional responses are NOT: they need per-orbit pr and evap "
                "from raw output that has been archived, so the numerator of "
                "this secant has no surviving input. "
                + HYDROLOGICAL_RESPONSE_SECANT["re_measured_by"]
                + " is what restores it.",
            "hydrological_response_note":
                "Fractional change in each land-mean quantity per kelvin of "
                "global-mean surface temperature. Precipitation and land "
                "evaporation are a secant between the two converged fluxes on "
                "this build; lake evaporation is a bracket, low end the land "
                "rate and high end the 6.7 %/K convexity of saturation vapour "
                "pressure, because the Penman response has not been measured "
                "here. See notes/audits/hydrological-sensitivity.md.",
            "runoff_percent_per_kelvin":
                None if carve_problems else round(runoff_per_kelvin(water), 2),
            "runoff_percent_per_kelvin_note":
                "The number the second currency turns on, and the one most "
                "easily got wrong. It is NOT d_runoff_per_d_precipitation times "
                "the precipitation response: a kelvin raises land evaporation "
                "slightly faster than precipitation, and the two amplified "
                "terms nearly cancel. The amplification applies to a "
                "precipitation-only perturbation, which a dust or circulation "
                "item can be and a temperature item never is.",
            "basins_is_a_range_because":
                "lake evaporation per kelvin is bracketed rather than measured. "
                "Runoff is P - E_land and contains no lake evaporation, so only "
                "the basin column carries the bracket.",
            "land_water_balance": water,
            "basin_response": basins,
        },
        "stopping_rule": "Refine an input when its plausible range exceeds the "
                         "effect of the thing you last refined.",
    }
    for label, delta, naive, atten, note in rows:
        runoff_pct, basin_range = carve(atten)
        payload["albedo_items"].append(
            {"item": label, "delta_land_albedo": delta,
             "kelvin_naive": round(naive, 3), "kelvin": round(atten, 3),
             "runoff_percent": runoff_pct, "basins": basin_range, "note": note})
    for label, (low, high), (k_low, k_high), note in forcings:
        runoff_pct, basin_range = carve(k_high)
        payload["forcing_items"].append(
            {"item": label, "w_m2": [low, high],
             "kelvin": [round(k_low, 3), round(k_high, 3)],
             "runoff_percent": runoff_pct, "basins": basin_range, "note": note})

    out = ROOT / "analysis" / "error_budget.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
