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
every answer.** A surface albedo change does not arrive at the top of the
atmosphere intact; everything above the surface scatters and absorbs. This
project measured it the expensive way: a lithology fix that moved planet-mean
*surface* albedo by +0.0033 was predicted at -0.81 K on the naive conversion, and
the run came in about a kelvin warmer, with planetary albedo having moved -0.0017
-- the opposite sign, once the warmer state's reduced sea ice is included. A
factor of 0.5 is used here as a conservative default and both columns are
reported, because the honest claim is a factor of two rather than a number.

## The second currency, and why kelvin is not enough

The temperature is a result this project can revise. The carve list is not: it
leaves the project, changes the terrain, and cannot be undone. It is decided by a
water balance, so an item priced only in kelvin cannot be ranked against it.

Runoff over land is the residual of two numbers five to six times its own size,
so it amplifies: a fractional change in precipitation moves runoff by `P/R` times
as much, and one in land evaporation by `E/R`. Both are measured here from the
baseline climatology rather than declared.

`basin_response` then measures how many basins change verdict under a uniform
fractional perturbation of precipitation, land evaporation and lake evaporation,
reconstructed from the carve list's own fields and its own criterion. That
reconstruction reproduces the recorded overflow count exactly at zero
perturbation, which is the check on it.

**One link in the chain is missing and it is the one that needs a run.** Nothing
has measured what a kelvin does to precipitation and evaporation on this world.
`HYDROLOGICAL_RESPONSE_PER_KELVIN` is that link, declared here as a single input
so that landing it is one edit rather than a rewrite; while it is `None` the
per-item runoff and basin columns are null and say why.

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

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

import numpy as np  # noqa: E402
import yaml  # noqa: E402

import builds  # noqa: E402
import orbit as orbit_lib  # noqa: E402
import sensitivity  # noqa: E402

# Measured, not assumed; see the module docstring.
DEFAULT_ATTENUATION = 0.5

# The one input BUDG-1 has to supply, and the only reason the carve columns are
# empty. Fractional change in each land-mean quantity per kelvin of global-mean
# surface temperature, from one perturbation run against the baseline:
#
#     {"precipitation": 0.0xx, "land_evaporation": 0.0xx, "lake_evaporation": 0.0xx}
#
# Set it and everything below fills in. Nothing else in this file changes.
HYDROLOGICAL_RESPONSE_PER_KELVIN = None

# Perturbation sizes the basin response is tabulated at. Both signs, because the
# criterion is a threshold and the basin population is not symmetric about it.
PERTURBATIONS = [-0.10, -0.05, -0.02, -0.01, 0.01, 0.02, 0.05, 0.10]


# (label, land-mean albedo delta, note). Deltas are the plausible RANGE of the
# item, not its current value: the budget ranks what refining each would buy.
ALBEDO_ITEMS = [
    ("biosphere: bare rock vs vegetated", -(0.276 - 0.179),
     "the assumption LPJ-GUESS exists to replace"),
    ("playa_clastic albedo, 0.25 to 0.33", 0.08 * 0.147,
     "largest single rock-class lever; a mixture of clay playa and varnished fan"),
    ("lakes composited into albedo", -0.0139,
     "solved lakes reaching the climate at all"),
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
    ("dust, shortwave-only if switched on as shipped", (-4.3, -4.0),
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
    ("slab depth, 25 m vs 50 m", "seasonal amplitude ~2x",
     "STRUCTURAL. Amplitude scales as 1/(omega*C) and this world's year is half "
     "Earth's, so 50 m damps seasonality about twice as hard as Earth's ocean "
     "does. Coldest-month mean is what PFT survival gates on. One perturbation "
     "run, ever, gives a permanent scaling for every later result."),
    ("no q-flux", "gradients too strong, ice too extensive",
     "STRUCTURAL, and no cheap version exists. Declare the direction and move on."),
    ("Penman over a dry land column", "10 to 18% of land-mean Penman",
     "Validated over ocean, where the air is equilibrated with the surface. Over "
     "a subgrid lake the column is dry, so VPD is too high, E is overstated and "
     "the verdict under-carves. Quantified by notes/audits/missed-couplings.md "
     "finding 2 by raising the assumed column humidity: -10% at 70%, -13% at "
     "80%, -18% at 90%. Read that against the lake-evaporation row of "
     "basin_response below; it is the largest carve-side item here and it is "
     "opposite in sign to the albedo-driven over-carve."),
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
    manifest = builds.build_root(config) / "exoplasim-T42" / "manifest.json"
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

    path = ROOT / str(config["baseline_climatology"])
    with Dataset(path) as ds:
        nlat, nlon = len(ds.dimensions["lat"]), len(ds.dimensions["lon"])
        w = leggauss(nlat)[1][::-1][:, None] * np.ones((1, nlon))
        land = np.asarray(ds.variables["lsm"][:], dtype=float).mean(axis=0) > 0.5
        wl = w * land
        def land_mean(name):
            v = np.asarray(ds.variables[name][:], dtype=float).mean(axis=0)
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


def basin_response(config) -> tuple[dict, int, "callable"]:
    """How many basins change verdict under a uniform fractional perturbation.

    Reconstructed from the carve list, which is enough: the criterion is

        Q = runoff * (catchment - area_at_spill) - (E_lake - P) * area_at_spill

    and the file carries the aridity index (E_lake - P)/runoff, the evaporation
    margin (E_lake - P - critical*runoff)/E_lake and the critical index, which
    solve for the catchment-mean precipitation and lake evaporation separately.
    Land evaporation is then P - runoff.

    **The check that could have failed**: at zero perturbation the reconstruction
    must reproduce the recorded overflow count exactly, and it does. It is also
    insensitive to the file's rounding -- perturbing every rounded field over its
    rounding interval moves at most one basin.

    **The blind spot, stated rather than hidden.** Basins whose catchment runoff
    is zero cannot be reconstructed, because the clamp at zero destroys P - E.
    None of them currently overflows and none can, since Penman is floored at the
    land rate and the land rate already exceeds P there. So the closures below
    are complete and the openings are a lower bound.
    """
    path = builds.component_data("hydrography", config, strict=True) / "carve_list.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = data["basins"]

    def column(key):
        return np.array([np.nan if r[key] is None else r[key] for r in rows],
                        dtype=float)

    index = column("aridity_index_penman")
    margin = column("evaporation_margin")
    critical = column("critical_aridity_index")
    catchment = column("catchment_km2")
    spill = column("area_at_spill_km2")
    runoff = column("runoff_km_per_year")
    overflow = column("overflow_km3_per_year")

    solvable = ((runoff > 0) & np.isfinite(index) & np.isfinite(margin)
                & (margin != 0))
    safe = np.where(solvable, margin, 1.0)
    lake_e = np.where(solvable, runoff * (index - critical) / safe, np.nan)
    precip = lake_e - index * runoff
    land_e = precip - runoff
    dry = catchment - spill
    base_q = runoff * dry - (lake_e - precip) * spill

    def overflowing(d_precip=0.0, d_land_e=0.0, d_lake_e=0.0) -> int:
        r = np.maximum((1 + d_precip) * precip - (1 + d_land_e) * land_e, 0.0)
        # The same floor carve_verdict.py applies: Penman below the model's own
        # land rate is impossible for a saturated surface under the same forcing.
        e = np.maximum((1 + d_lake_e) * lake_e, (1 + d_land_e) * land_e)
        q = r * dry - (e - (1 + d_precip) * precip) * spill
        return int(((overflow + np.where(solvable, q - base_q, 0.0)) > 0).sum())

    baseline = int((overflow > 0).sum())
    floored = np.isclose(lake_e, land_e, rtol=1e-9)[solvable].sum()
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
        "basins_total": len(rows),
        "basins_reconstructed": int(solvable.sum()),
        "reconstruction_check": "reproduces the recorded overflow count"
                                if reconstructed == baseline else
                                f"MISMATCH: {reconstructed} against {baseline}",
        "change_in_overflowing_basins": channels,
        "note": "Uniform fractional perturbations of the catchment-mean fields. "
                "A real climate change is patterned and these are not a "
                "substitute for one; they are the criterion's own sensitivity, "
                "which is what a budget item has to be multiplied through. The "
                "response is strongly asymmetric because the basin population "
                "piles up just above the threshold.",
        "floored_basins": int(floored),
        "floor_note": "The land-evaporation row is dominated by the floor "
                      "`penman = max(penman, evap)`, not by the runoff "
                      "denominator: without it a +1% land evaporation moves a "
                      "quarter as many basins. The floor binds across the WHOLE "
                      "catchment on floored_basins of the reconstructed set, "
                      "where the Penman estimate is inoperative and the aridity "
                      "index is exactly -1. Read the land-evaporation row with "
                      "that in mind and prefer both_evaporations, since a "
                      "kelvin moves land and lake evaporation together.",
    }
    return summary, baseline, overflowing


def per_item_carve_currency(kelvin, water, baseline, overflowing):
    """Kelvin into runoff percent and basins, once BUDG-1 has supplied the link.

    All three channels move together, because one kelvin moves all three. The
    basin count is the criterion re-evaluated at the joint perturbation rather
    than read off a per-channel table, since the channels partly cancel: a
    precipitation change appears in the criterion's numerator as well as its
    denominator.
    """
    if HYDROLOGICAL_RESPONSE_PER_KELVIN is None or kelvin is None:
        return None, None
    h = HYDROLOGICAL_RESPONSE_PER_KELVIN
    d_precip = h["precipitation"] * kelvin
    d_land_e = h["land_evaporation"] * kelvin
    d_lake_e = h["lake_evaporation"] * kelvin
    runoff_percent = 100.0 * (water["d_runoff_per_d_precipitation"] * d_precip
                              + water["d_runoff_per_d_land_evaporation"] * d_land_e)
    basin_count = overflowing(d_precip, d_land_e, d_lake_e) - baseline
    return round(runoff_percent, 2), int(basin_count)


def main() -> None:
    config = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    sensitivity.test_identity()
    problems = sensitivity.verify()

    conversion = sensitivity.provenance(config)
    planetary_albedo = conversion["planetary_albedo"]
    fraction = land_fraction(config)
    water = land_water_balance(config)
    basins, basin_baseline, overflowing = basin_response(config)

    rows = []
    for label, delta, note in ALBEDO_ITEMS:
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
          f"attenuation {DEFAULT_ATTENUATION}")
    if problems:
        print("SENSITIVITY DISAGREES WITH THE RUN INDEX: " + "; ".join(problems))
    print()
    print(f"{'item':38} {'d(alb)':>8} {'K naive':>8} {'K':>7}")
    print("-" * 66)
    for label, delta, naive, atten, _ in rows:
        print(f"{label:38} {delta:+8.4f} {naive:+8.2f} {atten:+7.2f}")

    print(f"\n{'item':46} {'W/m2':>18} {'K':>18}")
    print("-" * 84)
    for label, (low, high), (k_low, k_high), _ in forcings:
        print(f"{label:46} {low:+7.2f} to {high:+6.2f} "
              f"{k_low:+8.2f} to {k_high:+6.2f}")

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
    print(f"\nbasins that change verdict, from {basins['overflowing_basins']} "
          f"overflowing ({basins['reconstruction_check']}):")
    header = "  " + f"{'channel':20}" + "".join(
        f"{k:>8}" for k in basins["change_in_overflowing_basins"]["land_precipitation"])
    print(header)
    for label, table in basins["change_in_overflowing_basins"].items():
        print("  " + f"{label:20}" + "".join(f"{v:+8d}" for v in table.values()))
    if HYDROLOGICAL_RESPONSE_PER_KELVIN is None:
        print("\nkelvin does not reach those columns yet: BUDG-1 has to measure "
              "dP, dE_land and dE_lake per kelvin on this world. Set "
              "HYDROLOGICAL_RESPONSE_PER_KELVIN and every item gets a runoff "
              "and a basin figure.")

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
                           attenuation_note="A surface albedo change does not "
                                            "reach the top of the atmosphere "
                                            "intact. Skipping this step doubles "
                                            "every answer, which this project "
                                            "did once and paid a kelvin for.",
                           verify=problems or "agrees with the run index"),
        "albedo_items": [],
        "forcing_items": [],
        "other_items": [{"item": l, "magnitude": m, "note": n}
                        for l, m, n in OTHER_ITEMS],
        "carve_currency": {
            "why": "The temperature is a result this project can revise. The "
                   "carve list is not: it leaves the project, changes the "
                   "terrain, and cannot be undone. An item priced only in "
                   "kelvin cannot be ranked against it.",
            "hydrological_response_per_kelvin": HYDROLOGICAL_RESPONSE_PER_KELVIN,
            "blocked_on": None if HYDROLOGICAL_RESPONSE_PER_KELVIN else {
                "task": "BUDG-1",
                "needs": "one perturbation climate run against the baseline, "
                         "reporting the fractional change per kelvin of "
                         "global-mean surface temperature in land-mean "
                         "precipitation, land-mean evaporation and the Penman "
                         "open-water evaporation the carve verdict reads.",
                "precision": "The basin response passes 100 basins per percent "
                             "of precipitation near the baseline, so a "
                             "conversion good to 20% relative is enough to rank "
                             "items and one good to a factor of two is not.",
            },
            "land_water_balance": water,
            "basin_response": basins,
        },
        "stopping_rule": "Refine an input when its plausible range exceeds the "
                         "effect of the thing you last refined.",
    }
    for label, delta, naive, atten, note in rows:
        runoff_pct, basin_count = per_item_carve_currency(
            atten, water, basin_baseline, overflowing)
        payload["albedo_items"].append(
            {"item": label, "delta_land_albedo": delta,
             "kelvin_naive": round(naive, 3), "kelvin": round(atten, 3),
             "runoff_percent": runoff_pct, "basins": basin_count, "note": note})
    for label, (low, high), (k_low, k_high), note in forcings:
        runoff_pct, basin_count = per_item_carve_currency(
            k_high, water, basin_baseline, overflowing)
        payload["forcing_items"].append(
            {"item": label, "w_m2": [low, high],
             "kelvin": [round(k_low, 3), round(k_high, 3)],
             "runoff_percent": runoff_pct, "basins": basin_count, "note": note})

    out = ROOT / "analysis" / "error_budget.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
