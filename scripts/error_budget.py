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
fractional perturbation of precipitation, land evaporation and lake evaporation.
It READS the three water terms the carve list records against each basin and
re-evaluates the criterion on them. It used to invert the recorded aridity index
and evaporation margin to recover them instead, which put the verdict's formula
in a second place and broke the first time the criterion moved; BUDG-5.

**The kelvin reaches those columns through
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
sys.path.insert(0, str(ROOT / "lib"))

import numpy as np  # noqa: E402
import yaml  # noqa: E402

import builds  # noqa: E402
import orbit as orbit_lib  # noqa: E402
import climatology
import paths  # noqa: E402
import sensitivity  # noqa: E402

# Measured, not assumed; see the module docstring.
DEFAULT_ATTENUATION = 0.5

# BUDG-1. Fractional change in each land-mean quantity per kelvin of global-mean
# surface temperature. This is the whole of the kelvin-to-carve conversion; the
# derivation and its evidence are in `notes/audits/hydrological-sensitivity.md`.
#
# `precipitation` and `land_evaporation` are MEASURED, as a secant between the
# two converged fluxes that share a surface on this build, over the last ten
# orbits of each. No run was bought for them. They are a secant across about
# seven kelvin rather than a local slope, and design intent forbids reusing a
# sensitivity measured in one regime in another, so the converged points loop
# A's flux re-bracket produces are what localises them. That is a step, not a
# task.
#
# They must be measured with the SAME land mean `land_water_balance` uses, since
# the two are multiplied together: Gaussian quadrature weights, land from
# `lsm > 0.5`, annualised to Earth years. A response measured under a different
# land mean does not compose with an amplification measured under this one.
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


ALBEDO_REPORT = ROOT / "exoplasim" / "inputs" / "t42" / "albedo_report.json"


def albedo_items() -> list[tuple[str, float, str]]:
    """(label, land-mean albedo delta, note) for the albedo half of the budget.

    Two of these are MEASURED by `build_surface_albedo.py` and are read from the
    report it writes rather than transcribed. They were literals until 2026-08-18
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
    if "lakes" not in report:
        raise SystemExit(
            f"{ALBEDO_REPORT} has no lakes block, so it was built without "
            "--lakes and the lake item cannot be read from it."
        )
    return [
        ("biosphere: bare rock vs vegetated",
         -(float(ends["bare_rock"]) - float(ends["vegetated"])),
         "the assumption LPJ-GUESS exists to replace. Measured: "
         f"{ends['bare_rock']:.5f} bare against {ends['vegetated']:.5f} "
         "vegetated, both pre-lake and on the overridden rock table."),
        ("playa_clastic albedo, 0.25 to 0.33", 0.08 * 0.147,
         "largest single rock-class lever; a mixture of clay playa and varnished fan"),
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
     "broadcasts both, and acts on `nhdiff > 0` at line 844, so horizontal "
     "heat diffusion is one namelist key on the binary already built. That "
     "satisfies every one of `docs/src/pipeline/sequencing.md` A3's four A/B conditions with no further "
     "work, which makes this the CHEAPEST structural item to price rather than "
     "the impossible one. Read it as physics-is-not-a-knob: ocean heat "
     "transport exists, and the argument for switching it on is not that it "
     "improves an agreement. Bracket `hdiffk` and report the spread rather "
     "than tuning it -- a constant diffusivity is a BOUND on the missing "
     "transport, not the transport. CLIM-16."),
    ("ocean salinity, 20 to 40 psu", "freezing point 272.07 to 270.94 K",
     "DECLARED, not derived: config/planet.yaml carries Earth's 34.7 because no "
     "salt budget exists here, and the model's TFREEZE is now set from it "
     "rather than compiled in. OCN-9 corrected the old one-signed inference: "
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

    # Through the one resolver, which raises when no baseline is named; the
    # str() here turned a null into the literal path "None".
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
    """
    path = builds.component_data("hydrography", config, strict=True) / "carve_list.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = [r for r in data["basins"] if not r.get("carried_from_previous_pass")]
    carried = len(data["basins"]) - len(rows)

    required = ("precipitation_km_per_year", "lake_evaporation_km_per_year",
                "land_evaporation_km_per_year")
    absent = [k for k in required if k not in rows[0]]
    if absent:
        raise SystemExit(
            f"{path.relative_to(ROOT)} does not record {', '.join(absent)}, so "
            "the carve criterion cannot be re-evaluated from it. This budget no "
            "longer inverts the aridity index to recover them; see BUDG-5. "
            "Re-run hydrography/scripts/export_carve_list.py, which writes them.")

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
    """
    if kelvin is None:
        return None, None
    h = HYDROLOGICAL_RESPONSE_PER_KELVIN
    d_precip = h["precipitation"] * kelvin
    d_land_e = h["land_evaporation"] * kelvin
    runoff_percent = 100.0 * (water["d_runoff_per_d_precipitation"] * d_precip
                              + water["d_runoff_per_d_land_evaporation"] * d_land_e)
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

    def carve(kelvin):
        return per_item_carve_currency(kelvin, water, basin_baseline, overflowing)

    def carve_columns(runoff_pct, basin_range):
        """Both carve columns as strings, so a null prints as a null."""
        if runoff_pct is None:
            return "--", "--"
        low, high = basin_range
        return (f"{runoff_pct:+.1f}%",
                f"{low:+d}" if low == high else f"{low:+d} to {high:+d}")

    rows = []
    for label, delta, note in albedo_items():
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
        keys = CARVE_CHANNELS[channel]
        counts = sorted(overflowing(*[x * k for k in keys]) - basin_baseline
                        for x in (low, high))
        carve_rows.append((label, channel, (low, high), counts, note))
        print(f"{label:46} {channel} {low:+.0%} to {high:+.0%} "
              f"{f'{counts[0]:+d} to {counts[1]:+d}':>16}")

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
    h = HYDROLOGICAL_RESPONSE_PER_KELVIN
    lo, hi = h["lake_evaporation"]
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
            "P - E_land and does not contain it.",
        "carve_currency": {
            "why": "The temperature is a result this project can revise. The "
                   "carve list is not: it leaves the project, changes the "
                   "terrain, and cannot be undone. An item priced only in "
                   "kelvin cannot be ranked against it.",
            "hydrological_response_per_kelvin":
                {k: (list(v) if isinstance(v, tuple) else v)
                 for k, v in HYDROLOGICAL_RESPONSE_PER_KELVIN.items()},
            "hydrological_response_note":
                "Fractional change in each land-mean quantity per kelvin of "
                "global-mean surface temperature. Precipitation and land "
                "evaporation are a secant between the two converged fluxes on "
                "this build; lake evaporation is a bracket, low end the land "
                "rate and high end the 6.7 %/K convexity of saturation vapour "
                "pressure, because the Penman response has not been measured "
                "here. See notes/audits/hydrological-sensitivity.md.",
            "runoff_percent_per_kelvin": round(runoff_per_kelvin(water), 2),
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
