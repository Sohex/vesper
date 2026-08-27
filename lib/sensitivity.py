"""The flux-to-kelvin conversion, in one place, because it was in three.

Anything that turns a radiative quantity into a temperature reads this module.
Three sensitivities were in simultaneous use and they are not interchangeable:
one measured at T21 on a superseded terrain across a bracket that no longer
contains the baseline, one measured at T42 on a superseded build, and one taken
from a stellar sweep that crosses the ice transition. They differ by a factor of
2.2, and the error budget's ranking is computed with whichever one the caller
picked.

Two things are settled here, and they are separate errors.

## 1. The slope

`SLOPE_K_PER_FLUX_RATIO` is dT/df near the baseline flux, where f is the stellar
constant in Earth units. It is a LOCAL slope, measured across two converged
points on the ACTIVE build, and it must not be used outside the regime it was
measured in.

## WHAT RE-RUNS WHEN THIS CHANGES

The corollary in `docs/src/pipeline/loops.md`, asked of this file. Nothing here
drives a model run; what a change to the slope makes worthless is a list of
generated artifacts and one register of predictions, and it is the whole list:

    scripts/error_budget.py                      -> analysis/error_budget.json
    exoplasim/scripts/derive_design_flux.py      -> exoplasim/analysis/design_flux.json
    exoplasim/scripts/shortwave_band_weights.py  -> exoplasim/analysis/shortwave_band_weights.json
    exoplasim/scripts/cloud_optical_depth_bracket.py
                                                 -> exoplasim/analysis/cloud_optical_depth_bracket.json
    exoplasim/scripts/stephens_tables_vs_fits.py -> exoplasim/analysis/stephens_tables_vs_fits.json
    analysis/snow_albedo_zenith.py               -> analysis/snow_albedo_zenith.json
    exoplasim/scripts/predict_ocean_terms.py     -> prints; no artifact
    exoplasim/notes/forcing-bundle-predictions.md   every kelvin registered
                                                 through the conversion

A registered prediction that changes silently is worse than one that was wrong,
so an entry in that register is amended IN PLACE and the amendment says it moved.

## What `verify()` is for, and what it could not do

It recomputes the declared slope from the runs it names. That is a check on
INTERNAL CONSISTENCY. It was not a check on CURRENCY, and the difference cost the
whole measurement: both runs of the previous bracket were deleted, survived only
as archived identity under `archive/runs/`, and carried a build that is no longer
`config/planet.yaml`'s `source_build`. `verify()` read their archived numbers,
reproduced the constant they had produced, and agreed with itself for as long as
nobody asked which world it was describing.

So the currency checks are the load-bearing half. Each bracket run must be in the
LIVE `exoplasim/runs/INDEX.json` rather than only in the archive, and must carry
the `source_build` config names now, resolved through `lib/orogen.py`'s registry.
`test_currency_refuses_a_superseded_measurement()` drives those checks against
the superseded declaration and asserts they fire, so the failure is exercised
rather than asserted.

## 2. The denominator

A top-of-atmosphere forcing is an ABSORBED-flux change, so converting one to a
flux ratio divides by absorbed flux per unit flux ratio, not by incident:

    absorbed per unit f  =  S_earth / 4 * (1 - planetary_albedo)

Dividing by incident insolation instead understates every answer by
1/(1 - planetary_albedo). Dividing by the incident flux AT the baseline rather
than per unit f overstates it by 1/f, in the other direction and much smaller.
Both were being done at once, and the two do not cancel.

The same denominator makes an albedo change and a forcing the same quantity.
Absorbed shortwave is S/4 * (1 - alpha), so a change of `d` in PLANETARY albedo
is worth `-d / (1 - alpha)` in flux ratio. That identity is what
`planetary_albedo_to_kelvin` is; `test_identity()` checks the two routes against
each other, and they are the same arithmetic written twice, so a disagreement is
a coding error rather than a physical claim.

## What this module deliberately does not do

It does not convert a SURFACE albedo change into a planetary one. That
attenuation is a modelled quantity with a factor-of-two uncertainty and it
belongs with the budget that owns the assumption, not here.
"""

from __future__ import annotations

import argparse
import json
import copy
import math
from pathlib import Path

import numpy as np
import yaml

from paths import climatology_path, rel
import climatology as clim  # `climatology` is a parameter name below

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONFIG = PROJECT_ROOT / "config" / "planet.yaml"
RUN_INDEX = PROJECT_ROOT / "exoplasim" / "runs" / "INDEX.json"
ARCHIVED_RUNS = PROJECT_ROOT / "archive" / "runs"
CONVERGENCE_REPORTS = PROJECT_ROOT / "exoplasim" / "analysis" / "convergence"

# --------------------------------------------------------------------------
# The slope. Measured 2026-08-27, on the active build.
# --------------------------------------------------------------------------
#
# Two T21 runs on `canonical-10m-base`, same executable sha, same staged surface
# fields, and configs that differ in nothing but the flux the run was given:
#
#     f = 0.945   run_432e5e46adef   fitted asymptote 279.6805 K +/- 0.0203
#     f = 1.000   run_b45380e61f90   fitted asymptote 288.4654 K +/- 0.0465
#
# giving 159.7 K per unit flux ratio. Both pass all six convergence criteria.
#
# THE ASYMPTOTES ARE THE DECLARED ESTIMATOR and `verify()` recomputes exactly
# this secant, so the tolerance only has to absorb a re-assessment moving a fit
# inside its own stated uncertainty. The two half-widths propagate to
# sqrt(0.0203^2 + 0.0465^2) / 0.055 = 0.92 K per unit flux ratio, which is the
# floor any spread can honestly claim and is what the tolerance is rounded from.
#
# THE SPREAD IS THE ENVELOPE OF THE ESTIMATORS ON THIS BUILD, because it is
# wider than that floor and the estimators disagree by more than either fit:
#
#     159.7   the two fitted asymptotes                      DECLARED
#     159.2   the two settled-window means
#     158.2   the last ten orbits of each assessed window
#     155.3   the bootstrap climatology of the 0.945 run against the last ten
#             orbits of the 1.000 run, which is what `derive_design_flux.py`
#             reads and what `notes/audits/design-flux-two-point-response.md`
#             measures as an 8.54 K span
#
# The low end is the lowest of those; the high end is the declared value plus the
# propagated half-width, since no estimator sits above it. The 4 K width is the
# 0.945 run's last twelve orbits -- the ones its climatology was cut from, added
# after its assessment -- sitting about 0.2 K above its own fitted asymptote.
#
# THE BASELINE IS THE COLD ENDPOINT, NOT INSIDE THE BRACKET. This is a forward
# secant anchored at the design flux and it is honest about that. It is not the
# defect that retired 150.2: that bracket ran 0.95 to 1.00 and excluded 0.945
# entirely, so it described a regime the baseline was not in. What localises this
# one is a converged run BELOW 0.945 on this build, which loop A's flux
# re-bracket produces anyway.
#
# THIS IS A CHORD ACROSS 8.8 K AND IT CARRIES SEA-ICE RETREAT. The modelled sea
# ice mean fraction falls from 0.0799 to 0.0252 between the endpoints, so the
# ice-albedo feedback over that retreat is inside the number rather than outside
# it. Stated because it is the obvious explanation to reach for and it is the
# wrong one: the superseded 202.0 was measured across a pair whose ice barely
# moved and it is the HIGHER value, so the difference between the two is the
# build and the model source, not the ice.
#
# THIS REPLACES 202.0, measured at T42 on `precarve-craton` between f = 0.910 and
# f = 0.945. That build is superseded, both of its runs have been deleted, and
# the model source has moved under it. It is 26% high against this measurement.
#
# THIS IS NOT the stellar sweep's sensitivity. The 0.85-to-0.95 sweep spans 21
# W/m2 absorbed and 33 K, which is 330 K per unit flux ratio: 2.1x this, because
# it crosses the ice transition proper. Do not use one for the other.
SLOPE_K_PER_FLUX_RATIO = 159.7
SLOPE_SPREAD_K_PER_FLUX_RATIO = (155.3, 160.6)

# `report` is the file the asymptote is READ FROM, named rather than derived: a
# diagnostic assessment carries its mode in its filename and must never stand
# where a reader looks for a run's own verdict, so naming it is what keeps the
# two apart. The warm endpoint's orbits are all declared diagnostic, so its
# report is about the experiment and not about the planet's trajectory; the
# fitted asymptote of its temperature series is a property of the series either
# way, and its six criteria are recorded here because that is what makes it
# usable as an endpoint at all.
#
# `window` AND `io_regime` ARE PART OF THE MEASUREMENT, and leaving them out
# once cost this number. A run has more than one admissible window, and the
# asymptote is different in each: the cold endpoint reads 279.6805 K over
# orbits 37-69 and 279.8296 K over the twelve clean-I/O orbits after them, a
# difference of 0.15 K that moves this slope by 2.7. Naming only the run and
# the report file left the window as an argument nobody recorded, so when the
# default assessment moved to the clean block -- which it had to, once a verdict
# window was forbidden from spanning the I/O join -- the file under that name
# quietly became a reading of a different window, and `verify()` recomputed the
# slope from an unmatched pair without noticing.
#
# THE ARMS OF A DIFFERENCE HAVE TO BE READ ON ONE INSTRUMENT AND ONE WINDOW.
# PlaSim's low-I/O accumulation and the clean stream disagree by about 0.17 K on
# this run, which is 2 per cent of the 8.8 K chord and the whole of that 2.7.
# Both endpoints are therefore read over orbits 37-69 of their own 70-orbit
# low-I/O block, the same indices on both, and `_verify_bracket` refuses when a
# report's window is not the one recorded here or when the two arms differ in
# regime. `docs/src/practice/failure-modes.md` class 36.
SLOPE_BRACKET_RUNS = {
    "cold": {"run_id": "run_432e5e46adef", "flux_ratio": 0.945,
             "asymptote_k": 279.6805, "asymptote_half_width_k": 0.0203,
             "report": "run_432e5e46adef_convergence_through070.json",
             "window": (37, 69), "io_regime": "low_io",
             "status": "equilibrated_for_worldbuilding"},
    "warm": {"run_id": "run_b45380e61f90", "flux_ratio": 1.000,
             "asymptote_k": 288.4654, "asymptote_half_width_k": 0.0465,
             "report": "run_b45380e61f90_convergence_diagnostic.json",
             "window": (37, 69), "io_regime": "low_io",
             "status": "equilibrated_for_worldbuilding"},
}
# The propagated half-width, 0.92, rounded up. A recomputation that moves further
# than the two fits' own stated uncertainty is a real change and not fit noise.
SLOPE_TOLERANCE_K_PER_FLUX_RATIO = 1.0

# Both endpoints stage the same surface fields, which is what makes their
# difference a flux response rather than a surface change wearing a slope's
# units. Those fields are TERRAIN-ONLY: `config/planet.yaml`'s
# `baseline_climatology` is null and this build has no baseline, so the bootstrap
# surface is the most determined state that exists here and the slope is measured
# on it. It moves when the baseline surface exists, and this declaration is what
# has to be re-derived then.
SLOPE_GEOGRAPHY = "ecb13b14"

# The declaration this replaced, kept as a FIXTURE and not as a record: it is
# what `test_currency_refuses_a_superseded_measurement()` drives the currency
# checks with. Both runs are deleted and survive only under `archive/runs/`, and
# `precarve-craton` is not the active build, so every currency check has to fire
# on it. The numbers are the archived entries' own and reproduce 201.8.
SUPERSEDED_BRACKET_RUNS = {
    "cold": {"run_id": "run_bfa3f5269660", "flux_ratio": 0.910,
             "asymptote_k": 281.968, "asymptote_half_width_k": None,
             "report": "run_bfa3f5269660_convergence.json",
             "status": "quasi_equilibrated"},
    "warm": {"run_id": "run_524fbed77a9a", "flux_ratio": 0.945,
             "asymptote_k": 289.030, "asymptote_half_width_k": None,
             "report": "run_524fbed77a9a_convergence.json",
             "status": "equilibrated_for_worldbuilding"},
}
SUPERSEDED_SLOPE_K_PER_FLUX_RATIO = 202.0
SUPERSEDED_SPREAD_K_PER_FLUX_RATIO = (201.0, 204.0)
SUPERSEDED_GEOGRAPHY = "3a17c498"


def config(cfg: dict | None = None) -> dict:
    return cfg if cfg is not None else yaml.safe_load(
        CONFIG.read_text(encoding="utf-8"))


def incident_w_m2_per_flux_ratio(cfg: dict | None = None) -> float:
    """Global-mean insolation per unit of stellar flux ratio, S_earth / 4."""
    return float(config(cfg)["orbit"]["earth_solar_constant_w_m2"]) / 4.0


def _gaussian_weights(nlat: int, nlon: int) -> np.ndarray:
    from numpy.polynomial.legendre import leggauss
    return leggauss(nlat)[1][::-1][:, None] * np.ones((1, nlon))


def planetary_albedo_from_fluxes(net_down_w_m2: float,
                                 upward_w_m2: float) -> float:
    """Planetary albedo from globally meaned `rst` and `rsut`.

    Separate from `planetary_albedo` only in where the two fluxes come from, so
    a caller holding a run's own annual means -- rather than a climatology file
    -- gets the same arithmetic instead of writing it again. `rst` is NET
    downward at the top and `rsut` is the upward part, so their sum is the
    incident flux and no separate insolation figure is needed.
    """
    upward = abs(float(upward_w_m2))
    return upward / (float(net_down_w_m2) + upward)


def planetary_albedo(cfg: dict | None = None,
                     climatology: Path | None = None) -> tuple[float, dict]:
    """Measured, not assumed: reflected over incident at the top of atmosphere.

    `rst` is net downward shortwave at the top and `rsut` is the upward part, so
    their difference is the incident flux and no separate insolation figure is
    needed. Taking it from the climatology rather than declaring it means it
    moves when the surface does, which is the whole reason the budget's 0.266
    went stale.
    """
    from netCDF4 import Dataset

    cfg = config(cfg)
    # Through the one resolver, which RAISES when no baseline is named. Doing
    # it here as PROJECT_ROOT / str(cfg[...]) turned a null into the literal
    # path "None" and reported it as a missing file.
    path = climatology or climatology_path(root=PROJECT_ROOT)
    with Dataset(path) as ds:
        w = _gaussian_weights(len(ds.dimensions["lat"]),
                              len(ds.dimensions["lon"]))
        centres = np.asarray(ds.variables["time"][:], dtype=float)

        def mean(name: str) -> float:
            v = np.asarray(ds.variables[name][:], dtype=float)
            if v.ndim == 3:
                # Bins hold unequal numbers of raw records; weight by them
                # rather than equally. climatology.py, CLIM-13.
                v = clim.annual_mean(v, centres)
            return float((v * w).sum() / w.sum())
        net_down = mean("rst")
        upward = abs(mean("rsut"))
        run_id = getattr(ds, "vesper_run_id", None)
        flux = float(getattr(ds, "vesper_flux_ratio", float("nan")))
    incident = net_down + upward
    alpha = planetary_albedo_from_fluxes(net_down, upward)
    return alpha, {
        "source": rel(path),
        "run_id": run_id,
        "flux_ratio": flux,
        "incident_w_m2": round(incident, 3),
        "reflected_w_m2": round(upward, 3),
        "absorbed_w_m2": round(net_down, 3),
    }


def absorbed_w_m2_per_flux_ratio(alpha: float,
                                 cfg: dict | None = None) -> float:
    """The denominator. Absorbed shortwave per unit stellar flux ratio."""
    return incident_w_m2_per_flux_ratio(cfg) * (1.0 - alpha)


def kelvin_per_w_m2(alpha: float, cfg: dict | None = None,
                    slope: float | None = None) -> float:
    """Kelvin per W/m2 of top-of-atmosphere ABSORBED-flux forcing."""
    return (SLOPE_K_PER_FLUX_RATIO if slope is None else slope) \
        / absorbed_w_m2_per_flux_ratio(alpha, cfg)


def forcing_to_kelvin(w_m2: float, alpha: float, cfg: dict | None = None,
                      slope: float | None = None) -> float:
    """A top-of-atmosphere forcing in W/m2 to kelvin. Positive warms."""
    return w_m2 * kelvin_per_w_m2(alpha, cfg, slope)


def radiative_damping_w_m2_per_k(alpha: float, cfg: dict | None = None,
                                 slope: float | None = None) -> float:
    """W/m2 of top-of-atmosphere forcing per kelvin of response.

    The reciprocal of `kelvin_per_w_m2`, and therefore the SAME measurement seen
    from the other side rather than a second one. It is what an energy-balance
    relaxation time divides a heat capacity by: `tau = C / lambda`.

    It is here, and not carried as its own constant beside a slab model, because
    a private copy is a second sensitivity in a project that has already paid
    for having three. `assess_convergence.py` held 1.31 W/m2/K "measured, not
    assumed" with no pointer to the measurement, 11% above what this module's
    slope implies.
    """
    return 1.0 / kelvin_per_w_m2(alpha, cfg, slope)


def flux_ratio_to_kelvin(d_flux: float, slope: float | None = None) -> float:
    return d_flux * (SLOPE_K_PER_FLUX_RATIO if slope is None else slope)


def planetary_albedo_to_kelvin(d_alpha: float, alpha: float,
                               cfg: dict | None = None,
                               slope: float | None = None) -> float:
    """A change in PLANETARY albedo to kelvin. Brighter cools, so the sign flips."""
    return flux_ratio_to_kelvin(-d_alpha / (1.0 - alpha), slope)


# --------------------------------------------------------------------------
# Checks that can fail
# --------------------------------------------------------------------------

def test_identity(alpha: float = 0.3, cfg: dict | None = None) -> None:
    """The albedo route and the forcing route are one conversion, so they agree.

    Not a physical claim and not evidence about the slope: it is a definition,
    which is the only kind of thing worth asserting about arithmetic. It catches
    a sign flip or a dropped `1 - alpha` in either path, which is what it is for.
    The albedo is an argument with an arbitrary default because the identity holds
    for every value of it; using the measured one would suggest otherwise.
    """
    d_alpha = 0.01
    # The forcing that a planetary albedo change of d_alpha exerts on one unit of
    # flux ratio, which is the normalisation `slope` is stated in.
    forcing = -d_alpha * incident_w_m2_per_flux_ratio(cfg)
    a = planetary_albedo_to_kelvin(d_alpha, alpha, cfg)
    b = forcing_to_kelvin(forcing, alpha, cfg)
    if not math.isclose(a, b, rel_tol=1e-12):
        raise AssertionError(f"albedo route {a} != forcing route {b}")


def _live_entries() -> dict[str, dict]:
    """The runs that still EXIST, keyed by id. Empty if there is no index."""
    if not RUN_INDEX.is_file():
        return {}
    runs = json.loads(RUN_INDEX.read_text(encoding="utf-8"))["runs"]
    return {e.get("run_id"): e for e in runs}


def _archived_entries() -> dict[str, dict]:
    """The runs that survive only as IDENTITY, keyed by id.

    Kept separate from the live index rather than merged into it. Merging them is
    what let a bracket measured on two deleted runs go on reproducing its own
    constant: an archived entry is a record of what a run WAS, and reading a
    measurement out of one says nothing about the world the tree describes now.
    """
    entries = {}
    if ARCHIVED_RUNS.is_dir():
        for p in sorted(ARCHIVED_RUNS.glob("*/INDEX_ENTRY.json")):
            e = json.loads(p.read_text(encoding="utf-8"))
            entries.setdefault(e.get("run_id"), e)
    return entries


def _report_metrics(name: str) -> dict | None:
    """The metrics block of one named convergence report, or None if absent."""
    path = CONVERGENCE_REPORTS / name
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def active_build(cfg: dict | None = None) -> str:
    """The build `config/planet.yaml` names, refusing one the registry refuses."""
    import orogen  # local: the registry is only needed for the currency check

    name = str(config(cfg)["source_build"])
    refusal = orogen.activation_refusal(name=name)
    if refusal:
        raise RuntimeError(
            f"config names source_build {name!r}, which lib/orogen.py refuses: "
            f"{refusal}")
    return name


def _verify_bracket(bracket: dict, slope: float, spread: tuple[float, float],
                    tolerance: float, geography: str,
                    cfg: dict | None = None) -> list[str]:
    """The whole check, over any declaration. `verify()` passes this file's.

    Taking the declaration as an argument is what makes the currency checks
    testable: `test_currency_refuses_a_superseded_measurement()` hands it the
    superseded one and asserts the refusals come back.
    """
    live, archived = _live_entries(), _archived_entries()
    try:
        build = active_build(cfg)
    except Exception as exc:                       # noqa: BLE001 - reported, not raised
        return [f"the active build could not be resolved: {exc}"]

    problems, points, half_widths, regimes = [], {}, {}, {}
    for role, want in bracket.items():
        run_id = want["run_id"]
        entry = live.get(run_id)
        if entry is None:
            if run_id in archived:
                was = archived[run_id].get("source_build")
                problems.append(
                    f"{role} run {run_id} survives only as archived identity, on "
                    f"build {was}; the measurement it carries is superseded and "
                    f"cannot be recomputed")
            else:
                problems.append(f"{role} run {run_id} is in no run index at all")
            continue
        on = entry.get("source_build")
        if on != build:
            problems.append(f"{run_id} was run on build {on}, not the active "
                            f"{build}")
        flux = float(entry.get("physical", {}).get("flux_ratio", float("nan")))
        if not math.isclose(flux, want["flux_ratio"], abs_tol=1e-9):
            problems.append(f"{run_id} is at flux {flux}, "
                            f"not the {want['flux_ratio']} recorded here")
        geo = entry.get("physical", {}).get("geography")
        if geo != geography:
            problems.append(f"{run_id} is on geography {geo}, "
                            f"not the {geography} the bracket was measured on")

        report = _report_metrics(want["report"])
        if report is None:
            problems.append(f"{run_id} has no convergence report at "
                            f"{rel(CONVERGENCE_REPORTS / want['report'])}")
            continue
        if not report.get("sufficiently_equilibrated_for_worldbuilding"):
            problems.append(f"{run_id} no longer passes every convergence "
                            f"criterion: {report.get('failed_criteria')}")
        # THE WINDOW IS PART OF THE MEASUREMENT. A run has more than one
        # admissible window and the asymptote differs between them, so a report
        # that has been re-taken over a different one is a reading of something
        # else under the same filename. That is what happened here once: see the
        # comment on SLOPE_BRACKET_RUNS.
        want_window = want.get("window")
        if want_window is not None:
            got = (report.get("window_start_year_index"),
                   report.get("window_end_year_index"))
            if tuple(got) != tuple(want_window):
                problems.append(
                    f"{run_id} was measured over orbits "
                    f"{want_window[0]}-{want_window[1]} and its report now "
                    f"covers {got[0]}-{got[1]}; the asymptote is a property of "
                    f"the window, so re-take it over the recorded one or "
                    f"re-derive the whole bracket")
            regimes[role] = want.get("io_regime")
        metrics = report.get("metrics") or {}
        asymptote = metrics.get("temperature_asymptote_k")
        if asymptote is None:
            problems.append(f"{run_id} has no fitted asymptote")
            continue
        if abs(asymptote - want["asymptote_k"]) > 0.05:
            problems.append(f"{run_id} asymptote is {asymptote:.4f} K, "
                            f"not the {want['asymptote_k']:.4f} recorded here")
        # The index carries its own copy for a production assessment. Where it
        # does, the two records have to agree: a re-assessment that reached one
        # and not the other is exactly the drift this file is guarding.
        #
        # ONLY WHERE THIS ARM READS THE RUN'S OWN DEFAULT REPORT. An arm that
        # names a differently-scoped assessment -- the warm endpoint's
        # diagnostic one, the cold endpoint's low-I/O block -- is deliberately
        # reading a window the index does not describe, and requiring those to
        # agree would demand that two windows return one asymptote, which is
        # the same error this bracket was just repaired for.
        default_report = f"{run_id}_convergence.json"
        indexed = (entry.get("convergence_metrics") or {}).get(
            "temperature_asymptote_k")
        if want["report"] != default_report:
            indexed = None
        if indexed is not None and abs(indexed - asymptote) > 1e-6:
            problems.append(f"{run_id} reports asymptote {asymptote:.4f} K in "
                            f"its convergence report and {indexed:.4f} K in the "
                            f"run index")
        points[role] = (flux, asymptote)
        half = metrics.get("temperature_asymptote_half_width_k")
        if half is not None:
            half_widths[role] = float(half)
            declared_half = want.get("asymptote_half_width_k")
            if declared_half is not None and abs(half - declared_half) > 0.005:
                problems.append(
                    f"{run_id} asymptote half-width is {half:.4f} K, not the "
                    f"{declared_half:.4f} recorded here")

    # ONE INSTRUMENT ON BOTH ARMS. The two I/O regimes disagree by about 0.17 K
    # on this model, which is 2 per cent of the chord, so a difference taken
    # across them measures the instrument as much as the flux.
    if len(regimes) == 2 and len(set(regimes.values())) != 1:
        problems.append(
            "the two endpoints are declared on different I/O regimes ("
            + ", ".join(f"{role}: {regime}" for role, regime
                        in sorted(regimes.items()))
            + "); a difference across a change of instrument is not a flux "
              "response, so both arms have to be read on one of them")

    if len(points) == 2:
        (fc, tc), (fw, tw) = points["cold"], points["warm"]
        measured = (tw - tc) / (fw - fc)
        if abs(measured - slope) > tolerance:
            problems.append(
                f"the bracket now measures {measured:.1f} K per unit flux ratio "
                f"against the declared {slope}")
        if len(half_widths) == 2:
            floor = math.hypot(*half_widths.values()) / abs(fw - fc)
            lo, hi = spread
            if hi - lo < 2.0 * floor:
                problems.append(
                    f"the declared spread {lo}-{hi} is narrower than the "
                    f"{2 * floor:.2f} K per unit flux ratio the two fits' own "
                    f"uncertainties propagate to")
            if not lo <= measured <= hi:
                problems.append(
                    f"the bracket measures {measured:.1f} K per unit flux ratio, "
                    f"outside the declared spread {lo}-{hi}")
    return problems


def verify(cfg: dict | None = None) -> list[str]:
    """Recompute the declared slope, and check it is CURRENT. Empty means agreed.

    Two questions, and only the first was ever asked. Does the declaration match
    the runs it names -- and are those runs on the world this tree describes now?
    """
    return _verify_bracket(SLOPE_BRACKET_RUNS, SLOPE_K_PER_FLUX_RATIO,
                           SLOPE_SPREAD_K_PER_FLUX_RATIO,
                           SLOPE_TOLERANCE_K_PER_FLUX_RATIO, SLOPE_GEOGRAPHY,
                           cfg)


def test_currency_refuses_a_superseded_measurement(
        cfg: dict | None = None) -> None:
    """Drive the currency checks with the superseded declaration; they must fire.

    This is the test the check exists for, and it can fail: it asserts a named
    refusal rather than a difference. `SUPERSEDED_BRACKET_RUNS` names two runs
    that were deleted and archived, on a build that is not the active one, and
    the arithmetic between their archived asymptotes still reproduces 202.0
    exactly. A `verify()` that reads archived identity therefore passes on them,
    which is what it did. Anything that lets those two runs through again turns
    this red.
    """
    problems = _verify_bracket(SUPERSEDED_BRACKET_RUNS,
                               SUPERSEDED_SLOPE_K_PER_FLUX_RATIO,
                               SUPERSEDED_SPREAD_K_PER_FLUX_RATIO,
                               SLOPE_TOLERANCE_K_PER_FLUX_RATIO,
                               SUPERSEDED_GEOGRAPHY, cfg)
    for want in SUPERSEDED_BRACKET_RUNS.values():
        run_id = want["run_id"]
        if not any(run_id in p for p in problems):
            raise AssertionError(
                f"the currency check passed {run_id}, which is deleted and was "
                f"run on a build that is not active. verify() cannot tell a "
                f"superseded measurement from a current one.")



def test_window_and_regime_refuse_an_unmatched_pair(
        cfg: dict | None = None) -> None:
    """The window and I/O-regime guards must both be able to return no.

    Named refusals rather than differences, and each is driven by the exact
    substitution that once passed silently:

    THE WINDOW. Point the cold endpoint at the run's own default report. That
    file describes the SAME RUN and is a perfectly good assessment, which is
    why the swap went unnoticed -- it just covers the twelve clean-I/O orbits
    rather than the 37-69 the bracket was measured over, and the asymptote is
    0.15 K higher there. Nothing but the window guard distinguishes the two.

    THE REGIME. Relabel one arm and leave everything else alone, so the refusal
    can only come from the declared instruments differing.
    """
    swapped = copy.deepcopy(SLOPE_BRACKET_RUNS)
    swapped["cold"]["report"] = f"{swapped['cold']['run_id']}_convergence.json"
    problems = _verify_bracket(swapped, SLOPE_K_PER_FLUX_RATIO,
                               SLOPE_SPREAD_K_PER_FLUX_RATIO,
                               SLOPE_TOLERANCE_K_PER_FLUX_RATIO,
                               SLOPE_GEOGRAPHY, cfg)
    if not any("covers" in problem for problem in problems):
        raise AssertionError(
            "the cold endpoint's report was swapped for an assessment of a "
            "different window and the bracket accepted it. The asymptote is a "
            "property of the window, so nothing else in verify() can catch "
            f"this: {problems}")

    relabelled = copy.deepcopy(SLOPE_BRACKET_RUNS)
    relabelled["cold"]["io_regime"] = "clean_io"
    problems = _verify_bracket(relabelled, SLOPE_K_PER_FLUX_RATIO,
                               SLOPE_SPREAD_K_PER_FLUX_RATIO,
                               SLOPE_TOLERANCE_K_PER_FLUX_RATIO,
                               SLOPE_GEOGRAPHY, cfg)
    if not any("I/O regimes" in problem for problem in problems):
        raise AssertionError(
            "the two endpoints were declared on different I/O regimes and the "
            f"bracket accepted it: {problems}")

    # The paired positive, so a refusal that fires on everything is not read as
    # this guard working: the declaration as it stands passes both.
    problems = verify(cfg)
    if problems:
        raise AssertionError(
            f"the declared bracket does not pass its own guards: {problems}")

def provenance(cfg: dict | None = None) -> dict:
    """The conversion and where every number in it came from, for an artifact."""
    cfg = config(cfg)
    alpha, alpha_src = planetary_albedo(cfg)
    absorbed = absorbed_w_m2_per_flux_ratio(alpha, cfg)
    return {
        "module": "lib/sensitivity.py",
        "slope_k_per_flux_ratio": SLOPE_K_PER_FLUX_RATIO,
        "slope_spread_k_per_flux_ratio": list(SLOPE_SPREAD_K_PER_FLUX_RATIO),
        "slope_tolerance_k_per_flux_ratio": SLOPE_TOLERANCE_K_PER_FLUX_RATIO,
        "slope_bracket_runs": SLOPE_BRACKET_RUNS,
        "slope_geography": SLOPE_GEOGRAPHY,
        "slope_source_build": config(cfg)["source_build"],
        "slope_regime": "a forward secant from the baseline flux to 0.055 above "
                        "it, on the low-ice branch and carrying the sea-ice "
                        "retreat over that interval. The baseline is the cold "
                        "endpoint rather than inside the bracket; a converged run "
                        "below it on this build is what localises the slope. The "
                        "0.85-to-0.95 stellar sweep gives 330 K per unit flux "
                        "ratio because it crosses the ice transition proper; that "
                        "value is not interchangeable with this one.",
        "planetary_albedo": round(alpha, 5),
        "planetary_albedo_source": alpha_src,
        "incident_w_m2_per_flux_ratio": round(incident_w_m2_per_flux_ratio(cfg), 3),
        "absorbed_w_m2_per_flux_ratio": round(absorbed, 3),
        "kelvin_per_w_m2_absorbed": round(SLOPE_K_PER_FLUX_RATIO / absorbed, 4),
        "denominator_note": "A top-of-atmosphere forcing is an absorbed-flux "
                            "change, so it divides by absorbed flux per unit "
                            "flux ratio. Dividing by incident insolation "
                            "understates by 1/(1 - planetary_albedo).",
    }


def main() -> None:
    # Parsed before anything reads an artifact, so `--help` is answerable with
    # no climatology on disk. It took no arguments at all before, which meant
    # `--help` ran the whole computation and smoke_test's entry-point check was
    # passing on the data being present rather than on argparse being built.
    argparse.ArgumentParser(
        description="The project's one flux-to-kelvin conversion, with its "
                    "provenance. Takes no arguments. Prints the slope, the "
                    "planetary albedo it is evaluated at, and a verification "
                    "against the run index; needs the baseline climatology "
                    "named in config/planet.yaml.").parse_args()
    cfg = config()
    test_identity()
    test_currency_refuses_a_superseded_measurement(cfg)
    p = provenance(cfg)
    alpha = p["planetary_albedo"]
    print(f"slope          {SLOPE_K_PER_FLUX_RATIO} K per unit flux ratio "
          f"({SLOPE_SPREAD_K_PER_FLUX_RATIO[0]}-{SLOPE_SPREAD_K_PER_FLUX_RATIO[1]})")
    print(f"planet albedo  {alpha:.5f}  from {p['planetary_albedo_source']['source']}")
    print(f"incident       {p['incident_w_m2_per_flux_ratio']} W/m2 per unit flux ratio")
    print(f"absorbed       {p['absorbed_w_m2_per_flux_ratio']} W/m2 per unit flux ratio")
    print(f"conversion     {p['kelvin_per_w_m2_absorbed']} K per W/m2 absorbed")
    problems = verify(cfg)
    print("\nverify: " + ("current, and agrees with the runs it names"
                          if not problems else "; ".join(problems)))


if __name__ == "__main__":
    main()
