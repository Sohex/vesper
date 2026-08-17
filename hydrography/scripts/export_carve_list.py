#!/usr/bin/env python3
"""Write the carve verdict as a World Orogen `--preserve-basins` list.

Three outcomes per basin, and the middle one is the point:

  carved     The basin overflows, and with enough discharge to cut its sill
             through within the relaxation window. Retain 0, which Orogen treats
             as not preserved, so drainage enforcement carves the outlet as it
             would have without preservation at all.
  preserved  The basin does not overflow. Retain 1, rim intact.
  marginal   The basin overflows but only trickles over. Retain between 0 and 1,
             which cuts a notch at the saddle and tapers it over the divide band,
             producing a through-flowing valley with a residual lake.

Retain is what Orogen cuts with, so it answers a geomorphic question: how much of
the rim survives the water crossing it. Two separate things decide that, and they
are computed separately here.

**Whether the basin overflows** is the balance at spill level,

    Q = runoff * (catchment - area_at_spill) - (E - P) * area_at_spill

which is the same test as `(E - P) / runoff <= catchment / area_at_spill - 1`
wherever runoff is positive, and is defined where it is not. A basin with no
catchment at all still gains water when its own lake surface receives more
precipitation than it evaporates, and it then fills until it spills. Dividing by
runoff loses that case, which is how 228 dry pans came to be carved on one pass
and pinned shut on the next.

**How deeply the outlet is cut** is Q again, through stream power. Incision goes
as `K Q^m S^n`, so a basin that pours cuts its sill orders of magnitude faster
than one that trickles over it, and across the 1e4 to 1e6 years a landscape takes
to relax that is the difference between a drained depression and a lake with a
notch in its rim.

    retain_incision = 1 - (Q / Q_full) ** INCISION_EXPONENT, clipped to [0, 1]

`Q_full` is per basin, being the declared constant divided by the erodibility of
the rock at that basin's own sill: soft rock is cut through by less water. That
factor is the export's `erodibility`, which is already a relative stream-power
multiplier normalised to 1 over land, and it is the contrast a channel network
expresses rather than the contrast between intact rock samples. The two differ by
orders of magnitude and the expressed one is what a landscape model wants; see
`sill_erodibility` below.

No absolute time-to-cut is attempted even so. What this mapping carries is the
ordering, which the previous one did not have at all: a trickle and a torrent
both went out at retain 0.

**The uncertainty is a third thing, it is kept apart, and it turns out not to
decide anything.** Retain is the larger of the incision value and the margin
below, on the reasoning that either is a reason to leave a rim standing. But
taking the larger lets the margin raise a retain and never lower one, and the
margin is positive exactly when the basin does not overflow, which is exactly
when the incision term is already 1. Both read the sign of the same `Q`, so
`retain` equals `retain_incision` on every basin -- measured, 0 of 3,621 where it
does not. The margin is still computed and still written to the sidecar, because
how close a preserved basin sits to its threshold is worth reading; it is not a
second input to the cut, and a marginal count is a discharge statement.

A basin balances exactly at an evaporation of `E* = P + critical * runoff`. Under
the Penman estimate it evaporates `E_penman`. The fractional margin

    margin = (E_penman - E*) / E_penman
    retain_margin = min(1, margin / TOLERANCE)

says how much open-water evaporation would have to fall before the basin starts
overflowing. A basin needing a 3% change is genuinely marginal and keeps little
of its rim; one needing 40% is comfortably closed and keeps all of it.

TOLERANCE is 0.25, set from the uncertainty that actually dominates: the biosphere
is assumed rather than modelled, and that assumption is worth 3.7 to 7.1 K, which
moves evaporation by considerably more than the 1.7% Penman itself was validated
to.

An earlier mapping interpolated the critical index between the two evaporation
estimates. It is retained in the sidecar as `retain_span` but is not used, because
it measures the wrong thing: the land-evaporation end is a deliberately weak lower
bound, so the span is enormous and the critical value lands near its low end
almost regardless of the basin. It put 70% of marginal basins above retain 0.9,
which is to say it called them marginal and then declined to cut them.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
import json
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
import yaml

import carve_verdict as cv
from _paths import CONFIG, DATA, PROJECT_ROOT  # noqa: F401
from builds import component_data
from orbit import orbital_year_days
from paths import climatology_path
from orogen import Export
from lake_balance import BasinSet


# Declared before any verdict was taken with them, in the same way the albedo
# bracket's criteria were.
#
# Q_FULL_M3_S is the overflow above which the sill is treated as cut through
# within the relaxation window. 1 m3/s is the perennial-stream boundary: above
# it a channel flows year round and works on its bed every year, below it flow
# is seasonal to ephemeral and the same volume arrives as pulses that spend most
# of the window not cutting anything. It is a declared scale, not a fitted one,
# and the sidecar records it so a verdict can be re-read against a different
# choice. An order of magnitude either way moves retain by a factor of about
# three at fixed discharge, so report it whenever a marginal count is quoted.
#
# INCISION_EXPONENT is the discharge exponent m in `K Q^m S^n`, 0.5 being the
# middle of the 0.4 to 0.6 range detachment-limited bedrock studies fit. The
# mapping is far less sensitive to it than to Q_FULL.
Q_FULL_M3_S = 1.0
INCISION_EXPONENT = 0.5
KM3_PER_YEAR_TO_M3_PER_S = 1e9


def incision_retain(q_km3_per_year, year_s: float, q_full_m3_s=Q_FULL_M3_S):
    """Rim surviving the overflow. 1 keeps it, 0 cuts it.

    Takes the overflow at spill level in km3 per Vesper year. A basin that does
    not overflow gets 1 by construction, since Q is then zero or negative.

    `q_full_m3_s` is per basin once the sill's own rock is known, since a soft
    sill is cut through by less water than a hard one.
    """
    q = np.asarray(q_km3_per_year, dtype=float) * KM3_PER_YEAR_TO_M3_PER_S / year_s
    cut = np.power(np.clip(q, 0.0, None) / np.asarray(q_full_m3_s, dtype=float),
                   INCISION_EXPONENT)
    return np.clip(1.0 - cut, 0.0, 1.0)


def sill_erodibility(basins_path: Path, terrain_hash: str) -> np.ndarray:
    """The rock at each basin's outlet, as a relative stream-power multiplier.

    The export's `erodibility` is exactly this quantity and says so: a relative
    stream-power multiplier from the exposed rock, mean-normalised to 1 over
    land. So it enters as `K` does, and Q_full divides by it: a soft sill is cut
    through by less water than a hard one.

    **It is the expressed contrast, not the intact-rock one, and that is the
    number a landscape model wants.** Stock and Montgomery (1999) measure K
    across five orders of magnitude between lithologies, but Zondervan (2020)
    measures the contrast a real channel network expresses at about 4x, because
    channels adjust width and slope in response to the rock they are in. This
    field spans under 4x, which is the right order; using the intact-rock spread
    would make the rock the only term that mattered and the discharge decorative.

    The saddle is smaller than a mesh cell, so the two regions either side of it
    bracket the rock being cut rather than naming it. Their geometric mean is
    what is used, that being the right average for a multiplicative factor, and
    the two sides agree only loosely: correlation 0.34 on the current build.
    """
    export = Export()
    if export.terrain_hash != terrain_hash:
        raise SystemExit(
            f"basins.nc was built from terrain {terrain_hash[:16]} and the "
            f"configured export is {export.terrain_hash[:16]}. The sill lookup "
            "is by region index, which does not survive a terrain change."
        )
    try:
        ero = export.field("erodibility")
    except Exception:
        print("note: the export carries no erodibility field; sills are uniform")
        with Dataset(basins_path) as ds:
            return np.ones(ds.dimensions["basin"].size)

    with Dataset(basins_path) as ds:
        inner = np.asarray(ds["spill_region"][:])
        outer = np.asarray(ds["spill_exit_region"][:])
    if inner.size and max(int(inner.max()), int(outer.max())) >= ero.size:
        raise SystemExit("spill regions index past the export's mesh")

    # -1 marks a basin whose saddle was not resolved on one side or either.
    # Fall back through the side that exists to the land mean, which is 1.
    a = np.where(inner >= 0, ero[np.maximum(inner, 0)], np.nan)
    b = np.where(outer >= 0, ero[np.maximum(outer, 0)], np.nan)
    with np.errstate(invalid="ignore"):
        both = np.sqrt(a * b)
    out = np.where(np.isfinite(both), both,
                   np.where(np.isfinite(a), a, np.where(np.isfinite(b), b, 1.0)))
    return np.clip(out, 1e-3, None)


def main() -> None:
    ap = argparse.ArgumentParser()
    # Every per-build path below defaults to None and is resolved AFTER
    # parse_args. Resolving them while building the parser meant --help did
    # real work and died on a missing build, which is also why nothing caught
    # that the climatology default pointed at a superseded directory.
    #
    # The paths are per-build and strict. This script's output leaves the
    # project and changes the terrain, so a default that quietly reads or writes
    # the wrong build is the most expensive one here.
    ap.add_argument("--climatology", type=Path, default=None,
                    help="regular climatology; defaults to config's "
                         "baseline_climatology")
    ap.add_argument("--coupling", type=Path, default=None)
    # Hydrography is per-build now, so the basin set has to be selectable
    # alongside the coupling it was built with. Mixing a coupling matrix from one
    # terrain with a basin catalogue from another would misalign the rows in the
    # same silent way the longitude convention misaligned the columns.
    ap.add_argument("--basins", type=Path, default=None)
    # Which field stands for runoff generated over the catchment. This was
    # `mrro` and should not have been: `mrro` is not local runoff generation but
    # river-routed net divergence, because landmod.f90's roffstep calls mkradv,
    # which advects runoff downhill and modifies its argument in place. So it is
    # local generation minus river outflow plus river inflow, on ExoPlaSim's own
    # grid and its own downhill directions, which know nothing about our basins.
    #
    # P - E is the water balance the derivation actually wants: whatever falls on
    # the catchment and does not evaporate is what reaches the sink. Both of this
    # project's other consumers of catchment runoff already made this call --
    # pedology's `runoff_source: p_minus_e`, with the reasoning in
    # pedogenesis.yaml, and surface_water.py's lake solver. The carve verdict was
    # the only one left on `mrro`, and it is the one whose output changes the
    # terrain.
    ap.add_argument("--runoff-source", choices=("p_minus_e", "mrro"),
                    default="p_minus_e")
    # Iteration 2 onward. A verdict is taken on the basins a build still has,
    # but Orogen regenerates from the planet code and needs the whole catalogue,
    # including the basins an earlier pass already carved. Merging keeps the loop
    # monotone: a basin carved in a previous iteration stays carved, because the
    # terrain that justified re-examining it no longer exists. Without this the
    # list would silently re-preserve every basin the current build has already
    # lost, and the next build would undo the last one.
    ap.add_argument("--previous", type=Path, default=None,
                    help="carve_list.json from the pass that produced the "
                         "current build; its retain-0 basins are carried forward")
    ap.add_argument("--config", type=Path, default=CONFIG)
    ap.add_argument("--out-list", type=Path, default=None)
    ap.add_argument("--out-json", type=Path, default=None)
    args = ap.parse_args()

    _bd = component_data("hydrography", strict=True)
    if args.coupling is None:
        args.coupling = _bd / "coupling_exoplasim-T42.nc"
    if args.basins is None:
        args.basins = _bd / "basins.nc"
    if args.out_list is None:
        args.out_list = _bd / "carve_list.txt"
    if args.out_json is None:
        args.out_json = _bd / "carve_list.json"
    if args.climatology is None:
        args.climatology = climatology_path()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    basins = BasinSet(args.basins)
    resolution = str(config["model"]["resolution"]).upper()

    with Dataset(args.climatology) as ds:
        am = cv.annual_mean
        pr, evap, mrro = am(ds, "pr"), -am(ds, "evap"), am(ds, "mrro")
        ts, tas = am(ds, "ts"), am(ds, "tas")
        ps_pa, rss, rls = am(ds, "ps") * 100.0, am(ds, "rss"), am(ds, "rls")
        field_lon = np.asarray(ds["lon"][:])
    q_air, wind = cv.turbulent_forcing(args.climatology)

    land_albedo = cv.read_sra_field(
        PROJECT_ROOT / "exoplasim" / "inputs" / resolution.lower()
        / f"orogen_{resolution}_surf_0174.sra", *ps_pa.shape)
    penman = np.maximum(cv.penman_open_water(
        ts, tas, q_air, wind, ps_pa, rss, rls, land_albedo,
        float(config["planet"]["gravity_m_s2"])), evap)

    runoff_field = mrro if args.runoff_source == "mrro" else (pr - evap)
    means, _ = cv.basin_means(args.coupling,
                              {"pr": pr, "wet": evap, "pen": penman,
                               "ro": runoff_field},
                              basins.n, field_lon=field_lon)
    # Orbital period varies with flux, so take it from the config rather than
    # hardcoding. The literal here was 189.6145 d, the 0.90-flux year, while this
    # baseline runs at 0.96 and 180.655 d. It cancels out of the aridity index
    # and the equilibrium lake area, both ratios, but it was making the reported
    # runoff depth 5% high.
    year_s = orbital_year_days(config) * 86400.0
    # Clamped at zero for the same reason carve_verdict.py clamps it: a
    # catchment losing more to evaporation than it receives delivers nothing,
    # not a negative amount. It keeps the reported depth physical, and the
    # discharge below reads the same clamped field.
    runoff = np.maximum(means["ro"], 0.0) * year_s / 1000.0
    precip = means["pr"] * year_s / 1000.0
    crit = basins.catchment_km2 / np.maximum(basins.area_at_spill_km2, 1e-9) - 1.0

    # Reported, not decided on. The index is undefined where a basin has no
    # catchment runoff, which is exactly the set the discharge form handles.
    def index(field):
        e = means[field] * year_s / 1000.0
        return np.where(runoff > 0, (e - precip) / np.where(runoff > 0, runoff, 1.0),
                        np.inf)

    idx_wet, idx_pen = index("wet"), index("pen")

    # The verdict, as the water that has to leave at spill level. Runoff is
    # generated over the dry catchment; the lake surface itself gains P and
    # loses E, which is the term the ratio form divides away.
    dry_km2 = np.maximum(basins.catchment_km2 - basins.area_at_spill_km2, 0.0)

    def discharge(field):
        e = means[field] * year_s / 1000.0
        return runoff * dry_km2 - (e - precip) * basins.area_at_spill_km2

    q_wet, q_pen = discharge("wet"), discharge("pen")
    # E_wet is the model's own land evaporation and E_pen is Penman over open
    # water, floored at it, so q_wet >= q_pen everywhere and the sets nest.
    carved = q_pen > 0.0
    preserved = q_wet <= 0.0
    marginal = ~carved & ~preserved
    assert int(carved.sum() + preserved.sum() + marginal.sum()) == basins.n

    TOLERANCE = 0.25
    e_pen = means["pen"] * year_s / 1000.0
    e_balance = precip + crit * runoff              # evaporation that exactly balances
    with np.errstate(divide="ignore", invalid="ignore"):
        margin = np.where(e_pen > 0, (e_pen - e_balance) / np.where(e_pen > 0, e_pen, 1.0), 1.0)
    retain_margin = np.where(carved, 0.0, np.clip(margin / TOLERANCE, 0.0, 1.0))

    # What the water can actually cut, which is the half the margin never knew.
    # The sill's own rock sets how much water that takes.
    sill_ero = sill_erodibility(args.basins, basins.terrain_hash)
    q_full = Q_FULL_M3_S / sill_ero
    retain_incision = incision_retain(q_pen, year_s, q_full)

    # Either is a reason to leave a rim standing: that the basin may not overflow
    # at all, or that its overflow cannot cut. Taking the larger keeps both.
    retain = np.maximum(retain_incision, retain_margin)

    # The superseded mapping, kept for comparison in the sidecar only.
    span = idx_pen - idx_wet
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(np.isfinite(span) & (span > 1e-9),
                     (crit - idx_wet) / np.where(span > 1e-9, span, 1.0), 0.5)
    retain_span = np.where(carved, 0.0, np.where(preserved, 1.0, np.clip(1.0 - f, 0.0, 1.0)))

    # The verdict follows retain rather than the overflow test, because what
    # Orogen does to a basin is set by retain. A basin can overflow and still
    # keep most of its rim, if what crosses the sill is a trickle.
    verdict_name = np.where(retain <= 0.0, "carve",
                            np.where(retain >= 1.0, "preserve", "marginal"))
    n_carve = int((retain <= 0.0).sum())
    n_preserve = int((retain >= 1.0).sum())
    n_marginal = basins.n - n_carve - n_preserve

    with Dataset(args.basins) as ds:
        ids = [str(x) for x in ds["basin_id"][:]]

    carried = {}
    if args.previous is not None:
        prev = json.loads(args.previous.read_text(encoding="utf-8"))
        here = set(ids)
        for entry in prev["basins"]:
            if entry["id"] not in here:
                # Absent from this build because it was carved away. Its verdict
                # is not re-decidable and is carried at 0.
                carried[entry["id"]] = 0.0
            elif float(entry["retain"]) <= 0.0:
                carried[entry["id"]] = 0.0
        if carried and any(float(v) > 0 for v in carried.values()):
            raise RuntimeError("carried-forward entries must all be retain 0")

    pass_label = ("first pass" if not carried
                  else f"pass {len(carried) and 2}, {len(carried)} basins carried forward")
    n_carried = len(carried)
    n_zero = n_carve + n_carried
    n_total = basins.n + n_carried
    n_here = basins.n
    clim_name = args.climatology.name
    flux_earth = float(config["orbit"]["baseline_flux_earth"])
    with Dataset(args.climatology) as ds:
        _ts = np.asarray(ds["ts"][:]).mean(axis=0)
        _lat = np.asarray(ds["lat"][:])
        _w = np.cos(np.deg2rad(_lat))[:, None] * np.ones_like(_ts)
        mean_ts = float((_w * _ts).sum() / _w.sum())
    runoff_source = ("precipitation minus evaporation over the catchment"
                     if args.runoff_source == "p_minus_e"
                     else "the model's mrro field")
    land_surface = str(config["model"].get("land_albedo_source", "uniform"))
    glaciers = ("glaciers enabled" if (config["surface"].get("glaciers") or {}).get("enabled")
                else "glaciers off")
    header = f"""# Vesper carve verdict, {pass_label}
#
# Produced from a converged ExoPlaSim climatology: {resolution},
# {flux_earth:g} S-Earth, {land_surface} land surface, {glaciers},
# {mean_ts:.2f} K, climatology {clim_name}.
# Terrain {basins.terrain_hash[:16]}, basin catalogue unchanged.
#
# A basin overflows when more water arrives than its lake surface can evaporate,
#     Q = runoff * (catchment - area_at_spill) - (E - P) * area_at_spill  >  0
# and retain is then what that Q can cut, as 1 - (Q/Q_full)^{INCISION_EXPONENT:g}, where Q_full is
# {Q_FULL_M3_S:g} m3/s at land-mean rock and less where the sill itself is softer,
# floored by how uncertain the overflow test itself is. Open-water evaporation is
# the Penman combination equation with water's albedo and roughness, validated
# against the model over ocean cells to within 1.7%.
# Catchment runoff is {runoff_source}; see the sidecar.
#
#   retain 1.0   {n_preserve:4d} basins  closed, or spilling too little to cut
#   retain 0<r<1 {n_marginal:4d} basins  a notch: uncertain, or a trickle over the sill
#   retain 0.0   {n_zero:4d} basins  carved
#                            {n_carried:4d} of them carried forward, {n_carve:4d} decided here
#
# The counts above are of the whole {n_total:d}-entry catalogue. This pass could
# only decide the {n_here:d} basins the current build still has; the rest were
# carved by an earlier pass and are held at 0 to keep the loop monotone.
#
# Carved basins are listed explicitly at retain 0 rather than omitted, so this
# file is the complete verdict rather than a subset of it. If your parser would
# rather they were absent, say so and we will drop them.
#
# id                                    retain
"""
    lines = [header]
    for i in np.argsort(-retain):
        lines.append(f"{ids[i]:<38s} {retain[i]:.4f}")
    for bid in sorted(carried):
        lines.append(f"{bid:<38s} {0.0:.4f}")
    args.out_list.parent.mkdir(parents=True, exist_ok=True)
    args.out_list.write_text("\n".join(lines) + "\n", encoding="ascii")

    sidecar = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "terrain_hash": basins.terrain_hash,
        "climatology": str(args.climatology),
        "climate": {
            "resolution": resolution,
            "flux_earth": flux_earth,
            "stellar_spectrum": str(config["radiation"].get("stellar_spectrum")),
            "land_surface": land_surface,
            "mean_surface_temperature_k": round(mean_ts, 2),
            "note": None if land_surface == "modelled" else
                    "The land surface is assumed rather than modelled, and the "
                    "albedo bracket puts that assumption at 3.7 to 7.1 K.",
        },
        "method": {
            "test": "Q = runoff*(catchment - area_at_spill) - (E - P)*area_at_spill > 0",
            "test_note": "equivalent to (E - P) / runoff <= catchment / "
                         "area_at_spill - 1 wherever catchment runoff is "
                         "positive, and defined where it is not: a basin with no "
                         "catchment still fills if its own lake surface gains "
                         "more precipitation than it evaporates",
            "open_water_evaporation": "Penman combination, water albedo and "
                                      "roughness, floored at the model's land rate",
            "penman_ocean_validation_ratio": 1.017,
            "retain_mapping": "max(incision, margin), the larger of what the "
                              "overflow cannot cut and what the overflow test "
                              "cannot decide",
            "retain_incision_mapping": f"1 - (Q / Q_full)**{INCISION_EXPONENT:g}, clipped to [0, 1]",
            "retain_incision_q_full_m3_per_s": Q_FULL_M3_S,
            "retain_incision_q_full_note": "per basin, divided by the sill's own "
                "erodibility; the constant above is the value at land-mean rock",
            "retain_incision_exponent": INCISION_EXPONENT,
            "retain_incision_rationale": "stream power goes as K Q^m S^n, so "
                "overflow discharge and not distance from the threshold is what "
                "cuts a sill. Q_full is the perennial-stream scale, declared "
                "rather than fitted; no absolute time-to-cut is attempted "
                "because Stock and Montgomery (1999) measure K across five "
                "orders of magnitude by lithology",
            "sill_rock": "K from the export's erodibility field, a relative "
                "stream-power multiplier mean-normalised to 1 over land, taken "
                "as the geometric mean of the two regions either side of the "
                "saddle. It is the fluvially expressed contrast rather than the "
                "intact-rock one, which is what a landscape model wants; see "
                "Zondervan (2020) in references/INDEX.md",
            "sill_erodibility_range": [round(float(sill_ero.min()), 4),
                                       round(float(sill_ero.max()), 4)],
            "retain_margin_mapping": "min(1, ((E_penman - (P + critical*runoff)) / E_penman) / 0.25)",
            "retain_tolerance": TOLERANCE,
            "retain_tolerance_rationale": "fractional change in open-water "
                "evaporation that would flip the verdict; 0.25 is set by the "
                "assumed-biosphere uncertainty, which dominates Penman's own 1.7%",
        },
        "counts": {"carve": n_carve, "preserve": n_preserve, "marginal": n_marginal},
        "counts_by_penman_test_alone": {
            "carve": int(carved.sum()), "preserve": int(preserved.sum()),
            "disagreeing_with_land_evaporation": int(marginal.sum())},
        "basins": [
            {
                "id": ids[i],
                "verdict": str(verdict_name[i]),
                "retain": round(float(retain[i]), 4),
                "retain_incision": round(float(retain_incision[i]), 4),
                "retain_margin": round(float(retain_margin[i]), 4),
                "retain_span_superseded": round(float(retain_span[i]), 4),
                "overflow_km3_per_year": round(float(q_pen[i]), 6),
                "overflow_m3_per_s": round(
                    float(q_pen[i]) * KM3_PER_YEAR_TO_M3_PER_S / year_s, 6),
                "sill_erodibility": round(float(sill_ero[i]), 4),
                "evaporation_margin": None if not np.isfinite(margin[i])
                                      else round(float(margin[i]), 4),
                "critical_aridity_index": round(float(crit[i]), 4),
                "aridity_index_penman": None if not np.isfinite(idx_pen[i])
                                        else round(float(idx_pen[i]), 4),
                "aridity_index_land_evaporation": None if not np.isfinite(idx_wet[i])
                                                  else round(float(idx_wet[i]), 4),
                "catchment_km2": round(float(basins.catchment_km2[i]), 2),
                "area_at_spill_km2": round(float(basins.area_at_spill_km2[i]), 2),
                "runoff_km_per_year": round(float(runoff[i]), 8),
                "no_catchment_runoff": bool(runoff[i] <= 0),
            } for i in range(basins.n)
        ] + [
            {"id": bid, "verdict": "carve", "retain": 0.0,
             "carried_from_previous_pass": True}
            for bid in sorted(carried)
        ],
    }
    args.out_json.write_text(json.dumps(sidecar, indent=2) + "\n", encoding="utf-8")

    mid = (retain > 0.0) & (retain < 1.0)
    n_trickle = int((carved & (retain > 0.0)).sum())
    n_lake_fed = int((carved & (runoff <= 0)).sum())
    if carried:
        print(f"carried   {len(carried):5d}  retain 0.0, from the previous pass")
    print(f"carve     {n_carve:5d}  retain 0.0  (this pass, of {basins.n} remaining)")
    print(f"          {n_trickle:5d}  overflow, but too little to cut through")
    print(f"          {n_lake_fed:5d}  overflow fed by the lake surface, no catchment runoff")
    if n_marginal:
        print(f"marginal  {n_marginal:5d}  retain {retain[mid].min():.3f}"
              f"-{retain[mid].max():.3f}, median {np.median(retain[mid]):.3f}")
    else:
        print(f"marginal  {n_marginal:5d}")
    print(f"preserve  {n_preserve:5d}  retain 1.0")
    print(f"\nwrote {args.out_list}\n      {args.out_json}")


if __name__ == "__main__":
    main()
