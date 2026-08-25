#!/usr/bin/env python3
"""Solve endorheic lake levels from a climate forcing.

The machinery here is climate-independent; the forcing is not. Nothing in this
file decides what the climate is. It takes runoff and evaporation as arguments
and returns the lake each basin settles at, which is the step that has to wait
for ExoPlaSim.

Run directly to exercise the solver against a uniform placeholder forcing. That
is a smoke test of the machinery and a sensitivity sweep, not a result about
this world.

The balance
-----------

For a basin with catchment area C and lake area A, at equilibrium the water
arriving equals the water leaving::

    r (C - A)  +  p A  =  e A

where r is runoff depth generated per unit of dry catchment, p is precipitation
falling directly on the lake, and e is evaporation from open water, all as
depths per unit time. Solving for area::

    A = r C / (e - p + r)

so the lake is set by the ratio of catchment supply to the net evaporative
demand over water. If ``e - p + r <= 0`` the lake can never evaporate what it
receives and the basin fills to its spill and overflows regardless of size.

Overflow cascades. Most basins spill into another basin rather than to the
ocean (hydrography_report.json carries the count), so a full basin passes its
surplus downstream and the system has to be iterated to a fixed point.
"""

from __future__ import annotations

import argparse
import json

from netCDF4 import Dataset
import numpy as np

from _paths import ANALYSIS, DATA  # noqa: F401
from builds import component_data, grid_export
from orogen import Export

TERMINAL_OCEAN = -1


class BasinSet:
    """Basin hypsometry and topology, as built by `build_hydrography.py`."""

    def __init__(self, path=None):
        # Per-build, strict. Defaulting to the flat data/ paired one terrain's
        # basins with another's everything else; the flat directory currently
        # holds 2,107 basins against the active build's 3,629.
        path = path or (component_data("hydrography", strict=True) / "basins.nc")
        with Dataset(path) as ds:
            self.level_km = np.asarray(ds["level_km"][:])
            self.area_km2 = np.asarray(ds["flooded_area_km2"][:])
            self.volume_km3 = np.asarray(ds["volume_km3"][:])
            self.spill_km = np.asarray(ds["spill_km"][:])
            self.sink_elevation_km = np.asarray(ds["sink_elevation_km"][:])
            self.spill_target = np.asarray(ds["spill_target"][:])
            self.catchment_km2 = np.asarray(ds["catchment_km2"][:])
            self.capacity_km3 = np.asarray(ds["capacity_km3"][:])
            self.area_at_spill_km2 = np.asarray(ds["area_at_spill_km2"][:])
            self.terrain_hash = ds.terrain_hash
        self.n = self.level_km.shape[0]

    @property
    def has_impoundment(self):
        """Does the finished terrain still hold water at spill?

        PRESERVED IS NOT THE SAME AS STILL CLOSED. The catalogue's preserved
        set is the set the drainage conditioning was told not to breach, and
        the conditioning honours that for its carve passes only: erosion,
        ridge sharpening and soil creep run over protected cells by design and
        the sink is left unprotected so the basin floor can erode. A basin can
        therefore come out of a generation published as preserved, with
        `retain: 1` and a hypsometry curve measured on a surface that no longer
        exists, and impound nothing at all.

        Every quantity below that needs an impoundment -- a lake area, a level,
        a capacity, a depth to cut -- is undefined on such a basin. Test this
        rather than membership of the preserved set.
        """
        return self.capacity_km3 > 0.0

    @property
    def depth_at_spill_m(self):
        """Spill point down to the basin floor: what a sill has to be cut by.

        The depth of the lake when the basin is full, and so the length the
        overflow has to remove before the depression is gone.

        **Floored at 1 m, and the floor is a guard rather than a rounding
        correction.** The basins it catches are the ones the drainage
        conditioning flattened: their published `depthKm` runs up to a
        kilometre and a half on the pre-conditioning surface and their
        `capacity_km3` here is exactly zero. A 1 m floor makes such a basin
        carve at any coefficient, which is the right instruction and the wrong
        reason, so `has_impoundment` names the state and callers report it
        rather than letting the floor decide silently.

        notes/audits/basin-catalogue-floor.md carries the measurement.
        """
        return np.maximum((self.spill_km - self.sink_elevation_km) * 1000.0, 1.0)

    def level_for_area(self, basin: int, area_km2: float) -> float:
        """Invert the hypsometric curve. Area is monotone in level."""
        a, l = self.area_km2[basin], self.level_km[basin]
        return float(np.interp(area_km2, a, l))

    def volume_for_area(self, basin: int, area_km2: float) -> float:
        a, v = self.area_km2[basin], self.volume_km3[basin]
        return float(np.interp(area_km2, a, v))


def solve(
    basins: BasinSet,
    runoff_km_per_year: np.ndarray,
    lake_evaporation_km_per_year: np.ndarray,
    lake_precip_km_per_year: np.ndarray,
    *,
    max_iterations: int = 200,
    tolerance_km3: float = 1e-3,
):
    """Equilibrium lake area, level and volume for every basin.

    All three forcings are per-basin depths per year: `runoff_km_per_year` is
    generated over the dry catchment, the other two apply over open water.

    Returns a dict of arrays plus `converged`, which is False if the overflow
    cascade did not reach a fixed point. Do not quietly use a non-converged
    result; it means basins are exchanging overflow in a cycle.
    """
    n = basins.n
    supply = runoff_km_per_year * basins.catchment_km2      # km3/yr from land
    demand = lake_evaporation_km_per_year - lake_precip_km_per_year + runoff_km_per_year

    inflow_extra = np.zeros(n)      # overflow received from upstream basins
    area = np.zeros(n)
    overflow = np.zeros(n)
    converged = False

    for _ in range(max_iterations):
        total_supply = supply + inflow_extra
        with np.errstate(divide="ignore", invalid="ignore"):
            want = np.where(demand > 0, total_supply / demand, np.inf)
        area = np.minimum(np.nan_to_num(want, posinf=np.inf), basins.area_at_spill_km2)

        # A basin pinned at its spill passes on whatever it cannot evaporate.
        at_spill = area >= basins.area_at_spill_km2 - 1e-9
        overflow_new = np.where(
            at_spill,
            np.maximum(total_supply - basins.area_at_spill_km2 * demand, 0.0),
            0.0,
        )

        received = np.zeros(n)
        downstream = basins.spill_target
        routed = (downstream >= 0) & (overflow_new > 0)
        np.add.at(received, downstream[routed], overflow_new[routed])

        if np.max(np.abs(received - inflow_extra)) < tolerance_km3:
            inflow_extra = received
            overflow = overflow_new
            converged = True
            break
        inflow_extra = received
        overflow = overflow_new

    level = np.array([basins.level_for_area(b, area[b]) for b in range(n)])
    volume = np.array([basins.volume_for_area(b, area[b]) for b in range(n)])
    return {
        "area_km2": area,
        "level_km": level,
        "volume_km3": volume,
        "overflow_km3_per_year": overflow,
        "fills_to_spill": area >= basins.area_at_spill_km2 - 1e-9,
        "dry": area <= 0,
        "converged": converged,
    }


# ---------------------------------------------------------------------------
# The periodic steady state, and the criteria it is judged on
# ---------------------------------------------------------------------------
# DECLARED BEFORE ANY OF THIS WAS RUN, because the answer here is a size and not
# a yes. Each is a threshold on a quantity whose scale is set by something other
# than the result: the closure by the basin's own capacity, the publication
# threshold by the mesh the areas came off.

# The year must close on itself in the basin's own volume. Relative to capacity
# because capacities span orders of magnitude here, with an absolute floor so a
# small basin is not held to a bar below its own arithmetic.
PERIODIC_CLOSURE_RELATIVE = 1e-3
PERIODIC_CLOSURE_FLOOR_KM3 = 1e-6
# Years of the cycle to iterate before a basin that has not closed is refused.
PERIODIC_MAX_YEARS = 400
# No sub-step may move a basin by more than this share of its capacity. Sets the
# integration error, not the answer.
PERIODIC_MAX_STEP_FRACTION = 0.02
# Share of the explicit scheme's own stability limit a sub-step may take. Below
# 2 by the stability bound; 0.5 leaves margin for the slope changing inside the
# step, which it does at every break in a hypsometric curve.
PERIODIC_STABILITY_MARGIN = 0.5

# A basin HAS a seasonal cycle worth publishing when both hold. The fraction is
# a statement about the basin; the area floor is a statement about the
# instrument, and it is the caller's mean mesh cell area, because a lake area
# comes off a curve built from mesh cells and a swing smaller than one cell is
# below what the curve can express.
SEASONAL_AMPLITUDE_MIN_FRACTION = 0.10


def solve_periodic(
    basins: BasinSet,
    runoff_km_per_year: np.ndarray,
    lake_evaporation_km_per_year: np.ndarray,
    lake_precip_km_per_year: np.ndarray,
    bin_lengths_years: np.ndarray,
    *,
    area_resolution_km2: float,
    initial_volume_km3: np.ndarray | None = None,
    max_years: int = PERIODIC_MAX_YEARS,
    closure_relative: float = PERIODIC_CLOSURE_RELATIVE,
    closure_floor_km3: float = PERIODIC_CLOSURE_FLOOR_KM3,
):
    """Lake storage through the year, as a PERIODIC steady state.

    `solve()` above answers a different question and answers it exactly: where
    does a basin settle if the forcing never changes. This one integrates the
    same balance through the climatology's own time bins against the basin's own
    hypsometric curve, and asks for the year to close on itself rather than for
    the storage to stop moving. **The convergence criterion is the new thing
    here, not the solver.**

    The three forcings are `(nbin, n_basins)` in km per year, the same units and
    the same convention as `solve()`: runoff is generated over the dry catchment
    and the other two act on open water.

    `bin_lengths_years` is `(nbin,)` and ABSOLUTE, IN THE SAME YEAR THE FLUXES
    ARE PER, and must sum to one cycle of the forcing. It is absolute rather
    than normalised because a storage integration is the one calculation here
    that has a time in it: `solve()` finds a fixed point, which does not care
    how long a year is, but how far a lake gets through its cycle does.
    `surface_water.py` supplies fluxes per VESPER year against a climatology
    binned over one Vesper year, so its bin lengths are the normalised weights
    `lib/climatology.py` recovers. A caller mixing two year units gets an
    amplitude wrong by their ratio and a closure that still passes, which is why
    `_selftest` checks that the cycle length reaches the answer at all.

    `area_resolution_km2` is the mean mesh cell area the hypsometric curves were
    built from. It sets the publication threshold and nothing else.

    WHAT IS REFUSED RATHER THAN REPORTED. A basin with no impoundment on the
    finished terrain has no state to integrate and is refused outright; see
    `has_impoundment`. A basin whose year does not close within `max_years` is
    refused too, and `closed` says which. A caller that reads the arrays without
    reading `closed` is reading a transient.

    WHAT IS PUBLISHED. Per basin and per bin, the area, level and volume; and
    per basin, the peak-to-trough area range, that range as a fraction of the
    annual mean area, the residence time, and `seasonal` -- whether both
    published-cycle criteria hold. The honest product is a per-basin amplitude
    carried by the low-capacity basins, not a claim that every preserved basin
    has a season: a deep terminal lake holds years of supply and its surface
    barely moves inside one, while a shallow playa's area is almost all
    seasonal.
    """
    r = np.atleast_2d(np.asarray(runoff_km_per_year, dtype=float))
    e = np.atleast_2d(np.asarray(lake_evaporation_km_per_year, dtype=float))
    pcp = np.atleast_2d(np.asarray(lake_precip_km_per_year, dtype=float))
    dt_bin = np.asarray(bin_lengths_years, dtype=float)
    nbin, n = r.shape
    if e.shape != r.shape or pcp.shape != r.shape:
        raise ValueError(f"forcing shapes disagree: {r.shape}, {e.shape}, {pcp.shape}")
    if dt_bin.shape != (nbin,):
        raise ValueError(f"bin_lengths_years is {dt_bin.shape}, expected ({nbin},)")
    if n != basins.n:
        raise ValueError(f"forcing carries {n} basins against the set's {basins.n}")
    if not np.all(dt_bin > 0):
        raise ValueError("every bin length must be positive")

    capacity = basins.capacity_km3
    live = basins.has_impoundment
    cap_safe = np.where(live, capacity, 1.0)
    catchment = basins.catchment_km2

    def area_and_slope(vol, rows=None):
        """Invert each basin's volume-to-area curve, and give its local slope.

        Both curves are monotone. Vectorised over basins rather than looped:
        this is evaluated once per sub-step per bin per cycle, so a Python loop
        over the catalogue here is the whole cost of the routine. `rows`
        restricts it to a subset, which is what the sub-stepping uses once most
        basins have finished their bin.

        `dA/dV` comes back because the integrator needs it. It is what makes the
        balance stiff: a basin whose floor is flat gains a hundred mesh cells of
        surface over one level of the curve, and an explicit step taken without
        reference to that slope oscillates across it instead of settling on it.
        """
        vcurve = basins.volume_km3 if rows is None else basins.volume_km3[rows]
        acurve = basins.area_km2 if rows is None else basins.area_km2[rows]
        m = vcurve.shape[1]
        idx = np.clip((vcurve < vol[:, None]).sum(1) - 1, 0, m - 2)
        r_ = np.arange(vcurve.shape[0])
        v0, v1 = vcurve[r_, idx], vcurve[r_, idx + 1]
        a0, a1 = acurve[r_, idx], acurve[r_, idx + 1]
        span = v1 - v0
        good = span > 0
        frac = np.where(good, (vol - v0) / np.where(good, span, 1.0), 0.0)
        out = a0 + np.clip(frac, 0.0, 1.0) * (a1 - a0)
        slope = np.where(good, (a1 - a0) / np.where(good, span, 1.0), 0.0)
        keep = live if rows is None else live[rows]
        return np.where(keep, out, 0.0), np.where(keep, slope, 0.0)

    def area_of(vol, rows=None):
        return area_and_slope(vol, rows)[0]

    if initial_volume_km3 is None:
        # Seeded from the annual-mean equilibrium, which is the fixed point a
        # basin with no season has. A deep basin then starts AT its answer
        # instead of filling towards it over its residence time, which is what
        # would otherwise decide `max_years`.
        w = dt_bin / dt_bin.sum()
        eq = solve(basins, (w[:, None] * r).sum(0), (w[:, None] * e).sum(0),
                   (w[:, None] * pcp).sum(0))
        vol0 = np.clip(eq["volume_km3"], 0.0, capacity)
    else:
        vol0 = np.clip(np.asarray(initial_volume_km3, dtype=float), 0.0, capacity)
    vol0 = np.where(live, vol0, 0.0)

    inflow = np.zeros((nbin, n))        # km3/yr received from upstream, per bin
    downstream = basins.spill_target
    routes = downstream >= 0
    tol = np.maximum(closure_relative * capacity, closure_floor_km3)

    area_t = np.zeros((nbin, n))
    vol_t = np.zeros((nbin, n))
    overflow_t = np.zeros((nbin, n))    # km3 leaving over the sill, per bin
    closed = np.zeros(n, dtype=bool)
    years = 0

    for years in range(1, max_years + 1):
        vol = vol0.copy()
        overflow_t[:] = 0.0
        for k in range(nbin):
            target = dt_bin[k]
            todo = np.where(live, target, 0.0)
            # Per-basin adaptive sub-stepping, on the ACTIVE subset. A global
            # step would let one fast-emptying playa set the step for every deep
            # lake in the set; carrying the finished basins through the arrays
            # anyway would cost the same as never having finished them.
            rows = np.flatnonzero(todo > 0)
            for _ in range(256):
                if rows.size == 0:
                    break
                v = vol[rows]
                a, slope = area_and_slope(v, rows)
                rate = (r[k, rows] * np.maximum(catchment[rows] - a, 0.0)
                        + (pcp[k, rows] - e[k, rows]) * a + inflow[k, rows])
                big = np.abs(rate)
                with np.errstate(divide="ignore", invalid="ignore"):
                    limit = np.where(big > 0,
                                     PERIODIC_MAX_STEP_FRACTION * cap_safe[rows] / big,
                                     np.inf)
                    # The stability bound, and it is what the accuracy bound
                    # above cannot supply. Writing the balance as
                    # dV/dt = S - D*A(V) with D = runoff + E - P, the local
                    # relaxation rate is D * dA/dV, and an explicit step longer
                    # than 2/(D dA/dV) alternates about the fixed point instead
                    # of approaching it. On this terrain that is not an edge
                    # case: a basin with a flat floor takes a hundred mesh cells
                    # of surface in one level of its curve, and the alternation
                    # reads out as a seasonal area swing under a forcing that
                    # has no season at all.
                    lam = (r[k, rows] + e[k, rows] - pcp[k, rows]) * slope
                    stable = np.where(lam > 0,
                                      PERIODIC_STABILITY_MARGIN / np.where(lam > 0, lam, 1.0),
                                      np.inf)
                step = np.minimum(todo[rows], np.minimum(limit, stable))
                raw = v + rate * step
                overflow_t[k, rows] += np.maximum(raw - capacity[rows], 0.0)
                vol[rows] = np.clip(raw, 0.0, capacity[rows])
                todo[rows] = todo[rows] - step
                rows = rows[todo[rows] > 0]
            else:
                raise RuntimeError(
                    f"bin {k} did not integrate in 256 sub-steps; the step limit "
                    "or the hypsometric curve is degenerate")
            area_t[k] = area_of(vol)
            vol_t[k] = vol

        received = np.zeros((nbin, n))
        for k in range(nbin):
            rate_out = overflow_t[k] / dt_bin[k]
            np.add.at(received[k], downstream[routes], rate_out[routes])

        drift = np.abs(vol - vol0)
        inflow_moved = np.max(np.abs(received - inflow), axis=0)
        closed = live & (drift <= tol) & (inflow_moved <= tol)
        # CLOSURE PROPAGATES DOWN THE CASCADE. A basin whose own volume repeats
        # but whose upstream supply does not is not on a periodic steady state;
        # it is riding someone else's transient, and under a constant forcing
        # that is exactly how a basin comes to be reported with a season it
        # cannot have. Most basins here spill into another rather than to the
        # ocean, so this reaches a long way.
        for _ in range(n):
            bad = live & ~closed & routes
            if not np.any(bad):
                break
            hit = np.zeros(n, dtype=bool)
            hit[downstream[bad]] = True
            if not np.any(closed & hit):
                break
            closed &= ~hit
        if np.all(closed | ~live):
            inflow = received
            break
        # Carry the year forward: the state the cycle ended in is the state it
        # starts the next one in, which is what makes this a fixed point on the
        # cycle rather than a spin-up.
        vol0 = vol
        inflow = received

    mean_area = (area_t * dt_bin[:, None]).sum(0) / dt_bin.sum()
    swing = area_t.max(axis=0) - area_t.min(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        amplitude = np.where(mean_area > 0, swing / mean_area, 0.0)
    seasonal = (closed & (amplitude >= SEASONAL_AMPLITUDE_MIN_FRACTION)
                & (swing >= float(area_resolution_km2)))

    annual_inflow = ((r * np.maximum(catchment - area_t, 0.0) + inflow)
                     * dt_bin[:, None]).sum(0)
    mean_volume = (vol_t * dt_bin[:, None]).sum(0) / dt_bin.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        residence = np.where(annual_inflow > 0, mean_volume / annual_inflow, np.inf)

    level_t = np.zeros((nbin, n))
    m = basins.area_km2.shape[1]
    rows_all = np.arange(n)
    for k in range(nbin):
        idx = np.clip((basins.area_km2 < area_t[k][:, None]).sum(1) - 1, 0, m - 2)
        a0, a1 = basins.area_km2[rows_all, idx], basins.area_km2[rows_all, idx + 1]
        l0, l1 = basins.level_km[rows_all, idx], basins.level_km[rows_all, idx + 1]
        span = a1 - a0
        frac = np.where(span > 0, (area_t[k] - a0) / np.where(span > 0, span, 1.0), 0.0)
        level_t[k] = np.where(live, l0 + np.clip(frac, 0.0, 1.0) * (l1 - l0), 0.0)
    return {
        "area_km2": area_t,
        "level_km": level_t,
        "volume_km3": vol_t,
        "overflow_km3": overflow_t,
        "mean_area_km2": mean_area,
        "mean_volume_km3": mean_volume,
        "area_swing_km2": swing,
        "seasonal_amplitude": amplitude,
        "residence_time_years": residence,
        "seasonal": seasonal,
        "closed": closed,
        "refused_no_impoundment": ~live,
        "cycles": years,
        "all_closed": bool(np.all(closed | ~live)),
    }


def _sweep(basins: BasinSet, out_path) -> dict:
    """Uniform-forcing sensitivity sweep. Machinery check, not a result."""
    # From the export. See the note at the same denominator in `carve_verdict`:
    # the areas being divided are Orogen's, so the sphere has to be Orogen's,
    # and the literal this replaces was a config value copied into a script.
    planet_km2 = Export(grid_export()).surface_area_km2
    rows = []
    for runoff_mm in (10.0, 50.0, 200.0):
        for evap_mm in (400.0, 800.0, 1600.0):
            r = np.full(basins.n, runoff_mm * 1e-6)         # mm/yr -> km/yr
            e = np.full(basins.n, evap_mm * 1e-6)
            p = np.full(basins.n, runoff_mm * 1e-6)          # placeholder
            s = solve(basins, r, e, p)
            rows.append({
                "runoff_mm_per_year": runoff_mm,
                "lake_evaporation_mm_per_year": evap_mm,
                "converged": bool(s["converged"]),
                "lake_area_fraction_of_planet": float(s["area_km2"].sum() / planet_km2),
                "lake_volume_km3": float(s["volume_km3"].sum()),
                "basins_filling_to_spill": int(s["fills_to_spill"].sum()),
                "basins_dry": int(s["dry"].sum()),
            })
    result = {
        "note": ("Uniform placeholder forcing. This exercises the solver and shows "
                 "how lake extent scales; it is not a claim about this world's "
                 "climate. surface_water.py runs this solver under the real "
                 "climatology; this sweep only exercises the machinery."),
        "terrain_hash": basins.terrain_hash,
        "basins": basins.n,
        "capacity_km3": float(basins.capacity_km3.sum()),
        "sweep": rows,
    }
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def _synthetic(n_basins=6, levels=64, *, staircase=False) -> BasinSet:
    """A BasinSet with no file behind it, for checks that need no build.

    `staircase` reproduces the one feature of real terrain that breaks a naive
    integrator: a flat basin floor that takes a large share of the basin's
    surface in a single level of the curve. It is a fixture built to be hard in
    a named way, not a convenience.
    """
    b = BasinSet.__new__(BasinSet)
    lvl = np.linspace(0.0, 1.0, levels)
    if staircase:
        shape = np.where(lvl < 0.3, lvl / 3.0, 0.1 + 0.9 * (lvl - 0.3) / 0.7)
        shape = np.where(lvl < 0.32, shape, shape + 0.6)
        shape = shape / shape[-1]
    else:
        shape = lvl ** 2
    # The basins differ in DEPTH at one surface area and one catchment ratio, so
    # capacity spans a factor of a hundred while supply does not. That is what
    # makes residence time the axis these checks vary, which is the axis the
    # seasonal amplitude is a function of.
    depth_km = np.geomspace(0.02, 2.0, n_basins)[:, None]
    b.level_km = depth_km * lvl[None, :]
    b.area_km2 = np.tile(1e3 * shape, (n_basins, 1))
    b.volume_km3 = np.concatenate(
        [np.zeros((n_basins, 1)),
         np.cumsum(0.5 * (b.area_km2[:, 1:] + b.area_km2[:, :-1])
                   * np.diff(b.level_km, axis=1), axis=1)], axis=1)
    b.spill_km = depth_km[:, 0].copy()
    b.sink_elevation_km = np.zeros(n_basins)
    b.spill_target = np.full(n_basins, -1)
    # Chosen so the steady lake sits at about a third of the area at spill: a
    # set pinned at its own spill has no room to swing and would pass a
    # seasonal check by having no season available to it.
    b.catchment_km2 = 5.0 * b.area_km2[:, -1]
    b.capacity_km3 = b.volume_km3[:, -1].copy()
    b.area_at_spill_km2 = b.area_km2[:, -1].copy()
    b.terrain_hash = "synthetic"
    b.n = n_basins
    return b


def _selftest() -> int:
    """Checks with right answers, not outcomes that could only differ.

    Every one of these can fail. Three are identities the periodic solve must
    satisfy against the equilibrium solve it generalises, one is a monotonicity
    that a solver reading the wrong forcing would break, and one is the fixture
    that caught the integrator: a flat basin floor makes the balance stiff, and
    an explicit step taken without reference to the local slope alternates
    across it and reports a season under a forcing that has none.
    """
    problems: list[str] = []
    n_checks = 0

    def check(name: str, ok: bool, detail: str = "") -> None:
        nonlocal n_checks
        n_checks += 1
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name}{'' if ok else ': ' + detail}")
        if not ok:
            problems.append(name)

    nbin = 12
    dt = np.full(nbin, 1.0 / nbin)

    for label, staircase in (("smooth", False), ("flat-floored", True)):
        b = _synthetic(staircase=staircase)
        r = np.full((nbin, b.n), 60e-6)
        e = np.full((nbin, b.n), 900e-6)
        pc = np.full((nbin, b.n), 60e-6)
        eq = solve(b, r[0], e[0], pc[0])
        per = solve_periodic(b, r, e, pc, dt, area_resolution_km2=1.0)
        check(f"a {label} basin closes its year under a steady forcing",
              bool(per["all_closed"]), f"{int((~per['closed']).sum())} did not")
        mean_area = (per["area_km2"] * dt[:, None]).sum(0)
        rel = np.abs(mean_area - eq["area_km2"]) / np.maximum(eq["area_km2"], 1e-12)
        check(f"a {label} steady forcing reproduces the equilibrium solve",
              float(rel.max()) < 1e-6, f"max relative difference {rel.max():.3g}")
        check(f"a {label} steady forcing flags no basin as seasonal",
              not bool(per["seasonal"].any()),
              f"{int(per['seasonal'].sum())} flagged, max amplitude "
              f"{per['seasonal_amplitude'].max():.3g}")

    # A season in the forcing must show up as a season in the lake, and more of
    # it must show more. A solver that read the annual mean would return the
    # same amplitude for both.
    b = _synthetic()
    phase = np.cos(2 * np.pi * (np.arange(nbin) + 0.5) / nbin)[:, None]
    amps = []
    for swing in (0.0, 0.3, 0.6):
        e = 900e-6 * (1.0 + swing * phase) * np.ones((1, b.n))
        per = solve_periodic(b, np.full((nbin, b.n), 60e-6), e,
                             np.full((nbin, b.n), 60e-6), dt,
                             area_resolution_km2=1.0)
        amps.append(float(np.median(per["seasonal_amplitude"])))
    check("a bigger seasonal forcing gives a bigger seasonal lake",
          amps[0] < amps[1] < amps[2], f"amplitudes {amps}")

    # The cycle LENGTH has to reach the answer, and it is the one input a
    # closure test cannot check: a fixed point does not care how long a year is,
    # so a caller passing bin lengths in the wrong year unit would get a wrong
    # amplitude and a clean closure. Stretching the cycle gives a lake longer to
    # follow its forcing, so the swing must grow.
    b = _synthetic()
    e = 900e-6 * (1.0 + 0.5 * phase) * np.ones((1, b.n))
    stretched = []
    for cycle in (0.5, 1.0, 2.0):
        per = solve_periodic(b, np.full((nbin, b.n), 60e-6), e,
                             np.full((nbin, b.n), 60e-6),
                             np.full(nbin, cycle / nbin), area_resolution_km2=1.0)
        stretched.append(float(np.median(per["seasonal_amplitude"])))
    check("a longer cycle gives the lake more of its forcing to follow",
          stretched[0] < stretched[1] < stretched[2], f"amplitudes {stretched}")

    # Residence time is what sets the amplitude, which is the whole reason the
    # product is per basin. The synthetic set spans a factor of 100 in capacity
    # at one catchment ratio, so the deep end must damp.
    b = _synthetic(n_basins=8)
    e = 900e-6 * (1.0 + 0.5 * phase) * np.ones((1, b.n))
    per = solve_periodic(b, np.full((nbin, b.n), 60e-6), e,
                         np.full((nbin, b.n), 60e-6), dt, area_resolution_km2=1.0)
    amp, res = per["seasonal_amplitude"], per["residence_time_years"]
    fin = np.isfinite(res) & (amp > 0)
    corr = float(np.corrcoef(np.log(res[fin]), np.log(amp[fin]))[0, 1]) if fin.sum() > 2 else 0.0
    check("a longer residence time damps the seasonal swing", corr < -0.5,
          f"log-log correlation {corr:.3f}")

    # A basin the conditioning flattened has no state to integrate. WORLD-O3B6.
    b = _synthetic()
    b.capacity_km3 = b.capacity_km3.copy()
    b.capacity_km3[0] = 0.0
    b.volume_km3 = b.volume_km3.copy()
    b.volume_km3[0] = 0.0
    per = solve_periodic(b, np.full((nbin, b.n), 60e-6), np.full((nbin, b.n), 900e-6),
                         np.full((nbin, b.n), 60e-6), dt, area_resolution_km2=1.0)
    check("a basin with no impoundment is refused rather than integrated",
          bool(per["refused_no_impoundment"][0]) and not bool(per["closed"][0])
          and not bool(per["seasonal"][0]),
          "it was integrated")

    print(f"\n{n_checks} checks, {len(problems)} failed")
    return 1 if problems else 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--basins", default=None)
    ap.add_argument("--output", default=None)
    ap.add_argument("--selftest", action="store_true",
                    help="run the periodic solver's identity and fixture "
                         "checks; needs no build")
    args = ap.parse_args()
    if args.selftest:
        raise SystemExit(_selftest())
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    basins = BasinSet(args.basins)
    out = ANALYSIS / "lake_balance_sweep.json" if args.output is None else __import__(
        "pathlib").Path(args.output)
    res = _sweep(basins, out)
    print(f"{basins.n} basins, total capacity {res['capacity_km3']:,.0f} km3")
    print(f"{'runoff':>8} {'lake evap':>10} {'lake area':>11} {'volume km3':>13} "
          f"{'at spill':>9} {'dry':>6} {'conv':>5}")
    for r in res["sweep"]:
        print(f"{r['runoff_mm_per_year']:7.0f}m {r['lake_evaporation_mm_per_year']:9.0f}m "
              f"{r['lake_area_fraction_of_planet']*100:10.3f}% "
              f"{r['lake_volume_km3']:13,.0f} {r['basins_filling_to_spill']:9d} "
              f"{r['basins_dry']:6d} {str(r['converged']):>5}")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()


def carve_verdict(basins: BasinSet, aridity_index: np.ndarray, land_area_km2: np.ndarray):
    """Which basins overflow persistently, and so should not be basins at all.

    A basin that overflows year on year incises its outlet. Over the 1e4 to 1e6
    years a landscape needs to relax, that cuts the sill down and drains the
    lake, and the depression stops existing. So `solve()` returning "pinned at
    spill, overflowing forever" is a transient, not a landscape state.

    The verdict is climate-free per basin up to a single number. From the
    balance, a basin overflows exactly when

        (E - P) / runoff  <=  catchment / area_at_spill - 1

    and the right-hand side is pure geometry (the `critical_aridity_index`
    variable in basins.nc), recomputed here from catchment and area at spill.
    Pass the left-hand side per basin, from climate, and get back the basins
    whose outlets the terrain should have carved.

    Returns (carve, endorheic_land_km2). The area accounts for catchments that
    pass through a carved basin on their way somewhere else, so it is the land
    that still has nowhere to drain once the carving is done. Divide by total
    land yourself; passing only basin-draining land as `land_area_km2` and
    dividing by its own sum silently answers a different question.
    """
    crit = basins.catchment_km2 / np.maximum(basins.area_at_spill_km2, 1e-9) - 1.0
    carve = np.asarray(aridity_index) <= crit

    eff = np.arange(basins.n)
    for _ in range(basins.n):
        nxt = np.where(carve[eff] & (basins.spill_target[eff] >= 0),
                       basins.spill_target[eff], eff)
        nxt = np.where(carve[eff] & (basins.spill_target[eff] < 0), -1, nxt)
        nxt = np.where(eff < 0, -1, nxt)
        if np.array_equal(nxt, eff):
            break
        eff = nxt
    else:
        raise RuntimeError("carve cascade did not settle")

    return carve, float(np.sum(np.asarray(land_area_km2)[eff >= 0]))
