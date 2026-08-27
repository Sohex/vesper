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
            # The catalogue's finalPreserved.retainedFraction: final over
            # natural spill depth, both in the generator's dimensionless
            # elevation parameter. It is the ONE quantity that converts a
            # retain computed against the finished depression in metres into
            # the basis Orogen spends it against, and `export_carve_list.py`
            # is its only consumer. Optional because basins.nc predates the
            # field; a consumer that needs it says so by name rather than
            # every other consumer failing to open the file.
            self.retained_fraction = (
                np.asarray(ds["retained_fraction"][:])
                if "retained_fraction" in ds.variables else None)
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

# THERE IS NO STEP SIZE HERE, AND THAT IS THE POINT. Within one time bin the
# forcing is constant and the hypsometric curve is piecewise linear, so
# `dV/dt = S - D A(V)` is a LINEAR ordinary differential equation on every
# segment of the curve, and `_integrate_bin` walks the segments solving each in
# closed form. Nothing is discretised, so there is no accuracy fraction and no
# stability margin to choose; both used to stand here and both were numbers
# picked to be small rather than derived from anything.

# A basin HAS a seasonal cycle worth publishing when both hold. The fraction is
# a statement about the basin; the area floor is a statement about the
# instrument, and it is the caller's mean mesh cell area, because a lake area
# comes off a curve built from mesh cells and a swing smaller than one cell is
# below what the curve can express.
SEASONAL_AMPLITUDE_MIN_FRACTION = 0.10


# ---------------------------------------------------------------------------
# The exact bin integrator
# ---------------------------------------------------------------------------
# `phi1`, `phi2` and `psi` are the three entire functions the closed-form step
# is written in. Each is a ratio that is 0/0 at the origin, so each has a series
# branch, and the branch point is where the direct form starts losing digits to
# cancellation rather than where the ratio is undefined.


def _phi1(x):
    """(e^x - 1) / x, and 1 at x = 0."""
    x = np.asarray(x, dtype=float)
    small = np.abs(x) < 1e-8
    xs = np.where(small, 1.0, x)
    return np.where(small, 1.0 + x / 2.0 + x * x / 6.0, np.expm1(xs) / xs)


def _phi2(x):
    """(e^x - 1 - x) / x^2, and 1/2 at x = 0.

    The direct form cancels the leading `x` out of `expm1`, so it loses about
    `-log10|x|` digits; below 1e-4 the four-term series is good to 1e-19.
    """
    x = np.asarray(x, dtype=float)
    small = np.abs(x) < 1e-4
    xs = np.where(small, 1.0, x)
    series = 0.5 + x / 6.0 + x * x / 24.0 + x * x * x / 120.0
    return np.where(small, series, (np.expm1(xs) - xs) / (xs * xs))


def _psi(z):
    """-log(1 - z) / z, and 1 at z = 0. Defined for z < 1."""
    z = np.asarray(z, dtype=float)
    small = np.abs(z) < 1e-8
    zs = np.where(small, 1.0, z)
    return np.where(small, 1.0 + z / 2.0 + z * z / 3.0, -np.log1p(-zs) / zs)


def _sub_mesh_area(supply_km3_per_year, demand_km_per_year, area_at_zero_km2):
    """The lake a basin holding NO STORAGE still has, and why it is not zero.

    A hypsometric curve is built from mesh cells, so its first level already
    carries a whole cell of surface: `basins.nc` gives every basin a positive
    area at zero volume, up to a couple of hundred km2. A basin whose supply
    cannot fill that first cell therefore has a lake the curve cannot express,
    somewhere between nothing and one cell, and reading the curve's own first
    area for it charges open-water evaporation over water that is not there.
    Over this catalogue that is a tenth of a percent of everything arriving,
    destroyed rather than evaporated.

    What the balance can say without the curve is the whole of the answer: with
    no storage to draw down, the lake is exactly the area whose evaporative
    demand consumes the supply, `S / D`, and it is below one cell by the same
    condition that put the basin here. This is the same area `solve()` returns
    at the bottom of the curve, which reads `supply / demand` with no floor
    under it, so the two solves agree at the dry end instead of differing by a
    mesh cell.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        d = np.asarray(demand_km_per_year, dtype=float)
        return np.where(d > 0,
                        np.clip(np.asarray(supply_km3_per_year, dtype=float)
                                / np.where(d > 0, d, 1.0), 0.0, area_at_zero_km2),
                        area_at_zero_km2)


def _integrate_bin(vcurve, acurve, capacity, volume_km3, supply_km3_per_year,
                   demand_km_per_year, span_years):
    """Integrate one basin-storage balance across one time bin, in closed form.

    THE BALANCE. Over a bin the forcing does not change, so with `S` the volume
    arriving per year and `D` the net depth leaving per year over open water,

        dV/dt = S - D A(V)

    and `A(V)` is the basin's hypsometric curve, which `basins.nc` publishes as
    a piecewise linear function of volume. On one segment of that curve the
    right-hand side is LINEAR IN V, so the equation is a linear ODE and its
    solution is an exponential. This walks the segments the basin crosses and
    solves each one exactly. No step size is chosen anywhere: what is computed
    per segment is the TIME THE BASIN TAKES TO CROSS IT, and the bin ends inside
    whichever segment the remaining time runs out in.

    WHY THAT IS NOT AN OPTIMISATION. An explicit step has to be short enough to
    be stable against `D dA/dV`, and short enough that a basin taking a large
    throughflow does not move far in one step; on this terrain the second bound
    goes to zero for a small reservoir on a big river, because the bound is set
    by the RATE while the basin, pinned at its spill, does not move at all. The
    closed form has neither bound, and it charges the correct open-water demand
    while the basin is spilling instead of the demand at the level it started
    the step from.

    THE WORK IS BOUNDED BY THE CURVE. `A` is non-decreasing in `V`, so the
    right-hand side is monotone in `V`, so within a bin the storage moves in one
    direction only and crosses each of the curve's levels at most once. The
    iteration limit is therefore the curve's own length and cannot be tuned: if
    it is ever reached, either a curve is not monotone or a step made no
    progress, and both are defects rather than a bin that needed more steps.

    `vcurve` and `acurve` are `(nrows, nlevels)`, volume-ascending, and the
    remaining arguments are `(nrows,)` except `span_years`. Returns
    `(volume_end, overflow_km3, volume_integral, area_integral, area_end)`, the
    two integrals being over the bin, in km3 yr and km2 yr. A basin at its
    capacity with water still arriving holds there and the surplus goes to
    `overflow_km3`; a basin at zero with water still leaving holds there, with
    the lake `_sub_mesh_area` gives it rather than the curve's own first cell.
    """
    v = np.array(volume_km3, dtype=float)
    n, m = vcurve.shape
    span = float(span_years)
    overflow = np.zeros(n)
    int_v = np.zeros(n)
    int_a = np.zeros(n)
    todo = np.full(n, span)
    active = np.flatnonzero(todo > 0)

    for _ in range(m + 2):
        if active.size == 0:
            break
        vc, ac = vcurve[active], acurve[active]
        rows = np.arange(active.size)
        vi, ti = v[active], todo[active]
        s_i, d_i = supply_km3_per_year[active], demand_km_per_year[active]
        cap_i = capacity[active]

        # The level the basin sits at, and the area there. Continuous at a
        # level boundary, so which side it is read from does not matter.
        count_lt = (vc < vi[:, None]).sum(1)
        count_le = (vc <= vi[:, None]).sum(1)
        j = np.clip(count_le - 1, 0, m - 2)
        j0, j1 = vc[rows, j], vc[rows, j + 1]
        a0, a1 = ac[rows, j], ac[rows, j + 1]
        width = j1 - j0
        frac = np.where(width > 0, (vi - j0) / np.where(width > 0, width, 1.0), 0.0)
        a_here = a0 + np.clip(frac, 0.0, 1.0) * (a1 - a0)
        rate = s_i - d_i * a_here

        # Three states that hold for the rest of the bin, because the rate at
        # them cannot change while the basin does not move: spilling at
        # capacity, dry at zero, and sitting on a fixed point.
        spilling = (rate > 0) & (vi >= cap_i)
        drying = (rate < 0) & (vi <= vc[rows, 0])
        done = spilling | drying | (rate == 0)
        if done.any():
            w = active[done]
            hold = ti[done]
            int_v[w] += vi[done] * hold
            int_a[w] += np.where(drying[done], _sub_mesh_area(s_i, d_i, ac[:, 0])[done],
                                 a_here[done]) * hold
            overflow[w] += np.where(spilling[done], rate[done], 0.0) * hold
            todo[w] = 0.0

        move = ~done
        if not move.any():
            active = active[todo[active] > 0]
            continue

        # The next level in the direction of travel, taken STRICTLY so a
        # repeated volume in the curve cannot stall the walk.
        up = rate > 0
        ib = np.where(up, np.clip(count_le, 1, m - 1), np.clip(count_lt - 1, 0, m - 2))
        v_edge = vc[rows, ib]
        js = np.where(up, ib - 1, ib)
        e0, e1 = vc[rows, js], vc[rows, js + 1]
        b0, b1 = ac[rows, js], ac[rows, js + 1]
        ewidth = e1 - e0
        slope = np.where(ewidth > 0, (b1 - b0) / np.where(ewidth > 0, ewidth, 1.0), 0.0)

        # The relaxation rate on this segment, and the time to cross it. `z` is
        # how much of the way to the segment's own fixed point the far edge is:
        # at or past it the edge is never reached and the bin ends inside.
        lam = d_i * slope
        edge = v_edge - vi
        with np.errstate(divide="ignore", invalid="ignore"):
            safe_rate = np.where(rate != 0, rate, 1.0)
            z = np.where(rate != 0, lam * edge / safe_rate, 0.0)
            reachable = (z < 1.0) & (edge != 0)
            cross = np.where(reachable,
                             (edge / safe_rate) * _psi(np.where(reachable, z, 0.0)),
                             np.inf)
        step = np.where(move, np.minimum(cross, ti), 0.0)
        x = -lam * step
        reached = move & (cross <= ti)

        # TWO FORMS OF ONE ANSWER, and which is used is a conditioning choice
        # rather than an approximation. Written as a displacement from where the
        # storage started, the step subtracts two nearly equal numbers once the
        # basin has relaxed most of the way to the segment's fixed point, which
        # is the stiff case and so the normal one here. Written as a
        # displacement from that fixed point it does not, and the fixed point
        # itself comes off the segment's own coefficients rather than off the
        # current volume, so it carries no cancellation into the step either.
        relaxing = (np.abs(x) >= 1.0) & (slope != 0) & (d_i != 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            safe_slope = np.where(relaxing, slope, 1.0)
            safe_d = np.where(relaxing, d_i, 1.0)
            a_star = np.where(relaxing, s_i / safe_d, 0.0)
            v_star = np.where(relaxing, e0 + (a_star - b0) / safe_slope, 0.0)
        p1 = _phi1(x)
        p2 = _phi2(x)
        decay = np.exp(np.where(relaxing, x, 0.0))
        v_free = np.where(relaxing, v_star + (vi - v_star) * decay,
                          vi + rate * step * p1)
        v_new = np.where(reached, v_edge, v_free)
        dv_int = np.where(relaxing, step * (v_star + (vi - v_star) * p1),
                          vi * step + rate * step * step * p2)
        da_int = np.where(relaxing, step * (a_star + (a_here - a_star) * p1),
                          a_here * step + slope * rate * step * step * p2)
        w = active[move]
        int_v[w] += dv_int[move]
        int_a[w] += da_int[move]
        v[w] = np.clip(v_new[move], 0.0, cap_i[move])
        todo[w] = ti[move] - step[move]
        active = active[todo[active] > 0]
    else:
        if active.size:
            raise RuntimeError(
                f"{active.size} basins did not cross their hypsometric curve in "
                f"{m + 2} segment steps. The storage moves one way within a bin "
                "and each level is crossed at most once, so this means a curve "
                "whose volume is not ascending or a step that made no progress, "
                "not a bin that needed a finer integration")

    # The area the bin ENDS at, by the same rule the integral used, so the
    # sampled area and the time-mean area cannot disagree about a basin that
    # holds nothing.
    idx = np.clip((vcurve < v[:, None]).sum(1) - 1, 0, m - 2)
    rows = np.arange(n)
    c0, c1 = vcurve[rows, idx], vcurve[rows, idx + 1]
    d0, d1 = acurve[rows, idx], acurve[rows, idx + 1]
    cw = c1 - c0
    cf = np.where(cw > 0, (v - c0) / np.where(cw > 0, cw, 1.0), 0.0)
    area_end = np.where(
        v > 0.0,
        d0 + np.clip(cf, 0.0, 1.0) * (d1 - d0),
        _sub_mesh_area(supply_km3_per_year, demand_km_per_year, acurve[:, 0]))
    return v, overflow, int_v, int_a, area_end


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

    WHAT IS PUBLISHED. Per basin and per bin, the area, level and volume AT THE
    END OF THE BIN, and separately the bin's own TIME MEAN of area and volume;
    and per basin, the peak-to-trough area range, that range as a fraction of
    the cycle-mean area, the residence time, and `seasonal` -- whether both
    published-cycle criteria hold. The honest product is a per-basin amplitude
    carried by the low-capacity basins, not a claim that every preserved basin
    has a season: a deep terminal lake holds years of supply and its surface
    barely moves inside one, while a shallow playa's area is almost all
    seasonal.

    THE END-OF-BIN AND TIME-MEAN AREAS ARE DIFFERENT QUANTITIES and neither is
    the other rounded. A basin whose storage turns over inside a bin can sit at
    its capacity for most of the bin and end it far below, so the peak-to-trough
    range is taken on the samples the per-bin arrays publish, while everything
    that has to conserve water -- the cycle mean, the residence time, the annual
    balance -- is taken on the time mean.

    THE WATER BALANCE TRAVELS WITH THE ANSWER, per bin and per year, as
    `bin_balance_residual_km3` and `annual_water_residual_km3`. They are
    diagnostics rather than criteria: `closed` is still what says whether a
    basin may be read.
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
    catchment = basins.catchment_km2

    # THE BALANCE IS LINEAR IN THE LAKE'S AREA, and the integrator below is
    # built on that. Runoff is generated over the DRY catchment, so the supply
    # term is `r (C - A)` and folds into `S = r C` against `D = r + E - P` only
    # while the lake is no larger than the catchment it sits in. That is a
    # property of the catalogue rather than of the climate -- a catchment is
    # defined to contain its own basin floor -- so a violation is a defect in
    # `basins.nc` and is refused here rather than linearised away.
    over = live & (basins.area_at_spill_km2 > catchment)
    if np.any(over):
        b = int(np.flatnonzero(over)[np.argmax(
            (basins.area_at_spill_km2 - catchment)[over])])
        raise ValueError(
            f"{int(over.sum())} basins have more area at spill than catchment, "
            f"worst basin {b} at {basins.area_at_spill_km2[b]:,.1f} km2 against "
            f"{catchment[b]:,.1f} km2. A catchment contains its own basin floor, "
            "so this is a defect in the basin catalogue")

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
    mean_area_t = np.zeros((nbin, n))   # time-mean over the bin, not its end
    mean_vol_t = np.zeros((nbin, n))
    residual_t = np.zeros((nbin, n))    # the bin's own water balance, km3
    closed = np.zeros(n, dtype=bool)
    years = 0

    # The integration runs on the live subset only, and its curves are lifted
    # out of the cycle loop because they do not change.
    rows_live = np.flatnonzero(live)
    vcurve_live = basins.volume_km3[rows_live]
    acurve_live = basins.area_km2[rows_live]
    cap_live = capacity[rows_live]
    catch_live = catchment[rows_live]

    cycle_start = vol0
    for years in range(1, max_years + 1):
        cycle_start = vol0
        vol = vol0.copy()
        overflow_t[:] = 0.0
        for k in range(nbin):
            # S and D of the balance dV/dt = S - D A(V), constant across the
            # bin because the forcing is. Runoff is generated over the DRY
            # catchment, so the lake's own footprint comes out of the supply
            # and into the demand; the guard above is what licenses writing it
            # this way rather than clamping C - A at every evaluation.
            supply = r[k, rows_live] * catch_live + inflow[k, rows_live]
            demand = r[k, rows_live] + e[k, rows_live] - pcp[k, rows_live]
            v_end, ovf, iv, ia, a_end = _integrate_bin(
                vcurve_live, acurve_live, cap_live, vol[rows_live],
                supply, demand, dt_bin[k])
            # The bin's own water balance, and the two sides of it are computed
            # by different closed forms rather than by one from the other: the
            # storage change comes out of the exponential step and the loss out
            # of the exponential's integral. `_selftest` puts a tolerance on it.
            residual_t[k, rows_live] = (
                (v_end - vol[rows_live])
                - (supply * dt_bin[k] - demand * ia - ovf))
            vol[rows_live] = v_end
            overflow_t[k, rows_live] = ovf
            mean_vol_t[k, rows_live] = iv / dt_bin[k]
            mean_area_t[k, rows_live] = ia / dt_bin[k]
            area_t[k] = 0.0
            area_t[k, rows_live] = a_end
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

    # The cycle means are TIME means, from the integrator's own exact integral
    # over each bin, not the bin-end samples averaged. The two differ by the
    # curvature the storage has inside a bin, which is exactly what a basin
    # whose residence time is short compared with a bin is made of.
    mean_area = (mean_area_t * dt_bin[:, None]).sum(0) / dt_bin.sum()
    mean_volume = (mean_vol_t * dt_bin[:, None]).sum(0) / dt_bin.sum()
    swing = area_t.max(axis=0) - area_t.min(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        amplitude = np.where(mean_area > 0, swing / mean_area, 0.0)
    seasonal = (closed & (amplitude >= SEASONAL_AMPLITUDE_MIN_FRACTION)
                & (swing >= float(area_resolution_km2)))

    annual_inflow = ((r * np.maximum(catchment - mean_area_t, 0.0) + inflow)
                     * dt_bin[:, None]).sum(0)
    with np.errstate(divide="ignore", invalid="ignore"):
        residence = np.where(annual_inflow > 0, mean_volume / annual_inflow, np.inf)

    # The year's own water balance per basin: what arrived less what left less
    # what the storage gained. Zero to arithmetic on a closed cycle, and the
    # scale to read it against is what passed through, so both are returned.
    annual_supply = ((r * np.maximum(catchment - mean_area_t, 0.0)
                      + pcp * mean_area_t + inflow) * dt_bin[:, None]).sum(0)
    annual_loss = (e * mean_area_t * dt_bin[:, None]).sum(0)
    annual_overflow = overflow_t.sum(0)
    storage_change = vol - cycle_start
    water_residual = np.where(
        live, annual_supply - annual_loss - annual_overflow - storage_change, 0.0)

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
        "mean_area_by_bin_km2": mean_area_t,
        "mean_volume_by_bin_km3": mean_vol_t,
        "mean_area_km2": mean_area,
        "mean_volume_km3": mean_volume,
        # The conservation evidence, per basin: the bin's own balance residual
        # from the integrator's two closed forms, and the year's balance from
        # the published trajectory. Both are km3, and `annual_throughput_km3`
        # is the scale to read them against.
        "bin_balance_residual_km3": residual_t,
        "annual_water_residual_km3": water_residual,
        "annual_throughput_km3": annual_supply,
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


# The tolerances the integrator is judged on, DECLARED BEFORE IT WAS RUN. Each
# is a distance from an answer that is known independently, so each is set by
# what double precision can carry over the work involved and not by what came
# out. They are unrelated to PERIODIC_CLOSURE_RELATIVE, which is a physical
# statement about a cycle repeating.
#
# The closed form on one segment of a curve is the exponential itself, so it is
# right to the last few bits: 1e-13 leaves room for the handful of roundings
# between the two expressions and nothing else.
EXACTNESS_RELATIVE = 1e-13
# The bin's own water balance is two closed forms of the same integral, one for
# the storage change and one for the loss. Over at most one pass of a 128-level
# curve that is a few hundred roundings against the largest term in the balance.
BALANCE_RELATIVE = 1e-10
# An independent stiff integrator, at four orders tighter than this on its own
# tolerance, so what is being measured is this solver rather than the reference.
# Read against the basin's CAPACITY, which is the scale a lake volume has.
REFERENCE_TOLERANCE = 1e-6
REFERENCE_RTOL = 1e-10
REFERENCE_ATOL_KM3 = 1e-14


def _reference_bin(vcurve, acurve, capacity, v0, supply, demand, span_years):
    """One basin, one bin, by an independent stiff integrator.

    Radau on the same balance, with the hypsometric curve read through
    `np.interp` rather than through anything in this module. It exists to
    disagree with `_integrate_bin`, so it shares no code with it; the price is
    that it takes one basin at a time and cannot be what the solver uses.

    Returns the volume at the end of the bin. Water above capacity spills, so
    the area is held at the spill area there and the surplus is not carried.
    """
    from scipy.integrate import solve_ivp

    def rhs(_t, y):
        v = min(max(float(y[0]), 0.0), float(capacity))
        a = float(np.interp(v, vcurve, acurve))
        rate = supply - demand * a
        if y[0] >= capacity and rate > 0:
            return [0.0]
        if y[0] <= 0.0 and rate < 0:
            return [0.0]
        return [rate]

    sol = solve_ivp(rhs, (0.0, float(span_years)), [float(v0)], method="Radau",
                    rtol=REFERENCE_RTOL, atol=REFERENCE_ATOL_KM3, dense_output=False)
    if not sol.success:
        raise RuntimeError(f"the reference integrator failed: {sol.message}")
    return float(np.clip(sol.y[0, -1], 0.0, capacity))


def _selftest() -> int:
    """Checks with right answers, not outcomes that could only differ.

    Every one of these can fail. Three are identities the periodic solve must
    satisfy against the equilibrium solve it generalises, one is a monotonicity
    that a solver reading the wrong forcing would break, and the rest are on the
    bin integrator: the closed form against the exponential it claims to be, the
    water balance against itself, an independent stiff integrator against the
    two shapes of terrain that break a stepped one, and the throughflow case
    where a small reservoir sits on a large river.
    """
    problems: list[str] = []
    n_checks = 0

    def check(name: str, ok: bool, detail: str = "", margin: str = "") -> None:
        """`margin` is printed whether the check passes or not.

        A numeric check that only reports itself when it fails cannot be read
        for how much room it had, and the room is what says whether the next
        change to the integrator is safe.
        """
        nonlocal n_checks
        n_checks += 1
        tail = f" [{margin}]" if margin else ""
        print(f"[{'  ok  ' if ok else ' FAIL '}] {name}{tail}"
              f"{'' if ok else ': ' + detail}")
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

    # ---------------------------------------------------------------------
    # The bin integrator, against answers that do not come from it
    # ---------------------------------------------------------------------

    # ON ONE SEGMENT THE BALANCE IS A LINEAR ODE AND ITS ANSWER IS THE
    # EXPONENTIAL. A curve of two levels is exactly one segment, so the whole
    # bin has a closed form written out here from the coefficients rather than
    # from anything the module does. Nothing about stiffness is at issue: this
    # is whether the algebra is the algebra.
    v_top, a_top = 4.0, 200.0
    vcurve = np.array([[0.0, v_top]])
    acurve = np.array([[0.0, a_top]])
    cap = np.array([v_top])
    slope = a_top / v_top
    exact_worst = 0.0
    for S, D, v0, h in ((3.0, 0.05, 0.5, 0.25), (0.4, 0.4, 3.0, 0.08),
                        (12.0, 2.0, 0.1, 0.02), (0.05, 0.9, 2.0, 0.5)):
        lam = D * slope
        vstar = S / (D * slope)
        want = vstar + (v0 - vstar) * np.exp(-lam * h)
        want = min(max(want, 0.0), v_top)
        got, _, _, _, _ = _integrate_bin(vcurve, acurve, cap, np.array([v0]),
                                      np.array([S]), np.array([D]), h)
        exact_worst = max(exact_worst, abs(got[0] - want) / max(abs(want), 1e-30))
    check("on one segment the closed form is the exponential it claims to be",
          exact_worst < EXACTNESS_RELATIVE,
          f"worst relative difference {exact_worst:.3g} against {EXACTNESS_RELATIVE:g}",
          f"{exact_worst:.2g} of {EXACTNESS_RELATIVE:g}")

    # A SMALL RESERVOIR ON A LARGE RIVER, which is the case that broke the
    # stepped integrator: the storage does not move at all while a throughflow
    # many times its capacity passes over the sill each bin. The answers are
    # arithmetic -- the volume stays at capacity and the overflow is the
    # surplus over the whole bin at the spill area's own demand.
    S, D, h = 900.0, 0.02, 0.1
    got, ovf, iv, ia, _ = _integrate_bin(vcurve, acurve, cap, np.array([v_top]),
                                      np.array([S]), np.array([D]), h)
    want_ovf = (S - D * a_top) * h
    check("a basin pinned at its spill passes its throughflow in one step",
          abs(got[0] - v_top) < EXACTNESS_RELATIVE * v_top
          and abs(ovf[0] - want_ovf) < EXACTNESS_RELATIVE * want_ovf
          and abs(ia[0] - a_top * h) < EXACTNESS_RELATIVE * a_top * h,
          f"volume {got[0]:.12g} against {v_top}, overflow {ovf[0]:.12g} "
          f"against {want_ovf:.12g}",
          f"throughflow {S * h / cap[0]:.0f}x capacity in the bin")

    # AGAINST AN INDEPENDENT STIFF INTEGRATOR, on the two shapes of terrain
    # that break a stepped one: a flat basin floor, where the area runs away
    # with the volume, and a basin whose whole capacity turns over inside a bin.
    # Every case here is one an explicit scheme has to sub-step and this one
    # does not, so a shortcut that was only approximately right would show.
    fix = _synthetic(n_basins=6, staircase=True)
    worst_ref = 0.0
    worst_case = None
    cases = []
    for bi in range(fix.n):
        capb = float(fix.capacity_km3[bi])
        for supply_mult in (0.2, 1.0, 5.0, 60.0):
            for demand in (2e-4, 9e-4, 4e-3):
                cases.append((bi, capb * supply_mult / 0.1, demand,
                              0.35 * capb, 0.1))
    for bi, S, D, v0, h in cases:
        got, _, _, _, _ = _integrate_bin(
            fix.volume_km3[bi:bi + 1], fix.area_km2[bi:bi + 1],
            fix.capacity_km3[bi:bi + 1], np.array([v0]),
            np.array([S]), np.array([D]), h)
        want = _reference_bin(fix.volume_km3[bi], fix.area_km2[bi],
                              fix.capacity_km3[bi], v0, S, D, h)
        rel = abs(got[0] - want) / float(fix.capacity_km3[bi])
        if rel > worst_ref:
            worst_ref, worst_case = rel, (bi, S, D, v0)
    check(f"a stiff bin agrees with an independent integrator over {len(cases)} cases",
          worst_ref < REFERENCE_TOLERANCE,
          f"worst {worst_ref:.3g} of capacity against {REFERENCE_TOLERANCE:g}, "
          f"at basin {worst_case}",
          f"{worst_ref:.2g} of capacity, against {REFERENCE_TOLERANCE:g}")

    # THE BIN'S WATER BALANCE, over a real cycle rather than a single bin. The
    # storage change comes out of the exponential and the loss out of the
    # exponential's integral, so this is two closed forms of one integral
    # against each other and it can disagree.
    b = _synthetic(n_basins=8, staircase=True)
    e = 900e-6 * (1.0 + 0.7 * phase) * np.ones((1, b.n))
    per = solve_periodic(b, np.full((nbin, b.n), 60e-6), e,
                         np.full((nbin, b.n), 60e-6), dt, area_resolution_km2=1.0)
    scale = np.maximum(b.capacity_km3, np.abs(per["annual_throughput_km3"]))
    rel_bin = float((np.abs(per["bin_balance_residual_km3"]).max(axis=0)
                     / np.maximum(scale, 1e-30)).max())
    check("every bin's water balance closes on itself",
          rel_bin < BALANCE_RELATIVE,
          f"worst relative residual {rel_bin:.3g} against {BALANCE_RELATIVE:g}",
          f"{rel_bin:.2g} of {BALANCE_RELATIVE:g}")

    # AND THE YEAR'S. On a closed cycle the storage returns to where it started,
    # so what arrived over the year has to equal what evaporated plus what went
    # over the sill, per basin and summed over the set.
    rel_year = float((np.abs(per["annual_water_residual_km3"])
                      / np.maximum(scale, 1e-30)).max())
    total = float(np.abs(per["annual_water_residual_km3"]).sum()
                  / max(per["annual_throughput_km3"].sum(), 1e-30))
    check("the year's water closes, per basin and over the set",
          bool(per["all_closed"]) and rel_year < BALANCE_RELATIVE
          and total < BALANCE_RELATIVE,
          f"per basin {rel_year:.3g}, over the set {total:.3g}, "
          f"against {BALANCE_RELATIVE:g}",
          f"per basin {rel_year:.2g}, over the set {total:.2g}, "
          f"against {BALANCE_RELATIVE:g}")

    # A BASIN TOO DRY TO FILL ONE MESH CELL. Every hypsometric curve in this
    # project carries a whole cell of surface at zero volume, so a basin whose
    # supply cannot fill that cell has a lake the curve cannot express and the
    # curve's own first area is not it. The two things that must hold are that
    # the water still closes and that the periodic solve gives the same area as
    # the equilibrium solve, which reads supply over demand with no floor under
    # it. Reading the curve instead fails both, by one mesh cell per basin.
    b = _synthetic(n_basins=6)
    b.area_km2 = b.area_km2 + 150.0        # a cell of surface at zero volume
    b.volume_km3 = np.concatenate(
        [np.zeros((b.n, 1)),
         np.cumsum(0.5 * (b.area_km2[:, 1:] + b.area_km2[:, :-1])
                   * np.diff(b.level_km, axis=1), axis=1)], axis=1)
    b.capacity_km3 = b.volume_km3[:, -1].copy()
    b.area_at_spill_km2 = b.area_km2[:, -1].copy()
    r_dry = np.full((nbin, b.n), 2e-7)     # far too little to fill one cell
    e_dry = np.full((nbin, b.n), 2000e-6)
    p_dry = np.full((nbin, b.n), 2e-7)
    eq = solve(b, r_dry[0], e_dry[0], p_dry[0])
    per = solve_periodic(b, r_dry, e_dry, p_dry, dt, area_resolution_km2=1.0)
    scale = np.maximum(b.capacity_km3, np.abs(per["annual_throughput_km3"]))
    dry_bal = float((np.abs(per["annual_water_residual_km3"])
                     / np.maximum(scale, 1e-30)).max())
    dry_area = float(np.abs(per["mean_area_km2"] - eq["area_km2"]).max()
                     / max(eq["area_km2"].max(), 1e-30))
    check("a basin below one mesh cell of lake conserves water and matches solve()",
          bool((per["volume_km3"] <= 0).all()) and dry_bal < BALANCE_RELATIVE
          and dry_area < 1e-12,
          f"balance {dry_bal:.3g}, area difference {dry_area:.3g}",
          f"balance {dry_bal:.2g}, area {dry_area:.2g}, "
          f"lake {per['mean_area_km2'].max():.3g} km2 under a "
          f"{b.area_km2[0, 0]:.0f} km2 first cell")

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
