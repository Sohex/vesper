#!/usr/bin/env python3
"""Verify the design-flux band projection on fields whose answer is known in advance.

WORLDBUILDING CONTEXT: Vesper is a fictional planet. Everything below is a
property of `derive_design_flux.py`'s arithmetic, checked on synthetic band
fields rather than on any run. Nothing here reads a climatology, spawns a
process or touches the model.

WHAT THIS EXISTS TO PIN. `derive_design_flux.py` used to refuse when any band's
warm-season AMPLIFICATION -- the band's warmest-bin move divided by the global
annual one -- came out negative, on the ground that "the projection cannot be
trusted". Two things were wrong with that guard and both are checked here.

  * THE PROJECTION IS A STRAIGHT LINE through two measured points, so it is
    defined for either sign of its slope. A band whose warmest bin cools as the
    star brightens is a MEASUREMENT of this world, and the linear form carries
    it exactly. The sign test was not NECESSARY: case `negative_but_ordered`
    below is a projection with a negative warm response that stays a field
    across the whole declared candidate range, and the old guard would have
    refused it.
  * WHAT ACTUALLY BREAKS is the ordering. `warm` and `cold` are the maximum and
    the minimum of the SAME twelve output bins, so warm >= cold is an identity
    of the field. A band whose cold response is much larger than its warm one
    has a seasonal range that closes at some flux, and the sign of the warm
    response says nothing about where. The sign test was not SUFFICIENT: case
    `positive_but_inverts` below has every warm response POSITIVE and still puts
    a cell's warmest bin below its coldest inside the candidate range, and the
    old guard would have passed it.

  * AND THE LINE MUST PASS THROUGH ITS OWN POINTS. The superseded form converted
    flux to kelvin with `lib/sensitivity.py`'s canonical global slope while
    normalising each band's move by the span the two sources actually showed.
    The two are different numbers, so the fitted line missed the very point it
    was fitted to. `anchor_identity` below states the identity -- the projection
    at the bracket flux reproduces the bracket's own band means -- and its
    control reproduces the superseded form and requires it to FAIL, so the
    check is known to have teeth rather than assumed to.

Run it:

    python exoplasim/scripts/verify_design_flux_projection.py

Exits 0 when every case holds and 1 on the first that does not.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import CONFIG  # noqa: E402,F401  (puts lib/ on sys.path)

import derive_design_flux as ddf                       # noqa: E402
import sensitivity                                     # noqa: E402


F0 = 0.945          # the base point's flux; the fixtures' anchor
F1 = 1.000          # the bracket point's flux
NLON = 4


def fixture(warm_response: dict[int, float], cold_response: dict[int, float],
            base_warm_c: float = 20.0, base_cold_c: float = -10.0,
            d_global: float = -7.0) -> tuple[dict, dict, np.ndarray, int]:
    """A base and a bracket field on one row per 10-degree band.

    One latitude row per band and `NLON` cells in each, every cell land, so a
    band mean is the row mean and the arithmetic under test is not entangled
    with the land mask or the Gaussian weights. `warm_response` and
    `cold_response` are in KELVIN PER UNIT FLUX RATIO, keyed by band index; a
    band absent from either dict gets zero.
    """
    lat = np.array([b * 10.0 - 85.0 for b in range(18)])
    bands = ddf.band_index(lat, 10.0)
    nb = int(bands.max()) + 1
    shape = (len(lat), NLON)
    base_warm = np.full(shape, base_warm_c + 273.15)
    base_cold = np.full(shape, base_cold_c + 273.15)
    other_warm = np.empty(shape)
    other_cold = np.empty(shape)
    for b in range(nb):
        rows = bands == b
        # t_other = t_base - response * (f0 - f1), so the response the module
        # recovers is exactly the one asked for here.
        other_warm[rows, :] = base_warm[rows, :] - warm_response.get(b, 0.0) * (F0 - F1)
        other_cold[rows, :] = base_cold[rows, :] - cold_response.get(b, 0.0) * (F0 - F1)
    land = np.ones(shape, dtype=bool)
    area = np.full(shape, 1.0 / (len(lat) * NLON))
    base = {"warm": base_warm, "cold": base_cold, "lat": lat, "land": land,
            "area": area, "q_warm": None, "annual_global": 0.0}
    other = {"warm": other_warm, "cold": other_cold, "lat": lat, "land": land,
             "area": area, "q_warm": None, "annual_global": -d_global}
    return base, other, bands, nb


def superseded_projection(base: dict, other: dict, bands: np.ndarray, nb: int,
                          d_global: float):
    """The form this file exists to keep out: a ratio carried by a foreign slope.

    Reproduced here and nowhere else, so `anchor_identity`'s control is the
    arithmetic that was actually in the file rather than a description of it.
    """
    slope = sensitivity.SLOPE_K_PER_FLUX_RATIO
    amp = {"warm": np.zeros(nb), "cold": np.zeros(nb)}
    for b in range(nb):
        rows = bands == b
        for season in ("warm", "cold"):
            t_base = ddf.band_mean(base[season], rows, base["land"], base["area"])
            t_other = ddf.band_mean(other[season], rows, base["land"], base["area"])
            amp[season][b] = (t_base - t_other) / d_global
    warm_cells = amp["warm"][bands][:, None]
    cold_cells = amp["cold"][bands][:, None]

    def projected(f):
        dT = slope * (f - F0)
        return (base["warm"] + warm_cells * dT - 273.15,
                base["cold"] + cold_cells * dT - 273.15)
    return projected


def band_means_c(field: np.ndarray, bands: np.ndarray, nb: int,
                 land: np.ndarray, area: np.ndarray) -> np.ndarray:
    return np.array([ddf.band_mean(field, bands == b, land, area)
                     for b in range(nb)])


def case_anchor_identity() -> list[str]:
    """The projection reproduces the bracket point's band means at the bracket flux.

    A line fitted to two points passes through both of them. That is an
    identity, so the tolerance is floating-point and not a judgement.
    """
    problems = []
    d_global = -7.196
    warm = {b: -30.0 + 4.0 * b for b in range(18)}
    cold = {b: 90.0 - 3.0 * b for b in range(18)}
    base, other, bands, nb = fixture(warm, cold, d_global=d_global)
    response, _ = ddf.band_response(base, other, bands, nb, F0, F1, d_global)
    projected, _ = ddf.band_projector(base, bands, response, F0)
    got_warm, got_cold, _ = projected(F1)
    for season, got in (("warm", got_warm), ("cold", got_cold)):
        want = band_means_c(other[season], bands, nb, base["land"], base["area"]) - 273.15
        err = float(np.abs(band_means_c(got + 273.15, bands, nb, base["land"],
                                        base["area"]) - 273.15 - want).max())
        if err > 1e-9:
            problems.append(
                f"the projection at the bracket flux {F1} misses the bracket's "
                f"own {season}-season band means by up to {err:.4f} K. A line "
                f"through two points passes through both of them; this one "
                f"does not, so it is not the line it says it is.")

    # THE CONTROL. The superseded form on the same fixture must fail the same
    # identity, or this case cannot tell a right projection from a wrong one.
    # THE CONTROL, stated as its own identity rather than as a tolerance. The
    # superseded form carried each band's move by
    # slope * (f1 - f0) / (base - bracket global span) instead of by the move
    # itself, so it overshoots every band by exactly that ratio less one. If
    # that is not what it does, this case is not reproducing the form it claims
    # to keep out.
    old = superseded_projection(base, other, bands, nb, d_global)
    old_warm, _ = old(F1)
    got = band_means_c(old_warm + 273.15, bands, nb, base["land"], base["area"])
    want = band_means_c(other["warm"], bands, nb, base["land"], base["area"])
    at_base = band_means_c(base["warm"], bands, nb, base["land"], base["area"])
    overshoot = abs(sensitivity.SLOPE_K_PER_FLUX_RATIO * (F1 - F0) / d_global) - 1.0
    predicted = np.abs(at_base - want) * overshoot
    err = float(np.abs(np.abs(got - want) - predicted).max())
    if err > 1e-9:
        problems.append(
            f"the control does not reproduce the superseded form: its miss at "
            f"the anchor departs from {overshoot:.4f} of each band's own move "
            f"by up to {err:.6f} K. The identity above is only known to have "
            f"teeth while this reproduces the arithmetic it replaced.")
    if float(predicted.max()) < 0.1:
        problems.append(
            f"the fixture's largest band move is too small for the superseded "
            f"form's {overshoot:.1%} overshoot to show: it would miss by at "
            f"most {predicted.max():.4f} K, which is below anything a reader "
            f"would call a failure. Widen the fixture's responses.")
    return problems


def case_negative_but_ordered() -> list[str]:
    """A negative warm response the old sign test would have refused, and should not.

    Every band's warmest bin cools as the star brightens and every band's
    coldest bin warms, but the seasonal ranges are wide enough that none of them
    closes anywhere in the declared candidate range. The projection is a field
    at every candidate, so there is nothing here for a guard to refuse.
    """
    problems = []
    warm = {b: -5.0 for b in range(18)}
    cold = {b: 20.0 for b in range(18)}
    base, other, bands, nb = fixture(warm, cold, base_warm_c=25.0,
                                     base_cold_c=-45.0)
    response, rows = ddf.band_response(base, other, bands, nb, F0, F1, -7.0)
    if not all(r["warm_k_per_unit_flux"] < 0 for r in rows):
        problems.append("the fixture no longer has a negative warm response in "
                        "every band, so it does not test what it says it does")
    if not all(r["amp_warm"] < 0 for r in rows):
        problems.append("the reported amplification ratio no longer carries the "
                        "sign of the response; the old guard's input has moved")
    _, inverted = ddf.band_projector(base, bands, response, F0)
    bad = [f for f in ddf.DECLARED["candidates"] if inverted(f).any()]
    if bad:
        problems.append(
            f"the fixture inverts at {len(bad)} candidates ({bad[0]} to "
            f"{bad[-1]}); it is meant to stay a field across the whole declared "
            f"range so that the only thing the old sign test could have "
            f"objected to is the sign")
    return problems


def case_positive_but_inverts() -> list[str]:
    """A positive warm response the old sign test would have passed, and should not.

    Every band's warmest bin warms with the star, so every amplification is
    positive and the old guard is silent. The cold response is several times
    larger, which is what this world's high latitudes actually show, so the
    seasonal range closes and the projection stops being a field inside the
    declared candidate range.
    """
    problems = []
    warm = {b: 20.0 for b in range(18)}
    cold = {b: 300.0 for b in range(18)}
    base, other, bands, nb = fixture(warm, cold, base_warm_c=20.0,
                                     base_cold_c=-10.0)
    response, rows = ddf.band_response(base, other, bands, nb, F0, F1, -7.0)
    if not all(r["warm_k_per_unit_flux"] > 0 for r in rows):
        problems.append("the fixture no longer has a positive warm response in "
                        "every band, so the old sign test would have caught it "
                        "and the case proves nothing")
    if not all(r["amp_warm"] > 0 for r in rows):
        problems.append("the reported amplification ratio is not positive in "
                        "every band; the old guard would have fired here and "
                        "this case no longer separates the two tests")
    _, inverted = ddf.band_projector(base, bands, response, F0)
    bad = [f for f in ddf.DECLARED["candidates"] if inverted(f).any()]
    if not bad:
        problems.append(
            "the fixture never inverts inside the declared candidate range, so "
            "it does not show that a positive warm response can still leave the "
            "projection with a warmest bin below its coldest")
    return problems


def case_validity_window_is_an_interval() -> list[str]:
    """The candidates the projection is a field on form one contiguous block.

    Each cell's projected seasonal range is linear in flux, so the set where
    every cell's range is non-negative is an intersection of half-lines and
    therefore an interval. `derive()` asserts this before it scores anything;
    a gap would mean the arithmetic is not what the module says it is.
    """
    problems = []
    warm = {b: 20.0 - 2.0 * b for b in range(18)}
    cold = {b: 40.0 + 5.0 * b for b in range(18)}
    base, other, bands, nb = fixture(warm, cold, base_warm_c=20.0,
                                     base_cold_c=-30.0)
    response, _ = ddf.band_response(base, other, bands, nb, F0, F1, -7.0)
    _, inverted = ddf.band_projector(base, bands, response, F0)
    ok = [not inverted(f).any() for f in ddf.DECLARED["candidates"]]
    if not any(ok):
        problems.append("the fixture is a field at no candidate at all, so the "
                        "contiguity it is meant to demonstrate is vacuous")
        return problems
    if not all(ok):
        first, last = ok.index(True), len(ok) - 1 - ok[::-1].index(True)
        if not all(ok[first:last + 1]):
            problems.append(
                "the candidates the projection is a field on are not "
                "contiguous, which a response linear in flux cannot produce")
    return problems


def case_no_sign_guard_remains() -> list[str]:
    """No refusal in the module is spelled against the sign of a response.

    Named rather than inferred: a reinstated sign test would most likely come
    back under its old key, and a chosen flux standing on it would be a number
    with a retracted argument under it.
    """
    problems = []
    if "band_amplification_sign" in ddf.CHOOSABLE_REFUSALS:
        problems.append(
            "derive_design_flux.CHOOSABLE_REFUSALS still offers "
            "`band_amplification_sign` as a ground a chosen flux may stand on. "
            "The sign of a band's warm-season response is a measurement of this "
            "world, and a projection that refuses on it refuses on its data.")
    return problems


def main() -> None:
    cases = [
        ("the projection passes through both of its own points",
         case_anchor_identity),
        ("a negative warm response is projected, not refused",
         case_negative_but_ordered),
        ("a positive warm response does not certify the projection",
         case_positive_but_inverts),
        ("the validity window is one interval",
         case_validity_window_is_an_interval),
        ("no refusal is spelled against the sign of a response",
         case_no_sign_guard_remains),
    ]
    failed = 0
    for name, run in cases:
        problems = run()
        if problems:
            failed += 1
            print(f"[ FAIL ] {name}", flush=True)
            for p in problems:
                print(f"         {p}", flush=True)
        else:
            print(f"[  ok  ] {name}", flush=True)
    print(f"\n{len(cases)} cases, {failed} failed")
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
