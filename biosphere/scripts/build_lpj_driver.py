"""Build the binary driver file LPJ-GUESS reads, from climate and lithology.

One self-describing file carrying the gridlist, a soil code per cell derived from
World Orogen lithology, the land column property contract's soil hydraulic
states, and binned climate. `vesperinput.cpp` reads it.

Four things are worth knowing about what goes in.

**The soil's hydraulic states are READ, not derived.** Saturation, field
capacity, the wilting point, the retention closure's exponent and the share of
each physical layer the soil column carries all come from
`pedology/scripts/land_column_properties.py`, which is the one place in this
pipeline that evaluates the closure. `soilinput.cpp` inverted the Cosby texture
regressions for itself and `vesperinput.cpp` rescaled the profile by regolith
depth for itself; WORLD-OF6N removed both, so a driver without these states
leaves every cell on its LPJ soil code instead.

**Land comes from `surface_class`, never `land_mask`.** The two disagree by 1.9%
of the planet, all of it dry closed-basin floor below sea level that `land_mask`
would flood. That is the whole point of the Orogen fork.

**Insolation is supplied as net downward surface shortwave.** `driver.cpp`'s
`NETSWRAD_TS` path then applies no albedo correction of its own, which is what we
want for the equilibrium evapotranspiration it also drives: this project computes
surface albedo from lithology and ExoPlaSim has already used it, so the number is
better than driver.cpp's global 0.17 constant. It is NOT the incident photon
supply to a canopy, and PCAR-2 and BIO-25 own that end; see
`biosphere/notes/ecological-forcing-field-contract.md`.

**The producer's intervals travel, and are not remapped onto a calendar.**
Interval 0 starts at day 0 of the simulation year and `build_vesper_header.py`
fits the declination phase on that assumption, so the two must be regenerated
together after any orbit change. Everything else about the partition is carried:
each interval's start, end and duration in absolute seconds and its local solar
phase go in the file, taken from the climatology's own time axis, and the
consumer integrates them onto its own 24-hour step. Up to V7 the file demanded
twelve bins that were the model's months, so this script had to remap the
producer's intervals onto them first, and the producer's own partition was gone
before the consumer saw it.

**One climatology or many.** Pass several to `--climatology` and each becomes one
year of forcing, in the order given; `vesperinput` cycles through them, so the
spin-up sees the whole sequence rather than one arbitrary phase of it. One file
is a fixed climate. Several are how a stellar cycle reaches the biosphere, since
a single repeating year cannot represent a variable star at all.

    python biosphere/scripts/build_lpj_driver.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import struct
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import netCDF4 as nc
import numpy as np
import yaml

from _paths import CONFIG, GENERATED, PROJECT_ROOT, climatology_path
from paths import rel  # noqa: E402

import builds
import orbit
from gridding import land_fraction_of_class
from orogen import Export
from rootable import read_rootable

# V2 regolith depth, V3 bedrock water, V4 multiple years, V5 nitrogen deposition
# read per Earth year rather than per orbit, V6 the fourth climate array named
# for the variable it is actually the range of, V7 the land column property
# contract's retention states and per-layer usable shares in place of regolith
# depth and bedrock fraction. V5 carries the same bytes as V4 and V6 the same
# bytes as V5; each of those differs only in what a field MEANS, which is
# exactly the change a magic has to catch, since such a file read by the wrong
# binary looks perfectly valid. V7 is the first that changes the bytes: two
# doubles per cell become four plus one per physical layer. V5's change is
# argued in biosphere/notes/time-base-unit-contract.md, V6's in
# biosphere/notes/ecological-forcing-field-contract.md and V7's in
# pedology/notes/land-column-property-contract.md.
#
# V8 IS CHRONOLOGICAL AND CARRIES ITS OWN INTERVALS. Every version up to V7 had
# a fixed twelve bins per year and the reader knew what a bin meant, so the
# cadence lived in the reader and the interval length lived nowhere. This
# builder had to remap the producer's own intervals onto the model's months to
# fit that, which threw the producer's interval identity away before the
# consumer ever saw it. V8 carries a table of intervals with explicit start,
# end and duration in absolute seconds and the local solar phase of each, ANY
# COUNT, and the fields per interval as intensive quantities. The 24-hour
# hydrology and biogeochemistry boundary is then the CONSUMER's: vesperinput
# integrates the intervals onto its own day rather than the format assuming
# that boundary on its behalf. Twelve intervals a year drop in, and so do one
# per day and one per timestep, with no format change and no reader change.
#
# What V8 gives up is the smooth daily curve interp_monthly_means_conserve
# manufactured between bin centres. That curve had no source: the producer
# states an interval mean and says nothing about the shape inside the interval,
# so the smooth reconstruction was invented structure. V8's consumer-side
# integration invents nothing and is the identity when an interval is one
# absolute day. EFOR-1's contract and EFOR-4 argue it;
# biosphere/config/ecological_forcing_contract.yaml is the declaration.
MAGIC = b"VESPDRV8"

# Samples inside each interval. Zero today: the subdaily arm is EFOR-3's to
# deliver and PCAR-1's and FIRE-1's to consume, and vesperinput refuses a
# non-zero count by name rather than reading a field it has no operator for.
# Carried in the header so that adding them later is a value change and not a
# format change.
SUBDAILY_SAMPLES = 0

# Seconds in the absolute day of biosphere/notes/time-base-unit-contract.md.
# Not this world's rotation, which is 30 hours; the interval bounds below are
# absolute seconds and the consumer's step is the absolute day.
DAY_SECONDS = 86400.0

# The file header after the magic, as one struct format shared by the writer and
# by --self-test, so the round-trip check cannot pass against a layout the
# writer does not use. Six int32s and two doubles: cells, intervals per year,
# year length in days, years of climate, subdaily samples per interval, a pad
# that keeps the doubles eight-byte aligned, then CO2 and nitrogen deposition.
# vesperinput.cpp reads them in this order.
HEADER_FORMAT = "<iiiiiidd"

# One interval record: start, end and duration in absolute seconds from the
# start of the simulation year, and the local solar phase at the start as a
# fraction of a rotation. Four packed doubles, which vesperinput reads as one
# block and static_asserts the struct size of.
INTERVAL_FORMAT = "<dddd"

# Coordinate precision shared with pedology/scripts/build_soil.py, so the soil
# map keys match exactly. See where lon_signed is rounded.
COORD_DECIMALS = 4
PROVENANCE_BYTES = 64
KELVIN = 273.15

# LPJ soil codes, from the table in modules/soilinput.h:
#   0 ice, 1 coarse, 2 medium, 3 fine, 4 medium-coarse, 5 fine-coarse,
#   6 fine-medium, 7 fine-medium-coarse, 8 organic, 9 vertisols
#
# Mapping Orogen's rock classes onto them is a judgement call, so it is written
# down rather than buried. The reasoning is weathering product, not parent rock
# hardness: what matters to a plant is the texture of the regolith the rock
# breaks down into and the water it can hold.
#
# Nothing maps to 8 (organic), because organic soils are a product of the
# biosphere we are about to model and cannot be an input to it. Nothing maps to
# 0 (ice) either; glaciated ground is handled by the climate model, not here.
SOIL_CODE_BY_ROCK = {
    "water":               2,  # unused, land only, but keep the table total
    "morb":                4,  # basalt weathers to clay-rich but stony ground
    "oib":                 4,
    "flood_basalt":        6,  # deep basalt saprolite, finer than arc material
    "arc_basalt":          4,
    "arc_andesite":        4,
    "rift_bimodal":        4,
    "granite":             1,  # granite grus is coarse and sandy
    "granodiorite":        1,
    "gneiss":              4,  # gneiss weathers coarse but with more clay
    "schist":              6,  # phyllosilicates give fine, platy regolith
    "quartzite":           1,  # almost pure quartz sand
    "melange":             7,  # tectonic mixture, so a mixed texture
    "shelf_clastic":       7,  # interbedded sandstone and shale
    "carbonate":           3,  # karst residuum is clay-rich terra rossa
    "foreland_clastic":    7,
    "continental_clastic": 5,
    "pelagic":             3,  # abyssal clay
    "evaporite":           9,  # vertisol: the shrink-swell salt-affected case
    "playa_clastic":       3,  # playa mud, fine and poorly drained
}


def project_relative(path: Path) -> str:
    """Relative to the project when it is inside it, absolute when it is not.

    Provenance should read cleanly for tracked inputs without crashing on a
    scratch path outside the tree, which is exactly what /tmp climatologies are
    during a test.
    """
    path = Path(path).resolve()
    return rel(path)


def month_lengths(year_length: int) -> np.ndarray:
    """The model's own month lengths, from the generator that compiles them in.

    Imported rather than recomputed. `VESPER_MONTH_LENGTHS` sizes the C++
    `Date` and is what the model's monthly output tables are binned on, so a
    second copy of the rule here would be a second thing to keep in step.
    Nothing is DELIVERED on these any more: the driver carries the producer's
    own intervals and the consumer integrates them onto the absolute day. They
    survive here for the month-weighted summary in the provenance report.
    """
    from build_vesper_header import month_lengths as _lengths
    return np.asarray(_lengths(year_length), dtype=float)


def producer_interval_days(times: np.ndarray, year_length: int) -> np.ndarray:
    """How many absolute days each of the producer's time bins actually spans.

    Derived from the climatology's own time axis through `lib/climatology.py`,
    which is the only place pyburn's binning arithmetic belongs. The bins do NOT
    all hold the same number of raw records: in the clean I/O regime the counts
    are [15]*5 + [16] + [15]*5 + [16], so two bins in twelve are longer than the
    other ten and neither is where a calendar would put a long month.

    The returned spans sum to `year_length`. That is a STRETCH, not an identity:
    at the clean write interval the records cover 99.56% of the orbit and the
    remainder is in no bin at all (`lib/climatology.py`). Distributing the
    uncovered fraction in proportion is the only choice available from the
    product alone, and the alternative -- letting the intervals fall short of the
    year -- would leave a gap the consumer's calendar has no way to represent.
    """
    from climatology import bin_weights
    weights = bin_weights(np.asarray(times, dtype=float))
    return weights * float(year_length)


def interval_table(spans: np.ndarray, rotation_hours: float) -> np.ndarray:
    """The interval record for one year: start, end, duration, solar phase.

    `spans` is each interval's length in absolute days and they tile the year.
    Bounds are absolute seconds from the start of the year, exact and never
    inferred from a record index: a consumer that infers an interval from a
    record number cannot tell a missing interval from a short one, and that
    distinction is what the acceptance check is.

    The phase is the local solar phase at the interval's START, as a fraction of
    a rotation. Rotation is uniform, so it is exact; the ecological forcing
    contract declares it derived for that reason rather than carried from the
    model. It advances every interval because the absolute day and this world's
    rotation are not the same length, and a consumer that assumes a fixed phase
    is assuming a 24-hour rotator.
    """
    spans = np.asarray(spans, dtype=float)
    if spans.ndim != 1 or spans.size < 1:
        raise SystemExit("an interval table needs at least one interval")
    if not np.all(spans > 0.0):
        raise SystemExit(f"every interval must have positive length, got {spans}")
    edges = np.concatenate([[0.0], np.cumsum(spans)]) * DAY_SECONDS
    rotation_seconds = float(rotation_hours) * 3600.0
    return np.stack([
        edges[:-1],
        edges[1:],
        np.diff(edges),
        np.mod(edges[:-1] / rotation_seconds, 1.0),
    ], axis=1)


def integrate_onto_days(values: np.ndarray, spans: np.ndarray,
                        year_length: int) -> np.ndarray:
    """What the consumer will compute, computed here so it can be checked.

    This is vesperinput's own reconstruction: each absolute day takes the
    duration-weighted mean of the intervals overlapping it. It is the identity
    when an interval is one absolute day, it invents no structure inside an
    interval, and for an intensive quantity it conserves the duration-weighted
    total exactly, which is identity 7 of the forcing contract.

    Reproduced rather than trusted. The consumer is C++ in another tree and
    cannot be executed here, so the arithmetic it will run is written out and
    the conservation identity asserted against the artifact this script emits.
    """
    values = np.asarray(values, dtype=float)
    spans = np.asarray(spans, dtype=float)
    days = np.full(int(year_length), 1.0)
    if abs(spans.sum() - days.sum()) > 1e-6 * days.sum():
        raise SystemExit(
            f"the producer's intervals span {spans.sum():.6f} days and the "
            f"model's year {days.sum():.6f}; they must tile the same year")

    p_edges = np.concatenate([[0.0], np.cumsum(spans)])
    d_edges = np.concatenate([[0.0], np.cumsum(days)])
    overlap = np.clip(
        np.minimum(d_edges[1:, None], p_edges[None, 1:])
        - np.maximum(d_edges[:-1, None], p_edges[None, :-1]), 0.0, None)

    out = np.tensordot(overlap, values, axes=(1, 0)) / days.reshape(
        (-1,) + (1,) * (values.ndim - 1))

    total_in = float(np.tensordot(spans, values, axes=(0, 0)).sum())
    total_out = float(np.tensordot(days, out, axes=(0, 0)).sum())
    scale = max(abs(total_in), 1.0)
    if abs(total_out - total_in) > 1e-9 * scale:
        raise SystemExit(
            f"the consumer's day integration would lose mass: {total_in!r} in, "
            f"{total_out!r} out")
    return out


def area_weights(lat: np.ndarray, nlon: int) -> np.ndarray:
    w = np.cos(np.deg2rad(lat))
    return (w / w.sum())[:, None] * np.ones((1, nlon)) / nlon


def soil_codes(config: dict, land: np.ndarray) -> tuple[np.ndarray, dict]:
    """Dominant soil code per grid cell, integrated from the native mesh.

    Not sampled from the gridded export. That export resamples categorical fields
    by the region containing the cell centre, which at T42 throws away most of
    what is in a cell. Instead each soil code's share of a cell's *land* area is
    accumulated over the configured build's native mesh (the 10M export for new
    spatial-support work) and the largest share wins.

    Dominant rather than averaged, because a soil code is categorical: the mean
    of code 1 and code 3 is code 2, which is not what half granite and half
    carbonate behaves like.
    """
    mesh = Export(builds.mesh_export(config))
    grid_dir = builds.grid_export(config, str(config["model"]["resolution"]).upper())
    rock = mesh.substrate_class
    classes = {c["id"]: c["code"] for c in mesh.manifest["lithology"]["rockClasses"]}

    unknown = sorted({classes[i] for i in np.unique(rock) if i in classes}
                     - set(SOIL_CODE_BY_ROCK))
    if unknown:
        raise SystemExit(
            f"rock classes with no soil-code mapping: {unknown}. Add them to "
            f"SOIL_CODE_BY_ROCK and record why."
        )

    # Accumulate land-area share per soil code, then take the argmax per cell.
    shares: dict[int, np.ndarray] = {}
    for rock_id, code in classes.items():
        soil = SOIL_CODE_BY_ROCK[code]
        selected = rock == rock_id
        if not selected.any():
            continue
        share = land_fraction_of_class(mesh, grid_dir, selected)
        shares[soil] = shares.get(soil, np.zeros_like(share)) + share

    codes_present = sorted(shares)
    stacked = np.stack([shares[c] for c in codes_present])
    dominant = np.array(codes_present, dtype=np.int32)[np.argmax(stacked, axis=0)]

    # Cells with no land in the mesh but flagged land by the climatology's mask
    # would otherwise take whichever code sorts first. Medium is the neutral
    # default and the count is reported rather than hidden.
    empty = stacked.sum(axis=0) <= 0.0
    dominant[empty] = 2
    stranded = int(np.sum(empty & land))

    total = max(int(land.sum()), 1)
    summary = {
        "method": "dominant land-area share, integrated over the native mesh",
        "cells_with_no_mesh_land_defaulted_to_medium": stranded,
        "soil_code_share_of_land_cells": {
            int(c): round(float(np.sum(dominant[land] == c) / total), 4)
            for c in sorted(np.unique(dominant[land]))
        },
    }
    return dominant, summary


def self_test() -> int:
    """The interval arithmetic and the header layout, on fixtures.

    WHAT THIS COVERS AND WHAT IT DOES NOT. `integrate_onto_days` is the
    arithmetic `vesperinput.cpp:integrate_year` runs, written out here so it can
    be executed: LPJ-GUESS does not build on this tree, so the C++ is checked by
    `g++ -fsyntax-only` against a synthetic header and never run. What is
    executed here is therefore the OPERATOR and the FORMAT, not the reader. A
    fixture that does not get the verdict it was built for is a defect in this
    checker rather than in the builder.

        python biosphere/scripts/build_lpj_driver.py --self-test

    No climatology, no soil map, no model.
    """
    failures = 0

    def check(label: str, ok: bool, detail: str = "") -> None:
        nonlocal failures
        failures += 0 if ok else 1
        print(f"[{'  ok  ' if ok else ' FAIL '}] {label}"
              + (f"\n           {detail}" if detail and not ok else ""))

    def refuses(label: str, fn) -> None:
        try:
            fn()
        except SystemExit as error:
            check(label, True, str(error))
            return
        check(label, False, "it was accepted")

    year = 12
    rotation = 30.0

    # The identity case, and it is the one the whole format change is for: one
    # interval per absolute day passes through untouched, so a chronological
    # forcing is not reconstructed at all.
    spans = np.ones(year)
    values = np.arange(float(year))
    days = integrate_onto_days(values, spans, year)
    check("one interval per absolute day integrates to itself",
          bool(np.allclose(days, values)), f"{days} against {values}")

    # Uneven intervals, which is what the producer actually emits. The
    # duration-weighted total is what has to survive; the daily shape is a step
    # function and is not claimed to be anything else.
    # Fractional edges deliberately, because that is what the producer emits:
    # its bins are 15.08 days and change and no edge lands on a day boundary.
    spans = np.array([2.5, 5.0, 0.5, 4.0])
    values = np.array([10.0, -3.0, 7.0, 0.5])
    days = integrate_onto_days(values, spans, year)
    check("uneven intervals conserve the duration-weighted total",
          abs(float(days.sum()) - float((spans * values).sum())) < 1e-9,
          f"{days.sum()} against {(spans * values).sum()}")
    check("a day wholly inside one interval takes that interval's value",
          abs(days[3] - values[1]) < 1e-12, f"{days[3]} against {values[1]}")
    # Day 2 spans 2.0 to 3.0 and the first edge is at 2.5, so it is half in
    # interval 0 and half in interval 1. Half of 10 plus half of -3 is 3.5, and
    # a reader that snapped the edge to a day boundary would return 10 or -3.
    check("a day straddling two intervals takes their weighted mean",
          abs(days[2] - 3.5) < 1e-12, f"day 2 {days[2]}, expected 3.5")
    # Day 7 spans 7.0 to 8.0, and 7.5 is where interval 2 begins, so it is half
    # of -3 and half of 7.
    check("a day straddling the shortest interval takes its share too",
          abs(days[7] - 2.0) < 1e-12, f"day 7 {days[7]}, expected 2.0")

    # A field with a grid axis, which is the shape the builder actually passes.
    field = np.stack([values, values * 2.0], axis=-1)[:, None, :]
    days = integrate_onto_days(field, spans, year)
    check("the operator carries the grid axes through",
          days.shape == (year, 1, 2)
          and abs(float(days[:, 0, 1].sum()) - float((spans * values * 2.0).sum())) < 1e-9,
          f"shape {days.shape}")

    refuses("intervals that do not tile the year are refused",
            lambda: integrate_onto_days(np.array([1.0, 2.0]),
                                        np.array([1.0, 2.0]), year))

    table = interval_table(np.array([2.0, 5.0, 1.0, 4.0]), rotation)
    check("the interval table starts at zero and reaches the year",
          abs(table[0, 0]) < 1e-9
          and abs(table[-1, 1] - year * DAY_SECONDS) < 1e-6,
          f"{table[0, 0]} to {table[-1, 1]}")
    check("every interval begins where the previous one ends",
          bool(np.allclose(table[1:, 0], table[:-1, 1])))
    check("the declared duration is the bounds' own difference",
          bool(np.allclose(table[:, 2], table[:, 1] - table[:, 0])))
    check("the solar phase is in [0, 1) and wraps rather than accumulating",
          bool(np.all(table[:, 3] >= 0.0) and np.all(table[:, 3] < 1.0)),
          f"{table[:, 3]}")
    # 2 absolute days is 48 hours, which is 1.6 rotations of 30 hours, so the
    # phase at the second interval's start is 0.6 and not 1.6. A phase that
    # accumulated would pass every other check here.
    check("the phase after two absolute days is 1.6 rotations, wrapped to 0.6",
          abs(table[1, 3] - 0.6) < 1e-9, f"{table[1, 3]}")

    refuses("an interval of no length is refused",
            lambda: interval_table(np.array([2.0, 0.0, 1.0]), rotation))
    refuses("a negative interval is refused",
            lambda: interval_table(np.array([2.0, -1.0, 1.0]), rotation))

    packed = struct.pack(HEADER_FORMAT, 7, 12, 183, 3, 0, 0, 285.0, 0.5)
    check("the header packs to the size vesperinput reads",
          struct.calcsize(HEADER_FORMAT) == 6 * 4 + 2 * 8,
          f"{struct.calcsize(HEADER_FORMAT)} bytes")
    check("the header round-trips through its own format",
          struct.unpack(HEADER_FORMAT, packed) == (7, 12, 183, 3, 0, 0, 285.0, 0.5))
    check("an interval record is four packed doubles",
          struct.calcsize(INTERVAL_FORMAT) == 4 * 8,
          f"{struct.calcsize(INTERVAL_FORMAT)} bytes")
    check("the interval table's own bytes round-trip",
          struct.unpack(INTERVAL_FORMAT, table[1].astype("<f8").tobytes())
          == tuple(float(v) for v in table[1]))

    print(f"\n{failures} failed" if failures else "\nno failures")
    return failures


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--climatology", type=Path, nargs="+", default=None,
                        help="one or more climatologies, each one year of "
                             "forcing, cycled in the order given")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--states", type=Path, default=None,
                        help="the land column property contract's emitted "
                             "per-cell states. Without them every cell falls "
                             "back to its LPJ soil code, which carries its own "
                             "tabulated capacities.")
    parser.add_argument("--rootable", type=Path, default=None,
                        help="BIO-11 rootable-fraction artifact; default the "
                             "active build and rung")
    parser.add_argument("--ndep", type=float, default=0.5,
                        help="nitrogen deposition, kgN/ha per EARTH YEAR. "
                             "Absolute time, not per orbit: deposition is an "
                             "atmospheric flux and does not know how long this "
                             "world takes to go round its star. A declared "
                             "assumption -- this world has no deposition field "
                             "and no industry. Default is a low "
                             "pre-industrial-like value; report the "
                             "sensitivity, do not tune it.")
    parser.add_argument("--self-test", action="store_true",
                        help="run the interval arithmetic and header layout "
                             "fixtures and exit. No climatology, no model.")
    args = parser.parse_args()

    if args.self_test:
        raise SystemExit(1 if self_test() else 0)

    config = yaml.safe_load(CONFIG.read_text())
    climatologies = list(args.climatology) if args.climatology else [climatology_path()]
    # Each climatology becomes one year of forcing, so they must all describe
    # the same world as the soil map and the compiled header.
    import sys as _sys
    _sys.path.insert(0, str(PROJECT_ROOT / "lib"))
    from provenance import require_build
    for _c in climatologies:
        require_build(Path(_c), "climatology", config)
    for path in climatologies:
        if not Path(path).is_file():
            raise SystemExit(f"{path} does not exist")
    climatology = Path(climatologies[0])

    year_length = orbit.model_year_days(config)
    co2_ppm = float(config["atmosphere"]["pCO2_bar"]) / 1.0 * 1e6

    # The model's months. Nothing is remapped onto them any more -- the driver
    # carries the producer's own intervals -- and they survive here for the
    # month-weighted summary in the report, which is a human-readable digest and
    # not a partition anything is delivered on.
    months = month_lengths(year_length)

    # This world's rotation, for the local solar phase of each interval. The
    # phase is not a property of the interval bounds alone: it needs the
    # rotation, and the rotation is not the absolute day.
    rotation_hours = float(config["planet"]["rotation_hours"])

    # Each climatology contributes one year, stacked as [year][interval][lat][lon].
    tas_y, pr_y, rss_y, tsrange_y = [], [], [], []
    spans_y = []
    tables = []
    lat = lon = lsm = None
    for path in climatologies:
        with nc.Dataset(path) as data:
            this_lat = np.asarray(data["lat"][:], dtype=float)
            this_lon = np.asarray(data["lon"][:], dtype=float)
            if lat is None:
                lat, lon = this_lat, this_lon
                lsm = np.asarray(data["lsm"][0], dtype=float)
            elif not (np.allclose(lat, this_lat) and np.allclose(lon, this_lon)):
                raise SystemExit(
                    f"{path} is on a different grid from {climatologies[0]}; "
                    f"every year has to share one grid")
            spans = producer_interval_days(
                np.asarray(data["time"][:], dtype=float), year_length)
            spans_y.append(spans.tolist())

            tables.append(interval_table(spans, rotation_hours))

            # EVERY FIELD IS INTENSIVE, and the conversions happen here because
            # the producer converts into the contract unit and the consumer out
            # of it. Kelvin to C, and m/s to mm per absolute day, which is a
            # RATE and not a per-interval total: a total would need an interval
            # length to be chosen to form it over, and that choice is what put
            # 18 days of a rate measured over 16.09 into the last bin of V7.
            # 86400 is the absolute day of
            # biosphere/notes/time-base-unit-contract.md, not this world's
            # rotation.
            tas_y.append(np.asarray(data["tas"][:], dtype=float) - KELVIN)
            pr_y.append(np.asarray(data["pr"][:], dtype=float) * 1000.0 * 86400.0)
            rss_y.append(np.asarray(data["rss"][:], dtype=float))
            # The range of the SURFACE temperature. maxt and mint are extrema of
            # dt(:,NLEP) and bracket ts, not tas. The near-surface AIR extrema
            # are codes 201/202 and reach a product as tasmax/tasmin since
            # world-j0az; carrying them needs a driver array VESPDRV does not
            # have, which is why vesperinput still does not hand anything to
            # climate.dtr, whose one reader means an air-temperature range.
            # See the field contract.
            tsrange_y.append(np.maximum(
                np.asarray(data["maxt"][:], dtype=float)
                - np.asarray(data["mint"][:], dtype=float), 0.0))
    tas = np.stack(tas_y)   # [year][interval][lat][lon]
    pr = np.stack(pr_y)
    rss = np.stack(rss_y)
    tsrange = np.stack(tsrange_y)
    nyears = tas.shape[0]

    nintervals = tas.shape[1]
    if any(table.shape[0] != nintervals for table in tables):
        raise SystemExit(
            "the climatologies do not agree on how many intervals a year has: "
            f"{[int(table.shape[0]) for table in tables]}. One driver file "
            "carries one interval count, because a consumer indexes the years "
            "of a cycle through one array.")
    intervals = np.stack(tables)   # [year][interval][start, end, duration, phase]

    # WHAT THE CONSUMER WILL DO, DONE HERE SO IT CAN FAIL HERE. vesperinput
    # integrates the intervals onto its own absolute day; the same arithmetic
    # runs over every field below and raises if it would not conserve the
    # duration-weighted total. Identity 7 of the forcing contract, checked
    # against the artifact this script is about to write rather than asserted
    # about it.
    for year in range(nyears):
        for field in (tas[year], pr[year], rss[year], tsrange[year]):
            integrate_onto_days(field, np.asarray(spans_y[year]), year_length)

    land = lsm > 0.5
    rootable, rootable_provenance = read_rootable(
        config, lat, lon, args.rootable, land=land)
    # THE SIMULATED SET IS THE ROOTABLE GROUND, NOT THE OWNERSHIP MASK. `land`
    # is ExoPlaSim's binary 0.5 coastline rounding, and intersecting with it
    # dropped every cell the rounding gave to the ocean -- exactly the cells
    # SPAT-5's tile model exists for. BIO-11's f_rootable is integrated from
    # the native mesh and is positive only where there is rootable land, so it
    # already implies a positive land fraction and needs no mask beside it:
    # measured at T21 on canonical-10m-carve2, all 1617 cells with
    # f_rootable > 0 have code-1720 land fraction > 0, and the 620 this adds to
    # the previous 997 are precisely the sea-owned partial cells.
    #
    # `land` is kept for the support check above and for the report below,
    # where it still says what the ownership mask thinks.
    simulated_land = rootable > 0.0
    codes, soil_summary = soil_codes(config, simulated_land)

    # ExoPlaSim's longitudes run 0..360; LPJ-GUESS expects -180..180.
    #
    # Rounded to COORD_DECIMALS because LPJ-GUESS's soil map lookup keys on
    # std::pair<double,double> and compares it exactly. pedology/build_soil.py
    # writes its coordinates to the same precision, so the two round-trip to
    # identical doubles and the lookup needs no search radius. That agreement is
    # a contract between the two scripts, not a coincidence.
    lon_signed = np.round(np.where(lon > 180.0, lon - 360.0, lon), COORD_DECIMALS)
    lat = np.round(lat, COORD_DECIMALS)

    provenance = (f"{config.get('source_build')}|{climatology.parent.name}|"
                  f"root:{rootable_provenance['sha256'][:12]}").encode()
    provenance = provenance[:PROVENANCE_BYTES].ljust(PROVENANCE_BYTES, b"\0")

    output = args.output or (GENERATED / "vesper_driver.bin")
    output.parent.mkdir(parents=True, exist_ok=True)

    # THE SIMULATED SOIL'S RETENTION STATES, READ AND NOT DERIVED. The land
    # column property contract is the one description of this world's soil
    # hydraulics: `pedology/scripts/land_column_properties.py` evaluates the
    # named closure at this planet's gravity and emits, per cell, saturation,
    # field capacity, the wilting point, the closure's exponent and the share
    # of each physical layer the soil column carries. `soilinput.cpp` derived
    # its own from the Cosby regressions until WORLD-OF6N; it now refuses
    # without these.
    #
    # Keyed on the same rounded coordinates the soil map uses. Absent, the
    # sentinel goes in and `vesperinput` runs the driver's LPJ soil code path,
    # which is a texture class with its own tabulated capacities and no column
    # geometry.
    # The physical layer count is the contract's, declared once. `vesperinput`
    # checks what arrives against its compiled NSOILLAYER and refuses a
    # mismatch, so a driver built against another profile fails by name rather
    # than reading one array as another.
    contract = yaml.safe_load(
        (PROJECT_ROOT / "pedology" / "config" / "land_column_properties.yaml")
        .read_text(encoding="utf-8"))
    n_layers = int(contract["geometry"]["physical_layer_count"])

    states_by_coord: dict[tuple[float, float], tuple[float, ...]] = {}
    states_path = args.states
    if states_path is None:
        candidate = builds.land_column_states()
        states_path = candidate if candidate.is_file() else None
    if states_path is not None:
        lines = Path(states_path).read_text().splitlines()
        header = lines[0].split()
        usable_names = sorted(name for name in header
                              if name.startswith("u") and name[1:].isdigit())
        for needed in ("b", "theta_s", "theta_fc", "theta_wp"):
            if needed not in header:
                raise SystemExit(
                    f"{states_path} has no {needed} column; write it with "
                    "pedology/scripts/land_column_properties.py")
        if len(usable_names) != n_layers:
            raise SystemExit(
                f"{states_path} carries {len(usable_names)} per-layer usable "
                f"shares and the contract declares {n_layers} physical layers. "
                "They are the contract's weathered-bedrock rule and vesperinput "
                "applies them rather than deriving them; rewrite the states "
                "with pedology/scripts/land_column_properties.py.")
        index = [header.index(name) for name in
                 ("b", "theta_s", "theta_fc", "theta_wp")]
        index += [header.index(name) for name in usable_names]
        for line in lines[1:]:
            parts = line.split()
            key = (round(float(parts[0]), COORD_DECIMALS),
                   round(float(parts[1]), COORD_DECIMALS))
            states_by_coord[key] = tuple(float(parts[i]) for i in index)
        print(f"land column states from {Path(states_path).name}: "
              f"{len(states_by_coord)} cells, {n_layers} physical layers")
    else:
        print("no land column states found; every cell falls back to its LPJ "
              "soil code")

    rows = np.argwhere(simulated_land)
    missing_states = 0
    with output.open("wb") as handle:
        handle.write(MAGIC)
        # Six int32s and not four: the interval count replaces the bin count and
        # says nothing about what an interval is, the subdaily count is carried
        # so that adding samples later is a value change and not a format
        # change, and the pad keeps the doubles that follow eight-byte aligned.
        handle.write(struct.pack(HEADER_FORMAT, len(rows), nintervals,
                                 year_length, nyears, SUBDAILY_SAMPLES, 0,
                                 co2_ppm, args.ndep))
        handle.write(provenance)
        # THE INTERVAL TABLE, ONCE PER YEAR AND BEFORE ANY CELL. Every cell
        # shares one time axis, so carrying it per cell would be the same fact
        # written a hundred thousand times with nothing checking the copies
        # against each other.
        handle.write(intervals.astype("<f8").tobytes())
        for j, i in rows:
            handle.write(struct.pack("<dd", float(lon_signed[i]), float(lat[j])))
            handle.write(struct.pack("<ii", int(codes[j, i]), 0))
            key = (float(lon_signed[i]), float(lat[j]))
            states = states_by_coord.get(key)
            if states is None:
                # THE SENTINEL, and it is a saturation of zero. `vesperinput`
                # takes a non-positive saturation as "no contract states" and
                # runs the LPJ soil code path for that cell, which carries its
                # own tabulated capacities. A plausible-looking default here
                # would be a third derivation of the soil.
                states = (0.0, 0.0, 0.0, 0.0) + (1.0,) * n_layers
                missing_states += 1
            handle.write(struct.pack(f"<{4 + n_layers}d", *states))
            # Flattened [year][interval], matching what vesperinput indexes.
            # Precipitation goes in as the mean RATE in mm per absolute day, not
            # as a per-interval total: the contract carries no extensive field,
            # and the depth is the rate times a duration the interval table
            # already states.
            handle.write(tas[:, :, j, i].astype("<f8").tobytes())
            handle.write(pr[:, :, j, i].astype("<f8").tobytes())
            handle.write(rss[:, :, j, i].astype("<f8").tobytes())
            handle.write(tsrange[:, :, j, i].astype("<f8").tobytes())

    weights = area_weights(lat, len(lon))
    effective_weights = weights * rootable
    lw = effective_weights[simulated_land]

    def interval_mean(field: np.ndarray) -> np.ndarray:
        """Mean over years and intervals, weighting each interval by its length.

        Duration-weighted and not a plain mean over the interval axis: the
        producer's intervals are not equal, and averaging them flat over-weights
        the short ones. It is the one operator for every field here because
        every field here is intensive.
        """
        spans = np.asarray(spans_y, dtype=float)   # [year][interval]
        weighted = np.einsum("yi,yi...->...", spans, field)
        return weighted / spans.sum()

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "climatologies": [project_relative(p) for p in climatologies],
        "climatology_sha256": [hashlib.sha256(Path(p).read_bytes()).hexdigest()
                               for p in climatologies],
        "years_of_climate": nyears,
        "years_note": ("Cycled by vesperinput, so spin-up sees the whole "
                       "sequence. One year is a fixed climate; several are how a "
                       "stellar cycle reaches the biosphere."),
        "config_sha256": hashlib.sha256(CONFIG.read_bytes()).hexdigest(),
        # The parsed config, so that a later "is this artifact still current?"
        # can be answered per KEY. `config_sha256` alone cannot: it moves for an
        # edited comment exactly as it does for an edited parameter, and
        # `check_consistency.py` reported these artifacts stale on a comment
        # change until it had this to read. Same field, same purpose as
        # `source_config` in an ExoPlaSim run manifest.
        "source_config": config,
        "source_build": config.get("source_build"),
        "stellar_spectrum": config.get("radiation", {}).get("stellar_spectrum"),
        "flux_earth": float(config["orbit"]["baseline_flux_earth"]),
        "orbital_year_earth_days": orbit.orbital_year_days(config),
        "year_length_days": year_length,
        "driver_format": MAGIC.decode(),
        "intervals_per_year": nintervals,
        "subdaily_samples_per_interval": SUBDAILY_SAMPLES,
        # The producer's OWN partition, per climatology, from the climatology's
        # own time axis. It is what the driver carries: no remap happens any
        # more, and there is no second partition for it to disagree with.
        "producer_interval_days": spans_y,
        "model_month_days": months.astype(int).tolist(),
        "model_months_note": (
            "The model's months are still the calendar Date sizes with, and "
            "they are no longer a partition anything is delivered on. Up to V7 "
            "the driver remapped the producer's intervals onto them because the "
            "format demanded twelve bins; V8 carries the producer's intervals "
            "with explicit bounds and vesperinput integrates them onto its own "
            "absolute day."),
        "consumer_reconstruction": (
            "duration-weighted integration onto the absolute day, in "
            "vesperinput. Exact when an interval is one absolute day, conserves "
            "the duration-weighted total for any interval structure, and "
            "invents no shape inside an interval. The smooth curve "
            "interp_monthly_means_conserve manufactured between bin centres is "
            "gone, and it had no source: the producer states an interval mean "
            "and says nothing about the shape inside the interval."),
        "land_cells": int(land.sum()),
        "simulated_rootable_cells": int(len(rows)),
        "fully_nonrootable_land_cells_omitted": int(np.sum(land & ~simulated_land)),
        # The cells the binary coastline rounding gives to the ocean and which
        # carry rootable ground anyway. They were dropped until SPAT-5's tile
        # model needed a land tile on them; the count is the size of what the
        # ownership mask was hiding from the biosphere.
        "simulated_cells_outside_the_ownership_mask": int(
            np.sum(simulated_land & ~land)),
        "rootable_surface": rootable_provenance,
        "rootable_area_fraction_of_model_land": float(
            effective_weights[land].sum() / weights[land].sum()),
        "rootable_area_fraction_of_simulated_ground": float(
            effective_weights[simulated_land].sum()
            / weights[simulated_land].sum()),
        "land_column_states_source": (project_relative(Path(states_path))
                                      if states_path
                                      else "none, LPJ soil codes assumed"),
        "land_column_states_sha256": (
            hashlib.sha256(Path(states_path).read_bytes()).hexdigest()
            if states_path else None),
        "physical_layers": n_layers,
        "cells_without_land_column_states": missing_states,
        "soil_hydraulics": (
            "read from the land column property contract's emitted states and "
            "derived nowhere in this pipeline but "
            "pedology/scripts/land_column_properties.py. soilinput.cpp's Cosby "
            "inversion and vesperinput.cpp's regolith rescaling were removed "
            "under WORLD-OF6N."),
        "co2_ppm": co2_ppm,
        "ndep_kgn_ha_earth_year": args.ndep,
        "ndep_kgn_ha_absolute_day": args.ndep / orbit.EARTH_SIDEREAL_YEAR_DAYS,
        "ndep_note": (
            "Declared, not measured. This world has no deposition field and no "
            "industry, so any value is an assumption and NPP inherits it. The "
            "unit is per EARTH year; vesperinput divides by the Earth year to "
            "reach the per-absolute-day rate it hands the model."),
        "insolation": "NETSWRAD_TS, net downward surface shortwave (rss), W/m2",
        "fourth_array": (
            "maxt - mint, the range of the SURFACE temperature. NOT the "
            "near-surface air diurnal range: maxt/mint are extrema of "
            "dt(:,NLEP) and bracket ts, not tas. vesperinput does not assign it "
            "to climate.dtr. The air range is tasmax - tasmin, codes 201/202, "
            "which this format has no array for. "
            "biosphere/notes/ecological-forcing-field-contract.md"),
        "land_definition": "lsm from the climatology, itself built from surface_class",
        "simulation_population": (
            "BIO-11 f_rootable > 0, which is integrated from the native mesh "
            "and is NOT intersected with the binary ownership mask: the "
            "coastline rounding is not a statement about where ground is. LPJ "
            "outputs remain intensive over rootable ground and every extensive "
            "consumer applies the same fraction"),
        # Month-weighted, because the months are NOT equal: a plain mean over the
        # twelve would over-weight the eleven short ones.
        "land_mean_temperature_c": float(
            np.average(interval_mean(tas)[simulated_land], weights=lw)),
        "land_mean_precip_mm_per_earth_year": float(
            np.average(interval_mean(pr)[simulated_land], weights=lw)
            * orbit.EARTH_CALENDAR_YEAR_DAYS),
        "land_mean_net_sw_w_m2": float(
            np.average(interval_mean(rss)[simulated_land], weights=lw)),
        "soil": soil_summary,
        "soil_code_mapping": SOIL_CODE_BY_ROCK,
        "output_sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            cwd=PROJECT_ROOT).stdout.strip() or None,
    }
    # Name the generator in the artifact. Without it, a staleness check can
    # say WHICH file is out of date but not what to re-run, and it named the
    # wrong script for two of three artifacts for exactly that reason.
    report["generator"] = "biosphere/scripts/build_lpj_driver.py"
    report_path = output.with_name(output.stem + "_provenance.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")

    print(f"land cells     {int(land.sum())} by the ownership mask; "
          f"{len(rows)} simulated, of which "
          f"{int(np.sum(simulated_land & ~land))} the mask calls ocean")
    print(f"year length    {year_length} absolute days")
    print(f"intervals      {nintervals} per year, carried with explicit bounds: "
          + ", ".join(f"{s:.2f}" for s in spans_y[0]) + " days")
    print(f"               solar phase at each start "
          + ", ".join(f"{v:.3f}" for v in intervals[0, :, 3]))
    print(f"subdaily       {SUBDAILY_SAMPLES} samples per interval")
    print(f"climate years  {nyears} "
          f"({'cycled' if nyears > 1 else 'fixed climate, repeated'})")
    print(f"CO2            {co2_ppm:.0f} ppm     "
          f"N deposition {args.ndep} kgN/ha/Earth-yr "
          f"({args.ndep / orbit.EARTH_SIDEREAL_YEAR_DAYS:.3e} per absolute day)")
    print(f"land means     {report['land_mean_temperature_c']:.2f} C, "
          f"{report['land_mean_precip_mm_per_earth_year']:.0f} mm/Earth-yr, "
          f"{report['land_mean_net_sw_w_m2']:.1f} W/m2 net SW")
    print("soil codes     " + ", ".join(
        f"{k}:{v:.0%}" for k, v in soil_summary["soil_code_share_of_land_cells"].items()))
    print(f"\nwrote {project_relative(output)} "
          f"({output.stat().st_size / 1e6:.1f} MB)")
    print(f"      {report_path.name}")


if __name__ == "__main__":
    main()
