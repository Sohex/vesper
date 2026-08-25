"""Build the binary driver file LPJ-GUESS reads, from climate and lithology.

One self-describing file carrying the gridlist, a soil code per cell derived from
World Orogen lithology, and binned climate. `vesperinput.cpp` reads it.

Three things are worth knowing about what goes in.

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

**Bin order defines the calendar, and bin LENGTH does not.** Bin 0 starts at day
0 of the simulation year and `build_vesper_header.py` fits the declination phase
on that assumption, so the two must be regenerated together after any orbit
change. But a pyburn bin and a model month are two different partitions of the
same year, and the producer owns its own interval bounds. `interval_to_month`
below remaps between them conservatively rather than assuming they agree.

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

# V2 regolith depth, V3 bedrock water, V4 multiple years, V5 nitrogen deposition
# read per Earth year rather than per orbit, V6 the fourth climate array named
# for the variable it is actually the range of. V5 carries the same bytes as V4
# and V6 the same bytes as V5; each differs only in what a field MEANS, which is
# exactly the change a magic has to catch, since such a file read by the wrong
# binary looks perfectly valid. V5's change is argued in
# biosphere/notes/time-base-unit-contract.md and V6's in
# biosphere/notes/ecological-forcing-field-contract.md.
MAGIC = b"VESPDRV6"

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
    `Date`, and `interp_monthly_means_conserve` spreads a bin over exactly these
    days, so a second copy of the rule here would be a second thing to keep in
    step.
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


def interval_to_month(values: np.ndarray, spans: np.ndarray,
                      months: np.ndarray) -> np.ndarray:
    """Remap per-interval values onto the model's months, conserving the total.

    `values` is per-interval and INTENSIVE in time: a mean over the interval, or
    a rate per absolute day. Its leading axis is the interval axis. Both
    partitions tile the same year starting at day 0.

    One operator serves means and rates alike, which is the point of expressing
    precipitation as a rate before it gets here. For a mean the result is the
    duration-weighted mean over the month; for a rate it is the mean rate over
    the month, and multiplying it by the month length gives a total that
    conserves exactly:

        sum_k months[k] * out[k] = sum_b spans[b] * values[b]

    That identity is asserted below, because a remap that silently loses mass is
    worse than no remap: the annual total would still look plausible.
    """
    values = np.asarray(values, dtype=float)
    spans = np.asarray(spans, dtype=float)
    months = np.asarray(months, dtype=float)
    if abs(spans.sum() - months.sum()) > 1e-6 * months.sum():
        raise SystemExit(
            f"the producer's intervals span {spans.sum():.6f} days and the "
            f"model's year {months.sum():.6f}; they must tile the same year")

    p_edges = np.concatenate([[0.0], np.cumsum(spans)])
    m_edges = np.concatenate([[0.0], np.cumsum(months)])
    # Overlap in days between month k and interval b.
    overlap = np.clip(
        np.minimum(m_edges[1:, None], p_edges[None, 1:])
        - np.maximum(m_edges[:-1, None], p_edges[None, :-1]), 0.0, None)

    out = np.tensordot(overlap, values, axes=(1, 0)) / months.reshape(
        (-1,) + (1,) * (values.ndim - 1))

    total_in = float(np.tensordot(spans, values, axes=(0, 0)).sum())
    total_out = float(np.tensordot(months, out, axes=(0, 0)).sum())
    scale = max(abs(total_in), 1.0)
    if abs(total_out - total_in) > 1e-9 * scale:
        raise SystemExit(
            f"interval-to-month remap lost mass: {total_in!r} in, "
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--climatology", type=Path, nargs="+", default=None,
                        help="one or more climatologies, each one year of "
                             "forcing, cycled in the order given")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--soil-map", type=Path, default=None,
                        help="pedology soilmap.txt, for its regolith depth "
                             "column. Without it every cell is given the full "
                             "profile depth, which is LPJ-GUESS's own default.")
    parser.add_argument("--ndep", type=float, default=0.5,
                        help="nitrogen deposition, kgN/ha per EARTH YEAR. "
                             "Absolute time, not per orbit: deposition is an "
                             "atmospheric flux and does not know how long this "
                             "world takes to go round its star. A declared "
                             "assumption -- this world has no deposition field "
                             "and no industry. Default is a low "
                             "pre-industrial-like value; report the "
                             "sensitivity, do not tune it.")
    args = parser.parse_args()

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

    # The model's months, which are what interp_monthly_*_conserve spreads a
    # value over. Every field below is remapped onto these from the producer's
    # own intervals rather than being assumed to already be on them.
    months = month_lengths(year_length)

    # Each climatology contributes one year, stacked as [year][month][lat][lon].
    tas_y, pr_y, rss_y, tsrange_y = [], [], [], []
    spans_y = []
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

            def onto_months(values: np.ndarray) -> np.ndarray:
                return interval_to_month(values, spans, months)

            # Kelvin to C and m/s to mm per absolute day BEFORE the remap, so the
            # operator sees one intensive quantity per field and the conserved
            # total is the one that means something. 86400 is the absolute day of
            # biosphere/notes/time-base-unit-contract.md, not Vesper's rotation.
            tas_y.append(onto_months(
                np.asarray(data["tas"][:], dtype=float) - KELVIN))
            pr_y.append(onto_months(
                np.asarray(data["pr"][:], dtype=float) * 1000.0 * 86400.0))
            rss_y.append(onto_months(np.asarray(data["rss"][:], dtype=float)))
            # The range of the SURFACE temperature. maxt and mint are extrema of
            # dt(:,NLEP) and bracket ts, not tas; the near-surface AIR extrema are
            # codes 201/202, which no product carries yet. See the field
            # contract, and note that vesperinput does not hand this to
            # climate.dtr, whose one reader means an air-temperature range.
            tsrange_y.append(onto_months(np.maximum(
                np.asarray(data["maxt"][:], dtype=float)
                - np.asarray(data["mint"][:], dtype=float), 0.0)))
    tas = np.stack(tas_y)   # [year][month][lat][lon]
    pr = np.stack(pr_y)
    rss = np.stack(rss_y)
    tsrange = np.stack(tsrange_y)
    nyears = tas.shape[0]

    nbins = tas.shape[1]
    if nbins != 12:
        raise SystemExit(f"driver format expects 12 bins per year, got {nbins}")

    land = lsm > 0.5
    codes, soil_summary = soil_codes(config, land)

    # After the remap the driver's bins ARE the model's months, so a
    # precipitation total is the total for exactly the days it is spread over.
    bin_days = months

    # ExoPlaSim's longitudes run 0..360; LPJ-GUESS expects -180..180.
    #
    # Rounded to COORD_DECIMALS because LPJ-GUESS's soil map lookup keys on
    # std::pair<double,double> and compares it exactly. pedology/build_soil.py
    # writes its coordinates to the same precision, so the two round-trip to
    # identical doubles and the lookup needs no search radius. That agreement is
    # a contract between the two scripts, not a coincidence.
    lon_signed = np.round(np.where(lon > 180.0, lon - 360.0, lon), COORD_DECIMALS)
    lat = np.round(lat, COORD_DECIMALS)

    provenance = f"{config.get('source_build')}|{climatology.parent.name}".encode()
    provenance = provenance[:PROVENANCE_BYTES].ljust(PROVENANCE_BYTES, b"\0")

    output = args.output or (GENERATED / "vesper_driver.bin")
    output.parent.mkdir(parents=True, exist_ok=True)

    # Regolith depth, keyed on the same rounded coordinates the soil map uses.
    # Absent, every cell gets LPJ-GUESS's full 1.5 m profile, which is what the
    # unpatched model assumes anyway.
    default_depth_m = 1.5
    depth_by_coord: dict[tuple[float, float], float] = {}
    bedrock_by_coord: dict[tuple[float, float], float] = {}
    # Without a soil map, sub-bedrock layers keep the fresh-rock minimum. It is a
    # declared physical value, not a numerical guard; see pedogenesis.yaml.
    default_bedrock_fraction = 0.05
    soil_map = args.soil_map
    if soil_map is None:
        candidate = builds.soilmap()
        soil_map = candidate if candidate.is_file() else None
    if soil_map is not None:
        lines = soil_map.read_text().splitlines()
        header = lines[0].split()
        for needed in ("depth", "bedrockfrac"):
            if needed not in header:
                raise SystemExit(
                    f"{soil_map} has no {needed} column; rebuild it with "
                    f"build_soil.py")
        depth_column = header.index("depth")
        bedrock_column = header.index("bedrockfrac")
        for line in lines[1:]:
            parts = line.split()
            key = (round(float(parts[0]), COORD_DECIMALS),
                   round(float(parts[1]), COORD_DECIMALS))
            depth_by_coord[key] = float(parts[depth_column])
            bedrock_by_coord[key] = float(parts[bedrock_column])
        print(f"regolith depth from {soil_map.name}: {len(depth_by_coord)} cells")
    else:
        print("no soil map found; every cell gets the full 1.5 m profile")

    rows = np.argwhere(land)
    missing_depth = 0
    with output.open("wb") as handle:
        handle.write(MAGIC)
        handle.write(struct.pack("<iiii", len(rows), nbins, year_length, nyears))
        handle.write(struct.pack("<dd", co2_ppm, args.ndep))
        handle.write(provenance)
        for j, i in rows:
            handle.write(struct.pack("<dd", float(lon_signed[i]), float(lat[j])))
            handle.write(struct.pack("<ii", int(codes[j, i]), 0))
            key = (float(lon_signed[i]), float(lat[j]))
            depth = depth_by_coord.get(key)
            if depth is None:
                depth = default_depth_m
                missing_depth += 1
            handle.write(struct.pack("<d", depth))
            handle.write(struct.pack(
                "<d", bedrock_by_coord.get(key, default_bedrock_fraction)))
            # Flattened [year][bin], matching what vesperinput indexes.
            handle.write(tas[:, :, j, i].astype("<f8").tobytes())
            handle.write((pr[:, :, j, i] * bin_days[None, :]).astype("<f8").tobytes())
            handle.write(rss[:, :, j, i].astype("<f8").tobytes())
            handle.write(tsrange[:, :, j, i].astype("<f8").tobytes())

    weights = area_weights(lat, len(lon))
    lw = weights[land]

    def month_mean(field: np.ndarray) -> np.ndarray:
        """Mean over years and months, weighting each month by its own length."""
        per_year = np.tensordot(months, field, axes=(0, 1)) / months.sum()
        return per_year.mean(axis=0)

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
        "bins_per_year": nbins,
        "bin_days": bin_days.tolist(),
        # Both partitions, per climatology, so a later reader can see what was
        # remapped from what rather than having to re-derive it. The producer's
        # spans come from the climatology's own time axis; the driver's bins are
        # the model's months. They are NOT the same partition, and assuming they
        # were gave the last month 18 days of a rate measured over 16.09.
        "producer_interval_days": spans_y,
        "interval_remap": (
            "conservative overlap remap from the producer's time bins onto the "
            "model's months, biosphere/scripts/build_lpj_driver.py:"
            "interval_to_month. Total conserved exactly; the seasonal "
            "distribution is what moves."),
        "land_cells": int(len(rows)),
        "regolith_depth_source": (project_relative(soil_map)
                                  if soil_map else "none, full profile assumed"),
        "cells_without_depth": missing_depth,
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
            "to climate.dtr. biosphere/notes/ecological-forcing-field-contract.md"),
        "land_definition": "lsm from the climatology, itself built from surface_class",
        # Month-weighted, because the months are NOT equal: a plain mean over the
        # twelve would over-weight the eleven short ones.
        "land_mean_temperature_c": float(
            np.average(month_mean(tas)[land], weights=lw)),
        "land_mean_precip_mm_per_earth_year": float(
            np.average(month_mean(pr)[land], weights=lw)
            * orbit.EARTH_CALENDAR_YEAR_DAYS),
        "land_mean_net_sw_w_m2": float(
            np.average(month_mean(rss)[land], weights=lw)),
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

    print(f"land cells     {len(rows)}")
    print(f"year length    {year_length} days, months "
          f"{bin_days.astype(int).tolist()}")
    print("producer bins  " + ", ".join(f"{s:.2f}" for s in spans_y[0])
          + " days, remapped conservatively onto the months above")
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
