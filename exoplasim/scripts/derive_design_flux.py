#!/usr/bin/env python3
"""Derive the design flux from the comfort-band scoring `docs/src/pipeline/state.md` section 5b describes.

WORLDBUILDING CONTEXT: Vesper is a fictional planet and this script scores a
toy GCM's simulated seasons to pick a stellar flux for it. Nothing here refers
to the real world.

WHY IT EXISTS. The 2026-08-16 derivation that set `baseline_flux_earth` lived
only as prose: no script, no artifact, no named runs, no stated thresholds
(CLIM-24, inherited-earth-constants.md finding 5). docs/src/pipeline/sequencing.md loop C requires the flux
to be RE-DERIVED on every new terrain, so the criteria have to be re-runnable,
not reinvented. This codifies them. Every number below is fixed in this
docstring before the script is run against the data, per the project's
threshold rule.

THE METHOD, as 5b records it: seasonal temperature per latitude band, measured
on converged runs at two fluxes spanning the target, projected across
candidate means, scored on a warm-season comfort band. Concretely here:

- A "month" is one of the 12 output bins; per-cell warmest and coldest bin
  means of `tas` are the two seasons. Bin weights come from
  `lib/climatology.py` (they differ between I/O regimes).
- The projection is per 10-degree band: amplification = (band seasonal delta)
  / (global annual delta) between the two sources, applied to every land cell
  in the band, with the global mean moving along `lib/sensitivity.py`'s
  canonical slope. Bracket between two points, never extrapolate from one --
  the doctrine 5b itself was written under.
- The 0.91 source contributes `tas` only: its spin-up ran the cheap I/O
  regime, whose historical first-record defect touched wind and humidity, so
  humidity comes solely from the clean baseline climatology.

WHICH HALF A RADIATION CHANGE REACHES, because the shortwave cloud optics have
just moved the modelled climate by up to 20.6 K and the two halves of this
docstring answer to different things. The comfort band is a set of thresholds in
DEGREES of simulated surface air temperature: a statement about the kind of
world this project wants, not about the model that produces it. Changing the
radiation moves which stellar flux delivers a given temperature field, so it
moves the CANDIDATE RANGE and leaves the band where it is.
`exoplasim/notes/trace-gas-absorbers.md` section 4 records the same mechanism
from the other side: a radiation scheme short of greenhouse forcing reaches
these same thresholds at a higher stellar flux, and the search finds it. So the
range below is re-derived under the corrected optics and the band is not.

DECLARED THRESHOLDS, each with its status.

    comfort:  warmest-bin mean <= 33.0 C  AND  coldest-bin mean >= -25.0 C
    harsh:    warm side 33 to 38 C; cold side -25 down to -40 C
    extreme:  above 38 C, or below -40 C

    score = land-area fraction in the comfort band; the winner maximizes it.

`cold_floor_c`, `warm_extreme_c` and `cold_extreme_c` are DECLARED DESIGN
PREFERENCES and are KEPT. They say what kind of world this is meant to be, they
are fixed here ahead of the runs that will be judged against them, and nothing
in the model or in any run determines them. A design preference does not move
when the physics does, so the cloud optics do not reach them.

`warm_ceiling_c` is the one that was NOT a preference. Its recorded
justification was that the anchor flux "puts the tropics near +33 C in their
warmest month rather than +36", which is a number read off where the superseded
physics landed -- the same shape as the cold-extreme cap CLIM-30 exists to undo.
The corrected optics remove that justification outright: at the anchor flux the
modelled tropics are now several kelvin colder, so the sentence no longer picks
out 33. The value is RE-DECLARED here as a preference in its own right, at the
same number, and it now stands on what the other three stand on: the
warm-season monthly mean above which land on this world is not meant to read as
comfortable. A re-declaration at the old value is exactly where a tuning would
hide, so it is not left bare -- it is declared WITH A BRACKET THAT GETS SWEPT,
`THRESHOLD_BRACKET` below, and the report carries the design flux each bracket
point returns, so a reader can see what the answer owes to the preference. That
is the third of the four dispositions `docs/src/practice/conventions.md` allows
under "No tuned values".

`cold_extreme_cap` is declared by the caller and starts as None. The purged
artifact INFERRED it from the anchor's own row, which is a threshold fitted to
its own answer; this script now refuses to run until it is set. CLIM-30.

THE CANDIDATE RANGE, re-derived because the cloud optics moved under it. The
range's only job is to contain the winner, and a search that cannot reach its
own answer returns the edge of its range instead. Both bounds are built from
where the winner can be and then given margin, and the refusal below is what
makes the construction falsifiable rather than merely generous.

- The corrected optics COOL. world-f9ig put the tuned tswr1/tswr2/tswr3 fits
  aside for Stephens, Ackerman and Smith (1984)'s own tables, and
  `exoplasim/analysis/stephens_tables_vs_fits.json` brackets the flux ratio that
  offsets that at +0.026 to +0.102. That is a magnitude bound over a black
  surface with the whole of the model's cover on one layer, so the true offset
  is at most the top of it and may sit below the bottom. world-jgen's band-1
  correction pushes the other way and is currently unquantified: its artifact
  was computed under the deleted band-1 fit and is worthless, not stale
  (world-sii1).
- CEILING 1.100. The highest the winner can sit is the superseded-physics
  winner under the inferred cap, 0.945, plus the whole of the offset bound,
  which is 1.047. The margin above that is 0.053, or 10.7 K at
  `lib/sensitivity.py`'s slope and twice the width of the harsh-warm band. It
  covers a declared cap stricter than the inferred one, which moves the winner
  up, and the terrain moving between builds.
- FLOOR 0.790. The lowest the winner can sit is the superseded-physics
  UNCONSTRAINED winner, 0.8725, which the corrected optics can only push up.
  The margin below that is 0.0825, or 16.7 K, which clears the 15 K width of the
  harsh-cold band even at unit cold-season amplification, and 5b records the
  winter amplification as several times the summer's. It covers world-jgen's
  opposing warming and each carve iteration darkening and warming the world.
- STEP 0.0025, unchanged. Half a kelvin of global mean, finer than anything the
  projection resolves, so the winner is limited by the criteria and not by the
  grid.
- The range is deliberately wider than the bound requires, because the two
  costs are not symmetric. Too wide costs arithmetic on candidates that cannot
  win. Too narrow returns an edge, silently, and that is the failure this
  re-derivation exists to remove.

THE HUMIDITY-COUPLED VARIANT, because the audit showed the dry score stands in
for a quantity that couples temperature to humidity, and this world's aridity
is distributed by drainage rather than latitude. The variant scores the
isobaric equivalent temperature T_e = tas + (L/cp) q at the lowest level --
plain thermodynamics, no other content. Its ceiling is pinned by declaration:
chosen so the comfort fraction at the pin flux matches the dry score's, so the
two columns differ only in HOW the margin is distributed, and the measured
divergence is the winner shift and the per-band gap. The pin flux is a
convention and its only requirements are that it be fixed ahead of the run and
lie inside the candidate range; it is set to the anchor because that is a value
already written down, and its standing as a prior derivation plays no part.

THE ANCHOR IS A PRIOR, NOT A TARGET. `anchor_flux` is what
`config/planet.yaml` carries as PROVISIONAL. The report measures the distance to
it and nothing depends on that distance being small; under the corrected optics
agreement would be surprising, because the anchor was chosen under the tuned
cloud coefficients and the offset bound above says the answer has moved.

Checks that can fail, all of which REFUSE rather than report and carry on:
`cold_extreme_cap` must be declared; the anchor and the pin flux must be
interior to the candidate range; the two sources must span more than 5 K of
global mean; the amplification must be positive in every band for the warm
season; the declared cap must admit at least one candidate; and no winner --
dry or equivalent, capped or not -- may sit on either end of the candidate
range. The last of those was a reported flag in the artifact and is now a
refusal, because the range is derived and an edge winner means the derivation
was wrong.

Writes `exoplasim/analysis/design_flux.json`. Registered as step
`design_flux`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset
from numpy.polynomial.legendre import leggauss

from _paths import ANALYSIS, CONFIG, RUNS  # noqa: F401  (puts lib/ on sys.path)

import climatology  # noqa: E402
import sensitivity  # noqa: E402
from paths import climatology_path, rel  # noqa: E402

LATENT_OVER_CP = 2.501e6 / 1004.9   # K per unit specific humidity, plasim's constants

DECLARED = {
    "warm_ceiling_c": 33.0,      # design preference, re-declared, swept below
    "cold_floor_c": -25.0,       # design preference
    "warm_extreme_c": 38.0,      # design preference
    "cold_extreme_c": -40.0,     # design preference
    "cold_extreme_cap": None,    # declared by the caller; CLIM-30. None refuses
    "band_degrees": 10.0,
    "candidates": [round(0.790 + 0.0025 * i, 4) for i in range(125)],
    "anchor_flux": 0.945,
    "te_pin_flux": 0.945,
}

# The bracket the two comfort thresholds are declared with, swept so the report
# says what the design flux owes to the preference rather than asserting it does
# not. The warm ceiling's bracket is +-3 K, inside its own 5 K harsh band; the
# cold floor's is +-5 K, a third of its 15 K one.
THRESHOLD_BRACKET = {
    "warm_ceiling_c": [30.0, 31.5, 33.0, 34.5, 36.0],
    "cold_floor_c": [-30.0, -27.5, -25.0, -22.5, -20.0],
}


def seasonal_fields(path: Path) -> dict:
    """Per-cell warmest/coldest bin means of tas, plus masks and weights."""
    with Dataset(path) as ds:
        tas = np.asarray(ds["tas"][:], dtype=float)
        w = climatology.bin_weights(np.asarray(ds["time"][:], dtype=float))
        lat = np.asarray(ds["lat"][:], dtype=float)
        nlat, nlon = tas.shape[1:]
        land = np.asarray(ds["lsm"][:], dtype=float).mean(axis=0) > 0.5
        q = (np.asarray(ds["hus"][:], dtype=float)[:, -1]
             if "hus" in ds.variables and ds["hus"].ndim == 4 else None)
    gw = leggauss(nlat)[1][::-1]
    area = np.broadcast_to(gw[:, None], (nlat, nlon)) / (2.0 * nlon)
    annual = np.tensordot(w, tas, axes=(0, 0))
    return {
        "warm": tas.max(axis=0), "cold": tas.min(axis=0),
        "annual_global": float((annual * area).sum()),
        "q_warm": (None if q is None else
                   q[np.argmax(tas, axis=0),
                     np.arange(nlat)[:, None], np.arange(nlon)[None, :]]),
        "lat": lat, "land": land, "area": area,
    }


def tail_mean_fields(run_dir: Path, orbits: list[int]) -> Path:
    """Average the tail orbits' binned files into one 12-bin scratch file.

    Plain mean across orbits per bin, which is what build_climatology does;
    written to the analysis directory's scratch space so seasonal_fields can
    read one path either way.
    """
    import tempfile
    acc, count = {}, 0
    for orbit in orbits:
        p = run_dir / f"MOST.{orbit:05d}.nc"
        with Dataset(p) as ds:
            for name in ("tas", "lsm", "time"):
                a = np.asarray(ds[name][:], dtype=float)
                acc[name] = acc.get(name, 0.0) + a
            lat = np.asarray(ds["lat"][:], dtype=float)
            lon = np.asarray(ds["lon"][:], dtype=float)
        count += 1
    out = Path(tempfile.gettempdir()) / f"design_flux_tail_{run_dir.name}.nc"
    with Dataset(out, "w") as ds:
        ds.createDimension("time", acc["tas"].shape[0])
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))
        for name, dims in (("time", ("time",)), ("lat", ("lat",)), ("lon", ("lon",))):
            v = ds.createVariable(name, "f8", dims)
            v[...] = {"time": acc["time"] / count, "lat": lat, "lon": lon}[name]
        for name in ("tas", "lsm"):
            v = ds.createVariable(name, "f8", ("time", "lat", "lon"))
            v[...] = acc[name] / count
    return out


def band_index(lat: np.ndarray, width: float) -> np.ndarray:
    return np.floor((lat + 90.0) / width).astype(int)


def score(warm_c: np.ndarray, cold_c: np.ndarray, land: np.ndarray,
          area: np.ndarray, warm_ceiling: float,
          cold_floor: float | None = None) -> dict:
    d = DECLARED
    cold_floor = d["cold_floor_c"] if cold_floor is None else cold_floor
    la = float((area * land).sum())
    def frac(mask):
        return round(float((area * (land & mask)).sum() / la), 5)
    return {
        "comfort": frac((warm_c <= warm_ceiling) & (cold_c >= cold_floor)),
        "harsh_warm": frac((warm_c > warm_ceiling) & (warm_c <= d["warm_extreme_c"])),
        "extreme_warm": frac(warm_c > d["warm_extreme_c"]),
        "harsh_cold": frac((cold_c < cold_floor) & (cold_c >= d["cold_extreme_c"])),
        "extreme_cold": frac(cold_c < d["cold_extreme_c"]),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--bracket-run", default="run_bfa3f5269660",
                    help="run supplying the second flux point's tas tail")
    ap.add_argument("--tail-orbits", type=int, default=10)
    ap.add_argument("--output", type=Path, default=ANALYSIS / "design_flux.json")
    args = ap.parse_args()

    # Refuse before a single field is read, so nothing about these can be
    # contaminated by what the data turns out to say.
    cap = DECLARED["cold_extreme_cap"]
    if cap is None:
        raise SystemExit(
            "cold_extreme_cap is not declared. The comfort-maximizing rule alone "
            "does not settle this flux: the trade 5b records is tropics against "
            "poles, and 'polar margins severe but small' is a constraint on the "
            "cold-extreme land fraction that the prose never quantified. Set "
            "DECLARED['cold_extreme_cap'] to the largest fraction of land this "
            "world accepts below the cold-extreme threshold, in this file, before "
            "the runs it will be applied to exist. CLIM-30; it must not be "
            "inferred from the answer, which is what the purged artifact did.")
    edges = (DECLARED["candidates"][0], DECLARED["candidates"][-1])
    for name in ("anchor_flux", "te_pin_flux"):
        f = DECLARED[name]
        if f in edges or f not in DECLARED["candidates"]:
            raise SystemExit(f"{name} = {f} is not interior to the candidate "
                             f"range {edges[0]} to {edges[1]}; the pin and the "
                             "prior must both be scored on a row the search can "
                             "reach from either side")

    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    base_path = climatology_path()
    base = seasonal_fields(base_path)
    f0 = float(cfg["orbit"]["baseline_flux_earth"])

    run_dir = RUNS / args.bracket_run
    manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    f1 = float(manifest.get("stellar_flux_ratio",
               manifest.get("derived_parameters", {}).get("stellar_flux_ratio_earth")))
    last = max(int(p.stem.split(".")[1]) for p in run_dir.glob("MOST.[0-9]*.nc"))
    orbits = list(range(last - args.tail_orbits + 1, last + 1))
    other = seasonal_fields(tail_mean_fields(run_dir, orbits))

    d_global = base["annual_global"] - other["annual_global"]
    if abs(d_global) < 5.0:
        raise SystemExit(f"the two sources span only {d_global:.2f} K of global "
                         "mean; a bracket this weak cannot support a projection")

    bands = band_index(base["lat"], DECLARED["band_degrees"])
    nb = bands.max() + 1
    amp = {"warm": np.ones(nb), "cold": np.ones(nb)}
    band_rows = []
    for b in range(nb):
        rows = bands == b
        sel = base["land"][rows, :]
        if not sel.any():
            continue
        a_rows = base["area"][rows, :]
        wsum = float((a_rows * sel).sum())
        entry = {"band": f"{b*10-90:+d} to {b*10-80:+d}"}
        for season in ("warm", "cold"):
            t_base = float((base[season][rows, :] * a_rows * sel).sum() / wsum)
            t_other = float((other[season][rows, :] * a_rows * sel).sum() / wsum)
            amp[season][b] = (t_base - t_other) / d_global
            entry[f"{season}_c_at_{f0}"] = round(t_base - 273.15, 2)
            entry[f"amp_{season}"] = round(amp[season][b], 3)
        band_rows.append(entry)
    if any(r.get("amp_warm", 1) <= 0 for r in band_rows):
        raise SystemExit("a band's warm-season amplification is not positive; "
                         "the projection cannot be trusted")

    slope = sensitivity.SLOPE_K_PER_FLUX_RATIO
    amp_warm_cells = amp["warm"][bands][:, None]
    amp_cold_cells = amp["cold"][bands][:, None]

    te_warm0 = base["warm"] + LATENT_OVER_CP * base["q_warm"]

    def projected(f):
        """Warm, cold and equivalent-warm land fields at a candidate flux, in C."""
        dT = slope * (f - f0)
        return (base["warm"] + amp_warm_cells * dT - 273.15,
                base["cold"] + amp_cold_cells * dT - 273.15,
                te_warm0 + amp_warm_cells * dT - 273.15)

    # the humidity-coupled ceiling, pinned so the pin-flux comfort fractions match
    pin_warm, pin_cold, pin_te = projected(DECLARED["te_pin_flux"])
    dry_pin = score(pin_warm, pin_cold, base["land"], base["area"],
                    DECLARED["warm_ceiling_c"])
    lo, hi = 20.0, 90.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        s = score(pin_te, pin_cold, base["land"], base["area"], mid)
        if s["comfort"] < dry_pin["comfort"]:
            lo = mid
        else:
            hi = mid
    te_ceiling = round(0.5 * (lo + hi), 3)

    table = []
    for f in DECLARED["candidates"]:
        warm, cold, te_warm = projected(f)
        table.append({
            "flux": f,
            "dry": score(warm, cold, base["land"], base["area"],
                         DECLARED["warm_ceiling_c"]),
            "equivalent": score(te_warm, cold, base["land"], base["area"],
                                te_ceiling),
        })

    def winner(kind, cap=None):
        rows = [r for r in table if cap is None or r[kind]["extreme_cold"] <= cap]
        return max(rows, key=lambda r: r[kind]["comfort"])["flux"] if rows else None

    win_dry, win_te = winner("dry"), winner("equivalent")
    win_dry_capped, win_te_capped = winner("dry", cap=cap), winner("equivalent", cap=cap)
    if win_dry_capped is None or win_te_capped is None:
        raise SystemExit(
            f"the declared cold-extreme cap {cap} admits no candidate in "
            f"{edges[0]} to {edges[1]}: this world cannot be placed at any flux "
            "the search carries without exceeding the cap. That is a result about "
            "the cap and the terrain together, and it is not resolved by moving "
            "either one after the fact")

    # An edge winner is a refusal, not a flag. The candidate range is derived in
    # this file's docstring from where the corrected cloud optics can put the
    # balance; a winner sitting on an end means that derivation was wrong, and
    # the number the search would return is the end of the range rather than an
    # optimum. Widen the range and say why, in the docstring, before re-running.
    at_edge = sorted({n for n, f in (("dry", win_dry), ("equivalent", win_te),
                                     ("dry under the cap", win_dry_capped),
                                     ("equivalent under the cap", win_te_capped))
                      if f in edges})
    if at_edge:
        raise SystemExit(
            f"the winner is on the end of the candidate range for: "
            f"{', '.join(at_edge)}. The range {edges[0]} to {edges[1]} does not "
            "contain its own answer, so what the search would return is an edge "
            "and not an optimum. Re-derive the range in the docstring and re-run")

    # What the answer owes to the two comfort preferences, measured rather than
    # asserted: the design flux each bracket point of each threshold returns,
    # the other threshold and the cap held at their declared values.
    def winner_at(warm_ceiling, cold_floor):
        best_flux, best_comfort = None, -1.0
        for f in DECLARED["candidates"]:
            warm, cold, _ = projected(f)
            s = score(warm, cold, base["land"], base["area"], warm_ceiling, cold_floor)
            if s["extreme_cold"] > cap:
                continue
            if s["comfort"] > best_comfort:
                best_flux, best_comfort = f, s["comfort"]
        return best_flux

    sensitivity_rows = {
        "warm_ceiling_c": {str(v): winner_at(v, DECLARED["cold_floor_c"])
                           for v in THRESHOLD_BRACKET["warm_ceiling_c"]},
        "cold_floor_c": {str(v): winner_at(DECLARED["warm_ceiling_c"], v)
                         for v in THRESHOLD_BRACKET["cold_floor_c"]},
    }

    # Reported, never used as the cap: the smallest cold-extreme tolerance under
    # which the recorded prior would have been the constrained optimum. It is
    # what the purged artifact solved backwards, kept visible so the defect is
    # legible rather than repeated.
    anchor_row = next(r for r in table if r["flux"] == DECLARED["anchor_flux"])
    prior_implied_cap = anchor_row["dry"]["extreme_cold"]

    report = {
        "generated": datetime.now(timezone.utc).isoformat(),
        "generator": "exoplasim/scripts/derive_design_flux.py",
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "declared": DECLARED,
        "threshold_bracket": THRESHOLD_BRACKET,
        "sources": {
            "baseline": {"path": str(base_path), "flux": f0},
            "bracket": {"run": args.bracket_run, "flux": f1,
                        "tail_orbits": orbits,
                        "note": "tas only; its spin-up I/O regime historically "
                                "corrupted wind and humidity, so humidity comes "
                                "from the clean climatology alone"},
            "global_mean_span_k": round(d_global, 3),
            "slope_k_per_unit_flux": slope,
        },
        "bands": band_rows,
        "equivalent_ceiling_c": round(te_ceiling, 2),
        "equivalent_ceiling_note": "pinned so the pin-flux comfort fraction "
                                   "matches the dry score's; the divergence "
                                   "between columns is then distribution, not level",
        "candidates": table,
        "design_flux": win_dry_capped,
        "design_flux_equivalent": win_te_capped,
        "winner_dry_uncapped": win_dry,
        "winner_equivalent_uncapped": win_te,
        "projection_extrapolation_k": {
            "sources_span": [min(f0, f1), max(f0, f1)],
            "beyond_sources_at_design_flux": round(
                slope * max(0.0, min(f0, f1) - win_dry_capped,
                            win_dry_capped - max(f0, f1)), 2),
            "note": "how far past the two measured points the band "
                    "amplifications are carried to reach the design flux, in "
                    "global-mean kelvin; zero means the answer is interpolated",
        },
        "threshold_sensitivity": {
            "note": "the design flux each bracket point of a comfort threshold "
                    "returns, the other threshold and the cap held at their "
                    "declared values. This is the sweep the warm ceiling is "
                    "declared with, and it measures what the answer owes to a "
                    "preference instead of asserting it owes nothing",
            **sensitivity_rows,
        },
        "against_the_recorded_prior": {
            "prior": DECLARED["anchor_flux"],
            "distance": round(win_dry_capped - DECLARED["anchor_flux"], 4),
            "cap_the_prior_implies": prior_implied_cap,
            "note": "the prior is config/planet.yaml's PROVISIONAL "
                    "baseline_flux_earth, chosen under the tuned tswr1/tswr2/"
                    "tswr3 cloud coefficients world-f9ig deleted. Agreement with "
                    "it is not a success and disagreement is not a failure; the "
                    "distance is reported because the prior fixes a semi-major "
                    "axis that is compiled into the biosphere. cap_the_prior_"
                    "implies is the tolerance that would have made the prior "
                    "optimal, reported so the backwards-solved threshold stays "
                    "legible; it is not the declared cap and must never be used "
                    "as one.",
        },
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"sources: {f0} (baseline) and {f1} ({args.bracket_run}), "
          f"span {d_global:.2f} K")
    print(f"candidates {edges[0]} to {edges[1]}, cold-extreme cap {cap}")
    print(f"design flux: dry {win_dry_capped}  equivalent {win_te_capped}  "
          f"(uncapped: dry {win_dry}, equivalent {win_te})")
    print(f"recorded prior {DECLARED['anchor_flux']}, distance "
          f"{win_dry_capped - DECLARED['anchor_flux']:+.4f}")
    print(f"warm-ceiling bracket: {sensitivity_rows['warm_ceiling_c']}")
    print(f"cold-floor bracket:   {sensitivity_rows['cold_floor_c']}")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
