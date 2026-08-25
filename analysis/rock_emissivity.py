#!/usr/bin/env python3
"""Broadband thermal emissivity per Orogen rock class, from measured spectra.

    python analysis/rock_emissivity.py

Worldbuilding. Vesper is an invented planet; what follows is the derivation of
one surface property of the simulated land, and every number is a property of a
laboratory spectrum or of the model, not of anywhere real.

This is the derivation behind `surface.lithology_emissivity` in
`config/planet.yaml`, which `exoplasim/scripts/build_surface_emissivity.py`
turns into an ExoPlaSim boundary condition at surface code 177. It stands to
the longwave what `analysis/rock_albedo.py` stands to the shortwave, and it is
built out of the same library.

## What the model does with the number

`radmod.f90:lwr` carries a per-cell surface emissivity `zeps`, and until this
field existed it filled it from two scalars, ELWLAND and ELWSEA. The surface
net longwave is `-eps * (sigma*Ts^4 - LWdown)`, exactly linear in `eps`, so the
emissivity CONTRAST between two land surfaces is worth that contrast times the
same bracket the model's own surface longwave loss sits in. Nothing else in the
surface energy balance depends on it, which is what makes the size of the
effect computable in advance rather than only after a run.

## Kirchhoff, and the geometry that makes it usable

For an opaque sample in local thermodynamic equilibrium the directional
emissivity is one minus the directional-hemispherical reflectance at the same
wavelength. That identity is why only the hemispherical measurements in the
library are read here: a BIDIRECTIONAL reflectance samples one outgoing
direction and misses the rest of the hemisphere, so `1 - R` from it is an upper
bound on emissivity and not an estimate of it. Samples whose header does not
say hemispherical are dropped and counted.

## The two preparations, and why the answer is a bracket

`analysis/rock_albedo.py` documents the shortwave version of this trap: a
powder is not the slab it came from. In the thermal infrared the same
distinction runs the other way and is larger. The reststrahlen bands that pull
silicate emissivity down are surface-scattering features of a coherent
interface; break the material into particles comparable with the wavelength and
they weaken, so the same rock reads far closer to a blackbody as a powder than
as a slab.

A grid cell of this world's land is neither. It is some mixture of exposed
outcrop and the regolith weathered off it, and this project has no way to
measure that mixture. So the answer per class is a BRACKET, from the solid
samples to the particulate ones, and the value written into the config is the
midpoint. `--preparation solid` and `--preparation particulate` build the ends,
which is the arm to run if the mixture ever needs to be argued about.

## Where the spectra stop, and what that is worth

The library runs to about 14 or 15 micrometres for the hemispherical samples,
and a surface at this planet's land temperatures radiates roughly a third of
its energy beyond that. The in-band Planck-weighted emissivity is therefore
extended to the whole thermal spectrum, which assumes the out-of-band
emissivity equals the in-band one. The report carries the in-band Planck
fraction so the size of the assumption is visible, and the alternative bound --
a blackbody beyond the band edge -- is reported beside it, because it is the
one that would COMPRESS the between-class contrast and so the one that could
overturn the decision to carry a field at all.

## Soils are soils

`playa_clastic` is playa mud and alluvial fan fill, and `analysis/playa_albedo.py`
already argues that the right population for it is the aridisol soils of the
library rather than any rock. The same selection is used here, for the same
reason and from the same file list, so the two properties of that class describe
one material.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ECOSTRESS = ROOT / "references" / "ecospeclib-all"
REPORT = ROOT / "analysis" / "rock_emissivity.json"

H, C, KB = 6.626e-34, 2.998e8, 1.381e-23

# The Planck weight's temperature. A weighting temperature, not a result:
# emissivity is weighted by where the surface actually radiates. Bracketed
# below rather than asserted, and the bracket is reported.
WEIGHT_T_K = 290.0
WEIGHT_T_BRACKET_K = (250.0, 320.0)

# The narrowest band a sample must cover to be used. Below 8 um the Planck
# weight at these temperatures is small and the sample sets diverge; above
# 14 um only the bidirectional TIR series reaches, and Kirchhoff does not hold
# on it. Fixed before any spectrum was read.
BAND_UM = (8.0, 14.0)

# Preparations, by the library's own particle-size field.
SOLID = ("solid",)
PARTICULATE = ("coarse", "fine", "medium", "unsorted", "none")

# Orogen rock class -> the library populations that are that material.
# Each entry is (type, class, subclass) with None meaning "any". `water` is not
# a land surface and is absent; the ocean keeps ELWSEA.
SELECTORS: dict[str, list[tuple[str, str, str | None]]] = {
    "morb":                [("rock", "igneous", "mafic")],
    "oib":                 [("rock", "igneous", "mafic")],
    "flood_basalt":        [("rock", "igneous", "mafic")],
    "arc_basalt":          [("rock", "igneous", "mafic")],
    "arc_andesite":        [("rock", "igneous", "intermediate")],
    # Bimodal by definition: the mafic and felsic ends together, unweighted,
    # because the export carries no proportion for them.
    "rift_bimodal":        [("rock", "igneous", "mafic"),
                            ("rock", "igneous", "felsic")],
    "granite":             [("rock", "igneous", "felsic")],
    "granodiorite":        [("rock", "igneous", "intermediate")],
    "gneiss":              [("rock", "metamorphic", "gneis")],
    "schist":              [("rock", "metamorphic", "schist"),
                            ("rock", "metamorphic", "phyllite")],
    "quartzite":           [("rock", "metamorphic", "quartzite")],
    # Block-in-matrix: oceanic blocks in a sheared argillaceous matrix, so the
    # mafic and the pelitic-metamorphic populations together.
    "melange":             [("rock", "igneous", "mafic"),
                            ("rock", "metamorphic", "schist"),
                            ("rock", "metamorphic", "serpentinite")],
    "shelf_clastic":       [("rock", "sedimentary", "sandstone"),
                            ("rock", "sedimentary", "shale")],
    "carbonate":           [("rock", "sedimentary", "limestone"),
                            ("rock", "sedimentary", "dolomite")],
    "foreland_clastic":    [("rock", "sedimentary", "sandstone"),
                            ("rock", "sedimentary", "shale"),
                            ("rock", "sedimentary", "conglomerate")],
    "continental_clastic": [("rock", "sedimentary", "sandstone"),
                            ("rock", "sedimentary", "siltstone"),
                            ("rock", "sedimentary", "conglomerate")],
    "pelagic":             [("rock", "sedimentary", "shale"),
                            ("rock", "sedimentary", "limestone")],
    # Salt crust: the evaporite minerals themselves, not a rock population,
    # and selected BY MINERAL. `mineral.sulfate` and `mineral.halide` as whole
    # subclasses carry barite, celestite, jarosite, fluorite and cryolite,
    # which are not what precipitates out of a drying basin.
    "evaporite":           [("mineral", "sulfate", None),
                            ("mineral", "halide", None),
                            ("mineral", "chloride", None),
                            ("mineral", "carbonate", None)],
    # A soil, and the SAME aridisol suborders analysis/playa_albedo.py uses.
    "playa_clastic":       [("soil", "aridisol", None)],
}

# The suborders playa_albedo.py names. Kept here as one list so the two
# properties of playa_clastic are derived from one population.
PLAYA_SUBORDERS = ("salorthid", "gypsiorthid", "calciorthid", "camborthid",
                   "haplargid")

# THE EVAPORITE CLASS IS SELECTED BY MINERAL, AND HALITE IS HELD OUT.
#
# What precipitates out of a drying continental basin, in the order a brine
# gives it up: carbonate, then sulfate, then chloride. These are the minerals
# whose spectra are the class; the subclass as a whole is not, because it also
# carries barite, celestite, jarosite, fluorite, cryolite and the copper
# chlorides, none of which is an evaporite in this sense.
EVAPORITE_MINERALS = ("gypsum", "anhydrite", "mirabilite", "thenardite",
                      "glauberite", "epsomite", "hexahydrite", "polyhalite",
                      "kieserite", "bloedite", "aphthitalite", "calcite",
                      "aragonite", "dolomite", "magnesite")

# Halite is the endmember this library cannot settle. NaCl has no absorption
# band anywhere in the thermal window, so an optically thick pure halite powder
# is a poor emitter and the measurement says so -- and a natural salt crust is
# neither optically pure nor optically thick at these wavelengths, because it
# is a few millimetres of brine-cemented, dust-bearing crystal over damp playa
# mud whose radiation the crust transmits. The library carries no crust sample,
# so the pure-mineral value is the wrong PREPARATION in exactly the sense
# analysis/rock_albedo.py documents for powders against slabs, and further from
# the surface than a powder is. It is reported as a declared endmember and is
# not averaged into the class.
EVAPORITE_HELD_OUT = ("halite",)


def planck(microns: np.ndarray, t_k: float) -> np.ndarray:
    lam = microns * 1e-6
    return (2 * H * C ** 2 / lam ** 5) / (np.exp(H * C / (lam * KB * t_k)) - 1)


def planck_fraction(lo_um: float, hi_um: float, t_k: float) -> float:
    """Share of a blackbody's emission that falls inside a wavelength band."""
    full = np.geomspace(0.2, 500.0, 40000)
    band = np.geomspace(lo_um, hi_um, 4000)
    return float(np.trapezoid(planck(band, t_k), band)
                 / np.trapezoid(planck(full, t_k), full))


def header(path: str) -> dict[str, str]:
    out = {}
    with open(path, errors="replace") as fh:
        for line in fh:
            if ":" not in line:
                continue
            key, _, val = line.partition(":")
            key = key.strip().lower()
            if key in ("measurement", "name", "particle size", "wavelength range"):
                out[key] = val.strip()
            if len(out) == 4:
                break
    return out


def spectrum(path: str) -> tuple[np.ndarray, np.ndarray]:
    wl, refl = [], []
    with open(path, errors="replace") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) != 2:
                continue
            try:
                a, b = float(parts[0]), float(parts[1])
            except ValueError:
                continue
            wl.append(a)
            refl.append(b)
    w = np.asarray(wl)
    r = np.asarray(refl) / 100.0
    order = np.argsort(w)
    return w[order], r[order]


def emissivity(path: str, t_k: float) -> float | None:
    """Planck-weighted emissivity over BAND_UM, or None if the sample misses it."""
    wl, refl = spectrum(path)
    if wl.size < 50 or wl.min() > BAND_UM[0] or wl.max() < BAND_UM[1]:
        return None
    keep = (wl >= BAND_UM[0]) & (wl <= BAND_UM[1])
    wl, refl = wl[keep], refl[keep]
    if wl.size < 20:
        return None
    eps = 1.0 - np.clip(refl, 0.0, 1.0)
    w = planck(wl, t_k)
    return float(np.trapezoid(eps * w, wl) / np.trapezoid(w, wl))


def catalogue() -> list[dict]:
    """Every usable hemispherical spectrum, with the fields selection needs."""
    rows = []
    dropped_geometry = 0
    for path in sorted(glob.glob(os.path.join(ECOSTRESS, "*.spectrum.txt"))):
        parts = os.path.basename(path).split(".")
        if len(parts) < 5 or parts[0] not in ("rock", "soil", "mineral"):
            continue
        head = header(path)
        meas = head.get("measurement", "").lower()
        # Kirchhoff needs the WHOLE hemisphere. A bidirectional sample is an
        # upper bound on emissivity, not an estimate of it.
        if "hemispherical" not in meas or "bidirectional" in meas:
            dropped_geometry += 1
            continue
        rows.append({"path": path, "type": parts[0], "class": parts[1],
                     "subclass": parts[2], "size": parts[3],
                     "name": head.get("name", ""),
                     "stem": os.path.basename(path)})
    return rows, dropped_geometry


def matches(row: dict, sel: tuple[str, str, str | None]) -> bool:
    typ, cls, sub = sel
    if row["type"] != typ or row["class"] != cls:
        return False
    if sub is not None and row["subclass"] != sub:
        return False
    if typ == "soil" and cls == "aridisol":
        return any(s in row["stem"] for s in PLAYA_SUBORDERS)
    if typ == "mineral":
        name = row["name"].lower()
        return any(m in name for m in EVAPORITE_MINERALS)
    return True


def held_out(row: dict) -> bool:
    return (row["type"] == "mineral"
            and any(m in row["name"].lower() for m in EVAPORITE_HELD_OUT))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--weight-temperature", type=float, default=WEIGHT_T_K)
    ap.add_argument("--json", type=Path, default=REPORT)
    args = ap.parse_args()

    if not ECOSTRESS.is_dir():
        raise SystemExit(
            f"{ECOSTRESS} is absent. It is bulk reference data, excluded from "
            "git; references/INDEX.md records where the archive comes from.")

    rows, dropped = catalogue()
    t_k = args.weight_temperature

    cache: dict[str, float | None] = {}

    def eps_of(path: str) -> float | None:
        if path not in cache:
            cache[path] = emissivity(path, t_k)
        return cache[path]

    classes = {}
    endmembers = {}
    for code, sels in SELECTORS.items():
        got = {"solid": [], "particulate": []}
        for row in rows:
            if not any(matches(row, s) for s in sels) and not (
                    code == "evaporite" and held_out(row)):
                continue
            prep = ("solid" if row["size"] in SOLID
                    else "particulate" if row["size"] in PARTICULATE else None)
            if prep is None:
                continue
            value = eps_of(row["path"])
            if value is None:
                continue
            if code == "evaporite" and held_out(row):
                endmembers.setdefault(row["name"].split()[0].lower(),
                                      []).append(value)
                continue
            got[prep].append(value)
        entry = {}
        for prep in ("solid", "particulate"):
            arr = np.asarray(got[prep])
            entry[prep] = ({"n": 0} if arr.size == 0 else
                           {"n": int(arr.size), "mean": round(float(arr.mean()), 4),
                            "sd": round(float(arr.std(ddof=1)) if arr.size > 1 else 0.0, 4),
                            "min": round(float(arr.min()), 4),
                            "max": round(float(arr.max()), 4)})
        ends = [entry[p]["mean"] for p in ("solid", "particulate") if entry[p]["n"]]
        if not ends:
            raise SystemExit(f"no usable spectra for rock class {code}")
        entry["bracket"] = [round(min(ends), 4), round(max(ends), 4)]
        # A class the library measures in ONE preparation has no preparation
        # bracket, and saying so is not the same as saying its uncertainty is
        # zero. The soils and the salt crust are the two.
        entry["bracket_spans_preparations"] = all(
            entry[p]["n"] > 0 for p in ("solid", "particulate"))
        entry["emissivity"] = round(float(np.mean(ends)), 4)
        classes[code] = entry

    values = np.asarray([c["emissivity"] for c in classes.values()])
    in_band = planck_fraction(BAND_UM[0], BAND_UM[1], t_k)

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "quantity": "broadband thermal emissivity of the land surface, per "
                    "Orogen rock class, one minus directional-hemispherical "
                    "reflectance weighted by a Planck function",
        "library": "ECOSTRESS spectral library v1.0 (Meerdink et al. 2019), "
                   "which absorbs the ASTER library v2.0 (Baldridge et al. "
                   "2009); references/INDEX.md carries both",
        "band_um": list(BAND_UM),
        "weight_temperature_k": t_k,
        "weight_temperature_bracket_k": list(WEIGHT_T_BRACKET_K),
        "planck_fraction_in_band": round(in_band, 4),
        "out_of_band_assumption":
            "the in-band value is carried to the whole thermal spectrum. The "
            "opposing bound is a blackbody outside the band, which pulls every "
            "class toward 1 by (1 - in-band fraction) and is reported as "
            "blackbody_beyond_band because it is the bound that COMPRESSES the "
            "between-class contrast",
        "samples_dropped_non_hemispherical": dropped,
        "held_out_endmembers": {
            k: {"n": len(v), "mean": round(float(np.mean(v)), 4),
                "min": round(float(np.min(v)), 4),
                "max": round(float(np.max(v)), 4)}
            for k, v in sorted(endmembers.items())},
        "held_out_reason":
            "an optically thick pure halite powder has no absorption band in "
            "the thermal window and is measured as a poor emitter. A salt "
            "crust on a playa is neither optically pure nor optically thick "
            "there, and the library carries no crust sample, so the pure "
            "mineral is the wrong preparation rather than a low answer",
        "classes": classes,
        "spread": {
            "min": round(float(values.min()), 4),
            "max": round(float(values.max()), 4),
            "range": round(float(values.max() - values.min()), 4),
        },
        "blackbody_beyond_band": {
            "range": round(float((values.max() - values.min()) * in_band), 4),
        },
    }
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    width = max(len(c) for c in classes)
    for code, entry in sorted(classes.items(), key=lambda kv: kv[1]["emissivity"]):
        flag = " " if entry["bracket_spans_preparations"] else "*"
        print(f"{code:{width}s}  {entry['emissivity']:.4f} {flag}  bracket "
              f"{entry['bracket'][0]:.4f}-{entry['bracket'][1]:.4f}   "
              f"n solid {entry['solid']['n']:3d}  particulate "
              f"{entry['particulate']['n']:3d}")
    print("* one preparation only, so the bracket is not a preparation bracket")
    for name, e in report["held_out_endmembers"].items():
        print(f"held out: {name} {e['mean']:.4f} over {e['n']} samples")
    print(f"\nrange across classes {report['spread']['range']:.4f}; "
          f"Planck fraction inside {BAND_UM[0]}-{BAND_UM[1]} um at {t_k} K is "
          f"{in_band:.3f}")
    print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
