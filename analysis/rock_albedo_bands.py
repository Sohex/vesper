#!/usr/bin/env python3
"""Two-band shape factors for every Orogen rock class, under this world's star.

    python analysis/rock_albedo_bands.py

ExoPlaSim's shortwave scheme has two bands meeting at 0.75 um, and at
`NSIMPLEALBEDO = 0` the radiation reads the background albedo pair in codes 175
and 176 rather than the broadband 174. `vendor/orogen/js/lithology.js` carries
one broadband number per rock class and nothing carries a pair, so
`exoplasim/scripts/build_surface_albedo.py` had no split to write for the
substrate and wrote the same field to all three codes. That asserts of every
lithology on the planet the thing a measured reflectance spectrum says is false:
equal reflectance either side of 0.75 um. This script is the missing half.

## What is derived here, and what is NOT

Derived: the SHAPE, one number per class -- `band2_over_band1`, the ratio of the
class's stellar-flux-weighted reflectance above 0.75 um to its reflectance
below. Everything else in the output follows from that ratio and from the class's
existing broadband albedo.

Not derived: the LEVEL. The broadband albedo per class stays exactly what
`lithology.js` and `config/planet.yaml`'s `lithology_albedo_overrides` say it is.
Rebuilding those levels from spectra is a different question with a different
blocker -- laboratory preparation, argued at length in
`notes/audits/orogen-lithology.md` -- and it is not touched here.

That division is what makes this derivable at all, and it is the same argument
`analysis/playa_albedo.py` and `analysis/vegetation_albedo.py` both use. A
laboratory spectrum of a sieved powder is far brighter than the outcrop it came
from: Paragas et al. (2025) measure powder over slab at 2.19 to 5.11 on the same
rock. That offset is in the numerator and the denominator of a band ratio alike,
so it cancels from the ratio while dominating the level. What survives the
cancellation is where a mineral's absorption features sit, which is a property of
the mineral rather than of the sample jar.

The cancellation is not perfect and this script measures how imperfect rather
than asserting it: ECOSTRESS carries the same rock types at `solid`, `coarse`,
`medium` and `fine` preparation, so the ratio's sensitivity to preparation is
read off the library instead of being argued about. That sensitivity is one of
the two things the reported bracket is made of.

## Anchoring, and the residual `bio-18` asks for

The model recombines the pair with its own band weights: `radmod.f90` forms
`zsolars(1)*dsalb(1) + zsolars(2)*dsalb(2)`, and `lib/stellar.band_fractions`
reproduces `solarini` exactly, minwavel cut and band-edge interval included. The
pair written here is therefore anchored on THOSE weights,

    a1 = broadband / (z1 + z2 * rho),    a2 = rho * a1

so that `z1*a1 + z2*a2` returns the broadband albedo as an identity rather than
as an approximation. `build_surface_albedo.py` checks it cell by cell on the
fields it is about to write.

The alternative anchoring -- normalising on the flux fraction over the range the
spectra were MEASURED over -- does not satisfy that identity, because the
measured range is not the star's whole shortwave. The difference is reported per
class as `range_mismatch_residual`: it is what the recombined broadband would
come out as, minus what it should be, had the pair been anchored the naive way.
It is signed and it is not small enough to ignore, which is the point of
reporting it. `analysis/vegetation_albedo.py` carries the same number for the
canopy pair.

## Data

`references/ecospeclib-all`, the ECOSTRESS spectral library v1.0 (Meerdink et
al. 2019), 12,292 directional-hemispherical reflectance spectra -- the geometry
a model albedo wants, and the same library the canopy and playa derivations use.
Water is the one class ECOSTRESS does not cover over the full shortwave, and it
comes from the POSEIDON surface-albedo database instead.

`references/INDEX.md` records the library as "the basis for rebuilding Orogen's
rock_albedo table and for the two-band split at 0.75 um". This is the second
half of that sentence.
"""

from __future__ import annotations

import argparse
import datetime
import glob
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "lib"))
from lib.paths import rel  # noqa: E402
from lib import stellar  # noqa: E402

ECOSTRESS = ROOT / "references" / "ecospeclib-all"
POSEIDON = ROOT / "references" / "poseidon_surface_albedo"
OUTPUT = ROOT / "analysis" / "rock_albedo_bands.json"

BAND_SPLIT_UM = stellar.BAND_SPLIT_UM

# The range essentially every ECOSTRESS laboratory spectrum covers, and the
# range every spectrum is integrated over so the ratios are comparable. A
# spectrum that does not span it is dropped rather than extrapolated: padding a
# tail changes the ratio, which is the quantity being measured.
WORKING = (0.40, 2.49)

# How far from the working range a spectrum's own limits may sit and still be
# integrated over it. 0.401 admits the library's 0.4 um start against floating
# point; 2.49 admits both the 2.5 um vswir files and the combined 0.4-15 um ones.
QUALIFY = (0.401, 2.49)

# The check that can fail, fixed before any number below was computed.
#
# The pair is anchored on the model's own band weights, so the recombination is
# an identity and its only error is arithmetic. The bound is set from the
# INSTRUMENT and not from the size of any effect: `.sra` is written at
# `%12.5f`, so a value round-trips to the model with a quantum of 5e-6, and a
# recombination of two such values cannot be held tighter than that. Twice the
# quantum leaves the check able to fail on anything that is a construction error
# rather than a formatting one.
RECOMBINATION_TOLERANCE = 1.0e-5

# The reported residual is a different quantity and gets a different bound, in
# albedo rather than as a fraction. A class whose naive-anchoring residual
# exceeds this is flagged, because at that point the choice of anchor is worth
# more than the last digit the rock table carries: `lithology.js` quotes albedo
# to two decimals, so half its last digit is 0.005.
RESIDUAL_FLAG = 0.005

LITHOLOGY_JS = ROOT / "vendor" / "orogen" / "js" / "lithology.js"
PLANET_CONFIG = ROOT / "config" / "planet.yaml"

# Orogen rock class -> the ECOSTRESS material that stands in for it.
#
# One rule throughout: the proxy is chosen for its MINERALOGY, because that is
# what sets where the absorption features sit and therefore the only thing being
# transferred. Grain size, induration and weathering state are chosen for by the
# preparation spread instead, and end up in the bracket.
#
# `land_fraction` is the share of land the class holds on `precarve-craton`,
# carried here only so the table's own attention is auditable: the four classes
# above 0.10 are playa_clastic, continental_clastic, schist and granite, and
# those are the four the proxies had better be right about. Read the current
# numbers from `world_state.json`, never from here.
PROXIES: dict[str, dict] = {
    "water": {
        "why": "open water, and what a persistent lake is painted with. The "
               "thinnest class in the table: ECOSTRESS's seawater spectra are "
               "thermal-infrared only, so the shortwave comes from POSEIDON's "
               "USGS ocean spectrum with ECOSTRESS tap water beside it as the "
               "only independent measurement of the same liquid. Two "
               "measurements is not a population, and the bracket says so.",
        "eco": ("water.tapwater.*",),
        "poseidon": ("GoodisGordon2025-GG25/USGS_ocean_seawater_GG25.txt",),
    },
    "morb": {"why": "basalt", "eco": ("rock.igneous.mafic.*", "rock.igneous.basic.*")},
    "oib": {"why": "basalt", "eco": ("rock.igneous.mafic.*", "rock.igneous.basic.*")},
    "flood_basalt": {"why": "basalt",
                     "eco": ("rock.igneous.mafic.*", "rock.igneous.basic.*")},
    "arc_basalt": {
        "why": "basalt to basaltic andesite, so mafic and intermediate together",
        "eco": ("rock.igneous.mafic.*", "rock.igneous.basic.*",
                "rock.igneous.intermediate.*")},
    "arc_andesite": {"why": "andesite to dacite",
                     "eco": ("rock.igneous.intermediate.*",)},
    "rift_bimodal": {
        "why": "bimodal by definition: basalt and rhyolite, no intermediate term",
        "eco": ("rock.igneous.mafic.*", "rock.igneous.basic.*",
                "rock.igneous.felsic.*")},
    "granite": {"why": "granite", "eco": ("rock.igneous.felsic.*",)},
    "granodiorite": {"why": "granodiorite sits between granite and diorite",
                     "eco": ("rock.igneous.felsic.*", "rock.igneous.intermediate.*")},
    "gneiss": {"why": "gneiss", "eco": ("rock.metamorphic.gneis.*",)},
    "schist": {"why": "schist and phyllite, which is what the class is named",
               "eco": ("rock.metamorphic.schist.*", "rock.metamorphic.phyllite.*")},
    "quartzite": {"why": "quartzite", "eco": ("rock.metamorphic.quartzite.*",)},
    "melange": {
        "why": "block-in-matrix: mafic and cherty blocks in a sheared "
               "argillaceous matrix, so shale and slate with a mafic term",
        "eco": ("rock.sedimentary.shale.*", "rock.metamorphic.slate.*",
                "rock.igneous.mafic.*")},
    "shelf_clastic": {"why": "sandstone and shale",
                      "eco": ("rock.sedimentary.sandstone.*",
                              "rock.sedimentary.shale.*")},
    "carbonate": {"why": "limestone, dolomite and marble",
                  "eco": ("rock.sedimentary.limestone.*",
                          "rock.sedimentary.dolomite.*",
                          "rock.metamorphic.marble.*")},
    "foreland_clastic": {"why": "molasse and flysch: sandstone, shale, conglomerate",
                         "eco": ("rock.sedimentary.sandstone.*",
                                 "rock.sedimentary.shale.*",
                                 "rock.sedimentary.conglomerate.*")},
    "continental_clastic": {"why": "intracratonic sandstone, siltstone and shale",
                            "eco": ("rock.sedimentary.sandstone.*",
                                    "rock.sedimentary.siltstone.*",
                                    "rock.sedimentary.shale.*")},
    "pelagic": {"why": "abyssal clay and ooze: argillaceous, with a carbonate term "
                       "for the ooze",
                "eco": ("rock.sedimentary.argillaceou.*",
                        "rock.sedimentary.shale.*",
                        "rock.sedimentary.limestone.*")},
    "evaporite": {
        "why": "salt crust. Halite is the mineral the class is named for and the "
               "one whose spectrum is featureless through the near infrared; the "
               "hydrous sulfate is carried beside it because a real crust is not "
               "pure halite and gypsum's structural-water bands run the other "
               "way. references/INDEX.md, the Slater 1987 row.",
        "eco_named": (("mineral.halide.*vswir*", "Halite"),
                      ("mineral.sulfate.*vswir*", "Gypsum")),
    },
    "playa_clastic": {
        "why": "playa mud and fan fill, which is a soil rather than a rock. The "
               "same aridisol suborders analysis/playa_albedo.py selects, so the "
               "level and the shape describe one material.",
        "eco": ("soil.aridisol.salorthid.*", "soil.aridisol.gypsiorthid.*",
                "soil.aridisol.calciorthid.*", "soil.aridisol.camborthid.*",
                "soil.aridisol.haplargid.*"),
    },
}


def table_albedos() -> dict[str, float]:
    """The broadband level per class, from the table that defines it.

    Read from `lithology.js` rather than from an export, because the shapes are
    a property of the lithology and not of a build, and read from the source
    rather than from any prose copy of it. `config/planet.yaml`'s
    `lithology_albedo_overrides` is applied on top, since that is the level the
    climate actually receives.
    """
    src = LITHOLOGY_JS.read_text(encoding="utf-8")
    rows = re.findall(
        r"\{\s*id:\s*(\d+),\s*code:\s*'([a-z_]+)'.*?albedo:\s*([0-9.]+)\s*\}", src)
    if not rows:
        raise SystemExit(f"no ROCK_CLASSES parsed from {LITHOLOGY_JS}; the table "
                         "has been reshaped and this reader must follow it")
    table = {code: float(value) for _id, code, value in rows}
    import yaml
    config = yaml.safe_load(PLANET_CONFIG.read_text(encoding="utf-8"))
    overrides = (config.get("model") or {}).get("lithology_albedo_overrides") or {}
    applied = {}
    for code, spec in overrides.items():
        value = float(spec["albedo"] if isinstance(spec, dict) else spec)
        applied[code] = value
        table[code] = value
    return table, applied


def read_ecostress(path: Path):
    """(header, ascending Nx2 of wavelength um and reflectance percent)."""
    head, rows = {}, []
    for line in path.read_text(errors="ignore").splitlines():
        if re.match(r"^\s*[\d.]+\s+[-\d.]+\s*$", line):
            a, b = line.split()[:2]
            rows.append((float(a), float(b)))
        elif ":" in line:
            key, _, value = line.partition(":")
            head.setdefault(key.strip(), value.strip())
    if not rows:
        return head, None
    arr = np.array(rows)
    return head, arr[arr[:, 0].argsort()]


def band_integrals(wl, refl, star_wl, star_flux):
    """(band1, band2, broadband, band1 flux share) over WORKING, or None.

    All four against the same stellar spectrum and the same wavelength range, so
    the ratio of the first two is a property of the material.
    """
    if wl.min() > QUALIFY[0] or wl.max() < QUALIFY[1]:
        return None
    grid = star_wl[(star_wl >= WORKING[0]) & (star_wl <= WORKING[1])]
    if grid.size < 32:
        return None
    r = np.interp(grid, wl, refl)
    f = np.interp(grid, star_wl, star_flux)
    lo = grid <= BAND_SPLIT_UM
    if lo.sum() < 8 or (~lo).sum() < 8:
        return None

    def mean(sel):
        return float(np.trapezoid(r[sel] * f[sel], grid[sel])
                     / np.trapezoid(f[sel], grid[sel]))

    w1 = float(np.trapezoid(f[lo], grid[lo]))
    w2 = float(np.trapezoid(f[~lo], grid[~lo]))
    return mean(lo), mean(~lo), mean(np.ones_like(grid, dtype=bool)), w1 / (w1 + w2)


def preparation_of(name: str) -> str:
    """ECOSTRESS puts the sample preparation in the fourth dotted field."""
    parts = name.split(".")
    return parts[3] if len(parts) > 4 else "unknown"


def collect(spec: dict, star_wl, star_flux):
    """Every qualifying spectrum for one class, as (preparation, b1, b2, broad, g1)."""
    out = []
    for pattern in spec.get("eco", ()):
        for path in sorted(glob.glob(str(ECOSTRESS / (pattern + "spectrum.txt")))):
            _head, arr = read_ecostress(Path(path))
            if arr is None:
                continue
            got = band_integrals(arr[:, 0], arr[:, 1] / 100.0, star_wl, star_flux)
            if got:
                out.append((preparation_of(Path(path).name), *got))
    for pattern, mineral in spec.get("eco_named", ()):
        for path in sorted(glob.glob(str(ECOSTRESS / (pattern + "spectrum.txt")))):
            head, arr = read_ecostress(Path(path))
            if arr is None or not head.get("Name", "").lower().startswith(mineral.lower()):
                continue
            got = band_integrals(arr[:, 0], arr[:, 1] / 100.0, star_wl, star_flux)
            if got:
                out.append((preparation_of(Path(path).name), *got))
    for name in spec.get("poseidon", ()):
        path = POSEIDON / name
        if not path.is_file():
            continue
        table = np.loadtxt(path)
        order = np.argsort(table[:, 0])
        got = band_integrals(table[order, 0], table[order, 1], star_wl, star_flux)
        if got:
            out.append(("poseidon", *got))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--library", type=Path, default=ECOSTRESS)
    ap.add_argument("--output", type=Path, default=OUTPUT)
    args = ap.parse_args()
    if not args.library.is_dir():
        raise SystemExit(f"{args.library} is absent. It is bulk reference data, "
                         "excluded from git; unzip references/ecospeclib_*.zip "
                         "beside this path.")

    table, overridden = table_albedos()
    missing = set(table) - set(PROXIES)
    if missing:
        raise SystemExit(
            f"rock classes with no spectral proxy: {sorted(missing)}. A class "
            "with no entry here would be written spectrally flat, which is the "
            "defect this table exists to remove; add a proxy rather than "
            "letting it fall through.")

    low_path, hires = stellar.spectrum_paths()
    data = np.loadtxt(hires, skiprows=1)
    star_wl, star_flux = data[:, 0], data[:, 1]
    z1, z2 = (float(v) for v in stellar.band_fractions(hires))

    # How much of the star this library never saw. The anchoring hides it from
    # the broadband answer but not from the ratio, so it is stated rather than
    # left for someone to find.
    total = float(np.trapezoid(star_flux, star_wl))
    inside = (star_wl >= WORKING[0]) & (star_wl <= WORKING[1])
    covered = float(np.trapezoid(star_flux[inside], star_wl[inside])) / total

    classes: dict[str, dict] = {}
    for code, spec in PROXIES.items():
        rows = collect(spec, star_wl, star_flux)
        if not rows:
            raise SystemExit(
                f"no qualifying spectrum for {code!r}. The proxy selectors name "
                "material that this library does not carry over "
                f"{WORKING[0]}-{WORKING[1]} um; fix the selector rather than "
                "dropping the class, which would leave it spectrally flat.")
        rho = np.array([b2 / b1 for _p, b1, b2, _br, _g in rows if b1 > 0])
        g1 = float(np.median([g for *_x, g in rows]))
        preps: dict[str, list[float]] = {}
        for prep, b1, b2, _br, _g in rows:
            if b1 > 0:
                preps.setdefault(prep, []).append(b2 / b1)
        prep_median = {k: float(np.median(v)) for k, v in sorted(preps.items())}
        central = float(np.median(rho))
        # The bracket is measured twice and the wider answer is kept: the spread
        # between preparations, which is the offset this method claims cancels,
        # and the 10-90 spread within the class, which is how well one proxy
        # stands for a lithology at all.
        lo = min([central, *prep_median.values(), float(np.percentile(rho, 10))])
        hi = max([central, *prep_median.values(), float(np.percentile(rho, 90))])
        classes[code] = {
            "why": spec["why"],
            "spectra": len(rows),
            "band2_over_band1": round(central, 4),
            "band2_over_band1_bracket": [round(lo, 4), round(hi, 4)],
            "by_preparation": {k: round(v, 4) for k, v in prep_median.items()},
            "measured_band1_flux_share": round(g1, 4),
            # A class standing on fewer than five spectra has a bracket that is
            # a spread between samples rather than an estimate of one, and is
            # marked rather than quoted as though it were the latter.
            "thin": len(rows) < 5,
        }

    # The pair, per class, for whatever broadband level the caller holds. Two
    # shape factors rather than two albedos, because the level is not this
    # script's to set: `lithology_albedo_overrides` moves playa_clastic and a
    # lake moves whatever is under it.
    report_classes = {}
    for code, row in classes.items():
        rho = row["band2_over_band1"]
        s1 = 1.0 / (z1 + z2 * rho)
        s2 = rho * s1
        g1 = row["measured_band1_flux_share"]
        n1 = 1.0 / (g1 + (1 - g1) * rho)
        residual = z1 * n1 + z2 * rho * n1 - 1.0
        blo, bhi = row["band2_over_band1_bracket"]
        level = table[code]
        report_classes[code] = {
            **row,
            "broadband": level,
            "broadband_source": ("config model.lithology_albedo_overrides"
                                 if code in overridden
                                 else "vendor/orogen/js/lithology.js ROCK_CLASSES"),
            "band_pair": [round(level * s1, 4), round(level * s2, 4)],
            "range_mismatch_residual_albedo": round(residual * level, 5),
            "flagged": bool(abs(residual * level) > RESIDUAL_FLAG),
            "band1_shape": round(s1, 4),
            "band2_shape": round(s2, 4),
            "band1_shape_bracket": [round(1.0 / (z1 + z2 * bhi), 4),
                                    round(1.0 / (z1 + z2 * blo), 4)],
            "band2_shape_bracket": [round(blo / (z1 + z2 * blo), 4),
                                    round(bhi / (z1 + z2 * bhi), 4)],
            "range_mismatch_residual_relative": round(residual, 5),
        }

    report = {
        "generated": datetime.date.today().isoformat(),
        "generator": "analysis/rock_albedo_bands.py",
        "spectrum": rel(hires),
        "spectrum_sha256": hashlib.sha256(hires.read_bytes()).hexdigest(),
        "library": rel(args.library),
        "band_split_um": BAND_SPLIT_UM,
        "model_band_flux_fractions": [round(z1, 6), round(z2, 6)],
        "model_band_flux_source": "lib/stellar.band_fractions, which reproduces "
                                  "radmod.f90 solarini including the minwavel cut "
                                  "and the band-edge interval",
        "working_range_um": list(WORKING),
        "stellar_flux_fraction_inside_working_range": round(covered, 4),
        "recombination_tolerance": RECOMBINATION_TOLERANCE,
        "residual_flag_threshold_albedo": RESIDUAL_FLAG,
        "lithology_table": rel(LITHOLOGY_JS),
        "overridden_levels": overridden,
        "classes": report_classes,
        "how_to_use": (
            "band_i = broadband * band_i_shape, with broadband whatever the "
            "caller's rock table says. z1*band1 + z2*band2 = broadband is then "
            "an identity in the model's own weights. The shapes are keyed by "
            "rock CODE, not by id, because ids move when a class is added."),
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"{len(report_classes)} classes against {hires.name}, "
          f"z1 = {z1:.5f}  z2 = {z2:.5f}")
    print(f"working range {WORKING[0]}-{WORKING[1]} um carries "
          f"{covered:.1%} of the star's shortwave")
    print(f"{'class':22}{'n':>4}{'b2/b1':>8}{'bracket':>16}"
          f"{'broad':>7}{'band1':>8}{'band2':>8}{'resid':>9}")
    for code, row in report_classes.items():
        blo, bhi = row["band2_over_band1_bracket"]
        b1, b2 = row["band_pair"]
        print(f"{code:22}{row['spectra']:4d}{row['band2_over_band1']:8.3f}"
              f"{blo:8.3f}-{bhi:<7.3f}{row['broadband']:7.2f}{b1:8.4f}{b2:8.4f}"
              f"{row['range_mismatch_residual_albedo']:+9.4f}"
              f"{' *' if row['flagged'] else ''}")
    print(f"wrote {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
