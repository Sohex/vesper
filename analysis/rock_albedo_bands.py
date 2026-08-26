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

It is also a prediction that can fail. If the ratio really is a property of the
mineral, then within a class the spread BETWEEN preparations must be smaller
than the spread between samples, because the samples differ in mineralogy and
the preparations do not. `preparation_spread_below_sample_spread` in the output
is that comparison, class by class, and it was fixed before the first run.

## What the bracket is made of, and what is left in it

The spread between preparations is narrow and the spread between SAMPLES is two
to four times wider, so the residual uncertainty in the two-band substrate is a
mineralogy question rather than a laboratory one. The selector is what answers a
mineralogy question, and until the sample-level rules below existed the selector
was a whole ECOSTRESS taxon. Reading `lithology.js`'s definition of each class
against the library's own `Name` and `Description` per spectrum -- which no
version of this script had consulted -- separates two things that were not
separate before, and the numbers say which one mattered.

**What the taxon was carrying that is not a rock class at all.** `rock.
sedimentary.shale` files seven whole-rock chips of a phosphate ORE, plus one
shale whose description is oolitic collophane, and every class naming shale
inherits all eight. Removing them is every narrowing worth a decimal place in
this table -- the rest, granodiorite's 0.033 the largest of them, is noise beside
it: `melange` 0.848 to 0.493 in bracket width, `pelagic` 1.136 to 0.811,
`foreland_clastic` 1.056 to 0.804, `shelf_clastic` 1.069 to 0.886,
`continental_clastic` 1.056 to 0.807. It MOVES the central ratio too, by 0.15 in
`shelf_clastic` and 0.11 in `continental_clastic`, so it is a correction and not
only a tightening.

**What a better-matched proxy bought, which was nothing.** `lithology.js`'s
`basin_clastic` rule is explicit that an intracratonic basin fills with "quieter,
more mature clastics" where a foreland basin takes an orogen's debris, and
maturity is mineralogy: it is how much feldspar and lithic debris survives beside
the quartz. Taking the library's own immature end -- greywacke, arkose, arkosic
sandstone -- out of `continental_clastic` therefore follows from the class
definition. It moves the central ratio by 0.001 and widens the bracket from 0.807
to 0.928. The same holds elsewhere: the per-class rules leave `shelf_clastic` and
`foreland_clastic` untouched and widen `pelagic` and `melange` slightly, and
taking the alkaline syenites out of `arc_andesite` widens it from 0.313 to 0.486.
Those rules stay, because a class's proxy is chosen by what the class IS and a
rule dropped for widening a bracket would be a tuning with an extra step. What
they establish is a NEGATIVE result: compositional maturity is not what sets the
near-infrared slope, and choosing better-matched samples along that axis does not
narrow the substrate.

**What does set it, and why no rule here can use it.** The library's descriptions
record iron directly -- ferruginous staining, limonite specks, hematite bands --
and splitting the sandstone taxon on that record separates 1.444 for the
iron-stained samples from 1.253 for the rest, with the unstained subset's bracket
half the width of the whole taxon's. That is the axis, and it is reported as
`residual_after_selection` rather than used, because `lithology.js` does not say
which side of it any clastic class sits on. A red bed and a grey quartz arenite
are both intracratonic clastics, the class definition chooses between them
nowhere, and five unstained samples is not a population. **So the substrate's
remaining bracket is not waiting on more spectra or a finer taxonomy: it is
waiting on a decision about this world's clastic iron that the rock table has not
made.**

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
               "measurements is not a population, and the bracket says so. "
               "The preparation test does not apply here and is expected to "
               "miss: the two entries are two liquids from two libraries, not "
               "one material prepared two ways.",
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


# Sample-level selection, BENEATH the taxon.
#
# The taxon above is the selector for a whole ECOSTRESS class, and that is as
# fine as a taxon can cut: `rock.sedimentary.sandstone` runs from a quartz
# arenite to an arkose, and `rock.igneous.felsic` files a hydrothermally
# altered tuff beside a granite. The library carries a `Name` and a
# `Description` per spectrum -- the mineralogy it was logged with, in the
# library's own words -- and neither was consulted until now.
#
# THE RULE IS THE ONE THE TAXON IS CHOSEN BY, applied one level down: a sample
# stands for a class if its MINERALOGY is the class's, because that is what
# sets where the absorption features sit and it is the only thing being
# transferred. Grain size, induration and weathering state are still not
# grounds -- they are what the preparation axis measures and what the bracket
# already carries -- so a rhyolite is kept as a proxy for granite and a gabbro
# for basalt. What comes out is a sample whose MINERALS are not the class's.
#
# Every rule below is stated from `lithology.js`'s definition of the class and
# from the library's own description of the sample, and all of them were fixed
# before any ratio was recomputed. A rule that named a sample because dropping
# it tightened the bracket would be a tuning with an extra step.
#
# Matched case-insensitively as a regular expression against the header field
# named in the second element.

# Filed under a taxon whose mineralogy is not the taxon's, so wrong for every
# class that names that taxon.
NOT_THE_TAXON = (
    (r"phosphorite", "Name",
     "phosphate ore, filed under rock.sedimentary.shale: whole-rock chips of "
     "the Phosphoria Formation's Meade Peak Member whose spectrum is "
     "hydroxylapatite, buddingtonite and muscovite. No class here is an ore"),
    (r"collophane|phosphatic material", "Description",
     "oolitic collophane, the same apatite mineralogy under the name Shale "
     "(Phosphatic)"),
    (r"\baltered volcanic tuff", "Name",
     "an alteration assemblage rather than a rock class: what the spectrum is "
     "OF is the alteration. No class here is an alteration halo"),
    (r"obsidian", "Name",
     "volcanic glass, which has no mineral assemblage to transfer. The rule "
     "these proxies are chosen by has nothing to say about an amorphous solid"),
)

# Wrong for one class, by that class's own definition.
NOT_THE_CLASS: dict[str, tuple] = {
    # The basalt assemblage is plagioclase with pyroxene and olivine.
    # Anorthosite is a plagioclase cumulate, better than nine parts feldspar,
    # so it is a different rock rather than a coarse one.
    "morb": ((r"anorthosite", "Name", "plagioclase cumulate, not basalt"),),
    "oib": ((r"anorthosite", "Name", "plagioclase cumulate, not basalt"),),
    "flood_basalt": ((r"anorthosite", "Name",
                      "plagioclase cumulate, not basalt"),),
    "arc_basalt": ((r"anorthosite", "Name", "plagioclase cumulate, not basalt"),),
    "rift_bimodal": ((r"anorthosite", "Name",
                      "plagioclase cumulate, not basalt"),),
    "melange": ((r"anorthosite", "Name",
                 "plagioclase cumulate: the mafic blocks in a melange are "
                 "basalt and gabbro, not an anorthosite body"),),
    # 'Continental-arc andesite / dacite' and 'Arc-root granodiorite' are
    # calc-alkaline. A nepheline syenite is silica-undersaturated and carries a
    # feldspathoid, which is an assemblage neither class has.
    "arc_andesite": ((r"syenite", "Name",
                      "alkaline, feldspathoid-bearing where the class is "
                      "calc-alkaline andesite to dacite"),),
    "granodiorite": ((r"nepheline syenite", "Name",
                      "silica-undersaturated and feldspathoid-bearing, which "
                      "an arc-root granodiorite is not"),),
    # 'Carbonate platform'. Verde Antique and the serpentine marbles are
    # serpentinite veined with calcite, and the mineral doing the absorbing is
    # serpentine.
    "carbonate": ((r"serpentine marble|verde antique", "Name",
                   "serpentinite, not a carbonate assemblage"),),
    # lithology.js's basin_clastic rule is what separates these two: "Foreland
    # basins sit against an orogen and fill with its debris; intracratonic
    # basins fill with quieter, more mature clastics." Maturity IS mineralogy
    # -- it is how much feldspar and lithic debris survives beside the quartz
    # -- so the immature end of the sandstone taxon belongs to foreland_clastic
    # and not to continental_clastic. foreland_clastic keeps the whole range,
    # because molasse and flysch span it.
    "continental_clastic": (
        (r"greywacke|arkos", "Name",
         "the library's own names for the lithic and feldspathic end of the "
         "sandstone taxon, which is the debris a foreland basin takes and the "
         "opposite of the mature clastic an intracratonic basin fills with"),),
    # 'Pelagic ooze / abyssal clay' has no sand fraction by definition, and the
    # library names the samples that do.
    "pelagic": ((r"arenaceous", "Name",
                 "a sand fraction, which an abyssal clay or ooze does not have"),),
}


def sample_verdict(code: str, head: dict) -> str | None:
    """Why this sample is not a proxy for this class, or None to keep it."""
    for pattern, field, reason in NOT_THE_TAXON + NOT_THE_CLASS.get(code, ()):
        if re.search(pattern, head.get(field, ""), re.IGNORECASE):
            return reason
    return None


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


def collect(code: str, spec: dict, star_wl, star_flux):
    """(taxon rows, selected rows, what the selection dropped) for one class.

    A row is (preparation, b1, b2, broadband, band1 flux share). The taxon rows
    are every qualifying spectrum the class's ECOSTRESS selectors name; the
    selected rows are those whose mineralogy is the class's, by
    `sample_verdict`. Both are returned so that what the narrowing bought is a
    measurement rather than an assertion.
    """
    taxon, kept, dropped = [], [], {}

    def consider(head, row):
        taxon.append(row)
        why = sample_verdict(code, head)
        if why is None:
            kept.append(row)
        else:
            dropped[head.get("Sample No.", "?")] = {
                "name": head.get("Name", "?"), "why": why}

    for pattern in spec.get("eco", ()):
        for path in sorted(glob.glob(str(ECOSTRESS / (pattern + "spectrum.txt")))):
            head, arr = read_ecostress(Path(path))
            if arr is None:
                continue
            got = band_integrals(arr[:, 0], arr[:, 1] / 100.0, star_wl, star_flux)
            if got:
                consider(head, (preparation_of(Path(path).name), *got))
    for pattern, mineral in spec.get("eco_named", ()):
        for path in sorted(glob.glob(str(ECOSTRESS / (pattern + "spectrum.txt")))):
            head, arr = read_ecostress(Path(path))
            if arr is None or not head.get("Name", "").lower().startswith(mineral.lower()):
                continue
            got = band_integrals(arr[:, 0], arr[:, 1] / 100.0, star_wl, star_flux)
            if got:
                consider(head, (preparation_of(Path(path).name), *got))
    for name in spec.get("poseidon", ()):
        path = POSEIDON / name
        if not path.is_file():
            continue
        table = np.loadtxt(path)
        order = np.argsort(table[:, 0])
        got = band_integrals(table[order, 0], table[order, 1], star_wl, star_flux)
        if got:
            consider({"Name": name, "Sample No.": name}, ("poseidon", *got))
    return taxon, kept, dropped


def summarise(rows) -> dict:
    """The central ratio, the bracket and the two spreads, over one set of rows.

    Called twice per class -- once on the taxon and once on what survives the
    mineralogy selection -- so the two are computed by the same code and the
    difference between them is the selection and nothing else.
    """
    rho = np.array([b2 / b1 for _p, b1, b2, _br, _g in rows if b1 > 0])
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
    prep_spread = (max(prep_median.values()) - min(prep_median.values())
                   if len(prep_median) > 1 else None)
    sample_spread = float(np.percentile(rho, 90) - np.percentile(rho, 10))
    return {
        "spectra": len(rows),
        "band2_over_band1": round(central, 4),
        "band2_over_band1_bracket": [round(lo, 4), round(hi, 4)],
        "by_preparation": {k: round(v, 4) for k, v in prep_median.items()},
        "preparation_spread": (None if prep_spread is None
                               else round(prep_spread, 4)),
        "sample_spread_p10_p90": round(sample_spread, 4),
        # The method's own prediction: preparation must matter less than
        # mineralogy, because the ratio is claimed to be a property of the
        # mineral. Null where the class has only one preparation.
        "preparation_spread_below_sample_spread": (
            None if prep_spread is None else bool(prep_spread < sample_spread)),
        "measured_band1_flux_share": round(
            float(np.median([g for *_x, g in rows])), 4),
        # A class standing on fewer than five spectra has a bracket that is
        # a spread between samples rather than an estimate of one, and is
        # marked rather than quoted as though it were the latter.
        "thin": len(rows) < 5,
    }


# The axis the selection above cannot use, measured so that the reason it
# cannot is a number rather than an opinion. Iron content sets where the
# near-infrared slope goes, the library records it per sample, and no class
# definition in `lithology.js` says which side of it a clastic class sits on.
# This is a DIAGNOSTIC and never a selector: wiring it in would be choosing a
# mineralogy for the world from the shape of a bracket.
IRON_IN_DESCRIPTION = r"ferruginous|limonite|hematite|red |reddish|purple|brownstone"


def residual_after_selection(star_wl, star_flux) -> dict:
    """What splitting the sandstone taxon on the library's own record of iron
    is worth, against what one proxy per taxon is worth."""
    rows = []
    for path in sorted(glob.glob(str(ECOSTRESS / "rock.sedimentary.sandstone.*"
                                     "spectrum.txt"))):
        head, arr = read_ecostress(Path(path))
        if arr is None:
            continue
        got = band_integrals(arr[:, 0], arr[:, 1] / 100.0, star_wl, star_flux)
        if got:
            stained = bool(re.search(
                IRON_IN_DESCRIPTION,
                head.get("Name", "") + " " + head.get("Description", ""),
                re.IGNORECASE))
            rows.append((stained, (preparation_of(Path(path).name), *got)))

    def side(sel):
        chosen = [r for stained, r in rows if sel(stained)]
        if not chosen:
            return None
        got = summarise(chosen)
        lo, hi = got["band2_over_band1_bracket"]
        return {"spectra": got["spectra"],
                "band2_over_band1": got["band2_over_band1"],
                "bracket_width": round(hi - lo, 4)}

    return {
        "axis": "iron, as the ECOSTRESS description records it: ferruginous "
                "staining, limonite, hematite, and the colour names that "
                "follow from them",
        "taxon": "rock.sedimentary.sandstone",
        "why_not_a_selector": (
            "vendor/orogen/js/lithology.js does not say which side of this "
            "axis any clastic class sits on. A red bed and a grey quartz "
            "arenite are both intracratonic clastics, and choosing between "
            "them from the width of the bracket they produce would be a "
            "tuned value. This is what the remaining bracket is waiting on."),
        "all": side(lambda _s: True),
        "iron_stained": side(lambda st: st),
        "not_iron_stained": side(lambda st: not st),
    }


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
        taxon_rows, rows, dropped = collect(code, spec, star_wl, star_flux)
        if not taxon_rows:
            raise SystemExit(
                f"no qualifying spectrum for {code!r}. The proxy selectors name "
                "material that this library does not carry over "
                f"{WORKING[0]}-{WORKING[1]} um; fix the selector rather than "
                "dropping the class, which would leave it spectrally flat.")
        if not rows:
            raise SystemExit(
                f"the mineralogy selection left {code!r} with no sample at all. "
                "A rule in NOT_THE_TAXON or NOT_THE_CLASS is excluding the "
                "class's own rock; fix the rule rather than falling back to the "
                "taxon, which would hide it.")
        selected = summarise(rows)
        taxon = summarise(taxon_rows)
        blo, bhi = selected["band2_over_band1_bracket"]
        tlo, thi = taxon["band2_over_band1_bracket"]
        classes[code] = {
            "why": spec["why"],
            **selected,
            # What the sample-level selection did, so the narrowing is a
            # measurement and not a claim. `taxon` is the whole ECOSTRESS
            # taxon, which is what this table used before the library's own
            # sample descriptions were read.
            "excluded_samples": dropped,
            "taxon": taxon,
            "bracket_width": round(bhi - blo, 4),
            "taxon_bracket_width": round(thi - tlo, 4),
            "bracket_narrowed_by": round((thi - tlo) - (bhi - blo), 4),
            "central_moved_by": round(
                selected["band2_over_band1"] - taxon["band2_over_band1"], 4),
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
        "method_check": (
            "preparation_spread_below_sample_spread, per class: the ratio is "
            "claimed to be a property of the mineral, so within a class the "
            "spread between sample preparations must be smaller than the "
            "spread between samples. Fixed before the first run."),
        "recombination_tolerance": RECOMBINATION_TOLERANCE,
        "residual_flag_threshold_albedo": RESIDUAL_FLAG,
        "sample_selection": (
            "Each class's proxy is an ECOSTRESS taxon narrowed to the samples "
            "whose MINERALOGY is the class's, read from each spectrum's own "
            "Name and Description against the class definition in "
            "lithology.js. `taxon` per class is what the whole taxon gives, "
            "so what the narrowing bought is a measurement. Grain size, "
            "induration and weathering are not grounds for exclusion: they "
            "are what the preparation axis measures."),
        "residual_after_selection": residual_after_selection(star_wl, star_flux),
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
    tested = {k: r for k, r in report_classes.items()
              if r["preparation_spread_below_sample_spread"] is not None}
    missed = [k for k, r in tested.items()
              if not r["preparation_spread_below_sample_spread"]]
    print(f"preparation matters less than mineralogy in "
          f"{len(tested) - len(missed)} of {len(tested)} classes carrying more "
          "than one preparation"
          + ("" if not missed else
             ". The exceptions are " + ", ".join(missed) + ", where the "
             "cancellation this method rests on is weakest, so read their "
             "brackets and not their central values"))
    narrowed = {k: r["bracket_narrowed_by"] for k, r in report_classes.items()
                if r["bracket_narrowed_by"] > 0}
    widened = {k: -r["bracket_narrowed_by"] for k, r in report_classes.items()
               if r["bracket_narrowed_by"] < 0}
    print(f"the mineralogy selection dropped "
          f"{sum(len(r['excluded_samples']) for r in report_classes.values())} "
          f"class-sample pairs; it narrowed {len(narrowed)} brackets and "
          f"widened {len(widened)}")
    if narrowed:
        best = sorted(narrowed, key=narrowed.get, reverse=True)[:4]
        print("  narrowed most: "
              + ", ".join(f"{k} by {narrowed[k]:.3f}" for k in best))
    if widened:
        worst = sorted(widened, key=widened.get, reverse=True)[:4]
        print("  widened: " + ", ".join(f"{k} by {widened[k]:.3f}" for k in worst)
              + " -- kept, because a proxy is chosen by what the class is")
    res = report["residual_after_selection"]
    print(f"  what is left is iron: {res['taxon']} splits "
          f"{res['iron_stained']['band2_over_band1']:.3f} stained against "
          f"{res['not_iron_stained']['band2_over_band1']:.3f} unstained, and "
          "no class definition says which side a clastic class is on")
    print(f"wrote {args.output.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
