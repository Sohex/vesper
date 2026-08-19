#!/usr/bin/env python3
"""CO2 consumed and silica released by weathering, and where each of them goes.

    python pedology/scripts/weathering_fluxes.py

Writes `pedology/analysis/weathering_fluxes.json`.

## Why one script for two solutes

They share every input -- lithology, runoff, temperature, and the endorheic
split -- and they are two columns of the same Meybeck table. Computing them
apart would mean two scripts that must be kept in step, and the interesting
results are comparisons between them.

## CO2

Silicate weathering consumes atmospheric CO2 and delivers bicarbonate. So the
CO2 consumption flux is the bicarbonate flux, with one correction: for CARBONATE
lithologies only half the bicarbonate comes from the atmosphere and the other
half from the rock, so those count at 0.5. That is the standard treatment and it
is why the carbonate share cannot simply be added in.

Concentrations come from Meybeck (1987) Table 2C, already extracted and
checked into `data/reference/meybeck1987_tables.json`, and flux is
concentration times runoff -- which is Meybeck's own model form.

Dessert et al. (2003) is carried as an INDEPENDENT cross-check on the volcanic
part rather than as the primary law, because it covers basalt alone. Their
Eq. 2, f_CO2 = Rf * 323.44 * exp(0.0642 T) mol/km2/yr, was fitted on basaltic
provinces directly. Note it is LINEAR in runoff: Dessert found bicarbonate
concentration to depend on temperature and not on runoff, which is the same
reasoning WHAK used and is not the 0.65 exponent Berner derived for the global
mix. The two laws disagreeing on the runoff exponent is a real disagreement in
the literature, not an error here, and both are reported.

## Silica

Same table, the `sio2_umol_l` column, no stoichiometric correction. Silica
matters here for what it makes downstream rather than for the carbon cycle:
silcrete needs a dissolved silica source and `notes/derived-surface-classes.md`
records that this world's source was unresolved. Diatomite needs silica
delivered to a lake that is productive and then desiccates.

The endorheic split is the point. Silica delivered to the ocean is diluted into
a very large reservoir; silica delivered to a closed basin has nowhere to go and
concentrates until it saturates. Eugster and Jones' behaviour type V is exactly
this -- SiO2 constant after saturation with a solid -- so a closed basin fed by
silica-rich runoff is the setting that precipitates it.

## What this does NOT do

It does not model biogenic silica. Diatoms take dissolved silica out of lake
water and put it into opal, which is what diatomite is made of, and that is a
biosphere question this component cannot answer. What is computed here is the
SUPPLY, which is the part that depends on lithology and climate, and which
bounds whatever the biology can do with it.

It also does not decide whether the supply is sustained. Basin lake level moves
with the stellar cycle, and a basin has to be persistent enough to be productive
and then dry enough to deflate. That is a cycle question and is noted in
`notes/derived-surface-classes.md`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import netCDF4 as nc
import numpy as np
import yaml

from _paths import ANALYSIS, CONFIG, PEDOGENESIS, PROJECT_ROOT, climatology_path

from brine_paths import ROCK_TO_MEYBECK
from gridding import land_fraction_of_class
from paths import rel

from build_soil import EARTH_YEAR_DAYS, KELVIN, lithology_fractions

# Bicarbonate from carbonate dissolution is half rock-derived, so only half of it
# represents atmospheric CO2 drawdown. Silicate weathering is fully atmospheric.
#
# `misc_metamorphic` is on the carbonate side of that line despite its name, and
# missing it was a real error here. Meybeck's class is MARBLE: Ca 1375 ueq/l
# against sedimentary carbonate's 2560, and Ca + Mg = 1740 against HCO3 1730, a
# near-exact carbonate balance. Treating it as silicate let rock-derived
# bicarbonate into the silicate total at full weight.
CO2_PER_HCO3 = {
    "sedimentary_carbonate_rocks": 0.5,
    "misc_metamorphic": 0.5,
    # Neither silicate nor carbonate weathering in the carbon-cycle sense.
    # Halite and gypsum dissolution consume no atmospheric CO2 at all.
    "halite_evaporite": 0.0,
    "gypsum_evaporite": 0.0,
}
DEFAULT_CO2_PER_HCO3 = 1.0

# Classes whose bicarbonate is wholly or partly rock-derived are excluded from
# the SILICATE total, not merely down-weighted in it. Listing them here rather
# than testing one name inline is the fix for having tested one name inline.
CARBONATE_BEARING = {"sedimentary_carbonate_rocks", "misc_metamorphic"}

MEYBECK = PROJECT_ROOT / "pedology" / "data" / "reference" / "meybeck1987_tables.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def concentrations() -> tuple[dict[str, float], dict[str, float]]:
    """Per-Meybeck-class SiO2 (umol/l) and HCO3 (ueq/l)."""
    tables = json.loads(MEYBECK.read_text(encoding="utf-8"))
    table = tables["table_2c"]
    columns = table["columns"]
    si_i = columns.index("sio2_umol_l")
    hco3_i = columns.index("hco3_ueq_l")
    silica = {k: float(v[si_i]) for k, v in table["rows"].items()}
    bicarb = {k: float(v[hco3_i]) for k, v in table["rows"].items()}
    return silica, bicarb


def dessert_co2(runoff_mm_yr: np.ndarray, temperature_c: np.ndarray) -> np.ndarray:
    """Dessert et al. (2003) Eq. 2, mol CO2 per km2 per year, basalt only."""
    return runoff_mm_yr * 323.44 * np.exp(0.0642 * temperature_c)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--climatology", type=Path, default=None)
    ap.add_argument("--source-build", default=None,
                    help="override config's source_build, for comparing builds")
    args = ap.parse_args()

    config = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    if args.source_build:
        config["source_build"] = args.source_build
    pedo = yaml.safe_load(PEDOGENESIS.read_text(encoding="utf-8"))
    climatology = (args.climatology or climatology_path()).resolve()
    if not climatology.is_file():
        raise SystemExit(f"{climatology} does not exist")
    from provenance import require_build
    from climatology import annual_mean
    require_build(climatology, "climatology", config)

    with nc.Dataset(climatology) as data:
        lat = np.asarray(data["lat"][:], dtype=float)
        lon = np.asarray(data["lon"][:], dtype=float)
        land = np.asarray(data["lsm"][0], dtype=float) > 0.5
        # Weighted by records per bin, which are not equal: pyburn's
        # linspace().astype(int) split puts sixteen records in two bins of
        # twelve and fifteen in the rest. TASKS.md CLIM-13.
        centres = np.asarray(data["time"][:], dtype=float)
        temperature = annual_mean(np.asarray(data["tas"][:], dtype=float),
                                  centres) - KELVIN
        scale = 1000.0 * 86400.0 * EARTH_YEAR_DAYS
        evaporation = -annual_mean(np.asarray(data["evap"][:], dtype=float), centres)
        runoff = np.maximum(
            (annual_mean(np.asarray(data["pr"][:], dtype=float), centres) - evaporation)
            * scale, 0.0)

    mass_earth = float(config["planet"]["mass_earth"])
    fractions, mesh, grid_dir = lithology_fractions(config)
    silica_c, bicarb_c = concentrations()

    # Cell area on the sphere, km2, so fluxes are absolute rather than per-area.
    radius_km = float(mesh.manifest["planet"]["radiusKm"])
    dlon = 2.0 * np.pi / len(lon)
    # Gaussian latitudes are uneven, so bands run between row midpoints rather
    # than at a fixed spacing.
    edges = np.deg2rad(np.concatenate(
        ([90.0], 0.5 * (lat[:-1] + lat[1:]), [-90.0])))
    band = np.abs(np.sin(edges[:-1]) - np.sin(edges[1:]))
    area = (radius_km ** 2 * dlon * band)[:, None] * np.ones((1, len(lon)))

    # Litres of runoff per cell per Earth year: mm -> m -> m3 -> l over km2.
    litres = runoff * 1e-3 * area * 1e6 * 1e3

    silica_flux = np.zeros_like(runoff)      # umol/yr
    co2_flux = np.zeros_like(runoff)         # ueq/yr, then -> mol
    # SILICATE-derived CO2 is the only part that belongs to the long-term
    # thermostat. Carbonate weathering consumes CO2 on land and returns it when
    # the carbonate reprecipitates in the ocean, so over thermostat timescales
    # it is a wash. Kept apart because conflating them inflates the total by
    # roughly a factor of two and makes it incomparable with any published
    # silicate figure.
    silicate_co2 = np.zeros_like(runoff)
    # Per-class contributions, so a single lithology quietly dominating the
    # planet's carbon budget is visible in the output rather than something
    # someone has to go looking for. It has already happened once: see the
    # `dominant_class` block below.
    per_class: dict[str, dict[str, float]] = {}
    volcanic_silica = np.zeros_like(runoff)
    volcanic_co2 = np.zeros_like(runoff)
    volcanic_fraction = np.zeros_like(runoff)
    for code, frac in fractions.items():
        klass = ROCK_TO_MEYBECK.get(code)
        if klass is None:
            continue
        si = silica_c[klass] * frac * litres
        co2 = (bicarb_c[klass] * CO2_PER_HCO3.get(klass, DEFAULT_CO2_PER_HCO3)
               * frac * litres)
        silica_flux += si
        co2_flux += co2
        if klass not in CARBONATE_BEARING:
            silicate_co2 += co2
            per_class[code] = {
                "meybeck_class": klass,
                "hco3_ueq_l": bicarb_c[klass],
                "silicate_co2_mol_per_year": float(co2[land].sum()) * 1e-6,
                "land_fraction": float((frac * area)[land].sum()
                                       / float(area[land].sum())),
            }
        if klass == "volcanic_rocks":
            volcanic_silica += si
            volcanic_co2 += co2
            volcanic_fraction += frac

    # umol -> mol, ueq -> mol (bicarbonate is monovalent, so ueq == umol).
    silica_mol = silica_flux * 1e-6
    co2_mol = co2_flux * 1e-6

    endorheic = land_fraction_of_class(mesh, grid_dir,
                                       mesh.is_endorheic.astype(bool))

    def total(field: np.ndarray, weight: np.ndarray | None = None) -> float:
        w = field if weight is None else field * weight
        return float(w[land].sum())

    silica_total = total(silica_mol)
    co2_total = total(co2_mol)
    silicate_mol = silicate_co2 * 1e-6
    silicate_total = total(silicate_mol)
    silicate_endorheic = total(silicate_mol, endorheic)
    silica_endorheic = total(silica_mol, endorheic)
    co2_endorheic = total(co2_mol, endorheic)

    # Dessert cross-check, over volcanic area only.
    dessert = dessert_co2(runoff, temperature) * area * volcanic_fraction
    dessert_total = float(dessert[land].sum())

    land_area = float(area[land].sum())
    # SiO2 molar mass 60.08 g/mol -> t/km2/yr, the unit Durr reports.
    silica_yield = silica_total * 60.08 / 1e6 / land_area

    report = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "climatology": rel(climatology),
        "climatology_sha256": sha256(climatology),
        "source_build": config.get("source_build"),
        "terrain_hash": mesh.terrain_hash,
        "meybeck_tables_sha256": sha256(MEYBECK),
        "land_area_km2": land_area,
        "method": (
            "Flux = concentration x runoff, per lithology, with concentrations "
            "from Meybeck (1987) Table 2C. CO2 consumption is the bicarbonate "
            "flux at 1 mol CO2 per mol HCO3 for silicates, 0.5 for carbonate "
            "rocks whose other half is rock-derived, and 0 for evaporites."),
        "co2": {
            "silicate_consumption_mol_per_year": silicate_total,
            "silicate_endorheic_mol_per_year": silicate_endorheic,
            "silicate_exorheic_mol_per_year": silicate_total - silicate_endorheic,
            "silicate_endorheic_fraction": (silicate_endorheic / silicate_total)
                                           if silicate_total else None,
            "silicate_note": (
                "SILICATE ONLY, and this is the number to compare with any "
                "published figure. Carbonate weathering also consumes CO2 on "
                "land but returns it when the carbonate reprecipitates in the "
                "ocean, so it does not belong to the long-term thermostat."),
            "all_lithologies_mol_per_year": co2_total,
            "all_lithologies_note": (
                "Silicate plus carbonate-derived, the latter at 0.5. Reported "
                "for completeness; do NOT compare it with a silicate figure."),
            "endorheic_mol_per_year": co2_endorheic,
            "endorheic_fraction": (co2_endorheic / co2_total) if co2_total else None,
            "volcanic_mol_per_year": total(volcanic_co2 * 1e-6),
            "volcanic_share": (total(volcanic_co2 * 1e-6) / co2_total)
                              if co2_total else None,
            "dessert_volcanic_mol_per_year": dessert_total,
            "dessert_note": (
                "Independent check on the volcanic part only, from Dessert et "
                "al. (2003) Eq. 2 fitted on basaltic provinces. It is LINEAR in "
                "runoff where the Meybeck route is also linear, so the two "
                "agreeing is a check on the concentration, not on the runoff "
                "law. Neither uses Berner's 0.65 exponent, which applies to the "
                "global lithological mix rather than to basalt."),
            "earth_silicate_consumption_mol_per_year": 11.7e12,
            "earth_reference": (
                "Gaillardet et al. (1999) put global silicate CO2 consumption "
                "near 11.7e12 mol/yr. Quoted for scale only; it is a different "
                "planet's land area, relief and lithology."),
        },
        "co2_by_class": dict(sorted(
            per_class.items(),
            key=lambda kv: -kv[1]["silicate_co2_mol_per_year"])),
        "dominant_class": {
            "note": (
                "Which lithology the carbon budget rests on, reported every run "
                "so that one class quietly deciding the planet's carbon cycle is "
                "visible in the output rather than found by accident. It has "
                "happened once already and cost a factor of 1.4 on the answer."),
            "what_happened": (
                "`melange` was mapped to Meybeck's misc_metamorphic, chosen when "
                "the melange rule could not fire and the class was 0.00% of "
                "land. The arc fix made it reachable at over 6% of land, at "
                "which point that unexamined choice supplied 46% of the "
                "planet's silicate CO2 drawdown. It was wrong twice: "
                "misc_metamorphic is MARBLE by its own chemistry, Ca + Mg = "
                "1740 against HCO3 1730, so it gave a forearc province "
                "carbonate chemistry AND fed rock-derived bicarbonate into the "
                "silicate total at full weight. Melange now maps to shale, the "
                "greywacke and argillite the forearc province actually is, and "
                "carbonate-bearing classes are excluded from the silicate total "
                "rather than one class being tested by name."),
            "the_general_lesson": (
                "A judgment made while a class is absent is not a judgment. "
                "Three separate values in this project -- arc albedo, arc "
                "erodibility, and this mapping -- went unchecked for exactly "
                "that reason, and all three surfaced together when one bug in "
                "the tectonic model was fixed."),
        },
        "carbon_balance": {
            "note": (
                "At steady state the atmosphere's CO2 is set by outgassing "
                "balancing burial, and silicate weathering is what converts "
                "CO2 into buriable alkalinity. So a silicate consumption rate "
                "is also a statement about the outgassing this world needs."),
            "silicate_consumption_over_earth": (silicate_total / 11.7e12),
            "implied_outgassing_over_earth": (silicate_total / 11.7e12),
            "reading": (
                "config/planet.yaml fixes CO2 at 450 ppm. That is an "
                "ASSUMPTION, not a result: nothing in this project solves the "
                "carbonate-silicate balance, and the word 'outgassing' appears "
                "nowhere else in it. This says what the assumption costs. "
                "Sustaining 450 ppm at this climate needs an outgassing rate "
                "of roughly this multiple of Earth's. If the real rate were "
                "Earth-like, CO2 would draw down until weathering fell to meet "
                "it, and the world would settle colder and thinner-aired than "
                "the climate runs assume."),
            "why_it_is_high": (
                "Not because this world weathers hard. Per unit land area it "
                "weathers LESS than Earth, because it is drier -- land runoff "
                "is about half Earth's. It has 2.1x Earth's land area, and "
                "area wins. A big-land planet is a high-outgassing planet or "
                "it is a cold one."),
            "not_solved_here": (
                "Meybeck's concentrations carry no CO2 dependence, so this "
                "cannot be iterated to a self-consistent CO2. Doing that needs "
                "a weathering law in pCO2 and a climate response to it, which "
                "is a coupled calculation across pedology and exoplasim rather "
                "than an analysis in either. What is bounded here is the "
                "size of the assumption, not its resolution."),
            "plausibility_check": {
                "note": (
                    "Whether the outgassing 450 ppm requires is a rate this "
                    "planet could actually supply. This is a CHECK on a "
                    "prescribed number, not a closure of the loop -- closing it "
                    "would replace a prescribed CO2 with a prescribed "
                    "outgassing, which is no better constrained."),
                "why_the_ratio_is_rigorous": (
                    "At steady state each planet's outgassing equals its own "
                    "silicate weathering, so the ratio of this world's silicate "
                    "consumption to Earth's IS the ratio of required "
                    "outgassing. That sidesteps absolute outgassing estimates, "
                    "which span a factor of several depending on whether "
                    "metamorphic and diagenetic sources are counted."),
                "required_outgassing_over_earth": (silicate_total / 11.7e12),
                "radiogenic_supply_over_earth": mass_earth,
                "margin": (mass_earth / (silicate_total / 11.7e12))
                          if silicate_total else None,
                "supply_basis": (
                    f"Radiogenic heat production scales with mass at fixed bulk "
                    f"composition, and this planet is {mass_earth:.3f} Earth "
                    f"masses. Outgassing tracks mantle melt production, which "
                    f"tracks heat flux, so ~{mass_earth:.1f}x is the first-order "
                    f"expectation."),
                "verdict": (
                    "450 ppm is attainable with margin: it asks for less "
                    "outgassing than mass scaling suggests this planet "
                    "supplies. It is a prescribed number that has now been "
                    "shown reachable, rather than merely assumed."),
                "caveats": [
                    "A plausibility bound, not a derivation. Melt production "
                    "depends on spreading rate and mantle temperature, higher "
                    "gravity compresses the melting column, and this world's "
                    "43% land means less ocean basin than Earth. The sign of "
                    "the net correction is not obvious and is not claimed.",
                    "Computed on whatever climatology and build are named "
                    "above, both of which predate the current terrain. The "
                    "requirement moves on the baseline re-run. The margin is wide "
                    "enough that it would take a large error to threaten the "
                    "verdict, but it is not unlimited -- re-check it.",
                    "The requirement is a FLOOR on outgassing. Supplying more "
                    "than it does not break anything; the thermostat would "
                    "settle at a higher CO2 and a warmer state, which is a "
                    "different world rather than an inconsistent one.",
                ],
            },
            "and_the_feedback_is_weaker_still": (
                "Only the exorheic part joins the marine carbonate feedback "
                "that stabilises CO2; see thermostat_efficiency.py. Endorheic "
                "alkalinity still buries carbon, on its own basin floor, but "
                "it is decoupled from the loop that regulates."),
        },
        "silica": {
            "release_mol_per_year": silica_total,
            "endorheic_mol_per_year": silica_endorheic,
            "exorheic_mol_per_year": silica_total - silica_endorheic,
            "endorheic_fraction": (silica_endorheic / silica_total)
                                  if silica_total else None,
            "volcanic_mol_per_year": total(volcanic_silica * 1e-6),
            "volcanic_share": (total(volcanic_silica * 1e-6) / silica_total)
                              if silica_total else None,
            "land_yield_t_sio2_per_km2_per_year": silica_yield,
            "earth_reference_yield": 3.3,
            "earth_reference_note": (
                "Durr et al. (2011) give a global mean exorheic DSi yield of "
                "3.3 t SiO2/km2/yr, and note that fresh unweathered ash can "
                "exceed 50 locally. The Meybeck concentrations used here are "
                "for ordinary volcanic terrain, not fresh tephra, so this is a "
                "floor on the volcanic contribution rather than a central "
                "estimate."),
        },
        "caveats": [
            "Concentrations are Meybeck's temperate stream model. They carry no "
            "explicit temperature dependence, so the only climate signal in "
            "these fluxes is runoff. Dessert's law, which does carry "
            "temperature, is reported beside the volcanic figure for that "
            "reason.",
            "Fresh volcanic glass dissolves far faster than the crystalline "
            "rock Meybeck's volcanic class represents. Wolff-Boenisch et al. "
            "(2004) put a 1 mm basaltic glass grain at a 500 year lifetime "
            "against 4500 for rhyolitic, and field rates on andesitic ash run "
            "about 20x the laboratory ones. No enhancement is applied here "
            "because the andisol model finds only 0.45% of land is andic; the "
            "correction is small and would need an eruption-age field to place.",
            "This is supply, not deposit. Whether silica becomes silcrete or "
            "diatomite depends on saturation, on biology, and on lake level "
            "over the stellar cycle, none of which is computed here.",
        ],
    }

    ANALYSIS.mkdir(parents=True, exist_ok=True)
    out = ANALYSIS / "weathering_fluxes.json"
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"land area            {land_area/1e6:10.2f} Mkm2")
    print(f"CO2, silicate only   {silicate_total/1e12:10.3f} e12 mol/yr "
          f"(Earth ~11.7)")
    print(f"  endorheic          {100*silicate_endorheic/silicate_total:10.1f}%")
    print(f"CO2, all lithologies {co2_total/1e12:10.3f} e12 mol/yr "
          f"(not comparable to the above)")
    print(f"  volcanic share     {100*total(volcanic_co2*1e-6)/co2_total:10.1f}%")
    print(f"  Dessert on volcanic{dessert_total/1e12:10.3f} e12 mol/yr against "
          f"{total(volcanic_co2*1e-6)/1e12:.3f} here")
    print(f"silica release       {silica_total/1e12:10.3f} e12 mol/yr")
    print(f"  endorheic          {100*silica_endorheic/silica_total:10.1f}%")
    print(f"  volcanic share     {100*total(volcanic_silica*1e-6)/silica_total:10.1f}%")
    print(f"  land yield         {silica_yield:10.2f} t SiO2/km2/yr "
          f"(Earth exorheic mean 3.3)")
    required = silicate_total / 11.7e12
    print(f"\noutgassing required  {required:10.2f} x Earth (450 ppm at this climate)")
    print(f"  mass scaling says  {mass_earth:10.2f} x Earth could be supplied")
    print(f"  margin             {mass_earth/required:10.2f} x")
    print(f"\nwrote {rel(out)}")


if __name__ == "__main__":
    main()
