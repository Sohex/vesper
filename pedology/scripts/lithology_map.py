"""The one place Orogen's rock classes are mapped onto a weathering source.

Three tables in this component read the same twenty Orogen rock classes and
answer the same physical question about each of them -- how much acid-
neutralising capacity does this rock release into the water draining it -- for
three consumers:

  * `brine_paths.py` wants the solute chemistry of a basin's inflow, and reads
    Meybeck (1987) Table 2C directly through `ROCK_TO_MEYBECK`.
  * `build_soil.py:soil_ph` wants a base-cation supply relative to a
    calcite-saturated soil's, and reads `PH_GROUP` into
    `pedogenesis.yaml:ph.base_cation_supply_by_category`, whose values are that
    same Meybeck bicarbonate column normalised by the carbonate row.
  * `weathering_schemes.yaml:class_mapping` wants a rokgem lithology class, and
    the two rokgem schemes are fitted to river CO2 consumption, which is the
    same flux again.

**They are one question with three vocabularies, so they can be checked against
each other, and until 2026-08-30 they were not.** Two classes had drifted, and
in both the pH bin was the odd one out:

  * `melange` was metamorphic there and shale to the other two, a factor of 4.29
    in the supply on six per cent of this world's land area;
  * `playa_clastic` was an evaporite there, supply 1 by construction, and shale
    to the other two, a factor of 5.5 on fourteen per cent.

Nothing in the tree could see either. This module holds the tables so a future
divergence is a failing check rather than a discovery.

The check is `check()`; `main()` runs it. `build_soil.py:soil_ph` calls it on
every soil build, so a mapping edit cannot reach a soil map unexamined, and
`scripts/smoke_test.py` runs it per commit because it is a static read.

WHAT THE CHECK CANNOT DO is decide which reading is right. It compares the
supply a class's Meybeck row implies against the bracket its pH category
declares, and a class outside its category's bracket must appear in
`DECLARED_DISAGREEMENTS` with the argument for reading it two ways. A
declaration for a class that in fact agrees also fails, so a stale entry cannot
sit here unnoticed.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from _paths import PEDOGENESIS, WEATHERING_SCHEMES  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
MEYBECK = ROOT / "pedology" / "data" / "reference" / "meybeck1987_tables.json"


# ---------------------------------------------------------------------------
#  Orogen rock class -> Meybeck table 2C lithology
# ---------------------------------------------------------------------------
# Orogen has 20 classes against Meybeck's 10, so this is a judgment map and it is
# written out in full rather than derived, because a silent default here would
# put every unmatched class on one side of the divide. Each entry says why.
ROCK_TO_MEYBECK = {
    "morb":                "volcanic_rocks",
    "oib":                 "volcanic_rocks",
    "flood_basalt":        "volcanic_rocks",
    "arc_basalt":          "volcanic_rocks",
    "arc_andesite":        "volcanic_rocks",
    "rift_bimodal":        "volcanic_rocks",
    "granite":             "granite",
    "granodiorite":        "granite",
    "gneiss":              "gneiss",
    # Meybeck's class is "gneiss and mica-schists"; schist belongs with it.
    "schist":              "gneiss",
    # Quartzite is metamorphic by origin but chemically a quartz sand: nearly
    # inert. Meybeck's misc_metamorphic is serpentinite/marble/amphibolite and
    # releases 1375 ueq/l Ca, which quartzite emphatically does not.
    "quartzite":           "sandstone",
    # Melange -> SHALE, corrected 2026-08-16, and the correction mattered more
    # than any other entry in this map.
    #
    # It was `misc_metamorphic`, chosen when the melange rule could not fire and
    # the class was 0.00% of land, so nothing rested on it. The arc fix made the
    # rule reachable and melange became one of the largest classes on this
    # world's land, at which point the choice was supplying close to half of
    # this planet's silicate CO2 drawdown.
    #
    # It was wrong twice over. Meybeck's misc_metamorphic is MARBLE: 2.3% of
    # Earth's outcrop, Ca 1375 ueq/l against sedimentary carbonate's 2560, and
    # Ca + Mg = 1740 against HCO3 1730 -- a near-exact carbonate balance, which
    # is the signature of carbonate dissolution rather than silicate weathering.
    # So it gave a forearc province marble chemistry, AND that bicarbonate was
    # entering the silicate total at full weight when most of it is
    # rock-derived.
    #
    # What Orogen's class actually is decides the replacement. It is not the
    # accretionary prism alone: `vendor/orogen/js/lithology.js` assigns it to
    # the whole forearc province -- prism, forearc basin and serpentinite
    # together, named in that order in the comment above the rule -- because at
    # ~15 km cells they cannot be separated. That province is greywacke,
    # argillite and clastic basin fill by volume, which is Meybeck's shale, and
    # is what shelf, foreland and pelagic clastics already map to.
    #
    # `gneiss` is not the alternative it looks like. A forearc is not felsic
    # crystalline basement, and mapping it there would swap one wrong rock for
    # another while happening to give a smaller number.
    #
    # SERPENTINITE IS THE MIXTURE THIS ROW DOES NOT RESOLVE, and it is the
    # third thing the generator names. Melanges carry it and Meybeck has a
    # peridotite row for exactly that rock. Nothing in Orogen sizes the
    # serpentinite share of a forearc -- the class is one basement branch with
    # no internal division, and the generator's own comment says the prism
    # cannot be separated from the basin at its resolution -- so the fraction is
    # UNDECLARED here rather than assumed, and the pH block's `melange` bracket
    # spans its full range instead. What that costs is small in this column:
    # peridotite's bicarbonate is 450 ueq/l against shale's 580, so an
    # all-serpentinite melange supplies 0.1408 against an all-shale one's
    # 0.1815 and the whole unresolved mixture spans 22 per cent, against the
    # factor of 4.29 the metamorphic reading was wrong by.
    #
    # The bias directions, since intuition gets one of them backwards.
    # Meybeck's peridotite is LOWER in bicarbonate than shale, 450 against 580,
    # so reading the province as shale biases CO2 and Ca HIGH, not low. It
    # biases Mg (500 against 240) and silica (180 against 150) LOW. Ultramafic
    # rock being an efficient CO2 sink per unit area is a statement about
    # reaction rate, not about the solute concentration this table carries.
    "melange":             "shale",
    # "Shelf sandstone / shale" is a mixture. Shale is far more reactive and
    # dominates the solute load of any such mixture, so it is mapped there;
    # brine_paths' --clastic-as-sandstone tests the other choice.
    "shelf_clastic":       "shale",
    "carbonate":           "sedimentary_carbonate_rocks",
    "foreland_clastic":    "shale",
    "continental_clastic": "sandstone",
    "pelagic":             "shale",
    # Orogen's evaporite class is a SALT crust, so halite rather than gypsum.
    # There is no gypsum lithology in this world's rock table at all, which is
    # itself a result: see the note in brine_paths' output.
    "evaporite":           "halite_evaporite",
    # Playa mud is detrital fill, not an evaporite. It is the host that
    # evaporites grow in, and mapping it to an evaporite would beg the question.
    "playa_clastic":       "shale",
}


# ---------------------------------------------------------------------------
#  Orogen rock class -> pH supply category
# ---------------------------------------------------------------------------
# The bins `pedogenesis.yaml:ph.base_cation_supply_by_category` declares a
# supply for. Kept in code rather than in that config because it is a
# classification of the EXPORT's vocabulary and not a tunable: if Orogen adds a
# class this must fail loudly, which it does.
#
# They are CHEMICAL bins and not Orogen's own `category` field:
# quartzite is metamorphic by origin and sits in `igneous_felsic` here because
# it weathers like a quartz sand, and `melange` has a bin of its own because
# the mixture behind it is declared here rather than inherited from a source's
# outcrop weights. Reading Orogen's `category` instead is what put melange on
# gneiss chemistry until 2026-08-30.
PH_GROUP = {
    "morb": "igneous_mafic", "oib": "igneous_mafic",
    "flood_basalt": "igneous_mafic", "arc_basalt": "igneous_mafic",
    "arc_andesite": "igneous_felsic", "rift_bimodal": "igneous_felsic",
    "granite": "igneous_felsic", "granodiorite": "igneous_felsic",
    "gneiss": "metamorphic", "schist": "metamorphic",
    "quartzite": "igneous_felsic", "melange": "melange",
    "shelf_clastic": "sedimentary_clastic",
    "foreland_clastic": "sedimentary_clastic",
    "continental_clastic": "sedimentary_clastic",
    "pelagic": "sedimentary_clastic",
    "carbonate": "carbonate",
    # Playa mud is DETRITAL FILL and not an evaporite, which is the same call
    # ROCK_TO_MEYBECK and both rokgem class mappings already make. Orogen
    # zones the closed basin: `saltCrustMask` puts the salt crust in the SUMP,
    # "the part that repeatedly floods and dries", and gives the margins the
    # clastic load as playa mud and alluvial fans. So this class is by
    # definition the part of the basin where the dissolved load did NOT
    # precipitate, and putting it on the sump's buffer contradicts the
    # definition. The sump's own buffer is `ph.soda_buffer_ph`, carried by the
    # `evaporite` share of a cell in `build_soil.py:SUMP_GROUPS`.
    "evaporite": "evaporite", "playa_clastic": "sedimentary_clastic",
    "water": "sedimentary_clastic",
}


# ---------------------------------------------------------------------------
#  Where the two readings are allowed to differ, and why
# ---------------------------------------------------------------------------
# A class here is one whose Meybeck row implies a supply outside the bracket its
# pH category declares. Each entry is the argument for reading the same rock two
# ways. An entry whose class turns out to AGREE is a failure too: a declaration
# that has stopped being needed is a licence nobody is watching.
DECLARED_DISAGREEMENTS = {
    "arc_andesite":
        "Meybeck has ONE volcanic row and does not resolve andesite from "
        "basalt, so brine_paths reads andesite at basalt's solute chemistry "
        "because that is the only reading its source supports. The pH block's "
        "source set includes GEM-CO2, which HAS a felsic-volcanic class and "
        "puts andesite in it at 0.0725 -- which is `igneous_felsic`'s own "
        "bracket top. So this is not two answers to one question; it is one "
        "source resolving a distinction the other cannot. "
        "weathering_schemes.yaml makes the same split for the same reason.",
    "rift_bimodal":
        "As `arc_andesite`: bimodal rift volcanics are half felsic and Meybeck "
        "has no felsic-volcanic row to put them in.",
    "evaporite":
        "The pH block's evaporite supply is 1 BY CONSTRUCTION and not by "
        "measurement -- a parent that is itself the buffer mineral supplies "
        "what a soil buffered by it supplies -- so no Meybeck row is being "
        "read and none can disagree. Meybeck's halite row is a solute "
        "concentration on a different basis entirely: its cation sum is 6.19 "
        "TIMES the carbonate row's, which is not a fraction of anything.",
    "water": "Open water is not land and no soil is built on it. It has no "
             "Meybeck row; `PH_GROUP` carries it only so the mesh's own class "
             "legend can be walked without a lookup failing.",
}


# ---------------------------------------------------------------------------
#  rokgem lithology class -> Meybeck row
# ---------------------------------------------------------------------------
# So that `weathering_schemes.yaml:class_mapping` can be compared against
# ROCK_TO_MEYBECK at all. Both rokgem schemes carry a coarse class list;
# these are the Meybeck rows those classes stand for. `ice` and `undeclared`
# are not rocks and map to nothing.
ROKGEM_TO_MEYBECK = {
    "carb": "sedimentary_carbonate_rocks",
    "shale": "shale",
    "sand": "sandstone",
    "basalt": "volcanic_rocks",
    # GEM-CO2's crystalline-shield class. Meybeck's nearest row is gneiss and
    # mica-schist, which is what a shield is.
    "shield": "gneiss",
    # GEM-CO2's acid class is acid plutonic and volcanic together.
    "acid": "granite",
    # GKWM's granite row IS its crystalline row; there is no separate shield.
    "granite": "granite",
    "ice": None,
    "undeclared": None,
}

# Where a rokgem mapping is allowed to differ from ROCK_TO_MEYBECK, keyed
# (scheme, rock class). The reasons are argued in weathering_schemes.yaml.
DECLARED_SCHEME_DISAGREEMENTS = {
    ("gkwm", "gneiss"): "GKWM's granite row IS the crystalline shield row.",
    ("gkwm", "schist"): "GKWM's granite row IS the crystalline shield row.",
    # GKWM has no felsic-volcanic class either, so its andesite goes to basalt,
    # which is where Meybeck's one volcanic row already puts it. No exception is
    # needed; the two agree by both lacking the distinction.
    ("gem_co2", "arc_andesite"):
        "GEM-CO2 has a felsic-volcanic class and Meybeck does not.",
    ("gem_co2", "rift_bimodal"):
        "GEM-CO2 has a felsic-volcanic class and Meybeck does not.",
    ("gkwm", "evaporite"):
        "A halite crust is not a silicate weathering substrate and its "
        "dissolution consumes no CO2; no rokgem class represents it.",
    ("gem_co2", "evaporite"):
        "A halite crust is not a silicate weathering substrate and its "
        "dissolution consumes no CO2; no rokgem class represents it.",
}


def meybeck_supply_by_row(path: Path = MEYBECK) -> dict[str, float]:
    """Each Meybeck Table 2C row's bicarbonate, over the carbonate row's.

    The quantity `ph.base_cation_supply_by_category` declares. Bicarbonate and
    not the cation sum, because bicarbonate is the acid-neutralising capacity:
    base cations arriving paired with chloride or sulphate buffer nothing. At
    equal runoff a concentration ratio is a flux ratio, which is what makes a
    table of concentrations a table of supplies.
    """
    tables = json.loads(path.read_text(encoding="utf-8"))
    t2c = tables["table_2c"]
    column = t2c["columns"].index("hco3_ueq_l")
    rows = t2c["rows"]
    reference = float(rows["sedimentary_carbonate_rocks"][column])
    return {name: float(values[column]) / reference
            for name, values in rows.items()}


def check(ph_supplies: dict, ph_brackets: dict,
          class_mapping: dict | None = None) -> list[str]:
    """Every problem with the three mappings, as sentences. Empty means clean.

    `ph_supplies` and `ph_brackets` are `pedogenesis.yaml`'s
    `ph.base_cation_supply_by_category` and
    `ph.base_cation_supply_bracket_by_category`; `class_mapping` is
    `weathering_schemes.yaml`'s, read from that file when it is not passed.
    """
    problems: list[str] = []
    supply = meybeck_supply_by_row()

    missing = sorted(set(ROCK_TO_MEYBECK) - set(PH_GROUP))
    if missing:
        problems.append(
            "ROCK_TO_MEYBECK maps " + ", ".join(missing) + " and PH_GROUP does "
            "not; one rock table, two lists of it")
    unmapped = sorted(code for code in PH_GROUP
                      if code not in ROCK_TO_MEYBECK
                      and code not in DECLARED_DISAGREEMENTS)
    if unmapped:
        problems.append(
            "PH_GROUP carries " + ", ".join(unmapped) + " with no Meybeck row "
            "and no entry in DECLARED_DISAGREEMENTS, so nothing says what that "
            "rock releases")

    for code in sorted(PH_GROUP):
        category = PH_GROUP[code]
        row = ROCK_TO_MEYBECK.get(code)
        declared = DECLARED_DISAGREEMENTS.get(code)
        if category not in ph_supplies:
            problems.append(
                f"PH_GROUP puts {code} in {category!r}, which "
                "ph.base_cation_supply_by_category does not declare")
            continue
        if row is None:
            if declared is None:
                problems.append(f"{code} has no Meybeck row and no declaration")
            continue
        implied = supply[row]
        value = float(ph_supplies[category])
        ends = ph_brackets.get(category)
        if ends is None:
            agrees = abs(implied - value) <= 1e-9
            where = f"the category's declared supply {value:.4f}"
        else:
            agrees = float(ends[0]) <= implied <= float(ends[1])
            where = (f"the {category} bracket "
                     f"[{float(ends[0]):.4f}, {float(ends[1]):.4f}]")
        if agrees and declared is not None:
            problems.append(
                f"{code} is declared in DECLARED_DISAGREEMENTS but the two "
                f"mappings agree: its Meybeck row {row} implies "
                f"{implied:.4f}, inside {where}. Delete the declaration rather "
                "than leaving a licence nobody is watching")
        if not agrees and declared is None:
            problems.append(
                f"{code} reads as Meybeck {row}, supply {implied:.4f}, in "
                f"ROCK_TO_MEYBECK and as pH category {category!r} in PH_GROUP, "
                f"which is outside {where} -- a factor of "
                f"{implied / value:.2f}. One rock, two answers, neither aware "
                "of the other. Settle it, or record the argument for reading "
                "it two ways in DECLARED_DISAGREEMENTS")

    if class_mapping is None:
        import yaml
        class_mapping = yaml.safe_load(
            WEATHERING_SCHEMES.read_text(encoding="utf-8"))["class_mapping"]
    for scheme, mapping in sorted(class_mapping.items()):
        for code, rokgem_class in sorted(mapping.items()):
            if rokgem_class not in ROKGEM_TO_MEYBECK:
                problems.append(
                    f"weathering_schemes.yaml maps {code} to rokgem class "
                    f"{rokgem_class!r} under {scheme}, which "
                    "ROKGEM_TO_MEYBECK does not name")
                continue
            through = ROKGEM_TO_MEYBECK[rokgem_class]
            direct = ROCK_TO_MEYBECK.get(code)
            agrees = through == direct
            declared = DECLARED_SCHEME_DISAGREEMENTS.get((scheme, code))
            if agrees and declared is not None:
                problems.append(
                    f"({scheme}, {code}) is declared in "
                    "DECLARED_SCHEME_DISAGREEMENTS but the two mappings agree "
                    f"on {direct}")
            if not agrees and declared is None:
                problems.append(
                    f"weathering_schemes.yaml reads {code} as rokgem "
                    f"{rokgem_class!r} under {scheme}, which is Meybeck "
                    f"{through}, where ROCK_TO_MEYBECK reads it as {direct}. "
                    "That file says its reasoning follows ROCK_TO_MEYBECK; "
                    "either follow it or record the exception in "
                    "DECLARED_SCHEME_DISAGREEMENTS")
    return problems


def check_or_die(ph_supplies: dict, ph_brackets: dict) -> None:
    """`check`, raising. What `build_soil.py:soil_ph` calls."""
    problems = check(ph_supplies, ph_brackets)
    if problems:
        raise SystemExit(
            "the three mappings of Orogen's rock classes disagree; see "
            "pedology/scripts/lithology_map.py:\n  " + "\n  ".join(problems))


def main() -> None:
    import yaml
    ph = yaml.safe_load(PEDOGENESIS.read_text(encoding="utf-8"))["ph"]
    supplies = ph["base_cation_supply_by_category"]
    brackets = ph["base_cation_supply_bracket_by_category"]
    supply = meybeck_supply_by_row()
    print(f"{'rock class':22s} {'pH category':22s} {'Meybeck row':28s} "
          f"{'u_row':>7}  reading")
    for code in sorted(PH_GROUP):
        row = ROCK_TO_MEYBECK.get(code)
        implied = supply[row] if row else None
        note = "declared apart" if code in DECLARED_DISAGREEMENTS else "agrees"
        print(f"{code:22s} {PH_GROUP[code]:22s} {str(row):28s} "
              f"{'' if implied is None else f'{implied:7.4f}'}  {note}")
    problems = check(supplies, brackets)
    if problems:
        raise SystemExit("FAILED:\n  " + "\n  ".join(problems))
    print(f"\n{len(PH_GROUP)} rock classes; the three mappings agree or are "
          f"declared apart, {len(DECLARED_DISAGREEMENTS)} of them with a "
          "recorded argument")


if __name__ == "__main__":
    main()
