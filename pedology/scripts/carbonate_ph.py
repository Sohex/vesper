#!/usr/bin/env python3
"""The alkaline end of the soil pH block, solved rather than transcribed.

    python pedology/scripts/carbonate_ph.py            # validate, then tabulate
    python pedology/scripts/carbonate_ph.py --check    # validation only

Vesper is a fictional super-Earth; every number below is either a property of
water and carbonic acid or a property of this world's declared atmosphere.

`pedology/config/pedogenesis.yaml` used to carry the calcite-buffered pH, the
two parent brackets that hang off it and the endorheic alkalinity bracket as
literals. Each of those is a function of ONE input, `config/planet.yaml`'s
`pCO2_bar`, and the file that carried them states in its own comments that
planet.yaml owns that number and nothing there restates it. A pH written into
the config IS that restatement, one function application removed, and
`build_soil.py` consumed it on every run. This module is the function, so the
config states no number at all and the consumer reads what is emitted.

## The equilibrium

Slessarev, Lin, Bingham, Johnson, Dai, Schimel and Chadwick (2016), "Water
balance creates a threshold in soil pH at the global scale", Nature 540,
567-569, 10.1038/nature20139, Methods, "Calcite buffer", equation (6). It is the
charge balance of a solution exposed to calcite and open to an atmosphere at CO2
partial pressure `p`:

    2[Ca2+] + [H] = [HCO3-] + 2[CO3(2-)] + [OH-]

with the carbonate species set by Henry's law and the two carbonic acid
dissociations, and [Ca2+] set by the calcite solubility product through
[CO3(2-)]. Clearing denominators gives a quartic in the hydrogen ion activity:

    0 = H^4 * 2Ks/(K1 K2 KH p) + H^3 - H (K1 KH p + Kw) - 2 K1 K2 KH p

**The printed equation in the paper is not this, and the difference is a
typesetting loss rather than a disagreement.** As set, the last two terms read
`- H Kw K1 KH p - K1 K2 KH p`: the hydroxide and bicarbonate terms have been run
together into a product instead of a sum, and the factor of two that carbonate
carries into a charge balance has been dropped. Both are recovered by writing
the charge balance out.

**What the test below can and cannot settle between those forms**, said because
the difference matters and a check that only looks decisive is worse than none.
Solving as printed gives pH 8.76 at the paper's own pressure against the 8.2 the
paper publishes, so the test rejects it outright. The missing factor of two on
carbonate is worth 0.0012 pH at that pressure -- far inside the one decimal the
paper reports -- so the test is blind to it and the charge balance is the only
thing that settles it. The factor is here because two charges per carbonate ion
is arithmetic, not because a number came out better with it.

**The no-calcite limit is the same equation.** Drop the `H^4` term, which is the
one carrying `Ks`, and what is left is the charge balance of water and CO2 alone
with no source of alkalinity. That is the other end of the silicate bracket, and
it comes out of the same three lines rather than out of a second derivation.

## Solving it

By bisection on pH, not by a root finder over the quartic's coefficients. The
`H^4` coefficient runs to 1e10 while the constant term runs to 1e-21, so the
companion matrix of that polynomial is badly scaled; bisection has no such
problem and the root is guaranteed unique. Reading the coefficient signs left to
right gives `+ + 0 - -`, one sign change, so Descartes' rule says exactly one
positive real root exists. `--check` cross-checks the bisection against
`numpy.roots` anyway, because two methods agreeing is a test and one method
looking reasonable is not.

## Units

The equilibrium constants below are per ATMOSPHERE and `pCO2_bar` is in bar, so
the conversion is applied. It is worth 0.004 pH here, and it is done because the
constant's unit says to rather than because the answer needed it.
"""

from __future__ import annotations

import argparse
import math

# Equilibrium constants at 25 C, base-10 logarithms. Properties of water,
# carbonic acid and calcite; nothing in this project computes any of them, and
# they are the values Slessarev's Methods cite (Plummer and Busenberg for
# calcite, Sander's Henry's law compilation for KH).
LOG10_KH = -1.468        # Henry's constant for CO2, mol/l/atm
LOG10_K1 = -6.352        # first dissociation of carbonic acid, mol/l
LOG10_K2 = -10.329       # second dissociation of carbonic acid, mol/l
LOG10_KS = -8.48         # calcite solubility product, mol2/l2
KW = 1.0e-14             # dissociation constant of water, mol2/l2
TEMPERATURE_C = 25.0     # the temperature all five are stated at

# THE ONE PUBLISHED POINT THAT IS A TEST: Slessarev state pH 8.2 at 3.45e-4 atm,
# and they state both halves of it. The implementation has to land there before
# it is evaluated anywhere else.
PUBLISHED_PCO2_ATM = 3.45e-4
PUBLISHED_PH = 8.2
PUBLISHED_CONDITION = ("the laboratory pCO2 the measurements were made under, "
                       "the 1985 ambient mole fraction at standard pressure")
# The published value is given to one decimal, so agreeing to one decimal is the
# whole of what can be asked. Fixed on the paper's precision before the solver
# was run.
PUBLISHED_TOLERANCE = 0.05

# THE SECOND PUBLISHED FIGURE IS NOT A SECOND TEST POINT, and treating it as one
# is a mistake this file exists to have already made. The paper says the expected
# pH is 8.3 before 1977 and states NO pressure for it, so a test against 8.3 has
# to invent its input, and an invented input can be moved until the answer
# arrives. What the figure does say, at the one decimal it is printed to, is that
# the solved pH crosses 8.25 somewhere in the ambient CO2 of the decade before
# 1977. That IS checkable and it can fail: the crossing pressure must be below
# the 1985 laboratory value and above 3.0e-4 atm, which is where the ambient
# record sat through the 1960s and 1970s. The band is on the record and is fixed
# ahead of the solve.
ROUNDING_BOUNDARY_PH = 8.25
CROSSING_BAND_ATM = (3.0e-4, PUBLISHED_PCO2_ATM)

BAR_PER_ATM = 1.01325

# Soil air runs above atmospheric because roots and decomposers put CO2 into it
# faster than it diffuses out. The range is a DECLARED bracket on that
# enrichment rather than a measurement on this world, and it sets how far below
# the open-atmosphere value a field pH on carbonate can read.
SOIL_AIR_ENRICHMENT_BRACKET = (10.0, 100.0)

# Helvaci (2019) lists lake water at pH 8.5 to 11 among the conditions the
# Turkish borate deposits form under. The one measured closed-basin range this
# project holds, and the target the endorheic alkalinity bonus has to reach.
CLOSED_BASIN_PH_RANGE = (8.50, 11.00)


def _constants() -> tuple[float, float, float, float]:
    return (10.0 ** LOG10_KH, 10.0 ** LOG10_K1,
            10.0 ** LOG10_K2, 10.0 ** LOG10_KS)


def _residual(h: float, p_atm: float, with_calcite: bool) -> float:
    """The quartic, evaluated at hydrogen ion activity `h`.

    Strictly increasing in `h` over the positive reals for either value of
    `with_calcite`, which is what makes the bisection below safe.
    """
    kh, k1, k2, ks = _constants()
    quartic = 2.0 * ks / (k1 * k2 * kh * p_atm) * h ** 4 if with_calcite else 0.0
    return (quartic + h ** 3 - h * (k1 * kh * p_atm + KW)
            - 2.0 * k1 * k2 * kh * p_atm)


def _solve_ph(p_atm: float, with_calcite: bool) -> float:
    if not p_atm > 0.0:
        raise ValueError(f"CO2 partial pressure must be positive, got {p_atm}")
    lo, hi = 1.0e-14, 1.0          # pH 14 to pH 0; brackets any soil solution
    if _residual(lo, p_atm, with_calcite) > 0.0:
        raise ValueError(f"pH above 14 at pCO2 {p_atm} atm; outside the bracket")
    if _residual(hi, p_atm, with_calcite) < 0.0:
        raise ValueError(f"pH below 0 at pCO2 {p_atm} atm; outside the bracket")
    for _ in range(200):
        mid = math.sqrt(lo * hi)   # bisect in pH, which is where the tolerance is
        if _residual(mid, p_atm, with_calcite) > 0.0:
            hi = mid
        else:
            lo = mid
        if hi / lo < 1.0 + 1.0e-12:
            break
    return -math.log10(math.sqrt(lo * hi))


def calcite_equilibrium_ph(p_co2_atm: float) -> float:
    """pH of a solution exposed to calcite and open to an atmosphere at `p`."""
    return _solve_ph(p_co2_atm, with_calcite=True)


def co2_only_ph(p_co2_atm: float) -> float:
    """pH of water open to the same atmosphere with NO source of alkalinity."""
    return _solve_ph(p_co2_atm, with_calcite=False)


def atm_from_bar(p_co2_bar: float) -> float:
    """`config/planet.yaml`'s `pCO2_bar` in the units the constants are stated in."""
    return float(p_co2_bar) / BAR_PER_ATM


def rounding_boundary_pco2_atm() -> float:
    """The pCO2 at which the solved pH crosses `ROUNDING_BOUNDARY_PH`.

    Below this pressure the equilibrium pH prints as 8.3 at one decimal and
    above it as 8.2, which is the whole content of the paper's remark about
    pre-1977 measurements.
    """
    lo, hi = 1.0e-5, 1.0e-2
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if calcite_equilibrium_ph(mid) > ROUNDING_BOUNDARY_PH:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def check() -> list[str]:
    """Reproduce the paper's own published pH at the paper's own pressure.

    This is the test that can fail, and it does four things rather than one.
    Each is stated with what it catches, because a check whose failure mode is
    unnamed is not a check.
    """
    problems = []

    # 1. The published point. Catches the equation as printed, which lands on
    #    pH 8.76 here, and any constant that is wrong by more than a hundredth
    #    of a log unit.
    got = calcite_equilibrium_ph(PUBLISHED_PCO2_ATM)
    if abs(got - PUBLISHED_PH) > PUBLISHED_TOLERANCE:
        problems.append(
            f"calcite equilibrium at {PUBLISHED_PCO2_ATM:.3g} atm "
            f"({PUBLISHED_CONDITION}) gives pH {got:.4f}; Slessarev et al. "
            f"(2016) publish {PUBLISHED_PH}, and the tolerance is "
            f"{PUBLISHED_TOLERANCE}")

    # 2. Two methods on the same quartic. Catches the solver rather than the
    #    equation: a bracket that does not contain the root, a bisection that
    #    has not converged, a sign error in the residual. numpy is imported
    #    here rather than at module scope because the solver itself needs
    #    nothing but the standard library.
    import numpy as np
    kh, k1, k2, ks = _constants()
    p_atm = PUBLISHED_PCO2_ATM
    coefficients = [2.0 * ks / (k1 * k2 * kh * p_atm), 1.0, 0.0,
                    -(k1 * kh * p_atm + KW), -2.0 * k1 * k2 * kh * p_atm]
    positive = [r.real for r in np.roots(coefficients)
                if abs(r.imag) < 1.0e-25 and r.real > 0.0]
    if len(positive) != 1:
        problems.append(f"at {p_atm:.3g} atm the quartic has {len(positive)} "
                        "positive real roots, not one")
    else:
        by_roots = -math.log10(positive[0])
        if abs(by_roots - got) > 1.0e-6:
            problems.append(
                f"at {p_atm:.3g} atm bisection gives pH {got:.8f} and "
                f"numpy.roots gives {by_roots:.8f}; the two disagree")

    # 3. The rounding boundary the paper's pre-1977 remark implies. Catches a
    #    solver that is right at one pressure and wrong in its pressure
    #    dependence, which the single published point cannot see.
    low, high = CROSSING_BAND_ATM
    crossing = rounding_boundary_pco2_atm()
    if not low < crossing < high:
        problems.append(
            f"the solved pH crosses {ROUNDING_BOUNDARY_PH} at pCO2 "
            f"{crossing:.4e} atm, outside [{low:.3g}, {high:.3g}]. Slessarev "
            "et al. print 8.3 for pre-1977 measurements and 8.2 for the 1985 "
            "laboratory pressure, so the crossing belongs in the ambient CO2 "
            "of the decade before 1977")

    # 4. The bracket the config hangs on this is not inverted. Catches a sign
    #    convention that has flipped between the two limits.
    if not co2_only_ph(p_atm) < calcite_equilibrium_ph(p_atm):
        problems.append("the no-alkalinity end is not below the calcite end; "
                        "the bracket derived from these is inverted")
    return problems


def derived_ph_block(p_co2_bar: float) -> dict:
    """Every pH the carbonate system decides, at this world's pCO2.

    The one input is `config/planet.yaml`'s `pCO2_bar`. Returned to full
    precision; the caller rounds if it is writing the numbers down.
    """
    p_atm = atm_from_bar(p_co2_bar)
    low_enrichment, high_enrichment = SOIL_AIR_ENRICHMENT_BRACKET
    atmospheric = calcite_equilibrium_ph(p_atm)
    soil_air_low = calcite_equilibrium_ph(p_atm * high_enrichment)
    closed_low, closed_high = CLOSED_BASIN_PH_RANGE
    return {
        "pCO2_bar": float(p_co2_bar),
        "pCO2_atm": p_atm,
        "temperature_c": TEMPERATURE_C,
        "no_alkalinity_ph": co2_only_ph(p_atm),
        "calcite_atmospheric_ph": atmospheric,
        "calcite_soil_air_ph": {
            f"{int(low_enrichment)}x": calcite_equilibrium_ph(p_atm * low_enrichment),
            f"{int(high_enrichment)}x": soil_air_low,
        },
        "parent_carbonate": atmospheric,
        # The DRY end of the pH block, the buffer every parent relaxes onto as
        # drainage goes to zero. Numerically the carbonate parent, because both
        # are calcite saturation at this world's atmospheric pCO2, and named
        # separately because they are different roles: one is the pH of a
        # solution on carbonate rock, the other is where a soil that exports
        # nothing ends up whatever its rock was. The bracket is the same
        # equilibrium at the declared soil-air enrichment, so a field pH reads
        # below the open-atmosphere value and the declared value sits at the
        # top of its own bracket.
        "calcite_buffer_ph": atmospheric,
        "calcite_buffer_ph_bracket": [soil_air_low, atmospheric],
        "parent_bracket_silicate": [co2_only_ph(p_atm), atmospheric],
        "parent_bracket_carbonate": [soil_air_low, atmospheric],
        "endorheic_alkalinity_bonus_bracket": [closed_low - atmospheric,
                                               closed_high - atmospheric],
        "closed_basin_ph_range": list(CLOSED_BASIN_PH_RANGE),
        "soil_air_enrichment_bracket": list(SOIL_AIR_ENRICHMENT_BRACKET),
        "source": "Slessarev et al. (2016) Nature 540, 567-569, "
                  "10.1038/nature20139, Methods eq. (6), as the charge balance "
                  "rather than as printed; pedology/scripts/carbonate_ph.py",
    }


def resolve(ph_params: dict, p_co2_bar: float) -> tuple[dict, dict]:
    """Fill the `derived` sentinels in a `ph` block, and refuse a restatement.

    Returns `(params, provenance)`: a copy of the block with every derived key
    carrying a number, and the full-precision derivation for the report.

    **A literal where the sentinel belongs is an error rather than an
    override.** That is the whole point of the sentinel: a number written into
    the config is a second statement of this module's arithmetic, silently
    staleable against a change of `pCO2_bar`, and the config's own comments
    already say that `config/planet.yaml` owns that number and nothing here
    restates it. `config/planet.yaml`'s own `derived` keys work the same way.
    """
    problems = check()
    if problems:
        raise SystemExit("carbonate_ph.py does not reproduce Slessarev's "
                         "published pH; refusing to derive anything from it:\n  "
                         + "\n  ".join(problems))
    block = derived_ph_block(p_co2_bar)
    resolved = dict(ph_params)
    for key, value in (("parent_bracket_silicate", block["parent_bracket_silicate"]),
                       ("parent_bracket_carbonate", block["parent_bracket_carbonate"]),
                       ("calcite_buffer_ph_bracket",
                        block["calcite_buffer_ph_bracket"]),
                       ("endorheic_alkalinity_bonus_bracket",
                        block["endorheic_alkalinity_bonus_bracket"])):
        declared = ph_params.get(key)
        if declared != "derived":
            raise SystemExit(
                f"pedogenesis.yaml ph.{key} is {declared!r}; it must be the "
                "string `derived`. It is a function of config/planet.yaml's "
                "pCO2_bar through Slessarev's calcite equilibrium, and a "
                "number here is a restatement that cannot learn pCO2 moved. "
                "See pedology/scripts/carbonate_ph.py")
        resolved[key] = list(value)

    if ph_params.get("calcite_buffer_ph") != "derived":
        raise SystemExit(
            f"pedogenesis.yaml ph.calcite_buffer_ph is "
            f"{ph_params.get('calcite_buffer_ph')!r}; it must be the string "
            "`derived`. It is the dry end of the pH block, calcite saturation "
            "at this world's pCO2, and a number here is a restatement that "
            "cannot learn pCO2 moved. See pedology/scripts/carbonate_ph.py")
    resolved["calcite_buffer_ph"] = block["calcite_buffer_ph"]

    # Every entry but the evaporite is the calcite buffer: a soil that exports
    # nothing accumulates pedogenic calcite until the solution saturates, so
    # which buffer it lands on is not a property of the rock. Only the evaporite
    # sits above calcite saturation, and it is declared.
    parents = dict(ph_params["parent_by_category"])
    stated = [name for name, value in parents.items()
              if name != "evaporite" and value != "derived"]
    if stated:
        raise SystemExit(
            "pedogenesis.yaml ph.parent_by_category states a number for "
            + ", ".join(f"{name} ({parents[name]!r})" for name in sorted(stated))
            + "; every entry but `evaporite` must be the string `derived`. The "
              "buffer a soil reaches with nothing exported is the calcite "
              "equilibrium at this world's pCO2, and a number here is a "
              "restatement that cannot learn pCO2 moved. See "
              "pedology/scripts/carbonate_ph.py")
    for name in parents:
        if name != "evaporite":
            parents[name] = block["parent_carbonate"]
    resolved["parent_by_category"] = parents

    # The check that can fail on a value nobody edited. The declared quantity is
    # a base-cation supply and not a pH, so the bracket cannot be applied to it
    # directly; what the bracket bounds is the fresh solution that supply
    # implies. Bicarbonate carries the alkalinity, so a solution supplying a
    # fraction `u` of a calcite-saturated one sits at
    # calcite_ph + log10(u), which is above water in equilibrium with the
    # atmosphere and below calcite saturation exactly when `u` is a supply a
    # soil solution can have. A change of pCO2 moves both ends without moving
    # any declared supply.
    supplies = ph_params["base_cation_supply_by_category"]
    brackets = ph_params["base_cation_supply_bracket_by_category"]
    low, high = block["parent_bracket_silicate"]
    calcite = block["calcite_buffer_ph"]
    bad = []
    for name, value in sorted(supplies.items()):
        value = float(value)
        if not 0.0 < value <= 1.0:
            bad.append(f"{name} {value}, which is not a fraction of a "
                       "calcite-saturated soil's supply")
            continue
        implied = calcite + math.log10(value)
        if not low <= implied <= high:
            bad.append(f"{name} {value} implies a fresh solution at pH "
                       f"{implied:.4f}")
        ends = brackets.get(name)
        if ends is not None and not float(ends[0]) <= value <= float(ends[1]):
            bad.append(f"{name} {value} is outside its own declared bracket "
                       f"[{ends[0]}, {ends[1]}]")
    if bad:
        raise SystemExit(
            "pedogenesis.yaml ph.base_cation_supply_by_category does not "
            "describe soil solutions the carbonate system allows at pCO2 "
            f"{p_co2_bar} bar, where the silicate bracket is "
            f"[{low:.4f}, {high:.4f}]: " + "; ".join(bad))
    return resolved, block


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Solve the calcite equilibrium of Slessarev et al. (2016) "
                    "Methods eq. (6). Validates against their published pH at "
                    "their stated pressure, then tabulates this world's.")
    parser.add_argument("--check", action="store_true",
                        help="run the validation only and say nothing else")
    parser.add_argument("--pco2-bar", type=float, default=None,
                        help="override config/planet.yaml's pCO2_bar")
    args = parser.parse_args()

    problems = check()
    print(f"pCO2 {PUBLISHED_PCO2_ATM:.3g} atm  solved pH "
          f"{calcite_equilibrium_ph(PUBLISHED_PCO2_ATM):.4f}  published "
          f"{PUBLISHED_PH}   {PUBLISHED_CONDITION}")
    crossing = rounding_boundary_pco2_atm()
    print(f"pH {ROUNDING_BOUNDARY_PH} crossing at {crossing:.4e} atm "
          f"({crossing * 1e6:.1f} ppmv at one atmosphere total), so the "
          f"one-decimal figure reads 8.3 below it; band "
          f"[{CROSSING_BAND_ATM[0]:.3g}, {CROSSING_BAND_ATM[1]:.3g}]")
    if problems:
        raise SystemExit("VALIDATION FAILED:\n  " + "\n  ".join(problems))
    print(f"validation passes at {PUBLISHED_TOLERANCE} pH")
    if args.check:
        return

    from pathlib import Path
    import yaml
    p_bar = args.pco2_bar
    if p_bar is None:
        config = yaml.safe_load(
            (Path(__file__).resolve().parents[2] / "config" / "planet.yaml")
            .read_text(encoding="utf-8"))
        p_bar = config["atmosphere"]["pCO2_bar"]
    block = derived_ph_block(p_bar)
    print(f"\nat pCO2 {block['pCO2_bar']} bar ({block['pCO2_atm']:.6e} atm):")
    print(f"  no alkalinity, water and CO2 alone   {block['no_alkalinity_ph']:.4f}")
    print(f"  calcite, open to the atmosphere      {block['calcite_atmospheric_ph']:.4f}")
    for label, value in block["calcite_soil_air_ph"].items():
        print(f"  calcite, soil air at {label:>4}            {value:.4f}")
    lo, hi = block["parent_bracket_silicate"]
    print(f"  parent_bracket_silicate              [{lo:.4f}, {hi:.4f}]")
    lo, hi = block["parent_bracket_carbonate"]
    print(f"  parent_bracket_carbonate             [{lo:.4f}, {hi:.4f}]")
    lo, hi = block["endorheic_alkalinity_bonus_bracket"]
    print(f"  endorheic_alkalinity_bonus_bracket   [{lo:.4f}, {hi:.4f}]")


if __name__ == "__main__":
    main()
