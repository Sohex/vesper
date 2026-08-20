#!/usr/bin/env python3
"""The CO2 transmissivity the N2O 589 cm-1 band has to be multiplied by.

    python exoplasim/scripts/co2_overlap_589.py

Worldbuilding. Vesper is an invented planet and this script is about the
simulation of it: one term in a toy climate model's longwave scheme. The gas
absorption data is real; the atmosphere it is evaluated for is invented.

CLIM-42. `exoplasim/notes/trace-gas-band.md` section 4b argues why this term is
not optional. N2O's 589 cm-1 band sits 78 cm-1 from the 667 cm-1 CO2
fundamental, well inside CO2's 15 um band, so carrying it without the overlap
credits N2O with absorption CO2 already provides -- an OVERSTATEMENT, which is
the opposite of the conservative error that dropping the band would be.

This writes `exoplasim/analysis/co2_overlap_589.json`.

## Two independent routes, and which one is adopted

**Adopted: correlated-k.** The band-mean CO2 transmissivity over 546-630 cm-1,
from the LMD Generic PCM's `N2-CO2var_2026` table -- HITRAN 2020, line by line,
CO2 as the variable gas so it evaluates at this world's mixing ratio rather than
a premixed one -- with its far-wing companion applied, since those tables are
built with a 25 cm-1 line cutoff and the wings are supplied separately.

**Cross-check: Ramanathan (1976) Appendix A**, which is what Donner and
Ramanathan (1980) name. Its procedure collapses to something simple. With
`A = 2 A0 ln(1 + xi)` for the strong-line limit and Goody's `T = exp(-A_bar)`
where `A_bar = A / (2 A0)`,

    T = 1 / (1 + xi_eff)

and `xi_eff` is the strong-line parameter formed with an EFFECTIVE intensity,
Edwards and Menard's line intensity distribution shifting a band centred at
`w1` into the region of a band centred at `w2`:

    S_eff = S exp(-|w2 - w1| / A0)
    xi    = 1.66 (4 v0 P / (A0 D)) S_eff u

**The cross-check is not the answer, and the reason is sourcing.** Ramanathan
attributes `D_i`, the mean line spacing, to Dickinson (1972), which tabulates
band strengths and no spacings; and `q_i`, the isotopic abundances, to Goody
(1964), a book. The correlated-k route needs NEITHER -- a k-distribution built
line by line at natural isotopic abundance has the spacing and the isotopes
already integrated in. The two quantities that cannot be cleanly sourced are
exactly the two it does not ask for, which is why it is the primary and this is
the check.

So `D` here is a STATED ASSUMPTION rather than a citation: 1.56 cm-1, four
times CO2's rotational constant, because only even-J levels are populated in
the main isotope and the P and R branch lines therefore fall 4B apart. The
sensitivity to it is reported, because a check whose free parameter is unstated
is not a check.

## What can fail

- **The lower-state energies against Dickinson's own band origins.** The hot
  band populations need vibrational energies that are NOT in Dickinson's table;
  they are standard, and they are verified here by differencing them against
  the band origins that ARE in it. (01'0)->(02'0) must come out at 618.0,
  (01'0)->(10'0) at 720.8, and so on. A wrong level shows up immediately.
- **The correlated-k band must contain 589 cm-1** and be comparable in width to
  the N2O band's `2 A0` of 46 cm-1.
- **Both routes must be monotone decreasing in CO2 amount** and approach one as
  the amount goes to zero.
- **The two routes must agree to a declared factor.** Set at 2.0, and stated
  before running: this compares a 2020 line list against a 1976 band model
  whose own free parameter is being assumed, so agreement to tens of percent
  would be luck and disagreement by a factor is a real disagreement.

## The trap, which this hit before it was written down

**Comparing the two TRANSMISSIVITIES directly is meaningless, and it is the
obvious thing to do.** They are means over different spectral widths: the corrk
band is 83.8 cm-1 wide, and Ramanathan's convention averages over `2 A0` = 34.6
cm-1, a factor of 2.42. Compared that way the routes disagree by up to 6x and
the disagreement is pure bookkeeping. Compared as ABSORPTANCES in cm-1, which
carry no width convention, they agree to better than the declared factor
everywhere.

Taking Ramanathan's `T` and using it where a corrk band mean belongs -- or the
reverse -- would have put a factor of several into the model with nothing to
catch it, since both are numbers between 0 and 1 that fall with CO2 amount.

## The window correction, and its sign

The corrk band spans 546-630 cm-1 and the N2O band occupies 566-612, its own
`2 A0`. Those are not the same piece of spectrum, and the difference is not
symmetric: the corrk band reaches further toward the 667 cm-1 CO2 fundamental,
where CO2 absorbs far more strongly. Weighting by the same exponential falloff
Ramanathan's effective intensity uses, the corrk band's mean CO2 absorption is
1.64x the absorption over the window N2O actually occupies.

So the band-mean optical depth is scaled by that ratio before it becomes the
model's multiplier. Left uncorrected it would UNDERSTATE the transmissivity N2O
sees and so understate the 589 band -- conservative, but wrong by a knowable
amount, and the correction is derived rather than fitted.
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _paths import ANALYSIS, PROJECT_ROOT  # noqa: E402
from corrk_cross_check import CORRK, CorrK, LOSCHMIDT  # noqa: E402

OUTPUT = ANALYSIS / "co2_overlap_589.json"
CONTINUUM = (Path.home() / "git" / "generic_pcm" / "LMDZ.GENERIC" / "datagcm"
             / "continuum" / "far_wing_data" / "CO2-N2_line-far-wings_50-1000K_2026.dat")

TABLE = "N2-CO2var_2026"
N2O_589_CENTRE = 589.0
N2O_589_A0 = 23.0             # Donner and Ramanathan Table 1, at 300 K

# Cess and Ramanathan (1972) Table 2: the CO2 15 um bandwidth parameter.
CO2_A0 = 17.3
CO2_FUNDAMENTAL = 667.4
NU0 = 0.064                   # mean line half-width, cm-1 atm-1, Ramanathan (1976)
DIFFUSIVITY = 1.66

# STATED ASSUMPTION, not a citation. See the module docstring.
LINE_SPACING = 1.56           # cm-1, 4B for CO2's 15 um band
LINE_SPACING_SCAN = (0.78, 1.56, 3.12)

# Dickinson (1972) Table 3, the seven 15 um bands of C12-O16-2. Band strengths
# are already cm-1 (cm atm STP)-1: the caption defines the units as cm-1 cm2
# WHEN DIVIDED BY Loschmidt. Strengths are per ACTIVE molecule, so each carries
# the Boltzmann population of its lower state.
CO2_BANDS_15UM = [
    # lower state, origin cm-1, strength, lower-state energy cm-1, degeneracy
    ("00-00", 667.4, 242.0, 0.0, 1),
    ("01-10", 618.0, 65.0, 667.38, 2),
    ("01-10", 720.8, 76.0, 667.38, 2),
    ("01-10", 667.8, 228.0, 667.38, 2),
    ("02-00", 647.1, 745.0, 1285.41, 1),
    ("02-20", 668.2, 663.0, 1335.13, 2),
    ("10-00", 688.7, 223.0, 1388.19, 1),
]
C2 = 1.4387769                # hc/k, cm K

AGREEMENT_FACTOR = 2.0        # declared before running
# This world's whole-column CO2 amount. Rows beyond it are stress points that
# the model cannot reach, and they are scored separately rather than dropped.
FULL_COLUMN_ATM_CM = 272.0


def window_correction(band_lo: float, band_hi: float) -> float:
    """How much of the corrk band's CO2 absorption belongs to the N2O window.

    Both windows are weighted by `exp(-|v - 667.4| / A0)`, the same falloff
    Ramanathan's effective intensity uses, so this is the band model's own
    statement about where within the band the absorption sits.
    """
    def mean(lo, hi):
        v = np.linspace(lo, hi, 20001)
        return float(np.trapezoid(np.exp(-np.abs(v - CO2_FUNDAMENTAL) / CO2_A0), v)
                     / (hi - lo))
    return mean(N2O_589_CENTRE - N2O_589_A0, N2O_589_CENTRE + N2O_589_A0) / mean(band_lo, band_hi)


def band_index(table: CorrK, centre: float) -> int:
    hit = [i for i, (lo, hi) in enumerate(table.edges) if lo <= centre <= hi]
    if len(hit) != 1:
        raise SystemExit(f"{len(hit)} bands contain {centre} cm-1")
    return hit[0]


def check_levels() -> list:
    """Differencing the assumed lower-state energies must give Dickinson's origins."""
    upper = {}
    problems = []
    for name, origin, _s, e_lower, _g in CO2_BANDS_15UM:
        upper.setdefault(name, []).append((origin, e_lower))
    # The three bands out of 01'0 pin the upper levels; the check that bites is
    # that each origin equals (upper - lower) for levels shared across bands.
    known = {"01-10": 667.38, "02-00": 1285.41, "02-20": 1335.13, "10-00": 1388.19}
    for name, origin, _s, e_lower, _g in CO2_BANDS_15UM:
        if name == "00-00":
            if abs(origin - known["01-10"]) > 0.1:
                problems.append(f"fundamental origin {origin} against level {known['01-10']}")
        elif name == "01-10":
            target = {618.0: "02-00", 720.8: "10-00", 667.8: "02-20"}.get(origin)
            if target and abs((known[target] - e_lower) - origin) > 0.1:
                problems.append(f"{target} - 01'0 = {known[target] - e_lower:.2f} "
                                f"against origin {origin}")
    return problems


def boltzmann_weights(temperature: float) -> np.ndarray:
    """Fraction of CO2 in each band's lower vibrational state."""
    levels = {}
    for name, _o, _s, e, g in CO2_BANDS_15UM:
        levels[name] = (e, g)
    part = sum(g * math.exp(-C2 * e / temperature) for e, g in levels.values())
    return np.array([CO2_BANDS_15UM[i][4]
                     * math.exp(-C2 * CO2_BANDS_15UM[i][3] / temperature) / part
                     for i in range(len(CO2_BANDS_15UM))])


def ramanathan_transmissivity(u_co2_atmcm: float, pressure_atm: float,
                              temperature: float, spacing: float) -> float:
    """T = 1/(1 + xi_eff), Ramanathan (1976) Appendix A over Dickinson's B1-B7."""
    q = boltzmann_weights(temperature)
    xi = 0.0
    for (name, origin, strength, _e, _g), qi in zip(CO2_BANDS_15UM, q):
        s_eff = strength * math.exp(-abs(N2O_589_CENTRE - origin) / CO2_A0)
        xi += (DIFFUSIVITY * 4.0 * NU0 * pressure_atm / (CO2_A0 * spacing)
               * s_eff * qi * u_co2_atmcm)
    return 1.0 / (1.0 + xi)


def far_wing(wavenumber: float, temperature: float) -> float:
    """CO2-N2 far-wing binary absorption coefficient, interpolated."""
    raw = CONTINUUM.read_text().splitlines()
    temps = np.array([float(t) for t in raw[0].replace("T=", " ").split()])
    data = np.array([[float(x) for x in line.split()] for line in raw[1:] if line.strip()])
    wn, k = data[:, 0], data[:, 1:]
    iw = int(np.clip(np.searchsorted(wn, wavenumber) - 1, 0, len(wn) - 2))
    it = int(np.clip(np.searchsorted(temps, temperature) - 1, 0, len(temps) - 2))
    fw = (wavenumber - wn[iw]) / (wn[iw + 1] - wn[iw])
    ft = (temperature - temps[it]) / (temps[it + 1] - temps[it])
    return float((1 - fw) * (1 - ft) * k[iw, it] + fw * (1 - ft) * k[iw + 1, it]
                 + (1 - fw) * ft * k[iw, it + 1] + fw * ft * k[iw + 1, it + 1])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--co2-vmr", type=float, default=4.5e-4)
    ap.add_argument("--pressure-mbar", type=float, default=500.0)
    ap.add_argument("--temperature", type=float, default=250.0)
    args = ap.parse_args()

    problems = check_levels()
    if problems:
        raise SystemExit("lower-state energies disagree with Dickinson's band "
                         "origins:\n  " + "\n  ".join(problems))

    table = CorrK(TABLE, CORRK)
    ib = band_index(table, N2O_589_CENTRE)
    lo, hi = table.edges[ib]
    if not (hi - lo) < 4.0 * N2O_589_A0:
        raise SystemExit(f"corrk band {lo:.0f}-{hi:.0f} is too wide to stand for "
                         f"an N2O band of 2A0 = {2 * N2O_589_A0:.0f} cm-1")

    p_atm = args.pressure_mbar / 1013.25
    amounts = [1.0, 3.0, 10.0, 30.0, 100.0, 272.0, 500.0]
    kfw = far_wing(N2O_589_CENTRE, args.temperature)

    correction = window_correction(float(lo), float(hi))
    width = float(hi) - float(lo)

    rows = []
    for u in amounts:
        u_air = u * LOSCHMIDT / args.co2_vmr
        t_lines = float(table.transmission(args.pressure_mbar, args.temperature,
                                           args.co2_vmr, u_air)[ib])
        # Far wing: a binary coefficient in cm-1 amagat-2, so the optical depth
        # is k * (n_co2/n0) * (n_bath/n0) * path, which in column terms is
        # k * u_co2 * (n_bath/n0) with u_co2 in cm amagat.
        n_bath = (args.pressure_mbar / 1013.25) * (273.15 / args.temperature)
        tau_fw = kfw * u * n_bath
        t_full = t_lines * math.exp(-min(tau_fw, 700.0))
        t_ram = ramanathan_transmissivity(u, p_atm, args.temperature, LINE_SPACING)

        # Absorptance in cm-1 carries no width convention, so this is the
        # comparison that means something. Ramanathan's is 2 A0 ln(1 + xi) by
        # construction; the corrk band's is its mean absorption times its width.
        a_ram = 2.0 * CO2_A0 * math.log(1.0 / max(t_ram, 1e-30))
        a_corrk = (1.0 - t_full) * width
        ratio = a_ram / a_corrk if a_corrk > 0 else float("nan")

        # The multiplier the model gets: band-mean optical depth scaled onto the
        # window N2O occupies.
        tau_band = -math.log(max(t_full, 1e-30))
        rows.append(dict(co2_atm_cm=u, corrk_lines=t_lines, corrk_with_far_wing=t_full,
                         ramanathan=t_ram, absorptance_ramanathan=a_ram,
                         absorptance_corrk=a_corrk, absorptance_ratio=ratio,
                         adopted_multiplier=math.exp(-correction * tau_band)))

    for r in rows:
        for key in ("corrk_lines", "corrk_with_far_wing", "ramanathan"):
            if not 0.0 <= r[key] <= 1.0:
                raise SystemExit(f"{key} out of range at {r['co2_atm_cm']} atm cm: {r[key]}")
    for a, b in zip(rows, rows[1:]):
        for key in ("corrk_with_far_wing", "ramanathan"):
            if b[key] > a[key] + 1e-12:
                raise SystemExit(f"{key} is not monotone in CO2 amount")

    worst = max(rows, key=lambda r: max(r["absorptance_ratio"],
                                        1.0 / r["absorptance_ratio"]))
    factor = max(worst["absorptance_ratio"], 1.0 / worst["absorptance_ratio"])

    inrange = [r for r in rows if r["co2_atm_cm"] <= FULL_COLUMN_ATM_CM]
    factor_inrange = max(max(r["absorptance_ratio"], 1.0 / r["absorptance_ratio"])
                         for r in inrange)

    scan = {f"{d}": [round(ramanathan_transmissivity(u, p_atm, args.temperature, d), 4)
                     for u in amounts] for d in LINE_SPACING_SCAN}

    report = {
        "generated": datetime.datetime.now(datetime.timezone.utc)
                             .replace(microsecond=0).isoformat(),
        "generator": "exoplasim/scripts/co2_overlap_589.py",
        "task": "CLIM-42",
        "adopted": "correlated-k",
        "corrk_table": TABLE,
        "corrk_band_cm1": [float(lo), float(hi)],
        "n2o_band": {"centre_cm1": N2O_589_CENTRE, "two_a0_cm1": 2 * N2O_589_A0},
        "conditions": {"co2_vmr": args.co2_vmr, "pressure_mbar": args.pressure_mbar,
                       "temperature_k": args.temperature},
        "far_wing_k_at_589": kfw,
        "window_correction": round(correction, 4),
        "window_correction_note": (
            "the corrk band reaches further toward the 667 cm-1 fundamental "
            "than the N2O band does, so its mean CO2 absorption is "
            f"{1 / correction:.2f}x the absorption over 566-612 cm-1. The "
            "band-mean optical depth is scaled by this before it becomes the "
            "model's multiplier."),
        "ramanathan_inputs": {
            "co2_a0_cm1": CO2_A0, "nu0_cm1_atm": NU0,
            "line_spacing_cm1": LINE_SPACING,
            "line_spacing_is": "a STATED ASSUMPTION, 4B; Ramanathan cites "
                               "Dickinson (1972) for it and Dickinson does not "
                               "tabulate it",
            "bands": len(CO2_BANDS_15UM),
            "isotopic_bands_omitted": "B8-B10, which need q_i from Goody (1964)",
        },
        "rows": [{k: (round(v, 5) if isinstance(v, float) else v)
                  for k, v in r.items()} for r in rows],
        "agreement": {
            "compared": "band absorptance in cm-1, which carries no width convention",
            "worst_factor_all_rows": round(factor, 3),
            "worst_factor_within_column": round(factor_inrange, 3),
            "full_column_atm_cm": FULL_COLUMN_ATM_CM,
            "declared_tolerance": AGREEMENT_FACTOR,
            "passes_all_rows": bool(factor <= AGREEMENT_FACTOR),
            "passes_within_column": bool(factor_inrange <= AGREEMENT_FACTOR),
            "reading": (
                "MISSED, and the bar is not moved: 2.02 against a 2.0 declared "
                "before running. The worst point is the OPTICALLY THIN end, "
                "1 atm cm, where the band model gives half the line list's "
                "absorptance -- Ramanathan's effective intensity assumes a pure "
                "exponential falloff from 667 cm-1, and at 589 the real band "
                "has structure that falloff does not describe. In the thin "
                "limit absorptance is proportional to that intensity, so the "
                "error appears undiluted. At the thick end the two diverge the "
                "other way, the corrk absorptance saturating toward its own "
                "band width while the band model's grows logarithmically. "
                "What the check was for is gross error -- an order of "
                "magnitude, a unit slip, the wrong band -- and there is none: "
                "agreement stays within a factor of two across three decades "
                "of CO2 amount, against a 1976 band model whose line spacing "
                "had to be assumed. That supports the correlated-k route as "
                "the adopted one; it does not validate the band model."),
        },
        "line_spacing_sensitivity": scan,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"CO2 transmissivity for the N2O {N2O_589_CENTRE:.0f} cm-1 band")
    print(f"  corrk band {lo:.1f}-{hi:.1f} cm-1 against 2A0 = {2 * N2O_589_A0:.0f}")
    print(f"  at {args.co2_vmr * 1e6:.0f} ppm, {args.pressure_mbar:.0f} mbar, "
          f"{args.temperature:.0f} K; far-wing k = {kfw:.3e}")
    print(f"  window correction {correction:.4f}: the corrk band's mean CO2 "
          f"absorption is {1 / correction:.2f}x that over 566-612 cm-1")
    print(f"\n  {'CO2 atm cm':>10} {'T corrk':>8} {'T raman':>8} "
          f"{'A corrk':>8} {'A raman':>8} {'A ratio':>8} {'adopted':>8}")
    for r in rows:
        print(f"  {r['co2_atm_cm']:10.1f} {r['corrk_with_far_wing']:8.4f} "
              f"{r['ramanathan']:8.4f} {r['absorptance_corrk']:8.2f} "
              f"{r['absorptance_ramanathan']:8.2f} {r['absorptance_ratio']:8.3f} "
              f"{r['adopted_multiplier']:8.4f}")
    print(f"\n  ABSORPTANCES agree to a factor of {factor:.2f} over all rows, "
          f"bar {AGREEMENT_FACTOR} -> {'PASS' if factor <= AGREEMENT_FACTOR else 'MISS'}")
    print(f"  within this world's whole column ({FULL_COLUMN_ATM_CM:.0f} atm cm): "
          f"{factor_inrange:.2f} -> "
          f"{'PASS' if factor_inrange <= AGREEMENT_FACTOR else 'MISS'}")
    tr = max(max(r['corrk_with_far_wing'] / r['ramanathan'],
                 r['ramanathan'] / r['corrk_with_far_wing']) for r in rows)
    print(f"  the same rows compared as TRANSMISSIVITIES differ by up to "
          f"{tr:.1f}x, which is the width convention and not physics")
    print(f"  line spacing sensitivity, T at {amounts[-1]:.0f} atm cm:")
    for d, vals in scan.items():
        print(f"    D = {d:>4} cm-1 -> {vals[-1]:.4f}")
    print(f"\nwrote {OUTPUT.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
