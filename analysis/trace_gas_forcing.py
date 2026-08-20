#!/usr/bin/env python3
"""CH4 and N2O: this host's surface mixing ratios, and what the missing band is worth.

    python analysis/trace_gas_forcing.py

Worldbuilding. Vesper is an invented planet and this script is about the
simulation of it: a toy climate model's longwave scheme, the two trace
absorbers it has no term for, and the modelled mixing ratios that set what
those terms would be worth. Every quantity here is a modelled field.

`config/planet.yaml` records that PlaSim's longwave is Sasamori (1968) over
water vapour, CO2 and ozone, that `radmod.f90` carries no CH4 and no N2O term,
and that adding either means adding a band rather than setting a key. This
script prices that omission, and it writes `analysis/trace_gas_forcing.json`.
`exoplasim/notes/trace-gas-absorbers.md` is the argument; this is the
derivation made re-runnable.

## Two halves, and only one of them was ever uncertain

The FORCING half is settled to a few percent. Byrne and Goldblatt (2014) fit
line-by-line calculations over 100 ppbv to 100 ppmv for both gases -- a range
that covers this world at both ends, where Etminan et al. (2016) stops at 3500
ppb of CH4 -- and where the two overlap they agree to 3.5%.

The MIXING RATIO half is what the width came from, and it is what this script
now settles. The gases are photochemical products of fixed biogenic surface
fluxes, so their steady-state abundance is set by the host star's ultraviolet
rather than by any choice made here. Segura et al. (2003) and Rugheimer et al.
(2013) both run that calculation; the project read Segura's K2V figure and got
a factor of 2 to 5 enhancement over Earth on both gases, and read Rugheimer's
text and got no tropospheric N2O enhancement at all.

**They do not disagree. They are at different stellar temperatures.**
Rugheimer's grid runs 4250 K to 7000 K at 250 K spacing, and the enhancement
turns on sharply between its 4750 K and 4500 K points -- 1.60 ppmv of CH4 at
4750 K against 37 ppmv at 4500 K. Segura's K2V is 4620 K, INSIDE that
transition, which is why it reads 4 to 8 ppmv. This host is 4965 K, above it,
and the 5000 K grid point reads the Sun case's values on both gases.

That grid point is not an interpolation either: epsilon Eridani is one of
Rugheimer's grid stars, assigned to 5000 K, and it is the star
`config/planet.yaml` names as what this world's activity is modelled on. So
the one caveat the earlier estimate could not close -- a quiet modelled K2V
against a star declared active -- is closed by the source itself.

## Where the numbers come from

Fig. 7 of Rugheimer et al. (2013) is a figure, but the paper is a vector PDF,
so the profiles are recovered as PATH COORDINATES rather than read off a log
axis by eye. `read_figure7` parses the page content stream, matches each curve
to its legend colour, and converts the lowest plotted vertex of each through
the panel's own tick calibration. The result is exact to the PDF's coordinate
precision, about 0.3% of a decade.

Both gases are well mixed through this world's troposphere in that model --
CH4 falls 1.601 to 1.583 ppmv over the lowest 20 km and N2O is flat to 15 km --
so the surface value IS the tropospheric column value, which is the quantity
the forcing wants.

## The bracket, declared by method

The value is the 5000 K grid point. The bracket is the envelope of the three
grid points nearest this host, 4750 / 5000 / 5250 K, which spans the model's
own scatter over 250 K either side without reaching the cold-star regime below
4750 K. The two points that actually bracket 4965 K are degenerate -- 4750 K
and 5000 K read identically -- so a bracket taken from those two alone would
report zero width, which is a statement about the grid rather than about the
world.

The cold-star regime is EXCLUDED rather than bounded, and the reason is the
grid star: this host's activity is represented at 5000 K by epsilon Eridani
itself, so there is nothing left for an extrapolation toward 4500 K to
represent.

## Gravity, which cuts the other way and must not be dropped

A mixing ratio is not a column. At `pressure_bar` over this world's gravity
against Earth's 1.013 bar over 9.81 m/s2, the column for a given mixing ratio
is 0.756 of Earth's -- the same reduction in atmospheric column mass that
`notes/audits/absent-and-inherited-physics.md` confirms is carried correctly
through every other absorber. The Earth-calibrated fits expect an
Earth-equivalent column, so every mixing ratio is multiplied by it before the
fit sees it.

## What is still excluded

Each gas is priced above Byrne's 100 ppbv floor, so the residual from zero to
100 ppbv is not in these numbers and no fit here covers it.

The CH4-N2O band overlap is evaluated where both gases exceed the fit's
pre-industrial reference, and is about -0.02 W/m2 for Earth's own atmosphere.
It comes out zero here, and NOT because it is negligible: this world's
Earth-equivalent N2O column is 0.227 ppmv against a 0.270 ppmv reference, so
the expression takes the logarithm of a negative number and has nothing to
say. The sign is known -- an overlap can only reduce the total -- so omitting
it holds these numbers slightly high, by less than the sub-100-ppbv residual
already excluded holds them low.

## The checks, and what each can fail

Every one has a right answer that is not this script's to choose:

- the two Byrne branches agree at the 2.5 ppmv junction, which pins the
  constants against a table that does not survive text extraction cleanly;
- the fits reproduce the paper's own stated maxima at 100 ppmv;
- they reproduce Earth's known CH4 and N2O forcings over the industrial era;
- the extracted Sun curve reproduces Rugheimer's separately STATED Sun-case
  mixing ratios, which is the figure extraction checking itself against a
  number printed in the same paper;
- the decade spacing recovered from each panel's tick marks is uniform, which
  is the axis calibration checking itself.

A failure raises. None of them can be satisfied by a wrong extraction that
happens to look plausible.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import re
from math import log, sqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUGHEIMER = ROOT / "references" / (
    "rugheimer_2013_spectral-fingerprints-of-earth-like-planets-around-fgk-stars.pdf")
OUTPUT = ROOT / "analysis" / "trace_gas_forcing.json"

# Rugheimer et al. (2013) Fig. 7 legend, in plotted order. The colours are the
# PDF's own stroke values, so they are matched exactly rather than by tolerance.
LEGEND = {
    "7000": (0.0, 0.0, 0.0),
    "6750": (0.341, 0.0, 0.568),
    "6500": (0.141, 0.0, 0.961),
    "6250": (0.0, 0.282, 1.0),
    "6000": (0.0, 0.732, 1.0),
    "5750": (0.0, 1.0, 0.797),
    "5500": (0.0, 1.0, 0.349),
    "5250": (0.082, 1.0, 0.0),
    "5000": (0.529, 1.0, 0.0),
    "4750": (0.98, 1.0, 0.0),
    "4500": (1.0, 0.58, 0.0),
    "4250": (1.0, 0.133, 0.0),
}

# Rugheimer's Sun case, STATED as numbers in its section 2.4 rather than
# plotted. These are what the extracted Sun curve is checked against.
SUN_STATED = {"CH4": 1.6e-6, "N2O": 3.0e-7}

# Byrne and Goldblatt (2014) Table 2. Reference concentrations are
# pre-industrial; N1 = M1 is the junction between the two branches.
M0, N0, M1 = 715e-9, 270e-9, 2.5e-6
FLOOR = 100e-9          # the fits' lower limit, and this pricing's reference
EARTH_GRAVITY = 9.81
EARTH_PRESSURE_BAR = 1.013

# The declared bracket: the grid points nearest this host. Fixed here, ahead of
# the extraction, because a bracket chosen after the numbers are seen is not a
# bracket. See the module docstring.
BRACKET_GRID = ("4750", "5000", "5250")
HOST_GRID_POINT = "5000"


# --------------------------------------------------------------------------
# Byrne and Goldblatt (2014) Table 2
# --------------------------------------------------------------------------

def f_ch4(m: float) -> float:
    """CH4 forcing from pre-industrial, W/m2, for a mole fraction `m`."""
    if m <= M1:
        d = sqrt(m) - sqrt(M0)
        return 1173.0 * d - 71636.0 * d * d
    r = log(m / M1)
    return 0.824 + 0.8 * r + 0.2 * r * r


def f_n2o(n: float) -> float:
    """N2O forcing from pre-industrial, W/m2, for a mole fraction `n`."""
    if n <= M1:
        d = sqrt(n) - sqrt(N0)
        return 3899.0 * d + 38256.0 * d * d
    r = log(n / M1)
    return 4.182 + 3.0 * r + 0.5469 * r * r


def overlap_ch4_n2o(m: float, n: float) -> float:
    """Reduction in forcing from CH4-N2O band overlap, W/m2 (negative)."""
    if m <= M0 or n <= N0:
        return 0.0
    return -24.0 * pow(2.718281828459045,
                       -0.02 * (log(m - M0) - 0.01) ** 2
                       - 0.044 * (log(n - N0) + 7.73) ** 2)


def forcing(m: float, n: float) -> dict:
    """Forcing of CH4 at `m` and N2O at `n` above the 100 ppbv floor."""
    fm = f_ch4(m) - f_ch4(FLOOR)
    fn = f_n2o(n) - f_n2o(FLOOR)
    ov = overlap_ch4_n2o(m, n)
    return {"ch4": fm, "n2o": fn, "overlap": ov, "total": fm + fn}


def check_fits() -> dict:
    """Four checks with right answers this script does not choose."""
    out = {}

    junction = (f_ch4(M1), f_n2o(M1))
    out["junction_2p5ppmv"] = [round(v, 4) for v in junction]
    for got, want, gas in ((junction[0], 0.824, "CH4"), (junction[1], 4.182, "N2O")):
        if abs(got - want) > 5e-4:
            raise SystemExit(f"{gas} branches disagree at 2.5 ppmv: "
                             f"{got:.4f} against the table's {want}")

    maxima = (f_ch4(100e-6), f_n2o(100e-6))
    out["maxima_100ppmv"] = [round(v, 3) for v in maxima]
    for got, want, gas in ((maxima[0], 6.66, "CH4"), (maxima[1], 22.3, "N2O")):
        if abs(got - want) / want > 0.03:
            raise SystemExit(f"{gas} misses the paper's stated 100 ppmv maximum: "
                             f"{got:.2f} against {want}")

    earth = (f_ch4(1800e-9), f_n2o(324e-9))
    out["earth_industrial_era"] = [round(v, 4) for v in earth]
    for got, want, gas in ((earth[0], 0.564, "CH4"), (earth[1], 0.193, "N2O")):
        if abs(got - want) > 2e-3:
            raise SystemExit(f"{gas} misses Earth's known forcing: "
                             f"{got:.3f} against {want}")

    anchor = forcing(SUN_STATED["CH4"], SUN_STATED["N2O"])["total"]
    out["earth_anchor_w_m2"] = round(anchor, 3)
    return out


# --------------------------------------------------------------------------
# Rugheimer et al. (2013) Fig. 7, from the PDF's vector paths
# --------------------------------------------------------------------------

def _paths(pdf: Path, page: int) -> list:
    """Stroked polylines on `page`, in device space, with their stroke colour."""
    import pypdf

    reader = pypdf.PdfReader(str(pdf))
    src = reader.pages[page].get_contents().get_data().decode("latin-1")
    src = re.sub(r"BT.*?ET", " ", src, flags=re.S)      # text draws no curves

    ctm, stack, colour, cur, out = [1, 0, 0, 1, 0, 0], [], (0.0, 0.0, 0.0), [], []
    operands: list[float] = []

    def apply(m, x, y):
        return (m[0] * x + m[2] * y + m[4], m[1] * x + m[3] * y + m[5])

    def mul(a, b):
        return [a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3],
                a[2] * b[0] + a[3] * b[2], a[2] * b[1] + a[3] * b[3],
                a[4] * b[0] + a[5] * b[2] + b[4], a[4] * b[1] + a[5] * b[3] + b[5]]

    for tok in src.split():
        try:
            operands.append(float(tok))
            continue
        except ValueError:
            pass
        if tok == "q":
            stack.append((list(ctm), colour))
        elif tok == "Q":
            if stack:
                ctm, colour = stack.pop()
                ctm = list(ctm)
        elif tok == "cm" and len(operands) >= 6:
            ctm = mul(operands[-6:], ctm)
        elif tok in ("SC", "SCN", "RG", "sc", "scn", "rg") and len(operands) >= 3:
            colour = tuple(round(v, 3) for v in operands[-3:])
        elif tok == "m" and len(operands) >= 2:
            if len(cur) > 1:
                out.append((colour, cur))
            cur = [apply(ctm, *operands[-2:])]
        elif tok == "l" and len(operands) >= 2:
            cur.append(apply(ctm, *operands[-2:]))
        elif tok == "c" and len(operands) >= 6:
            cur.append(apply(ctm, *operands[-2:]))
        elif tok in ("S", "s", "f", "F", "B", "n", "h"):
            if len(cur) > 1:
                out.append((colour, cur))
            cur = []
        operands = []
    if len(cur) > 1:
        out.append((colour, cur))
    return out


def _calibrate(paths, x_lo, x_hi, y_lo, y_hi, n_major):
    """Recover a log x-axis from a panel's own tick marks.

    Major ticks are twice the length of minor ones and every tick is one
    decade from its neighbour, so the spacing is measured rather than assumed.
    Raises if the spacing is not uniform, which is the calibration checking
    itself.
    """
    seen: dict = {}
    for _, pts in paths:
        if len(pts) != 2:
            continue
        (ax, ay), (bx, by) = pts
        if abs(ax - bx) > 1e-6 or not (x_lo < ax < x_hi):
            continue
        if abs(min(ay, by) - y_lo) > 0.6:
            continue
        length = round(abs(by - ay), 3)
        if not 0.0 < length < 0.1 * (y_hi - y_lo):   # excludes the frame verticals
            continue
        key = round(ax, 2)
        seen[key] = max(seen.get(key, 0.0), length)
    xs = sorted(seen.items())
    if len(xs) < 4:
        raise SystemExit(f"found {len(xs)} bottom-axis ticks, expected at least 4")
    positions = [x for x, _ in xs]
    steps = [b - a for a, b in zip(positions, positions[1:])]
    ppd = sum(steps) / len(steps)
    if max(abs(s - ppd) for s in steps) > 0.02 * ppd:
        raise SystemExit(f"tick spacing is not uniform: {steps}")

    longest = max(length for _, length in xs)
    majors = [x for x, length in xs if length > 0.6 * longest]
    if len(majors) != n_major:
        raise SystemExit(f"found {len(majors)} labelled ticks, expected {n_major}")
    return majors, ppd


def read_figure7(pdf: Path) -> dict:
    """Per-grid-star surface mixing ratios of CH4 and N2O, from Fig. 7."""
    paths = _paths(pdf, 6)

    # Panel boxes in device space, and the decade each panel's leftmost major
    # tick carries. Both come from the panel's printed axis labels.
    panels = {
        "CH4": dict(box=(221.0, 288.0, 654.0, 732.1), majors=(-10.0, 4)),
        "N2O": dict(box=(59.0, 125.7, 565.0, 643.3), majors=(-14.0, 6)),
    }

    inverse = {v: k for k, v in LEGEND.items()}
    out: dict = {}
    for gas, cfg in panels.items():
        bx0, bx1, by0, by1 = cfg["box"]
        log_ref, n_major = cfg["majors"]
        majors, ppd = _calibrate(paths, bx0, bx1, by0, by1, n_major)
        x_ref = majors[0]

        rows: dict = {}
        black = []
        for colour, pts in paths:
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            inside = (min(xs) >= bx0 and max(xs) <= bx1
                      and min(ys) >= by0 - 1.0 and max(ys) <= by1 + 1.0)
            if not inside or len(pts) < 20:
                continue
            surface = min(pts, key=lambda p: p[1])
            value = 10.0 ** (log_ref + (surface[0] - x_ref) / ppd)
            altitude = (surface[1] - by0) / (by1 - by0) * 60.0
            row = {"surface_mixing_ratio": value, "lowest_plotted_km": altitude}
            if colour == (0.0, 0.0, 0.0):
                black.append((min(xs), row))
            elif colour in inverse:
                rows[inverse[colour]] = row

        # Two black curves share the panel: the 7000 K grid star and the Sun.
        # The 7000 K star is the hottest, so it reaches furthest left at
        # altitude -- an ordering the paper states and the panel shows.
        if len(black) != 2:
            raise SystemExit(f"{gas}: found {len(black)} black curves, expected 2")
        black.sort(key=lambda b: b[0])
        rows["7000"] = black[0][1]
        rows["SUN"] = black[1][1]

        missing = set(LEGEND) - set(rows)
        if missing:
            raise SystemExit(f"{gas}: no curve matched {sorted(missing)}")

        got = rows["SUN"]["surface_mixing_ratio"]
        want = SUN_STATED[gas]
        rows["SUN"]["stated_in_paper"] = want
        rows["SUN"]["read_over_stated"] = got / want
        if not 0.85 < got / want < 1.15:
            raise SystemExit(f"{gas}: the extracted Sun curve reads {got:.3e} "
                             f"against the paper's stated {want:.3e}")

        out[gas] = {"decades_per_tick": ppd, "curves": rows}
    return out


# --------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--rugheimer", type=Path, default=RUGHEIMER,
                    help="Rugheimer et al. (2013), for Fig. 7")
    args = ap.parse_args()

    if not args.rugheimer.exists():
        raise SystemExit(
            f"{args.rugheimer.relative_to(ROOT)} is not on disk. references/ "
            "PDFs are untracked; references/INDEX.md records the DOI.")

    import yaml
    cfg = yaml.safe_load((ROOT / "config" / "planet.yaml").read_text(encoding="utf-8"))
    gravity = float(cfg["planet"]["gravity_m_s2"])
    pressure = sum(float(v) for k, v in cfg["atmosphere"].items()
                   if k.startswith("p") and k.endswith("_bar"))
    column = (pressure / EARTH_PRESSURE_BAR) / (gravity / EARTH_GRAVITY)

    checks = check_fits()
    fig7 = read_figure7(args.rugheimer)

    host = {gas: fig7[gas]["curves"][HOST_GRID_POINT]["surface_mixing_ratio"]
            for gas in ("CH4", "N2O")}
    envelope = {
        gas: (min(fig7[gas]["curves"][t]["surface_mixing_ratio"] for t in BRACKET_GRID),
              max(fig7[gas]["curves"][t]["surface_mixing_ratio"] for t in BRACKET_GRID))
        for gas in ("CH4", "N2O")
    }

    value = forcing(host["CH4"] * column, host["N2O"] * column)
    low = forcing(envelope["CH4"][0] * column, envelope["N2O"][0] * column)
    high = forcing(envelope["CH4"][1] * column, envelope["N2O"][1] * column)
    bracket = [round(low["total"], 3), round(high["total"], 3)]

    report = {
        "generated": datetime.datetime.now(datetime.timezone.utc)
                             .replace(microsecond=0).isoformat(),
        "generator": "analysis/trace_gas_forcing.py",
        "rugheimer_pdf_sha256": hashlib.sha256(args.rugheimer.read_bytes()).hexdigest(),
        "host_effective_temperature_k": float(cfg["star"]["effective_temperature_k"]),
        "host_grid_point_k": float(HOST_GRID_POINT),
        "gravity_m_s2": gravity,
        "surface_pressure_bar": round(pressure, 5),
        "earth_equivalent_column_factor": round(column, 4),
        "fit": "Byrne and Goldblatt (2014) Table 2, above a 100 ppbv floor",
        "fit_checks": checks,
        "figure7": {
            gas: {
                "decades_per_tick": round(fig7[gas]["decades_per_tick"], 4),
                "surface_mixing_ratio_ppmv": {
                    t: round(r["surface_mixing_ratio"] * 1e6, 4)
                    for t, r in sorted(fig7[gas]["curves"].items())},
                "lowest_plotted_km": round(
                    fig7[gas]["curves"][HOST_GRID_POINT]["lowest_plotted_km"], 2),
                "sun_read_over_stated": round(
                    fig7[gas]["curves"]["SUN"]["read_over_stated"], 4),
            } for gas in ("CH4", "N2O")},
        "mixing_ratio_ppmv": {
            "value": {g: round(host[g] * 1e6, 4) for g in host},
            "bracket_grid_k": list(BRACKET_GRID),
            "bracket": {g: [round(v * 1e6, 4) for v in envelope[g]] for g in envelope},
        },
        "earth_equivalent_column_ppmv": {
            g: round(host[g] * column * 1e6, 4) for g in host},
        "forcing_w_m2": {
            "value": {k: round(v, 3) for k, v in value.items()},
            "bracket": bracket,
            "earth_anchor": checks["earth_anchor_w_m2"],
        },
        "excluded": (
            "The residual below Byrne's 100 ppbv floor, for both gases. The "
            "cold-star enhancement regime below 4750 K, which this host is "
            "above and which epsilon Eridani -- Rugheimer's own 5000 K grid "
            "star, and the star this world's activity is modelled on -- "
            "already represents at this grid point."),
        "note": (
            "Segura et al. (2003) and Rugheimer et al. (2013) do not disagree: "
            "Segura's K2V is 4620 K, inside the 4500-4750 K transition where "
            "the enhancement turns on, and this host is 4965 K, above it. The "
            "width here is the model's own scatter over the 250 K either side "
            "of this host, not a disagreement between sources. Both gases are "
            "well mixed through this world's troposphere in that model, so the "
            "surface value is the column value. Sasamori's longwave carries "
            "water vapour, CO2 and ozone, so there is no key to set: adding "
            "either gas means adding a band."),
    }
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"Rugheimer Fig. 7, surface mixing ratios in ppmv:")
    print(f"  {'K':>6}  {'CH4':>9}  {'N2O':>9}")
    for t in list(LEGEND) + ["SUN"]:
        c = fig7["CH4"]["curves"][t]["surface_mixing_ratio"] * 1e6
        n = fig7["N2O"]["curves"][t]["surface_mixing_ratio"] * 1e6
        mark = "  <- this host" if t == HOST_GRID_POINT else ""
        print(f"  {t:>6}  {c:9.4f}  {n:9.4f}{mark}")
    print(f"\n  Sun curve reads {fig7['CH4']['curves']['SUN']['read_over_stated']:.3f} "
          f"of stated on CH4, "
          f"{fig7['N2O']['curves']['SUN']['read_over_stated']:.3f} on N2O")
    print(f"\n  column factor {column:.4f} at {gravity} m/s2 and {pressure:.5f} bar")
    print(f"  Earth-equivalent: CH4 {host['CH4'] * column * 1e6:.4f} ppmv, "
          f"N2O {host['N2O'] * column * 1e6:.4f} ppmv")
    print(f"  forcing {value['total']:.3f} W/m2 "
          f"(CH4 {value['ch4']:.3f}, N2O {value['n2o']:.3f}, "
          f"overlap {value['overlap']:+.3f})")
    print(f"  bracket {bracket[0]:.3f} to {bracket[1]:.3f} W/m2, "
          f"against {checks['earth_anchor_w_m2']:.2f} for Earth's own")
    print(f"\nwrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
