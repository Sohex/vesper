"""Check the Howard-1956 shortwave absorptances against modern correlated-k data.

`shortwave_band_weights.py` derives `h2osww` and `co2sww`, and the CO2 patch's
closed form, from Howard, Burch and Williams (1956) laboratory band absorptions
run through Yamamoto's (1962) construction. That is 1950s data reached through
two closed-form fits, and nothing in this repository had checked it against a
line-by-line treatment. This script does that, against the LMD Generic PCM's
correlated-k tables (HITRAN2020 through SpeCT and exo-k, Chaverot et al. 2025).

WHAT IS COMPARED, AND WHY IT IS COMPARABLE
------------------------------------------
Both sides compute the same quantity: the fraction of a star's TOTAL incident
flux absorbed by one gas in a HOMOGENEOUS path at 760 mm Hg holding a stated
absorber amount. That is Yamamoto's definition, it is what Lacis and Hansen
Eq. 21 fits, and it is what `shortwave_band_weights.py` reconstructs. Only the
absorption data differs:

    project   Howard's total band absorption per band, divided by the band
              width for his Eq. 11 band-average, over 8 CO2 and 9 H2O bands
    here      1 - sum_g w_g exp(-k_g u) per correlated-k band, over 78 bands
              from 10 to 30000 cm-1

The SPECTRA ARE THE SAME OBJECTS. This module imports `blend` and `Spectrum`
from `shortwave_band_weights`, so the BT-Settl 4965 K and 5772 K blends are
byte-identical to the ones the derivation used and a model systematic in the
spectrum cannot show up as a disagreement about the absorption.

WHAT CANNOT BE MADE COMPARABLE, AND IS THEREFORE REPORTED RATHER THAN HIDDEN
----------------------------------------------------------------------------
- The tables are PREMIXED N2+CO2+H2O, so no table holds H2O without CO2. H2O is
  isolated by dividing the mixture transmission by the dry-slice transmission,
  which is the random-overlap assumption -- the same one the project makes when
  it charges CO2 with the fraction of its interval water vapour leaves.
  `--checks` proves the division works: it must return the same H2O absorptance
  from the 376 ppm and the 1000 ppm table, and it does to 0.2%.
- The tables carry LINE CENTRES ONLY, +-25 cm-1, with the plinth removed on the
  MT_CKD convention. The H2O self and foreign continuum is absent and the bundle
  has no H2O far-wing file to add it back. It absorbs in the windows BETWEEN the
  bands, which is where this star's flux boost is largest, so its omission
  pushes the corrk H2O weight DOWN. The disagreement reported here is therefore
  an upper bound on the shortfall, not a central estimate.
- Neither side is Earth's real column. Both evaluate a homogeneous 760 mm Hg
  path, because that is what Howard's constants describe and what Eq. 21 is
  stated at; the model reaches it by scaling the amount, not the pressure.

THE DATA
--------
Read-only from `~/git/generic_pcm`, which is CeCILL licensed: it is a DESIGN
REFERENCE and no code crosses. `notes/external-data.md` records the route.
The reading of the file layout is checked rather than assumed: `--checks`
verifies that k scales exactly with the CO2 mixing ratio between the 376 ppm and
1000 ppm tables in the bands where CO2 dominates, which it must, because CO2 is
a trace N2-broadened gas in both.

This writes nothing. The finding is `exoplasim/notes/corrk-cross-check.md`.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from _paths import ANALYSIS
from shortwave_band_weights import (CO2_BANDS, CO2_FIT_RANGE, Spectrum, blend,
                                    co2_closed_form, SOLAR_TEFF)

# The Generic PCM bundle, outside this repo and read-only.
CORRK = Path.home() / "git" / "generic_pcm" / "LMDZ.GENERIC" / "datagcm" / "corrk_data"
TABLE_376 = "N2-0.000376CO2-H2Ovar_2026"
TABLE_1000 = "N2-0.001CO2-H2Ovar_2026"
BANDS = "40x38"

LOSCHMIDT = 2.6867811e19      # molecules cm-3 at 273.15 K and 101325 Pa
H2O_PER_PRCM = 3.3428e22      # molecules cm-2 in one precipitable cm
P_STANDARD_MBAR = 1013.25     # the 760 mm Hg the Howard fits are stated at
DRY = 1e-8                    # the tables' lowest H2O mixing ratio

# The paths `shortwave_band_weights.py` evaluates its absorptances at, so the two
# calculations are quoted at the same amount. Recompute them there rather than
# trusting these if the climatology moves; they are arguments, not state.
DEFAULT_CO2_PLANET = 225.5626   # atmos-cm, column * pressure reduction * zbetta
DEFAULT_CO2_EARTH = 298.5466    # the same at Earth's gravity
DEFAULT_WATER_CM = 2.7891       # precipitable cm, magnified by zbetta
DEFAULT_WATER_VMR = 1e-2        # broadening partner fraction; --checks scans it
DEFAULT_TEMPERATURE = 290.0     # Howard's laboratory was room temperature
EARTH_MEAN_INSOLATION = 1361.0 / 4.0

# The derivation's own CO2 band intervals, taken from it rather than restated,
# so --bands cannot drift from the thing it is comparing against.
HOWARD_CO2_INTERVALS = {name: (b["lo"], b["hi"]) for name, b in CO2_BANDS.items()}


class CorrK:
    """One premixed correlated-k table, IR and VI bands joined without overlap."""

    def __init__(self, name: str, root: Path = CORRK, bands: str = BANDS):
        d = root / name
        if not d.is_dir():
            raise SystemExit(
                f"{d} not found. This check needs the LMD Generic PCM datagcm bundle; "
                "notes/external-data.md has the route."
            )
        self.name = name
        self.p = np.loadtxt(d / "p.dat", skiprows=1)          # log10 mbar
        self.T = np.loadtxt(d / "T.dat", skiprows=1)          # K
        token = (d / "Q.dat").read_text().split()
        ngas = int(token[0])
        self.gases = token[1:1 + ngas]
        nq = int(token[1 + ngas])
        self.Q = np.array([float(x) for x in token[2 + ngas:2 + ngas + nq]])
        gw = np.loadtxt(d / "g.dat", skiprows=1)
        self.g = gw[:-1]                                      # the last weight is 0

        edges, blocks = [], []
        for tag in ("IR", "VI"):
            token = (d / bands / f"narrowbands_{tag}.in").read_text().split()
            nb = int(token[0])
            e = np.array([float(x) for x in token[1:]]).reshape(nb, 2)
            raw = np.fromstring((d / bands / f"corrk_gcm_{tag}.dat").read_text(), sep=" ")
            shape = (len(self.T), len(self.p), len(self.Q), nb, len(gw))
            if raw.size != int(np.prod(shape)):
                raise SystemExit(f"{tag} table is {raw.size} values, expected {np.prod(shape)}")
            # gasv8(L_NTREF, L_NPREF, L_REFVAR, L_NSPECT, L_NGAUSS), Fortran order.
            k = raw.reshape(shape, order="F")[..., :len(self.g)]
            # IR runs 10-3000 cm-1 and VI 2000-30000; keep IR only below the join.
            keep = e[:, 1] <= 2000.0 + 1e-6 if tag == "IR" else np.ones(nb, bool)
            edges.append(e[keep])
            blocks.append(k[:, :, :, keep])
        self.edges = np.vstack(edges)
        self.k = np.concatenate(blocks, axis=3)
        self.logk = np.log10(np.maximum(self.k, 1e-200))

    def kvec(self, p_mbar: float, t_k: float, q: float) -> np.ndarray:
        """k per air molecule, (band, g), trilinear in log10 p, T and log10 q."""
        def bracket(grid, value):
            i = int(np.clip(np.searchsorted(grid, value) - 1, 0, len(grid) - 2))
            return i, (value - grid[i]) / (grid[i + 1] - grid[i])

        it, wt = bracket(self.T, t_k)
        ip, wp = bracket(self.p, math.log10(p_mbar))
        iq, wq = bracket(np.log10(self.Q), math.log10(q))
        out = np.zeros(self.k.shape[3:])
        for a, wa in ((it, 1 - wt), (it + 1, wt)):
            for b, wb in ((ip, 1 - wp), (ip + 1, wp)):
                for c, wc in ((iq, 1 - wq), (iq + 1, wq)):
                    out += wa * wb * wc * self.logk[a, b, c]
        return 10.0 ** out

    def transmission(self, p_mbar: float, t_k: float, q: float, u_air: float) -> np.ndarray:
        """sum_g w_g exp(-k_g u_air) per band. u_air in air molecules cm-2."""
        tau = np.clip(self.kvec(p_mbar, t_k, q) * u_air, 0.0, 700.0)
        return (self.g[None, :] * np.exp(-tau)).sum(axis=1)


def flux_fractions(table: CorrK, spectrum: Spectrum) -> np.ndarray:
    return np.array([spectrum.fraction_in_band(lo, hi) for lo, hi in table.edges])


def broadband(fraction: np.ndarray, band_absorption: np.ndarray) -> float:
    return float((fraction * band_absorption).sum())


def co2_vmr_of(name: str) -> float:
    """The CO2 mixing ratio the table is premixed at, from its own directory name."""
    return float(name.split("-")[1].replace("CO2", ""))


def air_column_for_co2(atmos_cm: float, vmr: float) -> float:
    return atmos_cm * LOSCHMIDT / vmr


def air_column_for_water(precipitable_cm: float, vmr: float) -> float:
    return precipitable_cm * H2O_PER_PRCM / vmr


def spectra(refresh: bool = False) -> tuple[Spectrum, Spectrum]:
    """The derivation's own spectra, so only the absorption data differs."""
    out = []
    for teff, lo, hi, label in ((SOLAR_TEFF, 5700, 5800, "sun 5772 K"),
                                (4965.0, 4900, 5000, "k25v 4965 K")):
        wave, flux, meta = blend(teff, lo, hi, refresh)
        out.append(Spectrum(wave, flux, label, meta))
    return out[0], out[1]


def run_checks(t376: CorrK, t1000: CorrK, sun: Spectrum, water_cm: float, t_k: float) -> None:
    """The three things that could have failed, before any number is quoted."""
    ratio_expected = co2_vmr_of(TABLE_1000) / co2_vmr_of(TABLE_376)
    it = int(np.argmin(abs(t376.T - t_k)))
    ip = int(np.argmin(abs(t376.p - math.log10(P_STANDARD_MBAR))))
    a, b = t376.k[it, ip, 0], t1000.k[it, ip, 0]
    # CO2-dominated bands only: where the ratio is meaningful the k must scale
    # EXACTLY with the mixing ratio, because CO2 is trace and N2-broadened in both.
    strong = a.max(axis=1) > 1e-27
    r = (b / np.maximum(a, 1e-300))[strong]
    co2ish = r > 2.0
    print(f"  k(1000 ppm)/k(376 ppm) where CO2 dominates: expected {ratio_expected:.4f}, "
          f"got {np.median(r[co2ish]):.4f} "
          f"[{np.percentile(r[co2ish], 2):.4f}, {np.percentile(r[co2ish], 98):.4f}] "
          f"over {co2ish.sum()} points")

    f = flux_fractions(t376, sun)
    print(f"  flux fraction inside the 10-30000 cm-1 table span, solar: {f.sum():.5f} "
          "(the rest is below 0.33 um, where neither gas absorbs)")

    print("  T_mix/T_dry must return the same H2O absorptance from both tables:")
    for q in (1e-3, 1e-2, 1e-1):
        u = air_column_for_water(water_cm, q)
        vals = []
        for tab in (t376, t1000):
            aa = 1.0 - tab.transmission(P_STANDARD_MBAR, t_k, q, u) / tab.transmission(
                P_STANDARD_MBAR, t_k, DRY, u)
            vals.append(broadband(flux_fractions(tab, sun), aa))
        print(f"    q={q:.0e}  {vals[0]:.6f} vs {vals[1]:.6f}  ({100*(vals[1]/vals[0]-1):+.2f}%)")


def fit_against_line_list(table: CorrK, sun: Spectrum, clear: np.ndarray,
                          vmr: float, t_k: float) -> dict:
    """Refit radmod's two-log CO2 closed form to the correlated-k absorptance.

    The form, the path range and the fit protocol are `shortwave_band_weights`'s
    own, imported rather than restated, so the only thing that changes between
    the two fits is where the absorptance comes from: Howard's 1956 band set
    there, HITRAN2020 through the correlated-k tables here.

    WHY THE WHOLE CURVE AND NOT TWO CONSTANTS. PHYS-10 was opened to correct the
    per-band H2O overlap at 2.7 and 2.0 um, and correcting only those makes the
    total WORSE. Summed over Howard's eight intervals the correlated-k
    absorptance is 0.004494 against the derivation's 0.005536, 18.8% low, while
    the correlated-k answer over ALL bands is 0.005098, only 7.9% low: about an
    eighth of the CO2 shortwave absorption falls outside every interval Howard
    measured, and the band-mean overlap error was partly standing in for it. So
    the honest fix is to take the level from the line list across the range,
    which needs no band bookkeeping at all.
    """
    from scipy.optimize import curve_fit

    lo, hi = CO2_FIT_RANGE
    u = np.logspace(math.log10(lo), math.log10(hi), 160)
    fs = flux_fractions(table, sun)
    y = np.array([
        broadband(fs, (1.0 - table.transmission(P_STANDARD_MBAR, t_k, DRY,
                                                air_column_for_co2(float(ui), vmr))) * clear)
        for ui in u
    ])
    popt, _ = curve_fit(
        co2_closed_form, u, y, p0=[1.0e-3, 1.0, 1.0e-3, 1.0e-2], sigma=y,
        bounds=([0.0, 0.0, 0.0, 0.0], [np.inf] * 4), maxfev=400000,
    )
    coefficients = [float(f"{v:.5g}") for v in popt]
    residual = co2_closed_form(u, *coefficients) / y - 1.0
    # The two-log form was chosen against Howard's SHAPE, so it is worth knowing
    # whether it fits the line list worse everywhere or only at the ends of a
    # range the model never visits. A T42 column sits near 200-300 atmos-cm.
    near = (u >= 100.0) & (u <= 1000.0)
    return {"a1": coefficients[0], "b1": coefficients[1],
            "a2": coefficients[2], "b2": coefficients[3],
            "max_relative_error": float(np.abs(residual).max()),
            "rms_relative_error": float(np.sqrt((residual ** 2).mean())),
            "max_relative_error_100_to_1000": float(np.abs(residual[near]).max()),
            "worst_at_atmos_cm": float(u[np.abs(residual).argmax()]),
            "u": u, "y": y}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--temperature", type=float, default=DEFAULT_TEMPERATURE)
    ap.add_argument("--water-cm", type=float, default=DEFAULT_WATER_CM)
    ap.add_argument("--water-vmr", type=float, default=DEFAULT_WATER_VMR)
    ap.add_argument("--co2-planet", type=float, default=DEFAULT_CO2_PLANET)
    ap.add_argument("--co2-earth", type=float, default=DEFAULT_CO2_EARTH)
    ap.add_argument("--checks", action="store_true", help="run the falsifiable checks and stop")
    ap.add_argument("--bands", action="store_true",
                    help="per-band CO2 against the derivation's own Howard intervals")
    ap.add_argument("--scan", action="store_true", help="report the sensitivity to every choice")
    ap.add_argument("--fit", action="store_true",
                    help="refit radmod's CO2 closed form to the line list (PHYS-10)")
    ap.add_argument("--refresh", action="store_true", help="re-fetch the BT-Settl grid points")
    args = ap.parse_args()

    sun, star = spectra(args.refresh)
    t376 = CorrK(TABLE_376)
    vmr = co2_vmr_of(TABLE_376)
    t_k, q = args.temperature, args.water_vmr

    if args.checks:
        print("CHECKS")
        run_checks(t376, CorrK(TABLE_1000), sun, args.water_cm, t_k)
        return

    fs, fr = flux_fractions(t376, sun), flux_fractions(t376, star)
    u_water = air_column_for_water(args.water_cm, q)
    t_dry_w = t376.transmission(P_STANDARD_MBAR, t_k, DRY, u_water)
    t_wet = t376.transmission(P_STANDARD_MBAR, t_k, q, u_water)
    clear = t_wet / t_dry_w                       # H2O-only transmission per band
    a_h2o = 1.0 - clear

    print(f"correlated-k: {TABLE_376}, {BANDS} bands, {t_k:.0f} K, {P_STANDARD_MBAR:.2f} mbar")
    print(f"water {args.water_cm:.4f} precipitable cm at vmr {q:.0e}; "
          f"CO2 paths {args.co2_planet:.2f} and {args.co2_earth:.2f} atmos-cm\n")

    a_sun, a_star = broadband(fs, a_h2o), broadband(fr, a_h2o)
    print(f"H2O   absorptance solar {a_sun:.6f}  star {a_star:.6f}  h2osww {a_star/a_sun:.4f}")

    for label, u_cm in (("planet", args.co2_planet), ("Earth ", args.co2_earth)):
        a_co2 = 1.0 - t376.transmission(P_STANDARD_MBAR, t_k, DRY, air_column_for_co2(u_cm, vmr))
        bare = broadband(fs, a_co2)
        net_s, net_r = broadband(fs, a_co2 * clear), broadband(fr, a_co2 * clear)
        print(f"CO2 {label} absorptance solar {net_s:.6f}  star {net_r:.6f}  "
              f"co2sww {net_r/net_s:.4f}   before the H2O overlap {bare:.6f} "
              f"(suppression {net_s/bare:.3f})   {net_s * EARTH_MEAN_INSOLATION:.2f} W/m2")

    if args.fit:
        report = json.loads((ANALYSIS / "shortwave_band_weights.json").read_text())
        old = report["co2"]["closed_form_fit"]
        fit = fit_against_line_list(t376, sun, clear, vmr, t_k)
        print("\nradmod.f90 CO2 shortwave closed form, A(u) = a1 ln(1+b1 u) + a2 ln(1+b2 u)")
        print(f"  fitted over {CO2_FIT_RANGE[0]:g} to {CO2_FIT_RANGE[1]:g} atmos-cm "
              f"at {args.water_cm:.4f} precipitable cm of water\n")
        print(f"  {'':>4} {'Howard 1956':>14} {'line list':>14} {'change':>10}")
        for k in ("a1", "b1", "a2", "b2"):
            print(f"  {k:>4} {old[k]:14.5g} {fit[k]:14.5g} "
                  f"{100*(fit[k]/old[k]-1):+9.1f}%")
        print(f"\n  fit quality: max {100*fit['max_relative_error']:.2f}% at "
              f"u = {fit['worst_at_atmos_cm']:.3g} atmos-cm, "
              f"rms {100*fit['rms_relative_error']:.2f}% relative "
              f"(Howard fit: max {100*old['max_relative_error']:.2f}%, "
              f"rms {100*old['rms_relative_error']:.2f}%)")
        print(f"  over 100-1000 atmos-cm, where a T42 column sits: max "
              f"{100*fit['max_relative_error_100_to_1000']:.2f}%")
        for u_at in (args.co2_planet, args.co2_earth):
            a_old = co2_closed_form(u_at, old["a1"], old["b1"], old["a2"], old["b2"])
            a_new = co2_closed_form(u_at, fit["a1"], fit["b1"], fit["a2"], fit["b2"])
            print(f"  at u = {u_at:8.2f} atmos-cm: {a_old:.6f} -> {a_new:.6f} "
                  f"({100*(a_new/a_old-1):+.1f}%), "
                  f"{(a_new-a_old)*EARTH_MEAN_INSOLATION:+.2f} W/m2 of insolation")
        print("\nThis prints. Put the four numbers in radmod.f90 and rebuild; TASKS.md PHYS-10.")
        return

    if args.bands:
        # Where the two disagree, band by band, on the intervals the derivation
        # itself uses. The corrk bands are apportioned onto Howard's intervals by
        # the fraction of each corrk band that falls inside, so nothing is
        # counted twice and the flux weighting stays the spectrum's own.
        report = json.loads((ANALYSIS / "shortwave_band_weights.json").read_text())
        per_band = report["co2"]["per_band"]
        a_co2 = 1.0 - t376.transmission(P_STANDARD_MBAR, t_k, DRY,
                                        air_column_for_co2(args.co2_planet, vmr))
        width = t376.edges[:, 1] - t376.edges[:, 0]
        print(f"\n{'band':>6} {'corrk clear':>12} {'derivation':>11} "
              f"{'corrk contrib':>14} {'derivation':>11} {'delta':>10}")
        total_c = total_d = 0.0
        for name, row in per_band.items():
            lo, hi = HOWARD_CO2_INTERVALS[name]
            share = np.clip(np.minimum(t376.edges[:, 1], hi)
                            - np.maximum(t376.edges[:, 0], lo), 0.0, None) / width
            got = broadband(fs * share, a_co2 * clear)
            bare = broadband(fs * share, a_co2)
            total_c += got
            total_d += row["contribution"]
            print(f"{name:>6} {got/max(bare, 1e-30):12.3f} {row['h2o_transmission_in_band']:11.3f} "
                  f"{got:14.6f} {row['contribution']:11.6f} {got - row['contribution']:+10.6f}")
        print(f"{'sum':>6} {'':12} {'':11} {total_c:14.6f} {total_d:11.6f} "
              f"{total_c - total_d:+10.6f}")
        print(f"corrk over ALL bands rather than only Howard's eight: "
              f"{broadband(fs, a_co2 * clear):.6f}")
        return

    if not args.scan:
        return

    print("\nEvery choice this comparison had to make, and what it is worth:")
    print("  H2O weight vs water path, at the stated vmr and temperature")
    for w in (0.01, 0.1, 1.0, args.water_cm, 5.0, 10.0):
        u = air_column_for_water(w, q)
        aa = 1.0 - t376.transmission(P_STANDARD_MBAR, t_k, q, u) / t376.transmission(
            P_STANDARD_MBAR, t_k, DRY, u)
        lh = 2.9 * w / ((1 + 141.5 * w) ** 0.635 + 5.925 * w)
        s = broadband(fs, aa)
        print(f"    w={w:8.4f} cm  A_sun={s:.6f}  weight={broadband(fr, aa)/s:.4f}  "
              f"corrk/LacisHansen-Eq21={s/lh:.3f}")
    print("  H2O weight vs the broadening partner fraction, at the stated path")
    for qq in (1e-3, 3e-3, 1e-2, 3e-2, 1e-1):
        u = air_column_for_water(args.water_cm, qq)
        aa = 1.0 - t376.transmission(P_STANDARD_MBAR, t_k, qq, u) / t376.transmission(
            P_STANDARD_MBAR, t_k, DRY, u)
        print(f"    vmr={qq:.0e}  weight={broadband(fr, aa)/broadband(fs, aa):.4f}")
    print("  H2O and CO2 weights vs temperature, at the stated paths")
    for tt in (230.0, 260.0, 290.0, 320.0):
        u = air_column_for_water(args.water_cm, q)
        aa = 1.0 - t376.transmission(P_STANDARD_MBAR, tt, q, u) / t376.transmission(
            P_STANDARD_MBAR, tt, DRY, u)
        ac = 1.0 - t376.transmission(P_STANDARD_MBAR, tt, DRY,
                                     air_column_for_co2(args.co2_planet, vmr))
        print(f"    T={tt:5.0f} K  h2osww={broadband(fr, aa)/broadband(fs, aa):.4f}  "
              f"co2sww={broadband(fr, ac)/broadband(fs, ac):.4f}  "
              f"A_CO2_solar(no overlap)={broadband(fs, ac):.6f}")
    print("  CO2 absorptance vs column, solar-weighted, no overlap")
    for u_cm in (1.0, 10.0, 100.0, 1000.0, 1.0e4):
        ac = 1.0 - t376.transmission(P_STANDARD_MBAR, t_k, DRY, air_column_for_co2(u_cm, vmr))
        print(f"    u={u_cm:9.1f} atmos-cm  A_sun={broadband(fs, ac):.6f}  "
              f"weight={broadband(fr, ac)/broadband(fs, ac):.4f}")


if __name__ == "__main__":
    main()
