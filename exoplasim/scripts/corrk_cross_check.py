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
    here      1 - sum_g w_g exp(-k_g u) per correlated-k band, over the
              bundle's IR and VI band sets from 10 to 30000 cm-1. They are
              joined without overlap and do not meet: `run_checks` measures the
              gap between them and the flux that falls in it.

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
  The checks prove the division works, and they RUN BEFORE ANY NUMBER IS
  QUOTED and raise: the two tables must return the same H2O absorptance, and
  the bound is the CO2 absorptance the division exists to remove, because a
  residual of removing something cannot exceed the thing removed. `--checks`
  runs them and stops there; every other invocation runs them first.
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
REFERENCE and no code crosses. `docs/src/reference/external-data.md` records the route.
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
                                    co2_closed_form, radmod_co2_fit, SOLAR_TEFF)

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
    """One premixed correlated-k table, its IR and VI band sets joined here.

    The join is this script's, not the bundle's: see `__init__`. `join` is the
    wavenumber the two sets meet at and `truncated` is the IR band that was cut
    back to reach it, which a consumer has to know before pricing anything.
    """

    def __init__(self, name: str, root: Path = CORRK, bands: str = BANDS):
        d = root / name
        if not d.is_dir():
            raise SystemExit(
                f"{d} not found. This check needs the LMD Generic PCM datagcm bundle; "
                "docs/src/reference/external-data.md has the route."
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

        half = {}
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
            # EACH SET IS CONTIGUOUS IN ITSELF, to the bit, and the join below is
            # then the only place this script can open a hole. Checking it here
            # is what makes that statement a fact about the files rather than an
            # assumption about them.
            step = float(np.abs(e[1:, 0] - e[:-1, 1]).max()) if nb > 1 else 0.0
            if step != 0.0:
                raise SystemExit(
                    f"the {tag} set of {name}/{bands} is not contiguous: adjacent "
                    f"edges differ by up to {step:g} cm-1. The join below assumes "
                    "each half partitions its own span")
            half[tag] = (e, k)

        # THE JOIN IS THIS SCRIPT'S AND NOT THE BUNDLE'S. Upstream the two sets
        # OVERLAP -- IR runs 10-3000 cm-1 and VI 2000-30000 -- and the Generic
        # PCM never concatenates them: `rad_correlatedk_read_opacity_tables.F90`
        # sets `IR_VI_wnlimit = 3000.` and hands `WNOI` and `WNOV` to two
        # independent solvers. One flux-weighted partition is what the checks
        # below need, so the two are joined here, and where that join is put is
        # a decision this file owns.
        #
        # IT IS THE VI SET'S OWN FIRST EDGE, so the halves meet to the bit, and
        # the IR band that straddles it is TRUNCATED rather than dropped. Cutting
        # at a round 2000 and dropping every IR band reaching past it left
        # 1974.95 to 2000 in neither half: 1.7e-4 of the solar flux inside the
        # table span and inside no band, which every `broadband` below then
        # priced as transparent. The truncated band keeps its own
        # k-distribution, which is the only opacity the bundle carries over that
        # sliver, and applying a band's distribution to a tenth of its width is
        # an approximation where treating it as clear is a hole. world-olt.
        #
        # NOTHING IS WIDENED. An edge moved outward would make the arithmetic in
        # `run_checks` telescope again while silencing the only instrument that
        # can see a hole; this narrows one band and leaves that check sharp.
        e_ir, k_ir = half["IR"]
        e_vi, k_vi = half["VI"]
        self.join = float(e_vi[0, 0])
        keep = e_ir[:, 0] < self.join
        e_ir, k_ir = e_ir[keep].copy(), k_ir[:, :, :, keep]
        self.truncated = (float(e_ir[-1, 0]), float(e_ir[-1, 1]), self.join)
        e_ir[-1, 1] = self.join
        self.edges = np.vstack([e_ir, e_vi])
        self.k = np.concatenate([k_ir, k_vi], axis=3)
        self.logk = np.log10(np.maximum(self.k, 1e-200))

    def kvec(self, p_mbar: float, t_k: float, q: float) -> np.ndarray:
        """k per air molecule, (band, g), trilinear in log10 p, T and log10 q."""
        def bracket(grid, value):
            # A ONE-POINT axis is not an error and must not interpolate. Several
            # tables in the bundle fix their variable gas -- every PCM Studio
            # set does -- so `Q` has a single entry, and the general form below
            # would divide by `grid[1] - grid[0]` with no `grid[1]` to read.
            # numpy clips the index to -1, the subtraction gives zero, and the
            # whole k-vector comes back NaN with only a RuntimeWarning. Found by
            # using those tables for CLIM-42; the transmissions were silently
            # NaN until the warning was read.
            if len(grid) == 1:
                return 0, 0.0
            i = int(np.clip(np.searchsorted(grid, value) - 1, 0, len(grid) - 2))
            return i, (value - grid[i]) / (grid[i + 1] - grid[i])

        it, wt = bracket(self.T, t_k)
        ip, wp = bracket(self.p, math.log10(p_mbar))
        iq, wq = bracket(np.log10(self.Q), math.log10(q))
        def pair(index, weight, axis):
            """The one or two nodes this axis contributes, with their weights."""
            if len(axis) == 1:
                return ((index, 1.0),)
            return ((index, 1 - weight), (index + 1, weight))

        out = np.zeros(self.k.shape[3:])
        for a, wa in pair(it, wt, self.T):
            for b, wb in pair(ip, wp, self.p):
                for c, wc in pair(iq, wq, self.Q):
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


def run_checks(t376: CorrK, t1000: CorrK, sun: Spectrum, star: Spectrum,
               water_cm: float, t_k: float, co2_atmos_cm: float) -> list[str]:
    """The three things that could have failed, before any number is quoted.

    EVERY ONE OF THESE HAS A RIGHT ANSWER AND A BOUND DERIVED FROM IT. They
    used to print and return nothing, so "the three things that could have
    failed" could not fail: this file was the one sibling of
    `co2_overlap_589.py`, `trace_gas_band_model.py` and `cloud_band_weight.py`
    that did not raise. Each bound below is an inequality the physics or the
    arithmetic forces, not a number chosen to fit what the tables happen to
    give, and none of them has a tolerance in it.

    1. THE MIXING-RATIO SCALING. CO2 is trace and N2-broadened in both tables,
       so `k` is exactly proportional to the CO2 mixing ratio and the ratio of
       the two tables' `k` is exactly the ratio of their mixing ratios. The
       selection admits bands where another absorber contributes, and that can
       only DILUTE the ratio toward 1 -- it cannot push it above. So the median
       must lie in [1, ratio_expected]. A misread file layout, a table read at
       the wrong index, or the same table read twice all land outside it.

    2. THE BAND BINNING CONSERVES FLUX, FOR BOTH SPECTRA. The per-band fractions
       and a direct integration of the same spectrum over the table's full span
       are the same sum in two groupings, so they differ only by float64
       reassociation. The bound is `N * eps` over the band count, which is a
       rigorous upper bound on any association order. It ran on the solar
       spectrum alone, and every number this script exists to quote -- `h2osww`
       and `co2sww` are both ratios of a stellar broadband to a solar one --
       rests on the STELLAR partition, which nothing looked at. A hole in the
       star grid, or a star spectrum whose support stops short of the table
       span, gave a wrong ratio under "all three hold". world-60x0.

    3. THE RANDOM-OVERLAP DIVISION. `T_mix / T_dry` is what isolates H2O, and
       the two tables must return the same H2O absorptance from it. The error
       the division can introduce is bounded by the quantity it removes: the
       CO2 absorptance itself over the same path. That is rigorous and needs no
       tolerance -- a disagreement larger than the CO2 signal is not a residual
       of removing CO2.

    Returns the failures; the caller raises.
    """
    failed = []
    ratio_expected = co2_vmr_of(TABLE_1000) / co2_vmr_of(TABLE_376)
    it = int(np.argmin(abs(t376.T - t_k)))
    ip = int(np.argmin(abs(t376.p - math.log10(P_STANDARD_MBAR))))
    a, b = t376.k[it, ip, 0], t1000.k[it, ip, 0]
    # CO2-dominated bands only: where the ratio is meaningful the k must scale
    # EXACTLY with the mixing ratio, because CO2 is trace and N2-broadened in both.
    strong = a.max(axis=1) > 1e-27
    r = (b / np.maximum(a, 1e-300))[strong]
    co2ish = r > 2.0
    median = float(np.median(r[co2ish]))
    print(f"  k(1000 ppm)/k(376 ppm) where CO2 dominates: expected {ratio_expected:.4f}, "
          f"got {median:.4f} "
          f"[{np.percentile(r[co2ish], 2):.4f}, {np.percentile(r[co2ish], 98):.4f}] "
          f"over {co2ish.sum()} points")
    if not 1.0 <= median <= ratio_expected:
        failed.append(
            f"the k ratio's median is {median:.4f}, outside [1, "
            f"{ratio_expected:.4f}]. Dilution by another absorber can only pull "
            "it toward 1, so nothing physical puts it outside that range: the "
            "table layout is being read wrong")

    edges = sorted((float(lo), float(hi)) for lo, hi in t376.edges)
    lo_cm1, hi_cm1 = float(edges[0][0]), float(edges[-1][1])
    # THE HOLE DETECTOR STAYS EVEN THOUGH THE JOIN NOW CLOSES. It is the only
    # instrument that can see flux inside the span and inside no band, and it is
    # what found the 1974.95 to 2000 cm-1 sliver the old round-number cut left.
    # Reporting the count rather than assuming zero is the whole point.
    holes = [(edges[i][1], edges[i + 1][0]) for i in range(len(edges) - 1)
             if edges[i + 1][0] > edges[i][1]]
    tlo, thi, tjoin = t376.truncated
    print(f"  IR and VI joined at {t376.join:.6f} cm-1, the VI set's own first "
          f"edge; the IR band straddling it is cut from {tlo:.6f}-{thi:.6f} to "
          f"{tlo:.6f}-{tjoin:.6f} and keeps its own k")
    # BOTH SPECTRA, because both are broadband-weighted downstream. The star is
    # the one the quoted ratios divide by.
    for spec in (sun, star):
        f = flux_fractions(t376, spec)
        direct = spec.fraction_in_band(lo_cm1, hi_cm1)
        hole_flux = sum(spec.fraction_in_band(a, b) for a, b in holes)
        slack = len(f) * float(np.finfo(np.float64).eps)
        print(f"  flux fraction inside the {lo_cm1:g}-{hi_cm1:g} cm-1 table span, "
              f"{spec.label}: {f.sum():.5f} in {len(f)} bands + {hole_flux:.5f} in "
              f"{len(holes)} gap(s) = {direct:.5f} integrated in one piece "
              "(the rest is below 0.33 um, where neither gas absorbs)")
        if abs(f.sum() + hole_flux - direct) > slack:
            failed.append(
                f"for {spec.label} the per-band flux fractions and the gaps sum "
                f"to {f.sum() + hole_flux:.9f} and the same spectrum integrated "
                f"over the whole span gives {direct:.9f}. The band fractions are "
                "differences of one cumulative integral, so a contiguous "
                "partition telescopes exactly and can only differ by "
                f"{slack:.1e} of float64 reassociation: the band edges overlap, "
                "or the span is misread")
        if holes:
            failed.append(
                f"{len(holes)} window(s) carrying {hole_flux:.6f} of {spec.label}'s "
                "flux sit inside the table span and inside no band, so the "
                "correlated-k tables price them as transparent in both gas sets. "
                "The join is made at the VI set's first edge and should leave "
                f"none: {holes}")

    print("  T_mix/T_dry must return the same H2O absorptance from both tables:")
    for q in (1e-3, 1e-2, 1e-1):
        u = air_column_for_water(water_cm, q)
        vals, co2_signal = [], 0.0
        for tab in (t376, t1000):
            frac = flux_fractions(tab, sun)
            aa = 1.0 - tab.transmission(P_STANDARD_MBAR, t_k, q, u) / tab.transmission(
                P_STANDARD_MBAR, t_k, DRY, u)
            vals.append(broadband(frac, aa))
            # What the division had to remove, on this table's own CO2 amount.
            a_co2 = 1.0 - tab.transmission(
                P_STANDARD_MBAR, t_k, DRY,
                air_column_for_co2(co2_atmos_cm, co2_vmr_of(tab.name)))
            co2_signal = max(co2_signal, broadband(frac, a_co2))
        gap = abs(vals[1] - vals[0])
        print(f"    q={q:.0e}  {vals[0]:.6f} vs {vals[1]:.6f}  "
              f"({100*(vals[1]/vals[0]-1):+.2f}%), against a CO2 signal of "
              f"{co2_signal:.6f}")
        if gap > co2_signal:
            failed.append(
                f"at q={q:.0e} the two tables' H2O absorptance differs by "
                f"{gap:.6f}, more than the {co2_signal:.6f} of CO2 absorptance "
                "the division exists to remove. A residual cannot exceed what "
                "was removed, so the random-overlap division is not working")
    return failed


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

    # THE CHECKS RUN BEFORE ANY NUMBER IS QUOTED, always, which is what
    # "before any number is quoted" was supposed to mean. `--checks` used to
    # return here, so the only path that ran them was the one that computed
    # nothing, and every quoted number came out of a run that had checked
    # nothing. They cost about a second.
    print("CHECKS")
    failed = run_checks(t376, CorrK(TABLE_1000), sun, star, args.water_cm, t_k,
                        args.co2_earth)
    if failed:
        for line in failed:
            print(f"  FAIL: {line}")
        raise SystemExit(
            "the correlated-k tables are not being read as this script "
            "assumes, so nothing computed from them means anything")
    print("  all three hold\n")
    if args.checks:
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
        # This refit LANDED: radmod.f90 carries it. Read the four back out of the
        # source and say whether they still match, so a later edit to either side
        # shows up here as a disagreement rather than as nothing.
        in_model = radmod_co2_fit()
        same = all(abs(in_model[k] - fit[k]) <= 5e-5 * abs(fit[k]) for k in ("a1", "b1", "a2", "b2"))
        print("\nradmod.f90 carries " + " ".join(
            f"{k}={in_model[k]:.6g}" for k in ("a1", "b1", "a2", "b2")))
        print("  " + ("matches this refit, as PHYS-10 landed it"
                      if same else
                      "DOES NOT match this refit -- one of the two has moved"))
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
