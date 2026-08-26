#!/usr/bin/env python
"""Stephens's tuned tables against the analytic fits radmod evaluates instead.

WORLDBUILDING FRAME. Vesper is a simulated super-Earth. Everything below is a
property of the model's shortwave scheme or of the papers it is ported from.

WHAT THIS SETTLES. `radmod`'s computed-cloud branch does not read Stephens's
tuned tables. It evaluates three analytic fits of its own,

    beta1  = tswr1 sqrt(mu0)
    beta2  = tswr2 sqrt(mu0) / ln(3 + 0.1 tau)
    1-om0  = tswr3 mu0**2 ln(1000/tau)

at tswr1 = 0.077, tswr2 = 0.065 and tswr3 = 0.0055, the last of which this
project re-weights for a K dwarf through `model.cloud_absorption_scale`. Those
three are NOT Stephens (1978) Eqs. (11a) and (11b) -- the surface polynomials in
ln(tau) and mu0 whose coefficients are that paper's Table 1, and which Stephens,
Ackerman and Smith (1984) p. 687 withdrew as containing errors. They are a
separate and much coarser parameterization, so the 1984 withdrawal does not
reach them directly. What it does reach is the reference they should be
compared against: 1984 Table 1(a) replaces the 1978 single-scattering albedo
table outright, and 1984 Tables 1(b) and 1(c) are the 1978 backscatter tables
with, in the paper's words, only a couple of points smoothed in each.

So the question is not whether radmod evaluates a withdrawn formula. It is how
far its own fits sit from the tables, and whether that matters at the optical
depths this model's cloud layers reach.

Both tables are transcribed below from the scanned papers, and the transcription
is checked rather than trusted: the two backscatter tables must agree except at
a few points, and the two albedo tables must not.

Run it as `python exoplasim/scripts/stephens_tables_vs_fits.py`.
"""

from __future__ import annotations

import argparse
import json
import math

import numpy as np

import _paths  # noqa: F401  anchors every path on this file and adds lib/

import sensitivity  # noqa: E402
from paths import rel  # noqa: E402

DEFAULT_BRACKET = _paths.ANALYSIS / "cloud_optical_depth_bracket.json"
DEFAULT_OUT = _paths.ANALYSIS / "stephens_tables_vs_fits.json"

TAU = np.array([1, 2, 5, 10, 16, 25, 40, 60, 80, 100, 200, 500], dtype=float)
MU0 = np.array([1.0, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1])

# Stephens (1978) Table 2, p. 2126. Average values of 1 - omega.
COALB_1978 = np.array([
    [0.0225, 0.0277, 0.0312, 0.0346, 0.0398, 0.0346, 0.0364, 0.0294, 0.0173],
    [0.0225, 0.0250, 0.0284, 0.0294, 0.0312, 0.0294, 0.0260, 0.0191, 0.0104],
    [0.0208, 0.0208, 0.0208, 0.0208, 0.0191, 0.0173, 0.0156, 0.0104, 0.0052],
    [0.0173, 0.0173, 0.0156, 0.0147, 0.0139, 0.0104, 0.0085, 0.0064, 0.0028],
    [0.0156, 0.0139, 0.0121, 0.0121, 0.0104, 0.0087, 0.0062, 0.0045, 0.0017],
    [0.0121, 0.0104, 0.0099, 0.0087, 0.0076, 0.0069, 0.0045, 0.0035, 0.0017],
    [0.0104, 0.0087, 0.0073, 0.0064, 0.0055, 0.0052, 0.0035, 0.0023, 0.0014],
    [0.0069, 0.0069, 0.0055, 0.0045, 0.0045, 0.0035, 0.0026, 0.0017, 0.00052],
    [0.0069, 0.0052, 0.0045, 0.0043, 0.0035, 0.0035, 0.0026, 0.0017, 0.00017],
    [0.0052, 0.0035, 0.0035, 0.0029, 0.0026, 0.0017, 0.0017, 0.0012, 0.0],
    [0.0035, 0.0035, 0.0026, 0.0017, 0.0017, 0.0017, 0.0016, 0.00035, 0.0],
    [0.0017, 0.0017, 0.0016, 0.0014, 0.0010, 0.0009, 0.00017, 0.0, 0.0],
])

# Stephens, Ackerman and Smith (1984) Table 1(a), p. 688. The REPLACEMENT.
COALB_1984 = np.array([
    [0.0225, 0.0222, 0.0218, 0.0208, 0.0199, 0.0155, 0.0109, 0.0059, 0.0017],
    [0.0213, 0.0200, 0.0179, 0.0176, 0.0156, 0.0118, 0.0078, 0.0038, 0.0010],
    [0.0195, 0.0166, 0.0146, 0.0125, 0.0096, 0.0069, 0.0043, 0.0021, 0.0005],
    [0.0173, 0.0138, 0.0114, 0.0093, 0.0070, 0.0049, 0.0026, 0.0013, 0.0003],
    [0.0156, 0.0111, 0.0090, 0.0073, 0.0052, 0.0035, 0.0019, 0.0009, 0.0002],
    [0.0115, 0.0088, 0.0069, 0.0052, 0.0038, 0.0026, 0.0014, 0.0007, 0.00014],
    [0.0104, 0.0055, 0.00425, 0.0032, 0.0023, 0.00145, 0.0008, 0.0003, 0.0001],
    [0.0083, 0.0050, 0.0038, 0.0028, 0.0020, 0.0013, 0.0007, 0.00034, 0.0001],
    [0.0069, 0.0043, 0.0035, 0.0022, 0.0018, 0.0011, 0.0006, 0.0003, 0.0000],
    [0.0060, 0.0043, 0.0035, 0.0022, 0.0018, 0.0011, 0.0006, 0.0003, 0.0000],
    [0.0044, 0.0031, 0.0025, 0.0016, 0.0011, 0.00072, 0.0004, 0.00019, 0.0000],
    [0.0026, 0.0018, 0.0014, 0.0010, 0.00072, 0.00048, 0.00029, 0.00015, 0.0],
])

# Stephens (1978) Table 2. Average values of beta1 and beta2.
BETA1_1978 = np.array([
    [0.04213, 0.05569, 0.06569, 0.07690, 0.09323, 0.11113, 0.12951, 0.14069, 0.11959],
    [0.04723, 0.06146, 0.07077, 0.08028, 0.09236, 0.10172, 0.10767, 0.10337, 0.07942],
    [0.05820, 0.06924, 0.07440, 0.07818, 0.08153, 0.08107, 0.07756, 0.06798, 0.04833],
    [0.06805, 0.07262, 0.07368, 0.07333, 0.07231, 0.06853, 0.06263, 0.05265, 0.03594],
    [0.07337, 0.07376, 0.07283, 0.07072, 0.06803, 0.06309, 0.05642, 0.04646, 0.03103],
    [0.07681, 0.07443, 0.07225, 0.06906, 0.06531, 0.05975, 0.05264, 0.04270, 0.02807],
    [0.07914, 0.07489, 0.07188, 0.06797, 0.06361, 0.05753, 0.05013, 0.04021, 0.02611],
    [0.08048, 0.07518, 0.07171, 0.06740, 0.06267, 0.05633, 0.04877, 0.03885, 0.02505],
    [0.08119, 0.07536, 0.07166, 0.06715, 0.06223, 0.05576, 0.04811, 0.03819, 0.02456],
    [0.08202, 0.07566, 0.07172, 0.06699, 0.06188, 0.05528, 0.04753, 0.03761, 0.02414],
    [0.08307, 0.07626, 0.07212, 0.06720, 0.06192, 0.05517, 0.04733, 0.03738, 0.02406],
    [0.08743, 0.08000, 0.07553, 0.07028, 0.06465, 0.05755, 0.04935, 0.03916, 0.02617],
])

BETA2_1978 = np.array([
    [0.04769, 0.06274, 0.07338, 0.08552, 0.10224, 0.12002, 0.13789, 0.14652, 0.12066],
    [0.05371, 0.06903, 0.07873, 0.08863, 0.10033, 0.10904, 0.11326, 0.10647, 0.07940],
    [0.06597, 0.07688, 0.08172, 0.08497, 0.08709, 0.08018, 0.08012, 0.06877, 0.04741],
    [0.07585, 0.07931, 0.07949, 0.07812, 0.07567, 0.07048, 0.06287, 0.05156, 0.03390],
    [0.08007, 0.07867, 0.07657, 0.07323, 0.06893, 0.06257, 0.05432, 0.04335, 0.02768],
    [0.08068, 0.07592, 0.07241, 0.06778, 0.06245, 0.05552, 0.04713, 0.03676, 0.02288],
    [0.07700, 0.07001, 0.06557, 0.06030, 0.05453, 0.04758, 0.03958, 0.03022, 0.01839],
    [0.06994, 0.06214, 0.05750, 0.05222, 0.04661, 0.04012, 0.03291, 0.02475, 0.01479],
    [0.06341, 0.05555, 0.05104, 0.04604, 0.04076, 0.03481, 0.02832, 0.02111, 0.01249],
    [0.05340, 0.04605, 0.04198, 0.03755, 0.03295, 0.02787, 0.02246, 0.01658, 0.00966],
    [0.04149, 0.03526, 0.03187, 0.02826, 0.02458, 0.02061, 0.01645, 0.01199, 0.00684],
    [0.02506, 0.02083, 0.01861, 0.01629, 0.01398, 0.01154, 0.00902, 0.00637, 0.00324],
])

# Stephens et al. (1984) Tables 1(b) and 1(c): the 1978 tables, smoothed at a
# couple of points each and printed to four decimals.
BETA1_1984 = np.array([
    [0.0421, 0.0557, 0.0657, 0.0769, 0.0932, 0.1111, 0.1295, 0.1407, 0.1196],
    [0.0472, 0.0615, 0.0708, 0.0803, 0.0924, 0.1017, 0.1077, 0.1034, 0.0794],
    [0.0582, 0.0692, 0.0744, 0.0782, 0.0815, 0.0812, 0.0776, 0.0680, 0.0483],
    [0.0682, 0.0726, 0.0737, 0.0733, 0.0723, 0.0685, 0.0626, 0.0527, 0.0359],
    [0.0734, 0.0738, 0.0728, 0.0707, 0.0680, 0.0631, 0.0564, 0.0465, 0.0310],
    [0.0768, 0.0744, 0.0723, 0.0691, 0.0653, 0.0598, 0.0526, 0.0427, 0.0281],
    [0.0791, 0.0749, 0.0719, 0.0680, 0.0636, 0.0575, 0.0501, 0.0402, 0.0261],
    [0.0805, 0.0752, 0.0717, 0.0674, 0.0627, 0.0563, 0.0488, 0.0389, 0.0251],
    [0.0812, 0.0754, 0.0717, 0.0672, 0.0622, 0.0558, 0.0481, 0.0382, 0.0246],
    [0.0820, 0.0757, 0.0717, 0.0670, 0.0619, 0.0553, 0.0475, 0.0376, 0.0241],
    [0.0831, 0.0763, 0.0721, 0.0672, 0.0619, 0.0552, 0.0473, 0.0374, 0.0241],
    [0.0874, 0.0800, 0.0755, 0.0703, 0.0647, 0.0576, 0.0494, 0.0392, 0.0262],
])

BETA2_1984 = np.array([
    [0.0477, 0.0627, 0.0734, 0.0855, 0.1022, 0.1200, 0.1379, 0.1465, 0.1207],
    [0.0537, 0.0690, 0.0788, 0.0886, 0.1003, 0.1090, 0.1133, 0.1065, 0.0794],
    [0.0660, 0.0769, 0.0817, 0.0850, 0.0871, 0.0864, 0.0801, 0.0688, 0.0474],
    [0.0759, 0.0793, 0.0795, 0.0781, 0.0757, 0.0705, 0.0629, 0.0516, 0.0339],
    [0.0801, 0.0787, 0.0766, 0.0732, 0.0689, 0.0626, 0.0543, 0.0434, 0.0277],
    [0.0807, 0.0759, 0.0724, 0.0678, 0.0625, 0.0555, 0.0471, 0.0368, 0.0229],
    [0.0770, 0.0700, 0.0656, 0.0603, 0.0545, 0.0476, 0.0396, 0.0302, 0.0184],
    [0.0699, 0.0621, 0.0575, 0.0522, 0.0466, 0.0401, 0.0329, 0.0248, 0.0148],
    [0.0634, 0.0556, 0.0510, 0.0460, 0.0408, 0.0348, 0.0283, 0.0211, 0.0125],
    [0.0534, 0.0461, 0.0420, 0.0376, 0.0330, 0.0279, 0.0225, 0.0166, 0.0097],
    [0.0415, 0.0353, 0.0319, 0.0283, 0.0246, 0.0206, 0.0165, 0.0120, 0.0068],
    [0.0251, 0.0208, 0.0186, 0.0163, 0.0140, 0.0115, 0.0090, 0.0064, 0.0032],
])


def check_transcription() -> dict:
    """Checks that can fail, on the transcription and on the paper's claims.

    Stephens et al. (1984) p. 688-689 say Tables 1(b) and 1(c) are the 1978
    backscatter tables "smoothed slightly", and that "this smoothing has altered
    only a couple of points in each table". Two independent readings of a scan
    can both be wrong the same way, but they cannot be wrong in a way that keeps
    that claim true: a slipped digit anywhere in either table shows up here as a
    third or fourth altered point. The albedo tables must fail the same test,
    because replacing them is what the 1984 paper is for.
    """
    b1 = np.abs(BETA1_1984 - BETA1_1978) > 0.00006   # beyond 4-decimal rounding
    b2 = np.abs(BETA2_1984 - BETA2_1978) > 0.00006
    om = np.abs(COALB_1984 - COALB_1978) > 1e-9
    findings = {
        "beta1_points_altered": int(b1.sum()),
        "beta2_points_altered": int(b2.sum()),
        "coalbedo_points_altered": int(om.sum()),
        "coalbedo_points_total": int(om.size),
    }
    if findings["beta1_points_altered"] > 4 or findings["beta2_points_altered"] > 4:
        raise SystemExit(
            "the two backscatter tables differ at more points than the 1984 "
            "paper claims; the transcription is wrong: " + json.dumps(findings))
    if findings["coalbedo_points_altered"] < om.size // 2:
        raise SystemExit(
            "the two single-scattering albedo tables agree too widely to be a "
            "replacement; the transcription is wrong: " + json.dumps(findings))
    # The 1984 albedo table is monotone in mu0 at every optical depth and the
    # 1978 one is not. That is the substance of the revision, so it is asserted
    # rather than described.
    mono84 = bool(np.all(np.diff(COALB_1984, axis=1) <= 1e-12))
    mono78 = bool(np.all(np.diff(COALB_1978, axis=1) <= 1e-12))
    findings["coalbedo_1984_monotone_in_mu0"] = mono84
    findings["coalbedo_1978_monotone_in_mu0"] = mono78
    if not mono84 or mono78:
        raise SystemExit(
            "the monotonicity the 1984 revision introduces is not present in "
            "the transcription: " + json.dumps(findings))
    return findings


def fits(tswr1: float, tswr2: float, tswr3: float):
    """radmod's three analytic fits on the tables' own (tau, mu0) grid."""
    tau = TAU[:, None]
    mu = MU0[None, :]
    beta1 = tswr1 * np.sqrt(mu) * np.ones_like(tau)
    beta2 = tswr2 * np.sqrt(mu) / np.log(3.0 + 0.1 * tau)
    coalb = tswr3 * mu * mu * np.log(1000.0 / tau)
    return beta1, beta2, coalb


def reflect_conservative(beta, tau, mu0):
    """Stephens Eq. (1). Band 1, where cloud droplet absorption is negligible."""
    x = beta * tau / mu0
    return x / (1.0 + x)


def reflect_absorbing(beta, coalb, tau, mu0):
    """Stephens Eqs. (2) and (3). Returns (reflectance, absorptance), band 2."""
    om = 1.0 - coalb
    un = np.maximum(1e-12, coalb)
    uz = un + 2.0 * beta * om
    u = np.sqrt(uz / un)
    teff = np.minimum(50.0, np.sqrt(un * uz) * tau / mu0)
    e = np.exp(teff)
    r = (u + 1.0) ** 2 * e - (u - 1.0) ** 2 / e
    re = (u * u - 1.0) / r * (e - 1.0 / e)
    tr = 4.0 * u / r
    return re, 1.0 - re - tr


def interp_table(table, tau, mu0):
    """Bilinear in (tau, mu0) on the tables' grid, as the 1984 paper prescribes.

    p. 689: "the appropriate values of beta1, beta2 and omega0 were interpolated
    from the tables using standard bilinear interpolation techniques". Outside
    the grid the edge value is held, which is what the fits also do in effect.
    """
    t = np.clip(np.asarray(tau, dtype=float), TAU[0], TAU[-1])
    m = np.clip(np.asarray(mu0, dtype=float), MU0[-1], MU0[0])
    mu_asc = MU0[::-1]
    tab = table[:, ::-1]
    it = np.clip(np.searchsorted(TAU, t) - 1, 0, len(TAU) - 2)
    im = np.clip(np.searchsorted(mu_asc, m) - 1, 0, len(mu_asc) - 2)
    ft = (t - TAU[it]) / (TAU[it + 1] - TAU[it])
    fm = (m - mu_asc[im]) / (mu_asc[im + 1] - mu_asc[im])
    return ((1 - ft) * (1 - fm) * tab[it, im]
            + ft * (1 - fm) * tab[it + 1, im]
            + (1 - ft) * fm * tab[it, im + 1]
            + ft * fm * tab[it + 1, im + 1])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bracket", type=_paths.Path, default=DEFAULT_BRACKET,
                    help="the world-jgen bracket JSON, for this model's own "
                         "layer optical depths and insolation")
    ap.add_argument("--tswr1", type=float, default=0.077)
    ap.add_argument("--tswr2", type=float, default=0.065)
    ap.add_argument("--tswr3", type=float, default=None,
                    help="default: the run's value, which this project scales "
                         "for a K dwarf through cloud_absorption_scale")
    ap.add_argument("--out", type=_paths.Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    checks = check_transcription()

    bracket = json.loads(args.bracket.read_text(encoding="utf-8"))
    tswr3 = args.tswr3
    if tswr3 is None:
        tswr3 = 0.0055 * 1.192  # cloud_absorption_scale, exoplasim/analysis
    beta1_fit, beta2_fit, coalb_fit = fits(args.tswr1, args.tswr2, tswr3)

    # Where the fits sit against the tables, over the tables' whole grid.
    grid = {
        "beta1_fit_over_table": (beta1_fit / BETA1_1984).tolist(),
        "beta2_fit_over_table": (beta2_fit / BETA2_1984).tolist(),
        "coalbedo_fit_over_table_1984": (coalb_fit / np.maximum(1e-6, COALB_1984)).tolist(),
    }

    # Restricted to what this model's layers actually carry. The band-2 optical
    # depths come from the world-jgen bracket, which built them from the model's
    # own temperature, humidity and sigma grid.
    layers = bracket["layers"]
    tau_model = np.array([L["tau_band2_eq10b"] for L in layers])
    tau1_model = np.array([L["tau_band1_eq10a"] for L in layers])
    # A flux-weighted set of zenith cosines: weight mu0 by mu0 itself, which is
    # what the incident flux does.
    mu_set = np.linspace(0.05, 1.0, 20)
    mu_w = mu_set / mu_set.sum()

    per_layer = []
    for j, L in enumerate(layers):
        t1 = tau1_model[j]
        t2 = tau_model[j]
        d_r1 = d_r2 = d_a2 = d_a2_78 = 0.0
        for m, w in zip(mu_set, mu_w):
            b1f = args.tswr1 * math.sqrt(m)
            b2f = args.tswr2 * math.sqrt(m) / math.log(3.0 + 0.1 * t2)
            cof = tswr3 * m * m * math.log(1000.0 / max(1e-10, t2))
            b1t = float(interp_table(BETA1_1984, t1, m))
            b2t = float(interp_table(BETA2_1984, t2, m))
            cot = float(interp_table(COALB_1984, t2, m))
            co78 = float(interp_table(COALB_1978, t2, m))
            b2_78 = float(interp_table(BETA2_1978, t2, m))
            d_r1 += w * (reflect_conservative(b1t, t1, m)
                         - reflect_conservative(b1f, t1, m))
            rf, af = reflect_absorbing(b2f, cof, t2, m)
            rt, at = reflect_absorbing(b2t, cot, t2, m)
            r78, a78 = reflect_absorbing(b2_78, co78, t2, m)
            d_r2 += w * (rt - rf)
            d_a2 += w * (at - af)
            d_a2_78 += w * (at - a78)
        per_layer.append({
            "level": L["level"],
            "water_path_g_m2_mean": L["water_path_g_m2_mean"],
            "tau_band1": t1, "tau_band2": t2,
            "d_reflectance_band1_table_minus_fit": d_r1,
            "d_reflectance_band2_table_minus_fit": d_r2,
            "d_absorptance_band2_table_minus_fit": d_a2,
            "d_absorptance_band2_1984_minus_1978": d_a2_78,
        })

    # The flux consequence, on the same cover and insolation the world-jgen
    # bracket measured, and bracketed the same way: the biggest and smallest a
    # single layer carrying all of the model's cover could be worth.
    #
    # ONLY THE REFLECTANCE ENTERS. A top-of-atmosphere forcing is the change in
    # what leaves the top, and over a black surface that is the cloud's own
    # reflectance alone; the absorptance change moves flux between the cloud and
    # the ground without leaving the column. A real surface returns part of what
    # the cloud transmits, which reduces the magnitude, so these figures are the
    # upper bound rather than the estimate. The absorptance columns stay in the
    # table because they say where the two schemes disagree.
    cover = bracket["arms"]["uniform"]["cloud_fraction_total_mean"]
    f1 = bracket["arms"]["uniform"]["band1_incident_w_m2"]
    zsolar1 = bracket["constants"]["zsolar1"]
    f2 = f1 / zsolar1 * (1.0 - zsolar1)
    alpha = bracket["planetary_albedo"]

    wet = [p for p in per_layer if p["water_path_g_m2_mean"] >= 1.0]
    band1 = [-cover * f1 * p["d_reflectance_band1_table_minus_fit"] for p in wet]
    band2 = [-cover * f2 * p["d_reflectance_band2_table_minus_fit"] for p in wet]
    both = [b1 + b2 for b1, b2 in zip(band1, band2)]
    lo, hi = min(both), max(both)

    out = {
        "what": "Stephens's tuned tables against radmod's own analytic fits",
        "issue": "world-f9ig",
        "transcription_checks": checks,
        "coefficients": {"tswr1": args.tswr1, "tswr2": args.tswr2,
                         "tswr3": tswr3},
        "bracket_source": rel(args.bracket),
        "grid_ratios": grid,
        "per_layer": per_layer,
        "flux_w_m2_if_tables_replaced_fits": {
            "cover": cover, "band1_incident_w_m2": f1, "band2_incident_w_m2": f2,
            "per_layer_band1": band1, "per_layer_band2": band2,
            "bracket_w_m2": [lo, hi],
            "bracket_kelvin": [sensitivity.forcing_to_kelvin(lo, alpha),
                               sensitivity.forcing_to_kelvin(hi, alpha)],
        },
        "settles_it": "nothing here. A T21 pair, and only after world-jgen's "
                      "own pair has run, because tswr1 was tuned against the "
                      "optical depth world-jgen corrects",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")

    print(f"transcription  beta1 altered at {checks['beta1_points_altered']} "
          f"points, beta2 at {checks['beta2_points_altered']}, "
          f"1-omega at {checks['coalbedo_points_altered']} of "
          f"{checks['coalbedo_points_total']}")
    print(f"coefficients   tswr1 {args.tswr1} tswr2 {args.tswr2} "
          f"tswr3 {tswr3:.6g}")
    print("")
    print("  lev  W g/m2   tau1    tau2   dRe1     dRe2     dA2   dA2(84-78)")
    for p in per_layer:
        print("  %3d %8.3f %7.3f %7.3f %+8.4f %+8.4f %+8.4f %+8.4f"
              % (p["level"], p["water_path_g_m2_mean"], p["tau_band1"],
                 p["tau_band2"], p["d_reflectance_band1_table_minus_fit"],
                 p["d_reflectance_band2_table_minus_fit"],
                 p["d_absorptance_band2_table_minus_fit"],
                 p["d_absorptance_band2_1984_minus_1978"]))
    print("")
    print(f"flux           {lo:+.2f} to {hi:+.2f} W/m2 if the tables replaced "
          f"the fits, over the layers carrying cloud water")
    print(f"               {out['flux_w_m2_if_tables_replaced_fits']['bracket_kelvin'][0]:+.2f} "
          f"to {out['flux_w_m2_if_tables_replaced_fits']['bracket_kelvin'][1]:+.2f} K")
    print(f"written        {rel(args.out)}")


if __name__ == "__main__":
    main()
