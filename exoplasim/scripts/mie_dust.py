"""Bohren & Huffman Mie code plus lognormal size integration.

Written to produce the six numbers ExoPlaSim's `aerofile` actually consumes
(Qext, Qsca, Qsca*beta per shortwave band), and the mass extinction efficiency
used to anchor them.

Validated two ways before use:
  1. Bohren & Huffman's own worked example, m = 1.55, x = 5.213 -> Qsca = 3.10543
  2. OPAC's published mineral-transported component at 0.55 um
Neither is a free parameter; both must reproduce or the code is wrong.
"""

import numpy as np


def bhmie(x, m, nang=2):
    """Qext, Qsca, g for a sphere of size parameter x and refractive index m = n + ik.

    Downward recurrence on the logarithmic derivative, per Bohren & Huffman
    (1983) appendix A. Series truncated at the standard Wiscombe criterion.
    """
    y = m * x
    nmx = int(max(np.ceil(abs(y)), np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0)) + 16)
    nstop = int(np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0))

    # logarithmic derivative D_n(y), downward recurrence from a zero seed
    d = np.zeros(nmx + 1, dtype=complex)
    for n in range(nmx, 0, -1):
        d[n - 1] = (n / y) - 1.0 / (d[n] + n / y)

    psi0, psi1 = np.cos(x), np.sin(x)
    chi0, chi1 = -np.sin(x), np.cos(x)
    xi1 = complex(psi1, -chi1)
    xi0 = complex(psi0, -chi0)

    qext_sum = 0.0
    qsca_sum = 0.0
    g_sum = 0.0
    a_prev = b_prev = 0.0 + 0.0j

    for n in range(1, nstop + 1):
        fn = (2.0 * n + 1.0)
        psi = (2.0 * n - 1.0) * psi1 / x - psi0
        chi = (2.0 * n - 1.0) * chi1 / x - chi0
        xi = complex(psi, -chi)

        dn = d[n]
        an = ((dn / m + n / x) * psi - psi1) / ((dn / m + n / x) * xi - xi1)
        bn = ((dn * m + n / x) * psi - psi1) / ((dn * m + n / x) * xi - xi1)

        qext_sum += fn * (an.real + bn.real)
        qsca_sum += fn * (abs(an) ** 2 + abs(bn) ** 2)

        # asymmetry parameter, Bohren & Huffman eq. 4.62
        if n > 1:
            nn = n - 1.0
            g_sum += (nn * (nn + 2.0) / (nn + 1.0)) * (
                a_prev.real * an.real + a_prev.imag * an.imag
                + b_prev.real * bn.real + b_prev.imag * bn.imag)
        g_sum += ((2.0 * n + 1.0) / (n * (n + 1.0))) * (
            an.real * bn.real + an.imag * bn.imag)

        a_prev, b_prev = an, bn
        psi0, psi1 = psi1, psi
        chi0, chi1 = chi1, chi
        xi1 = complex(psi1, -chi1)

    qext = (2.0 / x ** 2) * qext_sum
    qsca = (2.0 / x ** 2) * qsca_sum
    g = (4.0 / (x ** 2 * qsca)) * g_sum
    return qext, qsca, g


def lognormal_integrate(lam_um, n, k, r_mod_um, sigma_g,
                        r_min_um=0.005, r_max_um=60.0, rho_g_cm3=2.6,
                        n_r=600):
    """Distribution-averaged Qext, Qsca, g and mass extinction efficiency.

    r_mod is the NUMBER-median radius, sigma_g the geometric standard deviation.
    Returns MEE in m^2/g. Truncation limits matter: OPAC's components are
    truncated, and reproducing its numbers requires reproducing its limits.

    A thin wrapper on `distribution_integrate`, which takes any dN/dlnr. The two
    are kept apart because a lognormal is the only shape OPAC uses and this is
    the signature its published numbers are reproduced through, while sea salt
    is emitted as a sum of three lognormals and cannot use it.
    """
    r = np.logspace(np.log10(r_min_um), np.log10(r_max_um), n_r)
    lnsig = np.log(sigma_g)
    # dN/dlnr for a lognormal, normalised later so the constant drops out
    dndlnr = np.exp(-0.5 * (np.log(r / r_mod_um) / lnsig) ** 2)
    return distribution_integrate(lam_um, n, k, r, dndlnr, rho_g_cm3)


def distribution_integrate(lam_um, n, k, r_um, dndlnr, rho_g_cm3):
    """As `lognormal_integrate`, for an arbitrary number distribution.

    `r_um` is a radius grid and `dndlnr` the number per unit ln(radius) on it,
    to any normalisation: every quantity returned is a ratio of integrals over
    the same weight, so the constant divides out. The grid must be fine enough
    to resolve the Mie ripple structure it is averaging over, which for a
    logarithmic grid means a few hundred points per decade.
    """
    r = np.asarray(r_um, dtype=float)
    w = np.asarray(dndlnr, dtype=float)

    m = complex(n, k)
    qext = np.empty_like(r)
    qsca = np.empty_like(r)
    gg = np.empty_like(r)
    for i, rr in enumerate(r):
        x = 2.0 * np.pi * rr / lam_um
        qext[i], qsca[i], gg[i] = bhmie(x, m)

    area = np.pi * r ** 2
    lnr = np.log(r)
    # cross-sections per particle, integrated over the distribution
    c_ext = np.trapezoid(qext * area * w, lnr)
    c_sca = np.trapezoid(qsca * area * w, lnr)
    c_g = np.trapezoid(gg * qsca * area * w, lnr)
    vol = np.trapezoid((4.0 / 3.0) * np.pi * r ** 3 * w, lnr)

    qext_eff = c_ext / np.trapezoid(area * w, lnr)
    ssa = c_sca / c_ext
    g_eff = c_g / c_sca
    # Mass extinction efficiency. c_ext is um^2, vol is um^3, and a density in
    # g/cm^3 is 1e-12 g/um^3 and the um^2 -> m^2 conversion is another 1e-12;
    # the two cancel exactly, so the simplified quotient below is m^2/kg-ready
    # m^2/g as-is.
    mee = c_ext / (vol * rho_g_cm3)
    return qext_eff, ssa, g_eff, mee


if __name__ == "__main__":
    # --- test 1: Bohren & Huffman worked example -------------------------
    # r = 0.525 um at lam = 0.6328 um. Compute x rather than using the rounded
    # 5.213 the book prints: that rounding alone shifts Qsca by 4e-4 and reads
    # as a code error when it is a transcription of the input.
    x_bh = 2.0 * np.pi * 0.525 / 0.6328
    qe, qs, g = bhmie(x_bh, complex(1.55, 0.0))
    print(f"B&H m=1.55 x={x_bh:.6f}:  Qsca={qs:.5f}   (expect 3.10543)"
          f"   {'OK' if abs(qs - 3.10543) < 1e-4 else 'FAIL'}")
    # non-absorbing sphere: Qext must equal Qsca exactly
    print(f"  Qext-Qsca = {qe - qs:+.2e} (must be 0 for k=0)"
          f"   {'OK' if abs(qe - qs) < 1e-12 else 'FAIL'}")

    # --- test 2: OPAC mineral-transported at 0.55 um ---------------------
    # OPAC MITR: r_mod 0.50 um, sigma 2.20, r 0.02-5.0 um, rho 2.6 g/cm3
    # n = 1.53, k = 0.0055 at 0.55 um
    qe, ssa, g, mee = lognormal_integrate(0.55, 1.53, 0.0055, 0.50, 2.20,
                                          r_min_um=0.02, r_max_um=5.0)
    print(f"OPAC MITR 0.55um:  MEE={mee:.4f} m2/g (expect 0.3685)"
          f"   ssa={ssa:.4f} (expect 0.8374)   g={g:.4f} (expect 0.7736)")
