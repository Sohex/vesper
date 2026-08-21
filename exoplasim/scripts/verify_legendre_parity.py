"""Check the parity relations Phase 3 rests on, against legini's own recurrence.

Phase 3 halves the inverse transforms by computing a mirror pair of latitudes
from one pass over the spectral modes. That is only valid if the weight
matrices obey

    qi(w, -mu) = s(w) * qi(w, +mu)        s(w) = (-1)^(m+n)
    qj(w, -mu) = -s(w) * qj(w, +mu)

i.e. P and its mu-derivative have OPPOSITE parity. Calculus says so; this
reproduces legini's recurrence verbatim and checks it numerically, because the
recurrence is what the model actually uses and a sign convention in it is
exactly the kind of thing that makes a plausible-looking parity split wrong.
"""
import math

NTRU = 10
NCSP = (NTRU + 1) * (NTRU + 2) // 2


def legini_recurrence(zsin):
    """zpli (P) and zpld (dP/dmu-ish) for one latitude, as legini computes them."""
    zpli = [0.0] * (NCSP + 1)
    zpld = [0.0] * (NCSP + 1)
    zcsq = 1.0 - zsin * zsin
    f1m = math.sqrt(1.5)
    zpli[1] = math.sqrt(0.5)
    zpli[2] = f1m * zsin
    zpld[1] = 0.0
    lm = 2
    for m in range(0, NTRU + 1):
        if m > 0:
            lm += 1
            f2m = -f1m * math.sqrt(zcsq / (m + m))
            f1m = f2m * math.sqrt(m + m + 3.0)
            zpli[lm] = f2m
            if lm < NCSP:
                lm += 1
                zpli[lm] = f1m * zsin
                zpld[lm - 1] = -m * f2m * zsin
        amsq = m * m
        for n in range(m + 2, NTRU + 1):
            lm += 1
            z1 = math.sqrt(((n - 1) ** 2 - amsq) / (4 * (n - 1) ** 2 - 1))
            z2 = zsin * zpli[lm - 1] - z1 * zpli[lm - 2]
            zpli[lm] = z2 * math.sqrt((4 * n * n - 1) / (n * n - amsq))
            zpld[lm - 1] = (1 - n) * z2 + n * z1 * zpli[lm - 2]
        if lm < NCSP:
            z3 = math.sqrt((NTRU * NTRU - amsq) / (4 * NTRU * NTRU - 1))
            zpld[lm] = -NTRU * zsin * zpli[lm] + (NTRU + NTRU + 1) * zpli[lm - 1] * z3
        else:
            zpld[lm] = -NTRU * zsin * zpli[lm]
    return zpli, zpld


def mode_index():
    """w -> (m, n), in legini's own enumeration order."""
    out = {}
    lm = 0
    for m in range(0, NTRU + 1):
        for n in range(m, NTRU + 1):
            lm += 1
            out[lm] = (m, n)
    return out


mu = 0.37384  # arbitrary interior latitude
p_pos, q_pos = legini_recurrence(mu)
p_neg, q_neg = legini_recurrence(-mu)
modes = mode_index()

bad_p = bad_q = 0
worst_p = worst_q = 0.0
for w in range(1, NCSP + 1):
    m, n = modes[w]
    s = 1.0 if (m + n) % 2 == 0 else -1.0
    # P should satisfy  P(-mu) = s * P(mu)
    ref = s * p_pos[w]
    if abs(p_pos[w]) > 1e-12:
        rel = abs(p_neg[w] - ref) / max(abs(ref), 1e-300)
        worst_p = max(worst_p, rel)
        if rel > 1e-12:
            bad_p += 1
    # Q should satisfy  Q(-mu) = -s * Q(mu)
    ref = -s * q_pos[w]
    if abs(q_pos[w]) > 1e-12:
        rel = abs(q_neg[w] - ref) / max(abs(ref), 1e-300)
        worst_q = max(worst_q, rel)
        if rel > 1e-12:
            bad_q += 1

print(f"modes checked: {NCSP}   (NTRU={NTRU}, mu={mu})")
print(f"P(-mu) == (-1)^(m+n) P(mu)   : {bad_p} violations, worst rel {worst_p:.2e}")
print(f"Q(-mu) == -(-1)^(m+n) Q(mu)  : {bad_q} violations, worst rel {worst_q:.2e}")
print()
print("PASS" if bad_p == 0 and bad_q == 0 else "FAIL -- the parity premise does not hold")
