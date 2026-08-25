"""The standard error of a mean whose samples are not independent.

Worldbuilding frame: the series here are per-orbit diagnostics of the Vesper
climate model. Nothing in this module is about the real world.

WHY THIS EXISTS. A window mean's standard error is the scatter divided by the
root of the SAMPLE COUNT only when the samples are independent. Consecutive
orbits of this model are not: the difference series between two arms of one
configuration has a lag-1 autocorrelation well above zero, and the model's own
non-overlapping window means move by more than the differences its comparison
test was being asked to resolve. Two instruments in this project took their
standard error from the raw count -- `exoplasim/scripts/compare_equilibria.py`
from the orbit-to-orbit scatter inside a window, and
`exoplasim/scripts/assess_convergence.py` through `curve_fit`, which scales its
covariance as though the residuals were independent -- and both therefore
reported a confidence they had not earned. A test that understates its own
error rejects pairs that differ by less than the model's own drift, and does it
with a sigma figure that looks decisive.

WHAT REPLACES IT. For a stationary series with memory, the variance of a mean
over n consecutive samples is

    var(mean) = sigma^2 * tau / n

where sigma^2 is the marginal variance and tau the INTEGRATED AUTOCORRELATION
TIME, tau = 1 + 2 * sum_{k>=1} rho_k. tau is 1 for independent samples and
(1+r)/(1-r) for an AR(1) process with lag-1 correlation r. The effective sample
size is n / tau: the number of independent samples the series is worth.

THE ESTIMATE NEEDS A LONGER SERIES THAN THE WINDOW IT CORRECTS. tau cannot be
recovered from a span comparable to tau itself -- there are not enough
independent samples in it to estimate anything -- so `integrated_time` reports
whether the span it was given can support its own answer, and every caller is
expected to carry that verdict rather than quietly use the number.
`RELIABLE_SPAN_MULTIPLE` is where that bar is set, and it was fixed before any
of this was run against a model series.

WHAT IT IS NOT. tau is a property of a STATIONARY series. A run still
approaching equilibrium has a trend, a trend inflates every autocorrelation
toward one, and the resulting tau describes the approach rather than the
variability. `stationary_enough` is the guard for that, and it is the caller's
job to apply it to the span before asking for tau.
"""

from __future__ import annotations

import numpy as np

# A tau estimated from a span shorter than this multiple of itself is not an
# estimate. Ten is the practical floor: below it the sum over lags is running
# out of independent samples faster than it is accumulating signal. Declared
# here before it was applied to any model series.
RELIABLE_SPAN_MULTIPLE = 10.0


def autocorrelation(series, max_lag: int | None = None) -> np.ndarray:
    """rho_k for k = 0 .. max_lag, normalised so rho_0 is 1.

    The biased (divide by n at every lag) estimator, deliberately: it is the
    one whose lag sum has bounded variance, and the unbiased form's tail blows
    up exactly where the integrated time is accumulated.
    """
    x = np.asarray(series, dtype=float)
    n = x.size
    if n < 2:
        raise ValueError("an autocorrelation needs at least two samples")
    if max_lag is None:
        max_lag = n - 1
    max_lag = int(min(max_lag, n - 1))
    d = x - x.mean()
    denom = float(np.dot(d, d))
    if denom == 0.0:
        # A constant series has no variability to correlate; every mean of it
        # is exact. Report independence rather than dividing by zero.
        out = np.zeros(max_lag + 1)
        out[0] = 1.0
        return out
    return np.array([float(np.dot(d[: n - k], d[k:])) / denom
                     for k in range(max_lag + 1)])


def integrated_time(series) -> dict:
    """tau, the effective sample size, and whether the span can support them.

    Geyer's initial monotone positive sequence. The lag sum is taken in
    adjacent PAIRS, Gamma_m = rho_2m + rho_2m+1, which is provably positive and
    decreasing for a reversible process; the sum is truncated at the first pair
    that is not, and the pairs are forced non-increasing on the way. Truncating
    on the first negative rho instead keeps a run of noise that happens to
    start positive, which is how a short series returns a confident tau it has
    no basis for.

    Returns tau (never below 1), the effective sample size n/tau, the lag-1
    correlation for reporting, and `reliable`: whether the span is at least
    `RELIABLE_SPAN_MULTIPLE` times the tau it just produced.
    """
    x = np.asarray(series, dtype=float)
    n = x.size
    rho = autocorrelation(x)
    gammas = []
    for m in range(1, (n - 1) // 2 + 1):
        g = rho[2 * m - 1] + rho[2 * m]
        if g <= 0.0:
            break
        gammas.append(g if not gammas else min(g, gammas[-1]))
    tau = max(1.0, 1.0 + 2.0 * float(sum(gammas)))
    return {"tau": tau,
            "effective_sample_size": n / tau,
            "lag1": float(rho[1]) if rho.size > 1 else float("nan"),
            "samples": int(n),
            "reliable": bool(n >= RELIABLE_SPAN_MULTIPLE * tau)}


def mean_standard_error(series, tau: float | None = None) -> float:
    """sigma * sqrt(tau / n) for the mean of `series`.

    `tau` is passed in wherever it was estimated over a longer span than this
    one, which is the normal case and the whole point: a window is corrected by
    an autocorrelation time measured on something long enough to carry it.
    """
    x = np.asarray(series, dtype=float)
    if x.size < 2:
        return float("nan")
    if tau is None:
        tau = integrated_time(x)["tau"]
    return float(np.std(x, ddof=1) * np.sqrt(max(tau, 1.0) / x.size))


def samples_for_standard_error(sigma: float, tau: float, target: float) -> float:
    """How many consecutive samples a mean needs to reach `target` standard error.

    The inverse of the relation above, and the arithmetic that prices a window:
    n = tau * (sigma / target)^2. It is the answer to "how long must the window
    be", and the cost of the memory is the factor tau -- everything else is
    what an independent series would already have needed.
    """
    if target <= 0 or not np.isfinite(sigma) or not np.isfinite(tau):
        return float("nan")
    return float(max(tau, 1.0) * (sigma / target) ** 2)


def stationary_enough(series) -> dict:
    """Whether a span is flat enough for its tau to mean anything.

    THE CRITERION, fixed before it was applied: the linear trend fitted across
    the span must move the series by less than the span's own standard
    deviation. A trend larger than the scatter is a state that is still moving,
    every lag correlation on it is pushed toward one, and the tau that comes
    out describes the approach and not the variability around it.

    Returns the fitted drift per sample, the excursion it implies across the
    span, the span's standard deviation, and the verdict.
    """
    x = np.asarray(series, dtype=float)
    n = x.size
    if n < 3:
        return {"stationary": False, "reason": "fewer than three samples",
                "drift_per_sample": float("nan"), "excursion": float("nan"),
                "scatter": float("nan")}
    slope = float(np.polyfit(np.arange(n, dtype=float), x, 1)[0])
    excursion = abs(slope) * (n - 1)
    scatter = float(np.std(x, ddof=1))
    return {"stationary": bool(excursion < scatter),
            "reason": ("" if excursion < scatter else
                       f"the fitted trend moves the series {excursion:.4g} "
                       f"across the span against a scatter of {scatter:.4g}, "
                       "so this span is still approaching rather than varying"),
            "drift_per_sample": slope, "excursion": excursion,
            "scatter": scatter}
