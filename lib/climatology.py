"""Weights for the time bins of a pyburn climatology.

An annual mean taken from a binned `MOST.NNNNN.nc` is a mean over 12 bins, and
those bins do not all hold the same number of raw records. pyburn reduces
`ntimes` raw records to `nbin` bins with

    indices = np.linspace(0, ntimes, nbin+1).astype(int)
    counts  = np.diff(indices)
    binmean = np.add.reduceat(data, indices[:-1], axis=0) / counts

(`pyburn.py`, the `type(times)==int` branch). Each bin is divided by its own
count, so each bin mean is correct. What is not correct is averaging the bins
with EQUAL weight afterwards, because `linspace(...).astype(int)` distributes
the remainder of `ntimes / nbin` across the bins by truncation. Weight by the
counts and the annual mean is right.

## The two regimes, and which one is the awkward one

This is a property of the record count, so it differs between the model's two
I/O regimes, and the result is the opposite of what one would guess:

    NLOWIO = 1   36 raw records per orbit   36/12 = 3 exactly, bins all equal
    NLOWIO = 0  182 raw records per orbit   [15]*5 + [16] + [15]*5 + [16]

So the CHEAP spin-up regime bins exactly and the CLEAN regime a climatology is
required to use is the one whose bins are uneven. Two bins in twelve carry a
sixteenth record and are therefore under-weighted by 5.5% under equal weighting,
while the other ten are over-weighted by 1.1%.

Bin 0 is NOT the long one -- it holds the minimum. A derivation that assumed
eleven equal bins and took the twelfth as the residual reported a 1.19x excess
on bin 0 and a mechanism to explain it; both were artefacts of the assumption.
`TASKS.md` CLIM-13 carries the correction and the evidence.

## A separate, smaller truncation that weighting does not fix

At `NLOWIO = 0` the write interval is 32 timesteps and the orbit is 5850, so
182 records cover 5824 steps and the last 26 -- 0.44% of the orbit -- are never
written. No choice of weights recovers them; an annual mean from a binned file
is a mean over 99.56% of the orbit however it is taken.
"""

from __future__ import annotations

import numpy as np

# Candidate raw-record counts are bounded: fewer than one record per bin cannot
# produce the file, and the search only has to reach a few thousand.
MAX_NTIMES = 20000


def counts_for(ntimes: int, nbin: int) -> np.ndarray:
    """Records per bin, reproducing pyburn's own binning arithmetic."""
    return np.diff(np.linspace(0, ntimes, nbin + 1, True).astype(int))


def infer_ntimes(centres: np.ndarray) -> int:
    """How many raw records the binned file was reduced from.

    The count is not stored anywhere in the file, so it is RECOVERED from the
    bin centres and then checked, rather than assumed. pyburn sets each bin's
    timestamp to the mean of the raw timestamps in it, so with a fixed write
    interval the gap between consecutive centres is `(c_i + c_i+1)/2` records.
    Exactly one `ntimes` reproduces the observed pattern of gaps for a regular
    stream, and this returns it.

    This is a check that can fail, which is the point: if no candidate matches,
    the file was not produced by the binning this module models and its weights
    are not derivable here. Raising is then the honest answer.

    ONE CASE IS GENUINELY UNDER-DETERMINED and is reported as such rather than
    resolved by picking a fit. When the centres are evenly spaced, every
    multiple of `nbin` reproduces them and the record count cannot be recovered
    from the file. It does not matter: uniform spacing means `nbin` divides the
    record count, and every candidate then gives the SAME equal weights. So
    `nbin` is returned, and callers that want the true record count for its own
    sake have to get it from the run rather than from here.
    """
    centres = np.asarray(centres, dtype=float)
    nbin = len(centres)
    if nbin < 3:
        raise ValueError("need at least three bins to infer the record count")
    gaps = np.diff(centres)
    if np.allclose(gaps, gaps[0], rtol=1e-9, atol=1e-9):
        return nbin       # degenerate; see the docstring. Weights are equal.
    # Scale-free: compare the SHAPE of the gap pattern, since the write interval
    # in timesteps is the unknown common factor.
    shape = gaps / gaps.min()
    for ntimes in range(nbin, MAX_NTIMES + 1):
        counts = counts_for(ntimes, nbin)
        if counts.min() <= 0:
            continue
        predicted = 0.5 * (counts[:-1] + counts[1:]).astype(float)
        predicted /= predicted.min()
        if np.allclose(predicted, shape, rtol=1e-9, atol=1e-9):
            return ntimes
    raise ValueError(
        f"no record count in [{nbin}, {MAX_NTIMES}] reproduces the bin centres "
        f"{centres.tolist()}. This file was not written by pyburn's integer "
        "binning of a regular stream, so its bin weights cannot be derived here.")


def bin_weights(centres: np.ndarray) -> np.ndarray:
    """Normalised weights for an annual mean over a pyburn climatology's bins.

    Pass the file's `time` variable. Returns weights summing to 1, equal only
    when the record count divides the bin count.
    """
    counts = counts_for(infer_ntimes(centres), len(centres)).astype(float)
    return counts / counts.sum()


def annual_mean(field: np.ndarray, centres: np.ndarray) -> np.ndarray:
    """Time mean of `field` over its leading bin axis, correctly weighted."""
    w = bin_weights(centres)
    if field.shape[0] != len(w):
        raise ValueError(f"field has {field.shape[0]} bins, the time axis "
                         f"{len(w)}")
    return np.tensordot(w, np.asarray(field, dtype=float), axes=(0, 0))


def annual_mean_of(ds, name: str) -> np.ndarray:
    """`annual_mean` for a field read straight off an open climatology Dataset.

    The common case, and worth a name because it keeps the time axis and the
    field from being taken from different places. Raises through `bin_weights`
    if the file's bin centres are not something pyburn's integer binning could
    have produced, so a misapplied weighting fails rather than quietly scaling
    the wrong axis.
    """
    return annual_mean(np.asarray(ds[name][:], dtype=float),
                       np.asarray(ds["time"][:], dtype=float))


def masked_mean(field: np.ndarray, centres: np.ndarray,
                mask: np.ndarray) -> np.ndarray:
    """Mean over a SUBSET of the bins, weighted by records per bin.

    Seasonal composites take a few bins out of the twelve, and the weights have
    to be renormalised over the subset rather than reused whole. Separate from
    `annual_mean` because taking a subset of a normalised weight vector and not
    renormalising is the kind of thing that looks right.
    """
    w = bin_weights(centres)[mask]
    if w.sum() <= 0:
        raise ValueError("the mask selects no bins")
    w = w / w.sum()
    return np.tensordot(w, np.asarray(field, dtype=float)[mask], axes=(0, 0))
