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
CLIM-13 carries the correction and the evidence.

## Evenly spaced bin centres do not mean evenly filled bins

pyburn stamps each bin with the mean of the raw timestamps in it, so with a
regular stream the gap between two centres is `(c_i + c_i+1)/2` write intervals.
Evenly spaced centres therefore say that every ADJACENT PAIR of bins holds the
same total, which is weaker than saying the bins are equal.
`linspace(...).astype(int)` gives counts that take only the two values `k` and
`k+1`, so a constant pair sum has exactly two shapes: every bin holds `k`, or
the counts alternate `k, k+1, k, k+1, ...`. The alternating shape is what
`ntimes = nbin * (2m+1) / 2` produces, for whole `m`, whenever `nbin` is even:
18 raw records in 12 bins, then 30, then 42, and so on upward. Odd `nbin` never
reaches it.

The two shapes are not merely hard to tell apart on the file. They are
IDENTICAL on it. 36 records written every 160 timesteps and 18 records written
every 320 both stamp the twelve centres 319, 799, 1279, ... exactly, because
the write interval is the unknown common factor and halving the record count
while doubling the interval leaves every bin mean where it was. So a file whose
centres are evenly spaced does not carry its own record count, and it does not
carry its weights either: the alternating shape weights neighbouring bins 1 to
2 where the equal shape weights them 1 to 1.

The collision is not an exotic corner. The one evenly binned regime this
project runs, `NLOWIO = 1` at 36 records in 12 bins, is exactly the case that
collides, with 18 records in 12 bins. Refusing it costs nothing beyond the
declaration, because a product binned from that regime is separately tainted --
`exoplasim/scripts/build_climatology.py` stamps it and its consumers refuse it
-- so the refusal here is the same answer arriving earlier.

`infer_ntimes` therefore enumerates the record counts an evenly spaced axis is
consistent with and asks whether they agree about the WEIGHTS, which is the
only thing a caller needs. Where they agree the weights are determined and are
returned even though the record count is not; that is what an odd `nbin` gives.
Where they disagree, which for an even `nbin` is always, it REFUSES. A caller
that knows the record count from the RUN -- from the write interval and the
orbit length in its namelist, which is where the number actually lives --
passes it as `ntimes` to `bin_weights` and the rest, and the declaration is
checked against the centres before it is used. Guessing on the caller's behalf
is what returned a factor of two here.

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


def _checked_axis(centres: np.ndarray) -> np.ndarray:
    """The bin centres, or a refusal saying which axis was handed over instead.

    Two shapes are refused here rather than deeper in, because both of them
    reach the record-count recovery as an ordinary array and come back out of it
    as an ordinary weight vector.

    A BIN INDEX is the first. `np.arange(nbin)` asserts that the bins are evenly
    spaced instead of measuring them, and every integer axis that starts at the
    time origin has that same shape. A real bin centre is the mean stamp of the
    records already written into that bin, so it lies strictly after the origin;
    a first centre of exactly zero is an index, not a time.

    A non-increasing axis is the second: the gaps between centres are what the
    record count is read out of, so an axis that repeats or reverses is not a
    time coordinate and nothing below can be trusted to fail on it.
    """
    centres = np.asarray(centres, dtype=float)
    nbin = len(centres)
    if nbin < 3:
        raise ValueError("need at least three bins to infer the record count")
    if np.all(centres == np.round(centres)) and centres[0] == 0.0:
        raise ValueError(
            f"the bin axis {centres.tolist()} is integers from zero, which is "
            "the shape of a bin INDEX and not of a time coordinate: a bin "
            "centre is the mean timestamp of the records in it and lies after "
            "the time origin. The weighting reads the spacing of real centres, "
            "so an index asserts that the bins are even instead of measuring "
            "them. Pass the file's own `time` variable")
    if not np.all(np.diff(centres) > 0):
        raise ValueError(
            f"the bin centres {centres.tolist()} do not increase, so they are "
            "not the time coordinate of a binned climatology")
    return centres


def _gap_shape(gaps: np.ndarray) -> np.ndarray:
    """The gap pattern with the write interval divided out.

    Scale-free, because the write interval in timesteps is the unknown common
    factor between the observed gaps and the record counts that produced them.
    """
    return gaps / gaps.min()


def _predicted_shape(counts: np.ndarray) -> np.ndarray:
    """The gap shape `counts` would stamp on the centres of a regular stream."""
    predicted = 0.5 * (counts[:-1] + counts[1:]).astype(float)
    return predicted / predicted.min()


def _even_axis_witness(centres: np.ndarray):
    """`(ntimes, counts)` for the first record count that contradicts an evenly
    spaced axis, or None when nothing does.

    Evenly spaced centres say that every ADJACENT PAIR of bins holds the same
    total, not that the bins are equal, so the file is consistent with a whole
    family of record counts. The family is ENUMERATED here rather than reasoned
    about, and the only question put to it is the one the caller needs answered:
    do its members all weight the bins the same way? The witness is a member
    that weights them differently from the equal binning, so its existence is
    what says the weights are not in the file. None says every consistent count
    weights the bins equally, and the weights are then determined even though
    the record count is not.

    For even `nbin` the witness always exists: `nbin` records fill the bins
    equally and `3*nbin/2` records fill them alternating, and both stamp this
    same axis.
    """
    nbin = len(centres)
    equal = counts_for(nbin, nbin)
    equal = equal / equal.sum()
    for ntimes in range(nbin + 1, MAX_NTIMES + 1):
        counts = counts_for(ntimes, nbin)
        if counts.min() <= 0:
            continue
        pair = counts[:-1] + counts[1:]
        if not np.all(pair == pair[0]):
            continue                      # would not stamp an even axis at all
        if not np.allclose(counts / counts.sum(), equal, rtol=0, atol=1e-12):
            return ntimes, counts
    return None


def infer_ntimes(centres: np.ndarray) -> int:
    """How many raw records the binned file was reduced from.

    The count is not stored anywhere in the file, so it is RECOVERED from the
    bin centres and then checked, rather than assumed. pyburn sets each bin's
    timestamp to the mean of the raw timestamps in it, so with a fixed write
    interval the gap between consecutive centres is `(c_i + c_i+1)/2` records.
    Where the gaps are not all equal, exactly one `ntimes` reproduces the
    observed pattern for a regular stream, and this returns it.

    This is a check that can fail, which is the point: if no candidate matches,
    the file was not produced by the binning this module models and its weights
    are not derivable here. Raising is then the honest answer.

    EVENLY SPACED CENTRES ARE THE UNDER-DETERMINED CASE, and this REFUSES them
    whenever the counts they are consistent with disagree about the weights,
    which for an even bin count is always. The module docstring has the
    collision: an equal binning and an alternating one stamp the identical axis,
    and they weight neighbouring bins 1 to 1 and 1 to 2. Callers that hold the
    record count independently pass it to `bin_weights` as `ntimes` rather than
    coming through here. Where the consistent counts all agree on the weights,
    which happens for an odd bin count, the answer is determined and `nbin` is
    returned; the record count itself is still not knowable from the file, so a
    caller that wants it for its own sake takes it from the run.
    """
    centres = _checked_axis(centres)
    nbin = len(centres)
    gaps = np.diff(centres)
    if np.allclose(gaps, gaps[0], rtol=1e-9, atol=1e-9):
        witness = _even_axis_witness(centres)
        if witness is None:
            return nbin       # every consistent count weights the bins equally
        ntimes, counts = witness
        raise ValueError(
            f"the {nbin} bin centres {centres.tolist()} are evenly spaced, "
            f"which says every adjacent PAIR of bins holds the same total and "
            f"NOT that the bins are equal. This file is consistent with any "
            f"multiple of {nbin} raw records, which fill the bins equally, and "
            f"equally consistent with {ntimes} records, which fill them "
            f"{counts.tolist()}. Both stamp exactly these centres, because the "
            f"write interval is the unknown common factor and halving the "
            f"record count while doubling the interval moves no bin mean. The "
            f"two weight neighbouring bins 1 to 1 and 1 to 2, so the weights "
            f"are not in this file any more than the record count is. Pass the "
            f"run's own record count as `ntimes`, from its write interval and "
            f"orbit length")
    shape = _gap_shape(gaps)
    for ntimes in range(nbin, MAX_NTIMES + 1):
        counts = counts_for(ntimes, nbin)
        if counts.min() <= 0:
            continue
        if np.allclose(_predicted_shape(counts), shape, rtol=1e-9, atol=1e-9):
            return ntimes
    raise ValueError(
        f"no record count in [{nbin}, {MAX_NTIMES}] reproduces the bin centres "
        f"{centres.tolist()}. This file was not written by pyburn's integer "
        "binning of a regular stream, so its bin weights cannot be derived here.")


def bin_counts(centres: np.ndarray, ntimes: int | None = None) -> np.ndarray:
    """Records per bin for a climatology, from its centres or a declaration.

    `ntimes` is the raw record count, and passing it is how a caller supplies
    what an evenly spaced axis cannot: the write interval and the orbit length
    are in the run's namelist, so the count is knowable there and only there.
    A declaration is CHECKED against the centres before it is used, so a count
    from the wrong run or the wrong I/O regime fails rather than reweighting the
    file to fit. What the check cannot do is separate the two shapes an evenly
    spaced axis is consistent with, since neither leaves a trace on the file;
    the declaration is the whole of the evidence there.
    """
    centres = _checked_axis(centres)
    nbin = len(centres)
    if ntimes is None:
        return counts_for(infer_ntimes(centres), nbin)
    ntimes = int(ntimes)
    counts = counts_for(ntimes, nbin)
    if counts.min() <= 0:
        raise ValueError(f"{ntimes} raw records cannot fill {nbin} bins: "
                         f"pyburn's binning leaves {counts.tolist()}")
    gaps = np.diff(centres)
    if not np.allclose(_predicted_shape(counts), _gap_shape(gaps),
                       rtol=1e-9, atol=1e-9):
        raise ValueError(
            f"a declared {ntimes} raw records in {nbin} bins would hold "
            f"{counts.tolist()} and stamp centres spaced "
            f"{_predicted_shape(counts).tolist()}, but this file's centres are "
            f"spaced {_gap_shape(gaps).tolist()}. The declaration and the file "
            "describe different runs")
    return counts


def bin_weights(centres: np.ndarray, ntimes: int | None = None) -> np.ndarray:
    """Normalised weights for an annual mean over a pyburn climatology's bins.

    Pass the file's `time` variable. Returns weights summing to 1, equal only
    when the record count divides the bin count. An evenly spaced axis raises
    unless `ntimes` says what the run wrote, because equal weights and
    alternating ones are both consistent with it.
    """
    counts = bin_counts(centres, ntimes).astype(float)
    return counts / counts.sum()


def annual_mean(field: np.ndarray, centres: np.ndarray,
                ntimes: int | None = None) -> np.ndarray:
    """Time mean of `field` over its leading bin axis, correctly weighted."""
    w = bin_weights(centres, ntimes)
    if field.shape[0] != len(w):
        raise ValueError(f"field has {field.shape[0]} bins, the time axis "
                         f"{len(w)}")
    return np.tensordot(w, np.asarray(field, dtype=float), axes=(0, 0))


def annual_mean_of(ds, name: str, ntimes: int | None = None) -> np.ndarray:
    """`annual_mean` for a field read straight off an open climatology Dataset.

    The common case, and worth a name because it keeps the time axis and the
    field from being taken from different places. Raises through `bin_weights`
    if the file's bin centres are not something pyburn's integer binning could
    have produced, or if they are evenly spaced and no record count was
    declared, so a misapplied weighting fails rather than quietly scaling the
    wrong axis.
    """
    return annual_mean(np.asarray(ds[name][:], dtype=float),
                       np.asarray(ds["time"][:], dtype=float), ntimes)


def masked_mean(field: np.ndarray, centres: np.ndarray,
                mask: np.ndarray, ntimes: int | None = None) -> np.ndarray:
    """Mean over a SUBSET of the bins, weighted by records per bin.

    Seasonal composites take a few bins out of the twelve, and the weights have
    to be renormalised over the subset rather than reused whole. Separate from
    `annual_mean` because taking a subset of a normalised weight vector and not
    renormalising is the kind of thing that looks right.
    """
    w = bin_weights(centres, ntimes)[mask]
    if w.sum() <= 0:
        raise ValueError("the mask selects no bins")
    w = w / w.sum()
    return np.tensordot(w, np.asarray(field, dtype=float)[mask], axes=(0, 0))
