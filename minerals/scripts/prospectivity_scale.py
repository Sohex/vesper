"""Scaling shared by both prospectivity scripts: hosts in, 0-1 field out.

Ore prospectivity on Vesper is an index over an invented planet's rock map. Both
`build_prospectivity.py` and `build_downstream_prospectivity.py` build a raw
score per deposit type and then express it on 0-1, and this module is the one
copy of the two pieces they had each written for themselves.

## The divisor is a property of the RULE, not of a build

A field divided by `field[land].max()` is divided by ONE CELL. A build that
produces a single unusually favourable cell rescales every other cell downward,
so the same terrain under a different generation, or the same rule after a
carve, yields numbers that cannot be compared with the ones before it -- and
nothing in the artifact says that is what happened.

The divisor here is instead the rule's ATTAINABLE MAXIMUM: the host weight
multiplied by every modifier evaluated at the top of its declared input range.
It is computed in closed form from the config, so it moves when a weight is
edited and at no other time, and 0.5 means "half of what this rule can award"
on every build. Each config carries an `input_ranges` block giving the top of
each modifier's range and saying what makes it a bound.

The consequence is deliberate: a field no longer reaches 1.0 unless some cell
actually attains the rule's maximum. `orogenic_au` needs a schist cell that is
simultaneously a saturated fold belt, at saturated stress and deeply exhumed.
Where no cell does, the top of the field is below 1 and that is the information
the land-maximum divisor was destroying.

## Why the ceiling is checked rather than trusted

`normalise` raises when a score exceeds the ceiling it was handed. A closed form
that disagrees with the expression it is supposed to bound is a defect in one of
the two, and it has an unambiguous right answer, so it is a test rather than a
diagnostic. Each script's `ceiling` function mirrors its own score expression
term for term; adding a modifier to one without the other fails here.
"""

from __future__ import annotations

import numpy as np

# Tolerance on the ceiling check, for float32 round-off in the mesh fields. A
# real disagreement between a closed form and its expression is a factor, not a
# last-place error.
CEILING_TOLERANCE = 1e-6


def normalise(field: np.ndarray, land: np.ndarray, ceiling: float,
              label: str = "") -> np.ndarray:
    """Express a raw score on 0-1 against the rule's attainable maximum.

    Ocean is left at zero. `ceiling` comes from the calling script's closed form
    over its config, never from this build's data.
    """
    ceiling = float(ceiling)
    if not ceiling > 0:
        raise ValueError(
            f"{label or 'rule'} has an attainable maximum of {ceiling}, so it "
            "can never award anything. A rule that cannot fire is a defect in "
            "the config, not a field of zeros.")
    if land.any():
        top = float(field[land].max())
        if top > ceiling * (1.0 + CEILING_TOLERANCE):
            raise ValueError(
                f"{label or 'rule'} scored {top:.6g} against a closed-form "
                f"attainable maximum of {ceiling:.6g}. The ceiling function and "
                "the score expression have diverged; fix whichever is wrong "
                "rather than raising the ceiling to match.")
    out = np.zeros_like(field, dtype=np.float32)
    if land.any():
        out[land] = np.clip(field[land] / ceiling, 0.0, 1.0).astype(np.float32)
    return out


def host_ceiling(hosts: dict[str, float] | None, label: str = "") -> float:
    """The largest weight any host can contribute.

    Taken over the DECLARED hosts rather than over the ones this export happens
    to carry a code for. A rock class missing from a build's table lowers what
    that build attains; it does not lower what the rule can award.
    """
    if not hosts:
        raise ValueError(
            f"{label or 'rule'} declares no hosts, so it has no host weight to "
            "scale. Give it hosts or give it a rule branch of its own.")
    return float(max(hosts.values()))


def host_weight(spec: dict, substrate: np.ndarray, basement: np.ndarray,
                codes: dict[str, int]) -> np.ndarray:
    """Best host weight per region, from substrate or basement.

    Both count, because a deposit hosted in the basement is still there when
    cover survives above it. Classes absent from this export's rock table are
    skipped: a build that does not carry a class cannot place a deposit in it.
    """
    score = np.zeros(substrate.shape, dtype=float)
    for code, weight in (spec.get("hosts") or {}).items():
        if code not in codes:
            continue
        score = np.maximum(score, weight * ((substrate == codes[code])
                                            | (basement == codes[code])))
    return score
