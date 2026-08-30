#!/usr/bin/env python3
"""Derive the design flux from the comfort-band scoring `docs/src/pipeline/state.md` section 5b describes.

WORLDBUILDING CONTEXT: Vesper is a fictional planet and this script scores a
toy GCM's simulated seasons to pick a stellar flux for it. Nothing here refers
to the real world.

WHY IT EXISTS. The 2026-08-16 derivation that set `baseline_flux_earth` lived
only as prose: no script, no artifact, no named runs, no stated thresholds
(CLIM-24, inherited-earth-constants.md finding 5). docs/src/pipeline/sequencing.md loop C requires the flux
to be RE-DERIVED on every new terrain, so the criteria have to be re-runnable,
not reinvented. This codifies them. Every number below is fixed in this
docstring before the script is run against the data, per the project's
threshold rule.

THE METHOD, as 5b records it: seasonal temperature per latitude band, measured
on converged runs at two fluxes spanning the target, projected across
candidate means, scored on a warm-season comfort band. Concretely here:

- A "month" is one of the 12 output bins; per-cell warmest and coldest bin
  means of `tas` are the two seasons. Bin weights come from
  `lib/climatology.py` (they differ between I/O regimes).
- The projection is per 10-degree band: the band's seasonal RESPONSE, in
  kelvin per unit flux ratio, taken between the two sources and applied to
  every land cell in the band. Bracket between two points, never extrapolate
  from one -- the doctrine 5b itself was written under.

RESPONSE, NOT AMPLIFICATION, and the difference is the whole of what the
projection assumes. The response is a signed regression coefficient through two
measured points; an amplification is a ratio to the global-mean move, and
calling it that smuggles in the premise that a band moves WITH the global mean.
This world breaks that premise: at 32 degrees of obliquity, high-latitude
summer land is warmed by melting out of snow under clear skies, and brightening
the star adds moisture and cloud and retreats the sea ice, coupling that land to
an ocean whose heat capacity damps the summer peak. Bands whose winter warms by
tens of kelvin have warmest bins that COOL. The ratio is still reported per
band, because 5b argues in it and it is the readable statement of the trade;
nothing projects with it, and its sign is a measurement rather than a fault.

The projection therefore reproduces both of its own anchor points at the band
level, which the ratio form did not: it converted flux to kelvin with
`lib/sensitivity.py`'s canonical slope while normalising the band responses by
the span the two sources actually showed, and the two disagreed, so the fitted
line missed one of the points it was fitted to. Measured 2026-08-30 on the
baseline climatology of run_67323a923013 against run_b45380e61f90: the sources
span -7.196 K over -0.055 in flux, a local 130.8 K per unit flux ratio, against
the canonical 159.7, so every band mean at the bracket point came out 22 per
cent too far. The canonical slope is still what the candidate range's kelvin
margins are converted through -- that is a different use and it is the right
one there -- and the disagreement between the two estimates of the same global
slope is reported rather than averaged away.

THE ONE PROPERTY A PROJECTED FIELD MUST HAVE. `warm` and `cold` are the maximum
and the minimum of the same twelve bins, so warm >= cold is an identity of the
field. A per-band response applied across a wide flux range can break it: a band
whose winter response is several times its summer one has a seasonal range that
contracts with flux and eventually closes. That is the property the old sign
test was standing in for, and the sign test was neither necessary nor sufficient
for it. Measured on the same pair: three bands with a NEGATIVE warm response
stay ordered across the whole declared candidate range, while `-60 to -50` and
`+50 to +60`, whose warm responses are POSITIVE, close inside it. So the scored
candidates are the ones the projection is a field on, that window is reported,
and an inversion at one of the two MEASURED fluxes is a refusal -- there the
band-mean response is being asked to describe cells whose seasonal range it does
not.
- The bracket run contributes `tas` only: a spin-up under the cheap I/O
  regime carries the historical first-record defect in wind and humidity, so
  humidity comes solely from the clean climatology.

BOTH POINTS ARE ON THE VEGETATED BRANCH, and the branch is the condition
rather than the stage. What the method needs is two converged climatologies at
two fluxes with the same `model.land_albedo_source`; having soil and lakes is a
different property and not one it reads. So the first point is the BEST
AVAILABLE climatology whose run was vegetated, which on a first pass is the
bootstrap: the ordering there is bootstrap run, bootstrap climatology, bracket
run, design flux, adopt the flux, derived surface fields, baseline run, and
REQUIRING a baseline would be circular -- the baseline is the run on the full
surface fields and it still has to run at some flux, so demanding it first
means buying it at a provisional number and buying it again at the derived one.
Once a baseline exists it is the better determined of the two points and this
takes it, checking its branch off its own run manifest exactly as the bracket
run's is checked.

WHICH HALF A RADIATION CHANGE REACHES, because the shortwave cloud optics have
just moved the modelled climate by up to 20.6 K and the two halves of this
docstring answer to different things. The comfort band is a set of thresholds in
DEGREES of simulated surface air temperature: a statement about the kind of
world this project wants, not about the model that produces it. Changing the
radiation moves which stellar flux delivers a given temperature field, so it
moves the CANDIDATE RANGE and leaves the band where it is.
`exoplasim/notes/trace-gas-absorbers.md` section 4 records the same mechanism
from the other side: a radiation scheme short of greenhouse forcing reaches
these same thresholds at a higher stellar flux, and the search finds it. So the
range below is re-derived under the corrected optics and the band is not.

DECLARED THRESHOLDS, each with its status.

    comfort:  warmest-bin mean <= 33.0 C  AND  coldest-bin mean >= -25.0 C
    harsh:    warm side 33 to 38 C; cold side -25 down to -40 C
    extreme:  above 38 C, or below -40 C

    score = land-area fraction in the comfort band; the winner maximizes it.

`cold_floor_c`, `warm_extreme_c` and `cold_extreme_c` are DECLARED DESIGN
PREFERENCES and are KEPT. They say what kind of world this is meant to be, they
are fixed here ahead of the runs that will be judged against them, and nothing
in the model or in any run determines them. A design preference does not move
when the physics does, so the cloud optics do not reach them.

`warm_ceiling_c` is the one that was NOT a preference. Its recorded
justification was that the anchor flux "puts the tropics near +33 C in their
warmest month rather than +36", which is a number read off where the superseded
physics landed -- the same shape as the cold-extreme cap CLIM-30 exists to undo.
The corrected optics remove that justification outright: at the anchor flux the
modelled tropics are now several kelvin colder, so the sentence no longer picks
out 33. The value is RE-DECLARED here as a preference in its own right, at the
same number, and it now stands on what the other three stand on: the
warm-season monthly mean above which land on this world is not meant to read as
comfortable. A re-declaration at the old value is exactly where a tuning would
hide, so it is not left bare -- it is declared WITH A BRACKET THAT GETS SWEPT,
`THRESHOLD_BRACKET` below, and the report carries the design flux each bracket
point returns, so a reader can see what the answer owes to the preference. That
is the third of the four dispositions `docs/src/practice/conventions.md` allows
under "No tuned values".

`cold_extreme_cap` is DECLARED, at 0.05, on 2026-08-26 and ahead of the
re-derivation it constrains. It is a design preference of the same kind as the
three above it -- the largest share of land this world accepts below the
cold-extreme threshold -- and nothing in the model or in any run determines it.
CLIM-30.

It is the one threshold here with a HISTORY of being fitted to its own answer:
the purged artifact inferred 0.0642 from the anchor's own row, which reproduced
0.945 from an unconstrained winner of 0.8725. So it does not stand bare. It is
swept in `THRESHOLD_BRACKET` alongside the two comfort thresholds, and the
report carries the design flux each bracket point returns, which is what makes
the answer's debt to the preference readable instead of asserted. `main` still
refuses on None, so a future caller cannot get back to the inferred state by
deleting the value; and `cap_the_prior_implies` stays in the report as a
diagnostic that must never be used as the cap.

THE CANDIDATE RANGE, re-derived because the cloud optics moved under it. The
range's only job is to contain the winner, and a search that cannot reach its
own answer returns the edge of its range instead. Both bounds are built from
where the winner can be and then given margin, and the refusal below is what
makes the construction falsifiable rather than merely generous.

- The corrected optics COOL. world-f9ig put the tuned tswr1/tswr2/tswr3 fits
  aside for Stephens, Ackerman and Smith (1984)'s own tables, and
  `exoplasim/analysis/stephens_tables_vs_fits.json` brackets the flux ratio that
  offsets that at +0.026 to +0.102. That is a magnitude bound over a black
  surface with the whole of the model's cover on one layer, so the true offset
  is at most the top of it and may sit below the bottom. world-jgen's band-1
  correction pushes the other way and is currently unquantified: its artifact
  was computed under the deleted band-1 fit and is worthless, not stale
  (world-sii1).
- CEILING. The highest the winner can sit is the superseded-physics winner
  under the inferred cap, 0.945, plus the whole of the offset bound, which is
  1.047. The margin above that is declared as TWICE THE WIDTH OF THE HARSH-WARM
  BAND, 10 K, and the flux offset follows from `lib/sensitivity.py`'s slope. It
  covers a declared cap stricter than the inferred one, which moves the winner
  up, and the terrain moving between builds.
- FLOOR. The lowest the winner can sit is the superseded-physics UNCONSTRAINED
  winner, 0.8725, which the corrected optics can only push up. The margin below
  that is declared as THE WHOLE WIDTH OF THE HARSH-COLD BAND, 15 K, at unit
  cold-season amplification -- and 5b records the winter amplification as
  several times the summer's, so that is conservative. It covers world-jgen's
  opposing warming and each carve iteration darkening and warming the world.
- BOTH MARGINS ARE KELVIN AND THE FLUX ENDS ARE DERIVED FROM THEM, because the
  argument for each is made against a band width in degrees. They were once two
  flux literals, and when the slope was re-derived from 202.0 to 159.7 the
  floor's margin silently fell to 13.2 K against the 15 K band it claims to
  clear, while the ceiling's fell to 1.7 times a width it claims to double.
  Neither raised. A range that stops meaning what it says is exactly the
  failure the construction below is meant to make impossible.
- STEP 0.0025, unchanged. Half a kelvin of global mean, finer than anything the
  projection resolves, so the winner is limited by the criteria and not by the
  grid.
- The range is deliberately wider than the bound requires, because the two
  costs are not symmetric. Too wide costs arithmetic on candidates that cannot
  win. Too narrow returns an edge, silently, and that is the failure this
  re-derivation exists to remove.

THE HUMIDITY-COUPLED VARIANT, because the audit showed the dry score stands in
for a quantity that couples temperature to humidity, and this world's aridity
is distributed by drainage rather than latitude. The variant scores the
isobaric equivalent temperature T_e = tas + (L/cp) q at the lowest level --
plain thermodynamics, no other content. Its ceiling is pinned by declaration:
chosen so the comfort fraction at the pin flux matches the dry score's, so the
two columns differ only in HOW the margin is distributed, and the measured
divergence is the winner shift and the per-band gap. The pin flux is a
convention and its only requirements are that it be fixed ahead of the run and
lie inside the candidate range; it is set to the anchor because that is a value
already written down, and its standing as a prior derivation plays no part.

THE ANCHOR IS A PRIOR, NOT A TARGET. `anchor_flux` is what
`config/planet.yaml` carries as PROVISIONAL. The report measures the distance to
it and nothing depends on that distance being small; under the corrected optics
agreement would be surprising, because the anchor was chosen under the tuned
cloud coefficients and the offset bound above says the answer has moved.

Checks that can fail, all of which REFUSE rather than report and carry on:
`cold_extreme_cap` must be declared; the anchor and the pin flux must be
interior to the candidate range; the two sources must span more than 5 K of
global mean and must sit at different fluxes; the projection must not put a
land cell's warmest bin below its coldest at either MEASURED flux; at least one
candidate must survive that same test; the recorded prior must be one of the
candidates that does; the declared cap must admit at least one candidate; and
no winner -- dry or equivalent, capped or not -- may sit on either end of the
SCORED range. The last of those was a reported flag in the artifact and is now
a refusal, because the range is derived and an edge winner means the derivation
was wrong.

A FLUX THAT WAS CHOSEN, and why it is recorded by THIS script. Most of the
refusals above are about the inputs -- a missing run, a bare-rock bracket point,
a bracket too weak to project across -- and each is fixed by supplying what is
missing. A refusal that no run at any flux would lift is a different thing: the
instrument has looked at the evidence and cannot answer, and on that ground
alone the flux may be CHOSEN. `--chosen` records the choice.

`CHOOSABLE_REFUSALS` CARRIES ONE GROUND: `cold_extreme_cap_admits_no_candidate`.
Its predecessor `band_amplification_sign` is gone -- a negative warm-season
response is a measurement the linear form represents exactly, and the refusal
that replaced it, the projection inverting a cell's seasons at a measured flux,
is answered by a better projection or by bracketing points rather than by a
design decision. The cap ground passes the membership test because its minimum
is taken over the whole scored range: no run at any flux moves it, so only a
different terrain or a re-declared preference can, and both are design
decisions. Decided 2026-08-30 under world-dnrr, on
`notes/audits/design-flux-two-point-response.md`.

    --chosen FLUX --evidence notes/audits/<the finding>.md

WHAT MAKES THAT A RECORD AND NOT A RUBBER STAMP. `--chosen` runs the whole
derivation first and writes nothing unless it REFUSES, on a ground declared in
`CHOOSABLE_REFUSALS` below and nowhere else. A derivation that answers is
adopted, not overridden: the mode refuses and says to take the winner. A
derivation that refuses on a missing input refuses here too, because the fix is
to supply the input. So a choice cannot be reached by pointing the script at an
absent run, by narrowing the bracket, or by disliking the answer, and the
refusal the record carries is the exception the run actually raised rather than
a sentence someone typed.

The record carries the three things a later reader needs and a number alone does
not give them: WHY the derivation refused, as the verbatim refusal with the
bands that caused it; WHAT the choice rests on, as a path under `notes/` whose
content is hashed into the record; and WHAT WOULD REOPEN IT, which is declared
against the refusal in this file and is not the caller's to write. `basis` is
`derived` or `chosen` on every artifact this script writes, so nothing reading
one has to infer which it is holding, and `scripts/check_consistency.py` refuses
a chosen record that is missing any of the three.

Writes `exoplasim/analysis/design_flux.json`. Registered as step
`design_flux`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml
from netCDF4 import Dataset
from numpy.polynomial.legendre import leggauss

from _paths import ANALYSIS, CONFIG, RUNS  # noqa: F401  (puts lib/ on sys.path)

import climatology  # noqa: E402
import sensitivity  # noqa: E402
from paths import best_available_climatology, rel  # noqa: E402

LATENT_OVER_CP = 2.501e6 / 1004.9   # K per unit specific humidity, plasim's constants

# The two ends the candidate range is built from, in the units their argument
# is made in. `CEILING_ANCHOR` is the superseded-physics winner under the
# inferred cap plus the whole of the Stephens offset bound; `FLOOR_ANCHOR` is
# the superseded-physics UNCONSTRAINED winner, which the corrected optics can
# only push up. Both are FLUX because both are where a winner was observed.
CEILING_ANCHOR_FLUX = 1.047
FLOOR_ANCHOR_FLUX = 0.8725

# The margins, in KELVIN, which is the unit the docstring argues them in. The
# harsh bands are 33 to 38 C and -25 to -40 C, so 5 K and 15 K wide.
CEILING_MARGIN_K = 10.0      # twice the harsh-warm band
FLOOR_MARGIN_K = 15.0        # the whole harsh-cold band, at unit amplification
CANDIDATE_STEP = 0.0025      # half a kelvin of global mean, finer than the projection resolves


def _candidate_grid() -> list[float]:
    """The flux grid the winner is searched over, from the kelvin margins.

    Rounded OUTWARD to the step so a slope change can only widen the range,
    never narrow it below what the margins ask for. A range that is too wide
    costs arithmetic on candidates that cannot win; one that is too narrow
    returns its own edge, silently, which is the failure this construction
    exists to remove.
    """
    import math
    import sensitivity                                  # noqa: E402
    slope = sensitivity.SLOPE_K_PER_FLUX_RATIO
    lo = FLOOR_ANCHOR_FLUX - FLOOR_MARGIN_K / slope
    hi = CEILING_ANCHOR_FLUX + CEILING_MARGIN_K / slope
    lo = math.floor(lo / CANDIDATE_STEP) * CANDIDATE_STEP
    hi = math.ceil(hi / CANDIDATE_STEP) * CANDIDATE_STEP
    n = int(round((hi - lo) / CANDIDATE_STEP)) + 1
    return [round(lo + CANDIDATE_STEP * i, 4) for i in range(n)]


DECLARED = {
    "warm_ceiling_c": 33.0,      # design preference, re-declared, swept below
    "cold_floor_c": -25.0,       # design preference
    "warm_extreme_c": 38.0,      # design preference
    "cold_extreme_c": -40.0,     # design preference
    "cold_extreme_cap": 0.05,    # design preference, declared 2026-08-26, swept below
    "band_degrees": 10.0,
    # THE CANDIDATE RANGE IS DERIVED, not declared, and the docstring above
    # says why in kelvin rather than in flux. Its two margins are stated
    # against the BAND WIDTHS -- the ceiling clears twice the harsh-warm band,
    # the floor clears the whole harsh-cold band -- so the decision is a
    # temperature and the flux offset is that temperature over the slope.
    #
    # It was written as two flux literals and both arguments broke silently
    # when `lib/sensitivity.py`'s slope was re-derived from 202.0 to 159.7: the
    # floor's margin fell to 13.2 K against a 15 K band it claims to clear, and
    # the ceiling's to 1.7 times a width it claims to double. Neither failed.
    # The range simply stopped meaning what it said, which is the whole of
    # `docs/src/pipeline/loops.md`'s corollary in one line.
    "candidates": _candidate_grid(),
    "anchor_flux": 0.945,
    "te_pin_flux": 0.945,
}

# The bracket the two comfort thresholds are declared with, swept so the report
# says what the design flux owes to the preference rather than asserting it does
# not. The warm ceiling's bracket is +-3 K, inside its own 5 K harsh band; the
# cold floor's is +-5 K, a third of its 15 K one.
THRESHOLD_BRACKET = {
    "warm_ceiling_c": [30.0, 31.5, 33.0, 34.5, 36.0],
    "cold_floor_c": [-30.0, -27.5, -25.0, -22.5, -20.0],
    # The cap is bracketed to a factor of two either side of the declaration and
    # not to a fixed width in land fraction, because it is a SHARE and the
    # quantity a reader wants bounded is how much stricter or looser a defensible
    # alternative could have been. The bottom end is below the 0.0642 the purged
    # artifact inferred, so the sweep spans the value this declaration replaces
    # and a reader can see what that inference was worth.
    "cold_extreme_cap": [0.025, 0.05, 0.10],
}


class Refusal(SystemExit):
    """A guard that fired on the EVIDENCE rather than on a missing input.

    A `SystemExit` subclass, so it reaches a command line exactly as every other
    refusal in this script does and the message is the same message. What the
    type adds is that `--chosen` can tell the two classes apart: a missing run
    or a bare-rock bracket point is fixed by supplying what is missing and is
    never grounds for choosing a flux, while a guard that has read the fields
    and found the method inapplicable to this world is.
    """

    def __init__(self, key: str, message: str, context: dict | None = None):
        super().__init__(message)
        self.key = key
        self.context = context or {}


# THE REFUSALS A CHOSEN FLUX MAY STAND ON, and what would reopen each. Declared
# here, ahead of any run, for the same reason the thresholds are: the question
# "is this refusal one a choice can rest on" must not be answered by whoever is
# looking at the answer they wanted. Membership is the whole test -- `--chosen`
# writes nothing for a refusal outside this dict, and adding one is an edit to
# this file with the argument beside it.
#
# THE TEST FOR MEMBERSHIP IS WHETHER MORE INPUT WOULD ANSWER. A bracket that
# spans too little global mean, an absent run, a point on the wrong land-albedo
# branch: each is a refusal the project can BUY its way out of, and choosing a
# flux instead of buying it is laundering. A refusal that no run at any flux
# would lift is a statement about the method against this world, and a design
# decision is then the only thing left that can settle the number.
#
# The value is the reopen condition, which travels into the record. It is stated
# here rather than taken from the caller because it is a property of the guard:
# the caller supplies the evidence for the number, never the terms on which the
# derivation could resume.
#
# ONE MEMBER. `band_amplification_sign` was the previous entry, and it was
# wrong on its own terms: the projection is a straight line through two
# measured points and is defined for either sign of its slope, so a negative
# warm-season response was never a refusal the derivation had to take. What
# replaced it, `projection_inverts_at_a_measured_point`, fails the membership
# test above -- a per-cell response, or two points that bracket the candidate
# instead of sitting inside it, would answer it -- so it is a plain input
# refusal and not a ground for a choice.
#
# `cold_extreme_cap_admits_no_candidate` passes the test, and the argument is
# a measurement rather than a preference: the guard compares the declared cap
# against the MINIMUM cold-extreme land fraction over the whole scored range,
# so the refusal is a property of the terrain under this star at this
# obliquity, and no run at any flux moves the minimum. More input cannot
# answer it; only a different world or a different preference can, and both
# of those are design decisions. Decided 2026-08-30, world-dnrr;
# notes/audits/design-flux-two-point-response.md carries the decision and the
# measurement it rests on.
CHOOSABLE_REFUSALS: dict[str, str] = {
    "cold_extreme_cap_admits_no_candidate": (
        "reopen when the smallest cold-extreme land fraction any scored "
        "candidate reaches falls to the declared cap or below -- which is a "
        "property of the terrain under this star, so in practice a future "
        "carve that removes cold-extreme land, or a re-declaration of the "
        "cap with a new argument. Until one of those happens, re-running the "
        "derivation at any flux returns this same refusal."
    ),
}


def manifest_branch(manifest: dict) -> str | None:
    """`model.land_albedo_source` out of a run manifest, or None if it has none."""
    return manifest.get("source_config", {}).get("model", {}).get(
        "land_albedo_source")


def albedo_branch_of(clim_file: Path) -> str:
    """Which land-albedo branch the run behind a climatology was on.

    THE BRANCH IS THE CONDITION HERE, not the stage. This projects between two
    flux points and the projection is only a flux difference if both points sat
    on the same land albedo; a bare-rock point makes it an albedo difference
    wearing a flux label. The bracket run is already refused on this ground, and
    the first point earns the same check once it can be a baseline rather than
    always the bootstrap.

    IT REFUSES RATHER THAN ASSUMING when it cannot tell. A climatology carries
    `vesper_run_id` and the run carries the manifest, so the answer is on disk
    whenever the run is; where it is not, the honest report is that the branch
    is unknown, because the alternative is a projection between two worlds that
    nothing in the output would mark as such.
    """
    with Dataset(clim_file) as ds:
        run_id = getattr(ds, "vesper_run_id", None)
    if not run_id:
        raise SystemExit(
            f"{rel(clim_file)} carries no `vesper_run_id`, so the land "
            "albedo branch it was run on cannot be read. Rebuild it with "
            "exoplasim/scripts/build_climatology.py, which stamps the run id.")
    manifest = RUNS / str(run_id) / "run_manifest.json"
    if not manifest.is_file():
        raise SystemExit(
            f"{rel(clim_file)} names {run_id}, whose manifest is not on "
            f"disk at {rel(manifest)}. A run archived to identity cannot say "
            "which land albedo branch it was on, and this projection is only "
            "valid between two points on the same branch.")
    branch = manifest_branch(json.loads(manifest.read_text(encoding="utf-8")))
    if not branch:
        raise SystemExit(
            f"{rel(manifest)} records no model.land_albedo_source. The design "
            "flux is defined on the vegetated branch and an unrecorded branch "
            "cannot be checked against it.")
    return str(branch)


def seasonal_fields(path: Path) -> dict:
    """Per-cell warmest/coldest bin means of tas, plus masks and weights."""
    with Dataset(path) as ds:
        tas = np.asarray(ds["tas"][:], dtype=float)
        w = climatology.bin_weights(np.asarray(ds["time"][:], dtype=float))
        lat = np.asarray(ds["lat"][:], dtype=float)
        nlat, nlon = tas.shape[1:]
        land = np.asarray(ds["lsm"][:], dtype=float).mean(axis=0) > 0.5
        q = (np.asarray(ds["hus"][:], dtype=float)[:, -1]
             if "hus" in ds.variables and ds["hus"].ndim == 4 else None)
    gw = leggauss(nlat)[1][::-1]
    area = np.broadcast_to(gw[:, None], (nlat, nlon)) / (2.0 * nlon)
    annual = np.tensordot(w, tas, axes=(0, 0))
    return {
        "warm": tas.max(axis=0), "cold": tas.min(axis=0),
        "annual_global": float((annual * area).sum()),
        "q_warm": (None if q is None else
                   q[np.argmax(tas, axis=0),
                     np.arange(nlat)[:, None], np.arange(nlon)[None, :]]),
        "lat": lat, "land": land, "area": area,
    }


def tail_mean_fields(run_dir: Path, orbits: list[int]) -> Path:
    """Average the tail orbits' binned files into one 12-bin scratch file.

    Plain mean across orbits per bin, which is what build_climatology does;
    written to the analysis directory's scratch space so seasonal_fields can
    read one path either way.
    """
    import tempfile
    acc, count = {}, 0
    for orbit in orbits:
        p = run_dir / f"MOST.{orbit:05d}.nc"
        with Dataset(p) as ds:
            for name in ("tas", "lsm", "time"):
                a = np.asarray(ds[name][:], dtype=float)
                acc[name] = acc.get(name, 0.0) + a
            lat = np.asarray(ds["lat"][:], dtype=float)
            lon = np.asarray(ds["lon"][:], dtype=float)
        count += 1
    out = Path(tempfile.gettempdir()) / f"design_flux_tail_{run_dir.name}.nc"
    with Dataset(out, "w") as ds:
        ds.createDimension("time", acc["tas"].shape[0])
        ds.createDimension("lat", len(lat))
        ds.createDimension("lon", len(lon))
        for name, dims in (("time", ("time",)), ("lat", ("lat",)), ("lon", ("lon",))):
            v = ds.createVariable(name, "f8", dims)
            v[...] = {"time": acc["time"] / count, "lat": lat, "lon": lon}[name]
        for name in ("tas", "lsm"):
            v = ds.createVariable(name, "f8", ("time", "lat", "lon"))
            v[...] = acc[name] / count
    return out


def band_index(lat: np.ndarray, width: float) -> np.ndarray:
    return np.floor((lat + 90.0) / width).astype(int)


def score(warm_c: np.ndarray, cold_c: np.ndarray, land: np.ndarray,
          area: np.ndarray, warm_ceiling: float,
          cold_floor: float | None = None) -> dict:
    d = DECLARED
    cold_floor = d["cold_floor_c"] if cold_floor is None else cold_floor
    la = float((area * land).sum())
    def frac(mask):
        return round(float((area * (land & mask)).sum() / la), 5)
    return {
        "comfort": frac((warm_c <= warm_ceiling) & (cold_c >= cold_floor)),
        "harsh_warm": frac((warm_c > warm_ceiling) & (warm_c <= d["warm_extreme_c"])),
        "extreme_warm": frac(warm_c > d["warm_extreme_c"]),
        "harsh_cold": frac((cold_c < cold_floor) & (cold_c >= d["cold_extreme_c"])),
        "extreme_cold": frac(cold_c < d["cold_extreme_c"]),
    }


def band_mean(field: np.ndarray, rows: np.ndarray, land: np.ndarray,
              area: np.ndarray) -> float:
    """The land-area weighted mean of a field over one latitude band."""
    sel = land[rows, :]
    a = area[rows, :]
    return float((field[rows, :] * a * sel).sum() / float((a * sel).sum()))


def band_response(base: dict, other: dict, bands: np.ndarray, nb: int,
                  f0: float, f1: float, d_global: float) -> tuple[dict, list]:
    """Each band's seasonal RESPONSE in kelvin per unit flux ratio, and the rows.

    Signed, and taken between the two points directly, so the projection built
    from it reproduces both of its own anchors at the band level. See RESPONSE,
    NOT AMPLIFICATION in this module's docstring: the ratio to the global-mean
    move is reported per band and is never what anything projects with, because
    a band whose winter response is several times its summer one has a warmest
    bin that can move the OTHER WAY and the ratio's sign then reads as a fault
    rather than as the measurement it is.
    """
    d_flux = f0 - f1
    response = {"warm": np.zeros(nb), "cold": np.zeros(nb)}
    rows_out = []
    for b in range(nb):
        rows = bands == b
        if not base["land"][rows, :].any():
            continue
        entry = {"band": f"{b*10-90:+d} to {b*10-80:+d}"}
        for season in ("warm", "cold"):
            t_base = band_mean(base[season], rows, base["land"], base["area"])
            t_other = band_mean(other[season], rows, base["land"], base["area"])
            response[season][b] = (t_base - t_other) / d_flux
            entry[f"{season}_c_at_{f0}"] = round(t_base - 273.15, 2)
            entry[f"{season}_c_at_{f1}"] = round(t_other - 273.15, 2)
            entry[f"{season}_k_per_unit_flux"] = round(response[season][b], 2)
            entry[f"amp_{season}"] = round((t_base - t_other) / d_global, 3)
        rows_out.append(entry)
    return response, rows_out


def band_projector(base: dict, bands: np.ndarray, response: dict, f0: float):
    """`projected(f)` and `inverted_land(f)` for one base field and response.

    `projected` returns the warm, cold and equivalent-warm land fields at a
    candidate flux in degrees Celsius. `inverted_land` is the identity guard:
    `warm` and `cold` are the maximum and the minimum of the SAME twelve bins,
    so a cell whose projected warmest bin sits below its coldest is not a cell
    the comfort score can read.
    """
    warm_response_cells = response["warm"][bands][:, None]
    cold_response_cells = response["cold"][bands][:, None]
    q_warm = base.get("q_warm")
    te_warm0 = (base["warm"] if q_warm is None
                else base["warm"] + LATENT_OVER_CP * q_warm)

    def projected(f):
        df = f - f0
        return (base["warm"] + warm_response_cells * df - 273.15,
                base["cold"] + cold_response_cells * df - 273.15,
                te_warm0 + warm_response_cells * df - 273.15)

    def inverted_land(f):
        warm, cold, _ = projected(f)
        return base["land"] & (warm < cold)

    return projected, inverted_land


def derive(bracket_run: str, tail_orbits: int = 10) -> dict:
    """The derivation, or the refusal that says it cannot be made on this world.

    Returns the report with `basis` set to `derived`. Raises `Refusal` when a
    guard fires on the FIELDS -- the class of refusal a chosen flux may stand on
    -- and a plain `SystemExit` when a declaration is wrong or an input is
    missing, which is the class that is fixed by supplying what is missing.
    """
    # Refuse before a single field is read, so nothing about these can be
    # contaminated by what the data turns out to say.
    cap = DECLARED["cold_extreme_cap"]
    if cap is None:
        raise SystemExit(
            "cold_extreme_cap is not declared. The comfort-maximizing rule alone "
            "does not settle this flux: the trade 5b records is tropics against "
            "poles, and 'polar margins severe but small' is a constraint on the "
            "cold-extreme land fraction that the prose never quantified. Set "
            "DECLARED['cold_extreme_cap'] to the largest fraction of land this "
            "world accepts below the cold-extreme threshold, in this file, before "
            "the runs it will be applied to exist. CLIM-30; it must not be "
            "inferred from the answer, which is what the purged artifact did.")
    edges = (DECLARED["candidates"][0], DECLARED["candidates"][-1])
    for name in ("anchor_flux", "te_pin_flux"):
        f = DECLARED[name]
        if f in edges or f not in DECLARED["candidates"]:
            raise SystemExit(f"{name} = {f} is not interior to the candidate "
                             f"range {edges[0]} to {edges[1]}; the pin and the "
                             "prior must both be scored on a row the search can "
                             "reach from either side")

    cfg = yaml.safe_load(CONFIG.read_text(encoding="utf-8"))
    # THE BEST AVAILABLE, and the branch decides it rather than the stage.
    # What this needs is two converged climatologies ON THE VEGETATED BRANCH,
    # which is a statement about `model.land_albedo_source` and not about
    # having soil and lakes. Requiring a baseline would be circular: the design
    # flux IS the mean flux target, and a baseline run is by this project's
    # vocabulary the run on the full surface fields, which still has to be run
    # at SOME flux, so that ordering buys the expensive run at a provisional
    # number and then buys it again at the right one. But once a baseline
    # EXISTS it is the better determined of the two points and pinning to the
    # bootstrap keeps the projection anchored on the earlier world. So: take
    # the best available and check its branch the same way the bracket run's
    # is checked below, because the branch is the condition and the stage
    # never was.
    base_path, base_stage = best_available_climatology()
    base_branch = albedo_branch_of(base_path)
    if base_branch != "vegetated":
        raise SystemExit(
            f"{rel(base_path)} is the {base_stage} climatology and its run ran "
            f"with model.land_albedo_source = {base_branch!r}. The design flux "
            "is defined on the VEGETATED branch and both flux points have to "
            "be on it: the bare-rock arm is a BOUND, not a world, and "
            "projecting between the two measures the albedo difference rather "
            "than the flux difference.")
    base = seasonal_fields(base_path)
    f0 = float(cfg["orbit"]["baseline_flux_earth"])

    run_dir = RUNS / bracket_run
    if not run_dir.is_dir():
        raise SystemExit(
            f"{rel(run_dir)} does not exist. --bracket-run names a run by its "
            "UUID; ask exoplasim/runs/INDEX.json what is on disk. A run whose "
            "payload has been archived to identity cannot supply a tas tail.")
    manifest = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    # THE BRANCH IS CHECKED, because the whole argument for reading the
    # bootstrap rather than the baseline is that both points are vegetated. A
    # bare-rock bracket point would make the projection a comparison between two
    # different worlds and the script would not notice.
    branch = manifest_branch(manifest)
    if branch != "vegetated":
        raise SystemExit(
            f"{bracket_run} ran with model.land_albedo_source = {branch!r}. "
            "The design flux is defined on the VEGETATED branch and both flux "
            "points have to be on it: the bare-rock arm is a BOUND, not a "
            "world, and projecting between the two measures the albedo "
            "difference rather than the flux difference.")
    f1 = float(manifest.get("stellar_flux_ratio",
               manifest.get("derived_parameters", {}).get("stellar_flux_ratio_earth")))
    last = max(int(p.stem.split(".")[1]) for p in run_dir.glob("MOST.[0-9]*.nc"))
    orbits = list(range(last - tail_orbits + 1, last + 1))
    other = seasonal_fields(tail_mean_fields(run_dir, orbits))

    d_global = base["annual_global"] - other["annual_global"]
    if abs(d_global) < 5.0:
        raise SystemExit(f"the two sources span only {d_global:.2f} K of global "
                         "mean; a bracket this weak cannot support a projection")

    d_flux = f0 - f1
    if abs(d_flux) < 1e-9:
        raise SystemExit(f"both sources sit at flux {f0}; a response cannot be "
                         "taken between two points that are the same point")

    bands = band_index(base["lat"], DECLARED["band_degrees"])
    nb = bands.max() + 1
    response, band_rows = band_response(base, other, bands, nb, f0, f1, d_global)

    slope = sensitivity.SLOPE_K_PER_FLUX_RATIO
    projected, inverted_land = band_projector(base, bands, response, f0)

    def inverted_fraction(f):
        return float((base["area"] * inverted_land(f)).sum()
                     / (base["area"] * base["land"]).sum())

    # THE IDENTITY THE PROJECTION HAS TO KEEP, and the one the sign test was
    # standing in for. `warm` and `cold` are the maximum and the minimum of the
    # SAME twelve bins, so warm >= cold is a property of the field and not a
    # preference about it. A candidate at which the projection inverts a cell's
    # seasons is not a candidate whose comfort fraction means anything, and the
    # score would return one anyway.
    #
    # It fires at a MEASURED point or it does not fire at that point at all,
    # which is what makes it a check rather than a threshold: at f0 the
    # projection IS the base field, and at f1 it is the base field carried by a
    # band-mean response, so an inversion there is the band-uniform assumption
    # failing on a cell whose seasonal range the band mean does not describe.
    for name, f in (("bracket", f1),):
        if inverted_land(f).any():
            raise Refusal(
                "projection_inverts_at_a_measured_point",
                f"the projection puts {inverted_land(f).sum()} land cells' "
                f"warmest bin below their coldest at the {name} point's own "
                f"flux {f}, which is {inverted_fraction(f):.3%} of land. A "
                "band-mean response cannot describe those cells' seasons, so "
                "the per-band form is wrong for this world rather than "
                "under-determined.",
                {"all_bands": band_rows,
                 "inverted_land_fraction_at_bracket": round(inverted_fraction(f), 5),
                 "sources": {
                     "base": {"path": rel(base_path), "flux": f0,
                              "climatology_stage": base_stage,
                              "land_albedo_source": base_branch},
                     "bracket": {"run": bracket_run, "flux": f1,
                                 "tail_orbits": orbits,
                                 "land_albedo_source": branch},
                     "global_mean_span_k": round(d_global, 3),
                     "slope_k_per_unit_flux": slope}})

    # The candidates the projection is still a field on. Linear in f per cell,
    # so the valid set is an interval; asserted rather than assumed, because a
    # gap would mean the arithmetic is not what this comment says it is.
    scored = [f for f in DECLARED["candidates"] if not inverted_land(f).any()]
    if not scored:
        raise SystemExit(
            "the projection inverts a land cell's seasons at every candidate "
            f"in {DECLARED['candidates'][0]} to {DECLARED['candidates'][-1]}. "
            "There is no flux at which this projection is a field, so there is "
            "nothing for the comfort score to read")
    first = DECLARED["candidates"].index(scored[0])
    if DECLARED["candidates"][first:first + len(scored)] != scored:
        raise SystemExit(
            "the set of candidates the projection is a field on is not "
            "contiguous, which a linear response cannot produce. The "
            "projection is not the arithmetic this script says it is")
    scored_edges = (scored[0], scored[-1])

    # the humidity-coupled ceiling, pinned so the pin-flux comfort fractions match
    pin_warm, pin_cold, pin_te = projected(DECLARED["te_pin_flux"])
    dry_pin = score(pin_warm, pin_cold, base["land"], base["area"],
                    DECLARED["warm_ceiling_c"])
    lo, hi = 20.0, 90.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        s = score(pin_te, pin_cold, base["land"], base["area"], mid)
        if s["comfort"] < dry_pin["comfort"]:
            lo = mid
        else:
            hi = mid
    te_ceiling = round(0.5 * (lo + hi), 3)

    table = []
    for f in scored:
        warm, cold, te_warm = projected(f)
        table.append({
            "flux": f,
            "dry": score(warm, cold, base["land"], base["area"],
                         DECLARED["warm_ceiling_c"]),
            "equivalent": score(te_warm, cold, base["land"], base["area"],
                                te_ceiling),
        })

    def winner(kind, cap=None):
        rows = [r for r in table if cap is None or r[kind]["extreme_cold"] <= cap]
        return max(rows, key=lambda r: r[kind]["comfort"])["flux"] if rows else None

    win_dry, win_te = winner("dry"), winner("equivalent")
    win_dry_capped, win_te_capped = winner("dry", cap=cap), winner("equivalent", cap=cap)
    if win_dry_capped is None or win_te_capped is None:
        # A TYPED refusal, because it is a ground a choice may stand on: the
        # minimum below is over the whole scored range, so no run at any flux
        # lifts it. CHOOSABLE_REFUSALS carries the membership argument and the
        # reopen condition.
        raise Refusal(
            "cold_extreme_cap_admits_no_candidate",
            f"the declared cold-extreme cap {cap} admits no candidate in "
            f"{scored_edges[0]} to {scored_edges[1]}, the range the projection "
            f"is a field on inside the declared {edges[0]} to {edges[1]}. The "
            f"smallest cold-extreme land fraction any scored candidate reaches "
            f"is {min(r['dry']['extreme_cold'] for r in table)}. This world "
            "cannot be placed at any flux the search carries without exceeding "
            "the cap; that is a result about the cap and the terrain together, "
            "and it is not resolved by moving either one after the fact",
            {"cap": cap,
             "scored_range": list(scored_edges),
             "minimum_cold_extreme_land_fraction": float(
                 min(r["dry"]["extreme_cold"] for r in table)),
             "cold_extreme_land_fraction_by_flux": {
                 str(r["flux"]): r["dry"]["extreme_cold"] for r in table}})

    # An edge winner is a refusal, not a flag. The candidate range is derived in
    # this file's docstring from where the corrected cloud optics can put the
    # balance; a winner sitting on an end means that derivation was wrong, and
    # the number the search would return is the end of the range rather than an
    # optimum. Widen the range and say why, in the docstring, before re-running.
    at_edge = sorted({n for n, f in (("dry", win_dry), ("equivalent", win_te),
                                     ("dry under the cap", win_dry_capped),
                                     ("equivalent under the cap", win_te_capped))
                      if f in scored_edges})
    if at_edge:
        raise SystemExit(
            f"the winner is on the end of the scored range for: "
            f"{', '.join(at_edge)}. The candidates the projection is a field "
            f"on, {scored_edges[0]} to {scored_edges[1]} inside the declared "
            f"{edges[0]} to {edges[1]}, do not contain their own answer, so "
            "what the search would return is an edge and not an optimum. "
            "Widen the declared range in the docstring, or bring the "
            "projection's validity window to where the answer is, and re-run")

    # What the answer owes to the two comfort preferences, measured rather than
    # asserted: the design flux each bracket point of each threshold returns,
    # the other threshold and the cap held at their declared values.
    def winner_at(warm_ceiling, cold_floor, cap_value):
        best_flux, best_comfort = None, -1.0
        for f in scored:
            warm, cold, _ = projected(f)
            s = score(warm, cold, base["land"], base["area"], warm_ceiling, cold_floor)
            if s["extreme_cold"] > cap_value:
                continue
            if s["comfort"] > best_comfort:
                best_flux, best_comfort = f, s["comfort"]
        return best_flux

    # A bracket point that admits no candidate returns None rather than raising:
    # the DECLARED cap admitting nothing is a refusal above, but a bracket point
    # doing so is a result about that point, and reporting it is the whole
    # purpose of sweeping. A reader must be able to tell "no candidate" from
    # "the same winner", so the two are not both spelled as a missing row.
    sensitivity_rows = {
        "warm_ceiling_c": {str(v): winner_at(v, DECLARED["cold_floor_c"], cap)
                           for v in THRESHOLD_BRACKET["warm_ceiling_c"]},
        "cold_floor_c": {str(v): winner_at(DECLARED["warm_ceiling_c"], v, cap)
                         for v in THRESHOLD_BRACKET["cold_floor_c"]},
        "cold_extreme_cap": {str(v): winner_at(DECLARED["warm_ceiling_c"],
                                               DECLARED["cold_floor_c"], v)
                             for v in THRESHOLD_BRACKET["cold_extreme_cap"]},
    }

    # Reported, never used as the cap: the smallest cold-extreme tolerance under
    # which the recorded prior would have been the constrained optimum. It is
    # what the purged artifact solved backwards, kept visible so the defect is
    # legible rather than repeated.
    anchor_row = next((r for r in table if r["flux"] == DECLARED["anchor_flux"]),
                      None)
    if anchor_row is None:
        raise SystemExit(
            f"the recorded prior {DECLARED['anchor_flux']} is not a flux the "
            "projection is a field on, so the distance to it cannot be "
            "measured on the same instrument as the winner. That is a result "
            "about the projection and the prior together and it is not fixed "
            "by dropping the comparison")
    prior_implied_cap = anchor_row["dry"]["extreme_cold"]

    report = {
        "generated": datetime.now(timezone.utc).isoformat(),
        # DERIVED OR CHOSEN, on every artifact this script writes. A consumer
        # must never have to infer which of the two it is holding, and
        # `check_consistency.py` refuses an artifact that does not say.
        "basis": "derived",
        "generator": "exoplasim/scripts/derive_design_flux.py",
        "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "declared": DECLARED,
        "threshold_bracket": THRESHOLD_BRACKET,
        "sources": {
            "base": {"path": rel(base_path), "flux": f0,
                     # WHICH STAGE the first flux point came from, so a
                     # design flux derived before any baseline existed is
                     # distinguishable from one re-derived after. lib/paths.py.
                     "climatology_stage": base_stage,
                     "land_albedo_source": base_branch,
                     "note": "the best available climatology on the vegetated "
                             "branch; the design flux can be derived before a "
                             "baseline run exists and is re-derived on one "
                             "once it does"},
            "bracket": {"run": bracket_run, "flux": f1,
                        "tail_orbits": orbits,
                        "land_albedo_source": branch,
                        "note": "tas only; its spin-up I/O regime historically "
                                "corrupted wind and humidity, so humidity comes "
                                "from the clean climatology alone"},
            "global_mean_span_k": round(d_global, 3),
            "slope_k_per_unit_flux": slope,
        },
        "bands": band_rows,
        "equivalent_ceiling_c": round(te_ceiling, 2),
        "equivalent_ceiling_note": "pinned so the pin-flux comfort fraction "
                                   "matches the dry score's; the divergence "
                                   "between columns is then distribution, not level",
        "candidates": table,
        # WHERE THE PROJECTION IS STILL A FIELD, derived and not declared: the
        # candidates over which no land cell's projected warmest bin sits below
        # its coldest. Reported because the scored set is narrower than the
        # declared range and a reader must not have to infer that from a row
        # count.
        "projection_validity": {
            "declared_candidates": [edges[0], edges[1]],
            "scored_candidates": [scored_edges[0], scored_edges[1]],
            "dropped": len(DECLARED["candidates"]) - len(scored),
            "inverted_land_fraction_at_declared_edges": [
                round(inverted_fraction(edges[0]), 5),
                round(inverted_fraction(edges[1]), 5)],
            "note": "warm and cold are the maximum and the minimum of the same "
                    "twelve bins, so warm >= cold is an identity of the field. "
                    "A candidate that breaks it is not scored. The window is "
                    "necessary and not sufficient: a projection can be far "
                    "outside the regime it was measured in and still ordered, "
                    "which is what the extrapolation figures below are for",
        },
        "design_flux": win_dry_capped,
        "design_flux_equivalent": win_te_capped,
        "winner_dry_uncapped": win_dry,
        "winner_equivalent_uncapped": win_te,
        "projection_extrapolation_k": {
            "sources_span": [min(f0, f1), max(f0, f1)],
            "beyond_sources_at_design_flux": round(
                slope * max(0.0, min(f0, f1) - win_dry_capped,
                            win_dry_capped - max(f0, f1)), 2),
            # In units of the span the responses were measured over, because
            # that is the number that says whether the line is being read or
            # merely continued. `lib/sensitivity.py` says of its own slope that
            # a local response must not be used outside the regime it was
            # measured in, and the per-band responses here are the same kind of
            # object.
            "spans_beyond_sources_at_design_flux": round(
                max(0.0, min(f0, f1) - win_dry_capped,
                    win_dry_capped - max(f0, f1)) / abs(f0 - f1), 2),
            "spans_beyond_sources_at_scored_edges": [
                round(max(0.0, min(f0, f1) - scored_edges[0],
                          scored_edges[0] - max(f0, f1)) / abs(f0 - f1), 2),
                round(max(0.0, min(f0, f1) - scored_edges[1],
                          scored_edges[1] - max(f0, f1)) / abs(f0 - f1), 2)],
            "note": "how far past the two measured points the band responses "
                    "are carried to reach the design flux, in global-mean "
                    "kelvin and in multiples of the measured flux span; zero "
                    "means the answer is interpolated",
        },
        "threshold_sensitivity": {
            "note": "the design flux each bracket point of a comfort threshold "
                    "returns, the other threshold and the cap held at their "
                    "declared values. This is the sweep the warm ceiling is "
                    "declared with, and it measures what the answer owes to a "
                    "preference instead of asserting it owes nothing",
            **sensitivity_rows,
        },
        "against_the_recorded_prior": {
            "prior": DECLARED["anchor_flux"],
            "distance": round(win_dry_capped - DECLARED["anchor_flux"], 4),
            "cap_the_prior_implies": prior_implied_cap,
            "note": "the prior is config/planet.yaml's PROVISIONAL "
                    "baseline_flux_earth, chosen under the tuned tswr1/tswr2/"
                    "tswr3 cloud coefficients world-f9ig deleted. Agreement with "
                    "it is not a success and disagreement is not a failure; the "
                    "distance is reported because the prior fixes a semi-major "
                    "axis that is compiled into the biosphere. cap_the_prior_"
                    "implies is the tolerance that would have made the prior "
                    "optimal, reported so the backwards-solved threshold stays "
                    "legible; it is not the declared cap and must never be used "
                    "as one.",
        },
    }
    return report


def chosen_record(flux: float, evidence: Path, bracket_run: str,
                  tail_orbits: int = 10) -> dict:
    """The record of a flux that was CHOSEN because the derivation refused.

    Runs the derivation first and writes nothing unless it raises a `Refusal`
    whose key is in `CHOOSABLE_REFUSALS`. Three things therefore cannot happen.
    A derivation that ANSWERS cannot be overridden -- its winner is the flux and
    this raises. A refusal about a missing or unsuitable INPUT cannot be
    converted into a choice, because the fix is to supply the input. And the
    refusal in the record is the exception the run actually raised, with the
    bands that caused it, rather than a sentence someone typed.

    The caller supplies the number and the evidence for it. The reopen condition
    is NOT the caller's: it is declared against the refusal in this file,
    because it is a property of the guard rather than of the decision.
    """
    if flux not in DECLARED["candidates"]:
        raise SystemExit(
            f"the chosen flux {flux} is not one of this derivation's "
            f"candidates, {DECLARED['candidates'][0]} to "
            f"{DECLARED['candidates'][-1]} in steps of 0.0025. A chosen flux "
            "still has to sit inside the search space the refusal was taken "
            "in; outside it the record cites an argument it is not part of")
    if not evidence.is_file() or not evidence.read_bytes().strip():
        raise SystemExit(
            f"--evidence {evidence} is not a file with content in it. A chosen "
            "flux rests on a finding that a later reader can go and read, and "
            "the record carries the path and the hash of what was read")
    if "notes" not in evidence.resolve().parts:
        raise SystemExit(
            f"--evidence {evidence} is not under a `notes/` directory. A "
            "finding with its evidence lives there; a chosen flux that cites "
            "anything else is citing a document that is not one")

    # Bound outside the handler because Python deletes the `as` name when the
    # except block ends, and the record is built from what was caught.
    caught: Refusal | None = None
    try:
        derived = derive(bracket_run, tail_orbits)
    except Refusal as refusal:
        caught = refusal
        if refusal.key not in CHOOSABLE_REFUSALS:
            raise SystemExit(
                f"the derivation refused at {refusal.key!r}, which is not a "
                f"ground a chosen flux may stand on. The choosable refusals "
                f"are {sorted(CHOOSABLE_REFUSALS)}, and what makes a refusal "
                f"one of them is that no run at any flux would lift it. This "
                f"one is: {refusal}") from refusal
    else:
        raise SystemExit(
            f"the derivation ANSWERED: {derived['design_flux']}. A flux is "
            "chosen when the instrument refuses, never when it disagrees, so "
            "there is nothing to record here. Adopt the derived winner into "
            "config/planet.yaml, or change a declared threshold in this file "
            "and say why")

    return {
        "generated": datetime.now(timezone.utc).isoformat(),
        "basis": "chosen",
        "generator": "exoplasim/scripts/derive_design_flux.py",
        "generator_sha256": hashlib.sha256(
            Path(__file__).read_bytes()).hexdigest(),
        "design_flux": flux,
        "chosen": {
            "flux": flux,
            "refused_at": caught.key,
            "refusal": str(caught),
            "reopens_when": CHOOSABLE_REFUSALS[caught.key],
            "evidence": rel(evidence),
            "evidence_sha256": hashlib.sha256(
                evidence.read_bytes()).hexdigest(),
            "measured": caught.context,
            "note": "the derivation was RUN and refused; this record exists "
                    "because it did. `refusal` is the exception it raised and "
                    "`measured` is what it had read when it raised. Re-run "
                    "this script without --chosen to reproduce it",
        },
        "declared": DECLARED,
        "threshold_bracket": THRESHOLD_BRACKET,
        "design_flux_equivalent": None,
        "winner_dry_uncapped": None,
        "winner_equivalent_uncapped": None,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--bracket-run", required=True,
                    help="run id supplying the second flux point's tas tail, "
                         "e.g. run_1a2b3c4d5e6f. REQUIRED, and there is no "
                         "default because there cannot be one: runs are named "
                         "by UUID, so any default is the id of one historical "
                         "run and nothing else. `exoplasim/runs/INDEX.json` "
                         "says what exists; the `bracket_run` step in "
                         "config/pipeline.yaml says what makes a run usable "
                         "here")
    ap.add_argument("--tail-orbits", type=int, default=10)
    ap.add_argument("--output", type=Path, default=ANALYSIS / "design_flux.json")
    ap.add_argument("--chosen", type=float, default=None,
                    help="record a flux that was CHOSEN rather than derived. "
                         "The derivation is run first and this writes nothing "
                         "unless it refuses on a ground listed in "
                         "CHOOSABLE_REFUSALS; a derivation that answers is "
                         "adopted, not overridden. Requires --evidence")
    ap.add_argument("--evidence", type=Path, default=None,
                    help="the finding under notes/ a chosen flux rests on. Its "
                         "path and hash go into the record")
    args = ap.parse_args()

    if args.chosen is None:
        if args.evidence is not None:
            raise SystemExit(
                "--evidence is only meaningful with --chosen. A derived flux "
                "carries its own evidence: the two climatologies, the declared "
                "thresholds and the candidate table are all in the report")
        report = derive(args.bracket_run, args.tail_orbits)
    else:
        if args.evidence is None:
            raise SystemExit(
                "--chosen requires --evidence. A number with no argument "
                "beside it is the thing this mode exists to prevent, not the "
                "thing it exists to write")
        report = chosen_record(args.chosen, args.evidence, args.bracket_run,
                               args.tail_orbits)

    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if report["basis"] == "chosen":
        chosen = report["chosen"]
        print(f"basis: CHOSEN, not derived. design flux {report['design_flux']}")
        print(f"the derivation refused at {chosen['refused_at']}: "
              f"{chosen['refusal']}")
        print(f"evidence: {chosen['evidence']}")
        print(f"reopens when: {chosen['reopens_when']}")
        print(f"wrote {rel(args.output)}")
        return

    edges = (DECLARED["candidates"][0], DECLARED["candidates"][-1])
    src = report["sources"]
    print(f"sources: {src['base']['flux']} ({src['base']['path']}) and "
          f"{src['bracket']['flux']} ({src['bracket']['run']}), "
          f"span {src['global_mean_span_k']:.2f} K")
    print(f"candidates {edges[0]} to {edges[1]}, cold-extreme cap "
          f"{DECLARED['cold_extreme_cap']}")
    print(f"design flux: dry {report['design_flux']}  equivalent "
          f"{report['design_flux_equivalent']}  (uncapped: dry "
          f"{report['winner_dry_uncapped']}, equivalent "
          f"{report['winner_equivalent_uncapped']})")
    print(f"recorded prior {DECLARED['anchor_flux']}, distance "
          f"{report['against_the_recorded_prior']['distance']:+.4f}")
    print(f"warm-ceiling bracket: "
          f"{report['threshold_sensitivity']['warm_ceiling_c']}")
    print(f"cold-floor bracket:   "
          f"{report['threshold_sensitivity']['cold_floor_c']}")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
