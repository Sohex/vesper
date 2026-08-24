#!/usr/bin/env python3
"""Scale the atmospheric state in a restart, leaving its mean and its boundary alone.

    python exoplasim/scripts/scale_restart.py IN OUT --factor 0.5

Worldbuilding frame: a utility on the Vesper project's climate model. The
"atmosphere" here is a prognostic array in a checkpoint file.

WHY THIS EXISTS. `world-bxr` measures an energy sink in the dry adiabatic
dynamical core and needs to know which term carries it. Every candidate has a
different ORDER in the amplitude of the flow, so multiplying the state by a
factor and re-measuring separates them by an exponent that is predicted in
advance and can come out wrong:

    dE/dt ~ eps^1   a term linear in the state against the FIXED boundary,
                    which here means the orographic geopotential
    dE/dt ~ eps^2   the linear scheme -- the semi-implicit split, the time
                    scheme, the vertical discretisation -- or a product of the
                    state with the unscaled reference profile t0
    dE/dt ~ eps^3   the quadratic nonlinear terms, whose truncation or aliasing
                    is the hypothesis on world-bxr
    dE/dt ~ eps^4   the cubic terms, which this grid does NOT dealias:
                    NLON = 3*NTRU+1 dealiases products of two fields and the
                    sigma-coordinate tendencies carry products of three

The spectral SHAPE is untouched, so the tail fraction near the truncation is
identical across arms. That is deliberate: it is the one control that separates
"how much energy is up there" from "how much flow there is", and the two moved
together in every arm run before this.

WHAT IS SCALED, AND WHAT IS NOT.

  scaled whole      sz sd szm sdm      vorticity and divergence, whose (0,0)
                    mode is zero by construction; planetary vorticity is added
                    at runtime and is not in this file
  scaled anomaly    st sq stm sqm    the (0,0) coefficient of every level is
                    preserved, so the global mean temperature and the global mean
                    humidity do not move. Measured against `so`, neither field
                    carries any terrain signature at all -- R^2 below 0.005 at
                    every level -- so the global mean is the whole of what has to
                    be held.
  scaled residual   sp spm   ln(ps) is 96 percent OROGRAPHIC: regressed on `so`
                    across its coefficients it gives R^2 = 0.96, because surface
                    pressure is mostly a statement about how much atmosphere
                    stands above the terrain. Scaling it whole removes kilometres
                    of that atmosphere from over unchanged mountains, and the
                    model does not survive it: at 0.5 and at 0.25 the surface
                    layer reached negative absolute temperatures within seconds
                    and the run took SIGFPE on `log(z/z0)`. So the component
                    proportional to `so` is projected out and HELD, and only the
                    residual scales. The (0,0) coefficient is held with it.
  untouched         so   spectral orography: the boundary condition, and the
                    reference the ln(ps) projection is taken against.
                    sr   PUMA's restoration temperature, which PlaSim does not use.
                    everything else, including the surface, the ocean and the
                    accumulators. A dry adiabatic run exchanges nothing with them.

THE FORMAT is the one `compare_restarts.py` documents: sequential unformatted
Fortran, a 16-character name record then a data record, reals 8 bytes because
the model is compiled with -fdefault-real-8. Records are rewritten in place at
identical length, so nothing else in the file moves.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _paths  # noqa: F401
import restart_format  # noqa: E402

# The whole field scales: the (0,0) mode of vorticity and divergence is zero.
SCALE_WHOLE = ("sz", "sd", "szm", "sdm")
# The anomaly scales and the (0,0) mode of each level is held.
SCALE_ANOMALY = ("st", "sq", "stm", "sqm")
# The component proportional to the orography is held and the residual scales.
SCALE_RESIDUAL = ("sp", "spm")


def scale_record(body: bytes, factor: float, nrsp: int, hold_mean: bool,
                 hold_along: np.ndarray | None = None) -> bytes:
    if len(body) % 8:
        raise SystemExit("expected an 8-byte real record")
    vals = np.frombuffer(body, dtype="<f8").copy()
    if len(vals) % nrsp:
        raise SystemExit(f"record of {len(vals)} reals is not a multiple of "
                         f"NRSP={nrsp}; the truncation in this file is not what "
                         f"the header says")
    coefficients = []
    for base in range(0, len(vals), nrsp):
        # index 0 is the real part of (n=0,m=0) and index 1 its imaginary part,
        # which is identically zero. Holding both is holding the global mean.
        start = 2 if hold_mean else 0
        level = vals[base: base + nrsp]
        held = np.zeros(nrsp - start)
        if hold_along is not None:
            x = hold_along[start:]
            c = float(x @ level[start:] / (x @ x))
            held = c * x
            coefficients.append(c)
        # Written as an interpolation rather than as held + factor*(level - held)
        # so that factor 1.0 reproduces the input BIT FOR BIT: 1.0*level plus
        # 0.0*held is exact, where the round trip through the projection is not.
        # The held leading coefficients are never touched at all, for the same
        # reason.
        level[start:] = factor * level[start:] + (1.0 - factor) * held
    return vals.tobytes(), coefficients


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path)
    ap.add_argument("dest", type=Path)
    ap.add_argument("--factor", type=float, required=True,
                    help="multiplier on the flow; 1.0 rewrites the file unchanged")
    ap.add_argument("--only", default=None,
                    help="scale only these records, comma separated, out of the "
                         "groups above. For isolating which field a failure "
                         "comes from; the default scales the whole state")
    ap.add_argument("--report", type=Path, default=None,
                    help="write a JSON provenance record here")
    args = ap.parse_args()

    if args.factor <= 0:
        raise SystemExit("--factor must be positive; a sign flip is a different "
                         "experiment and this script does not claim it")

    raw = args.source.read_bytes()
    recs = restart_format.decode(raw, args.source)
    at = {rec.name: i for i, rec in enumerate(recs)}

    if "nrsp" not in at:
        raise SystemExit("no nrsp record: this is not a plasim restart")
    (nrsp,) = struct.unpack("<i", recs[at["nrsp"]].payload)

    if "so" not in at:
        raise SystemExit("no so record: the orographic component of ln(ps) "
                         "cannot be held, and scaling it whole destroys the run")
    so = np.frombuffer(recs[at["so"]].payload, dtype="<f8")

    touched = {}
    projections = {}
    wanted = None if args.only is None else {n.strip() for n in args.only.split(",")}
    for name in SCALE_WHOLE + SCALE_ANOMALY + SCALE_RESIDUAL:
        if wanted is not None and name not in wanted:
            continue
        if name not in at:
            raise SystemExit(f"record {name!r} is absent; refusing to scale a "
                             f"state this script does not recognise")
        rec = recs[at[name]]
        new, coeffs = scale_record(
            rec.payload, args.factor, nrsp,
            hold_mean=name in SCALE_ANOMALY + SCALE_RESIDUAL,
            hold_along=so if name in SCALE_RESIDUAL else None)
        recs[at[name]] = restart_format.Record(name=name, payload=new,
                                               offset=rec.offset)
        touched[name] = len(rec.payload) // 8
        if coeffs:
            projections[name] = coeffs

    args.dest.parent.mkdir(parents=True, exist_ok=True)
    out = restart_format.encode(recs)
    args.dest.write_bytes(out)

    digest = hashlib.sha256(out).hexdigest()
    payload = {
        "note": "amplitude-scaled restart for the world-bxr order test; "
                "exoplasim/scripts/scale_restart.py",
        "generated": datetime.now(timezone.utc).isoformat(),
        "source": str(args.source),
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "dest": str(args.dest),
        "dest_sha256": digest,
        "factor": args.factor,
        "only": args.only,
        "nrsp": nrsp,
        "scaled_whole": list(SCALE_WHOLE),
        "scaled_anomaly_mean_held": list(SCALE_ANOMALY),
        "scaled_residual_orography_held": list(SCALE_RESIDUAL),
        "orographic_projection": projections,
        "reals_per_record": touched,
    }
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"factor {args.factor:g}: {len(touched)} spectral records rewritten, "
          f"{sum(touched.values())} reals")
    print(f"  {args.dest}  sha256 {digest[:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
