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
  scaled anomaly    st sq sp stm sqm spm    the (0,0) coefficient of every level
                    is preserved, so the global mean temperature, the global
                    mean humidity and the global mean ln(ps) do not move. Scaling
                    the mean ln(ps) would change the atmosphere's MASS, and the
                    sink is reported per square metre of it.
  untouched         so   spectral orography: the boundary condition. Holding it
                    fixed is what makes the eps^1 arm mean something.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _paths  # noqa: F401

# The whole field scales: the (0,0) mode of vorticity and divergence is zero.
SCALE_WHOLE = ("sz", "sd", "szm", "sdm")
# The anomaly scales and the (0,0) mode of each level is held.
SCALE_ANOMALY = ("st", "sq", "sp", "stm", "sqm", "spm")


def records(path: Path):
    raw = path.read_bytes()
    out, pos = [], 0
    while pos < len(raw):
        (n,) = struct.unpack_from("<i", raw, pos)
        body = raw[pos + 4: pos + 4 + n]
        (m,) = struct.unpack_from("<i", raw, pos + 4 + n)
        if m != n:
            raise SystemExit(f"{path}: record markers disagree at byte {pos}")
        out.append((pos, n, body))
        pos += 8 + n
    return raw, out


def scale_record(body: bytes, factor: float, nrsp: int, hold_mean: bool) -> bytes:
    if len(body) % 8:
        raise SystemExit("expected an 8-byte real record")
    vals = list(struct.unpack(f"<{len(body)//8}d", body))
    if len(vals) % nrsp:
        raise SystemExit(f"record of {len(vals)} reals is not a multiple of "
                         f"NRSP={nrsp}; the truncation in this file is not what "
                         f"the header says")
    for base in range(0, len(vals), nrsp):
        # index 0 is the real part of (n=0,m=0) and index 1 its imaginary part,
        # which is identically zero. Holding both is holding the global mean.
        start = base + (2 if hold_mean else 0)
        for i in range(start, base + nrsp):
            vals[i] *= factor
    return struct.pack(f"<{len(vals)}d", *vals)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path)
    ap.add_argument("dest", type=Path)
    ap.add_argument("--factor", type=float, required=True,
                    help="multiplier on the flow; 1.0 rewrites the file unchanged")
    ap.add_argument("--report", type=Path, default=None,
                    help="write a JSON provenance record here")
    args = ap.parse_args()

    if args.factor <= 0:
        raise SystemExit("--factor must be positive; a sign flip is a different "
                         "experiment and this script does not claim it")

    raw, recs = records(args.source)
    names = {}
    for idx, (pos, n, body) in enumerate(recs):
        if n == 16:
            names[body.decode("ascii", "replace").strip()] = idx + 1

    nrsp_idx = names.get("nrsp")
    if nrsp_idx is None:
        raise SystemExit("no nrsp record: this is not a plasim restart")
    (nrsp,) = struct.unpack("<i", recs[nrsp_idx][2])

    out = bytearray(raw)
    touched = {}
    for name in SCALE_WHOLE + SCALE_ANOMALY:
        idx = names.get(name)
        if idx is None:
            raise SystemExit(f"record {name!r} is absent; refusing to scale a "
                             f"state this script does not recognise")
        pos, n, body = recs[idx]
        new = scale_record(body, args.factor, nrsp,
                           hold_mean=name in SCALE_ANOMALY)
        out[pos + 4: pos + 4 + n] = new
        touched[name] = len(body) // 8

    args.dest.parent.mkdir(parents=True, exist_ok=True)
    args.dest.write_bytes(bytes(out))

    digest = hashlib.sha256(bytes(out)).hexdigest()
    payload = {
        "note": "amplitude-scaled restart for the world-bxr order test; "
                "exoplasim/scripts/scale_restart.py",
        "generated": datetime.now(timezone.utc).isoformat(),
        "source": str(args.source),
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "dest": str(args.dest),
        "dest_sha256": digest,
        "factor": args.factor,
        "nrsp": nrsp,
        "scaled_whole": list(SCALE_WHOLE),
        "scaled_anomaly_mean_held": list(SCALE_ANOMALY),
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
