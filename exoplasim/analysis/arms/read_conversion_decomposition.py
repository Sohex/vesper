#!/usr/bin/env python
"""Read `CONVDECOMP` and `ENERGY FIXER` out of a run's `plasim_diag`. world-0ov.

WORLDBUILDING CONTEXT: Vesper is a fictional super-Earth and every number here
is a diagnostic of the climate model that simulates it, in W/m2 of modelled
atmosphere.

WHAT IT READS. `plasim.f90` writes one `CONVDECOMP` line per `ndiag` interval
carrying the reference conversion's decomposition -- `d02`, `d26`, `d27`, `cimp`,
`cvadv`, `ct`, `dt`, `ctm`, `ctp` -- and, when the energy fixer is on, one
`ENERGY FIXER applied W/m2, K/day` line carrying the correction it is applying
and the residual imbalance it is correcting against.

WHY IT READS THE DIAG AND NOT THE RESTART. `epilog` deallocates `adenergy` and
then writes it into the restart 213 lines later, so every run with `nenergy > 0`
segfaults after its last timestep. The diagnostics are already on disk when that
happens, so a crashed run's `plasim_diag` is complete and is the record these
arms are read from. The defect is reported in
`exoplasim/notes/energy-diagnostic-restart-defect.md`.

THE FIRST THREE PRINTS ARE DROPPED, which is the instrument every earlier
world-0ov arm used: the first step out of a restart carries a start-up transient
of order 250 W/m2 and the decomposition needs the transient out before its terms
mean anything.
"""
import argparse
import hashlib
import json
import pathlib
import re
import subprocess
import sys

COLUMNS = ["d02", "d26", "d27", "cimp", "cvadv", "ct", "dt", "ctm", "ctp"]
DROP_PRINTS = 3

CONV = re.compile(r"^\s*CONVDECOMP\s+\S+(?:\s+\S+)*?\s+(\d+)\s+(\d+)\s+" +
                  r"\s+".join([r"(-?[\d.]+E[+-]\d+)"] * len(COLUMNS)))
FIXER = re.compile(r"^\s*ENERGY FIXER applied W/m2, K/day\s+(\d+)\s+"
                   r"(-?[\d.]+E[+-]\d+)\s+(-?[\d.]+E[+-]\d+)\s+(-?[\d.]+E[+-]\d+)")


def parse(diag: pathlib.Path) -> tuple[list[dict], list[dict]]:
    conv, fixer = [], []
    for line in diag.read_text(encoding="utf-8", errors="replace").splitlines():
        if "CONVDECOMP" in line:
            m = CONV.match(line)
            if m:
                row = {"nstep": int(m.group(1)), "naccu": int(m.group(2))}
                row.update({c: float(m.group(3 + i)) for i, c in enumerate(COLUMNS)})
                conv.append(row)
        elif "ENERGY FIXER applied" in line:
            m = FIXER.match(line)
            if m:
                fixer.append({"nstep": int(m.group(1)),
                              "applied_w_m2": float(m.group(2)),
                              "accumulated_k_per_day": float(m.group(3)),
                              "residual_w_m2": float(m.group(4))})
    return conv, fixer


def summarise(conv: list[dict], fixer: list[dict]) -> dict:
    kept = conv[DROP_PRINTS:]
    if not kept:
        return {"prints": len(conv), "kept": 0}
    last = kept[-1]
    out = {
        "prints": len(conv),
        "kept_after_dropping_first_%d" % DROP_PRINTS: len(kept),
        "final_print_nstep": last["nstep"],
        # The running means the model itself carries: each print is already an
        # accumulation, so the LAST print is the arm's answer and a mean over
        # prints would average a decaying transient into it.
        "final": {c: last[c] for c in COLUMNS},
        "final_sink_d26_minus_d27": last["d26"] - last["d27"],
        "final_displacement_cimp_minus_ct": last["cimp"] - last["ct"],
        # The identity every earlier arm used to certify the arrays:
        # (ctp + ctm)/2 must equal cimp, which is sdp = 2 sdt - adm holding.
        "identity_ctp_plus_ctm_over_2": 0.5 * (last["ctp"] + last["ctm"]),
        "identity_residual_vs_cimp": 0.5 * (last["ctp"] + last["ctm"]) - last["cimp"],
        "trajectory": [{"nstep": r["nstep"],
                        "sink": r["d26"] - r["d27"],
                        "displacement": r["cimp"] - r["ct"]} for r in kept],
    }
    if fixer:
        out["fixer_prints"] = len(fixer)
        out["fixer_final"] = fixer[-1]
        out["fixer_trajectory"] = fixer[-8:]
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=pathlib.Path, action="append", required=True,
                        help="run directory carrying a plasim_diag; repeatable")
    parser.add_argument("--label", action="append", required=True,
                        help="what each --run is, in the same order")
    parser.add_argument("--out", type=pathlib.Path, required=True)
    args = parser.parse_args()
    if len(args.run) != len(args.label):
        raise SystemExit("--run and --label must come in pairs")

    arms = {}
    for run_dir, label in zip(args.run, args.label):
        diag = run_dir.resolve() / "plasim_diag"
        if not diag.exists():
            raise SystemExit(f"no plasim_diag in {run_dir}")
        conv, fixer = parse(diag)
        arms[label] = {
            "run_directory": run_dir.resolve().name,
            "plasim_diag_sha256": hashlib.sha256(diag.read_bytes()).hexdigest(),
            **summarise(conv, fixer),
        }

    provenance = {
        "generator": "exoplasim/analysis/arms/read_conversion_decomposition.py",
        "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                                     text=True).stdout.strip(),
        "python": sys.version.split()[0],
        "prints_dropped": DROP_PRINTS,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"provenance": provenance, "arms": arms}, indent=2) + "\n",
                        encoding="utf-8")
    print(json.dumps(arms, indent=2)[:4000])


if __name__ == "__main__":
    main()
