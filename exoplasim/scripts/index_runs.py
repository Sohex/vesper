#!/usr/bin/env python3
"""Index every run directory by what it is, since its name says nothing.

    python exoplasim/scripts/index_runs.py           # print and write INDEX.json
    python exoplasim/scripts/index_runs.py --find flux_ratio=0.968

Run ids are UUIDs. That is deliberate -- see `run_id` in `run_exoplasim.py` --
and it moves the entire burden of "what was this run" onto the manifest each run
already writes. This reads those manifests and produces
`exoplasim/runs/INDEX.json`, which is the human- and machine-readable answer.

Runs written before the UUID change carry their parameters in their directory
name instead, and most carry no `physical` block. Those are indexed from their
manifest where possible and flagged `legacy` where not, rather than being parsed
back out of the name: a name is not a record.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
import sys as _sys
if str(ROOT / "lib") not in _sys.path:
    _sys.path.insert(0, str(ROOT / "lib"))
from paths import rel  # noqa: E402
RUNS = ROOT / "exoplasim" / "runs"
INDEX = RUNS / "INDEX.json"

FIELDS = ["resolution", "layers", "ranks", "flux_ratio", "co2_ppm",
          "rotation_hours", "obliquity_degrees", "eccentricity", "glaciers",
          "stellar_spectrum", "geography"]


def read(run_dir: Path) -> dict | None:
    mf = run_dir / "run_manifest.json"
    if not mf.is_file():
        return None
    try:
        m = json.loads(mf.read_text(encoding="utf-8"))
    except Exception:
        return None

    phys = m.get("physical")
    legacy = phys is None
    if legacy:
        # Recover what the manifest does hold. source_config is a full copy of
        # the configuration the run was prepared with, so most of it survives
        # even for runs that predate the fingerprint block.
        cfg = m.get("source_config") or {}
        p, a, mo = (cfg.get("planet", {}), cfg.get("atmosphere", {}),
                    cfg.get("model", {}))
        phys = {
            "resolution": mo.get("resolution"),
            "layers": mo.get("layers"),
            "ranks": mo.get("ncpus"),
            # The run's ACTUAL flux, from derived_parameters -- not
            # orbit.baseline_flux_earth, which is the config's baseline at
            # prepare time and is simply a different number whenever a run was
            # launched with --flux-ratio. Reading the baseline here reported the
            # 0.968 bootstrap as 0.945.
            "flux_ratio": (m.get("derived_parameters", {}) or {})
                          .get("stellar_flux_ratio_earth"),
            "co2_ppm": (round(1e6 * float(a["pCO2_bar"]), 3)
                        if a.get("pCO2_bar") is not None else None),
            "rotation_hours": p.get("rotation_hours"),
            "obliquity_degrees": p.get("obliquity_degrees"),
            "eccentricity": p.get("eccentricity"),
            "glaciers": ((cfg.get("surface", {}) or {})
                         .get("glaciers", {}) or {}).get("enabled"),
            "stellar_spectrum": (cfg.get("radiation", {}) or {})
                                .get("stellar_spectrum"),
            "geography": None,
        }

    orbits = len(list(run_dir.glob("MOST.*.nc")))
    conv = m.get("convergence_assessment") or {}
    return {
        "run_id": m.get("run_id", run_dir.name),
        "directory": run_dir.name,
        "status": m.get("status"),
        "legacy_name": legacy,
        "orbits_on_disk": orbits,
        "snapshots_on_disk": len(list((run_dir / "snapshots").glob("*.nc")))
                             if (run_dir / "snapshots").is_dir() else 0,
        "size_gb": round(sum(f.stat().st_size for f in run_dir.rglob("*")
                             if f.is_file()) / 1e9, 2),
        "executable_sha256": ((m.get("executable") or {}).get("sha256")
                              or m.get("executable_sha256")),
        # WHICH ARM, if the run was one. Run ids are UUIDs, so this index is the
        # only record of what each run WAS, and an arm run -- another precision,
        # another flag line, a patched model source -- is not the shipped model
        # and must not read as it. Null for every ordinary run. world-u5pf.
        "arm": ((m.get("executable") or {}).get("arm") or {}).get("build_tag"),
        "config_sha256": m.get("config_sha256"),
        "physical": phys,
        "source_build": (m.get("source_config") or {}).get("source_build"),
        # Every climatology this run produced, keyed by label, so a consumer can
        # resolve a climatology back to the run and build that made it WITHOUT
        # globbing the analysis directory and taking whatever sorts last.
        "climatologies": {
            label: {"regular": Path(e.get("regular") or "").name,
                    "snapshots": Path(e.get("snapshots") or "").name,
                    "orbit_count": e.get("orbit_count")}
            for label, e in (m.get("climatologies") or {}).items()},
        "converged": conv.get("sufficiently_equilibrated_for_worldbuilding",
                              conv.get("pass")),
        "convergence_metrics": conv.get("metrics") or {},
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--find", action="append", default=[],
                    help="filter as key=value, repeatable, e.g. resolution=T42")
    ap.add_argument("--output", type=Path, default=INDEX)
    args = ap.parse_args()

    rows = [r for r in (read(d) for d in sorted(RUNS.iterdir()) if d.is_dir())
            if r is not None]

    shown = rows
    for f in args.find:
        k, _, v = f.partition("=")
        shown = [r for r in shown if str(r["physical"].get(k)).lower() == v.lower()]

    print(f"{'directory':66}{'res':5}{'flux':>7}{'orb':>5}{'GB':>7}  conv")
    for r in shown:
        p = r["physical"]
        flux = p.get("flux_ratio")
        print(f"{r['directory'][:64]:66}{str(p.get('resolution') or '-'):5}"
              f"{(f'{flux:.3f}' if isinstance(flux,(int,float)) else '-'):>7}"
              f"{r['orbits_on_disk']:>5}{r['size_gb']:>7.2f}  "
              f"{'yes' if r['converged'] else ''}")
    print(f"\n{len(shown)} of {len(rows)} runs, "
          f"{sum(r['size_gb'] for r in rows):.1f} GB total")

    payload = {
        "note": "What each run directory is. Run ids are UUIDs and carry no "
                "meaning; this is where the physics lives. Generated by "
                "exoplasim/scripts/index_runs.py; do not edit.",
        "generated": datetime.now(timezone.utc).isoformat(),
        "runs": rows,
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {rel(args.output)}")


if __name__ == "__main__":
    main()
