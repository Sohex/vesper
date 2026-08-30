#!/usr/bin/env python3
"""BIO-17 regression fixtures for rootable/lake/canopy compositing."""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT.parents[2]
sys.path.insert(0, str(SCRIPT.parent))

from build_surface_albedo import composite_rootable

REPORT = PROJECT_ROOT / "exoplasim/analysis/rootable_albedo_gate_report.json"


def main() -> None:
    checks: list[dict] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"name": name, "passed": bool(passed), "detail": detail})

    # Ordinary substrate, wholly barren ground, a half-water cell, and a pure
    # persistent lake. The LPJ cover values are deliberately nonzero over the
    # last three: only the rootable fraction is allowed to decide where they act.
    total = np.array([0.30, 0.50, 0.18, 0.06])
    root_contribution = np.array([0.30, 0.00, 0.15, 0.00])
    rootable = np.array([1.00, 0.00, 0.50, 0.00])
    tree_fpc = np.array([0.50, 1.00, 0.40, 1.00])
    grass_fpc = np.zeros(4)
    result, tree, grass, _ = composite_rootable(
        total, root_contribution, rootable, tree_fpc, grass_fpc, 0.10, 0.20)
    check("ordinary rootable cell",
          np.isclose(result[0], 0.20) and np.isclose(tree[0], 0.50),
          f"albedo {result[0]:.3f}, forest {tree[0]:.3f}")
    check("barren cell remains barren",
          result[1] == total[1] and tree[1] == 0.0 and grass[1] == 0.0,
          f"albedo {result[1]:.3f}, forest {tree[1]:.3f}")
    # Half rootable substrate at 0.30 and half water at 0.06. Forty percent of
    # rootable ground becomes 0.10 canopy: 0.03 water + 0.09 bare + 0.02 canopy.
    check("partial lake contribution preserved",
          np.isclose(result[2], 0.14) and np.isclose(tree[2], 0.20),
          f"albedo {result[2]:.3f}, forest {tree[2]:.3f}; water contributes 0.030")
    check("persistent lake remains exact water",
          result[3] == 0.06 and tree[3] == 0.0,
          f"albedo {result[3]:.3f}, forest {tree[3]:.3f}")

    # The same extensive-share algebra must commute with the model's two-band
    # recombination, or fixing broadband lakes would still leave their spectrum
    # repainted.
    z1, z2 = 0.4, 0.6
    total1, total2 = total * 0.7, total * 1.2
    root1, root2 = root_contribution * 0.7, root_contribution * 1.2
    tree1, tree2 = 0.10 * 0.5, 0.10 * (4.0 / 3.0)
    band1 = composite_rootable(
        total1, root1, rootable, tree_fpc, grass_fpc, tree1, 0.20)[0]
    band2 = composite_rootable(
        total2, root2, rootable, tree_fpc, grass_fpc, tree2, 0.20)[0]
    broadband_total = z1 * total1 + z2 * total2
    broadband_root = z1 * root1 + z2 * root2
    broadband = composite_rootable(
        broadband_total, broadband_root, rootable, tree_fpc, grass_fpc,
        z1 * tree1 + z2 * tree2, 0.20)[0]
    check("two-band recombination",
          np.allclose(z1 * band1 + z2 * band2, broadband),
          f"maximum residual {np.max(np.abs(z1 * band1 + z2 * band2 - broadband)):.3e}")

    source = (PROJECT_ROOT / "exoplasim/scripts/build_surface_albedo.py").read_text()
    check("BIO-11 partition is consumed",
          "read_rootable_partition(" in source,
          "modelled mode reads rootable, water and barren fractions")
    check("forest uses effective tree share",
          "forest = np.where(land_cells, np.clip(tree_cover" in source,
          "code 212 receives f_rootable times conditional tree FPC")

    failed = [item for item in checks if not item["passed"]]
    report = {
        "contract_version": "vesper-rootable-albedo/1",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "checks": checks,
        "summary": {"passed": len(checks) - len(failed), "failed": len(failed)},
        "verdict": "PASS" if not failed else "FAIL",
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    for item in checks:
        print(f"{'PASS' if item['passed'] else 'FAIL'}  {item['name']}: "
              f"{item['detail']}")
    print(f"\n{report['verdict']}: {len(checks) - len(failed)} passed, "
          f"{len(failed)} failed; {REPORT.relative_to(PROJECT_ROOT)}")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
