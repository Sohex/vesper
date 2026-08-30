#!/usr/bin/env python3
"""Build the exact stock LPJ-GUESS 4.1.1 ntransform arm.

The comparison changes one compilation unit in a copy of the current Vesper
tree.  Everything else, including the CENTURY accelerator and Vesper input
module, therefore remains byte-identical to the active arm.  The stock source
is read from the original parent of the vendored subtree merge and held to its
known digest before it can enter a build.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from _paths import GENERATED, GUESS_SOURCE, PROJECT_ROOT

STOCK_OBJECT = "b489c42f4^2"
STOCK_PATH = "modules/ntransform.cpp"
STOCK_SHA256 = "77ecbdb5a68f447ef72685be6de082a765a53bf0765178c0dea36a1b7a0a0d39"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def source_fingerprint(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if not path.is_file() or "build" in relative.parts:
            continue
        digest.update(str(relative).encode("utf-8") + b"\0")
        digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> None:
    stock = subprocess.check_output(
        ["git", "show", f"{STOCK_OBJECT}:{STOCK_PATH}"], cwd=PROJECT_ROOT)
    if sha256_bytes(stock) != STOCK_SHA256:
        raise SystemExit("the recorded stock LPJ-GUESS 4.1.1 source has drifted")

    live_sha = sha256(GUESS_SOURCE / STOCK_PATH)
    base_fingerprint = source_fingerprint(GUESS_SOURCE)
    arm_root = GENERATED / f"ntransform_stock_4_1_1_{base_fingerprint[:12]}"
    source = arm_root / "source"
    build = arm_root / "build"
    binary = build / "guess"
    provenance_path = binary.with_suffix(".provenance.json")

    if binary.is_file() and provenance_path.is_file():
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        if (provenance.get("base_source_fingerprint") == base_fingerprint and
                provenance.get("binary_sha256") == sha256(binary)):
            print(binary)
            return
        raise SystemExit(f"existing arm at {arm_root} does not match its provenance")
    if arm_root.exists():
        raise SystemExit(f"incomplete stock arm already exists at {arm_root}")

    arm_root.mkdir(parents=True)
    shutil.copytree(GUESS_SOURCE, source, ignore=shutil.ignore_patterns("build"))
    (source / STOCK_PATH).write_bytes(stock)

    configure = ["cmake", "-S", str(source), "-B", str(build),
                 "-DCMAKE_BUILD_TYPE=Release", "-DUNIT_TESTS=OFF"]
    compile_command = ["cmake", "--build", str(build), "--parallel", "16"]
    subprocess.run(configure, check=True)
    subprocess.run(compile_command, check=True)

    provenance = {
        "contract_version": "vesper-ntransform-comparison-binary/1",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "profile": "stock-4.1.1",
        "release": "LPJ-GUESS 4.1.1, Zenodo 8065737, SVN r10118",
        "stock_git_object": STOCK_OBJECT,
        "substituted_path": STOCK_PATH,
        "stock_source_sha256": STOCK_SHA256,
        "live_source_sha256": live_sha,
        "base_source_fingerprint": base_fingerprint,
        "binary_sha256": sha256(binary),
        "invariant": (
            "copy of the active vendor/lpj-guess tree with only "
            "modules/ntransform.cpp replaced by the stock release source"
        ),
        "configure_command": configure,
        "compile_command": compile_command,
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n",
                               encoding="utf-8")
    print(binary)


if __name__ == "__main__":
    main()
