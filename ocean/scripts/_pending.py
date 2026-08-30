"""Shared fail-closed entry point for registered ocean work not yet accepted."""

from __future__ import annotations

import argparse


def pending_main(issue: str, purpose: str) -> None:
    parser = argparse.ArgumentParser(description=purpose)
    parser.add_argument(
        "--explain",
        action="store_true",
        help="print the blocking contract and exit without producing an artifact",
    )
    parser.parse_args()
    parser.exit(
        2,
        f"REFUSED: {issue} is not accepted; this registered step cannot yet "
        "produce a valid artifact.\n",
    )
