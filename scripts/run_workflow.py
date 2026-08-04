#!/usr/bin/env python3
"""Canonical dispatcher for the packaged AA analysis entry points.

This dispatcher provides one discoverable entry to list, verify or invoke a
single workflow.  External and licensed inputs are never downloaded
automatically.  Workflow-specific arguments are forwarded unchanged after
``--``. Use a new output directory for every run; output handling is governed
by the selected workflow.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ORDER_FILE = ROOT / "WORKFLOW_ORDER.tsv"


def load_steps() -> list[dict[str, str]]:
    with ORDER_FILE.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def command_for(entry: Path, forwarded: list[str]) -> list[str]:
    suffix = entry.suffix.lower()
    if suffix == ".py":
        return [sys.executable, str(entry), *forwarded]
    if suffix == ".r":
        executable = shutil.which("Rscript")
        if executable is None:
            raise RuntimeError("Rscript is not available on PATH")
        return [executable, str(entry), *forwarded]
    if suffix in {".mjs", ".js"}:
        executable = shutil.which("node")
        if executable is None:
            raise RuntimeError("node is not available on PATH")
        return [executable, str(entry), *forwarded]
    raise RuntimeError(f"Unsupported entry-point type: {entry}")


def main() -> None:
    parser = argparse.ArgumentParser()
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--list", action="store_true")
    action.add_argument("--check", action="store_true")
    action.add_argument("--run", metavar="WORKFLOW")
    parser.add_argument("forwarded", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    steps = load_steps()
    by_name = {step["workflow"]: step for step in steps}

    if args.list:
        for step in steps:
            print(
                f'{step["order"]}\t{step["workflow"]}\t'
                f'{step["canonical_entry"]}\t{step["current_scope"]}'
            )
        return

    if args.check:
        missing = []
        for step in steps:
            entry = ROOT / step["canonical_entry"]
            state = "present" if entry.is_file() else "missing"
            print(f'{state}\t{step["workflow"]}\t{entry}')
            if state == "missing":
                missing.append(entry)
        if missing:
            raise SystemExit(1)
        return

    step = by_name.get(args.run)
    if step is None:
        choices = ", ".join(by_name)
        raise SystemExit(f"Unknown workflow {args.run!r}. Choices: {choices}")
    entry = ROOT / step["canonical_entry"]
    if not entry.is_file():
        raise FileNotFoundError(entry)
    forwarded = args.forwarded
    if forwarded and forwarded[0] == "--":
        forwarded = forwarded[1:]
    subprocess.run(command_for(entry, forwarded), check=True)


if __name__ == "__main__":
    main()
