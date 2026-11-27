#!/usr/bin/env python3
"""Scaffold a tool-backed module: `python scripts/new_module.py my_tool` creates src/wgsreport/my_tool.py and tests/test_my_tool.py.

The skeleton follows docs/EXTENDING.md: a thin run(), a parse(), a pure summarize(), and a report helper you wire in by hand.
"""
import re
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
if len(sys.argv) != 2 or not re.fullmatch(r"[a-z][a-z0-9_]*", sys.argv[1]):
    sys.exit("usage: python scripts/new_module.py <snake_case_name>")
name = sys.argv[1]
mod, test = root / "src" / "wgsreport" / f"{name}.py", root / "tests" / f"test_{name}.py"
if mod.exists() or test.exists():
    sys.exit(f"{mod.name} or {test.name} already exists")

mod.write_text(f'''"""TODO: one line on what {name} reports.

Source of truth for thresholds / lists: TODO (cite it and put the file in data/ with a provenance header).
Tool (not bundled): TODO name, license, where to download (add the URL to sources.py).
Limits: TODO what this cannot see and what needs confirmation.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from .utils import log


def run(tool_dir: Path, input_file: Path, out_dir: Path) -> Path:
    """Start the tool. Resolve paths first (the tool may run from another directory); delete personal intermediates in a finally block."""
    tool_dir, input_file, out_dir = Path(tool_dir).resolve(), Path(input_file).resolve(), Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    log("{name}: running")
    subprocess.run([str(tool_dir / "TODO-executable"), str(input_file), str(out_dir)], check=True)
    return out_dir / "TODO-output"


def parse(path: Path) -> list[dict]:
    """Read the tool's output into plain rows. Build this from one real output, then shrink it into a synthetic test fixture."""
    rows = []
    for line in Path(path).read_text().splitlines():
        f = line.split("\\t")
        if len(f) >= 2 and not line.startswith("#"):
            rows.append({{"name": f[0], "value": float(f[1])}})
    return rows


def summarize(rows: list[dict]) -> dict:
    """Pure function: categories and flags from parsed rows. Say 'not interpreted' when a row does not match the thresholds' definition,
    and report ties or ambiguity as such."""
    flagged = [r for r in rows if r["value"] > 0]       # TODO real thresholds from the cited data file
    return {{"rows": rows, "flagged": flagged}}
''')
test.write_text(f'''"""Tests for {name} on a synthetic fixture (never commit real tool output)."""
from wgsreport import {name}


def test_parse_and_summarize(tmp_path):
    f = tmp_path / "out.tsv"
    f.write_text("# comment\\nA\\t1.5\\nB\\t0\\n")
    rows = {name}.parse(f)
    assert [r["name"] for r in rows] == ["A", "B"]
    s = {name}.summarize(rows)
    assert [r["name"] for r in s["flagged"]] == ["A"]
''')
print(f"created {{mod.relative_to(root)}} and {{test.relative_to(root)}}\\nNext: wire cli.py, analyze.py, report.py (see docs/EXTENDING.md)")
