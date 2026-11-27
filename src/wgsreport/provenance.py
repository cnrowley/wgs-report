"""Execution manifest: what software and reference data produced a report, so a result can be reproduced or compared.

Nothing personal goes in: tool versions, score identifiers with file checksums, reference-data releases. No paths, sample names or genotypes.
"""
from __future__ import annotations

import hashlib
import platform
import re
import subprocess
from datetime import date
from pathlib import Path

from . import __version__

TOOLS = {"bcftools": (["bcftools", "--version"], r"bcftools (\S+)"), "samtools": (["samtools", "--version"], r"samtools (\S+)"),
         "plink2": (["plink2", "--version"], r"PLINK v(\S+ \S+)")}
PACKAGES = ("numpy", "pandas", "scipy", "matplotlib", "yaml")


def _tool_version(cmd: list[str], rx: str) -> str | None:
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, errors="replace", timeout=30)
    except (OSError, subprocess.TimeoutExpired):
        return None
    m = re.search(rx, out.stdout + out.stderr)
    return m.group(1) if m else None


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def collect(score_ids: list[str] | None = None, pgs_cache: Path | None = None, work: Path | None = None, extra: dict | None = None) -> dict:
    import importlib
    pkgs = {}
    for p in PACKAGES:
        try:
            pkgs[p] = getattr(importlib.import_module(p), "__version__", "?")
        except ImportError:
            pkgs[p] = None
    scores = {}
    for sid in score_ids or []:
        f = Path(pgs_cache) / f"{sid}_hmPOS_GRCh38.txt.gz" if pgs_cache else None
        scores[sid] = _sha256(f)[:16] if f and f.exists() else None
    ref = {}
    if work and (Path(work) / "panel" / "samples.tsv").exists():
        ref["reference_panel_samples"] = sum(1 for _ in open(Path(work) / "panel" / "samples.tsv")) - 1
    return {"wgsreport": __version__, "date": date.today().isoformat(), "python": platform.python_version(),
            "platform": f"{platform.system()} {platform.machine()}", "packages": pkgs,
            "tools": {n: _tool_version(*v) for n, v in TOOLS.items()}, "pgs_score_sha256_16": scores, **ref, **(extra or {})}
