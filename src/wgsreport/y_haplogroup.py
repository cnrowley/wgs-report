"""Y-chromosome haplogroup with Yleaf (Ralf et al., 2018; https://github.com/genid/Yleaf, MIT). Not bundled: download the standalone
release, point its `config.txt` at your GRCh38 FASTA, and pass the directory that holds the `yleaf` executable.

The haplogroup is deterministic from known markers (it reads thousands of Y SNPs), which makes it one of the most reliable results in a
report. It describes the direct paternal line only, is not a health result, and can be absent in females.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from .utils import log


def run_yleaf(yleaf_exe: Path, cram: Path, fasta: Path, out_dir: Path, threads: int = 4) -> Path:
    exe = Path(yleaf_exe).resolve()
    if not exe.exists():
        raise SystemExit(f"{exe} not found. Download the Yleaf standalone release (sources.YLEAF_URL) and pass the `yleaf` executable.")
    out_dir = Path(out_dir).resolve()
    log("Y haplogroup: running Yleaf")
    # -force: Yleaf asks an interactive question when the output directory exists; --no-update-check keeps it offline.
    subprocess.run([str(exe), "-cram", str(Path(cram).resolve()), "-cr", str(Path(fasta).resolve()), "-rg", "hg38", "-o", str(out_dir),
                    "-t", str(threads), "--no-update-check", "-force"], check=True)
    return out_dir / "hg_prediction.hg"


def parse(hg_prediction: Path) -> dict:
    """First data row of Yleaf's hg_prediction.hg -> {haplogroup, marker, valid_markers, qc}."""
    lines = [l for l in Path(hg_prediction).read_text().splitlines() if l.strip()]
    hdr = lines[0].split("\t")
    row = dict(zip(hdr, lines[1].split("\t")))
    return {"haplogroup": row.get("Hg", ""), "marker": row.get("Hg_marker", ""), "valid_markers": int(float(row.get("Valid_markers", 0) or 0)),
            "qc": float(row.get("QC-score", 0) or 0)}
