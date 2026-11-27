"""PGS Catalog scoring-file download and parsing (GRCh38-harmonized files)."""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from .utils import download

from .sources import PGS_FTP as FTP
AUTOSOMES = [str(i) for i in range(1, 23)]


def _gzip_ok(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as f:
            while f.read(1 << 24):
                pass
        return True
    except (OSError, EOFError):
        return False


def fetch_score(pgs_id: str, cache_dir: Path) -> Path:
    """Download (once) and integrity-check a harmonized scoring file."""
    dest = Path(cache_dir) / f"{pgs_id}_hmPOS_GRCh38.txt.gz"
    for _ in range(3):
        download(FTP.format(id=pgs_id), dest)
        if _gzip_ok(dest):
            return dest
        dest.unlink(missing_ok=True)
    raise RuntimeError(f"could not obtain a valid scoring file for {pgs_id}")


def read_score(path: Path) -> tuple[pd.DataFrame, dict]:
    """Parse a harmonized PGS Catalog file.

    Returns (df, header) with df columns chrom, pos, effect, other, weight for
    autosomal single-nucleotide variants. Weights are converted to log scale
    when the file declares an OR/HR weight type.
    """
    header: dict[str, str] = {}
    with gzip.open(path, "rt") as f:
        skip = 0
        for line in f:
            if line.startswith("#"):
                skip += 1
                if "=" in line:
                    k, v = line.lstrip("#").strip().split("=", 1)
                    header[k] = v
            else:
                break
    df = pd.read_csv(path, sep="\t", skiprows=skip, comment=None, dtype=str, low_memory=False)
    need = {"hm_chr", "hm_pos", "effect_allele", "effect_weight"}
    if not need.issubset(df.columns):
        raise ValueError(f"{path.name}: missing harmonized columns {need - set(df.columns)}")
    other = df["other_allele"] if "other_allele" in df else pd.Series(np.nan, index=df.index)
    if "hm_inferOtherAllele" in df:
        other = other.fillna(df["hm_inferOtherAllele"])
    out = pd.DataFrame({
        "chrom": df["hm_chr"],
        "pos": pd.to_numeric(df["hm_pos"], errors="coerce"),
        "effect": df["effect_allele"].str.upper(),
        "other": other.str.upper() if other.dtype == object else other,
        "weight": pd.to_numeric(df["effect_weight"], errors="coerce"),
        "rsid": (df["hm_rsID"] if "hm_rsID" in df else df["rsID"] if "rsID" in df else pd.Series("", index=df.index)).fillna(""),
    })
    n_raw = len(out)
    out = out.dropna(subset=["chrom", "pos", "effect", "weight"])
    out = out[out["chrom"].isin(AUTOSOMES)]
    out["pos"] = out["pos"].astype(int)
    snv = out["effect"].str.fullmatch("[ACGT]")
    other_ok = out["other"].isna() | out["other"].astype(str).str.fullmatch("[ACGT]")
    out = out[snv & other_ok].reset_index(drop=True)
    wt = header.get("weight_type", "NR")
    if wt in ("OR", "HR"):
        out["weight"] = np.log(out["weight"])
    header["n_rows_raw"] = str(n_raw)
    header["n_rows_usable"] = str(len(out))
    return out, header
