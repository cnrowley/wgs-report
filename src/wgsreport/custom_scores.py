"""Build polygenic scores locally from GWAS summary statistics that cannot be redistributed as weights.

Currently implemented: PGC3 schizophrenia (Trubetskoy et al. 2022, Nature, doi:10.1038/s41586-022-04434-5),
using the same clumping-and-thresholding (P+T) recipe as the paper: MAF > 1%, LD clumping at r^2 > 0.1 within
500 kb (1000 Genomes EUR as the LD reference), P < 0.05 (the threshold giving the paper's median liability R^2).

The result is written in the same table layout as a parsed PGS Catalog score (chrom,pos,effect,other,weight; GRCh38)
under <work>/custom_scores/, referenced from traits.yaml as `custom:<name>`.
"""
from __future__ import annotations

import gzip
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from .refpanel import DEFAULT_PATTERN, SUPERPOPS, load_panel_meta  # noqa: F401
from .utils import log, require_tools, run

SUMSTATS_URL = "https://ndownloader.figshare.com/files/34517828"   # PGC3 SCZ European autosomal (figshare 19426775, CC BY 4.0)
README_URL = "https://ndownloader.figshare.com/files/34865091"


# Column maps for supported summary-statistics files. maf_from: how to get a minor-allele frequency.
DATASETS = {
    "pgc3_scz": dict(rsid="ID", a1="A1", a2="A2", beta="BETA", p="PVAL", sep="\t", comment="#"),
}


def read_sumstats(name: str, path: Path, p_max: float, min_info: float = 0.6, min_maf: float = 0.01) -> pd.DataFrame:
    m = DATASETS[name]
    df = pd.read_csv(path, sep=m["sep"], comment=m.get("comment"), low_memory=False)
    d = pd.DataFrame({"rsid": df[m["rsid"]].astype(str), "A1": df[m["a1"]].astype(str).str.upper(),
                      "A2": df[m["a2"]].astype(str).str.upper(), "BETA": pd.to_numeric(df[m["beta"]], errors="coerce"),
                      "PVAL": pd.to_numeric(df[m["p"]], errors="coerce")})
    if name == "pgc3_scz":
        freq = (df["FCAS"] * df["NCAS"] + df["FCON"] * df["NCON"]) / (df["NCAS"] + df["NCON"])
        d["maf_ok"] = np.minimum(freq, 1 - freq) > min_maf
        d["q_ok"] = df["IMPINFO"] >= 0.8
    else:
        f = pd.to_numeric(df[m["freq"]], errors="coerce")
        n = pd.to_numeric(df[m["n"]], errors="coerce")
        d["maf_ok"] = np.minimum(f, 1 - f) > min_maf
        d["q_ok"] = n >= 0.5 * n.max()
        if "info" in m:
            d["q_ok"] &= pd.to_numeric(df[m["info"]], errors="coerce") >= min_info
    d = d[d["rsid"].str.startswith("rs") & d["A1"].str.fullmatch("[ACGT]") & d["A2"].str.fullmatch("[ACGT]")]
    d = d[d["maf_ok"] & d["q_ok"] & (d["PVAL"] < p_max) & d["BETA"].notna()]
    return d.drop(columns=["maf_ok", "q_ok"]).drop_duplicates("rsid").reset_index(drop=True)


def _extract_eur(chrom: str, ids_file: Path, kg_dir: Path, eur_keep: Path, out_dir: Path, threads: int) -> Path:
    prefix = out_dir / f"eur.chr{chrom}"
    if prefix.with_suffix(".pgen").exists():
        return prefix
    vcf = Path(kg_dir) / DEFAULT_PATTERN.format(chrom=chrom)
    tmp = out_dir / f"tmp.eur.chr{chrom}.vcf.gz"
    run(f"bcftools view -S {eur_keep} -i 'ID=@{ids_file}' -m2 -M2 -v snps {vcf} -Ou "
        f"| bcftools annotate -x INFO,^FMT/GT -Oz -o {tmp}")
    run(["plink2", "--vcf", str(tmp), "--vcf-half-call", "missing", "--rm-dup", "force-first",
         "--maf", "0.01", "--make-pgen", "--threads", str(threads), "--out", str(prefix)])
    tmp.unlink()
    return prefix


def build_pt_score(name: str, sumstats: Path, kg_dir: Path, work_dir: Path, p_max: float = 0.05,
                   clump_r2: float = 0.1, clump_kb: int = 500, jobs: int = 8, threads: int = 4) -> Path:
    require_tools("bcftools", "plink2")
    work = Path(work_dir)
    out_dir = work / "custom_scores" / f"{name}_build"
    out_dir.mkdir(parents=True, exist_ok=True)
    final = work / "custom_scores" / f"{name}_pt.tsv"
    if final.exists():
        log(f"already built: {final}")
        return final

    ss = read_sumstats(name, sumstats, p_max)
    log(f"{name}: {len(ss):,} SNPs with P<{p_max}, MAF>1%, passing quality filters")
    ids_file = out_dir / "ids.txt"
    ids_file.write_text("\n".join(ss["rsid"]) + "\n")
    meta = load_panel_meta(kg_dir)
    eur_keep = out_dir / "eur.keep"
    eur_keep.write_text("\n".join(meta.loc[meta.super_pop == "EUR", "IID"]) + "\n")

    log("extracting 1000G EUR genotypes at those rsIDs (22 chromosomes)")
    with ThreadPoolExecutor(jobs) as ex:
        prefixes = list(ex.map(lambda c: _extract_eur(str(c), ids_file, kg_dir, eur_keep, out_dir, threads), range(1, 23)))
    lst = out_dir / "merge.list"
    lst.write_text("\n".join(str(p) for p in prefixes) + "\n")
    merged = out_dir / "eur"
    run(["plink2", "--pmerge-list", str(lst), "--make-pgen", "--out", str(merged)])

    pv = pd.read_csv(f"{merged}.pvar", sep="\t", comment="#", header=None, usecols=range(5),
                     names=["chrom", "pos", "rsid", "ref", "alt"], dtype={"chrom": str})
    ss = ss[ss["rsid"].isin(pv["rsid"])]
    ss[["rsid", "PVAL"]].rename(columns={"rsid": "SNP", "PVAL": "P"}).to_csv(out_dir / "clump_in.tsv", sep="\t", index=False)
    log(f"{len(ss):,} SNPs present in the 1000G EUR panel; clumping (r2>{clump_r2}, {clump_kb} kb)")
    run(["plink2", "--pfile", str(merged), "--clump", str(out_dir / "clump_in.tsv"), "--clump-id-field", "SNP",
         "--clump-p-field", "P", "--clump-p1", str(p_max), "--clump-p2", str(p_max), "--clump-r2", str(clump_r2),
         "--clump-kb", str(clump_kb), "--threads", str(threads), "--out", str(out_dir / "clumped")])
    cl = pd.read_csv(out_dir / "clumped.clumps", sep=r"\s+")
    idx = cl["ID"] if "ID" in cl.columns else cl["SNP"]
    keep = ss[ss["rsid"].isin(set(idx))].merge(pv, on="rsid")
    log(f"{len(keep):,} independent index SNPs retained")
    # effect allele A1 (log-OR BETA) -> GRCh38 coordinates from the 1000G panel
    ok = ((keep["A1"] == keep["alt"]) & (keep["A2"] == keep["ref"])) | ((keep["A1"] == keep["ref"]) & (keep["A2"] == keep["alt"]))
    n_bad = int((~ok).sum())
    keep = keep[ok]
    out = pd.DataFrame({"chrom": keep["chrom"], "pos": keep["pos"], "effect": keep["A1"], "other": keep["A2"],
                        "weight": keep["BETA"]})
    final.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(final, sep="\t", index=False)
    log(f"wrote {len(out):,} weights to {final} ({n_bad} dropped for allele mismatch)")
    for p in prefixes:
        for ext in (".pgen", ".pvar", ".psam", ".log"):
            Path(str(p) + ext).unlink(missing_ok=True)
    return final


def build_pgc3_scz(sumstats, kg_dir, work_dir, **kw):
    return build_pt_score("pgc3_scz", sumstats, kg_dir, work_dir, p_max=0.05, **kw)
