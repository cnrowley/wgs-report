"""Harmonize PGS weights to the reference panel, score the panel and score the user.

Coding convention used everywhere in this module
-------------------------------------------------
Every panel variant is a biallelic SNV with REF/ALT alleles. For each polygenic score we
store a weight on the ALT allele,

    w_alt = +w   if the score's effect allele is ALT
    w_alt = -w   if the score's effect allele is REF   (the constant 2*w is dropped; it shifts
                                                       every individual equally and so cancels
                                                       in z-scores)

so that  score = sum_v w_alt[v] * dosage_alt[v]  for panel individuals and for the user alike.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from .utils import COMP, log, run

PALINDROMIC_MAF_CUTOFF = 0.42  # strand of A/T and C/G SNPs is only trusted below this MAF


def harmonize(score: pd.DataFrame, pvar: pd.DataFrame, af: pd.Series) -> tuple[pd.Series, dict]:
    """Map one score to panel ALT-allele weights.

    score: columns chrom,pos,effect,other,weight (see pgs.read_score)
    pvar : panel variants (chrom,pos,id,ref,alt)
    af   : panel ALT allele frequency indexed by variant id (used only for palindromic SNPs)
    Returns (weights indexed by variant id, stats dict).
    """
    n_in = len(score)
    m = score.merge(pvar, on=["chrom", "pos"], how="inner")
    stats = {"n_variants": n_in, "n_position_match": len(m)}
    e, o, r, a = m["effect"], m["other"], m["ref"], m["alt"]
    other_ok_alt = o.isna() | (o == r)       # 'other' consistent with effect==alt
    other_ok_ref = o.isna() | (o == a)
    direct_alt = (e == a) & other_ok_alt
    direct_ref = (e == r) & other_ok_ref
    ce = e.map(COMP)
    co = o.map(COMP)
    flip_alt = ~direct_alt & ~direct_ref & (ce == a) & (co.isna() | (co == r))
    flip_ref = ~direct_alt & ~direct_ref & (ce == r) & (co.isna() | (co == a))
    sign = np.where(direct_alt | flip_alt, 1.0, np.where(direct_ref | flip_ref, -1.0, np.nan))
    m = m.assign(w_alt=m["weight"].to_numpy() * sign, flipped=flip_alt | flip_ref)
    stats["n_strand_flipped"] = int(m["flipped"].sum())
    stats["n_allele_mismatch"] = int(np.isnan(sign).sum())
    m = m[~np.isnan(m["w_alt"])]
    # palindromic SNPs: drop when MAF is too close to 0.5 for the strand to be trusted
    pal = ((m["ref"].map(COMP) == m["alt"]))
    maf = af.reindex(m["id"]).to_numpy()
    maf = np.minimum(maf, 1 - maf)
    drop = pal.to_numpy() & ~(maf < PALINDROMIC_MAF_CUTOFF)  # NaN maf -> dropped as well
    stats["n_palindromic_dropped"] = int(drop.sum())
    m = m[~drop]
    w = m.groupby("id")["w_alt"].sum()   # duplicate rows for the same variant are summed
    stats["n_scored"] = int(len(w))
    stats["coverage"] = stats["n_scored"] / n_in if n_in else 0.0
    return w, stats


def write_weight_file(weights: pd.DataFrame, alt: pd.Series, path: Path) -> None:
    """plink2 --score input: ID, A1 (ALT), one column per score."""
    df = weights.copy()
    df.insert(0, "A1", alt.reindex(df.index).to_numpy())
    df.index.name = "ID"
    df.to_csv(path, sep="\t", float_format="%.8g")


def score_panel(pfile: Path, weights: pd.DataFrame, alt: pd.Series, workdir: Path,
                threads: int = 8) -> pd.DataFrame:
    """Score every panel sample with plink2. Returns DataFrame indexed by IID, one column per score."""
    workdir = Path(workdir)
    wfile = workdir / "weights.tsv"
    write_weight_file(weights, alt, wfile)
    n = weights.shape[1]
    out = workdir / "panel_scores"
    run(["plink2", "--pfile", str(pfile), "--score", str(wfile), "1", "2", "header-read",
         "cols=scoresums", "--score-col-nums", f"3-{2 + n}", "--threads", str(threads),
         "--out", str(out)])
    df = pd.read_csv(f"{out}.sscore", sep="\t")
    df = df.rename(columns={df.columns[0]: "IID"}).set_index("IID")
    df.columns = [c[:-4] if c.endswith("_SUM") else c for c in df.columns]
    return df[list(weights.columns)]


# ---------------------------------------------------------------------------------------
# User genotypes
# ---------------------------------------------------------------------------------------

def vcf_has_chr_prefix(vcf: Path) -> bool:
    out = run(["bcftools", "index", "-s", str(vcf)]).stdout if Path(str(vcf) + ".tbi").exists() \
        or Path(str(vcf) + ".csi").exists() else run(["bcftools", "view", "-h", str(vcf)]).stdout
    for line in out.splitlines():
        if line.startswith("##contig=<ID="):
            return line.split("ID=")[1].startswith("chr")
        if line and not line.startswith("#"):
            return line.split("\t")[0].startswith("chr")
    return False


def read_user_sites(vcf: Path, positions: set[tuple[str, int]]) -> dict:
    """Stream the user's VCF once and keep records starting at `positions` (non-'chr' naming).

    Returns {(chrom,pos): [(ref, [alts], gt_alleles(list[int|None]), filter_ok(bool)), ...]}.
    Several records can start at one position (e.g. an SNV and a deletion).
    """
    cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT]\n", str(vcf)]
    out: dict = {}
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20) as p:
        for line in p.stdout:  # type: ignore[union-attr]
            c, pos, ref, alt, flt, gt = line.rstrip("\n").split("\t")
            key = (c[3:] if c.startswith("chr") else c, int(pos))
            if key not in positions:
                continue
            gts = gt.replace("|", "/").split("/")
            alleles = [None if g in (".", "") else int(g) for g in gts]
            out.setdefault(key, []).append((ref, alt.split(","), alleles, flt in ("PASS", ".")))
        if p.wait() != 0:
            raise RuntimeError("bcftools query failed on user VCF")
    return out


def user_dosages(sites: dict, pvar: pd.DataFrame, pass_only: bool = True) -> tuple[np.ndarray, dict, np.ndarray]:
    """ALT-allele dosage (0/1/2, NaN = missing) for every panel variant, in pvar order.

    A panel variant absent from the user's VCF is assumed homozygous REF (the VCF only lists
    variants). Records failing filters, or with no-call genotypes, become NaN (mean-imputed later).
    """
    n = len(pvar)
    dos = np.zeros(n, dtype=np.float32)
    assumed = np.zeros(n, dtype=bool)
    stats = {"n_variants": n, "from_vcf_het": 0, "from_vcf_hom_alt": 0, "assumed_hom_ref": 0,
             "missing_filtered_or_nocall": 0, "third_allele": 0, "diploid_mismatch_ref": 0}
    chroms = pvar["chrom"].to_numpy()
    poss = pvar["pos"].to_numpy()
    refs = pvar["ref"].to_numpy()
    alts = pvar["alt"].to_numpy()
    for i in range(n):
        recs = [r for r in sites.get((chroms[i], int(poss[i])), []) if len(r[0]) == 1]
        rec = recs[0] if recs else None
        if rec is None:
            stats["assumed_hom_ref"] += 1
            assumed[i] = True
            continue
        uref, ualts, gt, ok = rec
        if pass_only and not ok:
            dos[i] = np.nan
            stats["missing_filtered_or_nocall"] += 1
            continue
        if any(g is None for g in gt):
            dos[i] = np.nan
            stats["missing_filtered_or_nocall"] += 1
            continue
        if uref != refs[i]:
            stats["diploid_mismatch_ref"] += 1
            dos[i] = np.nan
            continue
        seq = [uref] + ualts
        called = [seq[g] if g < len(seq) else "?" for g in gt]
        d = sum(1 for x in called if x == alts[i])
        if any(x not in (refs[i], alts[i]) for x in called):
            stats["third_allele"] += 1  # carries an allele absent from the panel; counts toward neither
        dos[i] = d
        if d == 1:
            stats["from_vcf_het"] += 1
        elif d == 2:
            stats["from_vcf_hom_alt"] += 1
    return dos, stats, assumed


def user_scores(dos: np.ndarray, weights: pd.DataFrame, mean_alt_freq: np.ndarray) -> pd.Series:
    """Weighted sum of ALT dosages; missing dosages are replaced by 2*panel ALT frequency."""
    d = np.where(np.isnan(dos), 2.0 * mean_alt_freq, dos).astype(np.float64)
    return pd.Series(d @ weights.to_numpy(dtype=np.float64), index=weights.columns)


def top_contributors(dos: np.ndarray, weights: pd.Series, pvar: pd.DataFrame,
                     mean_alt_freq: np.ndarray, n: int = 8) -> pd.DataFrame:
    """Variants moving the user's score furthest from the panel average (|w * (d - 2p)|)."""
    d = np.where(np.isnan(dos), 2.0 * mean_alt_freq, dos)
    contrib = weights.to_numpy() * (d - 2.0 * mean_alt_freq)
    idx = np.argsort(-np.abs(contrib))[:n]
    t = pvar.iloc[idx][["chrom", "pos", "ref", "alt"]].copy()
    t["genotype_alt_count"] = d[idx]
    t["weight_on_alt"] = weights.to_numpy()[idx]
    t["contribution"] = contrib[idx]
    t["rsid_hint"] = ""
    return t.reset_index(drop=True)
