"""Screen the user's VCF for ClinVar pathogenic / likely pathogenic variants."""
from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from .registry import data_path
from .utils import download, log, run

from .sources import CLINVAR_URL as URL
STARS = {"practice_guideline": 4, "reviewed_by_expert_panel": 3,
         "criteria_provided,_multiple_submitters,_no_conflicts": 2,
         "criteria_provided,_single_submitter": 1}


def acmg_genes() -> set[str]:
    """Gene symbols on the ACMG secondary-findings list (data/acmg_sf_v3.3.tsv, with provenance in its header)."""
    rows = data_path("acmg_sf_v3.3.tsv").read_text().splitlines()
    return {l.split("\t")[0] for l in rows if l and not l.startswith("#") and not l.startswith("gene\t")}


def carrier_genes() -> dict[str, str]:
    """ACMG tier-3 carrier-screening genes (Gregg 2021, Tables 1-6): gene -> 'condition (AR|XL)'."""
    out = {}
    for l in data_path("acmg_carrier_tier3.tsv").read_text().splitlines():
        if l and not l.startswith("#") and not l.startswith("gene\t"):
            g, cond, inh, _ = l.split("\t")
            out[g] = f"{cond} ({inh})" if cond else inh
    return out


# Genes added by hand because the literature supports them but ClinGen lists no Definitive/Strong cancer curation under that name.
EXTRA_CANCER_GENES = {"HOXB13": "AD (prostate cancer; G84E founder variant, added from the literature, not from the ClinGen list)"}


def cancer_genes() -> dict[str, str]:
    """ClinGen Definitive/Strong cancer-predisposition genes (data/cancer_predisposition_clingen.tsv): gene -> inheritance modes."""
    out = {}
    for l in data_path("cancer_predisposition_clingen.tsv").read_text().splitlines():
        if l and not l.startswith("#") and not l.startswith("gene\t"):
            g, moi, _ = l.split("\t")
            out[g] = moi
    return {**EXTRA_CANCER_GENES, **out}


def fetch(dest: Path) -> Path:
    download(URL, dest)
    return dest


def load_plp(clinvar_vcf: Path, min_stars: int = 2) -> pd.DataFrame:
    """P/LP ClinVar records (one row per ALT) with at least `min_stars` review stars."""
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\t%INFO/GENEINFO\t%INFO/MC\t%INFO/AF_TGP\t%INFO/RS\n"
    out = run(["bcftools", "query", "-f", fmt, str(clinvar_vcf)]).stdout
    rows = []
    for line in out.splitlines():
        c, pos, ref, alt, vid, sig, rev, dn, gi, mc, af, rs = line.split("\t")
        if "athogenic" not in sig or "Conflicting" in sig or sig.startswith("Benign"):
            continue
        if not (sig.startswith("Pathogenic") or sig.startswith("Likely_pathogenic")):
            continue
        stars = STARS.get(rev.replace("|", ","), 0)
        if stars < min_stars:
            continue
        rows.append((c, int(pos), ref, alt, vid, sig, stars, dn.replace("_", " "),
                     gi.split(":")[0].split("|")[0], mc.split("|")[-1].split(",")[-1] if mc != "." else "",
                     float(af) if af not in (".", "") else None, f"rs{rs}" if rs != "." else ""))
    return pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "clinvar_id", "significance", "stars",
                                       "condition", "gene", "consequence", "af_1000g", "rsid"])


def screen(plp: pd.DataFrame, sites: dict) -> pd.DataFrame:
    """Match the user's genotypes (sites) against the P/LP table."""
    key = {(r.chrom, r.pos, r.ref, r.alt): r for r in plp.itertuples()}
    acmg = acmg_genes()
    carrier = carrier_genes()
    cancer = cancer_genes()
    hits = []
    for (chrom, pos), recs in sites.items():
        for ref, alts, gt, ok in recs:
            if not ok or any(g is None for g in gt):
                continue
            for j, alt in enumerate(alts, start=1):
                n = sum(1 for g in gt if g == j)
                if n and (chrom, pos, ref, alt) in key:
                    r = key[(chrom, pos, ref, alt)]
                    hits.append({**r._asdict(), "zygosity": "homozygous" if n == 2 or len(gt) == 1 and n == 1 and chrom == "X" else "heterozygous",
                                 "acmg_sf_gene": r.gene in acmg, "carrier_panel_gene": r.gene in carrier,
                                 "cancer_gene": cancer.get(r.gene, "")})
    df = pd.DataFrame(hits).drop(columns=["Index"], errors="ignore")
    return add_pattern(df)


def add_pattern(df: pd.DataFrame) -> pd.DataFrame:
    """Label each hit: homozygous, possible compound heterozygote (>=2 distinct heterozygous P/LP variants in one gene,
    phase unknown without parents or long reads) or single heterozygous (carrier status if the condition is recessive,
    possible disease risk if dominant: inheritance is not assessed here)."""
    if df.empty:
        return df
    df = df.copy()
    het_n = df[df.zygosity == "heterozygous"].groupby("gene")["clinvar_id"].nunique()
    def lab(r):
        if r.zygosity == "homozygous":
            return "homozygous"
        return "possible compound heterozygote (phase unknown)" if het_n.get(r.gene, 0) >= 2 else "single heterozygous"
    df["pattern"] = df.apply(lab, axis=1)
    return df
