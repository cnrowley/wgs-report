"""Genotype the curated SNP panel and decode APOE."""
from __future__ import annotations

import pandas as pd

from .registry import data_path


def load_panel(include_sensitive: bool) -> pd.DataFrame:
    df = pd.read_csv(data_path("snp_panel.tsv"), sep="\t", dtype={"chrom": str})
    if not include_sensitive:
        df = df[df["category"] != "sensitive"]
    return df.reset_index(drop=True)


def panel_positions(panel: pd.DataFrame) -> set[tuple[str, int]]:
    return set(zip(panel["chrom"], panel["pos"].astype(int)))


def genotype(panel: pd.DataFrame, sites: dict, cram_calls: dict | None = None) -> pd.DataFrame:
    """Genotype the panel from the VCF; where read-level calls exist (see cram.py) they confirm or replace assumptions."""
    rows = []
    for r in panel.itertuples():
        recs = [x for x in sites.get((r.chrom, int(r.pos)), []) if len(x[0]) == 1]
        status, alleles = "called", None
        if not recs:
            status, alleles = "assumed hom-ref", [r.ref, r.ref]
        else:
            ref, alts, gt, ok = recs[0]
            seq = [ref] + alts
            if not ok or any(g is None for g in gt):
                status = "filtered/no-call"
            else:
                alleles = [seq[g] for g in gt]
                if len(alleles) == 1:        # haploid call (male X): treat as hemizygous
                    alleles = alleles * 2
        cc = (cram_calls or {}).get(r.rsid)
        if cc and cc["alleles"]:
            if status == "assumed hom-ref":
                if sorted(cc["alleles"]) == [r.ref, r.ref]:
                    status = f"confirmed from reads ({cc['depth']}x)"
                else:
                    status, alleles = f"reads disagree with VCF ({cc['note']})", cc["alleles"]
            elif alleles is not None and sorted(alleles) != sorted(cc["alleles"]):
                status = f"VCF/reads discordant (reads: {'/'.join(cc['alleles'])})"
            elif alleles is not None:
                status = f"called; confirmed from reads ({cc['depth']}x)"
            elif alleles is None:
                status, alleles = f"from reads ({cc['note']})", cc["alleles"]
        if alleles is None:
            n_eff, gt_str, note = None, "./.", "Genotype not reliably called"
        else:
            n_eff = sum(a == r.effect for a in alleles)
            gt_str = "/".join(sorted(alleles))
            note = [r.note0, r.note1, r.note2][n_eff]
        rows.append({"rsid": r.rsid, "gene": r.gene, "category": r.category, "trait": r.trait,
                     "chrom": r.chrom, "pos": int(r.pos), "genotype": gt_str, "effect_allele": r.effect,
                     "n_effect": n_eff, "interpretation": note, "status": status})
    out = pd.DataFrame(rows)
    ev = pd.read_csv(data_path("snp_evidence.tsv"), sep="\t", comment="#").set_index("rsid")["level"]
    out["evidence"] = out["rsid"].map(ev).fillna("C")
    return out


def decode_apoe(geno: pd.DataFrame) -> str | None:
    g = geno.set_index("rsid")
    if "rs429358" not in g.index or "rs7412" not in g.index:
        return None
    c4, t2 = g.loc["rs429358", "n_effect"], g.loc["rs7412", "n_effect"]
    if pd.isna(c4) or pd.isna(t2):
        return "undetermined"
    c4, t2 = int(c4), int(t2)
    table = {(0, 0): "e3/e3", (0, 1): "e2/e3", (0, 2): "e2/e2", (1, 0): "e3/e4",
             (1, 1): "e2/e4 (or e1/e3; e2/e4 assumed)", (2, 0): "e4/e4"}
    return table.get((c4, t2), "unusual combination")
