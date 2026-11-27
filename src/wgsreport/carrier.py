"""Carrier status for family planning: ClinVar pathogenic / likely pathogenic variants in recessive and X-linked genes.

Which genes: ClinGen Definitive/Strong/Moderate AR or XL curations (data/gene_inheritance.tsv) plus the ACMG tier-3 carrier list
(Gregg 2021). Variants are NOT filtered by population frequency: common carrier alleles matter most for family planning.

Only results are reported. One sentence (LIMIT) states what a "nothing found" does not cover; do not add tables of untested genes.
"""
from __future__ import annotations

from .clinvar import carrier_genes
from .registry import data_path

LIMIT = ("Known ClinVar SNVs and short indels only. This is not a laboratory carrier screen: finding nothing does not rule out being a "
         "carrier, and it misses deletions and duplications, including the usual causes of spinal muscular atrophy and alpha-thalassemia.")


def inheritance_table() -> dict[str, dict]:
    lines = [x for x in data_path("gene_inheritance.tsv").read_text().splitlines() if x and not x.startswith("#")]
    head = lines[0].split("\t")
    out = {}
    for line in lines[1:]:
        d = dict(zip(head, line.split("\t")))
        out[d["gene"]] = d
    return out


def recessive_genes() -> set[str]:
    return set(inheritance_table()) | set(carrier_genes())


def _is_risk_allele(significance: str) -> bool:
    s = str(significance).lower()
    return "risk" in s or "low_penetrance" in s


def _status(moi: str, pattern: str, chrom: str, sex: str | None) -> str:
    if "XL" in moi and str(chrom).replace("chr", "") == "X":
        if sex == "male":
            return "one X copy with the variant (hemizygous): this is the affected genotype in males, review it"
        if pattern == "homozygous":
            return "two copies: review"
        return "one copy of an X-linked variant: can be passed to children; some X-linked conditions also affect female carriers"
    if pattern == "homozygous":
        return "two copies: this is the affected genotype for a recessive condition, review it"
    if pattern.startswith("possible compound"):
        return "two different variants in this gene: carrier if on the same chromosome, affected genotype if on opposite ones (phase unknown)"
    return "carrier (one copy)"


def findings(cv_hits, sex: str | None = None) -> dict:
    """Split ClinVar hits into carrier results, risk alleles in recessive genes, and hits in genes with no recessive curation."""
    inh, acmg = inheritance_table(), carrier_genes()
    out = {"carrier": [], "risk_allele": [], "uncurated": []}
    if cv_hits is None or len(cv_hits) == 0:
        return out
    for h in cv_hits.to_dict("records"):
        gene = h["gene"]
        curated = inh.get(gene)
        if curated is None and gene not in acmg:
            out["uncurated"].append(h)
            continue
        moi = curated["moi"] if curated else ("XL" if "(XL)" in acmg[gene] else "AR")
        row = dict(h, moi=moi, also_dominant=bool(curated and curated["also_dominant"]), acmg_panel=gene in acmg,
                   diseases=(curated["diseases"] if curated else acmg[gene]),
                   status=_status(moi, h.get("pattern", ""), h.get("chrom", ""), sex))
        risk = _is_risk_allele(h.get("significance", ""))
        if risk:
            copies = "two copies" if h.get("pattern") == "homozygous" else "one copy"
            row["status"] = f"{copies} of a risk allele: can raise your own risk, and can be passed to children"
        out["risk_allele" if risk else "carrier"].append(row)
    out["carrier"].sort(key=lambda r: (r["status"] == "carrier (one copy)", not r["acmg_panel"], r["gene"]))
    return out
