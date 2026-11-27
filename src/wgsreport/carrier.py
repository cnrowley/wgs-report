"""Carrier status for family planning: ClinVar pathogenic / likely pathogenic variants in recessive and X-linked genes.

Which genes: ClinGen Definitive/Strong/Moderate AR or XL curations (data/gene_inheritance.tsv) plus the ACMG tier-3 carrier list
(Gregg 2021). Variants are NOT filtered by population frequency: common carrier alleles matter most for family planning.

This is not a carrier screen of laboratory standard. It matches known ClinVar SNVs and short indels only. For several of the most
important carrier conditions the usual cause is a deletion, duplication, rearrangement or repeat that this match cannot see; those are
listed in GAPS and shown with every report, whatever the result.
"""
from __future__ import annotations

from .clinvar import carrier_genes
from .registry import data_path

# Genes whose carrier status this pipeline cannot establish, or only partly, and why. Qualitative only.
GAPS = {
    "SMN1": ("spinal muscular atrophy", "not assessed", "carriers usually have a deleted or converted SMN1 copy; the SMN1+SMN2 depth total cannot show this"),
    "HBA1": ("alpha-thalassemia", "not assessed", "usually caused by deletions; the depth estimate is not a carrier test"),
    "HBA2": ("alpha-thalassemia", "not assessed", "usually caused by deletions; the depth estimate is not a carrier test"),
    "FMR1": ("fragile X syndrome", "repeat module only", "a CGG repeat: only the repeat-expansion module sizes it, and AGG interruptions are not read"),
    "FXN": ("Friedreich ataxia", "partial", "mostly a GAA repeat (repeat-expansion module); point variants are matched here"),
    "DMD": ("Duchenne / Becker muscular dystrophy", "partial", "exon deletions and duplications, a main cause, are not detected"),
    "F8": ("hemophilia A", "partial", "inversions and large deletions are not detected"),
    "CYP21A2": ("congenital adrenal hyperplasia", "unreliable", "a near-identical pseudogene causes gene conversions that short reads mis-map"),
    "GBA1": ("Gaucher disease", "unreliable", "a pseudogene next door causes recombinant alleles that short reads mis-map"),
    "GBA": ("Gaucher disease", "unreliable", "a pseudogene next door causes recombinant alleles that short reads mis-map"),
    "TNXB": ("TNXB-related Ehlers-Danlos syndrome", "unreliable", "overlaps a pseudogene region that short reads mis-map"),
    "RPGR": ("X-linked retinitis pigmentosa", "partial", "a repetitive region (ORF15) is poorly read by short reads"),
    "CFTR": ("cystic fibrosis", "partial", "known SNVs and small indels are matched; exon deletions and duplications are not detected"),
    "ARX": ("ARX-related disorders", "partial", "tract expansions and structural changes may be missed"),
    "AFF2": ("fragile XE", "partial", "a repeat expansion; only sized if the repeat module covers it"),
}


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


def panel_coverage(cv_hits) -> list[dict]:
    """One row per ACMG tier-3 gene: what this pipeline can say about it."""
    found = set() if cv_hits is None or len(cv_hits) == 0 else set(cv_hits["gene"])
    rows = []
    for gene, condition in carrier_genes().items():
        gap = GAPS.get(gene)
        if gene in found:
            result = "ClinVar variant found (see above)"
        elif gap and gap[1] in ("not assessed", "unreliable", "repeat module only"):
            result = f"{gap[1]}: {gap[2]}"
        else:
            result = "no known ClinVar SNV/short indel found" + (f"; {gap[2]}" if gap else "; deletions and duplications not assessed")
        rows.append({"gene": gene, "condition": condition, "result": result})
    return rows
