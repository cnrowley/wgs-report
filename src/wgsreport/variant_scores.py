"""Explicit-weight variant sets: a short list of (often rare, large-effect) variants whose weights are in the trait's own units.

Why this exists next to the genome-wide polygenic scores: some traits (Lp(a) above all) are dominated by a handful of variants that
a HapMap3-based genome-wide score barely weights. A published score such as PGS000667 (43 variants, weights in nmol/L) is not
standardised against a reference panel, because rare variants make the reference distribution degenerate (most people score 0).
Instead the carried variants and their summed weight are reported directly, in the score's own units.

Dosage rule (the same additive model as everywhere else): copies of the effect allele, 0/1/2. A heterozygote carries one copy and
contributes one weight; a homozygote contributes two. A site missing from a variant-only VCF is assumed homozygous reference.
Limitation: a score whose effect allele is the GRCh38 reference allele is not handled (none of the curated sets has one); such a
site would be counted 0 when absent from the VCF.
"""
from __future__ import annotations

import pandas as pd


def positions(df: pd.DataFrame) -> set[tuple[str, int]]:
    return set(zip(df["chrom"].astype(str), df["pos"].astype(int)))


def score_variant_set(df: pd.DataFrame, sites: dict, pass_only: bool = True) -> dict:
    """df: output of pgs.read_score (chrom, pos, effect, other, weight). sites: output of score.read_user_sites.

    Returns total weight carried, the carried variants (copies and contribution) and counts of what could not be assessed."""
    carried, assumed, unusable = [], 0, 0
    total = 0.0
    for r in df.itertuples(index=False):
        recs = [x for x in sites.get((str(r.chrom), int(r.pos)), []) if len(x[0]) == 1]
        if not recs:
            assumed += 1
            continue
        ref, alts, gt, ok = recs[0]
        if (pass_only and not ok) or any(g is None for g in gt):
            unusable += 1
            continue
        seq = [ref] + alts
        copies = sum(1 for g in gt if g < len(seq) and seq[g] == r.effect)
        if copies:
            total += copies * float(r.weight)
            carried.append({"rsid": getattr(r, "rsid", ""), "chrom": str(r.chrom), "pos": int(r.pos), "effect_allele": r.effect, "copies": copies,
                            "weight": float(r.weight), "contribution": copies * float(r.weight)})
    carried.sort(key=lambda x: -abs(x["contribution"]))
    return {"n_variants": len(df), "total": float(total), "carried": carried,
            "assumed_reference": assumed, "unusable": unusable}
