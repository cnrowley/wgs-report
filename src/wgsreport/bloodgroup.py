"""ABO blood group inferred from three variants (GRCh38; the ABO gene is on the reverse strand).

  rs8176719  chr9:133257521 T>TC   reference allele (no C) = O1 (c.261delG); TC = intact allele (A or B)
  rs8176746  chr9:133255935 G>T    T = B-defining c.796A (together with rs8176747)
  rs8176747  chr9:133255928 C>G    G = B-defining c.803C

Orientation was checked against 1000 Genomes: the TC allele has AF 0.34 globally and 0.39 in Europeans (non-O), the B
markers 0.15 / 0.08. The inference ignores rarer alleles (A2, cis-AB, O2, weak subgroups, O1v), so it is an estimate to
confirm with a blood test. Rh (D) status is estimated separately from RHD copy number (copynumber.py).
"""
from __future__ import annotations

POSITIONS = {("9", 133257521), ("9", 133255935), ("9", 133255928)}


def _alt_count(sites: dict, key, ref: str, alt: str):
    """Number of copies of `alt` (0-2) or None if filtered/no-call. Absent record = homozygous reference."""
    recs = sites.get(key, [])
    for r, alts, gt, ok in recs:
        if r == ref and alt in alts:
            if not ok or any(g is None for g in gt):
                return None, "filtered/no-call"
            return sum(1 for g in gt if g == alts.index(alt) + 1), "called"
    if recs:
        return None, "different alleles at the site"
    return 0, "assumed homozygous reference (absent from VCF)"


def infer_abo(sites: dict) -> dict:
    n_intact, s1 = _alt_count(sites, ("9", 133257521), "T", "TC")
    n_b1, s2 = _alt_count(sites, ("9", 133255935), "G", "T")
    n_b2, s3 = _alt_count(sites, ("9", 133255928), "C", "G")
    if None in (n_intact, n_b1, n_b2):
        return {"group": "undetermined", "genotype": "", "note": f"{s1}; {s2}; {s3}"}
    n_o = 2 - n_intact
    n_b = max(n_b1, n_b2)                      # both markers define B; take the larger as a guard against one miscall
    n_a = n_intact - n_b
    if n_a < 0:
        return {"group": "undetermined", "genotype": "", "note": "B markers exceed non-O alleles (rare allele or miscall)"}
    alleles = sorted(["O"] * n_o + ["A"] * n_a + ["B"] * n_b)
    geno = "/".join(alleles)
    group = "AB" if n_a and n_b else "B" if n_b else "A" if n_a else "O"
    note = "; ".join(sorted({s1, s2, s3}))
    phase_note = ""
    if n_o == 1 and n_b == 1 and n_a == 0:
        phase_note = "B marker assumed to lie on the non-O chromosome."
    return {"group": group, "genotype": geno, "note": (note + ". " + phase_note).strip(". ") + ".",
            "disclaimer": "Inferred from three variants; rare ABO alleles are not considered. Confirm with a blood test."}
