"""Reference-build checks. Everything here assumes GRCh38; a GRCh37/hg19 VCF would be scored at the wrong positions without any
error, so it is refused before any scoring (override: --skip-build-check).

Two independent checks:
  1. header_build(): contig lengths (and ##reference text) in the VCF header. chr1 is 248,956,422 bp in GRCh38 and 249,250,621 in GRCh37.
  2. ref_concordance(): after genotypes are read at the panel positions, the VCF's REF allele must match the panel's REF allele.
     In the right build almost every record agrees; in another build about a quarter do by chance."""
from __future__ import annotations

import re

LENGTHS = {"GRCh38": {"1": 248956422, "2": 242193529, "X": 156040895}, "GRCh37": {"1": 249250621, "2": 243199373, "X": 155270560}}
NAMES = [("GRCh38", re.compile(r"grch38|hg38|hs38|GCA_000001405\.1[5-9]|GCA_000001405\.2\d", re.I)),
         ("GRCh37", re.compile(r"grch37|hg19|b37|hs37|human_g1k_v37|GCA_000001405\.1\b", re.I))]
MAX_REF_MISMATCH = 0.05      # a correct GRCh38 VCF disagrees with the panel REF at well under 1% of sites


def header_build(header: str) -> tuple[str, str]:
    """Return (build, evidence); build is 'GRCh38', 'GRCh37' or 'unknown'. Contig lengths decide; names are used only when no length is present."""
    found = {}
    for m in re.finditer(r"##contig=<[^>]*ID=(?:chr)?(1|2|X)\b[^>]*length=(\d+)", header):
        found[m.group(1)] = int(m.group(2))
    matched = {b: [c for c, n in found.items() if ls.get(c) == n] for b, ls in LENGTHS.items()}
    if matched["GRCh38"] and matched["GRCh37"]:
        return "unknown", "header contigs mix GRCh37 and GRCh38 lengths"
    for b, hits in matched.items():
        if hits:
            others = [c for c in found if c not in hits]
            return (b, f"contig length of chr{hits[0]} is {found[hits[0]]:,}") if not others else \
                   ("unknown", f"contig lengths partly match {b} (chr{hits[0]}) but chr{others[0]} does not")
    if found:
        return "unknown", "contig lengths match neither GRCh37 nor GRCh38: " + ", ".join(f"{c}={n:,}" for c, n in found.items())
    for line in header.splitlines():
        if line.startswith("##reference") or line.startswith("##source"):
            for b, rx in NAMES:
                if rx.search(line):
                    return b, f"header says: {line[:100]}"
    return "unknown", "no contig lengths or reference line in the header"


def ref_concordance(stats: dict) -> tuple[float, int]:
    """(share of VCF records at panel positions whose REF disagrees with the panel REF, number of records compared)."""
    mism = stats.get("diploid_mismatch_ref", 0)
    seen = mism + stats.get("from_vcf_het", 0) + stats.get("from_vcf_hom_alt", 0)
    return (mism / seen if seen else 0.0), seen


def check_build(header: str, stats: dict | None = None) -> str | None:
    """Return an error message if the VCF is not GRCh38, else None. `stats` (from score.user_dosages) enables the REF-concordance check."""
    b, why = header_build(header)
    if b == "GRCh37":
        return (f"This VCF looks like GRCh37/hg19 ({why}). wgs-report needs GRCh38: lift it over (e.g. with CrossMap or bcftools +liftover) "
                "or re-call it on GRCh38. Use --skip-build-check to override at your own risk.")
    if stats is not None:
        rate, n = ref_concordance(stats)
        if n >= 200 and rate > MAX_REF_MISMATCH:
            return (f"{rate:.0%} of {n:,} VCF records disagree with the GRCh38 reference allele at the scored positions "
                    f"(header check: {b}, {why}). The VCF is probably on a different reference build. Use --skip-build-check to override.")
    return None
