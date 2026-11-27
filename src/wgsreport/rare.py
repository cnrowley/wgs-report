"""Rare predicted loss-of-function (pLoF) variants in medically actionable genes.

This catches variants that ClinVar has not (yet) classified: for dominant, loss-of-function-mechanism genes (e.g.
BRCA1/2, LDLR, MLH1) a novel nonsense/frameshift/canonical-splice variant is often pathogenic, which is the
main place whole-genome sequencing is clearly informative for an individual. Many ACMG genes act by other
mechanisms (gain of function, dominant-negative), and some LoF variants are benign or in the last exon, so every
hit here is a *candidate for review*, not a diagnosis.

Annotation uses bcftools csq with the Ensembl GRCh38 GFF3 (chromosome names prefixed with 'chr').
"""
from __future__ import annotations

import gzip
import io
from pathlib import Path

import pandas as pd

from .utils import download, log, require_tools, run

from .sources import ENSEMBL_GFF_URL as GFF_URL
LOF = ("stop_gained", "frameshift", "splice_acceptor", "splice_donor", "start_lost")


# Genes on the cancer-predisposition list whose disease mechanism is gain of function, activation or a duplication, so a truncating
# variant is not the mechanism (SCG5 is listed only because of the SCG5-GREM1 duplication). Written from the gene-disease mechanisms;
# extend it when a scan shows a false lead of this kind.
NO_LOF_MECHANISM = frozenset({"ALK", "EGFR", "KIT", "MET", "PDGFRA", "RET", "CDK4", "SCG5", "GREM1", "PRSS1"})


def drop_clinvar_matches(lof, clinvar_hits):
    """Remove scan rows that are already reported as ClinVar P/LP hits (chromosome names may differ: '15' vs 'chr15')."""
    if lof is None or len(lof) == 0 or clinvar_hits is None or len(clinvar_hits) == 0:
        return lof
    norm = lambda c: str(c).replace("chr", "")
    keys = {(norm(r.chrom), int(r.pos), r.ref, r.alt) for r in clinvar_hits.itertuples()}
    keep = [(norm(r.chrom), int(r.pos), r.ref, r.alt) not in keys for r in lof.itertuples()]
    return lof[keep].reset_index(drop=True)


def load_af_genes() -> dict:
    """AF-associated rare-LoF genes with their odds ratios (see data/af_lof_genes.tsv for the source)."""
    from .registry import data_path
    rows = [l.split("\t") for l in data_path("af_lof_genes.tsv").read_text().splitlines()
            if l and not l.startswith("#") and not l.startswith("gene\t")]
    return {r[0]: float(r[1]) for r in rows}


def prepare_gff(work: Path) -> Path:
    """Download Ensembl GFF3 and rewrite contig names to match a UCSC-style (chr1) VCF."""
    out = Path(work) / "ensembl.chr.gff3.gz"
    if out.exists():
        return out
    raw = Path(work) / "ensembl.gff3.gz"
    download(GFF_URL, raw)
    with gzip.open(raw, "rt") as fi, gzip.open(out, "wt") as fo:
        for line in fi:
            if line.startswith("#"):
                if not line.startswith("##sequence-region"):
                    fo.write(line)
            elif line.split("\t", 1)[0].isdigit() or line.split("\t", 1)[0] in ("X", "Y"):
                fo.write("chr" + line)
            elif line.startswith("MT\t"):
                fo.write("chrM" + line[2:])
    raw.unlink()
    return out


def lof_scan(vcf: Path, fasta: Path, gff: Path, genes: set[str], kg_dir: Path | None = None,
             max_af: float = 0.01) -> pd.DataFrame:
    """PASS variants with LoF consequences in `genes`, annotated with 1000 Genomes allele frequency when available."""
    require_tools("bcftools")
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%FILTER\t%BCSQ\n"
    out = run(f"bcftools view -f PASS,. {vcf} | bcftools csq -f {fasta} -g {gff} -l -p a --force -Ou "
              f"| bcftools query -f '{fmt}' -i 'INFO/BCSQ!=\".\"'").stdout
    rows = []
    for line in out.splitlines():
        chrom, pos, ref, alt, gt, flt, bcsq = line.split("\t")
        for c in bcsq.split(","):
            if c.startswith("@"):
                continue
            f = c.split("|")
            if len(f) < 3:
                continue
            cons, gene = f[0], f[1]
            if "NMD_transcript" in cons:
                continue
            if gene in genes and any(k in cons for k in LOF):
                rows.append({"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt, "gene": gene,
                             "consequence": cons, "transcript": f[2], "genotype": gt,
                             "change": f[5] if len(f) > 5 else "",
                             "flag": ("multi-allelic indel: probable repeat/alignment artefact; check reads"
                                      if "," in alt and gt in ("1/2", "1|2", "2/1", "2|1") else "")})
    df = pd.DataFrame(rows).drop_duplicates(["chrom", "pos", "alt", "gene"]) if rows else pd.DataFrame(
        columns=["chrom", "pos", "ref", "alt", "gene", "consequence", "transcript", "genotype", "change", "flag"])
    if len(df) and kg_dir is not None:
        df["af_1000g"] = [_kg_af(kg_dir, r.chrom, r.pos, r.alt) for r in df.itertuples()]
        df = df[df["af_1000g"].isna() | (df["af_1000g"] < max_af)]
    return df.reset_index(drop=True)


def _kg_af(kg_dir: Path, chrom: str, pos: int, alt: str):
    from .refpanel import DEFAULT_PATTERN
    f = Path(kg_dir) / DEFAULT_PATTERN.format(chrom=chrom.replace("chr", ""))
    if not f.exists():
        return None
    out = run(["bcftools", "query", "-r", f"{chrom.replace('chr', '')}:{pos}-{pos}", "-f", "%ALT\t%INFO/AF\n", str(f)]).stdout
    for line in out.splitlines():
        alts, afs = line.split("\t")
        for a, af in zip(alts.split(","), afs.split(",")):
            if a == alt:
                return float(af)
    return None
