"""Classical HLA typing from a CRAM with T1K (Song et al., Genome Res 2023, doi:10.1101/gr.277585.122), plus a table of the
HLA alleles with guideline-level clinical associations.

T1K (GPL-3.0) is not bundled. Build it yourself (`make`; see sources.T1K_URL) and create its HLA index
(`perl t1k-build.pl -o hlaidx --download IPD-IMGT/HLA`); pass the T1K directory and `hlaidx_dna_seq.fa`.

Reads from the MHC region and all unmapped reads are extracted with samtools and genotyped against the IPD-IMGT/HLA alleles.
Typing from 30x short reads is usually right at the 2-field level (protein) for the common classical genes but is
less reliable for rare alleles and for DRB3/4/5 and the DQ/DP genes; phase between genes is not known. Treat a result that
matters (a drug hypersensitivity allele above all) as a lead to confirm with a clinical HLA lab.
"""
from __future__ import annotations

import re
import shlex
import subprocess
from pathlib import Path

from .utils import log, require_tools

MHC = "chr6:28500000-33500000"          # GRCh38 MHC, with margin

# (gene, allele prefix) -> (condition or drug, source). Two-field matches (e.g. B*57:01 matches B*57:01:01).
ASSOCIATIONS = {
    ("B", "57:01"): ("Abacavir hypersensitivity (avoid abacavir)", "CPIC guideline, strong"),
    ("B", "15:02"): ("Carbamazepine / oxcarbazepine / phenytoin: Stevens-Johnson syndrome and toxic epidermal necrolysis, "
                     "highest in people of Asian ancestry", "CPIC guideline, strong"),
    ("A", "31:01"): ("Carbamazepine hypersensitivity (several reaction types)", "CPIC guideline, strong"),
    ("B", "58:01"): ("Allopurinol severe cutaneous adverse reactions (avoid or use with caution)", "CPIC guideline, strong"),
    ("B", "27:05"): ("Ankylosing spondylitis and related spondyloarthritis (risk allele; most carriers never develop it)", "association"),
}
# Celiac disease risk haplotypes are heterodimers (DQA1 + DQB1); the pairing (cis vs trans) cannot be read from short reads.
CELIAC = {"DQ2.5": ("05:01", "02:01"), "DQ2.2": ("02:01", "02:02"), "DQ8": ("03:01", "03:02")}


def extract_reads(cram: Path, fasta: Path, out_prefix: Path, threads: int = 4) -> tuple[Path, Path]:
    """FASTQ pair with MHC-region reads and unmapped reads (reads from HLA alt contigs need an alt-aware reference)."""
    require_tools("samtools")
    out_prefix = Path(out_prefix)
    r1, r2 = Path(f"{out_prefix}_1.fq"), Path(f"{out_prefix}_2.fq")
    q = shlex.quote   # paths are quoted so file names with spaces or shell metacharacters cannot alter the command
    sh = (f"set -e; t=$(mktemp -d); trap 'rm -rf $t' EXIT; "
          f"samtools view -@{threads} -u -T {q(str(fasta))} {q(str(cram))} {MHC} > $t/a.bam; "
          f"samtools view -@{threads} -u -f 4 -T {q(str(fasta))} {q(str(cram))} '*' > $t/b.bam; "
          f"samtools merge -@{threads} -u -f $t/m.bam $t/a.bam $t/b.bam; "
          f"samtools collate -@{threads} -u -O $t/m.bam $t/c | samtools fastq -@{threads} -n -1 {q(str(r1))} -2 {q(str(r2))} -0 /dev/null -s /dev/null -")
    log("HLA: extracting MHC and unmapped reads from the CRAM")
    subprocess.run(["bash", "-c", sh], check=True)
    return r1, r2


def run_t1k(t1k_dir: Path, index_fa: Path, r1: Path, r2: Path, out_dir: Path, threads: int = 4, prefix: str = "t1k") -> Path:
    exe = (Path(t1k_dir) / "run-t1k").resolve()
    index_fa, r1, r2 = Path(index_fa).resolve(), Path(r1).resolve(), Path(r2).resolve()
    if not exe.exists():
        raise SystemExit(f"{exe} not found. Build T1K (https://github.com/mourisl/T1K) and pass its directory.")
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    log("HLA: genotyping with T1K")
    subprocess.run([str(exe), "-1", str(r1), "-2", str(r2), "-f", str(index_fa), "--preset", "hla-wgs", "-t", str(threads),
                    "-o", prefix], check=True, cwd=str(out_dir))
    return out_dir / f"{prefix}_genotype.tsv"


CLASSICAL = ("A", "B", "C", "DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1")   # shown in the report; all genes stay in the JSON


def _slot(field: str) -> list[str]:
    """One allele slot of the T1K table: 'A*01:01:01' or tied candidates 'DQB1*02:02:01,DQB1*02:02:15,DQB1*03:03:02'."""
    return [x.strip().replace("HLA-", "") for x in field.split(",") if x.strip() and x.strip() != "."]


def parse_genotype(path: Path) -> list[dict]:
    """T1K genotype table -> [{gene, alleles, slots, tied, zygosity, quality, low_confidence}].
    `slots` keeps every tied candidate per allele slot, `alleles` the first of each, `tied` is True when a slot has candidates that
    differ at the two-field (protein) level. Alleles with quality <= 0 are ignored, as T1K advises."""
    rows = []
    for line in Path(path).read_text().splitlines():
        f = line.split("\t")
        if len(f) < 8 or not f[0] or f[0].startswith("#"):
            continue
        gene = f[0].replace("HLA-", "")
        slots, quals = [], []
        for name, q in ((f[2], f[4]), (f[5], f[7])):
            cands = _slot(name)
            if cands and float(q) > 0:
                slots.append(cands)
                quals.append(float(q))
        if not slots:
            rows.append({"gene": gene, "alleles": [], "slots": [], "tied": False, "zygosity": "not typed", "quality": None, "low_confidence": True})
            continue
        q = min(quals)
        firsts = [c[0] for c in slots]
        zyg = "heterozygous" if len(slots) == 2 and firsts[0] != firsts[1] else "homozygous (or second allele undetected)"
        tied = any(len({_two_field(x) for x in c}) > 1 for c in slots)
        rows.append({"gene": gene, "alleles": firsts, "slots": slots, "tied": tied, "zygosity": zyg, "quality": q, "low_confidence": q < 10})
    return rows


def _two_field(allele: str) -> str:
    m = re.match(r"^[A-Za-z0-9]+\*(\d+:\d+)", allele)
    return m.group(1) if m else ""


def _copies(call: dict, fld: str) -> tuple[int, bool]:
    """(slots containing a candidate with this two-field allele, whether any such slot is a tie that also has other alleles)."""
    n, amb = 0, False
    for slot in call.get("slots") or [[a] for a in call.get("alleles", [])]:
        hit = [x for x in slot if _two_field(x) == fld]
        if hit:
            n += 1
            amb = amb or len({_two_field(x) for x in slot}) > 1
    if n == 1 and call["zygosity"].startswith("homozygous"):
        n = 2
    return n, amb


def interpret(calls: list[dict]) -> dict:
    """Flag the alleles in ASSOCIATIONS and the celiac DQ haplotype pairs (a tied candidate list counts as compatible)."""
    by_gene = {c["gene"]: c for c in calls}
    flagged = []
    for (gene, fld), (what, src) in ASSOCIATIONS.items():
        c = by_gene.get(gene)
        if not c:
            continue
        n, amb = _copies(c, fld)
        if n:
            flagged.append({"allele": f"{gene}*{fld}", "copies": n, "meaning": what, "source": src,
                            "low_confidence": c["low_confidence"], "tied_with_other_alleles": amb})
    celiac = {k: bool(_copies(by_gene.get("DQA1", {"slots": [], "zygosity": ""}), a)[0] and _copies(by_gene.get("DQB1", {"slots": [], "zygosity": ""}), b)[0])
              for k, (a, b) in CELIAC.items()}
    return {"flagged": flagged, "celiac_compatible": celiac, "celiac_typed": bool(by_gene.get("DQA1", {}).get("alleles") and by_gene.get("DQB1", {}).get("alleles")),
            "any_low_confidence": any(c["low_confidence"] for c in calls if c["gene"] in CLASSICAL)}
