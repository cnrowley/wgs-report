"""Checks that read the alignment file (CRAM/BAM) directly.

A variant-only VCF cannot distinguish "homozygous reference" from "not callable". For the curated SNP panel we
therefore count reads at each site. We also call mitochondrial variants and build a primary-contig GRCh38 FASTA
from the public EBI CRAM reference registry (sequences addressed by the MD5 in the CRAM header).
"""
from __future__ import annotations

import json
import re
import urllib.request
from pathlib import Path

import pandas as pd

from .utils import log, require_tools, run

from .sources import EBI_MD5_URL as EBI_MD5
PRIMARY = re.compile(r"^chr([0-9]+|X|Y|M)$")


def fetch_reference(cram: Path, out_fasta: Path) -> Path:
    """Download the primary contigs named in the CRAM header (verified against the header MD5) into one FASTA."""
    import hashlib
    require_tools("samtools")
    out_fasta = Path(out_fasta)
    if out_fasta.exists():
        return out_fasta
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    hdr = run(["samtools", "view", "-H", str(cram)]).stdout
    tmp = out_fasta.with_suffix(".part")
    with open(tmp, "w") as fh:
        for line in hdr.splitlines():
            if not line.startswith("@SQ"):
                continue
            f = dict(x.split(":", 1) for x in line.split("\t")[1:])
            if not PRIMARY.match(f["SN"]):
                continue
            log(f"  fetching {f['SN']}")
            with urllib.request.urlopen(EBI_MD5.format(md5=f["M5"]), timeout=300) as r:
                seq = r.read().decode().strip().replace("\n", "").upper()
            if hashlib.md5(seq.encode()).hexdigest() != f["M5"] or len(seq) != int(f["LN"]):
                raise RuntimeError(f"checksum mismatch for {f['SN']}")
            fh.write(f">{f['SN']}\n")
            fh.writelines(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80))
    tmp.rename(out_fasta)
    run(["samtools", "faidx", str(out_fasta)])
    return out_fasta


def _count_bases(ref: str, s: str) -> dict:
    """Count read bases in a samtools mpileup column (indels are skipped, not counted)."""
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    i = 0
    while i < len(s):
        c = s[i]
        if c == "^":
            i += 2
        elif c == "$":
            i += 1
        elif c in "+-":
            m = re.match(r"[+-](\d+)", s[i:])
            i += len(m.group(0)) + int(m.group(1))
        else:
            b = ref if c in ".," else c.upper()
            if b in counts:
                counts[b] += 1
            i += 1
    return counts


def call_from_counts(counts: dict, ref: str, alt: list[str], haploid: bool = False,
                     min_depth: int = 8) -> tuple[list[str] | None, str]:
    """Simple allele-fraction caller for a known site: returns (alleles, note)."""
    depth = sum(counts.values())
    if depth < min_depth:
        return None, f"depth {depth} too low"
    top = sorted(counts.items(), key=lambda kv: -kv[1])
    (a1, n1), (a2, n2) = top[0], top[1]
    f2 = n2 / depth
    if haploid:
        return [a1] * 2, f"depth {depth}"
    if f2 < 0.15:
        return [a1, a1], f"depth {depth}"
    if f2 >= 0.25 and n1 / depth >= 0.25:
        return sorted([a1, a2]), f"depth {depth}, {a1}:{n1} {a2}:{n2}"
    return None, f"ambiguous allele fractions ({a1}:{n1}, {a2}:{n2})"


def panel_calls(cram: Path, fasta: Path, panel: pd.DataFrame, sex: str = "male") -> dict:
    """Genotype every panel SNP from reads. Returns {rsid: {alleles, depth, counts, note}}."""
    require_tools("samtools")
    from concurrent.futures import ThreadPoolExecutor

    def one(r):
        reg = f"chr{r.chrom}:{int(r.pos)}-{int(r.pos)}"
        out = run(["samtools", "mpileup", "-f", str(fasta), "-r", reg, "-Q", "20", "-q", "20", "-B",
                   "-d", "100000", "-a", str(cram)]).stdout
        return out.splitlines()[-1] if out.strip() else ""

    with ThreadPoolExecutor(6) as ex:
        lines = list(ex.map(one, panel.itertuples()))
    out = "\n".join(lines)
    cols = {}
    for line in out.splitlines():
        f = line.split("\t")
        cols[(f[0], int(f[1]))] = (f[2].upper(), f[4] if len(f) > 4 else "")
    res = {}
    for r in panel.itertuples():
        ref, s = cols.get((f"chr{r.chrom}", int(r.pos)), ("N", ""))
        counts = _count_bases(ref, s)
        haploid = r.chrom in ("X", "Y") and sex == "male"
        alleles, note = call_from_counts(counts, ref, str(r.alt).split(","), haploid)
        res[r.rsid] = {"alleles": alleles, "depth": sum(counts.values()), "counts": counts, "note": note,
                       "chrom": r.chrom, "pos": int(r.pos)}
    return res


def mito_variants(cram: Path, fasta: Path, min_qual: float = 50) -> dict:
    """Haploid chrM variant calls with allele fractions (heteroplasmy-aware depth, not a haplogroup call)."""
    require_tools("bcftools", "samtools")
    mp = run(f"bcftools mpileup -f {fasta} -r chrM -Q 20 -q 20 -a AD -d 100000 -Ou {cram} "
             f"| bcftools call -mv --ploidy 1 -Ov").stdout
    depth = run(["samtools", "depth", "-r", "chrM", "-a", str(cram)]).stdout.split("\n")
    d = [int(x.split("\t")[2]) for x in depth if x]
    rows = []
    for line in mp.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if float(f[5]) < min_qual:
            continue
        ad = dict(zip(f[8].split(":"), f[9].split(":"))).get("AD", "")
        nums = [int(x) for x in ad.split(",") if x != "."]
        af = sum(nums[1:]) / sum(nums) if nums and sum(nums) else None
        rows.append({"pos": int(f[1]), "ref": f[3], "alt": f[4], "qual": float(f[5]), "alt_fraction": af})
    return {"mean_depth": sum(d) / len(d) if d else None, "variants": rows}


def run_checks(cram: Path, fasta: Path, panel: pd.DataFrame, out_json: Path, sex: str = "male") -> dict:
    log("genotyping the SNP panel from reads")
    calls = panel_calls(cram, fasta, panel, sex)
    log("calling mitochondrial variants")
    mt = mito_variants(cram, fasta)
    out = {"cram": Path(cram).name, "panel": calls, "mito": mt}
    Path(out_json).write_text(json.dumps(out, indent=1))
    return out
