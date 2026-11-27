"""Whole-genome sanity checks that are deterministic from the data: sex-chromosome dosage, whole-chromosome gains/losses,
runs of homozygosity (parental relatedness). No calibration cohort is needed, and the answers are usually unambiguous.

Chromosome dosage uses read depth in sampled windows (`samtools bedcov`), each chromosome relative to the median autosome. It detects whole-chromosome changes (for example XXY, X0, trisomy 21). It does not detect mosaicism below roughly 20%, small
copy-number changes or balanced rearrangements, and read-count ratios shift a few percent with GC content and mappability.

Runs of homozygosity come from `bcftools roh` on autosomal SNVs. Long runs (>= 1.5 Mb) mean the two parental haplotypes were
inherited from a shared ancestor; the autosomal fraction F_ROH estimates the inbreeding coefficient (McQuillan et al., AJHG 2008,
doi:10.1016/j.ajhg.2008.02.004).
"""
from __future__ import annotations

import statistics
import subprocess
from pathlib import Path

from .utils import require_tools, run

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
PAR_X = ((10001, 2781479), (155701383, 156030895))          # GRCh38 pseudoautosomal regions on X
X_LEN = 156040895
# Expected inbreeding coefficients for the offspring of relatives (theory): used only to put F_ROH in context.
F_EXPECTED = (("first cousins", 1 / 16), ("second cousins", 1 / 64), ("third cousins", 1 / 256))


def _windows(chrom: str, start: int, end: int, n: int = 40, width: int = 50_000):
    step = (end - start - width) // max(n - 1, 1)
    return [(chrom, start + i * step, start + i * step + width) for i in range(n)]


def classify(ratio: dict) -> tuple[str, str, list]:
    """(sex-chromosome call, note, flagged autosomes) from depth ratios relative to the autosomal median."""
    x, y = ratio["chrX"], ratio["chrY"]
    flags = [f"{c}: {ratio[c]:.2f}x the autosomal median" for c in AUTOSOMES if c in ratio and (ratio[c] > 1.25 or ratio[c] < 0.65)]
    x_copies = 1 if x < 0.75 else (2 if x < 1.25 else 3)
    y_present = y > 0.2
    names = {(2, False): "XX", (1, True): "XY", (1, False): "X0 (Turner pattern)", (2, True): "XXY pattern",
             (3, False): "XXX pattern", (3, True): "XXXY pattern"}
    note = ("X depth suggests about %d cop%s; Y reads %s. Extra Y copies (XYY) cannot be told apart from XY by this screen."
            % (x_copies, "y" if x_copies == 1 else "ies", "present" if y_present else "absent"))
    return names[(x_copies, y_present)], note, flags


def dosage(cram: Path, fai: Path) -> dict:
    """Median depth of 40 evenly spaced 50 kb windows per chromosome (MAPQ >= 20), relative to the autosomal median.
    Medians ignore windows in gaps or repeats, so acrocentric short arms and heterochromatin do not distort the ratios.
    chrX windows skip the pseudoautosomal regions and the centromere; chrY windows cover the euchromatic arm."""
    lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in Path(fai).read_text().splitlines()}
    spec = {c: _windows(c, int(lens[c] * 0.05), int(lens[c] * 0.95)) for c in AUTOSOMES if c in lens}
    spec["chrX"] = _windows("chrX", 3_000_000, 155_000_000)
    spec["chrY"] = _windows("chrY", 3_000_000, 25_000_000)
    from .copynumber import bedcov
    med = {}
    for c, regs in spec.items():
        d = bedcov(cram, regs, min_mapq=20)
        med[c] = statistics.median(d)
    base = statistics.median(med[c] for c in AUTOSOMES if c in med)
    ratio = {c: med[c] / base for c in med}
    call, note, flags = classify(ratio)
    return {"ratios": {c: round(r, 3) for c, r in ratio.items()}, "sex_chromosome_call": call, "note": note,
            "autosome_flags": flags, "baseline_depth": round(base, 1), "method": "median of 40 windows x 50 kb per chromosome, MAPQ>=20"}


def x_het_fraction(vcf: Path) -> dict:
    """Heterozygous fraction of non-PAR chrX SNV calls: near zero in a person with one X, about a third with two."""
    out = run(["bcftools", "view", "-H", "-v", "snps", "-f", "PASS,.", "-r", f"chrX:{PAR_X[0][1] + 1}-{PAR_X[1][0] - 1}", str(vcf)]).stdout
    n = het = 0
    for line in out.splitlines():
        gt = line.split("\t")[9].split(":")[0].replace("|", "/")
        if gt in ("0/1", "1/0", "1/2"):
            het += 1
            n += 1
        elif gt in ("1/1", "0/0", "1"):
            n += 1
    return {"n_calls": n, "het_fraction": round(het / n, 4) if n else None}


def roh(vcf: Path, fai: Path, min_len: int = 1_500_000, min_quality: float = 30.0) -> dict:
    """Autosomal runs of homozygosity from genotypes only (bcftools roh, default allele frequency 0.4)."""
    require_tools("bcftools")
    lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in Path(fai).read_text().splitlines()}
    auto_len = sum(lens[c] for c in AUTOSOMES if c in lens)
    regions = ",".join(AUTOSOMES)
    cmd = (f"bcftools view -f PASS,. -v snps -m2 -M2 -r {regions} '{vcf}' -Ou | "
           f"bcftools roh --GTs-only 30 --AF-dflt 0.4 -O r -")
    out = subprocess.run(["bash", "-c", "set -o pipefail; " + cmd], capture_output=True, text=True, check=True).stdout
    return parse_roh(out, auto_len, min_len, min_quality)


def parse_roh(text: str, auto_len: int, min_len: int = 1_500_000, min_quality: float = 30.0) -> dict:
    """Summarise `bcftools roh -O r` output (RG lines: sample, chrom, start, end, length, markers, quality)."""
    runs = []
    for line in text.splitlines():
        f = line.split("\t")
        if f[0] == "RG":
            length, q = int(f[5]), float(f[7])
            if length >= min_len and q >= min_quality:
                runs.append({"chrom": f[2], "start": int(f[3]), "end": int(f[4]), "length": length, "markers": int(f[6]), "quality": q})
    total = sum(r["length"] for r in runs)
    froh = total / auto_len if auto_len else None
    return {"n_runs": len(runs), "total_bp": total, "f_roh": froh, "longest_bp": max((r["length"] for r in runs), default=0),
            "runs": sorted(runs, key=lambda r: -r["length"])[:15], "min_run_bp": min_len,
            "context": [{"relationship": n, "expected_f": f} for n, f in F_EXPECTED]}


def run_qc(cram: Path | None, vcf: Path, fai: Path) -> dict:
    out = {}
    if cram:
        require_tools("samtools")
        out["dosage"] = dosage(cram, fai)
    out["x_het"] = x_het_fraction(vcf)
    out["roh"] = roh(vcf, fai)
    return out
