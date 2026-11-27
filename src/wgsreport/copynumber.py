"""Copy-number estimates for repeat-dependent loci, from read depth in the alignment file.

Identifiability matters more than precision: for SMN1/SMN2, C4A/C4B, HBA1/HBA2, FCGR3A/B and LPA the biologically meaningful
quantity (gene-specific dose, isotype, allele size, which chromosome carries what) generally cannot be recovered from total depth,
however deep the sequencing. Treat every value as a hypothesis for a dedicated assay.

Depth is summed over *all* reads in the region (mapping quality 0 included), because reads from near-identical
copies cannot be placed uniquely: they land on one of the reference copies at random. The sum of depth over all
reference copies of a repeat is therefore proportional to the total number of copies in the person's genome.

    total copies (both chromosomes) = 2 x (sum over reference copies of mean depth) / (single-copy baseline depth)

These are estimates with a typical error of about +-1 copy. They give a *diploid total* and cannot say how copies
are split between your two chromosomes. Lp(a) in particular depends on the size of the shorter apo(a) isoform,
which this method cannot read.
"""
from __future__ import annotations

import statistics
import tempfile
from pathlib import Path

from .utils import require_tools, run

# GRCh38 gene spans (Ensembl). "expected" is the typical diploid total in people of European ancestry.
LOCI = {
    "AMY1 (A+B+C)": dict(regions=[("chr1", 103650592, 103664629), ("chr1", 103687415, 103701064),
                                  ("chr1", 103745126, 103758726)], expected="2-20 (median ~6)", confidence="low-moderate",
                         meaning="Salivary amylase gene copies. Depth gives only a diploid total of the AMY1 genes (they are near-identical "
                                 "to each other and partly to AMY2); associations with diet and BMI are inconsistent. Research interest only."),
    "SMN1 + SMN2": dict(regions=[("chr5", 70049635, 70078522), ("chr5", 70925030, 70953942)], expected="4 (2+2)", confidence="low for SMN1 alone",
                        meaning="SMN1 and SMN2 are near-identical, so this total cannot tell how many copies are SMN1 (the spinal muscular atrophy "
                                "gene). SMN1 dose needs the exon-7 base ratio or a clinical assay (MLPA/ddPCR); a silent 2+0 carrier cannot be "
                                "seen at all. Research interest only."),
    "HBA1 + HBA2": dict(regions=[("chr16", 172021, 173910), ("chr16", 176660, 177527)], expected="4", confidence="not reliable",
                        meaning="Alpha-globin genes (near-identical, GC-rich, subtelomeric). This depth estimate is NOT a carrier test: on five public "
                                "genomes with no deletion listed in the 1000 Genomes callset it read 4.5-5.3 copies (typical 4), and depth in this "
                                "region depends strongly on the sequencing library. It can neither show nor exclude an alpha-thalassemia deletion. "
                                "If carrier status matters (for example when planning children), use clinical deletion/duplication testing with sequencing."),
    "C4A + C4B": dict(regions=[("chr6", 31981991, 32005247), ("chr6", 32014725, 32035975)], expected="4 (2-8)", confidence="low-moderate",
                      meaning="Complement C4 total copies. C4A versus C4B and the long/short (HERV insertion) forms matter biologically and are not "
                              "separated here; an association of C4A expression with schizophrenia is small per copy. Research interest only."),
    "GSTM1": dict(regions=[("chr1", 109687309, 109709552)], expected="2 (0 in ~50% of Europeans)", confidence="moderate",
                  meaning="GSTM1 gene deletion is common; a halved depth suggests one copy deleted, while a homozygous deletion would drop "
                          "to near zero. Weak, inconsistent health associations."),
    "UGT2B17": dict(regions=[("chr4", 68516301, 68627039)], expected="2 (0-2)", confidence="moderate",
                    meaning="Deletion lowers testosterone glucuronidation (affects some steroid tests and drug metabolism); copy number is "
                            "usually well captured when the whole gene is deleted."),
    "RHD + RHCE": dict(regions=[("chr1", 25272393, 25330445), ("chr1", 25362249, 25430192)], expected="4 (3 = one RHD deleted; 2 = both deleted)", confidence="low",
                       meaning="Rh(D) blood group: the common Rh-negative state is a homozygous RHD deletion (about 2 total copies of RHD+RHCE). "
                               "RHD and RHCE are 97% identical, so this is a depth estimate only; confirm by serology."),
    "FCGR3B": dict(regions=[("chr1", 161623193, 161631963)], expected="2 (1-4)", confidence="low-moderate",
                   meaning="Neutrophil Fc receptor copy number (near-identical to FCGR3A); associations with autoimmune disease are weak."),
}

# LPA kringle IV exon pairs (160 + 182 bp) in the GRCh38 reference; unique exons at the 3' end (kringle V / protease).
LPA_UNIQUE = [("chr6", 160531482, 160531890), ("chr6", 160532531, 160532649), ("chr6", 160537855, 160537961),
              ("chr6", 160540043, 160540183), ("chr6", 160541107, 160541181), ("chr6", 160542688, 160542808),
              ("chr6", 160545440, 160545533), ("chr6", 160547789, 160547937)]
LPA_REPEAT_REGION = ("chr6", 160577000, 160651000)      # KIV-2 repeat block (all copies, plus a few flanking KIV types)
LPA_FLANK = ("chr6", 160651000, 160664276)

# Single-copy baseline windows: spread across autosomes, median taken to ignore gaps and repeats.
def _baseline_windows(fai: Path, n_per_chrom: int = 4, width: int = 50000):
    out = []
    for line in Path(fai).read_text().splitlines():
        c, ln = line.split("\t")[:2]
        if not c.startswith("chr") or not c[3:].isdigit():
            continue
        for k in range(1, n_per_chrom + 1):
            s = int(int(ln) * k / (n_per_chrom + 1))
            out.append((c, s, s + width))
    return out


def bedcov(cram: Path, regions: list[tuple[str, int, int]], min_mapq: int = 0) -> list[float]:
    """Mean per-base depth for each region (all reads when min_mapq=0)."""
    with tempfile.TemporaryDirectory() as td:
        bed = Path(td) / "r.bed"
        bed.write_text("".join(f"{c}\t{s - 1}\t{e}\n" for c, s, e in regions))
        out = run(["samtools", "bedcov", "-Q", str(min_mapq), str(bed), str(cram)]).stdout
    res = []
    for line in out.splitlines():
        f = line.split("\t")
        res.append(float(f[3]) / (int(f[2]) - int(f[1])))
    return res


def baseline_depth(cram: Path, fai: Path) -> float:
    d = bedcov(cram, _baseline_windows(fai))
    d = [x for x in d if x > 5]
    return statistics.median(d)


def estimate(cram: Path, fai: Path) -> dict:
    require_tools("samtools")
    base = baseline_depth(cram, fai)
    out = {"baseline_depth": base, "loci": {}}
    for name, spec in LOCI.items():
        depths = bedcov(cram, spec["regions"])
        # parts of the ref copies that are unique flanks would inflate, so use gene spans as given
        n = 2 * sum(depths) / base
        out["loci"][name] = {"estimated_total_copies": n, "expected": spec["expected"], "meaning": spec["meaning"],
                             "confidence": spec["confidence"],
                             "depths": depths}
    # LPA: kringle IV repeat block vs unique 3' exons
    uniq = bedcov(cram, LPA_UNIQUE)
    # weight by length
    ulen = [e - s + 1 for _, s, e in LPA_UNIQUE]
    u = sum(d * l for d, l in zip(uniq, ulen)) / sum(ulen)
    rep = bedcov(cram, [LPA_REPEAT_REGION])[0]
    flank = bedcov(cram, [LPA_FLANK])[0]
    block = LPA_REPEAT_REGION[2] - LPA_REPEAT_REGION[1]
    out["lpa"] = {"unique_exon_depth": u, "repeat_block_depth": rep, "flank_depth": flank,
                  "repeat_block_length": block, "baseline_depth": base,
                  "repeat_block_to_unique_ratio": rep / u}
    return out
