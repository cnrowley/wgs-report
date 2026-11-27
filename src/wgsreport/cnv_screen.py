"""Read-depth screen for a short list of large, individually actionable copy-number variants (CNVs).

Which CNVs: data/cnv_actionable.tsv (rule in its header): strong evidence of disease, a substantial effect on a carrier, and a documented
follow-up, with the strength of that guidance shown. Weaker or non-actionable CNVs are deliberately not reported; for research use,
`wgsreport query copy-number --region` reads any region.

Sources
  Region coordinates (GRCh38): ClinGen recurrent-region curation list.
  22q11.2 deletion psychiatric and developmental estimates: Kirov, Rees & Walters, BJPsych Advances 2015;21:157-163, and Mollon et al.,
  Mol Psychiatry 2023;28:1480. Other disease links: data/cnv_disease.tsv (inclusion rule and sources in its header).

Method
  Mean depth in 10 kb bins, reads with MAPQ >= 30, relative to the autosomal median = 2 copies. These regions are flanked by segmental
  duplications, so a bin is used only if at least 90% of its total depth survives the MAPQ filter. A heterozygous deletion leaves
  that fraction unchanged, so the filter does not remove the signal.

What it cannot do: it is not a CNV caller and not a genome-wide search. It cannot give breakpoints, and can miss events shorter than
about 200 kb, mosaic events and exon-level changes in disease genes. Bins without reads are not used, so loss of both copies shows as
'not assessable'. Any call needs an independent assay (chromosomal microarray or MLPA).
"""
from __future__ import annotations

import statistics
from pathlib import Path

from .copynumber import _baseline_windows, bedcov
from .registry import data_path
from .utils import require_tools

BIN = 10_000
MIN_MAPQ = 30
MIN_UNIQUE_FRACTION = 0.9
LOW, HIGH = 1.5, 2.5                  # per-bin copy thresholds
SEGMENT_BINS = 20                     # 200 kb of contiguous usable bins
MIN_USABLE_FRACTION = 0.7             # fewer usable bins than this: missing data or loss of both copies, never 'not detected'
MIN_BASELINE = 15.0
WORD = {"del": "deletion", "dup": "duplication"}
PSYCHIATRIC = {
    "22q11.2 (A-D; DiGeorge / velocardiofacial)": (
        "Published estimates: lifetime schizophrenia in carriers about 12% (95% CI 7-18; Kirov 2015) to 25% (Mollon 2023), against about 1% "
        "in the general population; intellectual disability, autism or congenital malformations 88% (53-100)."),
}


def _table(name: str) -> list[dict]:
    lines = [x.split("\t") for x in data_path(name).read_text().splitlines() if x and not x.startswith("#")]
    return [dict(zip(lines[0], x)) for x in lines[1:]]


def load_actionable() -> list[dict]:
    return [dict(r, start=int(r["start"]), end=int(r["end"])) for r in _table("cnv_actionable.tsv")]


def load_diseases() -> list[dict]:
    return _table("cnv_disease.tsv")


ACTIONABLE = load_actionable()
DISEASES = load_diseases()


def regions() -> list[tuple[str, int, int]]:
    """Distinct regions to measure (deletion and duplication entries share one)."""
    return sorted({(r["chrom"], r["start"], r["end"]) for r in ACTIONABLE}, key=lambda x: (int(x[0][3:]), x[1]))


def key(chrom: str, start: int, end: int) -> str:
    return f"{chrom}:{start}-{end}"


def bins(chrom: str, start: int, end: int, width: int = BIN) -> list[tuple[str, int, int]]:
    return [(chrom, s, min(s + width - 1, end)) for s in range(start, end + 1, width)]


def usable_copies(all_depth: list[float], q_depth: list[float], baseline: float) -> list[float]:
    """Copy number of every window that is uniquely mappable; the others are dropped (order is kept)."""
    out = []
    for a, q in zip(all_depth, q_depth):
        if a >= 1 and q / a >= MIN_UNIQUE_FRACTION:
            out.append(2 * q / baseline)
    return out


def _longest_run(flags: list[bool]) -> tuple[int, int]:
    best, at, n = 0, 0, 0
    for i, f in enumerate(flags):
        n = n + 1 if f else 0
        if n > best:
            best, at = n, i - n + 1
    return best, at


def classify(copies: list[float], n_total: int | None = None) -> dict:
    """Call one region from the per-bin copy numbers; n_total is the number of bins before filtering."""
    n = len(copies)
    if n < 10 or (n_total and n < MIN_USABLE_FRACTION * n_total):
        return {"call": "not assessable", "n_bins": n}
    med = statistics.median(copies)
    res = {"n_bins": n, "median_copies": med, "fraction_low": sum(c < LOW for c in copies) / n,
           "fraction_high": sum(c > HIGH for c in copies) / n, "call": "not detected"}
    for kind, frac, test, whole in (("deletion", "fraction_low", lambda c: c < LOW, med <= 1.35),
                                    ("duplication", "fraction_high", lambda c: c > HIGH, med >= 2.65)):
        if whole and res[frac] >= 0.8:
            res["call"] = kind
            return res
        run, at = _longest_run([test(c) for c in copies])
        if run >= min(SEGMENT_BINS, max(10, n // 2)):
            res.update(call=f"partial {kind}", segment_bins=run, segment_median_copies=statistics.median(copies[at:at + run]))
            return res
    return res


def _x_windows() -> list[tuple[str, int, int]]:
    """10 kb windows across non-pseudoautosomal chrX: a built-in one-copy control in XY genomes."""
    return [("chrX", s, s + BIN - 1) for s in range(10_000_000, 150_000_000, 1_000_000)]


def estimate(cram: Path, fai: Path) -> dict:
    require_tools("samtools")

    def depth(windows):
        return bedcov(cram, windows, 0), bedcov(cram, windows, MIN_MAPQ)

    base_windows = [w for c, s, e in _baseline_windows(fai) for w in bins(c, s, e - 1)]
    a, q = depth(base_windows)
    kept = [qq for aa, qq in zip(a, q) if aa > 5 and qq / aa >= MIN_UNIQUE_FRACTION]
    if not kept or statistics.median(kept) < MIN_BASELINE:
        return {"error": f"median autosomal depth below {MIN_BASELINE:.0f}x: not assessed"}
    base = statistics.median(kept)
    out = {"baseline_depth": base, "regions": {}}
    for chrom, start, end in regions():
        w = bins(chrom, start, end)
        out["regions"][key(chrom, start, end)] = classify(usable_copies(*depth(w), base), len(w))
    xc = usable_copies(*depth(_x_windows()), base)
    out["x_control_copies"] = statistics.median(xc) if xc else None
    return out


def findings(res: dict) -> list[dict]:
    """One row per actionable CNV. A CNV is flagged only when the observed direction is the actionable one.

    A whole-region call in that direction is 'confirmable'; a partial segment is flagged but may not involve the critical genes.
    Psychiatric estimates are withheld when the result carries psychiatric_estimates = False."""
    rows = []
    for spec in ACTIONABLE:
        r = res.get("regions", {}).get(key(spec["chrom"], spec["start"], spec["end"]))
        if r is None:
            r = {"call": "not assessable"}
        call, word = r["call"], WORD[spec["change"]]
        flag = call.endswith(word)
        diseases = [d for d in DISEASES if flag and d["region"] == spec["name"] and spec["change"] in d["change"].split(",")]
        rows.append(dict(spec, call=call if flag or call == "not assessable" else "not detected", copies=r.get("median_copies"), flag=flag,
                         partial=flag and call.startswith("partial"), diseases=diseases,
                         psychiatric=PSYCHIATRIC.get(spec["name"], "") if flag and res.get("psychiatric_estimates", True) else ""))
    return rows
