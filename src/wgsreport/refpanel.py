"""Build a compact PLINK2 reference panel from 1000 Genomes (GRCh38) at the variants we score.

The panel serves three purposes:
  * defines the allele coding (REF/ALT) used to resolve "absent from the VCF => homozygous reference";
  * provides the population distribution of every polygenic score (for z-scores / percentiles);
  * provides the PCA space for the ancestry estimate.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from .utils import log, require_tools, run

# 1000G phase-3 GRCh38 lift-over release (EBI). Files are named per chromosome.
DEFAULT_PATTERN = "ALL.chr{chrom}_GRCh38.genotypes.20170504.vcf.gz"
PANEL_FILE = "integrated_call_samples_v3.20130502.ALL.panel"
SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]


def load_panel_meta(kg_dir: Path) -> pd.DataFrame:
    cands = [Path(kg_dir) / PANEL_FILE, Path(kg_dir).parent / PANEL_FILE]
    p = next((c for c in cands if c.exists()), cands[0])
    if not p.exists():
        raise SystemExit(f"Sample panel file not found: {p}\n"
                         "Download it from the 1000 Genomes FTP (see README).")
    df = pd.read_csv(p, sep="\t", usecols=["sample", "pop", "super_pop", "gender"])
    return df.rename(columns={"sample": "IID"})


def _build_chrom(chrom: str, regions: pd.DataFrame, kg_dir: Path, pattern: str,
                 out_dir: Path, threads: int) -> Path:
    vcf = Path(kg_dir) / pattern.format(chrom=chrom)
    if not vcf.exists():
        raise SystemExit(f"Missing reference VCF: {vcf}")
    reg = out_dir / f"regions.chr{chrom}.tsv"
    regions.to_csv(reg, sep="\t", header=False, index=False)
    tmp = out_dir / f"tmp.chr{chrom}.vcf.gz"
    prefix = out_dir / f"chr{chrom}"
    if prefix.with_suffix(".pgen").exists():
        return prefix
    run(f"bcftools view -R {reg} -m2 -M2 -v snps {vcf} -Ou "
        f"| bcftools annotate -x INFO,^FMT/GT -Oz -o {tmp}")
    run(["plink2", "--vcf", str(tmp), "--set-all-var-ids", "@:#:$r:$a",
         "--rm-dup", "force-first", "--vcf-half-call", "missing",
         "--make-pgen", "--threads", str(threads), "--out", str(prefix)])
    tmp.unlink()
    reg.unlink()
    return prefix


def build_panel(universe: pd.DataFrame, kg_dir: Path, out_dir: Path,
                pattern: str = DEFAULT_PATTERN, jobs: int = 6, threads: int = 4) -> Path:
    """Extract biallelic SNVs at `universe` (chrom,pos) positions into out_dir/panel.{pgen,pvar,psam}."""
    require_tools("bcftools", "plink2")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    final = out_dir / "panel"
    import hashlib
    fp = hashlib.sha256(universe.sort_values(["chrom", "pos"]).to_csv(index=False).encode()).hexdigest()
    fpfile = out_dir / "universe.sha256"
    if final.with_suffix(".pgen").exists() and (not fpfile.exists() or fpfile.read_text() != fp):
        log("score set changed; rebuilding reference panel")
        for f in out_dir.glob("*"):
            f.unlink()
    fpfile.write_text(fp)
    if not final.with_suffix(".pgen").exists():
        log(f"extracting {len(universe):,} positions from 1000 Genomes ({jobs} chromosomes in parallel)")
        by_chrom = {c: g[["chrom", "pos"]].drop_duplicates().sort_values("pos")
                    for c, g in universe.groupby("chrom")}
        with ThreadPoolExecutor(jobs) as ex:
            prefixes = list(ex.map(
                lambda c: _build_chrom(c, by_chrom[c], kg_dir, pattern, out_dir, threads),
                sorted(by_chrom, key=int)))
        lst = out_dir / "merge.list"
        lst.write_text("\n".join(str(p) for p in prefixes) + "\n")
        run(["plink2", "--pmerge-list", str(lst), "--make-pgen", "--out", str(final)])
        for p in prefixes:
            for ext in (".pgen", ".pvar", ".psam", ".log"):
                Path(str(p) + ext).unlink(missing_ok=True)
    meta = load_panel_meta(kg_dir)
    meta.to_csv(out_dir / "samples.tsv", sep="\t", index=False)
    for sp in SUPERPOPS:
        if (out_dir / f"freq_{sp}.afreq").exists():
            continue
        keep = out_dir / f"keep_{sp}.txt"
        pd.DataFrame({"#IID": meta.loc[meta.super_pop == sp, "IID"]}).to_csv(keep, sep="\t", index=False)
        run(["plink2", "--pfile", str(final), "--keep", str(keep), "--freq",
             "--out", str(out_dir / f"freq_{sp}")])
    if not (out_dir / "freq_ALL.afreq").exists():
        run(["plink2", "--pfile", str(final), "--freq", "--out", str(out_dir / "freq_ALL")])
    return final


def read_pvar(panel: Path) -> pd.DataFrame:
    df = pd.read_csv(str(panel) + ".pvar", sep="\t", comment="#", header=None,
                     names=["chrom", "pos", "id", "ref", "alt"], usecols=range(5), dtype={"chrom": str})
    return df


def read_afreq(out_dir: Path, pop: str) -> pd.Series:
    f = pd.read_csv(Path(out_dir) / f"freq_{pop}.afreq", sep="\t")
    return f.set_index("ID")["ALT_FREQS"]
