#!/usr/bin/env python3
"""Supervised ADMIXTURE of one genome against 1000 Genomes reference groups.

Needs: `wgsreport prepare` output (work dir), your VCF, the ADMIXTURE binary (Alexander et al. 2009; free for
academic use, not redistributable, so it is not bundled), plink2. Nothing personal is written outside --out.

Two runs:
  global   K=4 continental groups (AFR, EUR, EAS, SAS) fixed by the reference labels;
  europe   K=5 European 1000G populations (CEU, GBR, FIN, TSI, IBS) fixed, over EUR samples only.
Proportions are model-based estimates against a few reference populations, not statements about origin.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from wgsreport import score
from wgsreport.utils import log, run


def user_genotypes(vcf: Path, work: Path, ids: list[str]) -> pd.DataFrame:
    pvar = pd.read_pickle(work / "pvar.pkl")
    pvar = pvar[pvar["id"].isin(ids)].reset_index(drop=True)
    sites = score.read_user_sites(vcf, set(zip(pvar["chrom"], pvar["pos"].astype(int))))
    dos, _, _ = score.user_dosages(sites, pvar)
    return pvar, dos


def write_user_vcf(pvar, dos, out: Path):
    """Single-sample VCF at the panel variants (ref/alt as in the panel); NaN dosage -> missing."""
    gts = {0: "0/0", 1: "0/1", 2: "1/1"}
    with open(f"{out}.vcf", "w") as f:
        f.write("##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        for c in sorted(set(pvar["chrom"]), key=str):
            f.write(f"##contig=<ID={c}>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tUSER\n")
        for d, r in zip(dos, pvar.itertuples()):
            g = "./." if np.isnan(d) else gts[int(round(d))]
            f.write(f"{r.chrom}\t{int(r.pos)}\t{r.id}\t{r.ref}\t{r.alt}\t.\t.\t.\tGT\t{g}\n")
    run(f"bgzip -f {out}.vcf && tabix -p vcf {out}.vcf.gz")


def run_admixture(binary: Path, bed: Path, pop_labels: list[str], k: int, out: Path, threads: int):
    Path(f"{bed}.pop").write_text("\n".join(pop_labels) + "\n")
    subprocess.run([str(binary), "--supervised", f"-j{threads}", f"{bed}.bed", str(k)], cwd=out, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, type=Path)
    ap.add_argument("--work-dir", required=True, type=Path)
    ap.add_argument("--admixture", required=True, type=Path, help="path to the ADMIXTURE binary")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    work = a.work_dir
    meta = pd.read_csv(work / "panel" / "samples.tsv", sep="\t")
    keep = [l.strip() for l in open(work / "prune.prune.in")]
    log(f"{len(keep):,} LD-pruned variants")

    # reference genotypes at pruned variants
    run(["plink2", "--pfile", str(work / "panel" / "panel-merge"), "--extract", str(work / "prune.prune.in"),
         "--snps-only", "just-acgt", "--max-alleles", "2", "--make-pgen", "--out", str(out / "kg")])
    pvar, dos = user_genotypes(a.vcf, work, keep)
    write_user_vcf(pvar, dos, out / "user")
    run(["plink2", "--pfile", str(out / "kg"), "--export", "vcf", "bgz", "id-paste=iid", "--out", str(out / "kg_export")])
    run(f"tabix -f -p vcf {out}/kg_export.vcf.gz")
    run(f"bcftools merge -m none {out}/kg_export.vcf.gz {out}/user.vcf.gz -Oz -o {out}/all.vcf.gz")
    run(["plink2", "--vcf", str(out / "all.vcf.gz"), "--vcf-half-call", "missing", "--make-bed", "--out", str(out / "all")])
    fam = pd.read_csv(out / "all.fam", sep=r"\s+", header=None)
    iid = fam[1].tolist()
    lab = meta.set_index("IID")
    results = {}
    # --- global
    sup = {"AFR": "AFR", "EUR": "EUR", "EAS": "EAS", "SAS": "SAS"}
    labels = [sup.get(lab.loc[i, "super_pop"], "-") if i in lab.index else "-" for i in iid]
    import shutil
    shutil.copy(out / "all.bed", out / "global.bed"); shutil.copy(out / "all.bim", out / "global.bim")
    shutil.copy(out / "all.fam", out / "global.fam")
    run_admixture(a.admixture.resolve(), out / "global", labels, 4, out, a.threads)
    q = pd.read_csv(out / "global.4.Q", sep=r"\s+", header=None)
    # supervised column order = order of first appearance in .pop file
    order = [x for x in dict.fromkeys(labels) if x != "-"]
    q.columns = order
    results["global"] = q.iloc[iid.index("USER")].round(4).to_dict()
    # --- europe
    eu = [i for i in iid if i == "USER" or (i in lab.index and lab.loc[i, "super_pop"] == "EUR")]
    (out / "eur.keep").write_text("IID\n" + "\n".join(eu) + "\n")
    run(["plink2", "--bfile", str(out / "all"), "--keep", str(out / "eur.keep"), "--make-bed", "--out", str(out / "europe")])
    fam = pd.read_csv(out / "europe.fam", sep=r"\s+", header=None)
    ids = fam[1].tolist()
    labels = [lab.loc[i, "pop"] if i in lab.index else "-" for i in ids]
    run_admixture(a.admixture.resolve(), out / "europe", labels, 5, out, a.threads)
    q = pd.read_csv(out / "europe.5.Q", sep=r"\s+", header=None)
    q.columns = [x for x in dict.fromkeys(labels) if x != "-"]
    results["europe"] = q.iloc[ids.index("USER")].round(4).to_dict()
    import json
    (out / "admixture.json").write_text(json.dumps(results, indent=1))
    print(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
