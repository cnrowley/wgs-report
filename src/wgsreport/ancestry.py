"""Genetic ancestry estimate: project the user into the 1000 Genomes PCA space."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .utils import log, run

N_PCS = 10


def build_pca(panel: Path, workdir: Path, threads: int = 8) -> dict:
    """LD-prune common panel variants, run PCA on all 1000G samples, return file locations."""
    workdir = Path(workdir)
    pr = workdir / "prune"
    run(["plink2", "--pfile", str(panel), "--maf", "0.05", "--geno", "0.02", "--indep-pairwise",
         "1000", "50", "0.1", "--threads", str(threads), "--out", str(pr)])
    n_in = sum(1 for _ in open(f"{pr}.prune.in"))
    log(f"PCA on {n_in:,} LD-pruned common variants")
    pca = workdir / "pca"
    run(["plink2", "--pfile", str(panel), "--extract", f"{pr}.prune.in", "--freq",
         "--pca", str(N_PCS), "allele-wts", "--threads", str(threads), "--out", str(pca)])
    return {"eigenvec": f"{pca}.eigenvec", "allele": f"{pca}.eigenvec.allele",
            "afreq": f"{pca}.afreq", "eigenval": f"{pca}.eigenval"}


def project_user(dos: np.ndarray, pvar: pd.DataFrame, files: dict) -> np.ndarray:
    """Project the user's dosages onto the panel PCs (variance-standardized allele loadings)."""
    load = pd.read_csv(files["allele"], sep="\t")
    # plink2 writes one row per (variant, allele); A1 is the allele the loading refers to.
    pcs = [c for c in load.columns if c.startswith("PC")]
    freq = pd.read_csv(files["afreq"], sep="\t").set_index("ID")["ALT_FREQS"]
    idx = pd.Series(np.arange(len(pvar)), index=pvar["id"])
    load = load[load["ID"].isin(idx.index)]
    ii = idx.reindex(load["ID"]).to_numpy()
    d = dos[ii].astype(np.float64)
    p = freq.reindex(load["ID"]).to_numpy()
    a1_is_alt = (load["A1"].to_numpy() == pvar["alt"].to_numpy()[ii])
    d = np.where(a1_is_alt, d, 2.0 - d)
    p = np.where(a1_is_alt, p, 1.0 - p)
    d = np.where(np.isnan(d), 2.0 * p, d)
    z = (d - 2.0 * p) / np.sqrt(2.0 * p * (1.0 - p))
    return z @ load[pcs].to_numpy(dtype=np.float64)


def calibrate_scale(panel_pcs: pd.DataFrame, panel: Path, files: dict, workdir: Path,
                    threads: int = 8) -> np.ndarray:
    """Per-PC scale between our projection formula and plink2's sample PCs.

    plink2's eigenvec are unit-normalized sample vectors, so applying the allele loadings to a
    standardized genotype vector gives PCs up to a per-PC constant. We estimate that constant by
    projecting the panel samples with plink2 itself (--score on the same loadings).
    """
    workdir = Path(workdir)
    load = pd.read_csv(files["allele"], sep="\t")
    pcs = [c for c in load.columns if c.startswith("PC")]
    wfile = workdir / "pca_loadings.tsv"
    load[["ID", "A1"] + pcs].to_csv(wfile, sep="\t", index=False)
    out = workdir / "pca_proj"
    run(["plink2", "--pfile", str(panel), "--read-freq", files["afreq"], "--score", str(wfile),
         "1", "2", "header-read", "cols=scoresums", "variance-standardize",
         "--score-col-nums", f"3-{2 + len(pcs)}", "--threads", str(threads), "--out", str(out)])
    proj = pd.read_csv(f"{out}.sscore", sep="\t")
    proj = proj.rename(columns={proj.columns[0]: "IID"}).set_index("IID")
    proj = proj[[c for c in proj.columns if c.endswith("_SUM")]]
    proj.columns = pcs
    ev = panel_pcs.set_index("IID")[pcs].reindex(proj.index)
    scale = np.array([np.linalg.lstsq(proj[c].to_numpy()[:, None], ev[c].to_numpy(), rcond=None)[0][0]
                      for c in pcs])
    r2 = [np.corrcoef(proj[c], ev[c])[0, 1] ** 2 for c in pcs]
    if min(r2) < 0.999:
        raise RuntimeError(f"PCA projection calibration failed (min R2={min(r2):.4f})")
    return scale


def classify(user_pcs: np.ndarray, panel_pcs: pd.DataFrame, meta: pd.DataFrame,
             k: int = 25, use_pcs: int = 6) -> dict:
    """k-nearest-neighbour ancestry call from the first PCs; also nearest sub-populations."""
    pcs = [f"PC{i}" for i in range(1, use_pcs + 1)]
    df = panel_pcs.merge(meta, on="IID")
    X = df[pcs].to_numpy()
    dist = np.sqrt(((X - user_pcs[:use_pcs]) ** 2).sum(1))
    nn = np.argsort(dist)[:k]
    sp = df["super_pop"].iloc[nn].value_counts(normalize=True)
    pops = df["pop"].iloc[nn].value_counts()
    best = sp.index[0]
    sd = df.groupby("super_pop")[pcs].std().loc[best].to_numpy()
    centroid = df[df.super_pop == best][pcs].mean().to_numpy()
    z_dist = float(np.sqrt((((user_pcs[:use_pcs] - centroid) / sd) ** 2).sum()))
    return {"superpop": best, "superpop_fractions": sp.to_dict(), "nearest_pops": pops.to_dict(),
            "distance_from_centroid_sd": z_dist, "min_distance": float(dist[nn[0]])}
