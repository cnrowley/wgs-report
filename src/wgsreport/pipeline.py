"""High-level pipeline stages: `prepare` (one-off, no personal data) and `analyze`."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import ancestry, pgs, refpanel, score
from .registry import load_registry, score_ids
from .utils import log, require_tools


MIN_SCORE_COVERAGE = 0.2


def check_coverage(pid: str, st: dict) -> None:
    """Fail loudly when a score matches almost none of the panel (typical of a genome-build or coordinate mismatch)."""
    if st["n_variants"] and st["coverage"] < MIN_SCORE_COVERAGE:
        raise SystemExit(f"{pid}: only {100 * st['coverage']:.1f}% of its variants matched the GRCh38 reference panel. "
                         "This usually means the weights are on a different genome build or use unexpected coordinates; "
                         "refusing to continue.")


def prepare(kg_dir: Path, work_dir: Path, jobs: int = 6, threads: int = 8,
            kg_pattern: str = refpanel.DEFAULT_PATTERN) -> None:
    """Download PGS weights, build the reference panel, harmonize weights, PCA, score the panel.

    Contains no personal data; the output directory can be reused for any number of individuals.
    """
    require_tools("bcftools", "plink2")
    work = Path(work_dir)
    (work / "pgs_cache").mkdir(parents=True, exist_ok=True)
    pan_dir = work / "panel"
    pan_dir.mkdir(exist_ok=True)

    reg = load_registry()
    ids = score_ids(reg)
    log(f"fetching {len(ids)} PGS Catalog scoring files (GRCh38-harmonized)")
    scores, headers = {}, {}
    for pid in ids:
        if pid.startswith("custom:"):
            f = work / "custom_scores" / f"{pid.split(':', 1)[1]}_pt.tsv"
            if not f.exists():
                log(f"  {pid}: not built yet (see `wgsreport build-score`); skipping")
                continue
            scores[pid] = pd.read_csv(f, sep="\t", dtype={"chrom": str})
            headers[pid] = {"trait_reported": "locally built from summary statistics"}
        else:
            scores[pid], headers[pid] = pgs.read_score(pgs.fetch_score(pid, work / "pgs_cache"))
        log(f"  {pid}: {len(scores[pid]):,} usable autosomal SNVs ({headers[pid].get('trait_reported', '?')})")

    universe = (pd.concat([s[["chrom", "pos"]] for s in scores.values()])
                .drop_duplicates().reset_index(drop=True))
    panel = refpanel.build_panel(universe, kg_dir, pan_dir, pattern=kg_pattern, jobs=jobs, threads=threads)

    pvar = refpanel.read_pvar(panel)
    af_eur = refpanel.read_afreq(pan_dir, "EUR")
    af_all = refpanel.read_afreq(pan_dir, "ALL")
    log(f"panel: {len(pvar):,} biallelic SNVs")

    cols, stats = {}, {}
    for pid, s in scores.items():
        w, st = score.harmonize(s, pvar, af_eur)
        cols[pid], stats[pid] = w, st
        check_coverage(pid, st)
        log(f"  {pid}: matched {st['n_scored']:,}/{st['n_variants']:,} variants "
            f"({100 * st['coverage']:.1f}%), flipped {st['n_strand_flipped']:,}, "
            f"mismatch {st['n_allele_mismatch']:,}, palindromic-dropped {st['n_palindromic_dropped']:,}")
    weights = pd.DataFrame(cols).reindex(pvar["id"]).fillna(0.0)
    weights = weights.loc[(weights != 0).any(axis=1)]
    weights.to_pickle(work / "weights.pkl")
    pvar.to_pickle(work / "pvar.pkl")
    af_all.to_pickle(work / "af_all.pkl")
    json.dump({"stats": stats, "headers": headers}, open(work / "score_stats.json", "w"), indent=1)

    log("scoring the 1000 Genomes panel")
    ps = score.score_panel(panel, weights, pvar.set_index("id")["alt"], work, threads=threads)
    ps.to_csv(work / "panel_scores.tsv", sep="\t")

    log("ancestry PCA")
    files = ancestry.build_pca(panel, work, threads=threads)
    pcs = pd.read_csv(files["eigenvec"], sep="\t")
    pcs = pcs.rename(columns={pcs.columns[0]: "IID"})
    pcs = pcs[["IID"] + [c for c in pcs.columns if c.startswith("PC")]]
    scale = ancestry.calibrate_scale(pcs, panel, files, work, threads=threads)
    np.save(work / "pc_scale.npy", scale)
    pcs.to_csv(work / "panel_pcs.tsv", sep="\t", index=False)
    json.dump(files, open(work / "pca_files.json", "w"))
    log("prepare complete")
