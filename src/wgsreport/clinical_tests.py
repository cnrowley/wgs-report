"""Common clinical genetic tests that whole-genome data can address, with this genome's result for each (see data/clinical_tests.tsv).

Gene-group rows screen ClinVar P/LP hits already computed (exact SNV/indel matches), SNP rows read the curated panel, module rows point to an optional module
or say it was not run. It is a map from familiar test names to what this report did, not a substitute for the clinical test."""
from __future__ import annotations

import pandas as pd

from .registry import data_path


def load() -> pd.DataFrame:
    return pd.read_csv(data_path("clinical_tests.tsv"), sep="\t", comment="#", dtype=str, keep_default_na=False,
                       names=["area", "test", "kind", "items", "coverage", "limits"], header=None)


def _genes(items: str, cv) -> str:
    gs = set(items.split(";"))
    if cv is None or not len(cv):
        return f"no ClinVar pathogenic/likely-pathogenic variant found in {len(gs)} genes"
    h = cv[cv["gene"].isin(gs)]
    if not len(h):
        return f"no ClinVar pathogenic/likely-pathogenic variant found in {len(gs)} genes"
    return "; ".join(f"{r.gene} {r.rsid or r.pos} ({r.zygosity})" for r in h.itertuples())


def _snps(items: str, geno) -> str:
    out = []
    for rs in items.split(";"):
        g = geno[geno["rsid"] == rs]
        if not len(g):
            out.append(f"{rs}: not assessed")
            continue
        r = g.iloc[0]
        out.append(f"{r['gene']} {rs} {r['genotype']} ({r['n_effect']} risk/effect copies; {str(r['status']).split(' (')[0]})")
    return "; ".join(out)


def _module(items: str, res: dict) -> str:
    if items.startswith("copynumber:"):
        key = items.split(":", 1)[1]
        loci = ((res.get("copynumber") or {}).get("loci")) or {}
        hit = {k: v for k, v in loci.items() if k.upper().startswith(key.upper())}
        if not hit:
            return "not run (needs a CRAM and `wgsreport copy-number`)"
        return "; ".join(f"{k}: about {v['estimated_total_copies']:.1f} copies" for k, v in hit.items())
    if items == "cnvscreen":
        nc = (res.get("copynumber") or {}).get("cnv_screen")
        if not nc or nc.get("error"):
            return "not run (needs a CRAM and `wgsreport copy-number`)"
        from .cnv_screen import findings
        rows = findings(nc)
        hit = [f"{f['name']}: {f['call']}" for f in rows if f["flag"]]
        unread = sum(f["call"] == "not assessable" for f in rows)
        base = "; ".join(hit) if hit else f"none of {len(rows) - unread} actionable changes detected"
        return base + (f" ({unread} not assessable)" if unread else "")
    if items == "lpa":
        return "see the lipoprotein(a) result in the trait section"
    if res.get(items):
        return "run: see the matching section above"
    return "not run (optional module; see docs/USAGE.md)"


def summarize(res: dict, geno, cv_hits) -> list[dict]:
    rows = []
    for r in load().itertuples():
        if r.kind == "genes":
            result = _genes(r.items, cv_hits)
        elif r.kind == "snps":
            result = _snps(r.items, geno)
        elif r.kind == "module":
            result = _module(r.items, res)
        else:
            result = "not covered"
        rows.append({"area": r.area, "test": r.test, "coverage": r.coverage or "none", "result": result, "limits": r.limits})
    return rows
