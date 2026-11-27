#!/usr/bin/env python3
"""Rebuild src/wgsreport/data/cancer_predisposition_clingen.tsv from ClinGen's gene-disease validity curations.

Keeps genes whose curation by one of the four ClinGen cancer expert panels is Definitive or Strong. Run it to refresh the list;
review the diff before committing. Needs network access to ClinGen only (no personal data is sent).
"""
import csv
import io
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgsreport import sources  # noqa: E402

PANELS = ("Hereditary Cancer", "Breast/Ovarian Cancer", "Colon Cancer", "Childhood, Adolescent and Young Adult Cancer Predisposition")
OUT = Path(__file__).resolve().parents[1] / "src" / "wgsreport" / "data" / "cancer_predisposition_clingen.tsv"

raw = urllib.request.urlopen(sources.CLINGEN_VALIDITY_URL, timeout=120).read().decode("utf-8")
rows = list(csv.reader(io.StringIO(raw)))
created = next((r[0] for r in rows if r and r[0].startswith("FILE CREATED")), "FILE CREATED: unknown").split(": ")[-1]
h = next(i for i, r in enumerate(rows) if r and r[0] == "GENE SYMBOL")
cols = rows[h]
best = {}
for r in rows[h + 2:]:
    if len(r) < len(cols):
        continue
    d = dict(zip(cols, r))
    if d["CLASSIFICATION"] in ("Definitive", "Strong") and any(p in d["GCEP"] for p in PANELS):
        g = best.setdefault(d["GENE SYMBOL"], {"moi": set(), "dis": []})
        g["moi"].add(d["MOI"])
        entry = f"{d['DISEASE LABEL']} ({d['MOI']}, {d['CLASSIFICATION']})"
        if entry not in g["dis"]:
            g["dis"].append(entry)
lines = ["# Cancer-predisposition genes with Definitive or Strong ClinGen gene-disease validity from the ClinGen cancer expert panels.",
         f"# Source: ClinGen gene-disease validity curations, file created {created} ({sources.CLINGEN_VALIDITY_URL}); ClinGen data are CC0.",
         "# Regenerate with scripts/build_cancer_genes.py. A gene listed here is a research-grade flag, not a clinical test panel.",
         "# moi lists every inheritance mode curated for the gene (a gene can be recessive for one disease and dominant for cancer risk).",
         "gene\tmoi\tdiseases"]
lines += [f"{g}\t{'/'.join(sorted(v['moi']))}\t{'; '.join(v['dis'])}" for g, v in sorted(best.items())]
OUT.write_text("\n".join(lines) + "\n")
print(f"{len(best)} genes -> {OUT}")
