#!/usr/bin/env python3
"""Rebuild src/wgsreport/data/gene_inheritance.tsv from ClinGen's gene-disease validity curations (all expert panels).

Keeps genes with at least one Definitive, Strong or Moderate recessive (AR), X-linked (XL) curation; notes whether the same gene also has a
dominant (AD/SD) curation. Needs network access to ClinGen only, or pass a downloaded CSV as the first argument.
"""
import csv
import io
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgsreport import sources  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "src" / "wgsreport" / "data" / "gene_inheritance.tsv"
KEEP = ("Definitive", "Strong", "Moderate")

raw = Path(sys.argv[1]).read_text() if len(sys.argv) > 1 else urllib.request.urlopen(sources.CLINGEN_VALIDITY_URL, timeout=120).read().decode("utf-8")
rows = list(csv.reader(io.StringIO(raw)))
created = next((r[0] for r in rows if r and r[0].startswith("FILE CREATED")), "FILE CREATED: unknown").split(": ")[-1]
h = next(i for i, r in enumerate(rows) if r and r[0] == "GENE SYMBOL")
genes = {}
for r in rows[h + 2:]:
    if len(r) < len(rows[h]):
        continue
    d = dict(zip(rows[h], r))
    if d["CLASSIFICATION"] not in KEEP:
        continue
    g = genes.setdefault(d["GENE SYMBOL"], {"rec": set(), "dis": [], "dom": False})
    if d["MOI"] in ("AR", "XL"):
        g["rec"].add(d["MOI"])
        if d["DISEASE LABEL"] not in g["dis"]:
            g["dis"].append(d["DISEASE LABEL"])
    elif d["MOI"] in ("AD", "SD"):
        g["dom"] = True
keep = {k: v for k, v in genes.items() if v["rec"]}
lines = ["# Genes with a Definitive, Strong or Moderate ClinGen gene-disease validity curation inherited as autosomal recessive (AR) or X-linked (XL).",
         f"# Source: ClinGen gene-disease validity curations, file created {created} ({sources.CLINGEN_VALIDITY_URL}); ClinGen data are CC0.",
         "# Regenerate with scripts/build_gene_inheritance.py. A gene that is absent has not been curated as recessive; that is not evidence that it is not.",
         "gene\tmoi\talso_dominant\tdiseases"]
lines += [f"{g}\t{'/'.join(sorted(v['rec']))}\t{'yes' if v['dom'] else ''}\t{'; '.join(sorted(v['dis']))}" for g, v in sorted(keep.items())]
OUT.write_text("\n".join(lines) + "\n")
print(f"{len(keep)} genes -> {OUT}")
