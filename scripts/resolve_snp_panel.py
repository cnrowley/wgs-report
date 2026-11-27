#!/usr/bin/env python3
"""Resolve rsIDs in snp_panel_src.yaml to GRCh38 coordinates/alleles via the Ensembl REST API.

Writes src/wgsreport/data/snp_panel.tsv. Only needed when the curated list changes; the resolved
table is shipped with the package. Uses only public identifiers (no personal data).
"""
import json
import sys
import time
import urllib.request
from pathlib import Path

import yaml

root = Path(__file__).resolve().parents[1] / "src" / "wgsreport" / "data"
src = yaml.safe_load(open(root / "snp_panel_src.yaml"))["snps"]
ids = [s["rsid"] for s in src]
req = urllib.request.Request(
    "https://rest.ensembl.org/variation/homo_sapiens",
    data=json.dumps({"ids": ids}).encode(),
    headers={"Content-Type": "application/json", "Accept": "application/json"},
)
for attempt in range(5):
    try:
        resp = json.load(urllib.request.urlopen(req, timeout=120))
        break
    except Exception as e:  # noqa: BLE001
        print("retry", e, file=sys.stderr)
        time.sleep(5 * (attempt + 1))
else:
    sys.exit("Ensembl request failed")

comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
rows, problems = [], []
for s in src:
    r = resp.get(s["rsid"])
    if not r:
        problems.append(f"{s['rsid']}: not found")
        continue
    maps = [m for m in r["mappings"] if m["assembly_name"] == "GRCh38" and m["seq_region_name"] in
            [str(i) for i in range(1, 23)] + ["X"]]
    if not maps:
        problems.append(f"{s['rsid']}: no primary GRCh38 mapping")
        continue
    m = maps[0]
    alleles = m["allele_string"].split("/")
    ref, alts = alleles[0], alleles[1:]
    if len(ref) != 1 or any(len(a) != 1 for a in alts):
        problems.append(f"{s['rsid']}: not a SNV ({m['allele_string']})")
        continue
    eff = s["effect"]
    if eff == "alt":
        if len(alts) != 1:
            problems.append(f"{s['rsid']}: 'alt' ambiguous with alleles {m['allele_string']}")
            continue
        eff = alts[0]
    if eff not in alleles:
        problems.append(f"{s['rsid']}: effect allele {eff} not in Ensembl alleles {m['allele_string']} "
                        f"(strand mix-up?)")
        continue
    rows.append((s["rsid"], m["seq_region_name"], m["start"], ref, ",".join(alts), eff,
                 s["gene"], s["category"], s["trait"], *s["notes"]))
for p in problems:
    print("PROBLEM:", p, file=sys.stderr)
with open(root / "snp_panel.tsv", "w") as f:
    f.write("rsid\tchrom\tpos\tref\talt\teffect\tgene\tcategory\ttrait\tnote0\tnote1\tnote2\n")
    for r in rows:
        f.write("\t".join(str(x) for x in r) + "\n")
print(f"resolved {len(rows)}/{len(src)} SNPs; {len(problems)} problems")
