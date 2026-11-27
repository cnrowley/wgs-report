#!/usr/bin/env python3
"""Print held-out accuracy by ancestry for candidate PGS Catalog scores (public metadata only; no personal data).

usage: python scripts/portability_table.py PGS001889 PGS002102 ...
R2 is the reported R2, or the square of a reported partial correlation. Use it to fill `r2_by_ancestry` / `alt_scores` in traits.yaml."""
import json
import sys
import urllib.request

API = "https://www.pgscatalog.org/rest"


def get(u):
    return json.load(urllib.request.urlopen(u, timeout=60))


for pid in sys.argv[1:]:
    s = get(f"{API}/score/{pid}")
    print(f"{pid}  {s['variants_number']} variants  {s['publication']['firstauthor']} {s['publication']['date_publication'][:4]}  {s['trait_reported']}")
    for r in get(f"{API}/performance/search?pgs_id={pid}&limit=100")["results"]:
        ss = r["sampleset"]["samples"][0]
        m = {x["name_short"]: x["estimate"] for k in ("effect_sizes", "class_acc", "othermetrics") for x in r["performance_metrics"].get(k) or []}
        r2 = m.get("R²") or m.get("partial-R2") or (m["Partial Correlation (partial-r)"] ** 2 if "Partial Correlation (partial-r)" in m else None)
        print(f"   {str(ss.get('ancestry_broad'))[:34]:34} n={ss.get('sample_number')}  R2={'' if r2 is None else round(r2, 3)}")
