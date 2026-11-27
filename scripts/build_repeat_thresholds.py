#!/usr/bin/env python3
"""Build src/wgsreport/data/repeat_thresholds.tsv: pathogenic repeat-size thresholds for the loci in ExpansionHunter's GRCh38 catalog.

Thresholds come from STRchive (Dashnow lab, MIT license; https://strchive.org). Usage:
    python scripts/build_repeat_thresholds.py path/to/ExpansionHunter/variant_catalog/hg38/variant_catalog.json
Review the diff before committing.
"""
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from wgsreport import sources  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "src" / "wgsreport" / "data" / "repeat_thresholds.tsv"
catalog = json.load(open(sys.argv[1]))
loci = json.loads(urllib.request.urlopen(sources.STRCHIVE_LOCI_URL, timeout=120).read())
lines = ["# Pathogenic repeat-size thresholds (repeat units) for ExpansionHunter catalog loci. Source: STRchive (https://strchive.org, MIT license), "
         "data/STRchive-loci.json.", "# Regenerate with scripts/build_repeat_thresholds.py. 'intermediate' covers STRchive's intermediate range (reduced-penetrance, "
         "premutation or gray-zone alleles).",
         "locus\tgene\tdisease\tinheritance\tmotif\tbenign_max\tintermediate_min\tintermediate_max\tpathogenic_min\tstrchive_id"]
for c in sorted(catalog, key=lambda x: x["LocusId"]):
    reg = c["ReferenceRegion"][0] if isinstance(c["ReferenceRegion"], list) else c["ReferenceRegion"]
    chrom, rng = reg.split(":")
    a, b = map(int, rng.split("-"))
    hit = [y for y in loci if (y.get("gene") or "").upper() == c["LocusId"].upper() and y.get("chrom") == chrom
           and y.get("start_hg38") and y["start_hg38"] <= b + 50 and y["stop_hg38"] >= a - 50]
    if not hit:
        continue
    y = hit[0]
    motif = ",".join(y.get("pathogenic_motif_reference_orientation") or y.get("reference_motif_reference_orientation") or [])
    f = lambda k: "" if y.get(k) is None else str(y[k])
    lines.append("\t".join([c["LocusId"], y["gene"], y.get("disease_id") or "", "/".join(y.get("inheritance") or []), motif,
                            f("benign_max"), f("intermediate_min"), f("intermediate_max"), f("pathogenic_min"), y["id"]]))
OUT.write_text("\n".join(lines) + "\n")
print(f"{len(lines) - 4} loci -> {OUT}")
