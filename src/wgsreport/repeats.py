"""Repeat-expansion screen with ExpansionHunter (Illumina, Apache-2.0; Dolzhenko et al., Genome Res 2020, doi:10.1101/gr.265927.120)
and pathogenic-size thresholds from STRchive (data/repeat_thresholds.tsv).

ExpansionHunter is not bundled: download a release (sources.EXPANSIONHUNTER_URL) and pass its directory. It genotypes about 30 known
loci from the CRAM. It sizes repeats accurately up to roughly the read length (150 bp, about 50 CAG units); longer expansions show up as
a lower bound. A call in the intermediate or pathogenic range needs confirmation in a clinical laboratory, most of all near a threshold.
Reads whose mate maps to a contig missing from the FASTA cannot be paired (use a FASTA with the alt contigs for the best sensitivity).
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

from .registry import data_path
from .utils import log

_COMP = str.maketrans("ACGTN", "TGCAN")


def canonical(motif: str) -> str:
    """Rotation- and strand-independent form of a repeat unit (so AAG, GAA, CTT and TTC all match)."""
    m = motif.upper()
    forms = []
    for s in (m, m.translate(_COMP)[::-1]):
        forms += [s[i:] + s[:i] for i in range(len(s))]
    return min(forms)


def load_thresholds() -> dict[str, dict]:
    out = {}
    for l in data_path("repeat_thresholds.tsv").read_text().splitlines():
        if not l or l.startswith("#") or l.startswith("locus\t"):
            continue
        f = l.split("\t")
        num = lambda x: int(x) if x else None
        out[f[0]] = {"gene": f[1], "disease": f[2], "inheritance": f[3], "motifs": [m for m in f[4].split(",") if m],
                     "benign_max": num(f[5]), "intermediate_min": num(f[6]), "intermediate_max": num(f[7]), "pathogenic_min": num(f[8]), "id": f[9]}
    return out


def run_expansionhunter(eh_dir: Path, cram: Path, fasta: Path, out_dir: Path, sex: str = "male", threads: int = 4) -> Path:
    eh_dir = Path(eh_dir).resolve()
    exe = next(iter(eh_dir.glob("**/bin/ExpansionHunter")), None)
    cat = next(iter(eh_dir.glob("**/variant_catalog/hg38/variant_catalog.json")), None)
    if not exe or not cat:
        raise SystemExit(f"ExpansionHunter binary or hg38 variant catalog not found under {eh_dir}")
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    log("repeat expansions: running ExpansionHunter")
    subprocess.run([str(exe), "--reads", str(Path(cram).resolve()), "--reference", str(Path(fasta).resolve()), "--variant-catalog", str(cat),
                    "--output-prefix", str(out_dir / "eh"), "--sex", sex, "--threads", str(threads)], check=True)
    return out_dir / "eh.json"


def _alleles(genotype: str) -> list[int]:
    return [int(x) for x in genotype.replace(">", "").split("/") if x.strip().isdigit()]


def category(n: int, th: dict) -> str:
    if th["pathogenic_min"] is not None and n >= th["pathogenic_min"]:
        return "pathogenic range"
    if th["intermediate_min"] is not None and n >= th["intermediate_min"]:
        return "intermediate range"
    if th["benign_max"] is not None and n <= th["benign_max"]:
        return "normal range"
    return "between normal and intermediate"


def parse(eh_json: Path) -> dict:
    """ExpansionHunter JSON -> per-locus rows with a category for each allele; loci whose counted unit is not the pathogenic motif are not interpreted."""
    d = json.loads(Path(eh_json).read_text())
    ths = load_thresholds()
    rows = []
    for lid, loc in sorted(d.get("LocusResults", {}).items()):
        th = ths.get(lid.upper()) or ths.get(lid)
        read_len = loc.get("ReadLength") or 150
        for vid, v in (loc.get("Variants") or {}).items():
            if v.get("VariantType") != "Repeat" or not v.get("Genotype"):
                continue
            unit = v.get("RepeatUnit") or ""
            counts = _alleles(v["Genotype"])
            row = {"locus": lid, "variant": vid, "gene": th["gene"] if th else lid, "disease": th["disease"] if th else "", "unit": unit,
                   "genotype": v["Genotype"], "ci": v.get("GenotypeConfidenceInterval", ""), "alleles": counts,
                   "inheritance": th["inheritance"] if th else "", "interpreted": False, "categories": [], "worst": "", "beyond_read_length": False}
            if th and counts and any(canonical(unit) == canonical(m) for m in th["motifs"]):
                row["interpreted"] = True
                row["categories"] = [category(n, th) for n in counts]
                order = ["normal range", "between normal and intermediate", "intermediate range", "pathogenic range"]
                row["worst"] = max(row["categories"], key=order.index)
                row["beyond_read_length"] = max(counts) * len(unit) > 0.8 * read_len
                row["thresholds"] = {k: th[k] for k in ("benign_max", "intermediate_min", "intermediate_max", "pathogenic_min")}
            rows.append(row)
    return {"loci": rows, "n_interpreted": sum(r["interpreted"] for r in rows),
            "flagged": [r for r in rows if r["worst"] in ("intermediate range", "pathogenic range")]}
