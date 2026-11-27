#!/usr/bin/env python3
"""Does a genome-wide polygenic score already capture a known large-effect locus, or should an explicit variant override it?

For each (trait, candidate variant) the script uses only public data: the PGS weights, the 1000 Genomes EUR genotypes and Ensembl coordinates.
It reports (a) whether the candidate is in the score and its weight, (b) the EUR allele frequency, (c) r^2 between the candidate's dosage and the
score computed from the OTHER score variants within +-WINDOW bp. Low r^2 means the score does not track the candidate through LD, so the candidate carries
information the score misses (an override candidate if its effect is large; effect sizes must come from the literature, this script does not estimate them).
usage: locus_override_check.py WORK_DIR KG_DIR EUR_SAMPLES OUT.tsv"""
import gzip
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
import yaml

work, kg, eur, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
WINDOW = 500_000
CAND = {
    "bilirubin": [("rs887829", "UGT1A1"), ("rs6742078", "UGT1A1"), ("rs4148323", "UGT1A1*6")],
    "urate": [("rs2231142", "ABCG2 Q141K"), ("rs12498742", "SLC2A9"), ("rs734553", "SLC2A9")],
    "ldl": [("rs11591147", "PCSK9 R46L"), ("rs429358", "APOE e4"), ("rs7412", "APOE e2"), ("rs5742904", "APOB R3527Q"), ("rs12740374", "SORT1")],
    "apob": [("rs11591147", "PCSK9 R46L"), ("rs429358", "APOE e4"), ("rs7412", "APOE e2")],
    "hdl": [("rs3764261", "CETP"), ("rs247616", "CETP"), ("rs1532085", "LIPC")],
    "tg": [("rs964184", "APOA5/ZPR1"), ("rs662799", "APOA5"), ("rs328", "LPL S447X"), ("rs1260326", "GCKR")],
    "apoa1": [("rs964184", "APOA5/ZPR1"), ("rs3764261", "CETP")],
    "bmi": [("rs9939609", "FTO"), ("rs17782313", "MC4R")],
    "bmd": [("rs3736228", "LRP5 A1330V"), ("rs4355801", "TNFRSF11B")],
    "height": [("rs143384", "GDF5"), ("rs1042725", "HMGA2")],
    "lpa": [("rs10455872", "LPA"), ("rs3798220", "LPA I4399M")],
    "cad": [("rs10455872", "LPA"), ("rs3798220", "LPA"), ("rs1333049", "9p21"), ("rs11591147", "PCSK9 R46L"), ("rs12740374", "SORT1")],
    "t2d": [("rs7903146", "TCF7L2"), ("rs13266634", "SLC30A8 R325W"), ("rs1801282", "PPARG P12A"), ("rs5219", "KCNJ11 E23K")],
    "af": [("rs2200733", "4q25/PITX2"), ("rs10033464", "4q25/PITX2"), ("rs6817105", "4q25/PITX2")],
    "prostate": [("rs138213197", "HOXB13 G84E"), ("rs1447295", "8q24"), ("rs16901979", "8q24")],
    "melanoma": [("rs1805007", "MC1R R151C"), ("rs1805008", "MC1R R160W"), ("rs1805009", "MC1R D294H"), ("rs1805006", "MC1R D84E"), ("rs1110400", "MC1R R163Q")],
    "hba1c": [("rs334", "HBB sickle"), ("rs1050828", "G6PD A-")],
    "hb": [("rs334", "HBB sickle"), ("rs855791", "TMPRSS6"), ("rs1800562", "HFE C282Y"), ("rs1050828", "G6PD A-")],
    "hct": [("rs334", "HBB sickle"), ("rs855791", "TMPRSS6"), ("rs1050828", "G6PD A-")],
    "rbc": [("rs334", "HBB sickle"), ("rs855791", "TMPRSS6")],
    "mcv": [("rs334", "HBB sickle"), ("rs855791", "TMPRSS6"), ("rs1800562", "HFE C282Y")],
    "mch": [("rs334", "HBB sickle"), ("rs855791", "TMPRSS6")],
    "wbc": [("rs2814778", "DARC/Duffy null")],
    "neut": [("rs2814778", "DARC/Duffy null")],
    "crp": [("rs1205", "CRP"), ("rs2794520", "CRP")],
    "shbg": [("rs1799941", "SHBG"), ("rs727428", "SHBG")],
    "cystc": [("rs911119", "CST3")],
    "alp": [("rs657152", "ABO")],
}
tr = {t["id"]: t for t in yaml.safe_load(open(Path(__file__).resolve().parents[1] / "src/wgsreport/data/traits.yaml"))["traits"]}
ids = sorted({r for v in CAND.values() for r, _ in v})
coords = {}
for i in range(0, len(ids), 40):
    req = urllib.request.Request("https://rest.ensembl.org/variation/homo_sapiens", data=json.dumps({"ids": ids[i:i + 40]}).encode(),
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
    resp = json.load(urllib.request.urlopen(req, timeout=120))
    for r, v in resp.items():
        m = [x for x in v["mappings"] if x["assembly_name"] == "GRCh38" and x["seq_region_name"] in [str(c) for c in range(1, 23)] + ["X"]]
        if m:
            coords[r] = (m[0]["seq_region_name"], int(m[0]["start"]), m[0]["allele_string"].split("/"))
samples = [l.strip() for l in open(eur) if l.strip()]
tmp = Path(out).with_suffix(".samples")
tmp.write_text("\n".join(samples))


def read_score(pid):
    rows = []
    with gzip.open(work / "pgs_cache" / f"{pid}_hmPOS_GRCh38.txt.gz", "rt") as f:
        cols = None
        for l in f:
            if l.startswith("#"):
                continue
            p = l.rstrip("\n").split("\t")
            if cols is None:
                cols = p
                continue
            d = dict(zip(cols, p))
            if d.get("hm_chr") and d.get("hm_pos") and d.get("effect_weight"):
                rows.append((d["hm_chr"], int(d["hm_pos"]), d["effect_allele"], float(d["effect_weight"])))
    return rows


def genos(chrom, lo, hi):
    """positions -> {(pos, ref, alt): dosage array of ALT} for biallelic SNPs, EUR samples."""
    vcf = kg / f"ALL.chr{chrom}_GRCh38.genotypes.20170504.vcf.gz"
    cmd = ["bcftools", "query", "-S", str(tmp), "-r", f"{chrom}:{lo}-{hi}", "-f", "%POS\t%REF\t%ALT[\t%GT]\n", str(vcf)]
    res = subprocess.run(cmd, capture_output=True, text=True).stdout.splitlines()
    d = {}
    for l in res:
        p = l.split("\t")
        if len(p[2]) == 1 and len(p[1]) == 1:
            d[(int(p[0]), p[1], p[2])] = np.array([sum(1 for a in g.replace("|", "/").split("/") if a == "1") for g in p[3:]], float)
    return d


rows = []
for trait, cands in CAND.items():
    pid = tr[trait]["pgs"]
    sc = read_score(pid)
    by_pos = {}
    for c, p, e, w in sc:
        by_pos.setdefault((c, p), []).append((e, w))
    for rsid, gene in cands:
        if rsid not in coords:
            rows.append((trait, pid, rsid, gene, "", "", "", "", "no GRCh38 position")); continue
        chrom, pos, al = coords[rsid]
        g = genos(chrom, max(1, pos - WINDOW), pos + WINDOW)
        gidx = {}
        for (gp, gr, ga), dos in g.items():
            gidx.setdefault(gp, []).append(((gr, ga), dos))
        cand = [(k, v) for k, v in g.items() if k[0] == pos]
        if not cand:
            rows.append((trait, pid, rsid, gene, chrom, pos, "", "", "not a biallelic SNP in 1000G")); continue
        (cp, cref, calt), cdos = cand[0]
        af = cdos.mean() / 2
        win = [(p, e, w) for (c, p, e, w) in sc if c == chrom and abs(p - pos) <= WINDOW]
        inscore = by_pos.get((chrom, pos))
        w_c = ""
        if inscore:
            e, w = inscore[0]
            w_c = w if e == calt else -w if e == cref else ""
        loc = np.zeros(len(samples))
        n = 0
        for p, e, w in win:
            if p == pos:
                continue
            for (gr, ga), dos in gidx.get(p, []):
                if e in (gr, ga):
                    loc += w * (dos if e == ga else 2 - dos)
                    n += 1
                    break
        r2 = float(np.corrcoef(loc, cdos)[0, 1] ** 2) if n > 2 and cdos.std() > 0 and loc.std() > 0 else float("nan")
        rows.append((trait, pid, rsid, gene, chrom, pos, f"{af:.3f}", w_c if w_c == "" else f"{w_c:.4g}", f"{r2:.3f}", n))
hdr = "trait\tpgs\trsid\tgene\tchrom\tpos\teur_af\tweight_alt_allele\tr2_with_rest_of_locus_score\tn_score_variants_in_window\tnote"
with open(out, "w") as f:
    f.write(hdr + "\n")
    for r in rows:
        r = list(r) + [""] * (11 - len(r))
        f.write("\t".join(str(x) for x in r) + "\n")
tmp.unlink()
print(f"wrote {out} ({len(rows)} rows)")
