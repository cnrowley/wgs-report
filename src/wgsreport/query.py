"""Ad-hoc queries against one genome: `wgsreport query snp | pgs | copy-number`.

Everything runs locally against your VCF/CRAM. Network use is limited to public identifiers: an rsID sent to Ensembl to learn its
GRCh38 coordinates (skip with --offline and --coords), or a PGS Catalog id to download a public scoring file. No genotype, file name
or sample identifier is ever sent.
"""
from __future__ import annotations

import json
import re
import statistics
import urllib.request
from pathlib import Path

import pandas as pd

from . import cram as cram_mod, copynumber, pgs, variant_scores
from .utils import require_tools, run

ENSEMBL = "https://rest.ensembl.org/variation/homo_sapiens/{rsid}?content-type=application/json"
POS_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|M|MT):(\d+)$", re.I)


def _norm_chrom(c: str) -> str:
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def resolve(queries: list[str], coords: Path | None = None, offline: bool = False) -> list[dict]:
    """rsIDs or chrom:pos (GRCh38) -> [{query, rsid, chrom, pos, ref, alt}]. Order: --coords file, bundled panel, Ensembl REST."""
    known: dict[str, dict] = {}
    if coords:
        d = pd.read_csv(coords, sep=None, engine="python", dtype=str)
        d.columns = [c.lower().lstrip("#") for c in d.columns]
        for r in d.itertuples():
            known[r.rsid] = {"chrom": _norm_chrom(r.chrom), "pos": int(r.pos), "ref": getattr(r, "ref", None), "alt": getattr(r, "alt", None)}
    from .snps import load_panel
    for r in load_panel(True).itertuples():
        known.setdefault(r.rsid, {"chrom": str(r.chrom), "pos": int(r.pos), "ref": r.ref, "alt": r.alt})
    out = []
    for q in queries:
        m = POS_RE.match(q)
        if m:
            out.append({"query": q, "rsid": "", "chrom": m.group(1).upper(), "pos": int(m.group(2)), "ref": None, "alt": None})
            continue
        if q not in known:
            if offline:
                raise SystemExit(f"{q}: coordinates unknown offline. Give chrom:pos (GRCh38) or a --coords file (rsid, chrom, pos).")
            known[q] = _ensembl(q)
        out.append({"query": q, "rsid": q, **known[q]})
    return out


def _ensembl(rsid: str) -> dict:
    if not re.fullmatch(r"rs\d+", rsid):
        raise SystemExit(f"'{rsid}' is neither an rsID (rs123) nor chrom:pos")
    req = urllib.request.Request(ENSEMBL.format(rsid=rsid), headers={"User-Agent": "wgsreport/0.1"})
    d = json.load(urllib.request.urlopen(req, timeout=60))
    maps = [m for m in d.get("mappings", []) if m.get("assembly_name") == "GRCh38" and m["seq_region_name"] in [str(i) for i in range(1, 23)] + ["X", "Y"]]
    if not maps:
        raise SystemExit(f"{rsid}: no primary GRCh38 mapping in Ensembl")
    m = maps[0]
    al = m["allele_string"].split("/")
    return {"chrom": m["seq_region_name"], "pos": int(m["start"]), "ref": al[0], "alt": ",".join(al[1:])}


def vcf_records(vcf: Path, chrom: str, pos: int) -> list[dict]:
    """Records overlapping a position, from an indexed VCF (tries chr1 and 1 naming)."""
    for name in (f"chr{chrom}", chrom):
        try:
            out = run(["bcftools", "query", "-r", f"{name}:{pos}-{pos}", "-f", "%POS\t%REF\t%ALT\t%FILTER\t[%GT]\n", str(vcf)]).stdout
        except RuntimeError as e:
            if "not found in the header" in str(e) or "unknown" in str(e).lower():
                continue
            raise
        recs = []
        for line in out.splitlines():
            p, ref, alt, flt, gt = line.split("\t")
            if int(p) != pos and len(ref) == 1:
                continue
            alleles = [ref] + alt.split(",")
            gts = [None if g in (".", "") else int(g) for g in re.split(r"[/|]", gt)]
            recs.append({"pos": int(p), "ref": ref, "alt": alt, "filter": flt, "gt": gt,
                         "genotype": "/".join("." if g is None else alleles[g] for g in gts)})
        if out or name == chrom:
            return recs
    return []


def query_snps(vcf: Path, queries: list[str], coords: Path | None = None, offline: bool = False,
               cram: Path | None = None, fasta: Path | None = None, sex: str = "male") -> list[dict]:
    require_tools("bcftools")
    rows = []
    for t in resolve(queries, coords, offline):
        recs = vcf_records(vcf, t["chrom"], t["pos"])
        row = {**t, "records": recs}
        if recs:
            row["status"] = "called" + ("" if all(r["filter"] in ("PASS", ".") for r in recs) else " (FILTER not PASS)")
        else:
            row["status"] = "no record: homozygous reference assumed (a variant-only VCF cannot prove it; use --cram to check)"
        if cram and fasta:
            panel = pd.DataFrame([{"rsid": t["query"], "chrom": t["chrom"], "pos": t["pos"], "alt": t["alt"] or ""}])
            row["reads"] = cram_mod.panel_calls(cram, fasta, panel, sex)[t["query"]]
        rows.append(row)
    return rows


def format_snps(rows: list[dict]) -> str:
    out = []
    for r in rows:
        head = f"{r['query']}  GRCh38 {r['chrom']}:{r['pos']}" + (f"  ref {r['ref']} alt {r['alt']}" if r.get("ref") else "")
        out.append(head)
        for x in r["records"]:
            out.append(f"  VCF: {x['ref']}>{x['alt']}  GT {x['gt']} = {x['genotype']}  FILTER {x['filter']}")
        out.append(f"  {r['status']}")
        if "reads" in r:
            rd = r["reads"]
            out.append(f"  reads: {'/'.join(rd['alleles']) if rd['alleles'] else 'no call'} ({rd['note']}); counts {rd['counts']}")
    return "\n".join(out)


def read_scoring_file(path: Path) -> pd.DataFrame:
    """PGS Catalog file (harmonized preferred). A non-harmonized file is read as GRCh38 chr_name/chr_position: the build is your responsibility."""
    try:
        return pgs.read_score(path)[0]
    except (ValueError, OSError):
        raw = pd.read_csv(path, sep="\t", comment="#", dtype=str)
        need = {"chr_name", "chr_position", "effect_allele", "effect_weight"}
        if not need <= set(raw.columns):
            raise SystemExit(f"{path}: need hm_chr/hm_pos (harmonized) or {sorted(need)} columns")
        df = pd.DataFrame({"chrom": raw["chr_name"].map(_norm_chrom), "pos": pd.to_numeric(raw["chr_position"], errors="coerce"),
                           "effect": raw["effect_allele"].str.upper(), "other": raw.get("other_allele"),
                           "weight": pd.to_numeric(raw["effect_weight"], errors="coerce"),
                           "rsid": raw["rsID"] if "rsID" in raw else ""}).dropna(subset=["chrom", "pos", "effect", "weight"])
        df["pos"] = df["pos"].astype(int)
        return df[df["effect"].str.fullmatch("[ACGT]")].reset_index(drop=True)


def query_pgs(vcf: Path, pgs_id: str | None, score_file: Path | None, cache: Path, top: int = 10) -> dict:
    """Raw additive score of one scoring file for this genome (sum of effect-allele copies x weight; absent sites = reference)."""
    from .score import read_user_sites
    require_tools("bcftools")
    if bool(pgs_id) == bool(score_file):
        raise SystemExit("give exactly one of --pgs-id (downloads the public file) or --score-file (offline)")
    path = Path(score_file) if score_file else pgs.fetch_score(pgs_id, cache)
    df = read_scoring_file(path)
    sites = read_user_sites(vcf, variant_scores.positions(df))
    res = variant_scores.score_variant_set(df, sites)
    res["carried"] = res["carried"][:top]
    res.update({"score": pgs_id or Path(score_file).name, "note": "raw weighted sum, not standardised: use `analyze` for percentiles against a reference "
                "panel. Weights are in the score's own units; sites missing from the VCF count as reference."})
    return res


def query_copy_number(cram: Path, fasta: Path, loci: list[str] | None = None, regions: list[str] | None = None) -> dict:
    """Depth-based diploid copy number for named loci (see copynumber.LOCI) and/or custom regions chrom:start-end (GRCh38)."""
    require_tools("samtools")
    fai = Path(str(fasta) + ".fai")
    base = copynumber.baseline_depth(cram, fai)
    out = {"baseline_depth": base, "loci": {}, "note": "total diploid copies from read depth against a median of autosomal windows; near-identical paralogues, "
           "GC and mappability bias give roughly +-1 copy of error (HBA: about +1 on NA12878). Not a diagnostic."}
    wanted = {k: v for k, v in copynumber.LOCI.items() if not loci or any(q.lower() in k.lower() for q in loci)} if (loci or not regions) else {}
    for q in loci or []:
        if not any(q.lower() in k.lower() for k in copynumber.LOCI):
            raise SystemExit(f"unknown locus '{q}'; available: {', '.join(copynumber.LOCI)}")
    for name, spec in wanted.items():
        d = copynumber.bedcov(cram, spec["regions"])
        out["loci"][name] = {"estimated_total_copies": 2 * statistics.mean(d) / base, "depths": d, "expected": spec["expected"], "confidence": spec["confidence"]}
    for rg in regions or []:
        m = re.fullmatch(r"(chr[0-9XYM]+):(\d+)-(\d+)", rg)
        if not m:
            raise SystemExit(f"region '{rg}' must look like chr6:160577000-160651000 (GRCh38)")
        d = copynumber.bedcov(cram, [(m.group(1), int(m.group(2)), int(m.group(3)))])
        out["loci"][rg] = {"estimated_total_copies": 2 * d[0] / base, "depths": d, "expected": "2 for single-copy autosomal sequence", "confidence": "unvalidated custom region"}
    return out
