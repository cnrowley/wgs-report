"""Consistency checks on the bundled data files (no network). A failure here means an edit broke a data contract."""
import re

import pandas as pd

from wgsreport import sources
from wgsreport.clinvar import acmg_genes
from wgsreport.registry import data_path, load_registry, score_ids
from wgsreport.rare import load_af_genes

TIERS = {"A", "B", "C"}


def test_traits_registry_contract():
    reg = load_registry()
    ids = [t["id"] for t in reg["traits"]]
    assert len(ids) == len(set(ids)), "duplicate trait ids"
    for t in reg["traits"]:
        for k in ("id", "label", "pgs", "kind", "tier", "cite"):
            assert k in t, f"{t.get('id')}: missing {k}"
        assert t["tier"] in TIERS and t["kind"] in ("quantitative", "binary")
        assert re.fullmatch(r"PGS\d{6}|custom:\w+", t["pgs"])
        if t["kind"] == "quantitative":
            assert 0 < t["r2"] < 1, f"{t['id']}: r2"
        else:
            assert t["or_per_sd"] > 1, f"{t['id']}: or_per_sd"
        if t.get("crosscheck"):
            assert re.fullmatch(r"PGS\d{6}", t["crosscheck"])
        if t["tier"] != "C" and t["kind"] == "quantitative" and t.get("unit"):
            assert set(t["mean"]) == set(t["sd"]), f"{t['id']}: mean/sd sexes differ"
    assert len(score_ids(reg)) == len(set(score_ids(reg)))


def test_related_rsids_exist_in_snp_panel():
    panel = set(pd.read_csv(data_path("snp_panel.tsv"), sep="\t")["rsid"])
    for t in load_registry()["traits"]:
        assert set(t.get("related_rsids", [])) <= panel, t["id"]


def test_acmg_table_is_complete_and_consistent():
    rows = [l.split("\t") for l in data_path("acmg_sf_v3.3.tsv").read_text().splitlines()
            if l and not l.startswith("#") and not l.startswith("gene\t")]
    genes = [r[0] for r in rows]
    assert len(genes) == len(set(genes)) == 84, "ACMG SF v3.3 has 84 genes; update the file and this number together"
    assert all(re.fullmatch(r"[A-Z0-9-]+", g) for g in genes)
    assert {r[2] for r in rows} <= {"AD", "AR", "SD", "XL"}
    assert {"BRCA1", "LDLR", "ABCD1", "CYP27A1", "PLN"} <= acmg_genes()


def test_af_gene_set():
    af = load_af_genes()
    assert len(af) == 21 and af["TTN"] > 1


def test_snp_panel_contract():
    df = pd.read_csv(data_path("snp_panel.tsv"), sep="\t", dtype={"chrom": str})
    assert df["rsid"].is_unique
    for r in df.itertuples():
        alleles = {r.ref, *r.alt.split(",")}
        assert r.effect in alleles, f"{r.rsid}: effect allele {r.effect} not among {alleles}"
        assert r.category in {"trait", "nutrition", "health", "pgx", "sensitive"}
        assert all(isinstance(x, str) and x for x in (r.note0, r.note1, r.note2)), r.rsid


def test_sources_are_templates():
    assert "{id}" in sources.PGS_FTP and "{md5}" in sources.EBI_MD5_URL
    assert str(sources.ENSEMBL_RELEASE) in sources.ENSEMBL_GFF_URL


def test_snp_panel_text_quality():
    df = pd.read_csv(data_path("snp_panel.tsv"), sep="\t", dtype=str)
    for r in df.itertuples():
        for col in ("trait", "note0", "note1", "note2"):
            v = getattr(r, col)
            assert v.count("(") == v.count(")"), f"{r.rsid} {col}: unbalanced parentheses: {v}"
    assert "very unlikely" not in " ".join(df["note0"]), "celiac proxy wording must not claim exclusion"


def test_carrier_panel_complete():
    from wgsreport.clinvar import carrier_genes
    g = carrier_genes()
    assert len(g) == 113                      # Gregg 2021 Tier 3: 97 autosomal recessive + 16 X-linked
    assert {"CFTR", "HBB", "SMN1", "BLM", "DMD", "FMR1"} <= set(g)


def test_cancer_gene_list():
    from wgsreport.clinvar import cancer_genes
    g = cancer_genes()
    assert len(g) >= 100 and {"BRCA2", "CHEK2", "ATM", "MUTYH", "TP53"} <= set(g)
    assert "AD" in g["ATM"] and "AR" in g["ATM"]           # both modes kept (cancer risk is dominant, ataxia telangiectasia recessive)


def test_ancestry_selection():
    from wgsreport.registry import for_ancestry, load_registry, score_ids, tier_from_r2
    reg = load_registry()
    lpa = next(t for t in reg["traits"] if t["id"] == "lpa")
    assert for_ancestry(lpa, "EUR", True)["pgs"] == lpa["pgs"]
    sas = for_ancestry(lpa, "SAS", True, {"PGS000689", lpa["pgs"]})
    assert sas["pgs"] == "PGS000689" and sas["tier"] == "B" and "switched" in sas["portability"]
    # alt score not in the panel -> falls back to the same score with its non-European accuracy, tier recomputed
    assert for_ancestry(lpa, "SAS", True, {lpa["pgs"]})["pgs"] == lpa["pgs"]
    apob = next(t for t in reg["traits"] if t["id"] == "apob")
    assert for_ancestry(apob, "SAS", True)["tier"] == "C"            # R2 0.075 < 0.09: percentile only
    assert "No published accuracy" in for_ancestry(apob, "AMR", True)["portability"]
    assert for_ancestry(apob, "SAS", False)["tier"] == "B"           # unclear ancestry: defaults, report warns separately
    assert tier_from_r2(0.25) == "A" and tier_from_r2(0.1) == "B" and tier_from_r2(0.05) == "C"
    assert {"PGS000689", "PGS001888"} <= set(score_ids(reg))


def test_related_effect_and_text_keys_are_consistent():
    reg = load_registry()
    for t in reg["traits"]:
        rel = set(t.get("related_rsids", []))
        assert set(t.get("related_effect", {})) <= rel, t["id"]
        assert set(t.get("related_text", {})) <= rel, t["id"]
        assert set(t.get("related_effect", {}).values()) <= {"up", "down"}, t["id"]
        if t.get("related_not_in_score"):
            assert rel, t["id"]
