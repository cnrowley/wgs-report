import numpy as np
import pandas as pd
import pytest

from wgsreport import score, snps, clinvar
from wgsreport.registry import load_registry


def pvar():
    return pd.DataFrame({"chrom": ["1", "1", "2", "2"], "pos": [100, 200, 300, 400],
                         "id": ["1:100:A:G", "1:200:C:T", "2:300:A:T", "2:400:G:C"],
                         "ref": ["A", "C", "A", "G"], "alt": ["G", "T", "T", "C"]})


def sc(rows):
    return pd.DataFrame(rows, columns=["chrom", "pos", "effect", "other", "weight"])


AF = pd.Series({"1:100:A:G": .3, "1:200:C:T": .3, "2:300:A:T": .5, "2:400:G:C": .1})


def test_harmonize_sign_and_flip():
    s = sc([("1", 100, "G", "A", 0.5),       # effect == ALT   -> +0.5
            ("1", 200, "C", "T", 0.2),       # effect == REF   -> -0.2
            ("2", 400, "C", None, 0.3)])     # effect == ALT, other unknown -> +0.3
    w, st = score.harmonize(s, pvar(), AF)
    assert w["1:100:A:G"] == 0.5 and w["1:200:C:T"] == -0.2 and w["2:400:G:C"] == 0.3
    # opposite-strand report of the 1:100 A/G SNP: effect C on the other strand == G
    s2 = sc([("1", 100, "C", "T", 0.5)])
    w2, st2 = score.harmonize(s2, pvar(), AF)
    assert w2["1:100:A:G"] == 0.5 and st2["n_strand_flipped"] == 1


def test_harmonize_drops_mismatch_and_ambiguous_palindromes():
    s = sc([("1", 100, "T", "C", 1.0),       # T/C is the opposite-strand form of A/G -> kept, flipped
            ("1", 200, "A", "C", 1.0),       # cannot match C/T on either strand
            ("2", 300, "T", "A", 1.0)])      # palindromic with MAF 0.5 -> dropped
    w, st = score.harmonize(s, pvar(), AF)
    assert "1:200:C:T" not in w.index and "2:300:A:T" not in w.index
    assert st["n_allele_mismatch"] == 1 and st["n_palindromic_dropped"] == 1


def test_user_dosages_absent_is_homref_and_filters():
    sites = {("1", 100): [("A", ["G"], [0, 1], True)],
             ("1", 200): [("C", ["T"], [1, 1], False)],           # filtered -> NaN
             ("2", 400): [("G", ["C", "T"], [2, 2], True)]}        # third allele only
    d, st, assumed = score.user_dosages(sites, pvar())
    assert assumed.tolist() == [False, False, True, False]
    assert d[0] == 1 and np.isnan(d[1]) and d[2] == 0 and d[3] == 0
    assert st["assumed_hom_ref"] == 1 and st["third_allele"] == 1


def test_user_scores_mean_imputes():
    w = pd.DataFrame({"s": [1.0, 2.0]})
    s = score.user_scores(np.array([np.nan, 1.0]), w, np.array([0.25, 0.5]))
    assert s["s"] == pytest.approx(0.5 * 1.0 + 1.0 * 2.0)


def test_snp_genotype_and_apoe():
    panel = snps.load_panel(include_sensitive=True)
    assert panel["rsid"].is_unique and {"rs429358", "rs7412"} <= set(panel["rsid"])
    sites = {}
    r = panel[panel.rsid == "rs429358"].iloc[0]
    sites[(r.chrom, int(r.pos))] = [(r.ref, r.alt.split(","), [0, 1], True)]
    g = snps.genotype(panel, sites)
    assert g.set_index("rsid").loc["rs429358", "n_effect"] == (1 if r.effect == r.alt.split(",")[0] else 1)
    assert snps.decode_apoe(g) in {"e3/e4", "undetermined", "e3/e3"}


def test_clinvar_screen_zygosity():
    plp = pd.DataFrame([{"chrom": "7", "pos": 10, "ref": "C", "alt": "T", "clinvar_id": "1", "significance": "Pathogenic",
                         "stars": 2, "condition": "x", "gene": "BRCA2", "consequence": "", "af_1000g": None, "rsid": ""}])
    hits = clinvar.screen(plp, {("7", 10): [("C", ["T"], [0, 1], True)]})
    assert len(hits) == 1 and hits.iloc[0]["zygosity"] == "heterozygous" and hits.iloc[0]["acmg_sf_gene"]


def test_registry_valid():
    reg = load_registry()
    assert len(reg["traits"]) >= 10
    for t in reg["traits"]:
        assert t["tier"] in "ABC" and (t["pgs"].startswith("PGS") or t["pgs"].startswith("custom:"))


def test_cyp2d6_parse_and_pileup_counts():
    from wgsreport.pgx import parse_cyp2d6
    from wgsreport.cram import _count_bases, call_from_counts
    assert parse_cyp2d6("*2/*68+*4")["phenotype"] == "intermediate metabolizer"
    assert parse_cyp2d6("*4/*4")["activity_score"] == 0
    c = _count_bases("T", ".,.,.^].G,,.+2AAg$")
    assert c["T"] == 9 and c["G"] == 2
    alleles, _ = call_from_counts({"A": 0, "C": 0, "G": 12, "T": 14}, "T", ["G"])
    assert alleles == ["G", "T"]


def test_clinvar_patterns():
    rows = [{"chrom": "1", "pos": p, "ref": "C", "alt": "T", "clinvar_id": str(p), "significance": "Pathogenic", "stars": 2,
             "condition": "x", "gene": g, "consequence": "", "af_1000g": None, "rsid": ""} for p, g in ((1, "CFTR"), (2, "CFTR"), (3, "GJB2"), (4, "HBB"))]
    plp = pd.DataFrame(rows)
    sites = {("1", 1): [("C", ["T"], [0, 1], True)], ("1", 2): [("C", ["T"], [0, 1], True)],
             ("1", 3): [("C", ["T"], [0, 1], True)], ("1", 4): [("C", ["T"], [1, 1], True)]}
    hits = clinvar.screen(plp, sites).set_index("gene")
    assert hits.loc["GJB2", "pattern"] == "single heterozygous"
    assert hits.loc["HBB", "pattern"] == "homozygous"
    assert set(hits.loc["CFTR", "pattern"]) == {"possible compound heterozygote (phase unknown)"}
