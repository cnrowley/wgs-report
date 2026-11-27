import json
"""Parsers for the optional external tools, on small synthetic fixtures (shapes taken from the real tool output)."""
import pytest
from wgsreport import hla_typing, pgx_pharmcat, report


def _allele(name):
    return {"gene": "CYP2C19", "name": name, "function": "x", "reference": False, "activityValue": "n/a"}


def _pharmcat_report():
    d = {"label": "*1/*2", "allele1": _allele("*1"), "allele2": _allele("*2"), "phenotypes": ["Intermediate Metabolizer"], "activityScore": None}
    n = {"label": "*1/*1", "allele1": _allele("*1"), "allele2": _allele("*1"), "phenotypes": ["Normal Metabolizer"], "activityScore": None}
    ann = lambda gene, rec: {"drugRecommendation": rec, "classification": "Strong", "implications": ["x"], "lookupKey": [{gene: "p"}]}
    return {"pharmcatVersion": "9.9.9", "dataVersion": "d",
            "genes": {"CYP2C19": {"recommendationDiplotypes": [d], "callSource": "MATCHER", "uncalledHaplotypes": ["*16"], "hasUndocumentedVariations": False},
                      "TPMT": {"recommendationDiplotypes": [n], "callSource": "MATCHER"},
                      "DPYD": {"recommendationDiplotypes": [], "callSource": "NONE"}},
            "drugs": {"CPIC Guideline Annotation": {
                "clopidogrel": {"guidelines": [{"annotations": [ann("CYP2C19", "Avoid standard dose clopidogrel")]}]},
                "azathioprine": {"guidelines": [{"annotations": [ann("TPMT", "Use standard dosing")]}]}}}}


def test_pharmcat_summary_hides_routine():
    s = pgx_pharmcat.summarize(_pharmcat_report())
    g = {x["gene"]: x for x in s["genes"]}
    assert set(g) == {"CYP2C19", "TPMT"}                       # empty DPYD call dropped
    assert g["CYP2C19"]["actionable"] and not g["TPMT"]["actionable"]
    assert [d["drug"] for d in s["drugs"]] == ["clopidogrel"] and s["routine_recommendations_hidden"] == 1


def test_t1k_parse_and_interpret(tmp_path):
    f = tmp_path / "g.tsv"
    f.write_text("B\t2\tB*57:01:01\t50\t60\tB*07:02:01\t40\t55\t\n"
                 "A\t1\tA*02:01:01\t90\t70\t.\t0\t-1\t\n"
                 "DQA1\t2\tDQA1*05:01:01\t30\t40\tDQA1*03:01:01\t20\t35\t\n"
                 "DQB1\t2\tDQB1*02:01:01\t30\t40\tDQB1*03:02:01\t20\t35\t\n"
                 "C\t1\t.\t0\t-1\t.\t0\t-1\t\n")
    calls = hla_typing.parse_genotype(f)
    by = {c["gene"]: c for c in calls}
    assert by["A"]["zygosity"].startswith("homozygous") and by["C"]["zygosity"] == "not typed"
    out = hla_typing.interpret(calls)
    assert [x["allele"] for x in out["flagged"]] == ["B*57:01"]
    assert out["celiac_compatible"] == {"DQ2.5": True, "DQ2.2": False, "DQ8": True}


def test_t1k_tied_candidates(tmp_path):
    f = tmp_path / "g.tsv"
    f.write_text("DQB1\t2\tDQB1*02:02:01,DQB1*02:02:15,DQB1*03:03:02\t30\t40\tDQB1*03:02:01\t20\t35\t\n"
                 "DQA1\t1\tDQA1*02:01:01\t58\t58\t.\t0\t-1\t\n"
                 "B\t2\tB*07:02:01,B*57:01:01\t5\t12\tB*44:02:01\t5\t12\t\n")
    calls = {c["gene"]: c for c in hla_typing.parse_genotype(f)}
    assert calls["DQB1"]["tied"] and len(calls["DQB1"]["slots"][0]) == 3
    out = hla_typing.interpret(list(calls.values()))
    assert out["celiac_compatible"] == {"DQ2.5": False, "DQ2.2": True, "DQ8": False}      # DQ8 needs DQA1*03:01, absent here
    fl = out["flagged"][0]
    assert fl["allele"] == "B*57:01" and fl["tied_with_other_alleles"]


def test_t1k_low_quality_ignored(tmp_path):
    f = tmp_path / "g.tsv"
    f.write_text("B\t1\tB*15:02:01\t5\t-1\t.\t0\t-1\t\n")
    c = hla_typing.parse_genotype(f)[0]
    assert c["alleles"] == [] and not hla_typing.interpret([c])["flagged"]


def test_report_sections_render():
    pc = pgx_pharmcat.summarize(_pharmcat_report())
    html = "".join(report._pharmcat_html(pc))
    assert "Intermediate Metabolizer" in html and "clopidogrel" in html and "standard dosing" in html
    h = {"calls": [{"gene": "B", "alleles": ["B*57:01:01"], "zygosity": "heterozygous", "quality": 60, "low_confidence": False}],
         "flagged": [{"allele": "B*57:01", "copies": 1, "meaning": "Abacavir", "source": "CPIC"}], "celiac_typed": False,
         "celiac_compatible": {}, "any_low_confidence": False}
    assert "Abacavir" in "".join(report._hla_html(h))
    assert report._pharmcat_html(None) == [] and report._hla_html(None) == []


def test_pharmcat_ambiguous_and_reference():
    rep = _pharmcat_report()
    mk = lambda lab, ph: {"label": lab, "phenotypes": [ph], "activityScore": None}
    rep["genes"]["CYP2C19"]["recommendationDiplotypes"] = [mk("*1/*4", "Intermediate Metabolizer"), mk("*1/*17", "Rapid Metabolizer")]
    rep["genes"]["RYR1"] = {"recommendationDiplotypes": [mk("Reference/Reference", "Uncertain Susceptibility")], "callSource": "MATCHER"}
    g = {x["gene"]: x for x in pgx_pharmcat.summarize(rep)["genes"]}
    assert g["CYP2C19"]["ambiguous"] and g["CYP2C19"]["diplotype"].startswith("ambiguous:") and g["CYP2C19"]["activity_score"] is None
    assert not g["RYR1"]["actionable"]


def test_cyp2d6_outside_normalised():
    assert pgx_pharmcat.normalize_outside({"CYP2D6": "*2/*68+*10"}) == {"CYP2D6": "*2/*10"}
    assert pgx_pharmcat.normalize_outside(None) is None


def test_genome_qc_classification_and_roh():
    from wgsreport import genome_qc as q
    auto = {f"chr{i}": 1.0 for i in range(1, 23)}
    assert q.classify({**auto, "chrX": 0.5, "chrY": 0.48})[0] == "XY"
    assert q.classify({**auto, "chrX": 1.0, "chrY": 0.01})[0] == "XX"
    assert q.classify({**auto, "chrX": 1.0, "chrY": 0.5})[0].startswith("XXY")
    assert q.classify({**auto, "chrX": 0.5, "chrY": 0.0})[0].startswith("X0")
    call, note, flags = q.classify({**auto, "chr21": 1.5, "chrX": 0.5, "chrY": 0.5})
    assert flags and "chr21" in flags[0]
    txt = ("# comment\nRG\ts\tchr1\t100\t3000100\t3000001\t900\t45\nRG\ts\tchr2\t100\t900100\t900001\t90\t60\n"
           "RG\ts\tchr3\t100\t4000100\t4000001\t500\t10\n")
    r = q.parse_roh(txt, auto_len=1_000_000_000)
    assert r["n_runs"] == 1 and r["total_bp"] == 3_000_001 and abs(r["f_roh"] - 0.003) < 1e-4    # short and low-quality runs excluded


def test_qc_report_section():
    from wgsreport import report
    qc = {"dosage": {"ratios": {"chrX": 0.5, "chrY": 0.5, "chr1": 1.0, "chr2": 1.01}, "sex_chromosome_call": "XY", "note": "n", "autosome_flags": []},
          "x_het": {"n_calls": 100, "het_fraction": 0.01},
          "roh": {"n_runs": 1, "total_bp": 3_000_000, "f_roh": 0.001, "longest_bp": 3_000_000, "min_run_bp": 1_500_000,
                  "context": [{"relationship": "first cousins", "expected_f": 0.0625}]}}
    html = "".join(report._qc_html(qc))
    assert "XY" in html and "none" in html and "first cousins" in html
    assert report._qc_html(None) == []


def _eh_json():
    v = lambda unit, gt: {"VariantType": "Repeat", "RepeatUnit": unit, "Genotype": gt, "GenotypeConfidenceInterval": "x-x/y-y"}
    return {"LocusResults": {
        "HTT": {"ReadLength": 150, "Variants": {"HTT": v("CAG", "19/30"), "HTT_CCG": v("CCG", "7/10")}},
        "FXN": {"ReadLength": 150, "Variants": {"FXN": v("GAA", "11/20")}},
        "DMPK": {"ReadLength": 150, "Variants": {"DMPK": v("CTG", "12/60")}},          # reverse-complement unit, pathogenic size
        "RFC1": {"ReadLength": 150, "Variants": {"RFC1": v("AARRG", "11/11")}}}}


def test_repeat_parse(tmp_path):
    from wgsreport import repeats
    p = tmp_path / "eh.json"
    p.write_text(__import__("json").dumps(_eh_json()))
    r = repeats.parse(p)
    by = {(x["locus"], x["variant"]): x for x in r["loci"]}
    assert by[("HTT", "HTT")]["worst"] == "intermediate range" and not by[("HTT", "HTT_CCG")]["interpreted"]
    assert by[("FXN", "FXN")]["worst"] == "normal range"
    assert by[("DMPK", "DMPK")]["worst"] == "pathogenic range"                        # CTG matches the CAG motif on the other strand
    assert not by[("RFC1", "RFC1")]["interpreted"]                                       # unit differs from the pathogenic motif
    assert {x["gene"] for x in r["flagged"]} == {"HTT", "DMPK"}
    assert repeats.canonical("AAG") == repeats.canonical("CTT") == repeats.canonical("GAA")


def test_repeat_report_section(tmp_path):
    from wgsreport import repeats, report
    p = tmp_path / "eh.json"
    p.write_text(__import__("json").dumps(_eh_json()))
    html = "".join(report._repeats_html(repeats.parse(p)))
    assert "intermediate range" in html and "paternal" in html and "lower bound" in html and "pathogenic range" in html and "Normal range: FXN" in html
    assert report._repeats_html(None) == []


def test_yleaf_parse(tmp_path):
    from wgsreport import y_haplogroup
    f = tmp_path / "hg_prediction.hg"
    f.write_text("Sample_name\tHg\tHg_marker\tTotal_reads\tValid_markers\tQC-score\tQC-1\tQC-2\tQC-3\nS1\tR-X1\tX1\t100\t1500\t0.97\t1\t1\t1\n")
    assert y_haplogroup.parse(f) == {"haplogroup": "R-X1", "marker": "X1", "valid_markers": 1500, "qc": 0.97}


def test_lof_scan_dedupe_and_mechanism_filter():
    import pandas as pd
    from wgsreport import rare
    lof = pd.DataFrame([{"chrom": "chr15", "pos": 10, "ref": "C", "alt": "CA", "gene": "BLM"},
                        {"chrom": "chr2", "pos": 5, "ref": "G", "alt": "T", "gene": "CALM2"}])
    hits = pd.DataFrame([{"chrom": "15", "pos": 10, "ref": "C", "alt": "CA"}])
    assert list(rare.drop_clinvar_matches(lof, hits)["gene"]) == ["CALM2"]
    assert "SCG5" in rare.NO_LOF_MECHANISM and "BRCA2" not in rare.NO_LOF_MECHANISM


def test_repeat_note_is_gene_specific():
    from wgsreport import report
    assert "paternal" in report._repeat_note("HTT", "intermediate range")
    n = report._repeat_note("ATXN1", "intermediate range")
    assert "paternal" not in n and "Interruptions" in n


def test_variant_set_additive_and_filtered():
    import pandas as pd
    from wgsreport import variant_scores
    df = pd.DataFrame({"chrom": ["6"] * 4, "pos": [10, 20, 30, 40], "effect": ["G", "T", "C", "A"], "other": ["A", "C", "T", "G"],
                       "weight": [100.0, 10.0, 5.0, 7.0], "rsid": ["rsA", "rsB", "rsC", "rsD"]})
    sites = {("6", 10): [("A", ["G"], [0, 1], True)],      # heterozygote: one copy
             ("6", 20): [("C", ["T"], [1, 1], True)],      # homozygote: two copies
             ("6", 30): [("T", ["C"], [0, 1], False)]}     # fails filter: not counted, reported as unusable
    r = variant_scores.score_variant_set(df, sites)       # ("6", 40) absent: assumed reference
    assert r["total"] == 120.0 and [c["rsid"] for c in r["carried"]] == ["rsA", "rsB"]
    assert r["carried"][0]["copies"] == 1 and r["carried"][1]["copies"] == 2
    assert r["unusable"] == 1 and r["assumed_reference"] == 1


def test_variant_reading_headline():
    r = {"label": "Lipoprotein(a)", "variant_set_primary": True, "variant_set_threshold": 125, "variant_set": {"unit": "nmol/L", "total": 130.0}}
    reading, value = report.variant_reading(r, [])
    assert "predict raised" in reading and "not a total concentration" in reading and "130" in value
    assert report.variant_reading({**r, "variant_set": {"unit": "nmol/L", "total": 0.0}}, ["rs1"])[1] == "raised"   # related variant alone is enough
    assert report.variant_reading({**r, "variant_set": {"unit": "nmol/L", "total": 0.0}}, []) is None               # nothing carried: polygenic score applies
    assert report.variant_reading({**r, "variant_set_primary": False}, ["rs1"]) is None


def test_variant_set_html_shows_reference_rank():
    v = {"pgs": "PGSX", "unit": "nmol/L", "n_variants": 5, "total": 100.0, "assumed_reference": 4, "unusable": 0,
         "carried": [{"rsid": "rs1", "chrom": "6", "pos": 1, "effect_allele": "T", "copies": 1, "weight": 100.0, "contribution": 100.0}],
         "reference": {"group": "EUR", "n": 503, "fraction_at_or_above": 0.02, "fraction_any_carried": 0.3}}
    html = report._variant_set_html(v)
    assert "2% have a total this high" in html and "30% carry any" in html


def test_pharmcat_format_change_and_version_warning():
    import pytest
    with pytest.raises(ValueError, match="format has probably changed"):
        pgx_pharmcat.summarize({"something": "else"})
    s = pgx_pharmcat.summarize({**_pharmcat_report(), "pharmcatVersion": "9.9.9"})
    assert "not tested" in s["version_warning"]
    assert pgx_pharmcat.summarize({**_pharmcat_report(), "pharmcatVersion": "3.4.0"})["version_warning"] is None


def test_pharmcat_install_check(monkeypatch, tmp_path, capsys):
    from wgsreport import pgx_pharmcat as ph
    monkeypatch.setattr(ph.shutil, "which", lambda x: None)
    with pytest.raises(SystemExit):
        ph.check_install(tmp_path / "pharmcat.jar")
    monkeypatch.setattr(ph.shutil, "which", lambda x: "/usr/bin/java")
    monkeypatch.setattr(ph, "installed_version", lambda jar: "9.9.9")
    assert ph.check_install(tmp_path / "pharmcat.jar") == "9.9.9"        # untested version warns, does not block


def test_provenance_has_no_personal_fields(tmp_path):
    from wgsreport import provenance
    pv = provenance.collect(["PGS000001"], tmp_path, tmp_path)
    assert {"wgsreport", "python", "packages", "tools"} <= set(pv)
    assert not any(k in pv for k in ("vcf", "sample", "path"))



def test_cancer_section_separates_dominant_recessive_and_low_review():
    import pandas as pd
    hits = pd.DataFrame([
        {"gene": "BRCA2", "cancer_gene": "AD", "rsid": "rs1", "zygosity": "heterozygous", "pattern": "single heterozygous", "condition": "Breast cancer",
         "significance": "Pathogenic", "stars": 3, "consequence": "frameshift_variant", "af_1000g": None},
        {"gene": "BLM", "cancer_gene": "AR", "rsid": "rs2", "zygosity": "heterozygous", "pattern": "single heterozygous", "condition": "Bloom syndrome",
         "significance": "Pathogenic", "stars": 2, "consequence": "frameshift_variant", "af_1000g": 0.0002}])
    low = [{**hits.iloc[0].to_dict(), "gene": "PALB2", "stars": 1}]
    html = "".join(report._cancer_html({"cancer_low_review": low, "traits": []}, hits))
    assert "BRCA2 (AD)" in html and "recessive-only" in html and "BLM (AR)" in html
    assert "Lower review status" in html and "PALB2" in html
    assert "No absolute lifetime risk is calculated" in html and "Not assessed" in html


def test_hoxb13_is_a_cancer_gene():
    from wgsreport import clinvar
    assert "HOXB13" in clinvar.cancer_genes() and "BRCA1" in clinvar.cancer_genes()


def test_pgx_consequences():
    from wgsreport import pgx
    assert pgx.consequences_cyp2c19(2, 0)[0]["action"] == pgx.ALT
    assert pgx.consequences_cyp2c19(0, 0)[0]["action"] == pgx.USE
    assert all(c["action"] == pgx.DOSE for c in pgx.consequences_dpyd(1))
    assert pgx.consequences_dpyd(2)[0]["action"] == pgx.ALT
    assert pgx.consequences_g6pd(1)[0]["action"] == pgx.ALT
    assert pgx.consequences_cyp2d6(0.0)[0]["action"] == pgx.ALT
    from wgsreport import report
    html = report._pgx_consequence_html([{"gene": "CYP2C19", "consequences": pgx.consequences_cyp2c19(2, 0)}])
    assert "avoid / choose alternative" in html and "clopidogrel" in html


def test_medication_report():
    from wgsreport import pgx, medication_report
    res = {"pgx": [{"gene": "CYP2C19", "result": "*2 x2", "phenotype": "poor", "test_type": "limited screen", "consequences": pgx.consequences_cyp2c19(2, 0)},
                   {"gene": "TPMT", "result": "0", "phenotype": "none", "test_type": "limited screen", "consequences": pgx.consequences_thiopurine(0, "TPMT")}],
           "provenance": {"wgsreport": "1.0", "date": "2026-01-01", "tools": {"bcftools": "1.20"}}}
    h = medication_report.build(res)
    assert "How this report was produced" in h and "clopidogrel" in h and "Not a clinical test" in h
    assert "RESEARCH USE ONLY" in h and "prescriber" not in h.lower() and "pharmacist" not in h.lower() and "doctor" not in h.lower()
    assert "SAMPLE" not in h


def test_cyp3a5_phenotype_matches_copies_of_star3():
    import pandas as pd
    from wgsreport import pgx
    for n, word in ((0, "expresser (*1/*1)"), (1, "intermediate"), (2, "non-expresser")):
        g = pd.DataFrame([{"rsid": "rs776746", "n_effect": n, "status": "called"}])
        row = [r for r in pgx.summarize(g) if r["gene"].startswith("CYP3A5")][0]
        assert word in row["phenotype"]


def test_clinical_tests_map():
    import pandas as pd
    from wgsreport import clinical_tests
    geno = pd.DataFrame([{"rsid": "rs6025", "gene": "F5", "genotype": "C/C", "n_effect": 0, "status": "called"}])
    cv = pd.DataFrame([{"gene": "LDLR", "rsid": "rs1", "pos": 5, "zygosity": "heterozygous"}])
    rows = {r["test"]: r for r in clinical_tests.summarize({}, geno, cv)}
    fh = rows["Familial hypercholesterolemia panel"]
    assert "LDLR" in fh["result"]
    assert "no ClinVar" in rows["Cardiomyopathy panel (HCM, DCM, ARVC)"]["result"]
    assert "F5" in rows["Thrombophilia (factor V Leiden, prothrombin G20210A)"]["result"]
    assert any(r["coverage"] == "none" for r in rows.values())
    assert "not run" in rows["Spinal muscular atrophy (SMN1 copy number)"]["result"]


def test_cnv_screen_classification():
    from wgsreport import cnv_screen as n
    assert n.classify([2.0] * 100)["call"] == "not detected"
    assert n.classify([1.0] * 100)["call"] == "deletion"
    assert n.classify([3.0] * 100)["call"] == "duplication"
    assert n.classify([2.0] * 40 + [1.0] * 40 + [2.0] * 20)["call"] == "partial deletion"
    assert n.classify([2.0, 1.0] * 50)["call"] == "not detected"          # scattered noisy bins are not a segment
    assert n.classify([1.0] * 5)["call"] == "not assessable"
    assert n.classify([2.0] * 40, n_total=100)["call"] == "not assessable"   # mostly missing data is never "not detected"
    # a heterozygous deletion keeps the unique-mapping fraction, so the filter must not drop it; multi-mapping windows go
    assert n.usable_copies([20, 40, 40], [20, 40, 10], 40.0) == [1.0, 2.0]


def test_cnv_list_states_guidance_strength_and_follow_up():
    from wgsreport import cnv_screen as n
    names = {r["name"] for r in n.ACTIONABLE}
    for r in n.ACTIONABLE:
        assert r["change"] in ("del", "dup") and r["guidance"] and r["follow_up"] and r["source"], r["name"]
        assert r["end"] - r["start"] >= 200_000
    for d in n.DISEASES:
        assert d["region"] in names, d
        assert d["evidence"].startswith(("independent replication", "UK Biobank only; known effect"))


def _screen(calls):
    from wgsreport import cnv_screen as n
    res = {n.key(*r): {"call": "not detected", "median_copies": 2.0} for r in n.regions()}
    for name, call in calls.items():
        spec = next(r for r in n.ACTIONABLE if name in r["name"])
        res[n.key(spec["chrom"], spec["start"], spec["end"])] = {"call": call, "median_copies": 1.0}
    return {"baseline_depth": 30.0, "regions": res, "x_control_copies": 1.0}


def test_cnv_screen_flags_only_the_actionable_direction():
    from wgsreport import cnv_screen as n
    from wgsreport.report import _cnv_screen_html
    flagged = lambda calls: [f for f in n.findings(_screen(calls)) if f["flag"]]
    assert flagged({"22q11.2": "duplication"}) == []                       # only the deletion is on the list
    assert [f["change"] for f in flagged({"PMP22), deletion": "duplication"})] == ["dup"]
    assert flagged({"17q12": "partial deletion"})[0]["partial"]
    page = _cnv_screen_html(_screen({"PMP22), deletion": "duplication"}))
    assert "Vincristine" in page and "OR 21.8" in page and "leg-crossing" not in page and "case reports and reviews only" in page
    assert "None of the 6" in _cnv_screen_html(_screen({}))
    assert "Not assessable" in _cnv_screen_html(_screen({"NF1": "not assessable"}))


def test_exclude_psychiatric_keeps_cnv_screen_without_psychiatric_estimates(tmp_path):
    from wgsreport.analyze import _load_copynumber
    from wgsreport.report import _cnv_screen_html
    p = tmp_path / "cn.json"
    p.write_text(json.dumps({"loci": {}, "cnv_screen": _screen({"22q11.2": "deletion"})}))
    page = _cnv_screen_html(_load_copynumber(p, False)["cnv_screen"])
    assert "vitamin D" in page and "schizophrenia" not in page.lower()
    assert "schizophrenia" in _cnv_screen_html(_load_copynumber(p, True)["cnv_screen"])


def _hit(gene, pattern="single heterozygous", chrom="7", sig="Pathogenic", af=0.02):
    return {"gene": gene, "pattern": pattern, "chrom": chrom, "significance": sig, "condition": "synthetic condition", "stars": 2,
            "rsid": "rs0", "af_1000g": af}


def test_carrier_findings_classify_by_inheritance_and_ignore_frequency():
    import pandas as pd
    from wgsreport import carrier
    hits = pd.DataFrame([_hit("CFTR"), _hit("DMD", chrom="X"), _hit("PAH", pattern="homozygous"),
                         _hit("HBB", sig="Pathogenic|risk_factor"), _hit("NOTAGENE")])
    f = carrier.findings(hits, sex="male")
    by = {r["gene"]: r for r in f["carrier"]}
    assert by["CFTR"]["status"] == "carrier (one copy)"                      # a 2% allele is not filtered out
    assert "affected genotype in males" in by["DMD"]["status"]
    assert by["PAH"]["status"].startswith("two copies")
    assert [r["gene"] for r in f["risk_allele"]] == ["HBB"]
    assert "raise your own risk" in f["risk_allele"][0]["status"]
    assert [h["gene"] for h in f["uncurated"]] == ["NOTAGENE"]
    assert "can be passed" in {r["gene"]: r for r in carrier.findings(hits, sex="female")["carrier"]}["DMD"]["status"]


def test_carrier_section_always_states_what_it_cannot_see():
    import pandas as pd
    from wgsreport import carrier
    from wgsreport.report import _carrier_html
    page = "".join(_carrier_html({"sex": "female"}, pd.DataFrame()))
    for gene in ("SMN1", "HBA1", "FMR1", "DMD", "CYP21A2", "F8"):
        assert gene in page
    assert "does not rule out being a carrier" in page and "Repeat expansions were not run" in page
    cov = {r["gene"]: r["result"] for r in carrier.panel_coverage(pd.DataFrame())}
    assert len(cov) == 113
    assert cov["SMN1"].startswith("not assessed") and cov["HBA2"].startswith("not assessed")
    assert "deletions and duplications not assessed" in cov["ASPA"]


def test_no_result_derived_numbers_are_hardcoded_in_report_text():
    from pathlib import Path
    src = Path(__file__).parents[1] / "src" / "wgsreport"
    assert "KIV-2 copies in total" not in (src / "report.py").read_text()
    assert "An estimate of" not in (src / "copynumber.py").read_text()
