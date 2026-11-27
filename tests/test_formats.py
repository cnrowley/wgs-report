"""Parsers against small synthetic files that mimic the upstream formats (no network)."""
import gzip
import shutil

import pandas as pd
import pytest

from wgsreport import clinvar, copynumber, pgs, pgx


def _pgs_file(tmp_path, weight_type="beta"):
    p = tmp_path / "PGS000001_hmPOS_GRCh38.txt.gz"
    text = (f"#format_version=2.0\n#weight_type={weight_type}\n#trait_reported=Test\n"
            "chr_name\tchr_position\teffect_allele\tother_allele\teffect_weight\thm_chr\thm_pos\thm_inferOtherAllele\n"
            "1\t100\tA\tG\t0.5\t1\t100\t\n"
            "2\t200\tC\t\t0.2\t2\t200\tT\n"          # other allele inferred
            "X\t300\tA\tG\t0.1\tX\t300\t\n"          # non-autosome dropped
            "3\t400\tAT\tG\t0.1\t3\t400\t\n"         # indel dropped
            "4\t500\tA\tG\t0.1\t\t\t\n")             # not harmonised: dropped
    with gzip.open(p, "wt") as f:
        f.write(text)
    return p


def test_pgs_read_score(tmp_path):
    df, hdr = pgs.read_score(_pgs_file(tmp_path))
    assert list(df["pos"]) == [100, 200] and list(df["other"]) == ["G", "T"]
    assert hdr["trait_reported"] == "Test" and hdr["n_rows_usable"] == "2"


def test_pgs_or_weights_become_log(tmp_path):
    df, _ = pgs.read_score(_pgs_file(tmp_path, "OR"))
    assert df["weight"].iloc[0] == pytest.approx(-0.6931, abs=1e-3)


def test_pgs_missing_columns_raises(tmp_path):
    p = tmp_path / "bad.txt.gz"
    with gzip.open(p, "wt") as f:
        f.write("#weight_type=beta\nchr_name\teffect_allele\teffect_weight\n1\tA\t0.1\n")
    with pytest.raises(ValueError, match="missing harmonized columns"):
        pgs.read_score(p)


@pytest.mark.skipif(shutil.which("bcftools") is None, reason="bcftools not installed")
def test_clinvar_parser(tmp_path):
    info = ("##fileformat=VCFv4.2\n"
            + "".join(f'##INFO=<ID={i},Number=.,Type=String,Description="x">\n'
                      for i in ("CLNSIG", "CLNREVSTAT", "CLNDN", "GENEINFO", "MC", "AF_TGP", "RS"))
            + "##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t1\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNDN=Test_disease;GENEINFO=BRCA1:672;MC=SO:0001587|nonsense;RS=1\n"
            "1\t200\t2\tA\tG\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;CLNDN=x;GENEINFO=G2:1;MC=SO:1|a;RS=2\n"
            "1\t300\t3\tA\tG\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=no_assertion_criteria_provided;CLNDN=y;GENEINFO=G3:1;MC=SO:1|b;RS=3\n")
    v = tmp_path / "cv.vcf"
    v.write_text(info)
    df = clinvar.load_plp(v)
    assert list(df["gene"]) == ["BRCA1"] and df["stars"].iloc[0] == 2


def test_copynumber_formula(monkeypatch):
    # single-copy baseline 40x, gene copies summing to 40x -> 2 copies; summing to 80x -> 4 copies
    monkeypatch.setattr(copynumber, "require_tools", lambda *a: None)     # CI has no samtools
    monkeypatch.setattr(copynumber, "baseline_depth", lambda cram, fai: 40.0)
    monkeypatch.setattr(copynumber, "bedcov", lambda cram, regions, min_mapq=0: [40.0 / len(regions)] * len(regions)
                        if regions != copynumber.LPA_UNIQUE else [40.0] * len(regions))
    out = copynumber.estimate("x.cram", "x.fai")
    assert out["loci"]["GSTM1"]["estimated_total_copies"] == pytest.approx(2.0)


def test_pgx_cyp2c19_phenotypes():
    def geno(n2, n3, n17):
        rows = [("rs4244285", n2), ("rs4986893", n3), ("rs12248560", n17)]
        return pd.DataFrame([{"rsid": r, "n_effect": n, "status": "called"} for r, n in rows])
    ph = lambda *a: pgx.summarize(geno(*a))[0]["phenotype"]
    assert ph(0, 0, 0).startswith("no tested") and ph(0, 0, 1) == "rapid metabolizer"
    assert ph(1, 0, 0).startswith("intermediate") and ph(2, 0, 0) == "poor metabolizer"


def test_abo_inference():
    from wgsreport import bloodgroup as bg

    def sites(intact, b1, b2):
        gt = lambda n: [1] * n + [0] * (2 - n)
        d = {}
        if intact: d[("9", 133257521)] = [("T", ["TC"], gt(intact), True)]
        if b1: d[("9", 133255935)] = [("G", ["T"], gt(b1), True)]
        if b2: d[("9", 133255928)] = [("C", ["G"], gt(b2), True)]
        return d
    assert bg.infer_abo(sites(0, 0, 0))["group"] == "O"           # absent = reference = O1/O1
    assert bg.infer_abo(sites(1, 0, 0))["genotype"] == "A/O"
    r = bg.infer_abo(sites(1, 1, 1)); assert r["group"] == "B" and r["genotype"] == "B/O"
    assert bg.infer_abo(sites(2, 1, 1))["group"] == "AB"
    assert bg.infer_abo(sites(2, 2, 2))["genotype"] == "B/B"
    assert bg.infer_abo(sites(0, 1, 1))["group"] == "undetermined"


def test_measured_compare_roundtrip(tmp_path):
    """Synthetic lab sheet (not anyone's data): parsing, treatment exclusion and notes must all work."""
    from wgsreport import measured
    csv = tmp_path / "labs.csv"
    csv.write_text("Group,,,,\nDate,LDL,HDL,Lpa,Hemoglobin\n2001-01-01,5.0,1.5,300,150\n2002-01-01,4.8,1.5,,150\n2003-01-01,2.0,1.5,,150\n")
    got = measured.load_wide(csv)
    assert got["ldl"] == [("2001-01-01", 5.0), ("2002-01-01", 4.8), ("2003-01-01", 2.0)]
    traits = [{"id": t, "kind": "quantitative", "label": t, "expected_sd": 0.5, "percentile": 70.0} for t in ("ldl", "hdl", "lpa")]
    rows = {r["id"]: r for r in measured.compare(traits, got, treated_from="2003-01-01")}
    assert rows["ldl"]["n_on_treatment"] == 1 and rows["ldl"]["mean"] == pytest.approx(4.9)
    assert "treatment" in rows["ldl"]["note"] and "KIV-2" in rows["lpa"]["note"]
    assert rows["lpa"]["verdict"] == "compare qualitatively"


def test_report_helpers():
    from wgsreport.report import _ordinal, _usable, _pct_range
    assert [_ordinal(n) for n in (1, 2, 3, 11, 12, 13, 21, 33, 42, 52, 62, 71, 92)] == \
        ["1st", "2nd", "3rd", "11th", "12th", "13th", "21st", "33rd", "42nd", "52nd", "62nd", "71st", "92nd"]
    assert _usable(100, 100) == "100%" and _usable(788056, 788116) == ">99.99%" and _usable(97, 100) == "97.0%"
    lo, hi = _pct_range(0.0, 503)
    assert lo < 50 < hi and hi - lo < 10


def test_measured_verdicts_and_flags():
    from wgsreport import measured
    t = lambda i, e: {"id": i, "kind": "quantitative", "label": i, "expected_sd": e, "percentile": 20.0}
    got = {"mcv": [("2001-01-01", 90.0), ("2002-01-01", 90.0), ("2003-01-01", 90.0)],
           "ldl": [("2001-01-01", 5.0), ("2002-01-01", 5.0)], "lpa": [("2001-01-01", 300.0)]}
    rows = {r["id"]: r for r in measured.compare([t("mcv", -0.5), t("ldl", 0.3), t("lpa", -0.5)], got)}
    assert "within the typical adult range" in rows["mcv"]["verdict"] and "microcyt" not in rows["mcv"]["verdict"]
    assert "3.5 mmol/L" in rows["ldl"]["clinical_flag"] and "without a confirmed unit" in rows["lpa"]["clinical_flag"]
    assert "50 mg/dL" not in measured.NOTES["lpa"].replace(">=50 mg/dL as separate", "")
