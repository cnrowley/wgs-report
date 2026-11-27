import pytest

from wgsreport import query


def test_resolve_positions_and_bundled_panel_offline():
    r = query.resolve(["chr10:94842866", "rs12913832"], offline=True)
    assert r[0]["chrom"] == "10" and r[0]["pos"] == 94842866
    assert r[1]["chrom"] == "15" and r[1]["pos"] == 28120472       # bundled curated panel, no network


def test_unknown_rsid_offline_is_refused_not_guessed():
    with pytest.raises(SystemExit):
        query.resolve(["rs999999999"], offline=True)


def test_coords_file_overrides(tmp_path):
    f = tmp_path / "c.tsv"
    f.write_text("rsid\tchrom\tpos\nrs1\tchr2\t500\n")
    assert query.resolve(["rs1"], coords=f, offline=True)[0]["pos"] == 500


def test_plain_text_scoring_file(tmp_path):
    f = tmp_path / "s.txt"
    f.write_text("#comment\nchr_name\tchr_position\teffect_allele\teffect_weight\n6\t160589086\tg\t2.5\nX\t5\tA\t1\n")
    df = query.read_scoring_file(f)
    assert list(df["effect"]) == ["G", "A"] and df["pos"].iloc[0] == 160589086


def test_pgs_requires_exactly_one_source(tmp_path):
    with pytest.raises(SystemExit):
        query.query_pgs(tmp_path / "x.vcf.gz", None, None, tmp_path)
    with pytest.raises(SystemExit):
        query.query_pgs(tmp_path / "x.vcf.gz", "PGS000001", tmp_path / "f", tmp_path)


def test_vcf_records_parse(monkeypatch):
    class R:
        stdout = "100\tA\tG,T\tPASS\t1|2\n"
    monkeypatch.setattr(query, "run", lambda *a, **k: R())
    rec = query.vcf_records("x.vcf.gz", "1", 100)[0]
    assert rec["genotype"] == "G/T" and rec["filter"] == "PASS"
