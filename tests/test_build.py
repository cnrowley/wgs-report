"""Wrong-reference-build protection: a GRCh37 VCF must be refused, a GRCh38 one accepted."""
from wgsreport import build


def hdr(chr1, prefix="chr", extra=""):
    chr2 = {248956422: 242193529, 249250621: 243199373}.get(chr1, 0)
    return f"##fileformat=VCFv4.2\n{extra}##contig=<ID={prefix}1,length={chr1}>\n##contig=<ID={prefix}2,length={chr2}>\n#CHROM\tPOS\n"


def test_header_grch38_and_grch37():
    assert build.header_build(hdr(248956422))[0] == "GRCh38"
    assert build.header_build(hdr(249250621, prefix=""))[0] == "GRCh37"          # b37 naming, no chr prefix
    assert build.header_build(hdr(249250621))[0] == "GRCh37"                    # hg19 naming
    assert build.header_build(hdr(12345))[0] == "unknown"                       # custom contig lengths


def test_header_falls_back_to_reference_line():
    assert build.header_build("##reference=file:///refs/hg19.fa\n#CHROM\n")[0] == "GRCh37"
    assert build.header_build("##reference=GRCh38_no_alt.fa\n#CHROM\n")[0] == "GRCh38"
    assert build.header_build("#CHROM\n")[0] == "unknown"


def test_wrong_build_is_refused_but_right_build_passes():
    assert "GRCh37" in build.check_build(hdr(249250621))
    assert build.check_build(hdr(248956422)) is None
    assert build.check_build(hdr(248956422), {"from_vcf_het": 5000, "from_vcf_hom_alt": 3000, "diploid_mismatch_ref": 20}) is None


def test_ref_allele_mismatch_catches_unlabelled_wrong_build():
    # header gives no clue, but a quarter of records disagree with the GRCh38 REF allele (what a shifted build looks like)
    stats = {"from_vcf_het": 4000, "from_vcf_hom_alt": 2000, "diploid_mismatch_ref": 2000}
    err = build.check_build(hdr(12345), stats)
    assert err and "different reference build" in err
    # too few records to judge: no false alarm
    assert build.check_build(hdr(12345), {"from_vcf_het": 20, "from_vcf_hom_alt": 10, "diploid_mismatch_ref": 20}) is None


def test_inconsistent_header_is_not_trusted_either_way():
    mixed = "##contig=<ID=chr1,length=249250621>\n##contig=<ID=chr2,length=242193529>\n"
    assert build.header_build(mixed)[0] == "unknown"


def test_sex_from_x_heterozygosity_handles_haploid_calls(monkeypatch, tmp_path):
    from wgsreport import analyze

    class R:
        def __init__(self, stdout): self.stdout = stdout

    def fake(cmd, *a, **k):
        if cmd[:3] == ["bcftools", "view", "-h"]:
            return R("")
        gts = ["1:35", "0:20", "1:30", "0/1:10"]       # haploid male-style calls with one stray diploid het
        return R("".join(f"chrX\t{i}\t.\tA\tG\t.\tPASS\t.\tGT:DP\t{g}\n" for i, g in enumerate(gts)))
    monkeypatch.setattr(analyze, "run", fake)
    assert analyze.detect_sex(tmp_path / "x.vcf", "auto")[0] == "male"
