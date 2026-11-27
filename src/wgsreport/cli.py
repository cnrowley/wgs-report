"""Command-line interface: `wgsreport prepare` and `wgsreport analyze`."""
from __future__ import annotations

import argparse
from pathlib import Path


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(prog="wgsreport", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("prepare", help="one-off: download scores, build reference panel (no personal data)")
    p.add_argument("--kg-dir", required=True, type=Path,
                   help="directory with 1000 Genomes GRCh38 per-chromosome VCFs and the sample .panel file")
    p.add_argument("--work-dir", required=True, type=Path)
    p.add_argument("--jobs", type=int, default=6, help="chromosomes processed in parallel")
    p.add_argument("--threads", type=int, default=8)

    b = sub.add_parser("build-score", help="build a score locally from GWAS summary statistics (not redistributable as weights)")
    b.add_argument("name", choices=["pgc3-scz"])
    b.add_argument("--sumstats", required=True, type=Path, help="summary-statistics file for the chosen score")
    b.add_argument("--kg-dir", required=True, type=Path)
    b.add_argument("--work-dir", required=True, type=Path)
    b.add_argument("--jobs", type=int, default=8)
    b.add_argument("--threads", type=int, default=4)

    d = sub.add_parser("doctor", help="check tools, bundled data and upstream endpoints/formats (use after upgrades or failures)")
    d.add_argument("--offline", action="store_true", help="skip network checks")

    f = sub.add_parser("fetch-reference", help="build a primary-contig GRCh38 FASTA matching your CRAM from the EBI reference registry")
    f.add_argument("--cram", required=True, type=Path)
    f.add_argument("--out", required=True, type=Path)

    c = sub.add_parser("cram-check", help="read-level checks: genotype the SNP panel from reads, call mitochondrial variants")
    c.add_argument("--cram", required=True, type=Path)
    c.add_argument("--fasta", required=True, type=Path)
    c.add_argument("--out", required=True, type=Path, help="JSON file to write (contains personal data; keep private)")
    c.add_argument("--sex", choices=["male", "female"], default="male")

    q = sub.add_parser("query", help="ad-hoc questions about one genome: SNP lookup, a PGS (download or file), copy number of a locus or region")
    qs = q.add_subparsers(dest="qcmd", required=True)
    qa = qs.add_parser("snp", help="genotype rsIDs or GRCh38 positions (chr10:94842866) from the VCF, optionally confirmed from reads")
    qa.add_argument("ids", nargs="+", help="rsIDs and/or chrom:pos")
    qa.add_argument("--vcf", required=True, type=Path, help="bgzipped, indexed VCF")
    qa.add_argument("--coords", type=Path, help="offline lookup table (columns rsid, chrom, pos[, ref, alt]); otherwise the bundled panel, then Ensembl")
    qa.add_argument("--offline", action="store_true", help="never contact Ensembl (rsIDs must be in --coords or the bundled panel)")
    qa.add_argument("--cram", type=Path); qa.add_argument("--fasta", type=Path, help="with --cram: also count the reads")
    qa.add_argument("--sex", choices=["male", "female"], default="male"); qa.add_argument("--json", action="store_true")
    qp = qs.add_parser("pgs", help="raw score of one PGS Catalog score (downloaded) or a scoring file you provide (offline)")
    qp.add_argument("--vcf", required=True, type=Path)
    qp.add_argument("--pgs-id", help="PGS Catalog id, e.g. PGS000667 (public file downloaded to --cache)")
    qp.add_argument("--score-file", type=Path, help="local scoring file (.txt/.gz); no network")
    qp.add_argument("--cache", type=Path, default=Path("pgs_cache")); qp.add_argument("--top", type=int, default=10)
    qc = qs.add_parser("copy-number", help="depth-based copy number of named loci and/or custom regions")
    qc.add_argument("--cram", required=True, type=Path); qc.add_argument("--fasta", required=True, type=Path, help="reference FASTA with .fai")
    qc.add_argument("--locus", action="append", help="name (substring) from the built-in list, e.g. SMN, HBA, GSTM1; repeatable")
    qc.add_argument("--region", action="append", help="custom GRCh38 region chrom:start-end; repeatable")

    n = sub.add_parser("copy-number", help="read-depth copy number (LPA KIV-2, SMN, HBA, AMY1, C4, GSTM1 ...) and a short list of actionable large CNVs")
    n.add_argument("--cram", required=True, type=Path)
    n.add_argument("--fasta", required=True, type=Path, help="FASTA (its .fai is used for baseline windows)")
    n.add_argument("--out", required=True, type=Path)

    ph = sub.add_parser("pharmcat", help="run PharmCAT (star alleles + CPIC recommendations) on your VCF; tool not bundled")
    ph.add_argument("--vcf", required=True, type=Path)
    ph.add_argument("--pharmcat-dir", required=True, type=Path, help="directory with pharmcat.jar and pharmcat_pipeline")
    ph.add_argument("--out", required=True, type=Path, help="output directory; the summary is written to <out>/pharmcat_summary.json")
    ph.add_argument("--no-absent-to-ref", action="store_true",
                    help="do not treat PGx positions missing from the VCF as reference (most genes then return 'No Result')")
    ph.add_argument("--cyp2d6", help="outside CYP2D6 diplotype from Cyrius/Aldy, e.g. '*1/*4' (PharmCAT cannot call CYP2D6 from a VCF)")

    q = sub.add_parser("genome-qc", help="sex-chromosome dosage, whole-chromosome gains/losses and runs of homozygosity")
    q.add_argument("--vcf", required=True, type=Path)
    q.add_argument("--cram", type=Path, help="alignment file (for chromosome dosage); optional")
    q.add_argument("--fasta", required=True, type=Path, help="FASTA whose .fai gives the contig lengths")
    q.add_argument("--out", required=True, type=Path)

    rp = sub.add_parser("repeat-expansions", help="repeat-expansion screen with ExpansionHunter (about 30 loci); tool not bundled")
    rp.add_argument("--cram", required=True, type=Path)
    rp.add_argument("--fasta", required=True, type=Path)
    rp.add_argument("--eh-dir", required=True, type=Path, help="unpacked ExpansionHunter release (contains bin/ExpansionHunter and variant_catalog/)")
    rp.add_argument("--out", required=True, type=Path, help="output directory; the summary is written to <out>/repeats.json")
    rp.add_argument("--sex", choices=["male", "female"], default="male")
    rp.add_argument("--threads", type=int, default=4)

    yh = sub.add_parser("y-haplogroup", help="Y-chromosome haplogroup with Yleaf (CRAM); tool not bundled")
    yh.add_argument("--cram", required=True, type=Path)
    yh.add_argument("--fasta", required=True, type=Path)
    yh.add_argument("--yleaf", required=True, type=Path, help="path to the yleaf executable")
    yh.add_argument("--out", required=True, type=Path, help="output directory (must not exist); the summary is written to <out>/y_haplogroup.json")
    yh.add_argument("--threads", type=int, default=4)

    hl = sub.add_parser("hla-type", help="classical HLA typing from a CRAM with T1K; tool not bundled")
    hl.add_argument("--cram", required=True, type=Path)
    hl.add_argument("--fasta", required=True, type=Path, help="GRCh38 FASTA matching the CRAM")
    hl.add_argument("--t1k-dir", required=True, type=Path, help="T1K directory containing run-t1k")
    hl.add_argument("--index", required=True, type=Path, help="T1K HLA index, e.g. hlaidx/hlaidx_dna_seq.fa")
    hl.add_argument("--out", required=True, type=Path, help="output directory; the summary is written to <out>/hla.json")
    hl.add_argument("--threads", type=int, default=4)

    a = sub.add_parser("analyze", help="score one personal VCF and write the report")
    a.add_argument("--vcf", required=True, type=Path, help="bgzipped GRCh38 VCF of your genome")
    a.add_argument("--work-dir", required=True, type=Path, help="output of `prepare`")
    a.add_argument("--out-dir", required=True, type=Path)
    a.add_argument("--skip-build-check", action="store_true", help="do not refuse a VCF that looks like GRCh37 or another non-GRCh38 build")
    a.add_argument("--sex", choices=["auto", "male", "female"], default="auto")
    a.add_argument("--sensitive", nargs="*", default=[], choices=["apoe", "psychiatric", "all"],
                   help="deprecated and ignored: APOE and the psychiatric score are now computed by default")
    a.add_argument("--exclude", nargs="*", default=[], choices=["apoe", "psychiatric"],
                   help="leave these modules out of the report (apoe: APOE genotype; psychiatric: schizophrenia score)")
    a.add_argument("--clinvar", type=Path, default=None, help="ClinVar GRCh38 VCF (downloaded if absent)")
    a.add_argument("--no-clinvar", action="store_true")
    a.add_argument("--clinvar-min-stars", type=int, default=2, choices=[1, 2, 3, 4],
                   help="minimum ClinVar review stars for the rare-disease screen (1 = also single-submitter records; more findings, more false leads)")
    a.add_argument("--treated-from", help="YYYY-MM-DD: lipid results on/after this date (e.g. statin start) are excluded from the comparison mean")
    a.add_argument("--measured", type=Path, help="spreadsheet/CSV of your own lab results (wide format) to compare with the genetic tendencies")
    a.add_argument("--height-mean", type=float, help="override reference mean height (cm) for your sex")
    a.add_argument("--height-sd", type=float, help="override reference height SD (cm) for your sex")
    a.add_argument("--cram-json", type=Path, help="output of `cram-check`: confirms/replaces VCF assumptions, adds mtDNA variants")
    a.add_argument("--cyp2d6", help="CYP2D6 diplotype from a star-allele caller, e.g. '*1/*4' (Cyrius, Aldy, DRAGEN)")
    a.add_argument("--lof-scan", action="store_true", help="scan ACMG genes for rare predicted loss-of-function variants (needs --fasta)")
    a.add_argument("--fasta", type=Path, help="GRCh38 FASTA matching the VCF (for --lof-scan)")
    a.add_argument("--kg-dir", type=Path, help="1000 Genomes GRCh38 VCF directory (allele-frequency filter for --lof-scan)")
    a.add_argument("--admixture-json", type=Path, help="output of scripts/admixture_run.py")
    a.add_argument("--copynumber-json", type=Path, help="output of `copy-number`")
    a.add_argument("--mt-haplogroup", help="mtDNA haplogroup label as printed by haplogrep")
    a.add_argument("--repeats-json", type=Path, help="repeats.json from `wgsreport repeat-expansions`")
    a.add_argument("--y-haplogroup-json", type=Path, help="y_haplogroup.json from `wgsreport y-haplogroup`")
    a.add_argument("--qc-json", type=Path, help="output of `genome-qc`")
    a.add_argument("--pharmcat-json", type=Path, help="pharmcat_summary.json from `wgsreport pharmcat`")
    a.add_argument("--hla-json", type=Path, help="hla.json from `wgsreport hla-type`")
    a.add_argument("--threads", type=int, default=8)

    args = ap.parse_args(argv)
    if args.cmd == "analyze":
        args.sensitive = {"apoe", "psychiatric"} - set(args.exclude)
    if args.cmd == "prepare":
        from .pipeline import prepare
        prepare(args.kg_dir, args.work_dir, jobs=args.jobs, threads=args.threads)
    elif args.cmd == "build-score":
        from .custom_scores import build_pt_score
        key = args.name.replace("-", "_")
        # PGC3 uses P<0.05 as in the paper
        build_pt_score(key, args.sumstats, args.kg_dir, args.work_dir, p_max=0.05,
                       jobs=args.jobs, threads=args.threads)
    elif args.cmd == "doctor":
        from .doctor import run_doctor
        raise SystemExit(run_doctor(args.offline))
    elif args.cmd == "fetch-reference":
        from .cram import fetch_reference
        fetch_reference(args.cram, args.out)
    elif args.cmd == "query":
        import json
        from . import query
        if args.qcmd == "snp":
            rows = query.query_snps(args.vcf, args.ids, args.coords, args.offline, args.cram, args.fasta, args.sex)
            print(json.dumps(rows, indent=1) if args.json else query.format_snps(rows))
        elif args.qcmd == "pgs":
            print(json.dumps(query.query_pgs(args.vcf, args.pgs_id, args.score_file, args.cache, args.top), indent=1))
        else:
            print(json.dumps(query.query_copy_number(args.cram, args.fasta, args.locus, args.region), indent=1))
    elif args.cmd == "copy-number":
        import json
        from . import cnv_screen
        from .copynumber import estimate
        fai = Path(str(args.fasta) + ".fai")
        cn = estimate(args.cram, fai)
        cn["cnv_screen"] = cnv_screen.estimate(args.cram, fai)
        Path(args.out).write_text(json.dumps(cn, indent=1))
    elif args.cmd == "repeat-expansions":
        import json
        from . import repeats
        eh = repeats.run_expansionhunter(args.eh_dir, args.cram, args.fasta, args.out, args.sex, args.threads)
        (args.out / "repeats.json").write_text(json.dumps(repeats.parse(eh), indent=1))
        for f in args.out.glob("eh_realigned.bam*"):
            f.unlink(missing_ok=True)      # realigned reads are personal data
    elif args.cmd == "y-haplogroup":
        import json
        from . import y_haplogroup
        hg = y_haplogroup.run_yleaf(args.yleaf, args.cram, args.fasta, args.out, args.threads)
        (args.out / "y_haplogroup.json").write_text(json.dumps(y_haplogroup.parse(hg), indent=1))
    elif args.cmd == "genome-qc":
        import json
        from .genome_qc import run_qc
        Path(args.out).write_text(json.dumps(run_qc(args.cram, args.vcf, Path(str(args.fasta) + ".fai")), indent=1))
    elif args.cmd == "pharmcat":
        import json
        from . import pgx_pharmcat
        calls = {"CYP2D6": args.cyp2d6} if args.cyp2d6 else None
        rep = pgx_pharmcat.run(args.vcf, args.out, args.pharmcat_dir, calls, absent_to_ref=not args.no_absent_to_ref)
        (args.out / "pharmcat_summary.json").write_text(json.dumps(pgx_pharmcat.load_summary(rep), indent=1))
    elif args.cmd == "hla-type":
        import json
        from . import hla_typing
        args.out.mkdir(parents=True, exist_ok=True)
        r1, r2 = hla_typing.extract_reads(args.cram, args.fasta, args.out / "mhc", args.threads)
        try:
            tsv = hla_typing.run_t1k(args.t1k_dir, args.index, r1, r2, args.out, args.threads)
            calls = hla_typing.parse_genotype(tsv)
            (args.out / "hla.json").write_text(json.dumps({"calls": calls, **hla_typing.interpret(calls)}, indent=1))
        finally:
            # reads are personal data and not needed afterwards, even after a failure (T1K keeps its own copies too)
            for f in [r1, r2, *args.out.glob("t1k_candidate*.f*"), *args.out.glob("t1k_aligned*.f*")]:
                f.unlink(missing_ok=True)
    elif args.cmd == "cram-check":
        from .cram import run_checks
        from .snps import load_panel
        run_checks(args.cram, args.fasta, load_panel(True), args.out, args.sex)
    else:
        from .analyze import analyze
        analyze(args)


if __name__ == "__main__":
    main()
