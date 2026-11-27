# What you provide and how to run it

## 1. Your input
| Input | Required? | Format | Notes |
|---|---|---|---|
| Whole-genome variant calls | **Yes** | bgzipped VCF (`.vcf.gz`) **with a `.tbi`/`.csi` index**, reference build **GRCh38** | One sample. Standard `GT` field and a `FILTER` column. `chr1` and `1` naming both work. The tool refuses GRCh37/hg19 (see below). |
| Alignments | No (unlocks more modules) | CRAM or BAM + index, the matching reference FASTA | Needed only for: read-level confirmation of the SNP panel, HLA typing, repeat expansions, copy number (LPA KIV-2, SMN, HBA...), chromosome dosage QC, mtDNA, Y haplogroup. |
| Your lab values | No | spreadsheet, see `--measured` | Compares measured traits with the genetic prediction. |

What the VCF must be: calls for **your whole genome**, typically from a sequencing provider (DRAGEN, GATK, DeepVariant output). It need not be a gVCF. Exome, array
or 23andMe-style files are not supported (too few positions for the scores). **CRAM alone is not enough:** the tool does not call variants. Call a VCF first
(`bcftools mpileup -f ref.fa x.cram | bcftools call -mv -Oz -o x.vcf.gz`, or GATK/DeepVariant) and index it with `tabix -p vcf x.vcf.gz`.

A variant-only VCF cannot show that a site is homozygous reference; the report says "assumed" wherever it relies on that. A CRAM lets the key checks use reads.

**Build check.** Before scoring, the contig lengths in the header (chr1 = 248,956,422 for GRCh38) and, after reading genotypes, the agreement between the VCF's REF
alleles and the panel's decide whether the file is GRCh38. A GRCh37/hg19 file stops with a message. If your provider gave GRCh37, lift it over (CrossMap or
`bcftools +liftover`) or ask for GRCh38. `--skip-build-check` overrides it.

## 2. Install
```bash
conda env create -f environment.yml && conda activate wgsreport && pip install -e .
wgsreport doctor            # tools, bundled data, upstream formats (add --offline to skip the network)
```
Linux is the tested platform (macOS should work for the core; Windows needs WSL). Needs bcftools, plink2, samtools (CRAM modules), about 70 GB for the 1000 Genomes
reference and 8 GB RAM.

## 3. One-off: build the reference (no personal data)
```bash
wgsreport prepare --kg-dir /path/to/1000G/grch38 --work-dir work     # 15-30 min; re-run after changing traits.yaml
```
`--kg-dir` holds `ALL.chr{1..22}_GRCh38.genotypes.20170504.vcf.gz` (+ `.tbi`) with `integrated_call_samples_v3.20130502.ALL.panel` beside it.

## 4. Per person
```bash
wgsreport analyze --vcf my.vcf.gz --work-dir work --out-dir results
open results/report.html
```
Outputs: `report.html`, `results.json`, `snp_panel_genotypes.tsv`, `clinvar_pathogenic_hits.tsv`. Everything stays on your machine.

Useful flags: `--sex male|female`, `--exclude apoe psychiatric`, `--measured labs.xlsx`, `--skip-build-check`, `--clinvar-min-stars 1`.

## 5. Optional modules (each adds a section; run, then pass the JSON to `analyze`)
| Module | Needs | Command | `analyze` flag |
|---|---|---|---|
| Read-level SNP checks, mtDNA | CRAM + FASTA | `wgsreport cram-check --cram x.cram --fasta ref.fa ...` | `--cram-json` |
| Copy number (LPA, SMN, HBA ...; actionable large CNVs) | CRAM + FASTA | `wgsreport copy-number --cram x.cram --fasta ref.fa` | `--copynumber-json` |
| Rare loss-of-function scan | VCF + FASTA | (inside analyze) | `--lof-scan --fasta ref.fa --kg-dir ...` |
| Pharmacogenomics | VCF, PharmCAT (download) | `wgsreport pharmcat --vcf my.vcf.gz ...` | `--pharmcat-json` |
| HLA typing | CRAM, T1K (download) | `wgsreport hla-type --cram x.cram ...` | `--hla-json` |
| Chromosome dosage, ROH | VCF (+CRAM) | `wgsreport genome-qc --vcf ... --fasta ref.fa` | `--qc-json` |
| Repeat expansions | CRAM, ExpansionHunter | `wgsreport repeat-expansions --cram x.cram ...` | `--repeats-json` |
| Y haplogroup | CRAM/VCF, Yleaf | `wgsreport y-haplogroup ...` | `--y-haplogroup-json` |

Each command prints its own `--help`. Third-party tools are downloaded by you into their own directory and never bundled (`docs/DATA_AND_LICENSES.md`).

## 6. If something goes wrong
`wgsreport doctor`, then `docs/EXTENDING.md` ("When an upstream dependency breaks"). The built-in CPIC-style pharmacogenomics summary works without PharmCAT.


## Quick queries: `wgsreport query`
Ask one question of one genome without running the full report. All of it runs locally.

```bash
# SNPs: rsIDs or GRCh38 positions; the VCF must be bgzipped and indexed (tabix -p vcf file.vcf.gz)
wgsreport query snp rs10455872 rs35705950 chr10:94842866 --vcf my.vcf.gz
wgsreport query snp rs4988235 --vcf my.vcf.gz --cram my.cram --fasta hg38.fa      # also counts the reads
wgsreport query snp rs4988235 --vcf my.vcf.gz --offline --coords my_coords.tsv     # no network (rsid, chrom, pos columns)

# A polygenic score: downloaded from the PGS Catalog, or a file you provide (no network)
wgsreport query pgs --vcf my.vcf.gz --pgs-id PGS000667 --cache pgs_cache
wgsreport query pgs --vcf my.vcf.gz --score-file my_score.txt.gz

# Copy number from read depth: named loci (SMN, HBA, GSTM1, AMY1, C4, UGT2B17, RHD, FCGR3B) and/or any region
wgsreport query copy-number --cram my.cram --fasta hg38.fa --locus SMN --region chr6:160577000-160651000
```
* rsID coordinates come from the bundled panel, then `--coords`, then Ensembl (only the rsID is sent; `--offline` forbids it). A genotype is never sent anywhere.
* A site missing from a variant-only VCF is reported as "no record: homozygous reference assumed"; add `--cram --fasta` to check it from reads.
* `query pgs` prints the raw additive sum (effect-allele copies x weight; the score's own units) and the variants carried. It is not standardised: use `analyze` for percentiles against a reference panel. A non-harmonized `--score-file` is read as GRCh38 `chr_name`/`chr_position`, so check its build.
* `query copy-number` gives a total diploid estimate with roughly plus or minus one copy of error; custom regions are unvalidated.


## Reading the output
[REPORT_GUIDE.md](REPORT_GUIDE.md) explains each section with screenshots; [MEDICATION_REPORT.md](MEDICATION_REPORT.md) describes the printable drug response summary.
