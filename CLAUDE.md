# CLAUDE.md: guide for AI assistants (and humans) working on wgs-report

`wgs-report` turns one person's whole-genome VCF (and optionally CRAM) into a private HTML report: polygenic scores with
evidence tiers, ancestry, a curated SNP panel, ClinVar, pharmacogenomics, rare loss-of-function candidates, copy-number
estimates. Everything runs locally. Read `README.md` for usage and `docs/METHODS.md` for the maths.

## Hard rules (do not relax these without the owner's explicit say-so)
1. **No personal data in the repo, ever.** VCF/CRAM/pgen, work directories, `results.json`, `cram_checks.json`, lab
   spreadsheets, report HTML, sample IDs, names, emails, lab values, "my genotype is ..." examples. Run
   `bash scripts/check_no_personal_data.sh` before every commit; put the user's own identifiers in the untracked
   `.personal_patterns` file so the check can grep for them. Use synthetic data in tests. Never write a number or interpretation learned from a user's result into report text, notes or docs: report text must be computed from the input.
2. **Nothing leaves the machine.** Network is for public reference data only (see `sources.py`). Never put a user's
   genotype, sample ID, results or lab values in a URL, header, query or issue text.
3. **Do not bundle data you may not redistribute** (GWAS summary statistics, PGC3 weights built from them,
   third-party tools such as ADMIXTURE). Scores are downloaded or built locally.
4. **Report only what was tested.** No tables or paragraphs about genes, tests or loci that were not run or cannot be detected: leave the
   row or section out. Where "nothing found" could be read as an all-clear, one short sentence says what was checked. Longer limits belong in `docs/`.
5. **Honesty over completeness.** Every reported prediction carries a tier and its published accuracy. Never convert a
   percentile to absolute risk; never present a polygenic score as a diagnosis; state each result's limits next to it.
   The intended users are well versed in genetics, so do **not** hide results behind consent gates or paternalistic wording:
   APOE and the psychiatric score are computed by default (`--exclude apoe psychiatric` leaves them out). Intelligence/education
   scores were removed by the maintainer's decision and should not be re-added. Do not implement behavioural or personality
   predictions, attractiveness or similar scores. Do not add a score only because one exists: it must meet the inclusion rule at the top of
   `data/traits.yaml` or be explicitly tier C ("percentile only", on request).
6. **Downloaded code is untrusted.** Run third-party tools from their own directory with `python -I` or as binaries; do
   not widen permissions to make a blocked run work. Ask the user.

## Repo map
```
src/wgsreport/
  cli.py          sub-commands: prepare, build-score, analyze, doctor, fetch-reference, cram-check, copy-number, query, pharmcat, hla-type, genome-qc, repeat-expansions, y-haplogroup
  pipeline.py     `prepare`: PGS weights -> 1000G panel -> harmonise -> score panel -> PCA (no personal data)
  analyze.py      `analyze`: reads the VCF, scores, ancestry, per-trait results, writes report + results.json
  score.py        harmonisation (strand/palindromes), panel + user scoring, user dosages, callability check
  pgs.py          PGS Catalog download + parsing        sources.py   ALL external URLs, one place
  refpanel.py     1000G subset builder                   ancestry.py  PCA projection, kNN
  clinvar.py      ClinVar P/LP screen; ACMG gene set     rare.py      bcftools-csq loss-of-function scan
  snps.py         curated SNP panel + APOE               pgx.py       CPIC-style summary, CYP2D6 activity score
  genome_qc.py    sex-chromosome dosage, chromosome gains/losses, ROH   repeats.py  ExpansionHunter + STRchive thresholds   y_haplogroup.py  Yleaf
  pgx_pharmcat.py PharmCAT wrapper (jar + pipeline script, not bundled)   hla_typing.py  T1K wrapper + HLA association table
  variant_scores.py  short explicit-weight variant sets in trait units (Lp(a) PGS000667)
  query.py        `query snp|pgs|copy-number`: ad-hoc single-genome questions
  medication_report.py  printable drug-response summary (no identifiers; never framed as something to give a clinician)
  clinical_tests.py  map of commonly ordered clinical tests to what the report covers (data/clinical_tests.tsv)
  disclaimer.py   the research-use notice shown at the top of every output
  provenance.py   execution manifest (versions, PGS checksums) written to results.json
  bloodgroup.py   ABO/Rh inference from SNPs                              utils.py     run(), download(), log(), tool checks
  cram.py         read-level genotyping, mtDNA, FASTA from the EBI registry
  copynumber.py   read-depth copy number (LPA KIV-2, SMN, HBA, AMY1, C4, GSTM1 ...)
  cnv_screen.py   read-depth screen of a short list of actionable large CNVs; add one only if it meets the rules in data/cnv_actionable.tsv / cnv_disease.tsv
  carrier.py      carrier status for family planning (ClinGen AR/XL + ACMG tier 3, no frequency filter; results plus one limit sentence)
  custom_scores.py  P+T scores built from GWAS summary statistics   measured.py  compare with the user's labs
  doctor.py       health check for tools, bundled data, upstream formats      report.py   HTML
  data/           traits.yaml  snp_panel_src.yaml -> snp_panel.tsv  acmg_sf_v3.3.tsv  af_lof_genes.tsv  acmg_carrier_tier3.tsv
scripts/          make_screenshots.py (docs images, public genome only)  resolve_snp_panel.py (Ensembl REST)  admixture_run.py  check_no_personal_data.sh
tests/            offline tests; tests/test_network.py is marked `network`
```

## Tailoring recipes
**Add or swap a polygenic score.** (1) Find candidates with the PGS Catalog REST API (`/rest/score/search?trait_id=...`
uses MONDO/EFO ids; fetch `/rest/performance/search?pgs_id=...` for held-out metrics). Prefer scores whose evaluation
cohort is independent of the GWAS used to build them, and note the evaluation ancestry. (2) Add an entry to
`data/traits.yaml`: `r2` (quantitative) or `or_per_sd` (binary), `tier`, `metric_note` (say whether any figure is
derived rather than reported), `cite`; use `sex:` for sex-specific scores and `crosscheck:` for a second score. (3) Run
`wgsreport prepare ...` (rebuilds the panel when the variant set changes; takes ~10-30 min), then `pytest`.
Million-variant scores are fine; 6-7M-variant scores make the panel and `weights.pkl` very heavy (>8 GB RAM).

**Add a curated SNP.** Edit `data/snp_panel_src.yaml` (rsid, gene, category, effect allele, three notes for 0/1/2 copies),
then run `python scripts/resolve_snp_panel.py` to regenerate `snp_panel.tsv` (it checks the effect allele is one of the
site's alleles). `tests/test_data.py` enforces the contract. Add genotype-to-phenotype logic for pharmacogenes in `pgx.py`.

**Update the ACMG list.** Edit `data/acmg_sf_v3.3.tsv` (rename the file and `acmg_genes()` for a new version), update
the count in `tests/test_data.py`. Source of the current file: Lee et al. 2025 (doi:10.1016/j.gim.2025.101454).

**Change the ancestry reference.** `refpanel.py` and `ancestry.py` assume 1000G phase 3 GRCh38 and its `.panel` file; to
use another reference keep the same output files (panel pgen, `samples.tsv` with IID/pop/super_pop, `pca_files.json`).

**Add a module that some readers may not want** (e.g. a new disease-risk result). Compute it by default and add its name to the
`--exclude` choices in `cli.py` so a user can leave it out; do not hide it behind opt-in.

## When an upstream dependency breaks
Run `wgsreport doctor` (add `--offline` to skip the network). The weekly `upstream-contracts` CI job runs the same checks.

| Dependency | Where it is used | Typical break | Fix |
|---|---|---|---|
| PGS Catalog harmonised files | `pgs.py`, `sources.PGS_FTP` | URL layout or column names change | update `sources.py` / `pgs.read_score`; `tests/test_formats.py` shows the expected columns |
| PGS Catalog REST | recipes only | pagination or field renames | adapt the query, not the pipeline |
| ClinVar VCF | `clinvar.load_plp` | INFO field renamed (CLNSIG, CLNREVSTAT, ...) | update `sources.CLINVAR_INFO_FIELDS` and the parser together |
| Ensembl GFF3 | `rare.py` | pinned release removed from FTP | bump `sources.ENSEMBL_RELEASE`; rerun `--lof-scan` on a known variant |
| EBI CRAM registry | `cram.fetch_reference` | endpoint moves | any source of the same MD5-addressed sequence works; keep the checksum check |
| bcftools / samtools / plink2 | everywhere | flag changes (e.g. plink2 `--pmerge` was unfinished in a6) | `doctor` prints versions; pin in `environment.yml` |
| 1000 Genomes | `refpanel.py` | file naming | `pipeline.prepare(kg_pattern=...)` (not yet a CLI flag); the `.panel` file must sit beside the VCF directory |
| Python libraries | everywhere | pandas/numpy behaviour changes | CI runs on a current Python; pin only if forced |

Debug order: `doctor` -> the failing parser's test in `tests/test_formats.py` -> fetch one real file by hand and diff its
header against the synthetic fixture -> fix the parser and extend the fixture so the regression is caught offline.

## Things learned the hard way
* PGS downloads must be atomic per process (`.part` files include the PID) and are integrity-checked with gzip.
* 1000G has half-calls: pass `--vcf-half-call missing` to plink2.
* A variant-only VCF cannot prove an absent site is reference. `score.callability_check` is a population-level diagnostic
  only; `cram-check` confirms the curated panel from reads. Say "assumed" where it is assumed.
* Palindromic SNPs with MAF >= 0.42 are dropped, others use the reference allele frequency to resolve strand.
* Percentiles from a ~500-person reference move by several points when the panel is rebuilt. Do not present decimals.
* PCs: project with the allele loadings and calibrate the scale against plink2's own projection.
* Copy-number from depth counts all reads (MAPQ 0 included); report diploid totals, not allele-specific counts.
* Do not infer sex-specific, ancestry-specific or clinical claims beyond what a score's own evaluation supports.

## External tools (PharmCAT, T1K): lessons
* They are downloaded by the user into their own directory, never bundled or vendored. The PharmCAT pipeline script imports its sibling
  `pcat` package, so `python -I` cannot be used: run it from its own directory with its own virtualenv (`pgx_pharmcat.run` does).
* PharmCAT matcher output can contain several tied diplotypes when data is missing. `summarize` must report "ambiguous", never take
  the first one (when data is missing, several diplotypes tie and the first can carry a misleading phenotype).
* PharmCAT does not accept tandem notation with a CYP2D6 *68 hybrid; pass the accompanying allele (`normalize_outside`). Outside calls go
  through the jar's `-po`, not the pipeline script.
* T1K `t1k-build.pl -o ../hlaidx` can write `.._dna_seq.fa`; rename to `hlaidx_dna_seq.fa`. Alleles with quality <= 0 are ignored.
* `samtools idxstats` scans a whole CRAM (a CRAM index has no read counts): use sampled windows (`genome_qc.dosage`). Acrocentric chromosomes
  look low by total depth because of reference gaps; medians of windows fix that.
* ExpansionHunter counts pure repeats, so interruptions (ATXN1 CAT, FMR1 AGG) are invisible: an 'intermediate' call at an interrupted locus can be benign.
  Only interpret a locus when the counted unit matches the pathogenic motif (`repeats.canonical`). Thresholds are from STRchive, never from memory.
* Yleaf asks an interactive question if its output directory exists: pass `-force`; it never needs the network with `--no-update-check`.
* Fixtures in `tests/test_tools.py` are hand-built from real output shapes. Never commit real PharmCAT/T1K output.

## Verification before finishing a change
`pytest` (offline, ~1 s) -> `wgsreport doctor` -> `bash scripts/check_no_personal_data.sh` -> for pipeline changes, run
`prepare`+`analyze` on a public 1000G genome (e.g. HG00096) rather than anyone's personal data.
