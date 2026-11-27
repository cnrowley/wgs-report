# Roadmap

Modules that are implemented: polygenic scores (PGS Catalog weights, plus scores built locally from GWAS summary
statistics), curated SNP panel, ClinVar P/LP screen, read-level genotyping from a CRAM, mtDNA variant calling,
CYP2D6 import and a pharmacogenomic summary, a rare loss-of-function scan (ACMG genes plus an AF gene set),
read-depth copy-number estimates, PharmCAT star-allele calls (`pharmcat`), HLA typing with T1K (`hla-type`), repeat expansions (`repeat-expansions`), Y haplogroup (`y-haplogroup`), sex-chromosome/ROH checks (`genome-qc`), supervised ADMIXTURE (`scripts/admixture_run.py`), comparison with your own lab results.

Planned, in rough priority order (none of these are implemented). Priorities were revised after two independent reviews: HLA typing,
a PharmCAT-based pharmacogenomics module, repeat-expansion screening, carrier screening (with explicit false-negative caveats for
SMN1/CFTR-type loci), and GIAB-based regression tests come first; telomere length, longevity scores and Y haplogroups come last.

1. **Callability model.** Per-site depth/mapping quality from a gVCF or the CRAM, so that "absent from the VCF" can be
   labelled hom-ref vs not callable per locus (today there is a population-level check only).
2. **ACMG secondary findings, properly.** Keep `acmg_sf_v3.3.tsv` current with each ACMG release; add variant
   classification beyond ClinVar P/LP (ACMG/AMP criteria) and de novo-style filtering of loss-of-function calls.
3. **Structural variants and repeat expansions** (Manta/ExpansionHunter or the aligner's own outputs).
4. **HLA typing:** implemented with T1K; still open: HLA-LA cross-check, alt-contig reads for alt-aware references, validation against GIAB truth sets.
5. **Haplogroups:** mtDNA via haplogrep (documented below), Y via Yleaf or similar.
6. **Carrier screening** over a curated recessive-gene list with consequence and population-frequency filters.
7. **Calibrated individual intervals** (PredInterval/CalPred-style) once a calibration cohort is available.
8. **More scores** only where published held-out accuracy meets the inclusion rule in `traits.yaml`.

## Ideas from related projects
Observed in the public README of `allelix` (a variant-level annotator, AGPL-3.0: implement the ideas independently, do not copy code):

* **Genome-build detection from the VCF itself** (contig lengths) and refusing a mismatch, instead of trusting the user. We only guard
  at `prepare` time.
* **Input quality summary in the report:** number of rows removed by the caller's FILTER, number of no-calls, and a plain statement
  when the VCF has no reference blocks. We count these internally (`score.py`) but do not surface them.
* **Explicit "not assessed" counts** for sites that could not be strand-resolved or where one position is split over several VCF
  records, and a "sources not used" line when a database is missing.
* **ClinVar conflict handling made visible:** list pathogenic records that were scored down because ClinVar's overall classification is
  conflicting or benign.
* **Rare-variant context:** gnomAD v4 allele frequencies (we use 1000 Genomes only) and a missense predictor such as AlphaMissense
  (check the licence first) for the rare-variant table.
* **Multi-sample VCF** support (`--sample`).
* **Per-variant pharmacogenomic evidence levels** (CPIC / ClinPGx) for genes without star-allele calling.
* Not adopted: SNPedia-based bulk annotation (its licence restricts reuse; magnitude scores invite over-interpretation).

Design rules for any new module: compute by default with an `--exclude` switch (users are expected to know genetics), report evidence tier and limitations next to each
finding, never convert a percentile to absolute risk, keep personal data out of the repository, do not bundle data
whose licence does not allow redistribution.

- Allele-specific LPA KIV-2 copy number from short reads (DRAGEN-style callers; ~50% of genomes phase), to combine with the explicit-variant set for Lp(a).

* Inherited cancer susceptibility: see docs/CANCER.md (current support vs investigation, acceptance criteria).

* CNVs and carrier status: the actionable-CNV screen needs a deletion positive control and its follow-up text checked against full sources; exon-level deletions/duplications, SMN1, alpha-globin and F8 need validated dedicated callers. Intelligence scores were removed by the maintainer and are not planned.
