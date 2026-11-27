# Data sources and licences

This repository contains **code and small curated tables only**. Everything below is downloaded on your machine at
run time or supplied by you; none of it is bundled.

| Item | Source | Notes |
|---|---|---|
| Polygenic score weights | PGS Catalog (https://www.pgscatalog.org), GRCh38-harmonised files | Each score has its own licence and publication; cite them. Terms: https://www.pgscatalog.org/about/ |
| Reference panel | 1000 Genomes phase 3, GRCh38 (user downloads) | Public data; you supply `--kg-dir` |
| GWAS summary statistics (optional local scores) | PGC3 schizophrenia (figshare 19426775, CC BY 4.0) | Read each dataset's terms before use; weights built from them are not redistributed here |
| ClinVar | NCBI | Public domain |
| Gene annotation (rare loss-of-function scan) | Ensembl GFF3 | See Ensembl terms |
| Reference FASTA for a CRAM | EBI CRAM reference registry (by MD5) | `wgsreport fetch-reference` |
| ADMIXTURE, haplogrep, Cyrius | Third-party tools you install yourself | Not bundled; check their licences (ADMIXTURE is free for academic use) |
| SNP panel (`snp_panel_src.yaml`) | Hand-curated from public literature | Coordinates resolved through the Ensembl REST API |
| AF gene set (`af_lof_genes.tsv`) | Supplement (eTable 4) to "Contributions of Common, Rare, and Somatic Genetic Variants to Incidence of Atrial Fibrillation", JAMA Cardiol 2025, doi:10.1001/jamacardio.2025.3664 | Carrier odds ratios are from that paper |

`acmg_sf_v3.3.tsv` is the ACMG SF v3.3 table (Lee et al. 2025) transcribed from the PMC author manuscript, 84 genes with phenotype and inheritance. Update it when ACMG publishes a new version.

## ACMG tier-3 carrier-screening genes
`data/acmg_carrier_tier3.tsv` lists the 113 genes of Gregg AR et al., Genet Med 2021;23:1793-1806 (doi:10.1038/s41436-021-01203-z), Tables 1-6, transcribed as a plain fact list (symbol, first listed condition, inheritance). It is used only to flag ClinVar hits in carrier-screening genes; the report states this is not a carrier test.

## External tools (not bundled)
PharmCAT (Apache-2.0, PharmGKB/CPIC) and T1K (GPL-3.0) are run as separate programs that the user downloads; this repository contains only wrappers
and parsers of their output, so their licenses do not apply to it. CPIC recommendation text appears in the user's own local report (CPIC data is
CC0). IPD-IMGT/HLA is downloaded by T1K's build script under the database's own terms.

## Added data and tools (second batch)
* `data/gene_inheritance.tsv`: recessive and X-linked genes from ClinGen gene-disease validity curations (CC0), rebuilt by `scripts/build_gene_inheritance.py`.
* `data/cnv_actionable.tsv`, `data/cnv_disease.tsv`: short facts and effect estimates with citations, summarised from the cited publications and guidance.
* `data/cancer_predisposition_clingen.tsv`: gene list derived from ClinGen gene-disease validity curations (CC0), rebuilt by `scripts/build_cancer_genes.py`.
* `data/repeat_thresholds.tsv`: repeat thresholds from STRchive (MIT), rebuilt by `scripts/build_repeat_thresholds.py`.
* ExpansionHunter (Apache-2.0) and Yleaf (MIT) are run as separate programs the user downloads; nothing of theirs is bundled.
