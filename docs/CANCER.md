# Inherited cancer susceptibility: assessment, current support, acceptance criteria

This page answers a proposal to add hereditary-cancer and polygenic cancer-risk analyses. It separates **what the tool does today** from
**what is under investigation**, and states the tests a new cancer feature must pass before it is shipped. Numbers marked *(unverified)* were
recalled, not checked against a source in this repository; do not copy them into the report without a citation.

## Verdict
* Monogenic screening **already exists** (ClinVar P/LP in ClinGen Definitive/Strong cancer genes, plus a rare loss-of-function scan). The proposal's
  Phase 1 is mostly an organisation and honesty problem, not a new capability. This release adds a dedicated section (below).
* Polygenic cancer scores must pass the project's inclusion rule (`data/traits.yaml`: twin heritability about 0.4 or more and held-out
  R2 >= 0.09 or OR/SD >= about 1.6). By that rule **prostate and melanoma are in; breast, colorectal, ovarian, pancreatic, lung, bladder and kidney are not**
  (published twin heritabilities (Mucci 2016, unverified here) are roughly 0.31, 0.15, 0.39, 0.36, 0.18, 0.31 and 0.38). Breast has an unusually well-validated score
  (Mavaddat 2018, PGS000004, 313 variants; OR/SD about 1.6, unverified), so it is a candidate for an explicit, owner-approved exception, not an automatic addition.
* Integrated absolute-risk models (CanRisk/BOADICEA) are **deferred**: they need family history and age that a VCF does not contain, and the hosted service
  would receive genotypes, which conflicts with the "nothing leaves the machine" rule. Local licensing is not established.

## Component by component
| Proposal component | Current support | Needs | Decision |
|---|---|---|---|
| P/LP variants in BRCA1/2, PALB2, CHEK2, ATM, RAD51C/D, BRIP1, Lynch genes, CDKN2A, STK11, BAP1, POT1, CDH1, CTNNA1, TP53, VHL, FH, FLCN, MET, SDHB/C/D | **Yes**: all are on the ClinGen Definitive/Strong list in `cancer_predisposition_clingen.tsv`; exact-match ClinVar screen on the VCF | none | shipped; now a dedicated report section |
| HOXB13 | **Added** by hand (not on the ClinGen list; G84E founder variant) | none | shipped, labelled as added from the literature |
| ClinVar review status | >= 2 stars in the main table; **1-star P/LP in cancer genes now shown in a separate "lower review status" table** | none | shipped |
| Recessive-only genes (BLM, MUTYH ...) | one copy now shown as carrier status, not as "predisposition" | none | shipped |
| VUS reporting | not shown | ClinVar VUS parsing, consequence/frequency filters, a clear "not evidence of risk" table | investigate |
| Penetrance numbers | not shown | a verified, cited number per gene, variant class and sex; ranges vary by study | investigate (no number without a citation) |
| Exon-level CNV / structural variants | not assessed (VCF); read-depth is crude | a validated per-exon caller; GeT-RM / GIAB truth; MLPA-confirmed samples | investigate; BRCA1, PMS2 and MSH2/EPCAM are the main loss |
| Per-gene coverage / callability | not assessed | CRAM depth over gene exons, or gVCF reference blocks | investigate; this is what separates "no variant" from "not assessed" |
| PMS2, CHEK2 1100delC region, pseudogene genes | flagged as poorly resolved | read-level confirmation | investigate |
| Prostate PGS | **Yes**, tier A (PGS000662) | none | shipped |
| Melanoma PGS | **Yes**, tier C (PGS019972) | none | shipped |
| Breast / colorectal / ovarian / pancreatic PGS | no | owner decision on an exception; independent validation, ancestry data, deposited weights | investigate; not added |
| Lung / bladder / kidney PGS | no | the cited gains come from clinical covariates (smoking, age), not the score | defer |
| CanRisk / BOADICEA | no | family history, age, licence, local execution | defer |

Notes on the evidence offered with the proposal. The 8-cancer UK Biobank study (JNCI Cancer Spectrum 2024, doi:10.1093/jncics/pkae008; the PDF is in `pdf/`)
trained its own scores from GWAS summary statistics with methods such as LDpred2 and PRS-CS inside UK Biobank, so those scores are **not** PGS Catalog
entries and have no deposited weights to use. Its text confirms AUCs of 0.681 (prostate) and 0.533 (ovarian); the other rows were not individually checked.
The study is also evidence that PGS adds little beyond clinical factors for several cancers, which argues for restraint.

## What the report now shows (section "Inherited cancer susceptibility")
A. P/LP variants in cancer-predisposition genes, split into dominant/X-linked genes, single copies in recessive-only genes (carrier status), and a
lower-review-status table. B. The cancer PGS that exist (prostate, melanoma) with tier and percentile, never as a probability. C. What was and was not assessed.
Monogenic and polygenic results are never merged into one number.

## Acceptance criteria for any new cancer feature
1. **Citation**: every penetrance figure, gene-cancer association and PGS metric has a DOI or source in the code or data file; none from memory.
2. **Inclusion rule**: a PGS meets the `traits.yaml` rule, or is an explicit tier C / owner-approved exception that says why.
3. **Independence**: the evaluation cohort is independent of the discovery GWAS, or the overlap is stated; the evaluated ancestry is stated.
4. **No merging**: monogenic findings and polygenic percentiles are not combined into one risk without a published, validated integrated model whose inputs are all present.
5. **No absolute risk** without the inputs the model needs (age, sex, family history); otherwise the genetic component only.
6. **Negative results**: a gene with no finding is never reported as "negative" unless its coverage was assessed; otherwise the report says "not assessed".
7. **Variant classes**: the report states which classes were assessed (SNV/indel, CNV, SV, splice, deep intronic) per gene or globally.
8. **VUS**: shown separately, labelled as not evidence of risk, never in a risk sentence.
9. **Truth data**: a benchmark on public data (GeT-RM, GIAB, 1000G, ClinVar-curated controls) is run and recorded in `docs/VALIDATION.md` before the feature is called supported.
10. **Offline tests** cover the logic with synthetic data; no personal data enters the repository; `check_no_personal_data.sh` passes.
11. **Privacy**: no genotype, family history or result leaves the machine; any service that would receive them is rejected or run locally.
12. **Maintainability**: the gene list is a data file with a regeneration script and a source date, and a test enforces its columns.

## Under investigation (not supported today)
VUS table; gene-level coverage and callability; exon-level CNV/SV detection and its truth set; penetrance data with verified citations; owner decision on
a breast (and possibly colorectal) PGS exception; local use of CanRisk/BOADICEA and its licence; ancestry portability of each cancer PGS.
