# Validation on a public genome (GIAB HG002)

HG002 (NA24385, an Ashkenazi Jewish male) is the Genome in a Bottle reference sample. It is public, so results can be published, unlike anyone's personal genome.
This run checks that the pipeline works end to end on an independently produced, high-confidence call set. **It tests robustness, not the accuracy of the scores or of the
optional external tools.**

Reproduce (about 160 MB download, then a few minutes):
```bash
B=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38
mkdir giab && cd giab && curl -O $B/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -O $B/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi && cd ..
wgsreport analyze --vcf giab/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz --work-dir work --out-dir giab/results
```

## What was observed (GIAB v4.2.1 GRCh38 benchmark, 4.05 million records)
| Check | Result | Reading |
|---|---|---|
| Build detection | accepted as GRCh38 | correct |
| Same file with a GRCh37 header (contig lengths swapped) | refused with a message | correct |
| Same file with positions shifted by 137 bp and no contig lengths | refused: 84% of 1,765 records disagree with the panel REF allele | correct; catches an unlabelled wrong build |
| REF-allele disagreements at scored positions in the real file | 0 of 650,000 records | the panel and the call set agree |
| Ancestry | EUR 100%, nearest 1000 Genomes populations TSI and IBS | consistent with Ashkenazi ancestry (Southern-European-like in 1000G PCA); the admixed-ancestry caveat applies |
| Genotype frequencies vs expectation | observed hom-ref 0.489 vs 0.483 expected | no sign of systematic genotype error |
| Score coverage | all scores >= 98% of variants usable | harmonisation works on an independent call set |
| Sex | "undetermined" | expected: the benchmark file has chromosomes 1-22 only; the report says so instead of guessing |
| Lp(a) explicit-variant set (PGS000667) | +19 nmol/L; 40% of the 503 EUR reference genomes have a total this high or higher | plausible for a non-carrier |

## More GIAB genomes (VCF only)
| Sample | Ancestry called | REF disagreements | Lp(a) score |
|---|---|---|---|
| HG001 (CEU) | EUR 100% | 0 of 644k | PGS002181 |
| HG003, HG004 (Ashkenazi) | EUR 100% | 0 of ~650k | PGS002181 |
| HG005-HG007 (Han Chinese trio) | EAS 100% | 0 of ~634k | PGS000689 (ancestry-aware switch; tier C) |

Open observation: the Ashkenazi trio all score high for height (z about +2 to +3) while HG001 does not. These genomes sit about 3.4 SD from the EUR centroid (the warning was at 4 SD and is now 3 SD), so ancestry-related bias in the EUR-standardised scores is plausible; measured heights were not checked.

## African-ancestry genomes (1000 Genomes NA18486 YRI, NA19625 ASW)
| Sample | Ancestry called | Warning | Reference group |
|---|---|---|---|
| NA18486 (YRI) | AFR 100%, 1.2 SD | none | AFR (661) |
| NA19625 (ASW, admixed) | AFR 72% / AMR 28%, 8.1 SD | ancestry warning fires | AFR (661) |

Both ran without errors. Per-trait portability notes appear (for example height, BMI and LDL have no published AFR accuracy; Lp(a) PGS has a published AFR R² of 0.00 against 0.34 in Europeans and the tier is recomputed). Duffy-null (rs2814778) was called with 2 copies in both.

Limits: both samples are members of the 1000G reference panel, so ancestry and percentiles include the sample itself; this tests the code paths, not accuracy. Only chr1-22 were extracted, so G6PD (chrX) was "assumed reference" (an artefact of the test file, the same as for any autosomes-only VCF). Neither sample carried HBB sickle or G6PD A-, so those paths have not fired on real African data. UGT1A1*6 and the locus-override check were not run in African ancestry.

## CRAM + VCF: NA12878 (HG001)
The 30x NYGC GRCh38 CRAM (ENA ERR3239334, 15.8 GB) was read remotely by `copy-number` (about 100 minutes; a full download was impractical). Baseline depth 32.8x.

| Locus | Estimate | Expected | Reading |
|---|---|---|---|
| SMN1+SMN2 | 4.3 | 4 | good |
| AMY1 | 6.4 | 2-20 | plausible |
| FCGR3B | 1.8 | 2 | good |
| C4A+C4B | 3.2 | 4 | about 0.8 low |
| HBA1+HBA2 | 5.0 | 4 | about 1 high; noisy locus (flanking depths 31x and 52x) |
| LPA KIV-2 block/unique depth ratio | 2.2 | n/a | no truth for repeat count |

Reading: the copy-number module is unbiased where loci are unique (SMN, FCGR3B) but the near-identical, subtelomeric HBA locus has an error of about +/-1 copy, so its calls are low-confidence (the report already says so). It is not systematically low, so a 2.9 estimate in another genome is not explained by a global low bias.

## Pharmacogenomics: PharmCAT against CDC GeT-RM consensus (60 public 1000 Genomes samples)
Public data only: the CDC consolidated GeT-RM table (363 samples) gives consensus diplotypes; genotypes come from the NYGC 30x GRCh38 1000 Genomes call set, cut to PharmCAT's regions, and were run through `pgx_pharmcat.run` with `--absent-to-ref` (the default). 60 samples with the most non-reference consensus calls were chosen (not a random sample). Script: `scripts/getrm_pgx_benchmark.py`.

| Gene | Exact diplotype agreement | Note |
|---|---|---|
| CYP2C9 | 59 / 60 | the one miss has consensus "*1/undefined" |
| CYP2C19 | 56 / 60 | all four misses are an extra `*38` called by PharmCAT (for example consensus *1/*1, PharmCAT *1/*38); the consensus predates that allele definition, so this is the metabolizer phenotype PharmCAT reports equals the phenotype of the consensus diplotype in all four (Normal, Rapid, Intermediate), so this is a nomenclature difference, not a wrong clinical call |
| TPMT | 60 / 60 | |
| CYP3A5 | 60 / 60 | |
| NUDT15 | 3 / 3 | only 3 samples with a consensus call |
| CYP2D6 | not called | PharmCAT needs an outside CYP2D6 call; the report says so rather than guessing |
| SLCO1B1, UGT1A1, CYP2B6, CYP4F2 | mostly disagree | not a fair comparison: the consensus uses older names (`*1A/*1B`, `*28+*60`, `*6`) that PharmCAT now expresses as `*37`, `*80+*28`, `*9` or `*1`. Needs a name-translation table before it can be scored |

Across CYP2C9, CYP2C19, TPMT, CYP3A5 and NUDT15, 238 of 243 comparisons agree exactly. This supports the wrapper and the `--absent-to-ref` assumption on a high-coverage call set; it does not test variant-only calls from other providers, CYP2D6, or copy-number-affected genes. Not tested: DPYD (the table's DPYD column was not compared).

## HLA typing (T1K) against GeT-RM: one sample only
The GeT-RM HLA reference set (108 DNA samples) cannot be benchmarked at scale from public sequence: only one of them, NA12273, is in the 1000 Genomes 30x set (the rest have no public WGS). T1K was run on that sample's MHC and unmapped reads (read remotely from the CRAM, about 2.5 minutes), two-field comparison with the consensus:

| Gene | Consensus | T1K | Reading |
|---|---|---|---|
| A | 02:01, 32:01 | 32:188 only (quality 2) | one allele wrong or missing; flagged low confidence |
| B | 35:01, 44:02 | 35:434 (tied), 44:02 (quality 2) | one allele replaced by a rare near-identical allele; flagged |
| C | 04:01, 07:04 | 07:04 | second allele not reported |
| DRB1 | 11:01, 13:05 | 11:01, 13:05 | correct (quality 8) |
| DRB3 / DQA1 | 02:02 / 05:05 | same | correct |
| DQB1 | 03:01 / 03:01 | 03:01 plus 03:276N | spurious rare null allele added |
| DPA1 / DPB1 | 01:03 / 02:01, 04:01 | same | correct |

Reading: roughly two thirds of alleles correct, the errors are rare near-identical alleles and dropped second alleles at low quality, and the report's low-confidence flag (quality below 10) marks most of them. One sample proves nothing about accuracy; this is a smoke test of the pipeline, and class I typing from 30x data should be treated as research-grade. A real benchmark needs public WGS for more GeT-RM samples or another truth set (for example 1000 Genomes HLA types from Gourraud 2014, not yet checked).

## Large-CNV depth method (`cnv_screen.py`)
The positive controls below test the depth method. Most of these duplications are no longer on the reported list (they are not actionable); of the reported changes, only the 17q12 region and the 22q11.2 region were exercised, both with duplications although the reported changes there are deletions.
**Positive controls (public 1000 Genomes, 32-36x).** Carriers were chosen from the 1000G 3,202-sample SV callset (`1KGP_3202.gatksv_svtools_novelins.freeze_V3`); only the carrier's region, two other regions, autosomal baseline windows (every second one) and chrX windows were read from the public CRAMs.

| Sample | Known duplication (callset) | Screen result | Median copies |
|---|---|---|---|
| NA20362 | 17q12 | duplication | 3.19 |
| HG03061 | 3q29 | duplication | 3.10 |
| HG00178 | 16p11.2 distal | duplication | 3.13 |
| NA20814 | 15q11.2 BP1-BP2 | duplication | 3.16 |
| HG00141 | 15q13.3 | duplication | 3.09 |
| HG00141 | 22q11.2 A-D | duplication | 3.11 |

6 of 6 called as whole-region duplications; in the same five genomes, 10 of 10 checks at two non-carrier regions (17q21.31, 7p22.1; medians 2.05-2.23) were not called. chrX read 0.97-0.99 copies in XY and 1.91-1.95 in XX samples. A separate run on NA12878 (XX) called none of the neuropsychiatric regions (medians 1.95-2.16; chrX 1.93).

**What this does not establish**
* **No deletion carrier**: the callset has no deletion of 200 kb or more in any screened region. One-copy depth is supported only by chrX in XY samples.
* Only whole-region events were tested; partial, atypical and NRXN1 exon-level events were not. NRXN1 exons vary (1.67-2.64 in a negative genome), so one low exon is inconclusive.
* The callset comes from the same short reads by other methods; it is not an independent wet-lab assay.
* The 17q11.2 (NF1) and 17p12 (PMP22) regions have no positive control.
* An early run on incompletely downloaded data reported "not detected" for regions with most bins missing; regions under 70% usable bins are now "not assessable".

## Alpha-globin depth estimate (`copynumber.py`): not usable for carrier status
On five public genomes with no HBA deletion in the 1000G SV callset, the HBA1+HBA2 estimate read 4.5-5.3 copies (expected 4). The report labels it "not reliable" and the clinical-test map lists alpha-thalassemia as not covered. Other paralog estimates in the same five also scattered (for example RHD+RHCE 3.3-4.3); treat them as research context.

## What this does not show
* No ground truth exists here for polygenic scores, so score values for HG002 are not validated (the benchmark VCF holds confident-region calls only; anything absent is assumed reference).
* PharmCAT, T1K (HLA), ExpansionHunter and Yleaf were not run on HG002, so their accuracy is untested here. Public truth sets for those (for example GeT-RM for pharmacogenes) would be the next step.
* One sample, one ancestry group, one caller. Other sequencing providers' VCFs are untested.

## Do large-effect loci override the genome-wide scores? (`scripts/locus_override_check.py`)
Lp(a) needed an explicit-variant cross-check because a HapMap3-style score barely weights its dominant variants. The same test was run for every tier A/B/C score on public
data only (the PGS weights, 1000 Genomes EUR genotypes, Ensembl coordinates). For each of 75 trait/variant pairs (candidates chosen from known biology, not from the scores) it asks: is the
variant in the score, and how well does the score computed from the *other* variants within 500 kb track it (r² of the candidate's dosage with that local score)? Low r² with no weight means the score
carries none of that variant's signal. The script does not estimate effect sizes: an override still needs a published effect and a unit-matched weight set.

| Verdict | Traits and loci (EUR, r² with the rest of the local score or weight in the score) |
|---|---|
| Captured, no override needed | bilirubin UGT1A1 (r² 0.96), urate ABCG2/SLC2A9 (in score; r² 0.7-0.8), triglycerides APOA5/GCKR/LPL, BMI FTO, CAD 9p21 and SORT1, T2D TCF7L2/SLC30A8/KCNJ11, height GDF5/HMGA2 (r² 0.4-0.7, in score), cystatin C CST3, HDL CETP/LIPC (in score) |
| **Missed: LDL and ApoB** | PCSK9 R46L (AF 2%, r² 0.19) and APOE e4 (AF 15%, r² 0.22) are in neither score. APOE is already computed; PCSK9 R46L is in the SNP panel. The LDL/ApoB results should say the score does not include them |
| **Missed: CAD** | LPA rs10455872 (AF 7%, not in score, r² 0.24) |
| **Partly missed** | AF 4q25 rs2200733 (r² 0.48, not in score); prostate 8q24 rs1447295 (r² 0.24, not in score) |
| **Weights look shrunken** | melanoma MC1R R151C/R160W are in the score at very small weights; D294H, D84E, R163Q are absent (AF about 1%). An explicit MC1R "R allele" count is the usual stratifier |
| Cannot be tested on EUR | HBB sickle (hba1c, hb, hct, rbc, mcv, mch), G6PD A- (hba1c, hb, hct), Duffy null (wbc, neut: r² 0.005, AF 0.6% in EUR), UGT1A1*6 (EAS). These matter in African or East Asian ancestry; rerun the script with the matching 1000 Genomes group |

Limits: LD is EUR only; candidates were chosen from known biology and the list is not exhaustive; r² is a redundancy measure, not a validation of the score; weights are in each score's own units.

**Acted on (report notes, not score replacements).** For the missed loci the report now shows the genome's genotype beside the score with the statement that the score does not include the variant: LDL and ApoB (APOE e4/e2, PCSK9 R46L, direction-aware), CAD (LPA rs10455872), AF (4q25 rs2200733), prostate (8q24 rs1447295) and melanoma (MC1R R alleles). The notes are driven by `related_rsids`, `related_effect`, `related_text` and `related_not_in_score` in `traits.yaml`. The ancestry-specific loci (HBB, G6PD, Duffy, UGT1A1*6) still need a check in the matching ancestry.
