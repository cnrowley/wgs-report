# Methods

1. **Scores.** GRCh38-harmonized PGS Catalog weights (`hmPOS_GRCh38`). SNVs only, autosomes only.
2. **Reference panel.** 1000 Genomes phase 3 (2,504 people) restricted to biallelic SNVs at score positions.
   The panel fixes REF/ALT coding. Weights are converted to ALT-allele weights; strand-flipped alleles are
   resolved; palindromic (A/T, C/G) SNVs with MAF >= 0.42 are dropped.
3. **Your genotypes.** Read from the VCF; sites absent are homozygous reference; filtered or no-call sites
   are mean-imputed. Deletions overlapping a site are ignored.
4. **Standardization.** The panel is scored with plink2. For the superpopulation your genome matches, the score
   is regressed on the first 4 PCs and your residual is expressed in SD. PCs are from LD-pruned common variants;
   your projection is calibrated against plink2's own projection (R^2 > 0.999 required).
5. **Prediction.** E[trait SD | score z] = sqrt(R^2) z with an approximate residual population interval +-1.96 sqrt(1-R^2)
   under a simplified bivariate-normal model. R^2 is the published held-out European value; uncertainty in R^2 and in z,
   transferability, age/sex effects, phenotype definitions, environment, within-family accuracy and calibration in an individual are
   ignored, so this is not a calibrated personal prediction interval. A plausible range for each percentile is shown from the
   sampling error of a ~500-person reference (se of z = sqrt(1/n + z^2/2n)). Binary traits: relative odds =
   OR_per_SD^z, not absolute risk.
6. **Selection of traits.** Heritability >= ~0.4 and held-out R^2 >= 0.09 (or OR/SD >= ~1.6, with AUROC and liability-scale R^2 reported as context because they depend on prevalence and cohort).
   Tier C traits are exceptions made on request. Tiers A/B/C rank published *predictive performance* (held-out R^2 or OR/SD) only; they
   say nothing about replication, confidence intervals or evidence quality, which differ between scores.
   `results.json` carries a `provenance` block (tool and package versions, PGS file checksums) and the report ends with a one-line summary of it.
7. **SNP panel.** Coordinates resolved with `scripts/resolve_snp_panel.py` (Ensembl REST); the script verifies that each
   stated effect allele is one of the site's alleles, catching strand mistakes.
8. **ClinVar.** P/LP, no conflicting interpretations, >= 2 review stars, exact chr/pos/ref/alt match.

9. **Read-level checks (`cram-check`).** The panel SNPs are counted from reads (samtools mpileup, base and mapping
   quality >= 20). A site is called homozygous when the second allele is below 15% of reads, heterozygous when both
   alleles are at least 25%, otherwise ambiguous; depth below 8 gives no call. Reads confirm, or contradict, "absent
   from the VCF means hom-ref". Mitochondrial variants are called haploid with bcftools.
10. **Pharmacogenomic summary.** Simplified CPIC-style phenotypes from the panel genotypes and an externally supplied
    CYP2D6 diplotype (activity-score method). It assumes no other star alleles and is not a clinical report.
11. **Rare loss-of-function scan.** `bcftools csq` on PASS variants with the Ensembl GFF3; stop-gain, frameshift, canonical
    splice and start-loss in the ACMG gene list and an AF gene set; 1000 Genomes AF < 1%. Transcripts flagged NMD-only
    are ignored; multi-allelic indels with 1/2 genotypes are flagged as probable artefacts. Output is candidates for
    review, not diagnoses.
12. **Copy number from depth.** Total copies = 2 x (sum of mean depth over all reference copies) / (single-copy baseline
    depth), counting reads of every mapping quality because near-identical copies cannot be placed uniquely. The
    baseline is the median depth of 80 windows spread over the autosomes. Estimates are diploid totals, +-1 copy.
13. **ADMIXTURE.** `scripts/admixture_run.py` runs supervised ADMIXTURE against 1000 Genomes labels (4 continental groups;
    5 European populations). Proportions are model-based and not statements of origin.
14. **Percentile stability.** With about 500 reference individuals, a rebuilt reference panel can move a percentile by
    several points. Read percentiles loosely, especially in the tails.

Known limitations: indel/CNV/structural variants and HLA/CYP2D6 are not scored; weights not in 1000G are lost
(coverage is reported per score); relatives in 1000G; 1000G EUR is only 503 people.

## Actionable large copy-number variants
`wgsreport copy-number` also screens a short list (`cnv_screen.py`, `data/cnv_actionable.tsv`). Inclusion needs all three: strong evidence that the change causes disease (ClinGen dosage score 3), a substantial effect on a carrier, and a follow-up documented in a published source. The source's strength differs (clinical recommendations, a surveillance guideline, a "limited" actionability rating, case reports, or only a patient-organisation guidebook) and is printed with each result. The follow-up text was summarised from search excerpts and must be checked against the sources. Low-penetrance and non-actionable CNVs, including most neuropsychiatric ones, are deliberately not reported.

Depth is measured in 10 kb bins at MAPQ >= 30 relative to the autosomal median; only bins keeping at least 90% of total depth after that filter are used. A whole-region call needs a median of <= 1.35 or >= 2.65 copies with 80% of bins agreeing; a partial call needs a contiguous run of abnormal bins (up to 200 kb). Under 70% usable bins is "not assessable". This is not a genome-wide CNV search.

Other disease links (`data/cnv_disease.tsv`) are shown only for a flagged change, and only if replicated in an independent cohort (Estonian Biobank, Auwerx et al. 2024) or Bonferroni-significant in UK Biobank and already known (Crawford et al. 2019). `analyze --exclude psychiatric` keeps the screen and drops psychiatric estimates.

## Carrier status
Section 4b (`carrier.py`) reports ClinVar P/LP SNVs and short indels in genes that ClinGen curates (Definitive, Strong or Moderate) as autosomal recessive or X-linked (`data/gene_inheritance.tsv`, rebuilt by `scripts/build_gene_inheritance.py`) or that are on the ACMG tier-3 carrier list (Gregg 2021). No population-frequency filter is applied. Lower-review records are listed separately. Genes where the common cause is invisible to this method (SMN1, HBA1/2, FMR1, CYP21A2, GBA1, DMD, F8 and others) are printed in every report with the reason. It is not a laboratory carrier screen, and a negative result does not exclude carrier status.

## Traits decided by a few large-effect variants (Lp(a))
A genome-wide score is the wrong tool when one locus plus a structural repeat decide most of the trait. For Lp(a), `traits.yaml` gives the trait
a `variant_set` (PGS000667, 43 conditionally independent LPA-region variants with weights in nmol/L, from independent Danish cohorts, R2 0.60 in
UK Biobank Europeans). With `variant_set_primary` the report headline comes from that set (additive per copy of the effect allele, so a heterozygote
counts once) or from a carried `related_rsids` variant (rs10455872, rs3798220), compared with the guideline threshold of 125 nmol/L. The genome-wide
percentile is shown as the diffuse background only. The KIV-2 repeat is not seen by any SNP list; allele-specific short-read callers exist
(see ROADMAP). The total has no baseline: it is the amount added to a non-carrier, not an absolute level. A measured value always outranks it.
