"""Every external endpoint the pipeline depends on, in one place.

When an upstream site changes, edit the constant here, then run `wgsreport doctor` (checks that each endpoint still
answers and still has the format the parsers expect) and `pytest -m network`.
"""
# PGS Catalog harmonised scoring files (GRCh38 positions). Format contract: see pgs.read_score and doctor.py.
PGS_FTP = "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{id}/ScoringFiles/Harmonized/{id}_hmPOS_GRCh38.txt.gz"
PGS_REST = "https://www.pgscatalog.org/rest"
# ClinVar GRCh38 VCF (always the latest weekly release). Required INFO fields: see CLINVAR_INFO_FIELDS.
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
CLINVAR_INFO_FIELDS = ("CLNSIG", "CLNREVSTAT", "CLNDN", "GENEINFO", "MC", "AF_TGP", "RS")
# Ensembl GFF3, pinned to a release (the 'current_gff3' directory renames its file every release and would break).
ENSEMBL_RELEASE = 116
ENSEMBL_GFF_URL = ("https://ftp.ensembl.org/pub/release-{r}/gff3/homo_sapiens/"
                   "Homo_sapiens.GRCh38.{r}.chr.gff3.gz").format(r=ENSEMBL_RELEASE)
# EBI CRAM reference registry (sequence by MD5).
EBI_MD5_URL = "https://www.ebi.ac.uk/ena/cram/md5/{md5}"
CHRM_MD5 = "c68f52674c9fb33aef52dcf399755519"            # GRCh38 chrM, 16,569 bp (used by doctor.py)
# Public GWAS summary statistics used by `wgsreport build-score` (the user supplies the files).
PGC3_SCZ_FIGSHARE = "https://ndownloader.figshare.com/files/34517828"
# 1000 Genomes (the user downloads these; listed for documentation and for doctor.py).
KG_FTP = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
# Optional external tools (the user downloads them; nothing is bundled). See pgx_pharmcat.py and hla_typing.py.
PHARMCAT_RELEASE_URL = "https://github.com/PharmGKB/PharmCAT/releases"      # jar + 'pipeline' tarball; Apache-2.0
T1K_URL = "https://github.com/mourisl/T1K"                                    # HLA/KIR genotyper; GPL-3.0
IMGT_HLA_URL = "https://github.com/ANHIG/IMGTHLA"                              # allele database T1K downloads
# ClinGen gene-disease validity curations (CSV; used by scripts/build_cancer_genes.py to build data/cancer_predisposition_clingen.tsv).
CLINGEN_VALIDITY_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
# STRchive: curated pathogenic repeat thresholds (MIT license), used by scripts/build_repeat_thresholds.py. ExpansionHunter: Illumina, Apache-2.0.
STRCHIVE_LOCI_URL = "https://raw.githubusercontent.com/dashnowlab/STRchive/main/data/STRchive-loci.json"
EXPANSIONHUNTER_URL = "https://github.com/Illumina/ExpansionHunter/releases"
YLEAF_URL = "https://github.com/genid/Yleaf/releases"                           # Y haplogroup caller, standalone Linux build; MIT
