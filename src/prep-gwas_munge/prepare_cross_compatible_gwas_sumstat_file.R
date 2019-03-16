############### SYNOPSIS ###################
# Prepare cross compatible GWAS summary stats file
# RolyPoly
# MAGMA
# DEPICT


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

### DESIGN DESCRIPTION:
# 1. read LDSC munged BMI sumstats (Hapmap3 SNPs only)
# 2. read SNPsnap mapping file
# 3. read BMI original GWAS file (needed for beta / se for RolyPoly)
# 4. add extra data cols to LDSC (merge data)

# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)

setwd(here("src/analysis-compare_genetic_prioritization_methods/"))


# gwas_name <- "BMI_UPDATE_Yengo2018"
gwas_name <- "BMI_UKBB_Loh2018"

# ======================================================================= #
# =============================== READ DATA ================================= #
# ======================================================================= #

if (gwas_name == "BMI_UPDATE_Yengo2018") {
  ### BMI_UPDATE_Yengo2018.sumstats.gz
  file.gwas_munge <- here("data/gwas_sumstats_ldsc/timshel-collection/BMI_UPDATE_Yengo2018.sumstats.gz")
  df.gwas_munge <- read_tsv(file.gwas_munge)
  file.gwas_orig <- here("data/gwas_sumstats_raw/BMI_HEIGHT_Yengo2018/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz")
  df.gwas_orig <- read_tsv(file.gwas_orig)
  df.gwas_orig <- df.gwas_orig %>% select(SNP, BETA, SE, P_ORIG=P)
  #### CHR     POS     SNP     Tested_Allele   Other_Allele    Freq_Tested_Allele_in_HRS       BETA    SE      P       N
  #### 7       92383888        rs10    A       C       0.06431  0.0013 0.0042   7.5e-01        598895
  #### 12      126890980       rs1000000       A       G       0.2219   0.0001 0.0021   9.6e-01        689928
} else if (gwas_name == "BMI_UKBB_Loh2018") {
  ### BMI_UKBB_Loh2018.sumstats.gz
  file.gwas_munge <- here("data/gwas_sumstats_ldsc/timshel-collection/BMI_UKBB_Loh2018.sumstats.gz")
  df.gwas_munge <- read_tsv(file.gwas_munge)
  file.gwas_orig <- here("data/gwas_sumstats_raw/alkesgroup_UKBB/UKBB_460k/body_BMIz.sumstats.gz")
  df.gwas_orig <- read_tsv(file.gwas_orig)
  df.gwas_orig <- df.gwas_orig %>% select(SNP, BETA=Beta, SE=se, P_ORIG=P)
  ## SNP     CHR     POS     A1      A2      REF     EAF     Beta    se      P       N       INFO
  ## rs10399793      1       49298   T       C       T       0.376191        0.00267751      0.00328604      3.1E-01 445921  0.342797
  ## rs2462492       1       54676   C       T       C       0.59937 0.000666627     0.00325463      8.5E-01 445921  0.340158
} else {
  stop("Wrong gwas_name")
}

file.snp_mapping <- here("data/snp_mapping/snpsnap_EUR_1KG_phase3-chrpos_mapping.tab.gz")
df.snp_mapping <- read_tsv(file.snp_mapping, col_types=cols(rsID="c", .default = col_guess())) # avoid read_tsv parsing rsID column as "col_time(format = "")"
bool.isin <- df.gwas_munge$SNP %in% df.snp_mapping$rsID
sum(bool.isin) # [1] 1177617
length(bool.isin) - sum(bool.isin) # missing = 39,694


# ======================================================================= #
# ======================= MERGE DATA (add cols) ========================= #
# ======================================================================= #

### drop all SNPs that have missing values (i.e HapMap3 SNPs not found in the GWAS summary stats)
df.gwas <- df.gwas_munge %>% drop_na()
nrow(df.gwas) # Yengo: 1019865

### add chr/pos
df.gwas <- df.gwas %>% left_join(df.snp_mapping, by=c("SNP"="rsID"))
df.gwas <- df.gwas %>% drop_na()
nrow(df.gwas) # Yengo: 1011703

### add beta/se/pval_orig
df.gwas <- df.gwas %>% left_join(df.gwas_orig, by="SNP")

### drop na and remove any duplicates (there were 2 duplicated rsIDs for Yengo)
df.gwas <- df.gwas %>% drop_na()
df.gwas <- df.gwas %>% filter(! duplicated(.$SNP))

### calculate P from Z
df.gwas <- df.gwas %>% mutate(P = pnorm(abs(Z), lower.tail=F)*2) # two-sided P-value as in GWAS

# ======================================================================= #
# ======================== Check for outliers and MHC ==================== #
# ======================================================================= #

### Report meta statistics:
df.gwas %>% summarise(N_max = max(N),
                      N_avg = mean(N),
                      Z_outlier_80 = sum(Z^2 > max(80, 0.0001*N)), # SNPs with unusual high chisq statistics 
                      Z_outlier_180 = sum(Z^2 > max(180, 0.0001*N)), # SNPs with unusual high chisq statistics 
                      Z_outlier_360 = sum(Z^2 > max(360, 0.0001*N)), # SNPs with unusual high chisq statistics 
                      Z_outlier_720 = sum(Z^2 > max(720, 0.0001*N)), # SNPs with unusual high chisq statistics 
                      N_MHC = sum((chr==6 & pos > 25e6) & (chr==6 & pos < 34e6)) # SNPs in MHC region [chr6:25Mb-34Mb]
                      ) %>% t()
### Yengo
# N_max         795640.0
# N_avg         686252.4
# Z_outlier_80    1980.0
# Z_outlier_180    373.0
# Z_outlier_360     79.0
# Z_outlier_720     11.0
# N_MHC           2380.0

### Lo
# N_max         457824
# N_avg         457824
# Z_outlier_80    1634
# Z_outlier_180    266
# Z_outlier_360     70
# Z_outlier_720     17
# N_MHC           2827

# Outlier : 16 53800387 --> FTO genes


# ======================================================================= #
# ============================= WRITE OUTPUT ============================= #
# ======================================================================= #


### filter out mhc
df.gwas %>% 
  filter(! ((chr==6 & pos > 25e6) & (chr==6 & pos < 34e6)) ) %>%
  write_tsv(sprintf("%s_no_mhc.sumstats.gz", gwas_name))


### filter out mhc + z_outliers
chisq_cutoffs <- c(80, 720)
for (chisq_cutoff in chisq_cutoffs) {
  print(chisq_cutoff)
  df.gwas %>% 
    filter(! ((Z^2 > max(chisq_cutoff, 0.0001*N)) | ((chr==6 & pos > 25e6) & (chr==6 & pos < 34e6))) ) %>%
    write_tsv(sprintf("%s_no_mhc_max_chisq_%s.sumstats.gz", gwas_name, chisq_cutoff))
}


# ======================================================================= #
# ========== APPENDIX: number of SNPs not mapped with chr and pos ====== #
# ======================================================================= #

### CONLCUSION
# We munged GWAS file include 1,217,311 HapMap3 SNPs.
# SNPsnap --->  missing = 39,694
# LDSC 1kg --> missing = 26,990
# Both mapping files result in more or the same mapping rate.

# ### SNPSNAP
# file.snp_mapping <- here("data/snp_mapping/snpsnap_EUR_1KG_phase3-chrpos_mapping.tab.gz")
# df.snp_mapping <- read_tsv(file.snp_mapping, col_types=cols(rsID="c", .default = col_guess())) # avoid read_tsv parsing rsID column as "col_time(format = "")"
# 
# bool.isin <- df.gwas_munge$SNP %in% df.snp_mapping$rsID
# sum(bool.isin) # [1] 1177617
# length(bool.isin) - sum(bool.isin) # missing = 39,694
# 
# ### LDSC 1KG data
# file.snp_mapping_ldsc_1kg <- "/raid5/projects/timshel/sc-genetics/ldsc/data/1000G_EUR_Phase3_plink/1000G.EUR.QC.CHR_1_22.bim"
# df.snp_mapping_ldsc_1kg <- read_tsv(file.snp_mapping_ldsc_1kg, col_names=c("chr", "rsID", "cm", "pos", "a1", "a2"))
# 
# bool.isin <- df.gwas_munge$SNP %in% df.snp_mapping_ldsc_1kg$rsID
# sum(bool.isin) # [1] 1190321
# length(bool.isin) - sum(bool.isin) # missing = 26,990

# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #




