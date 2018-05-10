############### SYNOPSIS ###################
# Join GWAS SNPs with baseline annotations

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

############################################



# ======================================================================================================= #
# ==========================================       SETUP      =========================================== #
# ======================================================================================================= #


# suppressMessages(library(rolypoly))
suppressMessages(library(tidyverse))
library(stringr)

dir.sc_genetics.data <- "/projects/timshel/sc-genetics/sc-genetics/data" # no trailing slash
dir.ldfiles <- file.path(dir.sc_genetics.data, "rolypoly/EUR_LD_FILTERED_NONAN_R") # Linkage disequilibrium files
# /projects/timshel/sc-genetics/sc-genetics/data/rolypoly/EUR_LD_FILTERED_NONAN_R

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


# ======================================================================================================= #
#=============================================== READ GWAS ================================================== #
# ======================================================================================================= #

# file.gwas <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_WHRadjBMIz.gwassumstats.rolypoly_fmt.tab.gz" ### body_WHRadjBMIz
# file.gwas <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMIz.gwassumstats.rolypoly_fmt.tab.gz" ### body_BMIz [MISSING]
# file.gwas <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/disease_T2D.gwassumstats.rolypoly_fmt.tab.gz" ### disease_T2D.gwassumstats
# file.gwas <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_HEIGHTz.gwassumstats.rolypoly_fmt.tab.gz" ### body_HEIGHTz
# file.gwas <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/mental_SCZ_Ripke2014.gwassumstats.rolypoly_fmt.tab.gz" ### mental_SCZ_Ripke2014
# file.gwas <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/lipids_TC.Willer2013.gwassumstats.rolypoly_fmt.tab.gz" ### lipids_TC.Willer2013
file.gwas <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz" ### body_BMI_Locke2015

df.gwas <- read_tsv(file.gwas)
names(df.gwas) # "rsID"    "beta"    "se"      "pval"    "snp_maf" "chr"     "pos"

# ======================================================================================================= #
#===================================== READ SNPsnap collection file ===================================== #
# ======================================================================================================= #

# 
# file.snpsnap <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/snp_mapping/snpsnap_EUR_1KG_phase3-chrpos_mapping.tab.gz"

### XXXX

# ======================================================================================================= #
#=============================================== READ ROLYPOLY LD data SNPs ================================================== #
# ======================================================================================================= #

df.rp_ld_snps <- read_csv("/raid5/projects/timshel/sc-genetics/sc-genetics/data/rolypoly/rolypoly_ld_snps.csv.gz")

# ======================================================================================================= #
#=============================================== READ LDSC Annotation ================================================== #
# ======================================================================================================= #



# ======================================================================================================= #
#=============================================== XXXXX ================================================== #
# ======================================================================================================= #



# ======================================================================================================= #
#=============================================== XXXXX ================================================== #
# ======================================================================================================= #







