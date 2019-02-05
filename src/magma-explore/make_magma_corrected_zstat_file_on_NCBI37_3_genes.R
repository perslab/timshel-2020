############### SYNOPSIS ###################
# Overall aim: create a MAGMA gene-based scores (.genes.out) file containing ALL MAGMA genes
# Step 1: create a 'dummy' covar file containing ALL MAGMA genes.
# Step 2: run magma with '--model correct=all' to get corrected ZSTAT
# We must first create a 'dummy' covar file with all genes because MAGMA only output gene scores for genes in the covar file

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(tidyverse)

setwd("/raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/")

# ======================================================================= #
# =========================== Make dummy covar file ===================== #
# ======================================================================= #

file.all_magma_genes <- "/tools/magma/1.07a/gene_loc.NCBI37.3/NCBI37.3.gene.loc" # get all genes in MAGMA universe (these genes are used for the annotation of GWASs)
### SNIPPET: file has no header
### First column is entrez_id
# 79501   1       69091   70008   +       OR4F5
# 100996442       1       142447  174392  -       LOC100996442
df.all_magma_genes <- read_tsv(file.all_magma_genes, col_names=F, col_types=cols_only(X1=col_character())) %>% rename(entrez_id=X1)
df.all_magma_genes <- df.all_magma_genes %>% mutate(dummy_covar=rnorm(n=n())) # variable must be 'realistic' otherwise MAGMA will throw it away.
df.all_magma_genes

### Write
# write_tsv(df.all_magma_genes, path="magma_dummy_covar_file.NCBI37_3.tab")

# ======================================================================= #
# =========================== Make MAGMA call =========================== #
# ======================================================================= #

### BMI gwas
# magma \
# --gene-results /scratch/tmp-magma_gwas/BMI_Yengo2018.txt.10UP.1.5DOWN.genes.raw \
# --gene-covar /raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/magma_dummy_covar_file.NCBI37_3.tab \
# --model correct=all direction-covar=greater \
# --settings abbreviate=0 gene-info \
# --out /raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/BMI_Yengo2018.resid.correct_all


