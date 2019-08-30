############### SYNOPSIS ###################
# Overall aim: create a MAGMA gene-based scores (.genes.out) file containing ALL MAGMA genes
# Step 0: run magma annotate to create .genes.annot file for a given window size
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
library(here)


source(here("src/lib/load_functions.R")) # load sc-genetics library


dir.wd <- "src/magma-gene_scores"
setwd(here(dir.wd))

# ======================================================================= #
# ================================ CONSTANTS ================================ #
# ======================================================================= #

RUN_MODE <- TRUE # if True, then MAGMA system calls will be made

magma_exec <- "/tools/magma/1.07a/magma"
dir.out_magma <- here(dir.wd, "/out_magma_BMI")
dir.create(dir.out_magma, showWarnings = FALSE)

# ======================================================================= #
# ============== [*RUN ONCE*] Annotate genes: create .genes.annot ======= #
# ======================================================================= #

### CMD
# magma --annotate window=<KB_upstream>,<KB_downstream> --snp-loc [SNPLOC_FILE] --gene-loc [GENELOC_FILE] --out [ANNOT_PREFIX]
# NB: window modifier must come directly after --annotate

# This will read SNP locations from the file [SNPLOC_FILE] and gene locations
# from the file [GENELOC_FILE], and produces the file [ANNOT_PREFIX].genes.annot
# containing the mapping of SNPs to genes. The [GENELOC_FILE] can be downloaded
# from the MAGMA website for different builds, the [SNPLOC_FILE] must be
# provided by the user. The .bim file of the genotype (reference) data can be used for this.
# NB: make sure that the gene and SNP locations are based on the same human genome reference build

### --snp-loc
# Requires a formatted data file, with rows corresponding to SNPs; a header is not allowed.
# The file must have three columns containing the SNP ID, chromosome and base pair position, in that order. 
# Allowed chromosome codes are 1-24, X and Y (where 23 and 24 correspond to X and Y respectively).
# Alternatively, a .bim file can be used as well.

### Params
file.gene_loc <- "/tools/magma/1.07a/gene_loc.NCBI37.3/NCBI37.3.gene.loc" # hg19
file.snp_loc <- "/tools/magma/1.07a/g1000_eur/g1000_eur.bim" # h19
file.annot_prefix <- file.path(dir.out_magma, "NCBI37_1kgp_up100kb_down100kb")

### Call
cmd_args <- sprintf("--annotate window=100,100 --snp-loc %s --gene-loc %s --out %s", file.snp_loc, file.gene_loc, file.annot_prefix)
cmd_args
if (RUN_MODE) {
  system2(magma_exec, cmd_args)
}


# ======================================================================= #
# ================= [*RUN ONCE*] Make dummy covar file ================= #
# ======================================================================= #

### file.gene_loc | get all genes in MAGMA universe (these genes are used for the annotation of GWASs)
### SNIPPET: file has no header
### First column is entrez_id
# 79501   1       69091   70008   +       OR4F5
# 100996442       1       142447  174392  -       LOC100996442
df.all_magma_genes <- read_tsv(file.gene_loc, col_names=F, col_types=cols_only(X1=col_character())) %>% rename(entrez_id=X1)
df.all_magma_genes <- df.all_magma_genes %>% mutate(dummy_gene_covar=rnorm(n=n())) # variable must be 'realistic' otherwise MAGMA will throw it away.
df.all_magma_genes

### Write
file.dummy_gene_covar <- file.path(dir.out_magma, "magma_dummy_gene_covar_file.NCBI37_3.tab")
write_tsv(df.all_magma_genes, path=file.dummy_gene_covar)


# ======================================================================= #
# ======================== [*RUN PER GWAS*] ============================= #
# ======================================================================= #

### Params
gwas_name <- "BMI_UKBB_Loh2018_no_mhc"
# gwas_name <- "BMI_UPDATE_Yengo2018_no_mhc"

# ======================================================================= #
# == Calc gene-based P-vals: create .genes.raw and .genes.out (uncorrected) == #
# ======================================================================= #
### NB: the .genes.out from this command are UNCORRECTED P-vals
# We can only obtain CORRECTED P-vals using --gene-covar/--model arguments
# However, we need the .genes.raw output files from this command to run the --gene-covar/--model command.

### CMD
# magma --bfile [DATA] --gene-annot [ANNOT].genes.annot --pval [PVAL_FILE] ncol=[N_COL] --out [GENE_PREFIX]

file.gwas <- file.path(here("/data/gwas_sumstats_ldsc/timshel-collection"), sprintf("%s.sumstats", gwas_name)) # MUST be an uncompressed file with "SNP" and "P" columns (by default).
file.annot <- sprintf("%s.genes.annot", file.annot_prefix) # add suffix
file.bfile <- "/tools/magma/1.07a/g1000_eur/g1000_eur" # no file extension, just prefix
file.out_genes_prefix <- file.path(dir.out_magma, gwas_name)

### Call
cmd_args <- sprintf("--bfile %s --gene-annot %s --pval %s ncol=N --out %s", file.bfile, file.annot, file.gwas, file.out_genes_prefix)
cat("magma", cmd_args)
if (RUN_MODE) {
  system2(magma_exec, cmd_args)
}

### OUTPUT files
# .genes.raw
# .genes.out

# ======================================================================= #
# ====== MAGMA covar: create gsa.genes.out [corrected gene-based vals] ======== #
# ======================================================================= #

### CMD
# magma --gene-results [GENE_RESULTS].genes.raw --gene-covar [COVAR_FILE] --out [OUTPUT_PREFIX]

### Params
file.gene_results <- sprintf("%s.genes.raw", file.out_genes_prefix)
file.out_gene_level_prefix <- file.path(dir.out_magma, sprintf("%s.resid_correct_all", gwas_name))
# file.dummy_gene_covar

### Call [correct=all]
cmd_args <- sprintf("--gene-results %s --gene-covar %s --model correct=all direction-covar=greater --settings abbreviate=0 gene-info --out %s", file.gene_results, file.dummy_gene_covar, file.out_gene_level_prefix)
cat("magma", cmd_args)
if (RUN_MODE) {
  system2(magma_exec, cmd_args)
}

### OUTPUT files
# .gsa.genes.out
# .gsa.out
# .log

# ======================================================================= #
# ============= Add more gene IDs: Enesbl + gene symbol IDs ============ #
# ======================================================================= #

### Read file from previous step
file.gsa_genes <- sprintf("%s.gsa.genes.out", file.out_gene_level_prefix) # output file from MAGMA --gene-covar --model cmd
df.magma <- read_table(file.gsa_genes, comment = "#")

### add ensembl ids
df.magma <- add_ensembl_ids_from_entrez(df.magma, colname_geneids_from="GENE", colname_geneids_to="ensembl_gene_id")
# [1] "Number of genes mapped: 17951"
# [1] "Number of genes not mapped: 193"
# [1] "Number of genes with a NON-unique mapping (genes with duplicated ensembl gene IDs after mapping): 41"
# [1] "Total mapping stats: 234 genes have no mapping (not mapped + duplicates) out of 18144 input genes."
# [1] "Total genes mapped (non NA genes): 17910"
# [1] "Returning tibble with the column 'ensembl_gene_id' added where all gene identifiers unique. Unmapped genes have NA values"

### gene symbol
df.magma <- hs_add_gene_symbol_from_ensembl_ids(df.magma, colname_geneids_from="ensembl_gene_id", colname_geneids_to="gene_symbol")
# [1] "Number of genes mapped: 17481"
# [1] "Number of genes not mapped: 663"
# [1] "Number of genes with a NON-unique mapping (genes with duplicated ensembl gene IDs after mapping): 0"
# [1] "Total mapping stats: 663 genes have no mapping (not mapped + duplicates) out of 18144 input genes."
# [1] "Total genes mapped (non NA genes): 17481"
# [1] "Returning tibble with the column 'gene_symbol' added where all gene identifiers unique. Unmapped genes have NA values"

# ======================================================================= #
# =========================== Add mouse ortholog gene =================== #
# ======================================================================= #

### add mouse orthologs
df.magma <- human_to_mouse_ortholog_gene_expression_mapping(df.magma, colname_geneids_from="ensembl_gene_id", colname_geneids_to="mmusculus_homolog_ensembl_gene")
# [1] "Number of genes with a NON-unique mapping (genes with duplicated gene IDs after mapping): 0"
# [1] "Number of genes mapped: 14972"
# [1] "Number of genes not mapped: 3172"


# ======================================================================= #
# =========================== Write output file =================== #
# ======================================================================= #

file.export <- sprintf("%s.gsa.genes.mapped.out", file.out_gene_level_prefix)
df.magma %>% write_tsv(file.export)

