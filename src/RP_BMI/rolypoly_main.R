





library(rolypoly)
library(tidyverse)
library(ggplot2)

#' 
## ------------------------------------------------------------------------
wd <- "/projects/timshel/sc-genetics/sc-genetics/src"
setwd(wd)

#' 
#' 
#' # PARAMETERS
#' 
## ------------------------------------------------------------------------
dir.ldfiles <- "../data/rolypoly/EUR_LD_FILTERED_NONAN_R" # Linkage disequilibrium files

#' 
#' 
#' # GWAS
#' 
#' REQUIRED COLUMNS IN DATA FRAME: c('chrom', 'pos', 'rsid', 'beta', 'se', 'maf')
#' 
## ------------------------------------------------------------------------
file.gwas <- "../data/gwas_sumstats/body_BMI_Locke2015.gwassumstats.rolypoly_fmt.tab.gz"
df.gwas <- read_delim(file.gwas, delim="\t") # suppressMessages or use "col_types=cols()"

# TODO: check that chr is always integer

#' 
#' Rename columns
## ------------------------------------------------------------------------
df.gwas.rp <- df.gwas %>% 
  rename(
    rsid = rsID,
    chrom = chr,
    maf = snp_maf) 
  
# must be data frame for RolyPoly to run

df.gwas.rp %>% head

#' 
#' 
#' 
#' # Gene annotation
#' 
## ------------------------------------------------------------------------
file.gene_annot <- "../data/gene_annotations/gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz"
df.gene_annot <- read_delim(file.gene_annot, col_types=cols(), delim="\t") # suppressMessages or use "col_types=cols()"

#' 
#' 
## ------------------------------------------------------------------------
df.gene_annot %>% head

#' 
#' ## Filter genes: filter out non-autosomal genes (discard X, Y, contigs etc)
## ------------------------------------------------------------------------
n.discarded <- df.gene_annot %>% filter(!chromosome_name %in% 1:22) %>% nrow()
print(sprintf("Number of non-autosomal genes discarded: %s", n.discarded))

df.gene_annot.filtered <- df.gene_annot %>% filter(chromosome_name %in% 1:22)

#' 
#' 
#' # Expression data
#' 
#' REQUIRED COLUMNS IN DATA FRAME: none (just cell-type/tissue names); rownames is genes
#' 
## ------------------------------------------------------------------------
file.expr <- "../data/expression/arc_lira/arc_lira-celltype_expr_ttest.hsapiens_orthologs.csv.gz"
df.expr <- read.csv(file.expr)

#' 
## ------------------------------------------------------------------------
df.expr %>% head

#' 
#' 
#' ## Filter genes (based on `df.gene_annot.filtered`, e.g. autosomal)
#' 
#' 
## ------------------------------------------------------------------------
bool.expr_genes_filtered <- rownames(df.expr) %in% df.gene_annot.filtered$ensembl_gene_id
df.expr.filtered <- df.expr[bool.expr_genes_filtered, ]

#' 
#' 
#' 
#' 
#' # Numerical transformation of expression data (square values)
#' 
#' 1. [only for non-ttest value]: Z-score (R scale)
#'     1. z-scores for each gene (across cells) "cell’s relative expression of each gene, in comparison to all other cells:
#' 2. Square to avoid negative values.
#' 
## ------------------------------------------------------------------------
### Square values
df.expr.filtered.rp <- df.expr.filtered^2

#' 
#' 
#' 
#' # Block annotation: window size, etc
#' 
#' REQUIRED COLUMNS IN DATA FRAME: c('chrom', 'start', 'end', 'label')
#' 
#' Set the parameter for the "gene window" size of the gene annotation by dplyr::transmute() the column names:
#'   1. label=ensembl_gene_id HUMAN
#'   2. chrom=chr
#'   3. start={gene_start-10kb; gene_start-10kb}
#'   4. end={gene_start+10kb; gene_end-10kb}
#' 
## ------------------------------------------------------------------------
WINDOW_SIZE <- 10*10^3 # 10kb

df.block_annotation <- df.gene_annot.filtered %>% transmute(
  label=ensembl_gene_id,
  chrom=chromosome_name,
  start=start_position,
  end=end_position
)

df.block_annotation %>% head

#' 
#' 
#' # Rolling rolypoly
#' 
#' We include all the previously described data into the main rolypoly function call. This function has many parameters to tinker with, however, should run fine with the defaults. Most importantly consider the number of bootstrap iterations to get accurate standard errors. We usually use at least 200.
#' 
## ------------------------------------------------------------------------

df.gwas.rp.test <- df.gwas.rp %>% slice(1:1000)



#' 
#' 
#' 
#' 
#' 
#' 
## ------------------------------------------------------------------------
library(doParallel)
registerDoParallel(40) # half of the number of cores 
getDoParWorkers()

#' 
#' 

#' 
#' 
#' 
## ---- results='hide'-----------------------------------------------------
# message = FALSE

# TODO: run with "highres" data

rm(rp)
x_timing <- system.time(rp <- rolypoly_roll(
  gwas_data = as.data.frame(df.gwas.rp),
  block_annotation = as.data.frame(df.block_annotation),
  block_data = as.data.frame (df.expr.filtered.rp),
  ld_folder = dir.ldfiles,
  bootstrap_iters = 50,
  gwas_link_parallel = T,
  bootstrap_parallel = T
))

print(x_timing)

save.image(file="out.rolypoly_main-SCRIPT.RData")


