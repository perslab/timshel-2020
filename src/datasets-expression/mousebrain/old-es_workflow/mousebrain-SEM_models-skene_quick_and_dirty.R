############### SYNOPSIS ###################
# *OBS* these SEMs use data from "mousebrain-agg_L5.csv.gz" ('CellTypeData_mousebrain.rda') that is *NOT* CPM and log-transformed as my default 'Seurat workflow' does
# 1: Export Mousebrain enrichment scores
# 2: Export Skene 'CTD' data: mean, specificity and rank



### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library


wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain"
setwd(wd)


# ======================================================================= #
# ==============================  FUNCTION  =============================== #
# ======================================================================= #

transform_to_quantiles <- function(vector_in, threshold_bin_zero=0){
  # Function to compute Skene quantile bins 
  # INPUT: a vector of numeric values to bin
  # OUTPUT: a numeric vector with values {0 ... numberOfBins-1?}
  # *REMARK*: note that vector_in < threshold_bin_zero will be put in bin 0.
  
  numberOfBins=41 # constant
  quantileValues = rep(0,length(vector_in))
  quantileValues[vector_in>threshold_bin_zero] = as.numeric(cut(vector_in[vector_in>threshold_bin_zero], # cut divides the range of x into intervals and codes the values in x according to which interval they fall. 
                                                 breaks=unique(quantile(vector_in[vector_in>threshold_bin_zero], probs=seq(0,1, by=1/numberOfBins), na.rm=TRUE)), 
                                                 include.lowest=TRUE))
  ### ORIGINAL
  # quantileValues = rep(0,length(vector_in))
  # quantileValues[vector_in>0] = as.numeric(cut(vector_in[vector_in>0], # cut divides the range of x into intervals and codes the values in x according to which interval they fall. 
  #                                              breaks=unique(quantile(vector_in[vector_in>0], probs=seq(0,1, by=1/numberOfBins), na.rm=TRUE)), 
  #                                              include.lowest=TRUE))
  
  return(quantileValues)
}

# ======================================================================= #
# ==============================  ***SKENE CTD EXPORT**  =============================== #
# ======================================================================= #

library(MAGMA.Celltyping) # loads EWCE package
library(tidyverse)

### Load the celltype data
# DATA_SET <- "MACA"
DATA_SET <- "mousebrain"
load(sprintf("/projects/timshel/sc-genetics/sc-genetics/src/nb-magma_celltyping/CellTypeData_%s.rda", DATA_SET)) # loads ctd file

### COPY FROM ewce-mousebrain.R This 'CellTypeData_mousebrain.rda' was loaded like this from the "mousebrain-agg_L5.csv.gz"
# file.data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-agg_L5.csv.gz"
# df.data_raw <- read_csv(file.data) %>% rename(gene=X1)
# # ^^ duplicates: ‘3110039M20Rik’, ‘4930556M19Rik’, ‘Apoc2’, ‘Atn1’, ‘C1s2’, ‘Ccdc142’, ‘Ccl19’, ‘Ccl21a’, ‘Ccl21b’, ‘Ccl21c’, ‘Ccl27a’, ‘Cd37’, ‘D130017N08Rik’, ‘Dancr’, ‘Fam205a2’, ‘Fbxw14’, ‘Flg’, ‘Gbp6’, ‘Gm15853’, ‘Gm16701’, ‘Gm2464’, ‘Gm3286’, ‘Hist2h2bb’, ‘Il11ra2’, ‘Itgam’, ‘l7Rn6’, ‘Ltbp4’, ‘Map2k7’, ‘Nova2’, ‘Ntn5’, ‘Olfr108’, ‘Olfr126’, ‘Olfr1284’, ‘Olfr1309’, ‘Olfr1316’, ‘Olfr1366’, ‘Olfr1396’, ‘Olfr1496’, ‘Olfr170’, ‘Olfr730’, ‘Olfr790’, ‘Olfr809’, ‘Pcdha11’, ‘Pcdhga8’, ‘Pik3c2g’, ‘Rp1’, ‘Schip1’, ‘Sgsm3’, ‘Smim20’, ‘Syngr4’, ‘Tead2’, ‘Tgfb1i1’, ‘Tulp2’, ‘U2af1l4’, ‘Umad1’, ‘Zfand4’ 
# df.data_raw <- df.data_raw %>% distinct(gene, .keep_all = TRUE) # remove duplicate gene names
# df.data_raw <- df.data_raw %>% column_to_rownames(var="gene")
# head(df.data_raw)


ctd = prepare.quantile.groups(ctd,specificity_species="mouse",numberOfBins=40) # contains mouse gene symbols (filtered on 1:1 orthologs)

### INFO
# the ctd data has been FILTERED on 1:1 orthologs by the EWCE package.
# the ctd data contains MOUSE GENE IDENTIFYERS, and hence has not been mapped yet.
# SEE "/raid5/projects/timshel/sc-genetics/sc-genetics/src/nb-magma_celltyping/ewce-mousebrain" for details on how the CTD RData file was generated.

### INFO #2:
# mouse_to_human_ortholog_gene_expression_mapping() MUST get as input a data frame with columnnames=annotations, rownames=mouse_gene_symbols

### mean_exp
df <- ctd[[1]][["mean_exp"]] %>% 
  as.data.frame %>%
  rename_all(funs(sapply(stringr::str_split(., "_"), '[[', 2))) # rename: drop the prefix for the ClusterName, e.g. drop 'Astrocytes_*'
write_csv(df %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.avg_expr.csv.gz")
df.human <- mouse_to_human_ortholog_gene_expression_mapping(df)
write_csv(df.human %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz")

### mean_exp --> z-score ---> [+MORE] --> quantiles
df <- ctd[[1]][["mean_exp"]] %>% 
  as.data.frame %>%
  rename_all(funs(sapply(stringr::str_split(., "_"), '[[', 2))) # rename: drop the prefix for the ClusterName, e.g. drop 'Astrocytes_*'
df.human <- mouse_to_human_ortholog_gene_expression_mapping(df)
df.expr <- df.human
df.expr.t <- t(as.data.frame(df.expr)) # transposing. returns a DATA FRAME with columns=genes, rownames=annotations
df.expr.t.scale <- scale(df.expr.t, center = TRUE, scale = TRUE) # returns a *MATRIX* with colnumns=genes and rownames=annotations
df.expr <- as.data.frame(t(df.expr.t.scale)) # transposing back. returns a DATA FRAME with columns=annotations, rownames=genes.
head(df.expr)
df.expr.na <- df.expr %>% filter_all(any_vars(is.na(.))) # we have NA values. Don't know why
df.expr.na
df.expr_no_na <- df.expr[complete.cases(df.expr), ] # *OBS* remove NA: they appently exists
write_csv(df.expr_no_na %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.z_score.hsapiens_orthologs.csv.gz")

# sqaure
df.expr_no_na.sq <- as.data.frame(df.expr_no_na^2) # *OBS*: x^2 or x**2 returns a matrix!
write_csv(df.expr_no_na.sq %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.z_score_sq.hsapiens_orthologs.csv.gz")
# log10(square+1)
df.expr_no_na.sq_log <- log10(df.expr_no_na.sq+1)
head(df.expr_no_na.sq_log)
write_csv(df.expr_no_na.sq_log %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.z_score_sq_log.hsapiens_orthologs.csv.gz")
# abs
df.expr_no_na.abs <- abs(df.expr_no_na) # returns data.frame
head(df.expr_no_na.abs)
write_csv(df.expr_no_na.abs %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.z_score_abs.hsapiens_orthologs.csv.gz")
# log10(abs+1)
df.expr_no_na.abs_log <- log10(df.expr_no_na.abs+1) 
head(df.expr_no_na.abs_log)
write_csv(df.expr_no_na.abs_log %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.z_score_abs_log.hsapiens_orthologs.csv.gz")
# pos_only
df.expr_no_na.pos_only <- replace(df.expr_no_na, df.expr_no_na < 0, 0) # replace negative values with 0
head(df.expr_no_na.pos_only)
write_csv(df.expr_no_na.pos_only %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.z_score_pos.hsapiens_orthologs.csv.gz")
# log10(pos_only+1)
df.expr_no_na.pos_only_log <- log10(df.expr_no_na.pos_only+1)
head(df.expr_no_na.pos_only_log)
write_csv(df.expr_no_na.pos_only_log %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.z_score_pos_log.hsapiens_orthologs.csv.gz")


df.expr_no_na.rank <- df.expr_no_na %>% transmute_all(funs(transform_to_quantiles), threshold_bin_zero=-9999) # *OBS* output tibble looses rownames
  # ^ We set threshold_bin_zero=-9999 because ~1/2 the z-scores will be negative.
rownames(df.expr_no_na.rank) <- rownames(df.expr_no_na) # set rownmaes *IMPORTANT*
head(df.expr_no_na.rank)
write_csv(df.expr_no_na.rank %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.z_score_quantiles.hsapiens_orthologs.csv.gz")



### specificity + quantiles
df <- ctd[[1]][["specificity"]] %>% 
  as.data.frame %>%
  rename_all(funs(sapply(stringr::str_split(., "_"), '[[', 2))) # rename: drop the prefix for the ClusterName, e.g. drop 'Astrocytes_*'
write_csv(df %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.specificity.csv.gz")
df.human <- mouse_to_human_ortholog_gene_expression_mapping(df)
write_csv(df.human %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.specificity.hsapiens_orthologs.csv.gz")
df.human.rank <- df.human %>% transmute_all(funs(transform_to_quantiles), threshold_bin_zero=0) # *OBS* output tibble looses rownames
rownames(df.human.rank) <- rownames(df.human) # set rownmaes *IMPORTANT*
write_csv(df.human.rank %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.specificity_quantiles.hsapiens_orthologs.csv.gz")


### quantiles
df <- ctd[[1]][["quantiles"]] %>% 
  as.data.frame %>%
  rename_all(funs(sapply(stringr::str_split(., "_"), '[[', 2))) # rename: drop the prefix for the ClusterName, e.g. drop 'Astrocytes_*'
write_csv(df %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.quantiles.csv.gz")
df.human <- mouse_to_human_ortholog_gene_expression_mapping(df)
write_csv(df.human %>% rownames_to_column(var="gene"), "mousebrain_skene.celltype_expr.quantiles.hsapiens_orthologs.csv.gz")

