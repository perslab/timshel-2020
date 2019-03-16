############### SYNOPSIS ###################
# EWCE generate cell-type data

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/nb-magma_celltyping"
setwd(wd)


library(tidyverse)
library(MAGMA.Celltyping) # also loads EWCE package

source("magma_celltyping_modified_functions.R")


# ======================================================================= #
# ====================  HYPOTHALAMUS SPECIFIC NOTES ==================== #
# ======================================================================= #

### Summary: for this dataset, we did not follow Skene's calculation of 'specificity' on non-CPM and non-log transformed data.
# 1) specificity was calculated on log(average): the input file "hypothalamus.celltype_expr.avg_expr.csv.gz" is already log-transformed
# 2) the cell-type expression 'average' was calculated from the *CPM normalized data* (@data slot) - see below
# TODO: calculate on non-log data by exp(..data..)

### THIS IS HOW THE the "hypothalamus.celltype_expr.avg_expr.csv.gz" file was calculated
### COPY FROM gene_expression_preprocessing.R
anno_specific_expr.avg_expr <- function(colname_ident, genes, seurat_obj) {
  # DESCRIPTION: calculate average gene expression across cell types. 
  # Returns average in log.space (this makes most sence, to dampen the gene expression variation)
  # RETURN: data frame with columnnames=annotations, rownames=genes 
  
  print(sprintf("Running average expression for colname_ident=%s", colname_ident))
  
  # AverageExpression(): uses `@data` slot (normalized and log-transformed single cell expression - NOT regressed)
  df.avg_expr <- log1p(AverageExpression(SetAllIdent(seurat_obj, id=colname_ident), genes.use=genes))
  # ---> RETURNS a DATA FRAME with columnnames=annotations, rownames=genes 
  
  return(df.avg_expr)
}




# ======================================================================= #
# ============================  LOAD DATA  =============================== #
# ======================================================================= #

file.data <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/hypothalamus/hypothalamus.celltype_expr.avg_expr.csv.gz" # file with MOUSE gene names
df.data_raw <- read_csv(file.data)
# df.data_raw <- df.data_raw %>% distinct(gene, .keep_all = TRUE) # remove duplicate gene names --> NO DUP
df.data_raw <- df.data_raw %>% column_to_rownames(var="gene")
head(df.data_raw)

# ======================================================================= #
# ============================  RUN  =============================== #
# ======================================================================= #

calculate.specificity.for.level <- function(df.mean_exp){
  normalised_meanExp = t(t(df.mean_exp)*(1/colSums(df.mean_exp)))
  specificity = normalised_meanExp/(apply(normalised_meanExp,1,sum)+0.000000000001)
  return(specificity)
}

ctd = list(list(
  "mean_exp"=as.matrix(df.data_raw), # matrix, genes x cells. with rownames and columnnames
  "annot"=colnames(df.data_raw),
  "specificity"=calculate.specificity.for.level(as.matrix(df.data_raw))
  ))

ctd
### STRUCTURE OF CTD FILES
# list[<data_set>] --> list of 3
# list[<data_set>]["mean_exp"] # data.frame. genes x cell_types
# list[<data_set>]["specificity"] # matrix genes x cell_types
# list[<data_set>]["annot"] # character. cell labels
# x <- ctd[[1]][["mean_exp"]]
# x <- ctd[[1]][["specificity"]]
# x <- ctd[[1]][["annot"]]

save(ctd,file="CellTypeData_hypothalamus_log_avg.rda")


### Load
# load("CellTypeData_mousebrain.rda") # loads ctd file
# str(ctd)

