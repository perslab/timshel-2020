############### SYNOPSIS ###################
# Compute cell-type specific expression for GTEx data (no ortholog gene mapping)

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ==============================  USAGE  =============================== #
# ======================================================================= #

# time Rscript gtex-cell_type_specific_expression.R |& tee gtex-cell_type_specific_expression.out.txt

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)
library(here)

library(Seurat)
library(Matrix)

library(parallel)


### Source custom scripts
dir.project_src <- "/projects/timshel/sc-arc_lira/src"
dir.pers_lab_sc_lib <- "/projects/timshel/git/perslab-sc-library"
source(sprintf("%s/constants-cell_type_annotations.R", dir.project_src)) # loads cell type annotations
source(sprintf("%s/seurat_functions/load_functions.R", dir.pers_lab_sc_lib)) # load Pers lab/Timshel single-cell library

source(here("src/lib/load_functions.R")) # load sc-genetics library


# ======================================================================= #
# ============================  PARAMS  ============================== #
# ======================================================================= #

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-gtex/"
setwd(wd)

N.CORES <- 5

# DATA_SET_NAME <- "gtex.sub_tissue.gene_tpm"
DATA_SET_NAME <- "gtex.tissue.gene_tpm"

# ======================================================================= #
# ============================  INITIALIZE  ============================== #
# ======================================================================= #

print("Reading Seurat data object...")

### Load data
#/data/rna-seq/gtex/v7-seurat_objs/gtex.seurat_obj.gene_median_tpm.RData
#/data/rna-seq/gtex/v7-seurat_objs/gtex.seurat_obj.gene_tpm.RData
file.RData.cell_atlas <- "/data/rna-seq/gtex/v7-seurat_objs/gtex.seurat_obj.gene_tpm.RData" # seurat_obj
load(file.RData.cell_atlas) 

# ======================================================================= #
# ========================  CREATE ANNOTATION COL  ====================== #
# ======================================================================= #

colname_ident <- "annotation"

# SMTS: tissue
# SMTSD: sub-tissue
seurat_obj@meta.data$annotation <- seurat_obj@meta.data$SMTS # tissue
# seurat_obj@meta.data$annotation <- seurat_obj@meta.data$SMTSD # sub-tissue


# ======================================================================= #
# ================================  PARAMs  ============================= #
# ======================================================================= #

genes <- rownames(seurat_obj@data) # get all genes

### TEST DATA
# genes <- rownames(seurat_obj@data)[1:10] # get all genes
# seurat_obj <- SubsetData(SetAllIdent(seurat_obj, id=colname_ident), max.cells.per.ident=10) # test data object


### Output files
# file.out.tstat <- sprintf("%s.celltype_expr.ttest.csv.gz", DATA_SET_NAME)
file.out.tstat.human <- sprintf("%s.celltype_expr.ttest.hsapiens_orthologs.csv.gz", DATA_SET_NAME)

# file.out.avg_expr <- sprintf("%s.celltype_expr.avg_expr.csv.gz", DATA_SET_NAME)
file.out.avg_expr.human <- sprintf("%s.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz", DATA_SET_NAME)


# ======================================================================= #
# =======================  Average expression  ========================== #
# ======================================================================= #

print("Running average expression")

### Average
df.avg.expr.human <- anno_specific_expr.avg_expr(colname_ident, genes, seurat_obj) # data frame with columnnames=annotations, rownames=genes 
write_csv(df.avg.expr.human %>% rownames_to_column(var="gene"), path=file.out.avg_expr.human) # write_* automatically deals with .gz, .bz2 or .xz output filenames



# ======================================================================= #
# ================================  t-test  ============================= #
# ======================================================================= #

# print("Running t-tests")
# 
# ### T-stat
# df.celltype_expr_ttest.human <- anno_specific_expr.tstat.all_idents_parallel(colname_ident, genes, seurat_obj, n.cores=N.CORES) # data frame with columnnames=annotations, rownames=genes 
# write_csv(df.celltype_expr_ttest.human %>% rownames_to_column(var="gene"), path=file.out.tstat.human) # write_* automatically deals with .gz, .bz2 or .xz output filenames
# 

# ======================================================================= #
# =======================      FINISH       ========================== #
# ======================================================================= #

print("Script done!")

