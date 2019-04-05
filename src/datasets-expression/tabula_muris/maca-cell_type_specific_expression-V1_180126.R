############### SYNOPSIS ###################
# Compute cell-type specific expression and ortholog mapping (map mouse genes to human)

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ==============================  USAGE  =============================== #
# ======================================================================= #

# time Rscript maca-cell_type_specific_expression.R |& tee maca-cell_type_specific_expression.out.txt

# time Rscript maca-cell_type_specific_expression.R maca.per_tissue_celltype |& tee maca-cell_type_specific_expression.per_tissue_celltype.out.txt
# time Rscript maca-cell_type_specific_expression.R maca.per_tissue |& tee maca-cell_type_specific_expression.per_tissue.out.txt
# time Rscript maca-cell_type_specific_expression.R maca.per_celltype |& tee maca-cell_type_specific_expression.per_celltype.out.txt


# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(here)

library(Seurat)
library(tidyverse)
library(Matrix)

library(parallel)


### Source custom scripts
# dir.project_src <- "/projects/timshel/sc-arc_lira/src"
# dir.pers_lab_sc_lib <- "/projects/timshel/git/perslab-sc-library"
# source(sprintf("%s/constants-cell_type_annotations.R", dir.project_src)) # loads cell type annotations
# source(sprintf("%s/seurat_functions/load_functions.R", dir.pers_lab_sc_lib)) # load Pers lab/Timshel single-cell library

source(here("src/lib/load_functions.R")) # load sc-genetics library


# ======================================================================= #
# ============================  PARAMS  ============================== #
# ======================================================================= #

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-maca/"
setwd(wd)

N.CORES <- 20

DATA_SET_NAME <- commandArgs(trailingOnly=TRUE)[1]
print(sprintf("RUNNING DATA_SET_NAME=%s", DATA_SET_NAME))

# DATA_SET_NAME <- "maca.per_tissue_celltype"
# DATA_SET_NAME <- "maca.per_tissue"
# DATA_SET_NAME <- "maca.per_celltype"

# ======================================================================= #
# ============================  INITIALIZE  ============================== #
# ======================================================================= #

print("Reading Seurat data object...")

### Load data
# file.RData.cell_atlas <- "/projects/timshel/maca/All_seurat_tiss.Robj" # Plate : 14 GB | tiss variable
# file.RData.cell_atlas <- "/projects/timshel/maca/10x_All_seurat_tiss.Robj" # 10x : 7.6G
file.RData.cell_atlas <- "/data/pub-others/tabula_muris/figshare/180126-facs/maca.seurat_obj.facs.figshare_180126.RData" # 6.1 GB | seurat_obj
load(file.RData.cell_atlas) 

# ======================================================================= #
# ========================  CREATE ANNOTATION COL  ====================== #
# ======================================================================= #

colname_ident <- "annotation"

if (DATA_SET_NAME == "maca.per_tissue_celltype") {
  seurat_obj@meta.data$annotation <- paste(seurat_obj@meta.data$tissue, seurat_obj@meta.data$cell_ontology_class, sep="-")
} else if (DATA_SET_NAME == "maca.per_tissue") {
  seurat_obj@meta.data$annotation <- seurat_obj@meta.data$tissue
} else if (DATA_SET_NAME == "maca.per_celltype") {
  seurat_obj@meta.data$annotation <- seurat_obj@meta.data$cell_ontology_class
} else {
  warning("Got the wrong argument")
}

# ======================================================================= #
# ================================  PARAMs  ============================= #
# ======================================================================= #

genes <- rownames(seurat_obj@data) # get all genes

### TEST DATA
# genes <- rownames(seurat_obj@data)[1:10] # get all genes
# seurat_obj <- SubsetData(SetAllIdent(seurat_obj, id=colname_ident), max.cells.per.ident=10) # test data object


### Output files
file.out.tstat <- sprintf("%s.celltype_expr.ttest.csv.gz", DATA_SET_NAME)
file.out.tstat.human <- sprintf("%s.celltype_expr.ttest.hsapiens_orthologs.csv.gz", DATA_SET_NAME)

file.out.avg_expr <- sprintf("%s.celltype_expr.avg_expr.csv.gz", DATA_SET_NAME)
file.out.avg_expr.human <- sprintf("%s.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz", DATA_SET_NAME)


# ======================================================================= #
# =======================  Average expression  ========================== #
# ======================================================================= #

print("Running average expression")

### Average
df.avg.expr <- anno_specific_expr.avg_expr(colname_ident, genes, seurat_obj) # data frame with columnnames=annotations, rownames=genes 
write_csv(df.avg.expr %>% rownames_to_column(var="gene"), path=file.out.avg_expr) # write_* automatically deals with .gz, .bz2 or .xz output filenames


### Ortholog
df.avg.expr.human <- mouse_to_human_ortholog_gene_expression_mapping(df.avg.expr)
write_csv(df.avg.expr.human %>% rownames_to_column(var="gene"), path=file.out.avg_expr.human) # write_* automatically deals with .gz, .bz2 or .xz output filenames


# ======================================================================= #
# ================================  t-test  ============================= #
# ======================================================================= #

print("Running t-tests")

### T-stat
df.celltype_expr_ttest <- anno_specific_expr.tstat.all_idents_parallel(colname_ident, genes, seurat_obj, n.cores=N.CORES) # data frame with columnnames=annotations, rownames=genes 
write_csv(df.celltype_expr_ttest %>% rownames_to_column(var="gene"), path=file.out.tstat) # write_* automatically deals with .gz, .bz2 or .xz output filenames


### Ortholog
df.celltype_expr_ttest.human <- mouse_to_human_ortholog_gene_expression_mapping(df.celltype_expr_ttest)
write_csv(df.celltype_expr_ttest.human %>% rownames_to_column(var="gene"), path=file.out.tstat.human) # write_* automatically deals with .gz, .bz2 or .xz output filenames


# ======================================================================= #
# =======================      FINISH       ========================== #
# ======================================================================= #

print("Script done!")

