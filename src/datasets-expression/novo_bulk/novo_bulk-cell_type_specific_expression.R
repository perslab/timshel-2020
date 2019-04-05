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

# time Rscript arc_lira-cell_type_specific_expression.R |& tee arc_lira-cell_type_specific_expression.out.txt

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(Seurat)
library(tidyverse)

library(parallel)

source(here("src/lib/load_functions.R")) # load sc-genetics library


# ======================================================================= #
# ============================  PARAMS  ============================== #
# ======================================================================= #

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-novo_bulk/"
setwd(wd)

N.CORES <- 10
DATA_SET_NAME <- "novo_bulk"

# =======================================

# ======================================================================= #
# ============================  LOAD DATA  ============================== #
# ======================================================================= #


load("/data/rna-seq/nn-lira-sema/171208/nn-lira-sema.seurat_obj.RData") # EDIT 20 Jan 2019

# ======================================================================= #
# ============================  PREP  ============================== #
# ======================================================================= #

### cell type annotation
colname_ident <- "annotation"

### Creating annotation
seurat_obj@meta.data$annotation <- paste0(seurat_obj@meta.data$brain_area, "_", seurat_obj@meta.data$treatment)
seurat_obj@meta.data %>% distinct(annotation)

# ======================================================================= #
# ================================  PARAMs  ============================= #
# ======================================================================= #

genes <- rownames(seurat_obj@data) # get all genes
# genes <- rownames(seurat_obj@data)[1:10] # get all genes

### Output files
file.out.avg_expr <- sprintf("%s.celltype_expr.avg_expr.csv.gz", DATA_SET_NAME)
file.out.avg_expr.human <- sprintf("%s.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz", DATA_SET_NAME)

# ======================================================================= #
# =======================  Average expression  ========================== #
# ======================================================================= #

### Average
df.avg.expr <- anno_specific_expr.avg_expr(colname_ident, genes, seurat_obj) # data frame with columnnames=annotations, rownames=genes 
write_csv(df.avg.expr %>% rownames_to_column(var="gene"), path=file.out.avg_expr) # write_* automatically deals with .gz, .bz2 or .xz output filenames


### Ortholog
df.avg.expr.human <- mouse_to_human_ortholog_gene_expression_mapping(df.avg.expr)
write_csv(df.avg.expr.human %>% rownames_to_column(var="gene"), path=file.out.avg_expr.human) # write_* automatically deals with .gz, .bz2 or .xz output filenames

# ======================================================================= #
# ======================= SEM models ========================== #
# ======================================================================= #

source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library

calculate_and_write_sems(df.avg_expr=df.avg.expr.human, name.dataset=DATA_SET_NAME) # run on human mapped genes
