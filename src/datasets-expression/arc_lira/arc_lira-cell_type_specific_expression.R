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

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-arc_lira/"
setwd(wd)

N.CORES <- 10
DATA_SET_NAME <- "arc_lira"

# ======================================================================= #
# ============================  INITIALIZE  ============================== #
# ======================================================================= #

### Load data
file.RData.cell_atlas <- "/projects/timshel/sc-arc_lira/src/out-generate_cell_atlas/arc_lira_cell_atlas-core_objs.RData"
load(file.RData.cell_atlas) # seurat_obj, df.cluster_markers

# ======================================================================= #
# ============================  p.res 0.8  ============================== #
# ======================================================================= #

### cell type annotation
p.res <- 0.8; p.res.string <- paste0("res.", p.res) # RES 0.8
colname_ident <- "annotation"

## Assigning cell type identity to clusters
df.cluster_annotation <- get_cell_atlas_annotation(resolution=p.res.string) # CALLS FUNCTIONS
seurat_obj <- SetAllIdent(seurat_obj, id=p.res.string)
cell_type_mapped <- plyr::mapvalues(x=seurat_obj@ident, from=df.cluster_annotation$cluster, to=df.cluster_annotation$annotation)
seurat_obj@ident <- cell_type_mapped # Set cell types ident
seurat_obj <- StashIdent(seurat_obj, save.name = "annotation") # saves ident in @meta.data


# ======================================================================= #
# ============================  p.res 4  ============================== #
# ======================================================================= #

# p.res <- 4; p.res.string <- paste0("res.", p.res) # RES 4
# colname_ident <- p.res.string

# ======================================================================= #
# ================================  PARAMs  ============================= #
# ======================================================================= #
genes <- rownames(seurat_obj@data) # get all genes
# genes <- rownames(seurat_obj@data)[1:10] # get all genes

### Output files
file.out.tstat <- sprintf("%s.celltype_expr.ttest.csv.gz", DATA_SET_NAME)
file.out.tstat.human <- sprintf("%s.celltype_expr.ttest.hsapiens_orthologs.csv.gz", DATA_SET_NAME)

file.out.avg_expr <- sprintf("%s.celltype_expr.avg_expr.csv.gz", DATA_SET_NAME)
file.out.avg_expr.human <- sprintf("%s.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz", DATA_SET_NAME)



# ======================================================================= #
# ================================  t-test  ============================= #
# ======================================================================= #

### T-stat
df.celltype_expr_ttest <- anno_specific_expr.tstat.all_idents_parallel(colname_ident, genes, seurat_obj, n.cores=N.CORES) # data frame with columnnames=annotations, rownames=genes 
write_csv(df.celltype_expr_ttest %>% rownames_to_column(var="gene"), path=file.out.tstat) # write_* automatically deals with .gz, .bz2 or .xz output filenames


### Ortholog
df.celltype_expr_ttest.human <- mouse_to_human_ortholog_gene_expression_mapping(df.celltype_expr_ttest)
write_csv(df.celltype_expr_ttest.human %>% rownames_to_column(var="gene"), path=file.out.tstat.human) # write_* automatically deals with .gz, .bz2 or .xz output filenames


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
# =======================      FINISH       ========================== #
# ======================================================================= #

print("Script done.")

