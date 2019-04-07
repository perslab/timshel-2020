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
library(parallel)

source(here("src/lib/load_functions.R")) # load sc-genetics library


# ======================================================================= #
# ============================  PARAMS  ============================== #
# ======================================================================= #

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-hypothalamus/"
setwd(wd)

N.CORES <- 10
DATA_SET_NAME <- "hypothalamus"

# ======================================================================= #
# ============================  INITIALIZE  ============================== #
# ======================================================================= #

### Load data
file.RData.cell_atlas <- "/projects/mludwig/Dataset_Alignment/Four_datasets/All_cell_types/Hypo_aligned_res_06.RData"
load(file.RData.cell_atlas) # hypo
seurat_obj <- hypo; rm(hypo) # copy and clean

### Normalize: we make sure that we know how the data has been normalized.
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)

### test
colSums(seurat_obj@data[,1:10])
colSums(exp(seurat_obj@data[,1:10])-1) # --> gives the same "scale.factor" value for every column.


# ======================================================================= #
# ============================  PREP  ============================== #
# ======================================================================= #

### cell type annotation
colname_ident <- "annotation"

### Creating annotation
seurat_obj@meta.data %>% distinct(protocol)
seurat_obj@meta.data %>% distinct(cell_type)
seurat_obj@meta.data$annotation <- paste0(seurat_obj@meta.data$protocol, "_", seurat_obj@meta.data$cell_type)
seurat_obj@meta.data %>% distinct(annotation)


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

# ### T-stat
# df.celltype_expr_ttest <- anno_specific_expr.tstat.all_idents_parallel(colname_ident, genes, seurat_obj, n.cores=N.CORES) # data frame with columnnames=annotations, rownames=genes 
# write_csv(df.celltype_expr_ttest %>% rownames_to_column(var="gene"), path=file.out.tstat) # write_* automatically deals with .gz, .bz2 or .xz output filenames
# 
# 
# ### Ortholog
# df.celltype_expr_ttest.human <- mouse_to_human_ortholog_gene_expression_mapping(df.celltype_expr_ttest)
# write_csv(df.celltype_expr_ttest.human %>% rownames_to_column(var="gene"), path=file.out.tstat.human) # write_* automatically deals with .gz, .bz2 or .xz output filenames


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
# =======================  Per dataset analysis (avg)  ========================== #
# ======================================================================= #

genes <- rownames(seurat_obj@data) # get all genes

### cell type annotation
colname_ident <- "annotation"

### Creating annotation
seurat_obj@meta.data %>% distinct(protocol)
seurat_obj@meta.data %>% distinct(cell_type)
seurat_obj@meta.data$annotation <- paste0(seurat_obj@meta.data$protocol, "_", seurat_obj@meta.data$cell_type)
seurat_obj@meta.data %>% distinct(annotation)


### Creating annotation
datasets <- seurat_obj@meta.data %>% distinct(protocol) %>% pull(protocol) # "chen", "campbell", "lira", "romanov"

for (dataset in datasets) {
  
  print(sprintf("dataset = %s", dataset))
  
  ### Output files
  DATA_SET_NAME <- sprintf("hypothalamus_%s", dataset) # setting file output name
  file.out.tstat <- sprintf("%s.celltype_expr.ttest.csv.gz", DATA_SET_NAME)
  file.out.tstat.human <- sprintf("%s.celltype_expr.ttest.hsapiens_orthologs.csv.gz", DATA_SET_NAME)
  file.out.avg_expr <- sprintf("%s.celltype_expr.avg_expr.csv.gz", DATA_SET_NAME)
  file.out.avg_expr.human <- sprintf("%s.celltype_expr.avg_expr.hsapiens_orthologs.csv.gz", DATA_SET_NAME)
  
  ### Get Seurat object
  print("Subsetting data...")
  seurat_obj.sub <- SubsetData(seurat_obj, subset.name="protocol", accept.value=dataset) #COULD TRY THIS: do.clean=T, subset.raw=T
  print("Done subsetting")
  
  ### Average
  df.avg.expr <- anno_specific_expr.avg_expr(colname_ident, genes, seurat_obj.sub) # data frame with columnnames=annotations, rownames=genes 
  # ^ OBS: input is *seurat_obj.sub*
  write_csv(df.avg.expr %>% rownames_to_column(var="gene"), path=file.out.avg_expr) # write_* automatically deals with .gz, .bz2 or .xz output filenames
  
  
  ### Ortholog
  df.avg.expr.human <- mouse_to_human_ortholog_gene_expression_mapping(df.avg.expr)
  write_csv(df.avg.expr.human %>% rownames_to_column(var="gene"), path=file.out.avg_expr.human) # write_* automatically deals with .gz, .bz2 or .xz output filenames
}

# ======================================================================= #
# =======================      FINISH       ========================== #
# ======================================================================= #

print("Script done.")

