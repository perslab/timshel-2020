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

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


# ======================================================================= #
# ============================  PARAMS  ============================== #
# ======================================================================= #

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-novo_bulk/"
setwd(wd)

N.CORES <- 10
DATA_SET_NAME <- "novo_bulk"

# ======================================================================= #
# ============================  LOAD DATA AND CREATE SEURAT OBJECT  ============================== #
# ======================================================================= #

# file.data <- "/data/rna-seq/nn-lira-sema/171208/GUS2016-142-NN_Tables-ARC_counts_only.csv"
# file.data <- "/data/rna-seq/nn-lira-sema/171208/GUS2016-142-NN_Tables-ARC_RPKM_only.csv"
file.data <- "/data/rna-seq/nn-lira-sema/171208/GUS2016-142-NN_Tables-counts_only.csv" # contains all data in counts values.

df.data.raw <- read_csv(file.data)
dim(df.data.raw) # genes are rows; samples are columns

### Let's tidy up data: 
# 1. remove duplicate gene symbols (they will cause trouble later) by *only keeping one of the entries*. (consider using Ensembl IDs instead to avoid duplicates)
# 2. set gene symbols as rownames of the data frame to be able to load smoothly into Seurat.

### Duplicated gene names
df.data.raw[duplicated(df.data.raw[,"Gene name"]), "Gene name"]

### Keep distinct
df.data <- df.data.raw %>% distinct(`Gene name`, .keep_all=T) 
# distinct: 
# If there are multiple rows for a given combination of inputs, ***only the first row will be preserved***.
# .keep_all: If TRUE, keep all variables in .data. If a combination of ... is not distinct, this keeps the first row of values.

# run once!
df.data <- df.data %>% column_to_rownames(var="Gene name") %>% select(-`Ensembl acc key`, -Description)

# IMPORTANT: Seurat does not accept a tibble as input. It will fail when doing e.g. NormalizeData
seurat_obj <- CreateSeuratObject(as.matrix(df.data), min.cells=0, min.genes=0, is.expr = 0)
seurat_obj

### Add metadata
brain_area <- as.character(sapply(strsplit(seurat_obj@cell.names, split = "_"), '[[', 1)) # e.g. "Arc"
tmp.string <- as.character(sapply(strsplit(seurat_obj@cell.names, split = "_"), '[[', 2)) # e.g. "Vehicle animal 20907"
treatment <- as.character(sapply(strsplit(tmp.string, split = " "), '[[', 1)) # e.g. "Vehicle"
animal <- as.character(sapply(strsplit(tmp.string, split = " "), '[[', 3)) # e.g. "20907"
df.metadata <- data.frame(brain_area_treatment=paste0(brain_area, "_", treatment),
                          brain_area=brain_area,
                          treatment=treatment,
                          animal=animal,
                          dummy_all_samples="dummy",
                          row.names=seurat_obj@cell.names)

seurat_obj <- AddMetaData(seurat_obj, df.metadata)

### Normalize data
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor=1e4)

##### Save
# save(seurat_obj, file = "nn-lira-sema.seurat_obj.RData") 
