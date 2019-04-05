
############### SYNOPSIS ###################
# Create Seurat object from top 10 genetic prioritized cell types
# cell_types_extract = ["MEINH7","HBINH9","DEGLU4","MEINH2","MEGLU6","DEINH3","HBINH3","HBGLU10","HBINH4","HBGLU3"]

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain"
setwd(wd)

library(Seurat)
library(tidyverse)


# ======================================================================= #
# ============================  LOAD DATA  =============================== #
# ======================================================================= #

### Data
file.data <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-L5.top10_cts.csv.gz"
df.data_raw <- read_csv(file.data) %>% rename(gene=X1)
dim(df.data_raw) # ---> 27998   1585
# ^^ duplicates: ???
df.data_raw <- df.data_raw %>% distinct(gene, .keep_all = TRUE) # remove duplicate gene names
df.data_raw <- df.data_raw %>% column_to_rownames(var="gene")
head(df.data_raw)
dim(df.data_raw) # ---> *27933*   1584

### Meta data
file.metadata <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-L5.top10_cts.metadata.csv"
df.metadata <- read_csv(file.metadata)


# ======================================================================= #
# ========================  Create Seurat Obj  ========================== #
# ======================================================================= #

# Seurat: create object and normalize


# seurat_obj <- CreateSeuratObject(df.expression, min.cells = 3, min.genes = 200, is.expr = 0)
seurat_obj <- CreateSeuratObject(as.data.frame(df.data_raw), min.cells = 0, min.genes = 0, is.expr = 0)
# seurat takes as input a data frame or matrix.
# the rownames() are used as gene names
# the colnames()/names() are used as cell id

# seurat_obj@meta.data is 'born with' the meta data columns: "nGene", "nUMI", "orig.ident"
# we could change the name of nUMI to "total_expr"...


seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)

# Add meta data to Seurat
df.tmp.metadata <- data.frame(df.metadata, row.names=df.metadata$CellID)
head(df.tmp.metadata)
seurat_obj <- AddMetaData(seurat_obj, df.tmp.metadata)
str(seurat_obj@meta.data)


# EXPORT: save Seurat object

save(seurat_obj, file="mousebrain-L5.top10_cts.seurat_obj.RData")


