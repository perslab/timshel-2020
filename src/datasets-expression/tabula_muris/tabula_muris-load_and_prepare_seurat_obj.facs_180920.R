############### SYNOPSIS ###################
# Prepare tabula muris figshare/180920-Robj/FACS_all.Robj data object: add annotations and save seurat object

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

# ======================================================================= #
# ============================  XXXXX  ============================== #
# ======================================================================= #

### 

### Load
file.RData.cell_atlas <- "/data/pub-others/tabula_muris/figshare/180920-Robj/FACS_all.Robj" # 8.1 GB | # tiss_FACS
load(file.RData.cell_atlas) 
seurat_obj <- tiss_FACS; 
rm(tiss_FACS) # rename

### Checks
any(is.na(seurat_obj@meta.data$tissue)) # ---> NO NA values in 'tissue'
# unique(seurat_obj@meta.data$subtissue)
sum(is.na(seurat_obj@meta.data$cell_ontology_class)) # ---> 170 cells with NA in cell_ontology_class
df.cell_ontology_class_with_na <- seurat_obj@meta.data %>% group_by(tissue) %>% summarise(n = sum(is.na(cell_ontology_class))) %>% arrange(desc(n))
df.cell_ontology_class_with_na
# tissue             n
# <chr>             <int>
#   1 Fat                 102
# 2 Lung                 40
# 3 Heart                28
# 4 Bladder               0
# 5 Brain_Myeloid         0
# 6 Brain_Non-Myeloid     0
df.cell_ontology_class_counts <- seurat_obj@meta.data %>% group_by(cell_ontology_class) %>% summarise(n = n()) %>% arrange(desc(n))
df.cell_ontology_class_counts # ---> there are no cells marked with 'unknown' or the like

# ======================================================================= #
# ============================  Set annotations  ============================== #
# ======================================================================= #

###  We use dot "." as seperated, because '-' and '_' are already used in some of the annotations.

### CELLTYPE
seurat_obj@meta.data$celltype <- with(seurat_obj@meta.data, cell_ontology_class)
### Replace "NA" in celltype with "unknown cell type"
seurat_obj@meta.data$celltype[is.na(seurat_obj@meta.data$celltype)] <- "unknown_cell_type"

### TISSUE CELLTYPE
seurat_obj@meta.data$tissue_celltype <- with(seurat_obj@meta.data, paste0(tissue, ".", celltype))

#### TISSUE SUBTISSUE CELLTYPE
seurat_obj@meta.data$subtissue_clean <- seurat_obj@meta.data$subtissue
seurat_obj@meta.data$subtissue_clean[seurat_obj@meta.data$subtissue_clean=="Unknown"] <- NA
seurat_obj@meta.data$subtissue_clean[seurat_obj@meta.data$subtissue_clean==""] <- NA
seurat_obj@meta.data$tissue_subtissue_celltype <- with(seurat_obj@meta.data, paste0(tissue, ".", subtissue_clean, ".", celltype)) # ---> this annotation is not very meaningful


### REPLACE ANY SPACES with UNDERSCORE (_). str_replace_all(): Replace all matched patterns in each string.
seurat_obj@meta.data$celltype <- stringr::str_replace_all(seurat_obj@meta.data$celltype, pattern="\\s+", replacement="_")
seurat_obj@meta.data$tissue_celltype <- stringr::str_replace_all(seurat_obj@meta.data$tissue_celltype, pattern="\\s+", replacement="_")
seurat_obj@meta.data$subtissue_clean <- stringr::str_replace_all(seurat_obj@meta.data$subtissue_clean, pattern="\\s+", replacement="_")
seurat_obj@meta.data$tissue_subtissue_celltype <- stringr::str_replace_all(seurat_obj@meta.data$tissue_subtissue_celltype, pattern="\\s+", replacement="_")

### EXTRA
df.n_tissue_celltype <- seurat_obj@meta.data %>% group_by(tissue_celltype) %>% summarise(n = n()) %>% arrange(n)
df.n_tissue_celltype
# tissue_celltype                                          n
# <chr>                                                <int>
#   1 Skin.leukocyte                                          15
# 2 Kidney.leukocyte                                        16
# 3 Marrow.pre-natural killer cell                          22
# 4 Lung.ciliated columnar cell of tracheobronchial tree    25
# 5 Marrow.basophil                                         25
# 6 Marrow.regulatory T cell                                27

# ======================================================================= #
# ============================ CLEAN OUT CELLS IN @raw.data SLOT  ============================== #
# ======================================================================= #

### @raw.data slot contains *MORE CELLS* than all the other slots. The is likely because Biohub did not use SubsetData(..., subset.raw=True)
dim(seurat_obj@data) # 23341 44949 | --> @data contains 44949 cells
dim(seurat_obj@meta.data) # 44949    23 | ---> @meta.data contains 44949 cells 
dim(seurat_obj@raw.data) # 23341 53760 | ---> @raw.data contains 53760 cells
# all(rownames(seurat_obj@meta.data)==seurat_obj@meta.data$cell) # True --> cell column and metadata rownames are identical. We just check to be sure.

### Remove cells from @raw.data. Keep only cells in the metadata
seurat_obj <- SubsetData(seurat_obj, cells.use=seurat_obj@meta.data$cell, subset.raw=T)

### Check that cell order are the same (just to make sure)
all(seurat_obj@cell.names==seurat_obj@raw.data@Dimnames[[2]]) # ok
all(seurat_obj@cell.names==seurat_obj@meta.data$cell) # ok

# ======================================================================= #
# ============================ EXPORT  ============================== #
# ======================================================================= #
           
### Normalize data (just to be sure)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)

### Save Robj
save(seurat_obj, file="/data/pub-others/tabula_muris/figshare/180920-Robj/tabula_muris.seurat_obj.facs.figshare_180920.RData")
