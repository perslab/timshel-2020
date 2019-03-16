############### SYNOPSIS ###################
# Get cell_type annotations

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)

# ======================================================================= #
# ============================  PARAMS  ============================== #
# ======================================================================= #

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-maca/"
setwd(wd)

# ======================================================================= #
# ============================  INITIALIZE  ============================== #
# ======================================================================= #

### Load data
file.RData.cell_atlas <- "/data/pub-others/tabula_muris/figshare/180126-facs/maca.seurat_obj.facs.figshare_180126.RData" # 6.1 GB | seurat_obj
load(file.RData.cell_atlas) 

# ======================================================================= #
# =========================  MAKE METADATA FILE  ========================== #
# ======================================================================= #


df.metadata <- seurat_obj@meta.data %>% group_by(tissue_celltype) %>% summarise(n_cells = n()) # count
df.metadata <- df.metadata %>% separate(tissue_celltype, into=c("tissue", "cell_type"), sep="\\.", remove=F, extra="merge")# split
df.metadata %>% write_csv(path="tabula_muris_facs.tissue_celltype.metadata.csv")
# /raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/tabula_muris/

# ======================================================================= #
# ============================  XXXX  ============================== #
# ======================================================================= #

# df.annotation <- read_csv("/data/pub-others/tabula_muris/figshare/180126-facs/metadata_FACS.csv")
# df.annotation <- read_csv("/data/pub-others/tabula_muris/figshare/180126-facs/annotations_FACS.csv")
# df.summary <- df.annotation %>% distinct(tissue, cell_ontology_class)







