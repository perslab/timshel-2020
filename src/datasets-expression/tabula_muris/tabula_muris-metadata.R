############### SYNOPSIS ###################
# Cell-type metadata

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
# ============================  PARAMS  ============================== #
# ======================================================================= #

wd <- here("src/datasets-expression/tabula_muris/")
setwd(wd)

# ======================================================================= #
# ============================  INITIALIZE  ============================== #
# ======================================================================= #

### Load data
# file.RData.cell_atlas <- "/data/pub-others/tabula_muris/figshare/180126-facs/maca.seurat_obj.facs.figshare_180126.RData" # 6.1 GB | seurat_obj
file.RData.cell_atlas <- "/data/pub-others/tabula_muris/figshare/180920-Robj/tabula_muris.seurat_obj.facs.figshare_180920.RData"
load(file.RData.cell_atlas) 

# ======================================================================= #
# =========================  MAKE METADATA FILE  ========================== #
# ======================================================================= #


df.metadata <- seurat_obj@meta.data %>% group_by(tissue_celltype) %>% summarise(n_cells = n()) # count
df.metadata <- df.metadata %>% separate(tissue_celltype, into=c("tissue", "cell_type"), sep="\\.", remove=F, extra="merge") # split

file.out <- "tabula_muris_facs.tissue_celltype.celltype_metadata.csv"
df.metadata %>% write_csv(path=file.out)


# ======================================================================= #
# ============================  XXXX  ============================== #
# ======================================================================= #

# df.annotation <- read_csv("/data/pub-others/tabula_muris/figshare/180126-facs/metadata_FACS.csv")
# df.annotation <- read_csv("/data/pub-others/tabula_muris/figshare/180126-facs/annotations_FACS.csv")
# df.summary <- df.annotation %>% distinct(tissue, cell_ontology_class)







