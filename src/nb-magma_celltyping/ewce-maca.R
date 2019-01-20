############### SYNOPSIS ###################
# EWCE generate cell-type data

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


### Runtime full MACA (n=10 cores) ==> ~15 min

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/nb-magma_celltyping"
setwd(wd)

library(Seurat) # only needed if fiddling with Seurat objects
library(MAGMA.Celltyping) # also loads EWCE package

source("magma_celltyping_modified_functions.R")

# ======================================================================= #
# ============================  MACA DATA  =============================== #
# ======================================================================= #

### Load data
file.RData.cell_atlas <- "/data/pub-others/tabula_muris/figshare/180126-facs/maca.seurat_obj.facs.figshare_180126.RData" # 6.1 GB | seurat_obj
load(file.RData.cell_atlas) 
### Set annotation
seurat_obj@meta.data$annotation <- paste(seurat_obj@meta.data$tissue, seurat_obj@meta.data$cell_ontology_class, sep="-")

### Subset/Switch
seurat_obj.use <- seurat_obj
# seurat_obj.use <- SubsetData(SetAllIdent(seurat_obj, id="annotation"), max.cells.per.ident=50, subset.raw=T)
# seurat_obj.use <- SubsetData(SetAllIdent(seurat_obj.use, id="annotation"), ident.use=c("Bladder-mesenchymal cell",
#                                                                                        "Bladder-bladder cell",
#                                                                                        "NA-NA",
#                                                                                        "Bladder-basal cell of urothelium",
#                                                                                        "Brain_Microglia-microglial cell"), subset.raw=T)
# 
# seurat_obj.use


### Setup
expData = as.matrix(seurat_obj.use@raw.data) # matrix, genes x cells [seurat_obj@raw.data is a data.frame]
annotLevels = list(l1=seurat_obj.use@meta.data$annotation)

### Runs in parallel
### https://github.com/NathanSkene/EWCE/blob/master/R/generate.celltype.data.r
filename.ctd = generate.celltype.data.PT(exp=expData,annotLevels=annotLevels,groupName="MACA", no_cores=10) # returns filename for saved celltype_data file
# ^^: saves file "CellTypeData_<groupName>.rda" using 'save(ctd, ..)' command

### Load
# load("CellTypeData_MACA.rda") # loads ctd file
# str(ctd)

