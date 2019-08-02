############### SYNOPSIS ###################
# Export Tabula Muris FACS (v. 180920) Seurat object to CSV

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(Seurat)
library(tidyverse)


wd <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-maca/"
setwd(wd)

# ======================================================================= #
# ============================  INITIALIZE  ============================== #
# ======================================================================= #

print("Reading Seurat data object...")
file.RData.cell_atlas <- "/data/pub-others/tabula_muris/figshare/180920-Robj/tabula_muris.seurat_obj.facs.figshare_180920.RData" # seurat_obj
load(file.RData.cell_atlas) 


# ======================================================================= #
# ================================ EXPORT TO CSV ============================= #
# ======================================================================= #

### Get expression data
df <- as.data.frame(as.matrix(seurat_obj@raw.data)) # need to convert to matrix first. Otherwise you will get the error "cannot coerce class "structure("dgCMatrix", package = "Matrix")" to a data.frame"
dim(df) # 23341 44949
colnames(df)[1:10]; rownames(df)[1:10]
df <- df %>% rownames_to_column(var="gene") %>% select(gene, everything()) %>% as.tibble() # set rownames as column
### df <- GetAssayData(object = seurat_obj, slot = "raw.data") # returns dgCMatrix
file.out.data <- "/scratch/data-for_fast_access/pub-others/tabula_muris_180920/tabula_muris.umi.csv"
data.table::fwrite(df, file=file.out.data,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch 
R.utils::gzip(file.out.data, overwrite=TRUE) # gzip

### Write meta-data
df.metadata <- seurat_obj@meta.data %>% as.tibble() %>% select(cell_id=cell,nGene,nReads,tissue,subtissue_clean,celltype,tissue_celltype,tissue_subtissue_celltype)
file.out.meta <- "/scratch/data-for_fast_access/pub-others/tabula_muris_180920/tabula_muris.metadata.csv"
data.table::fwrite(df.metadata, file=file.out.meta,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch 


# ======================================================================= #
# ================================ XXXXXXX ============================= #
# ======================================================================= #


# ======================================================================= #
# ================================ XXXXXXXX ============================= #
# ======================================================================= #

