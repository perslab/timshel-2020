############### SYNOPSIS ###################
# Export adipocyte Seurat object to CSV

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
library(here)

here("src/GE-adipocyte")

# ======================================================================= #
# ============================  INITIALIZE  ============================== #
# ======================================================================= #

print("Reading Seurat data object...")
# Here is the path to the Seurat object with gene ID's instead of gene symbols: 
file.RData.cell_atlas <- "/projects/pytrik/sc_adipose/analyze_10x_fluidigm/10x-adipocyte-analysis/output/10x-180831_gene-ids" # 
seurat_obj <- readRDS(file.RData.cell_atlas) 

# ======================================================================= #
# ============================ CLEAN OUT CELLS IN @raw.data SLOT  ============================== #
# ======================================================================= #

### @raw.data slot contains *MORE CELLS* than all the other slots. 
### Remove cells from @raw.data. Keep only cells in the metadata (or equivalently 'seurat_obj@cell.names')
seurat_obj <- SubsetData(seurat_obj, cells.use=row.names(seurat_obj@meta.data), subset.raw=T)

### Check that cell order are the same (just to make sure)
all(seurat_obj@cell.names==seurat_obj@raw.data@Dimnames[[2]]) # ok

# ======================================================================= #
# ================================ EXPORT TO CSV ============================= #
# ======================================================================= #
prefix_out <- "preadipocyte_developing_1808"
dir.out <- "/projects/timshel/sc-scheele_lab_adipose_fluidigm_c1/data-preadipocytes_developing"

### Get expression data
df <- as.data.frame(as.matrix(seurat_obj@raw.data)) # need to convert to matrix first. Otherwise you will get the error "cannot coerce class "structure("dgCMatrix", package = "Matrix")" to a data.frame"
dim(df) # 22979 88715
colnames(df)[1:10]; rownames(df)[1:10]
df <- df %>% rownames_to_column(var="gene") %>% select(gene, everything()) %>% as.tibble() # set rownames as column
### df <- GetAssayData(object = seurat_obj, slot = "raw.data") # returns dgCMatrix
file.out.data <- file.path(dir.out, sprintf("%s.umi.csv", prefix_out))
file.out.data
data.table::fwrite(df, file=file.out.data,  # fwrite cannot write gziped files
                  nThread=24, verbose=T) # write file ---> write to scratch
R.utils::gzip(file.out.data, overwrite=TRUE) # gzip

### Write meta-data
df.metadata <- seurat_obj@meta.data %>% 
  as.tibble() %>%
  select(cell_id=cellID, 
         branch_low_res, 
         branch_high_res, 
         everything()) # rearrange columns
head(df.metadata)
file.out.meta <- file.path(dir.out, sprintf("%s.metadata.csv", prefix_out))
data.table::fwrite(df.metadata, file=file.out.meta,  # fwrite cannot write gziped files
                  nThread=24, verbose=T) # write file ---> write to scratch


# ======================================================================= #
# ================================ XXXXXXX ============================= #
# ======================================================================= #


# ======================================================================= #
# ================================ XXXXXXXX ============================= #
# ======================================================================= #

