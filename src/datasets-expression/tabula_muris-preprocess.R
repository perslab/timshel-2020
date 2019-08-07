############### SYNOPSIS ###################
# Process Tabula Muris FACS (v. 180920, FACS_all.Robj) Seurat data object: add annotations and save seurat object
# Export Seurat object to CSV

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
# ============================  DOWNLOAD  ============================== #
# ======================================================================= #

# WEBSTE URL: https://figshare.com/articles/Robject_files_for_tissues_processed_by_Seurat/5821263/2
# Version 2,  Dataset posted on 21.09.2018, 00:01 by Tabula Muris Consortium
# DONWNLOAD DATE: Nov 5th 2018

file_download <- here("tmp-data/expression/tabula_muris-FACS_all.RData")
if (!file.exists(file_download)) {
  download_url <- "https://ndownloader.figshare.com/files/13092806"
  download.file(download_url, destfile=file_download)
}


# ======================================================================= #
# ============================  LOAD  ============================== #
# ======================================================================= #

### Load
load(file_download) 
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
# ==================  EXPORT CELL-TYPE METADATA FILE  =================== #
# ======================================================================= #


df.metadata <- seurat_obj@meta.data %>% group_by(tissue_celltype) %>% summarise(n_cells = n()) # count
df.metadata <- df.metadata %>% separate(tissue_celltype, into=c("tissue", "cell_type"), sep="\\.", remove=F, extra="merge") # split

file.out <- here("data/expression/tabula_muris/tabula_muris_facs.tissue_celltype.celltype_metadata.csv")
df.metadata %>% write_csv(path=file.out)


# ======================================================================= #
# ================================ EXPORT TO CSV ============================= #
# ======================================================================= #

### Get expression data
df <- as.data.frame(as.matrix(seurat_obj@raw.data)) # need to convert to matrix first. Otherwise you will get the error "cannot coerce class "structure("dgCMatrix", package = "Matrix")" to a data.frame"
dim(df) # 23341 44949
colnames(df)[1:10]; rownames(df)[1:10]
df <- df %>% rownames_to_column(var="gene") %>% select(gene, everything()) %>% as.tibble() # set rownames as column
### df <- GetAssayData(object = seurat_obj, slot = "raw.data") # returns dgCMatrix
file.out.data <- here("tmp-data/expression/tabula_muris.umi.csv"
data.table::fwrite(df, file=file.out.data,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch 
R.utils::gzip(file.out.data, overwrite=TRUE) # gzip

### Write meta-data
df.metadata <- seurat_obj@meta.data %>% as.tibble() %>% select(cell_id=cell,nGene,nReads,tissue,subtissue_clean,celltype,tissue_celltype,tissue_subtissue_celltype)
file.out.meta <- "tmp-data/expression/tabula_muris.metadata.csv"
data.table::fwrite(df.metadata, file=file.out.meta,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch 

