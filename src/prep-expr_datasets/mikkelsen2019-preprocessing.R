############### SYNOPSIS ###################
# Download and prepare mikkelsen-2019-natureneuroscience hypothalamus gene expression data

### OUTPUT:
# ....

### REMARKS:
# ....

### REFERENCE:

# Paper
# https://www.nature.com/articles/s41593-019-0349-8

# Supplementary materials
# 

# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library("tidyverse")
library("openxlsx")
library("here")
library("Seurat")
library("data.table")

# ======================================================================= #
# ========================== Download ================================ #
# ======================================================================= #
### You only need to run this step once

if (!dir.exists(here("tmp-data"))) dir.create(here("tmp-data"))
if (!dir.exists(here("tmp-data","expression"))) dir.create(here("tmp-data", "expression"))
if (!dir.exists(here("data","expression"))) dir.create(here("data", "expression"))
if (!dir.exists(here("data","expression", "mikkelsen2019"))) dir.create(here("data", "expression", "mikkelsen2019"))

# raw data
rawDataDownload <- here("tmp-data", "expression","mikkelsen2019-GSE125065_RAW.tar")

if (!file.exists(rawDataDownload)) {
  # Download UMI data
  downloadURLraw <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125065/suppl/GSE125065_RAW.tar"
  download.file(downloadURLraw, destfile=rawDataDownload)
}

# unbundle tarball contents:

# Archive/File	Name	Time	Size	Type
# Archive	GSE125065_RAW.tar	01/15/2019 12:55:11	61061120	TAR
# File	GSM3562050_AJ17001_barcodes.tsv.gz	01/14/2019 14:53:56	14793	TSV
# File	GSM3562050_AJ17001_genes.tsv.gz	01/14/2019 14:53:57	222705	TSV
# File	GSM3562050_AJ17001_matrix.mtx.gz	01/14/2019 14:54:03	27825957	MTX
# File	GSM3562051_AJ17002_barcodes.tsv.gz	01/14/2019 14:54:03	16075	TSV
# File	GSM3562051_AJ17002_genes.tsv.gz	01/14/2019 14:54:03	222705	TSV
# File	GSM3562051_AJ17002_matrix.mtx.gz	01/14/2019 14:54:10	32744552	MTX

# Need to set wd as tar will extract files to here
setwd(here("tmp-data", "expression"))

system2(command="tar", args= c("-xvf", rawDataDownload))

# sample GSM3562050

mat_GSM3562050 <- Matrix::readMM(here("tmp-data", "expression","GSM3562050_AJ17001_matrix.mtx.gz"))
vec_genes_GSM3562050 <- read.table(file=here("tmp-data", "expression","GSM3562050_AJ17001_genes.tsv.gz"))[[1]] %>% as.character 
vec_barcode_GSM3562050 <- read.table(file=here("tmp-data", "expression","GSM3562050_AJ17001_barcodes.tsv.gz"))[[1]] %>% as.character 

rownames(mat_GSM3562050) <- vec_genes_GSM3562050
colnames(mat_GSM3562050) <- vec_barcode_GSM3562050

dim(mat_GSM3562050)
# [1] 28692  3439

spmat_GSM3562050 <- as(object = mat_GSM3562050, Class="sparseMatrix")
rm(mat_GSM3562050)

# sample GSM3562051

mat_GSM3562051 <- Matrix::readMM(here("tmp-data", "expression","GSM3562051_AJ17002_matrix.mtx.gz"))
vec_genes_GSM3562051 <- read.table(file=here("tmp-data", "expression","GSM3562051_AJ17002_genes.tsv.gz"))[[1]] %>% as.character 
vec_barcode_GSM3562051 <- read.table(file=here("tmp-data", "expression","GSM3562051_AJ17002_barcodes.tsv.gz"))[[1]] %>% as.character 

rownames(mat_GSM3562051) <- vec_genes_GSM3562051
colnames(mat_GSM3562051) <- vec_barcode_GSM3562051

dim(mat_GSM3562051)
# [1] 28692  3793
spmat_GSM3562051 <- as(object = mat_GSM3562051, Class="sparseMatrix")
rm(mat_GSM3562051)

# metadata - received directly from authors
dirMetadata <- here("data", "expression", "mikkelsen2019")
filesMetadata <- dir(pattern="mapping", path = dirMetadata, full.names = T)

list_dt_metadata <- lapply(filesMetadata, fread)

# load normalized and filtered data for validation

filesdataNorm <- dir(pattern="normal", path = here("tmp-data", "expression"), full.names = T)
list_dt_dataNorm <- lapply(filesdataNorm,fread)

# ======================================================================= #
# ========================== Clean up metadata ================================ #
# ======================================================================= #

strsplit(x=dir(pattern="mapping", path = dirMetadata, full.names = F), split="-") %>%  sapply(., function(x) x[1]) -> vec_clustPrefix 

list_dt_metadata <- mapply(function(dt_metadata, prefix) {
  # Remove the "0" cluster which corresponds to unassigned (README.md)
  dt_metadata <- dt_metadata[dbCluster!="0"]
  # Add the GABA, Glut or Nonneuronal prefix
  dt_metadata[,dbCluster:=paste0(prefix, "_", dbCluster),]
}, 
prefix= vec_clustPrefix, 
dt_metadata = list_dt_metadata, 
SIMPLIFY = F)

# merge the data tables
dt_metadata <- Reduce(x = list_dt_metadata, f=rbind)

colnames(dt_metadata) <- c("cell_id","cell_type")

dim(dt_metadata)
# [1] 5973    2

# ======================================================================= #
# ========================== Filter cells ============================= #
# ======================================================================= #

# the initial 7,218 cells (3,784 male and 3,434 female), 89 cells with less than 500 UMIs or >40% of mitochondrial reads were discarded. 
# Gene expression of the remaining 7,129 cells was normalized .. 

# column-bind the count data
all.equal(rownames(spmat_GSM3562050),rownames(spmat_GSM3562051))
# [1] TRUE

spmat <- cbind(spmat_GSM3562050,spmat_GSM3562051)
#rm(spmat_GSM3562050,spmat_GSM3562051)

dim(spmat)
# [1] 28692  7232

# check that metadata cell barcodes correspond to those in expression data
table(colnames(spmat_GSM3562050) %in% dt_metadata$cell_id)
# FALSE  TRUE
# 580  2859

table(colnames(spmat_GSM3562051) %in% dt_metadata$cell_id)
# FALSE  TRUE
# 3779    14

# spmat_GSM3562051 barcode suffixes were changed when combining the matrices to avoid aliasing
# within male or female samples respectively the barcodes are unique
table(duplicated(colnames(spmat_GSM3562050)))
# FALSE
# 3439

table(duplicated(colnames(spmat_GSM3562051)))
# FALSE
# 3793

# but when merging:
table(duplicated(colnames(spmat)))
# FALSE  TRUE
# 7214    18

substr(x=colnames(spmat_GSM3562051),18,18) %>% as.numeric %>% summary
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 1       1       1       1       1       1

# change the names of spmat_GSM3562051
colnames(spmat_GSM3562051) <- gsub("-1","-2",colnames(spmat_GSM3562051))

# merge 
spmat <- cbind(spmat_GSM3562050,spmat_GSM3562051)

table(colnames(spmat) %in% dt_metadata$cell_id)
# FALSE  TRUE
# 1259  5973

# now we have matched a raw barcode to every metadata barcode
# let's filter on that 
spmat <- spmat[,colnames(spmat) %in% dt_metadata$cell_id]

# let#s reorder the metadata rows to match
all.equal(colnames(spmat), dt_metadata$cell_id)
# [1] 5972 string mismatches

dt_metadata <- dt_metadata[match(colnames(spmat),dt_metadata$cell_id),]

all.equal(colnames(spmat), dt_metadata$cell_id)
# [1] TRUE
# ======================================================================= #
# =============================== Filter genes =========================== #
# ======================================================================= #
# Genes with at least 2 counts in 5 cells were used for downstream analysis. 

spmat %>% as.matrix %>% '>='(2) %>% rowSums %>% '>='(5) -> vec_logicalGenesOK

table(vec_logicalGenesOK)
# FALSE  TRUE
# 15768 12924

spmat <- spmat[vec_logicalGenesOK,]

dim(spmat)
# [1] 12756  5973

spmat[0:4,0:3]
# [1] 4 x 3 sparse Matrix of class "dgCMatrix"
#                    AAACCTGAGAATAGGG-1 AAACCTGAGCCCTAAT-1 AAACCTGAGCCGCCTA-1
# ENSMUSG00000051951                  .                  .                  .
# ENSMUSG00000025902                  .                  .                  .
# ENSMUSG00000033845                  .                  .                  .
# ENSMUSG00000025903                  .                  1                  .
# ======================================================================= #
# ========================== Write to disk ================================ #
# ======================================================================= #

# n celltypes
dt_metadata %>% count(cell_type) %>% 
  write_csv(here("data","expression","mikkelsen2019","mikkelsen2019_cell_type.cell_type_metadata.csv"))

# counts
spmat %>% as.data.table -> dt_data # set rownames as column
dt_data <- data.table("gene"= rownames(spmat), dt_data)

dt_data[0:4,0:3]
#                  gene AAACCTGAGAATAGGG-1 AAACCTGAGCCCTAAT-1
# 1: ENSMUSG00000051951                  0                  0
# 2: ENSMUSG00000025902                  0                  0
# 3: ENSMUSG00000033845                  0                  0
# 4: ENSMUSG00000025903                  0                  1

file.out.data <- here("tmp-data","expression","mikkelsen2019.umi.csv")

data.table::fwrite(dt_data, file=file.out.data,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch 
R.utils::gzip(file.out.data, overwrite=TRUE) # gzip

# meta-data
file.out.meta <- here("tmp-data","expression","mikkelsen2019.metadata.csv")
data.table::fwrite(dt_metadata, file=file.out.meta,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch
