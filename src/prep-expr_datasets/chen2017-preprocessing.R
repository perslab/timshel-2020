############### SYNOPSIS ###################
# Download and pre=process chen-cellreports-2017 hypothalamus gene expression data

### OUTPUT:
# ....

### REMARKS:
# ....

### REFERENCE:

# https://www.cell.com/cell-reports/fulltext/S2211-1247(17)30321-2#secsectitle0150

# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library("here")
library("tidyverse")
library("Seurat")
library("openxlsx")

# ======================================================================= #
# ========================== Download ================================ #
# ======================================================================= #
### You only need to run this step once

if (!dir.exists(here("tmp-data"))) dir.create(here("tmp-data"))
if (!dir.exists(here("tmp-data","expression"))) dir.create(here("tmp-data", "expression"))
if (!dir.exists(here("data","expression"))) dir.create(here("data", "expression"))
if (!dir.exists(here("data","expression", "chen2017"))) dir.create(here("data", "expression", "chen2017"))

# raw data
rawDataDownload <- here("tmp-data", "expression","chen2017-GSE87544_Merged_17samples_14437cells_count.txt.gz")

if (!file.exists(rawDataDownload)) {
  # Download UMI data
  downloadURLraw <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87544/suppl/GSE87544_Merged_17samples_14437cells_count.txt.gz"
  download.file(downloadURLraw, destfile=rawDataDownload)
}

df.data_raw <- read.table(gzfile(rawDataDownload), header=T, sep="\t", stringsAsFactors=F, row.names = 1)

dim(df.data_raw)
# 23284 14437

# metadata - cluster assignments
metadataDownload <- here("data", "expression", "chen2017", "GSE87544_1443737Cells.SVM.cluster.identity.renamed.csv.gz")

if (!file.exists(metadataDownload)) {
  # Download UMI data
  downloadURLmeta <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE87nnn/GSE87544/suppl/GSE87544_1443737Cells.SVM.cluster.identity.renamed.csv.gz"
  download.file(downloadURLmeta, destfile=metadataDownload)
}

df.metadata <- read.csv(gzfile(metadataDownload), header=T, stringsAsFactors=F, row.names=1)
  
# get supp tab 2 to select cell clusters to use

suppTab1Download <- here("data", "expression", "chen2017", "mmc2.xlsx")

if (!file.exists(suppTab1Download)) {
  downloadURLsuppTab1 <- "https://www.cell.com/cms/10.1016/j.celrep.2017.03.004/attachment/86ebb5f9-2b5d-40ab-a4d9-b705b9e142dc/mmc2.xlsx"
  download.file(downloadURLsuppTab1, destfile=suppTab1Download)
}

df.suppTab1 <- openxlsx::read.xlsx(xlsxFile = suppTab1Download)

# ======================================================================= #
# ========================== Create Seurat object ============================= #
# ======================================================================= #

# Document S1 supp exp procedures
# https://www.cell.com/cell-reports/fulltext/S2211-1247(17)30321-2#secsectitle0150

seurat_obj <- CreateSeuratObject(counts = df.data_raw, 
                                 min.features = 0, 
                                 project = "chen-cellreports-2017", 
                                 meta.data = df.metadata)

seurat_obj
# An object of class Seurat 
# 23284 features across 14437 samples within 1 assay 
# Active assay: RNA (23284 features)

# ======================================================================= #
# =========================== Filter cells ============================== #
# ======================================================================= #

# https://www.cell.com/cell-reports/fulltext/S2211-1247(17)30321-2#secsectitle0150
# Based on the above clustering results, we re-assigned each of the 14,437 single cells (≥800 transcripts detected) to the 45 cell clusters. 

# All cell clusters were then pooled together and clusters with less than 10 cells, or
# representing double-droplets or without a marker identified or out of hypothalamus
# were filtered out. At the end, 45 cell clusters with distinct transcriptional features
# were identified.

seurat_obj_sub <- subset(x = seurat_obj, subset=SVM_clusterID %in% df.suppTab1$final_ClusterID & nCount_RNA>=800)

# this subset excludes two clusters 
# > table(seurat_obj$SVM_clusterID[!seurat_obj$SVM_clusterID %in% df.suppTab1$final_ClusterID])
# 
# SCO zothers
# 32    2348

dim(seurat_obj_sub)
# [1] 23284 12055

df.suppTab1$final_ClusterID %>% unique %>% length
# [1] 45

sum((seurat_obj_sub$SVM_clusterID %>% unique) %in% df.suppTab1$final_ClusterID)
# [1] 45

dim(seurat_obj_sub)
# [1] 23284 12055
# ======================================================================= #
# =============== EXPORT CELL-TYPE/ANNOTATION METADATA TO CSV =========== #
# ======================================================================= #

seurat_obj_sub@meta.data %>% count(SVM_clusterID) %>%
  write_csv(here("data","expression","chen2017","chen2017_SVM_clusterID.cell_type_metadata.csv"))

# ======================================================================= #
# ================================ EXPORT TO CSV ============================= #
# ======================================================================= #

### Write raw counts to temp
df <- as.data.frame(as.matrix(GetAssayData(seurat_obj_sub, slot="counts")))
dim(df) # 23284  2269
df <- df %>% rownames_to_column(var="gene") %>% select(gene, everything()) %>% as_tibble() # set rownames as column
file.out.data <- here("tmp-data","expression","chen2017.umi.csv")
data.table::fwrite(df, file=file.out.data,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch
R.utils::gzip(file.out.data, overwrite=TRUE) # gzip

### Write cell meta-data to temp
df.metadata <- seurat_obj_sub@meta.data %>% as_tibble()
df.metadata$cell_id <- rownames(seurat_obj_sub@meta.data)
file.out.meta <- here("tmp-data","expression","chen2017.metadata.csv")
data.table::fwrite(df.metadata, file=file.out.meta,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch
