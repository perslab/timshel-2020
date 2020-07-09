############### SYNOPSIS ###################
# Download and prepare moffitt-science-2018 hypothalamus gene expression data
# requires Seurat 3
### OUTPUT:
# ....

### REMARKS:
# ....

### REFERENCE:

# Supplementary materials
# https://science.sciencemag.org/highwire/filestream/717711/field_highwire_adjunct_files/0/aau5324-Moffitt-SM.pdf


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library("tidyverse")
library("openxlsx")
library("here")
library("Seurat")

# ======================================================================= #
# ========================== Download ================================ #
# ======================================================================= #
### You only need to run this step once

if (!dir.exists(here("tmp-data"))) dir.create(here("tmp-data"))
if (!dir.exists(here("tmp-data","expression"))) dir.create(here("tmp-data", "expression"))
if (!dir.exists(here("data","expression"))) dir.create(here("data", "expression"))
if (!dir.exists(here("data","expression", "moffitt2018"))) dir.create(here("data", "expression", "moffitt2018"))

# raw data
rawDataDownload <- here("tmp-data", "expression","moffitt2018-GSE113576_matrix.mtx.gz")

if (!file.exists(rawDataDownload)) {
  # Download UMI data
  downloadURLraw <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113576/suppl/GSE113576_matrix.mtx.gz"
  download.file(downloadURLraw, destfile=rawDataDownload)
}

spmat_data <- Matrix::readMM(rawDataDownload)
#spmat_data %>% as.matrix %>% as.data.frame -> df.data_raw

# barcodes
barcodesDownload <- here("tmp-data", "expression","moffitt2018-GSE113576_barcodes.tsv.gz")

if (!file.exists(barcodesDownload)) {
  downloadURLbarcodes <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113576/suppl/GSE113576_barcodes.tsv.gz"
  download.file(downloadURLbarcodes, destfile=barcodesDownload)
}

vec_barcodes <- read.delim(barcodesDownload, sep="\t", header= F) %>% '['(,1) %>% as.character

# genes 
genesDownload <- here("tmp-data", "expression","moffitt2018-GSE113576_genes.tsv.gz")

if (!file.exists(genesDownload)) {
  downloadURLgenes <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113576/suppl/GSE113576_genes.tsv.gz"
  download.file(downloadURLgenes, destfile=genesDownload)
}
vec_genes <- read.delim(genesDownload, sep = "\t", header = F) %>% '['(,1) %>% as.character

# metadata
metadataDownload <- here("data", "expression", "moffitt2018", "aau5324_Moffitt_Table-S1.xlsx")

if (!file.exists(metadataDownload)) {
  downloadURLmetadata <- "https://science.sciencemag.org/highwire/filestream/717711/field_highwire_adjunct_files/1/aau5324_Moffitt_Table-S1.xlsx"
  download.file(downloadURLmetadata, destfile=metadataDownload)
}

df.metadata <- openxlsx::read.xlsx(xlsxFile = metadataDownload, colNames = T, rowNames=F, startRow=2)

# add genes and barcodes to expression matrix 
colnames(spmat_data) <- vec_barcodes
rownames(spmat_data) <- vec_genes

# ======================================================================= #
# ========================== Clean up metadata ================================ #
# ======================================================================= #

colnames(df.metadata) <- gsub("\\(|\\)|-","", colnames(df.metadata))

df.metadata$cell_type <- ifelse(!is.na(df.metadata$Neuronal.cluster.determined.from.clustering.of.inhibitory.or.excitatory.neurons),
                               df.metadata$Neuronal.cluster.determined.from.clustering.of.inhibitory.or.excitatory.neurons,
                               df.metadata$Nonneuronal.cluster.determined.from.clustering.of.all.cells)
df.metadata$cell_type[is.na(df.metadata$cell_type)] <- df.metadata$Cell.class.determined.from.clustering.of.all.cells[is.na(df.metadata$cell_type)]

df.metadata$cell_type <- gsub(":|/|\\ ","_",df.metadata$cell_type)
df.metadata$cell_type <- gsub("__","_",df.metadata$cell_type)
rownames(df.metadata) <- df.metadata$Cell.name

all.equal(vec_barcodes, df.metadata$Cell.name)
#[1] TRUE

# ======================================================================= #
# ========================== Create Seurat object ============================= #
# ======================================================================= #

seurat_obj <- CreateSeuratObject(counts = spmat_data, 
                                 min.features = 0, 
                                 meta.data = df.metadata)

dim(seurat_obj)
#[1] 27998 31299
# ======================================================================= #
# =============================== Filter data =========================== #
# ======================================================================= #

# Supplementary materials
# https://science.sciencemag.org/highwire/filestream/717711/field_highwire_adjunct_files/0/aau5324-Moffitt-SM.pdf

# We found that the vast majority of clusters were not sensitive to
# small changes in this threshold. Cells that were members of unstable clusters or which clustered
# into small clusters (<10 cells) were marked as unstable and discarded from subsequent analysis.

# After clustering, we identified several clusters that expressed mixtures of markers for
# multiple cell types, e.g. neurons and oligodendrocytes, and that did not contain any genes that
# were uniquely enriched in these clusters (as judged by the largest z-score of individual genes for
#                                           the average expression profile of individual clusters). These clusters were marked as ambiguous,
# possibly representing doublets, and were excluded from subsequent analysis. Two inhibitory clusters (i33 and i34) were identified as putative doublets based on the
# strong expression of non-neuronal markers and were removed from subsequent analysis. We did
# not relabel the remaining clusters hence the cluster IDs e18, i33, and i34 are missing from our
# final data set.

seurat_obj_sub <- subset(x = seurat_obj, cells=colnames(seurat_obj)[!seurat_obj@meta.data$cell_type %in% c("Unstable", "Ambiguous")])

dim(seurat_obj_sub)
# [1] 27998 31211
# ======================================================================= #
# ========================== Write to disk ================================ #
# ======================================================================= #

seurat_obj_sub@meta.data %>% count(cell_type) %>% 
  write_csv(here("data","expression","moffitt2018","moffitt2018_cell_type.cell_type_metadata.csv"))

# ======================================================================= #
# ================================ EXPORT TO CSV ============================= #
# ======================================================================= #

### Get expression data
df <- as.data.frame(as.matrix(GetAssayData(seurat_obj_sub, slot= "counts")))
dim(df) # 26774 20921
df <- df %>% rownames_to_column(var="gene") %>% select(gene, everything()) %>% as_tibble() # set rownames as column
file.out.data <- here("tmp-data","expression","moffitt2018.umi.csv")
data.table::fwrite(df, file=file.out.data,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch 
R.utils::gzip(file.out.data, overwrite=TRUE) # gzip

### Write cell meta-data
df.metadata <- seurat_obj_sub@meta.data %>% as_tibble()
df.metadata$cell_id <- rownames(seurat_obj_sub@meta.data)
file.out.meta <- here("tmp-data","expression","moffitt2018.metadata.csv")
data.table::fwrite(df.metadata, file=file.out.meta,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch 
