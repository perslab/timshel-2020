############### SYNOPSIS ###################
# Download and prepare  romanov-natureneuroscience-2017 hypothalamus gene expression data

### OUTPUT:
# ....

### REMARKS:
# ....

### REFERENCE:
# https://www.nature.com/articles/nn.4462#methods
# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library("here")
library("tidyverse")
library("Seurat")

# ======================================================================= #
# ========================== Download ================================ #
# ======================================================================= #
### You only need to run this step once

if (!dir.exists(here("tmp-data"))) dir.create(here("tmp-data"))
if (!dir.exists(here("tmp-data","expression"))) dir.create(here("tmp-data", "expression"))
if (!dir.exists(here("data","expression"))) dir.create(here("data", "expression"))
if (!dir.exists(here("data","expression", "romanov2017"))) dir.create(here("data", "expression", "romanov2017"))

# ======================================================================= #
# ========================== Load files into memory ================================ #
# ======================================================================= #

# combined raw and metadatafile
combinedDataDownload <- here("data", "expression", "romanov2017", "romanov2017-GSE74672_expressed_mols_with_classes.xlsx.gz")

if (!file.exists(combinedDataDownload)) {
  # Download UMI data
  downloadURLcombined <-  "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74672/suppl/GSE74672_expressed_mols_with_classes.xlsx.gz"
  download.file(downloadURLcombined, destfile=combinedDataDownload)
  system2(command = "gunzip", args = combinedDataDownload)
  combinedDataDownload <- gsub("\\.gz","", combinedDataDownload)
}


#df.data_raw <- data.table::fread(rawDataDownload, nThread=24, showProgress=T)
mat.data_raw <- openxlsx::read.xlsx(combinedDataDownload, startRow=13,colNames=F,rowNames=T)

openxlsx::read.xlsx(combinedDataDownload, rows=1:12,colNames=T,rowNames=T) %>% 
  t %>% 
  as.data.frame -> df.metadata

vec_numericVars <- c("age (days postnatal)",
                     "sex (female=1,male=-1)",
                     "cell diameter",
                     "total molecules")

df.metadata[,vec_numericVars] <- apply(df.metadata[,vec_numericVars], 
                                       MARGIN=2, 
                                       FUN=as.numeric)

colnames(mat.data_raw) <- rownames(df.metadata)

# ======================================================================= #
# ========================= CLEAR UP METADATA =========================== #
# ======================================================================= #

colnames(df.metadata) <- gsub("\\ |,", "_",colnames(df.metadata))
colnames(df.metadata) <- gsub("\\(|\\)", "",colnames(df.metadata))
colnames(df.metadata) <- gsub("=", "",colnames(df.metadata))

df.metadata$level2_class_neurons_only<- gsub("\\ |,", "_",df.metadata$level2_class_neurons_only )
df.metadata$level2_class_neurons_only  <- gsub("\\(|\\)", "",df.metadata$level2_class_neurons_only )
df.metadata$level2_class_neurons_only  <- gsub("\\+", "pos",df.metadata$level2_class_neurons_only )
df.metadata$level2_class_neurons_only  <- gsub("-", "neg",df.metadata$level2_class_neurons_only )
df.metadata$level2_class_neurons_only  <- gsub("/", "_",df.metadata$level2_class_neurons_only )
df.metadata$level2_class_neurons_only  <- gsub("__", "_",df.metadata$level2_class_neurons_only )



# ======================================================================= #
# ====================== CREATE CELL_TYPE LABEL ========================== #
# ======================================================================= #
# include both neurons and glia

df.metadata$cell_type <- ifelse(!is.na(df.metadata$level2_class_neurons_only), df.metadata$level2_class_neurons_only, df.metadata$level1_class)

sum(is.na(df.metadata$cell_type))
#[1] 0

# ======================================================================= #
# ========================== Create Seurat object ============================= #
# ======================================================================= #


# Initialize a Seurat object and filter data
# Genes with <50 molecules in the whole data set [...] were excluded. 
## https://www.nature.com/articles/nn.4462#methods

seurat_obj <- CreateSeuratObject(counts = mat.data_raw, 
                              min.cells = 50, 
                              min.features = 0, 
                              meta.data=df.metadata)
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')


dim(seurat_obj)
#[1] 12884  2878
# ======================================================================= #
# =============================== FILTER GENES AND CELLS ========================== #
# ======================================================================= #

## https://www.nature.com/articles/nn.4462#methods

# Genes [...] expressed in >70% of the cells were excluded. 
GetAssayData(seurat_obj, slot= "counts") %>% as.matrix %>% '>'(0) %>% rowSums %>% '/'(ncol(seurat_obj)) -> vec_featsPctExpr

# Cells with more than 1,500 molecules/cell (excluding rRNA, mitochondrial RNA and repeats) were analyzed, resulting in a total of 3,131 cells. 

# We (Pers lab) also remove the "uc" class which appears to tag left overs (not included in  figures)
# We (Pers lab) also cell clusters with <= 3 cells
vec_smallClusters <- names(table(seurat_obj@meta.data$cell_type))[table(seurat_obj@meta.data$cell_type)<=3]

seurat_obj_sub <- subset(x = seurat_obj, 
                         features = rownames(seurat_obj)[vec_featsPctExpr<=0.7],
                         cells = colnames(seurat_obj)[!seurat_obj@meta.data$cell_type %in%  c("uc",vec_smallClusters)],
                         subset=nCount_RNA > 1500)

dim(seurat_obj_sub)
#[1] 12884  2730
# ======================================================================= #
# =============== EXPORT CELL-TYPE/ANNOTATION METADATA TO CSV =========== #
# ======================================================================= #

seurat_obj_sub@meta.data %>% count(cell_type) %>%
  write_csv(here("data","expression","romanov2017", "romanov2017_cell_type.cell_type_metadata.csv"))

# ======================================================================= #
# ================================ EXPORT TO CSV ============================= #
# ======================================================================= #

### Get expression data
df <- as.data.frame(as.matrix(GetAssayData(seurat_obj_sub, slot="counts")))
dim(df) # 26774 20921
df <- df %>% rownames_to_column(var="gene") %>% select(gene, everything()) %>% as_tibble() # set rownames as column
file.out.data <- here("tmp-data","expression", "romanov2017.umi.csv")
data.table::fwrite(df, file=file.out.data,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch
R.utils::gzip(file.out.data, overwrite=TRUE) # gzip

### Write cell meta-data
df.metadata <- seurat_obj_sub@meta.data %>% as_tibble()
df.metadata$cell_id <- rownames(seurat_obj_sub@meta.data)
file.out.meta <- here("tmp-data","expression", "romanov2017.metadata.csv")
data.table::fwrite(df.metadata, file=file.out.meta,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch
