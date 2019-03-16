############### SYNOPSIS ###################
# *SCRIPT NOT USED YET, BUT OK TO KEEP FOR INPIRATION FOR LOADING HYPOTHALAMUS DATA*
# Prepare hypothalamus gene epxression data
# Data from hypothalamus. All cells from Campbell et al., Chen et al., Romanov et al. and Pers Group
# COPY FROM /projects/mludwig/Dataset_Alignment/Data_preprocessing/Data_preprocessing_all_cells.Rmd
# DATE COPIED: December 10th, 2018

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(Seurat)
library(tidyverse)


wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-hypothalamus"
setwd(wd)

# ======================================================================= #
# ================================ Campbell ================================ #
# ======================================================================= #

## Create Seurat object
# Load data  
load(paste("/projects/timshel/sc-hypothalamus_atlas_app/data-SCGE/campbell2017/",
           "scge.umi.gene_names/GSE93374_Merged_all_020816_DGE.RData", sep=""))

# Initialize a Seurat object
# Keep all genes expressed in >= 3 cells. Keep all cells with at least 200 detected genes
campbell <- CreateSeuratObject(raw.data = df.dge.campbell2017, min.cells = 3, min.genes = 200)


## Add cell type information to Seurat object
# Load Campbell et al. meta data
campbell.meta <- read.csv(paste("/projects/timshel/sc-hypothalamus_atlas_app/",
                                "data-SCGE/metadata/Campbell2017.metavariables.csv", sep =""))

# Add meta data information
campbell.meta <- campbell.meta[, c("cell_id", "cell_type_all_lvl1")]
colnames(campbell.meta) <- c("cell_id", "cell_type")
campbell.meta <- data.frame(apply(campbell.meta, 2, as.character), stringsAsFactors=FALSE)
campbell@meta.data$cell_id <- rownames(campbell@meta.data)
campbell.meta <- left_join(campbell@meta.data, campbell.meta, by="cell_id")
campbell@meta.data$cell_type <- campbell.meta$cell_type


## Data preprocessing
# Compute the percent of mitochondrial RNA
mito.genes <- grep(pattern = "^MT-", x = rownames(x = campbell@data), value = TRUE,
                   ignore.case=TRUE)
percent.mito <- Matrix::colSums(campbell@raw.data[mito.genes,]) / Matrix::colSums(campbell@raw.data)

# Add data to Seurat object
campbell <- AddMetaData(object = campbell, metadata = percent.mito, col.name = "percent.mito")

# Filter cells
campbell <- FilterCells(object=campbell, subset.names=c("nGene", "percent.mito"), 
                        low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.2))

# Normalize the data
campbell <- NormalizeData(campbell, normalization.method = "LogNormalize", scale.factor = 1e4)

# Scale the data 
campbell <- ScaleData(campbell, vars.to.regress = c("nUMI", "percent.mito"))

## Save Seurat object 
save(campbell, file = "Campbell_all_cells_preprocessed.RData")




# ======================================================================= #
# ================================ Chen ================================ #
# ======================================================================= #

## Create Seurat object
# Load data  
chen.data <- read.table(gzfile(paste("/projects/timshel/sc-hypothalamus_atlas_app/data-SCGE/", 
                                     "chen2017/scge.umi.xiaoji_merged.gene_names/",
                                     "Merged_17samples_non_repetitive_name.dge.txt.gz", sep="")),
                        header = T)
rownames(chen.data) <- chen.data[,1]
chen.data <- chen.data[, -1]


# Initialize a Seurat object
# Keep all genes expressed in >= 3 cells. Keep all cells with at least 200 detected genes
chen <- CreateSeuratObject(raw.data = chen.data, min.cells = 3, min.genes = 200)


## Add cell type information to Seurat object
# Load Chen et al. meta data
chen.meta <- read.csv(paste("/projects/timshel/sc-hypothalamus_atlas_app/",
                            "data-SCGE/metadata/Chen2017.metavariables.csv", sep =""))

# Add meta data information
chen.meta <- chen.meta[, c("cell_id", "cell_type_all_lvl1")]
colnames(chen.meta) <- c("cell_id", "cell_type")
chen.meta <- data.frame(apply(chen.meta, 2, as.character), stringsAsFactors=FALSE)
chen@meta.data$cell_id <- rownames(chen@meta.data)
chen.meta <- left_join(chen@meta.data, chen.meta, by="cell_id")
chen@meta.data$cell_type <- chen.meta$cell_type


## Data preprocessing
# Filter cells
chen <- FilterCells(object = chen, subset.names = c("nGene"), 
                    low.thresholds = c(200), high.thresholds = c(5000))

# Normalize the data
chen <- NormalizeData(chen, normalization.method = "LogNormalize", scale.factor = 1e4)

# Scale the data
# Only nUMI is scaled out since mitochondrial RNAs have been removed from the Chen et al. data set  
chen <- ScaleData(chen, vars.to.regress = c("nUMI"))


## Save Seurat object 
save(chen, file = "Chen_all_cells_preprocessed.RData")


# ======================================================================= #
# ================================ Romanov ================================ #
# ======================================================================= #

## Create Seurat object
# Load data
romanov.data <- read.table(gzfile(paste("/projects/timshel/sc-hypothalamus_atlas_app/data-SCGE/",
                                        "romanov2016/scge.umi.gene_names/",
                                        "Romanov2017.dge_matrix.no_excel_errors.txt.gz", sep="")),
                           header = T)


# Initialize a Seurat object
# Keep all genes expressed in >= 3 cells. Keep all cells with at least 200 detected genes
romanov <- CreateSeuratObject(raw.data = romanov.data, min.cells = 3, min.genes = 200)


## Add cell type information to Seurat object
# Load Romanov et al. meta data
romanov.meta <- read.csv(paste("/projects/timshel/sc-hypothalamus_atlas_app/",
                               "data-SCGE/metadata/Romanov2017.metavariables.csv", sep =""))

# Add meta data information
romanov.meta <- romanov.meta[, c("cellID", "cell_type_all_lvl1")]
colnames(romanov.meta) <- c("cell_id", "cell_type")
romanov.meta <- romanov.meta[romanov.meta$cell_type != "", ]
romanov.meta$cell_id <- paste0("X", romanov.meta$cell_id)
romanov.meta <- data.frame(apply(romanov.meta, 2, as.character), stringsAsFactors=FALSE)


romanov@meta.data$cell_id <- rownames(romanov@meta.data)
romanov.meta <- left_join(romanov@meta.data, romanov.meta, by="cell_id")
romanov@meta.data$cell_type <- romanov.meta$cell_type

## Data preprocessing
# Filter cells
romanov <- FilterCells(object = romanov, subset.names = c("nGene"), 
                       low.thresholds = c(200), high.thresholds = c(5000))

# Normalize the data
romanov <- NormalizeData(romanov, normalization.method = "LogNormalize", scale.factor = 1e4)

# Scale the data
# Only nUMI is scaled out since mitochondrial RNAs have been removed from the Romanov et al. data set
romanov <- ScaleData(romanov, vars.to.regress = c("nUMI"))


## Save Seurat object 
save(romanov, file = "Romanov_all_cells_preprocessed.RData")

# ======================================================================= #
# ================================ Pers lab ================================ #
# ======================================================================= #

## Create Seurat object
# Load data  
pers.data.path <- paste("~/ygg-data/sc-10x/data-runs/170612-perslab-arc_lira",
                        "/agg-arc_lira-s1_12-no_norm/outs/filtered_gene_bc_matrices_mex/mm10", 
                        sep="")
pers.data <- Read10X(data.dir = pers.data.path)

length(na.omit(match(rownames(pers.data), rownames(df.dge.campbell2017))))


dim(pers.data)
colnames(pers.data)
# Initialize the Seurat object  
# Keep all genes expressed in >= 3 cells. Keep all cells with at least 200 detected genes
pers <- CreateSeuratObject(raw.data = pers.data, min.cells = 3, min.genes = 200)


## Add cell type information to Seurat object
# Load Seurat object with annotations 
load(paste("~/ygg-projects/timshel/sc-arc_lira/src/nb-cell_atlas-v2/",
           "arc_lira_cell_atlas-core_objs.annotated.RData", sep=""))

# Add meta data information
pers.meta <- cbind(rownames(seurat_obj@meta.data), seurat_obj@meta.data$annotation)
colnames(pers.meta) <- c("cell_id", "cell_type")
pers.meta <- data.frame(apply(pers.meta, 2, as.character), stringsAsFactors=FALSE)
pers@meta.data$cell_id <- rownames(pers@meta.data)
pers.meta <- left_join(pers@meta.data, pers.meta, by="cell_id")
pers@meta.data$cell_type <- pers.meta$cell_type


## Data preprocessing
# Compute the percent of mitochondrial RNA
mito.genes <- grep(pattern = "^mt-", x = rownames(x = pers@data), value = TRUE)
percent.mito <- Matrix::colSums(pers@raw.data[mito.genes, ])/Matrix::colSums(pers@raw.data)

# Add data to the Seurat object
pers <- AddMetaData(object = pers, metadata = percent.mito, col.name = "percent.mito")

# Filter outlier cells
pers <- FilterCells(object = pers, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(5000, 0.25))

# Normalize the data 
pers <- NormalizeData(object = pers, normalization.method = "LogNormalize", 
                      scale.factor = 1e4)

# Scale the data
pers <- ScaleData(pers, vars.to.regress = c("nUMI", "percent.mito"))


## Save Seurat object 
save(pers, file = "Pers_all_cells_preprocessed.RData")

