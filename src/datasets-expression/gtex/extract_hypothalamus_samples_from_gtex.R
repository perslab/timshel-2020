############### SYNOPSIS ###################
# Extact hypothalamus samples from GTEx data

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
library(Matrix)

# ======================================================================= #
# ============================  PARAMS  ============================== #
# ======================================================================= #

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-gtex/"
setwd(wd)

# ======================================================================= #
# ============================  XXXXXXXXX  ============================== #
# ======================================================================= #

### Load data
#/data/rna-seq/gtex/v7-seurat_objs/gtex.seurat_obj.gene_median_tpm.RData
#/data/rna-seq/gtex/v7-seurat_objs/gtex.seurat_obj.gene_tpm.RData
file.RData.cell_atlas <- "/data/rna-seq/gtex/v7-seurat_objs/gtex.seurat_obj.gene_tpm.RData" # seurat_obj
load(file.RData.cell_atlas) 

# ======================================================================= #
# ============================  XXXXXXXXX  ============================== #
# ======================================================================= #

# SMTS: tissue
# SMTSD: sub-tissue
# seurat_obj@meta.data$annotation <- seurat_obj@meta.data$SMTS # tissue
seurat_obj@meta.data$annotation <- seurat_obj@meta.data$SMTSD # sub-tissue

seurat_obj@meta.data %>% count(SMTS)
seurat_obj@meta.data %>% count(SMTSD) # 121 samples for 'Brain - Hypothalamus'

### Subset
seurat_obj.hypothalamus <- SubsetData(SetAllIdent(seurat_obj, id="annotation"), ident.use='Brain - Hypothalamus', subset.raw=T)

### Save object
save(seurat_obj.hypothalamus, file="gtex.hypothalamus.seurat_obj.gene_tpm.RData")


# ======================================================================= #
# ============================  XXXXXXXXX  ============================== #
# ======================================================================= #


















# ======================================================================= #
# ========================  NOT USED: READING RAW DATA  ====================== #
# ======================================================================= #


# Read expression data
dir.data <- "/data/rna-seq/gtex/v7"
file.data <- file.path(dir.data, "GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz")
df.data_raw <- read_delim(file.data, delim="\t", skip=2) # OBS: skip 2 first lines
head(df.data_raw)


### make data frame containing only gene expression values
df.expression <- df.data_raw %>% select(-Name, -Description) %>% as.data.frame()

## REMOVE ENSEMBL VERSION NUMBER FROM GENE ID.
gene_names <- sapply(str_split(df.data_raw$Name, pattern=fixed(".")), '[[', 1)
rownames(df.expression) <- gene_names # gene rownames. the rownames() are used by Seurat as gene names. 


# Read meta data:
file.annotations <- "/data/rna-seq/gtex/v7/GTEx_v7_Annotations_SampleAttributesDS.txt"
df.annot <- read_delim(file.annotations, delim="\t")
head(df.annot)


## view the meta data
# SMTS: tissue
# SMTSD: sub-tissue
df.annot %>% count(SMTSD) %>% arrange(desc(n))
df.annot %>% count(SMTS) %>% arrange(desc(n))


## Filtering sample annotations to keep only SAMPID in the expression data
df.annot.filtered <- df.annot %>% filter(SAMPID %in% colnames(df.expression))

