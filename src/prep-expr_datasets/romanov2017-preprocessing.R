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

if (!file.exists(combinedDataDownload) & !file.exists(gsub("\\.gz","", combinedDataDownload))) {
  # Download UMI data
  downloadURLcombined <-  "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74672/suppl/GSE74672_expressed_mols_with_classes.xlsx.gz"
  download.file(downloadURLcombined, destfile=combinedDataDownload)
  system2(command = "gunzip", args = combinedDataDownload)
}


#df.data_raw <- data.table::fread(rawDataDownload, nThread=24, showProgress=T)
mat.data_raw <- openxlsx::read.xlsx(gsub("\\.gz","", combinedDataDownload), startRow=13,colNames=F,rowNames=T)

openxlsx::read.xlsx(gsub("\\.gz","", combinedDataDownload), rows=1:12,colNames=T,rowNames=T) %>% 
  t %>% 
  as.data.frame -> df.metadata

vec_numericVars <- c("age (days postnatal)",
                     "sex (female=1,male=-1)",
                     "cell diameter",
                     "total molecules")

vec_charVars <- c("level1 class", 
                  "level2 class (neurons only)",
                  "level2 cluster number (neurons only)")

df.metadata[,vec_numericVars] <- apply(df.metadata[,vec_numericVars], 
                                       MARGIN=2, 
                                       FUN=as.numeric)

df.metadata[,vec_charVars] <- apply(df.metadata[,vec_charVars], 
                                       MARGIN=2, 
                                       FUN=as.character)

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

df.metadata$cell_type <- ifelse(!is.na(df.metadata$level2_class_neurons_only), 
                                df.metadata$level2_class_neurons_only, 
                                as.character(df.metadata$level1_class))

sum(is.na(df.metadata$cell_type))
#[1] 0

df.metadata$cell_type %>% table
# Adcyap1_1_Tac1                      Adcyap1_2
# 14                             20
# astrocytes                     Avp_1_high
# 267                              9
# Avp_2_high                   Avp_3_medium
# 8                              8
# circadian_1_Vip_Grppos_neg     circadian_2_Nms_VIPpos_neg
# 6                             12
# circadian_3_Per2                     Dopamine_1
# 13                             17
# Dopamine_2_low_VMAT2                     Dopamine_3
# 7                             12
# Dopamine_4                    endothelial
# 5                            240
# ependymal                        GABA_10
# 356                             20
# GABA_11_Nts_1                  GABA_12_Nts_2
# 25                             26
# GABA_13_Galanin               GABA_14_Npy_Agrp
# 23                             13
# GABA_15_Npynegmedium                 GABA_2_Gucy1a3
# 7                             19
# GABA_3_Crhpos_neg_Lhx6       GABA_4_Crhpos_neg_Pgr15l
# 20                             24
# GABA_5_Calcr_Lhx1               GABA_6_Otof_Lhx1
# 32                             24
# GABA_7_Pomcpos_neg                         GABA_8
# 27                             15
# GABA_9                          GABA1
# 22                             14
# Gadneglow_Gnrhneg_pos                           Ghrh
# 9                              3
# Hcrt                    Hmitpos_neg
# 14                             24
# microglia                           Npvf
# 48                              6
# oligos                     Oxytocin_1
# 1001                              9
# Oxytocin_2                     Oxytocin_3
# 7                             11
# Oxytocin_4                           Pmch
# 10                              4
# Qrfp                      Sst_1_low
# 3                             17
# Sst_2_high                   Sst_3_medium
# 9                             12
# Trh_1_low                   Trh_2_medium
# 26                             20
# Trh_3_high                             uc
# 9                            126
# Vglut2_1_Penk         Vglut2_10_Morn4_Prrc2a
# 2                              3
# Vglut2_11               Vglut2_12_Mgat4b
# 12                              4
# Vglut2_13Ninl_Rfx5_Zfp346               Vglut2_14_Col9a2
# 5                              5
# Vglut2_15_Hcn1_6430411K18Rik           Vglut2_16_Gm5595_Tnr
# 4                              3
# Vglut2_17_A930013F10Rik_Pou2f2      Vglut2_18_Zfp458_Ppp1r12b
# 3                              2
# Vglut2_2_Crhpos_neg     Vglut2_3_Crhnegpos_neg_low
# 15                             40
# Vglut2_4             Vglut2_5_Myt1_Lhx9
# 4                              3
# Vglut2_6_Prmt8_Ugdh            Vglut2_7_Pgam_Snx12
# 11                              6
# Vglut2_8                Vglut2_9_Gpr149
# 9                              6
# vsm
# 71


# ======================================================================= #
# ========================== Create Seurat object ============================= #
# ======================================================================= #
# Initialize a Seurat object and filter data

seurat_obj <- CreateSeuratObject(counts = mat.data_raw, 
                              min.cells = 0, 
                              min.features = 0, 
                              meta.data=df.metadata)
# Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')


dim(seurat_obj)
#[1] 24341  2881
# ======================================================================= #
# =============================== FILTER GENES AND CELLS ========================== #
# ======================================================================= #

## https://www.nature.com/articles/nn.4462#methods

# Genes [...] expressed in >70% of the cells were excluded. 
GetAssayData(seurat_obj, slot= "counts") %>% as.matrix %>% '>'(0) %>% rowSums %>% '/'(ncol(seurat_obj)) -> vec_featsPctExpr

summary(vec_featsPctExpr)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000000 0.001041 0.029851 0.122278 0.183270 1.000000

# Genes with <50 molecules in the whole data set [..] were excluded
GetAssayData(seurat_obj, slot= "counts") %>% as.matrix %>% rowSums %>% '>='(50) -> vec_logicalOver50Counts

table(vec_logicalOver50Counts)
# FALSE  TRUE
# 10568 13773

# Cells with more than 1,500 molecules/cell (excluding rRNA, mitochondrial RNA and repeats) were analyzed, resulting in a total of 3,131 cells. 
# We (Pers lab) also remove the "uc" class which appears to tag left overs (not included in  figures)
# We (Pers lab) also cell clusters with <= 3 cells
vec_smallClusters <- names(table(seurat_obj@meta.data$cell_type))[table(seurat_obj@meta.data$cell_type)<=3]

seurat_obj_sub <- subset(x = seurat_obj, 
                         features = rownames(seurat_obj)[vec_featsPctExpr<=0.7 & vec_logicalOver50Counts],
                         cells = colnames(seurat_obj)[!seurat_obj@meta.data$cell_type %in%  c("uc",vec_smallClusters)],
                         subset=nCount_RNA > 1500)

dim(seurat_obj_sub)
#[1] 13364  2733
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
