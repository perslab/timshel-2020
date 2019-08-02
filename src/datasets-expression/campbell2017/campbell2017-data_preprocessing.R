############### SYNOPSIS ###################
# Download and prepare Campbell2017 hypothalamus gene epxression data

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(here)
library(tidyverse)
library(Seurat)


# ======================================================================= #
# ========================== Download Campbell ================================ #
# ======================================================================= #
### You only need to run this step once

# Download UMI data  
download_url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93374/suppl/GSE93374_Merged_all_020816_DGE.txt.gz"
file_download <- here("tmp-data/expression/campbell2017-GSE93374_Merged_all_020816_DGE.txt.gz")
download.file(download_url, destfile=file_download)

# Read data
df.data_raw <- read.table(gzfile(file_download), header=T, sep="\t", stringsAsFactors=F)
# save(df.data_raw, file=here("tmp-data/expression/GSE93374_Merged_all_020816_DGE.RData"))


# ======================================================================= #
# ================================ Campbell ================================ #
# ======================================================================= #

# Initialize a Seurat object
seurat_obj <- CreateSeuratObject(raw.data = df.data_raw, min.cells = 0, min.genes = 0)
seurat_obj # 26774 genes across 21086 samples.

# ======================================================================= #
# ================================ META DATA ================================ #
# ======================================================================= #

df.taxon_mapping <- read_csv(here("data/expression/campbell2017/campbell_taxon_mapping.csv"))
## Add cell type information to Seurat object
df.metadata <- read_csv(here("data/expression/campbell2017/campbell2017.cell_metadata-181210.csv")) %>% 
  select(-cell_type_glia, -cell_type_neurons, -umi_count_per_cell,-n_genes_expressed_per_cell, -cell_barcode)



### Match cell-type and add to data frame
idx.match <- match(df.metadata$cell_type_all_lvl1, df.taxon_mapping$cell_type_id)
anyNA(idx.match) # FALSE --> All matches.
df.metadata <- df.metadata %>% mutate(
  taxonomy_lvl1 = df.taxon_mapping$taxonomy_lvl1[idx.match],
  taxonomy_lvl2 = df.taxon_mapping$taxonomy_lvl2[idx.match]
)
df.metadata

### Confirm that all meta data is in the Seurat object
all(df.metadata$cell_id %in% seurat_obj@cell.names) # --> TRUE

# Add data to Seurat object
df.metadata.seurat <- df.metadata %>% as.data.frame()
rownames(df.metadata.seurat) <- df.metadata.seurat$cell_id
head(df.metadata.seurat)
seurat_obj <- AddMetaData(object=seurat_obj, metadata=df.metadata.seurat)
# ^ row names must correspond exactly to the items in object@@cell.names (but does not need to be in the same order)

# ======================================================================= #
# ==================== Remove cells with missing annotations ============ #
# ======================================================================= #

cells_with_annotations <- seurat_obj@meta.data %>% filter(cell_type_all_lvl1!="miss") %>% pull(cell_id)
nrow(seurat_obj@meta.data) - length(cells_with_annotations) # 165 cells will be discarded
seurat_obj <- SubsetData(seurat_obj, cells.use=cells_with_annotations, subset.raw=T)

# ======================================================================= #
# ====================  ============ #
# ======================================================================= #

# REPLACE "/" with "-" and "\s+" with "_" in CELL TYPE NAMES
seurat_obj@meta.data$cell_type_all_lvl1 <- stringr::str_replace_all(seurat_obj@meta.data$cell_type_all_lvl1, c("/"="-", "\\s+"="_"))
seurat_obj@meta.data$cell_type_all_lvl2 <- stringr::str_replace_all(seurat_obj@meta.data$cell_type_all_lvl2, c("/"="-", "\\s+"="_"))
seurat_obj@meta.data$taxonomy_lvl2 <- stringr::str_replace_all(seurat_obj@meta.data$taxonomy_lvl2, c("/"="-", "\\s+"="_"))
seurat_obj@meta.data$taxonomy_lvl1 <- stringr::str_replace_all(seurat_obj@meta.data$taxonomy_lvl1, c("/"="-", "\\s+"="_"))

# ======================================================================= #
# ================================ COUNT =============================== #
# ======================================================================= #

seurat_obj@meta.data %>% count(cell_type_all_lvl1)
# a01.Oligodend3         392
# a02.Oligodend2         131
# a03.EndothelialCells   240
# a04.MuralCells          84
# a05.Oligodend1         224
# a06.NG2/OPC            627
# a07.PVMMicro           330
# a08.Fibroblast         172
# a09.Ependymocytes      467
# a10.Astrocyte           99
# a11.Tanycyte1         1184
# a12.Tanycyte2         3504
# a13.Neurons1            37
# a14.Neurons2           502
# a15.Neurons3           693
# a16.Neurons4           533
# a17.Neurons5           799
# a18.Neurons6         10515
# a19.ParsTuber1         150
# a20.ParsTuber2         238


### April 13th 2019:
# n03	30 ---> SHOULD BE n03.Th/Sst
# n07	279 ---> SHOULD BE n07.Arx/Nr5a2
# n08	287 ---> SHOULD BE n08.Th/Slc6a3
# n27	431 ---> SHOULD BE n27.Tbx19


seurat_obj@meta.data %>% count(cell_type_all_lvl2)
# n01.Hdc	50
# n02.Gm8773/Tac1	83
# n03	30 ---> SHOULD BE n03.Th/Sst
# n04.Sst/Nts	70
# n05.Nfix/Htr2c	35
# n06.Oxt	36
# n07	279 ---> SHOULD BE n07.Arx/Nr5a2
# n08	287 ---> SHOULD BE n08.Th/Slc6a3
# n09.Th/Slc6a3	202
# n10.Ghrh	504
# n11.Th/Cxcl12	327
# n12.Agrp/Sst	539
# n13.Agrp/Gm8773	1747
# n14.Pomc/Ttr	512
# n15.Pomc/Anxa2	369
# n16.Rgs16/Vip	474
# n17.Rgs16/Dlx1	106
# n18.Rgs16/Slc17a6	258
# n19.Gpr50	67
# n20.Kiss1/Tac2	657
# n21.Pomc/Glipr1	310
# n22.Tmem215	752
# n23.Sst/Unc13c	552
# n24.Sst/Pthlh	246
# n25.Th/Lef1	229
# n26.Htr3b	486
# n27	431 ---> SHOULD BE n27.Tbx19
# n28.Qrfp	66
# n29.Nr5a1/Adcyap1	745
# n30.Nr5a1/Nfib	316
# n31.Fam19a2	171
# n32.Slc17a6/Trhr	342
# n33.unassigned(1)	420
# n34.unassigned(2)	1381
# s01.Ependymocy1	37
# s02.Ependymocy2	430
# s03.Oligodendro1	23
# s04.Oligodendro2	12
# s05.Oligodendro3	96
# s06.Oligodendro4	392
# s07.Oligodendro5	39
# s08.Pars_Tuber1A	42
# s09.Pars_Tuber1B	15
# s10.Pars_Tuber1C	93
# s11.Endothelial	240
# s12.Mural_Cells1	29
# s13.Mural_Cells2	55
# s14.PVMs	26
# s15.Microglia	304
# s16.Fibroblasts1	14
# s17.Fibroblasts2	38
# s18.Fibroblasts3	120
# s19.Parstuber2A	230
# s20.Parstuber2B	8
# s27.Oligodendro6	185
# s28.NG2_OPC1	23
# s29.NG2_OPC2	604
# s30.b2_tanycytes1	592
# s31.b2_tanycytes2	592
# s32.Astrocytes	99
# s33.a1_tanycytes2	68
# s34.a1_tanycytes1	334
# s35.b1_tanycytes	1933
# s36.a2_tanycytes	1169

# ======================================================================= #
# ========================= EXPORT Seurat object ======================== #
# ======================================================================= #

### Optional
# save(seurat_obj, file=here("tmp-data/expression/campbell2017.annotated181210.seurat_obj.RData"))

# ======================================================================= #
# =============== EXPORT CELL-TYPE/ANNOTATION METADATA TO CSV =========== #
# ======================================================================= #

seurat_obj@meta.data %>% count(taxonomy_lvl1, taxonomy_lvl2, cell_type_all_lvl1) %>% 
  write_csv(here("data/expression/campbell2017/campbell_lvl1.cell_type_metadata.csv"))
seurat_obj@meta.data %>% count(taxonomy_lvl1, taxonomy_lvl2, cell_type_all_lvl2) %>%
  write_csv(here("data/expression/campbell2017/campbell_lvl2.cell_type_metadata.csv"))

# ======================================================================= #
# ================================ EXPORT TO CSV ============================= #
# ======================================================================= #

### Get expression data
df <- as.data.frame(as.matrix(seurat_obj@raw.data))
dim(df) # 26774 20921
df <- df %>% rownames_to_column(var="gene") %>% select(gene, everything()) %>% as.tibble() # set rownames as column
file.out.data <- here("tmp-data/expression/campbell2017.umi.csv")
data.table::fwrite(df, file=file.out.data,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch 
R.utils::gzip(file.out.data, overwrite=TRUE) # gzip

### Write cell meta-data
df.metadata <- seurat_obj@meta.data %>% as.tibble()
file.out.meta <- here("tmp-data/expression/campbell2017.cell_metadata.csv")
data.table::fwrite(df.metadata, file=file.out.meta,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch 







