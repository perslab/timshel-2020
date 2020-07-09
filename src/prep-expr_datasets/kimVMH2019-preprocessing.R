############### SYNOPSIS ###################
# Download and prepare  kim-cell-2019 ventrolateral subdivision of the ventromedial hypothalamus (VMHvl) gene expression data

### OUTPUT:
# ....

### REMARKS:
# ....

### REFERENCE:
# https://www.sciencedirect.com/science/article/pii/S0092867419310712?via%3Dihub#app2
# 
# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library("here")
library("tidyverse")
library("data.table")

# ======================================================================= #
# ========================== Download ================================ #
# ======================================================================= #
### You only need to run this step once

if (!dir.exists(here("tmp-data"))) dir.create(here("tmp-data"))
if (!dir.exists(here("tmp-data","expression"))) dir.create(here("tmp-data", "expression"))
if (!dir.exists(here("data","expression"))) dir.create(here("data", "expression"))
if (!dir.exists(here("data","expression", "kimVMH2019"))) dir.create(here("data", "expression", "kimVMH2019"))

# ======================================================================= #
# ========================== Load files into memory ================================ #
# ======================================================================= #

# combined raw and metadatafile
combinedDataDownload <- here("data", "expression", "kimVMH2019", "single-cell-rna-seq-data-for-_multimodal-analysis-of-cell-types-in-a-hypothalamic-node-controlling-social-behavior_.zip")

if (!file.exists(combinedDataDownload) & !file.exists(gsub("\\.gz","", combinedDataDownload))) {
  # Download UMI data
  downloadURLcombined <-  "https://data.mendeley.com/archiver/ypx3sw2f7c?version=2"
  download.file(downloadURLcombined, destfile=combinedDataDownload)
  system2(command = "unzip", args = c("-f", combinedDataDownload, paste0("-d ", here("data", "expression", "kimVMH2019"))))
}

# 10x files 

zip_archives_10x <- dir(pattern="10x", path=here("data", "expression", "kimVMH2019"), full.names = T)

for (zip_arch in zip_archives_10x) {
  system2(command = "unzip", args=c("-f", zip_arch, paste0("-d ", here("data", "expression", "kimVMH2019"))))
}

dirs_10x <- dir(pattern = "^2018_\\d{4}|^2019_\\d{4}", path=here("data", "expression", "kimVMH2019"), full.names = T)

list_spmat_10x <- lapply(dirs_10x, Seurat::Read10X)

# 10x metadata
dt_10x_meta <- fread(here("data", "expression", "kimVMH2019", "10x_VMH_metadata.csv"))
dim(dt_10x_meta)
#[1] 41580    16

# SMART-seq files
zip_archives_smartseq <- dir(pattern="SMART-seq_VMH_cpm.rda.zip", path=here("data", "expression", "kimVMH2019"), full.names = T)

system2(command = "unzip", args=c("-f", zip_archives_smartseq, paste0("-d ", here("data", "expression", "kimVMH2019"))))

smartseq_file <- here("data", "expression", "kimVMH2019","cpm.rda")

# create a temporary new environment to load a 'rda' file into
env <- new.env()
nm <- load(smartseq_file, env)[1]
spmat_smartseq <- env[[nm]] 
rm(env)

dim(spmat_smartseq)
# [1] 30862  4574

# load SMART-seq metadata 

dt_smartseq_meta <- fread(here("data", "expression", "kimVMH2019", "SMART-seq_VMH_metadata.csv"))
dim(dt_smartseq_meta)
# [1] 4574   22


# ======================================================================= #
# ========================== Clean up files ================================ #
# ======================================================================= #
# SMART-seq
all.equal(dt_smartseq_meta$sample_name, colnames(spmat_smartseq))
# [1] "4573 string mismatches"

all(dt_smartseq_meta$sample_name %in% colnames(spmat_smartseq))
#[1] TRUE

spmat_smartseq <- spmat_smartseq[,match(dt_smartseq_meta$sample_name, colnames(spmat_smartseq))]

spmat_smartseq[0:3,0:3]
# SM-GE4R1_S084_E1-50 SM-GE4R1_S150_E1-50 SM-GE4R1_S151_E1-50
# 0610005C13Rik               .                  .                    .
# 0610006L08Rik               .                  .                    .
# 0610007P14Rik             165.264            109.0218             136.043

dt_smartseq_meta[,cell_type:=gsub(" ","_",smart_seq_cluster_label),]
colnames(dt_smartseq_meta)[1] <- "cell_id"

dt_smartseq_meta$cell_type %>% table
# Astro_Aqp4           Clic4      Dlk1_1_Il1rapl2
# 18                   12                  210
# Dlk1_2_Nts          Dlk1_3_Maob          Dlk1_4_Maob
# 146                  258                  109
# Dlk1_5_Egflam         Dlk1_6_F2rl2            Endo_Gja5
# 139                   30                    6
# Esr1_1_Insm2   Esr1_2_Samd9l_Rorb  Esr1_3_Samd9l_Nmur2
# 106                   95                  100
# Esr1_4_Smoc2           Esr1_5_Ngf         Esr1_6_Ctxn3
# 150                  148                  122
# Irx5               Lhx1_1               Lhx1_2
# 13                   25                   11
# Lhx9        Micro_Siglech        Nfib_1_Prdm13
# 21                   48                   91
# Nfib_2_Car4       Nr5a1_1_Npffr2        Nr5a1_2_Ctxn3
# 26                  200                  170
# Nr5a1_3_Car8     Nr5a1_4|7_Glipr1       Nr5a1_8_Tcf7l2
# 155                  420                  175
# Nr5a1_9|11_Rorb Nr5a1_Foxp2_1_Prdm13  Nr5a1_Foxp2_2_C1ql2
# 74                  150                  216
# Nr5a1_Foxp2_3      Nr5a1_Smoc2_low              Nup62cl
# 14                   26                   91
# Oligo_Opalin_1       Oligo_Opalin_2             PVM_Mrc1
# 15                   12                    7
# Satb2_Calcrl           Satb2_Six3                 Scgn
# 303                   61                   53
# Slc18a3    Slc32a1_1_Nup62cl            Slc32a1_2
# 47                   11                  149
# SMC          Tbx3_1_Tac2         Tbx3_2_Anxa2
# 7                   14                   63
# Tbx3_3_Slc14a1   Tsix_Esr1_1_Gabrg3   Tsix_Esr1_2_Twist2
# 128                   65                   64

dt_smartseq_meta[,cell_type:=gsub("\\|","_",cell_type),] # get rid of REGEX special character


# 10x 
dt_10x_meta[,tv_cluster_label,] %>% table

# Dlk1_1        Dlk1_2        Dlk1_3        Dlk1_4        Dlk1_5
# 2046          1830          2124          2317          1487
# Dlk1_6      Doublets        Esr1_1      Esr1_2|3        Esr1_4
# 1304           195          1310          1390          1059
# Esr1_5        Esr1_6        Esr1_7       Nr5a1_1      Nr5a1_11
# 1070          1210           500          1223          1041
# Nr5a1_2       Nr5a1_3     Nr5a1_4|5       Nr5a1_6       Nr5a1_7
# 1794          1003          3827          1002          1050
# Nr5a1_8    Nr5a1_9|10 Nr5a1_Foxp2_1 Nr5a1_Foxp2_2       Nup62cl
# 838          1018          1500          1571          1201
# Satb2_1       Satb2_2       Satb2_3          Scgn     Tsix_Esr1
# 1567          2132           798           514          1659

dt_10x_meta[,tv_cluster_label:=gsub("\\|","_",tv_cluster_label),] # get rid of REGEX special character

colnames(dt_10x_meta)[colnames(dt_10x_meta)=="sample_name"] <- "cell_id"
colnames(dt_10x_meta)[colnames(dt_10x_meta)=="tv_cluster_label"] <- "cell_type"

# merge 10x matrices

# attach sample_specific prefix from metadata to expression matrix barcodes

dt_10x_meta$cell_id %>% gsub("[A-Z]","",.) %>% unique -> vec_10x_prefix 

list_spmat_10x <- lapply(list_spmat_10x, function(spmat_10x) {
  vec_propMatches <- sapply(vec_10x_prefix, function(prefix) {
    grep(prefix, dt_10x_meta$cell_id, value=T) %>% gsub(prefix,"",.) -> vec_barcodes 
    sum(vec_barcodes %in% colnames(spmat_10x))/length(vec_barcodes)
  })
  prefix_to_add <- names(vec_propMatches)[which.max(vec_propMatches)]
  colnames(spmat_10x) <- paste0(prefix_to_add, colnames(spmat_10x))
  spmat_10x
})

# column bind all 10x matrices
spmat_10x <- Reduce(x=list_spmat_10x, f = function(spmat1, spmat2) {
  vec_commonGenes <- intersect(rownames(spmat1), rownames(spmat2))
  cbind(spmat1[vec_commonGenes,], spmat2[vec_commonGenes,])
})

dim(spmat_10x)
# [1]  27998 183278

spmat_10x[0:3,0:3]
# 3 x 3 sparse Matrix of class "dgCMatrix"
# 1_AAACCTGAGCTGATAA 1_AAACCTGCACTGTCGG 1_AAACCTGCATAGTAAG
# Xkr4                     .                  5                  1
# Gm1992                   .                  .                  .
# Gm37381                  .                  .                  .

# ======================================================================= #
# =============================== FILTER CELLS ========================== #
# ======================================================================= #

# SMART-seq: already filtered

# 10x
# Cells that met any one of the following criteria were filtered out for downstream processing in each 10x run: 
 
# < 600 detected genes (for UMI count > 0), 
# > 30,000 UMI counts (potential multiplets), or 
# the proportion of the UMI count attributable to mitochondrial genes was greater than 15%. 
# Doublets were further removed by first classifying cells into broad cell classes 
# (neuronal versus non-neuronal) based on the co-expression of any pair of their marker genes 
# (Stmn2 for neurons; Cldn5 for endothelial cells; C1qc for microglia; Opalin for oligodendrocytes; 
# Gja1 for astrocytes; Pdgfra for OPCs; Mustn1 for mural cells; see Figure S4A).

condition_filter <- quote(dt_10x_meta$cell_type !="Doublets")
dt_10x_meta <- dt_10x_meta[eval(condition_filter),]
dim(dt_10x_meta)
#[1] 41385    16

spmat_10x <- spmat_10x[,match(dt_10x_meta$cell_id, colnames(spmat_10x))]
dim(spmat_10x)
#[1] 27998 41385
# ======================================================================= #
# =============== EXPORT CELL-TYPE/ANNOTATION COUNT TO CSV =========== #
# ======================================================================= #

dt_10x_meta %>% count(cell_type) %>%
  write_csv(here("data","expression","kimVMH2019", "kimVMH2019_10x_cell_type.cell_type_metadata.csv"))

dt_smartseq_meta %>% count(cell_type) %>%
  write_csv(here("data","expression","kimVMH2019", "kimVMH2019_smartseq_cell_type.cell_type_metadata.csv"))

# ======================================================================= #
# ================================ EXPORT TO CSV ============================= #
# ======================================================================= #

# 10x 
spmat_10x %>% as.matrix %>% as.data.table -> dt_10x
dt_10x <- data.table("gene"=rownames(spmat_10x), dt_10x)
file.out.data.10x <- here("tmp-data","expression", "kimVMH2019_10x.umi.csv")
data.table::fwrite(dt_10x, file=file.out.data.10x,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch
R.utils::gzip(file.out.data.10x, overwrite=TRUE) # gzip

### Write cell meta-data
file.out.meta.10x <- here("tmp-data","expression", "kimVMH2019_10x.metadata.csv")
data.table::fwrite(dt_10x_meta, file=file.out.meta.10x,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch

# SMART-seq
### Get expression data

spmat_smartseq %>% as.matrix %>% as.data.table -> dt_smartseq
dt_smartseq <- data.table("gene"=rownames(spmat_smartseq), dt_smartseq)
file.out.data.smartseq <- here("tmp-data","expression", "kimVMH2019_smartseq.umi.csv")
data.table::fwrite(dt_smartseq, file=file.out.data.smartseq,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch
R.utils::gzip(file.out.data.smartseq, overwrite=TRUE) # gzip

### Write cell meta-data
file.out.meta.smartseq <- here("tmp-data","expression", "kimVMH2019_smartseq.metadata.csv")
data.table::fwrite(dt_smartseq_meta, file=file.out.meta.smartseq,  # fwrite cannot write gziped files
                   nThread=24, verbose=T) # write file ---> write to scratch
