############### SYNOPSIS ###################
# AIM: show that the prioritization results are not confounded by scRNA-seq associated factors, in tabula muris
# N-cells
# N-UMI
# N-genes
# ...

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

.libPaths(c(.libPaths(), "~/R/x86_64-pc-linux-gnu-library/"))
library("vctrs", lib.loc = "~/R/x86_64-pc-linux-gnu-library/")
library("rlang", lib.loc = "~/R/x86_64-pc-linux-gnu-library/")
library("tidyverse")#, lib.loc = "~/R/x86_64-pc-linux-gnu-library/")
library(here)
library(corrr)
library(data.table)
# library(corrplot) # https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ========================== READ META DATA ============================= #
# ======================================================================= #

### Read
path.metadata.cell <- here("tmp-data", "expression", "tabula_muris.metadata.csv")
path.data <- here("tmp-data", "expression", "tabula_muris.umi.csv.gz")


df.data <- read_csv(path.data)
df.metadata.cell <- read_csv(path.metadata.cell)
df.metadata.clust <- get_metadata(dataset_prefix="tabula_muris")

all.equal(df.metadata.cell$cell_id, colnames(df.data)[-1])
#TRUE 

all(df.metadata.cell$tissue_celltype %in% df.metadata.clust$tissue_celltype)
#[1] TRUE

all(df.metadata.clust$tissue_celltype %in% df.metadata.cell$tissue_celltype)
#[1] TRUE

# ======================================================================= #
# ========================== count genes and  in data ============================= #
# ======================================================================= #


mat_genestats = sapply(df.metadata.clust$tissue_celltype, FUN = function(tissue_celltype) {
  mat_sub_raw <- df.data[,2:ncol(df.data)][,df.metadata.cell$tissue_celltype==tissue_celltype,]
  return(c("nCount_RNA"=median(colSums(mat_sub_raw)), "nCount_gene"=median(colSums(mat_sub_raw>0))))
})
  

dt_out = data.table::data.table("tissue_celltype"=df.metadata.clust$tissue_celltype,
                    "NCells"=df.metadata.clust$n_cells,
                    "nCount_RNA_median"=mat_genestats[1,],
                    "nCount_gene_median"=mat_genestats[2,])


# ======================================================================= #
# ========================== WRITE OUT ============================= #
# ======================================================================= #

file.out = here("out","qc_checks", "tabula_muris_cell_type_confounderstats.csv")
data.table::fwrite(dt_out, file = file.out)


# ======================================================================= #
# ========================== XXXXXXXXXXXXXX ============================= #
# ======================================================================= #



# ======================================================================= #
# ========================== XXXXXXXXXXXXXX ============================= #
# ======================================================================= #



# ======================================================================= #
# ========================== XXXXXXXXXXXXXX ============================= #
# ======================================================================= #