############### SYNOPSIS ###################
# Make .multi_geneset file from WGCNA output file for mousebrain all annotations combined

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/ldsc/"
setwd(wd)


library(tidyverse)

file.cluster_module_genes <- "/projects/jonatan/tmp-mousebrain/tables/mb_ClusterName_5_cell_cluster_module_genes.csv.gz" # Jon file, Dec 5th 2018
df.cluster_module_genes <- read_csv(file.cluster_module_genes)
df.cluster_module_genes <- df.cluster_module_genes %>% mutate(module_id = paste0(cell_cluster, ".", module))
n_distinct(df.cluster_module_genes$module) # [1] 14269
n_distinct(df.cluster_module_genes$module_id) # [1] 17108
n_distinct(df.cluster_module_genes$cell_cluster) # 160
n_distinct(df.cluster_module_genes$run) # 7
n_distinct(df.cluster_module_genes$ensembl) # 15206
