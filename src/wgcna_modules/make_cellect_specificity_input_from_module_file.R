############### SYNOPSIS ###################
# Generate specificity input for CELLECT from WGCNA module output table


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================ #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library

# ======================================================================= #
# =============================== CONSTANTS ================================ #
# ======================================================================= #

wd <- here("src/wgcna_modules")
setwd(wd)

# ======================================================================= #
# =========================== FILE ORIGIN INFO ========================== #
# ======================================================================= #

### Modules from FDR significant cell-types [v3, 190213] | NEW: WGCNA run on n=11 BMI_UKBB_Loh2018 FDR cell-types 
# Feb 13 [6:15 PM]
# Hej Pascal,
# Nu er den nye mb WGCNA analyse  endelig klar. Jeg gik tilbage til den tidligere udgave af pipelinen fra Januar.
# Koerslen hedder Neurons_sub_ClusterName_7.2_run1 og den relevante output fil er
# /projects/jonatan/mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_cell_cluster_module_genes.csv.gz
# Dodgerblue modulet gaar igen uden aendring, nu under navnet “lavenderblush”.
# (har i senere aendringer i scriptet soerget for at alle random seeds, inklusive navne, er reproducible) 
# dict_genomic_annot = {"wgcna.mousebrain-190213.fdr_sign_celltypes.continuous": 
# 					  	{"dataset":"mousebrain", 
# 					  	"file_multi_gene_set":"/projects/jonatan/mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_cell_cluster_module_genes.csv.gz"}
# 					 }
# New file path [August 2019] ----> /projects/jonatan/pub-perslab/18-mousebrain/18-mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_cell_cluster_module_genes.csv.gz

### Shell commands to generate out/wgcna/modules.mousebrain_bmicelltypes.rwgcna_table.csv
### NB: we filter out TEGLU4 (not significant after fix)
# FILE_M=/projects/jonatan/pub-perslab/18-mousebrain/18-mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_cell_cluster_module_genes.csv.gz
# zcat $FILE_M | grep -v TEGLU4 | gzip > /projects/timshel/sc-genetics/timshel-bmicelltypes2019/out/wgcna/modules.mousebrain_bmicelltypes.rwgcna_table.csv.gz

# ======================================================================= #
# =============================== WRANGE DATA ========================== #
# ======================================================================= #

file.table <- here("out/wgcna/modules.mousebrain_bmicelltypes.rwgcna_table.csv.gz")
df.long <- read_csv(file.table)

## long --> matrix
df <- df.long %>% select(module, ensembl, pkME) %>% spread(key='module', value='pkME')
df[is.na(df)] <- 0 # ALT:  df %>% mutate_all(~replace(., is.na(.), 0)) 

### map mouse to human orthologs
df <- df %>% as.data.frame() %>% column_to_rownames("ensembl")
df.human <- mouse_to_human_ortholog_gene_expression_mapping(df, type_mouse_gene_ids="ensembl")
# [1] "Reading gene mapping files..."
# [1] "Subsetting on orthologs..."
# [1] "Number of genes not mapped to human ortholog: 629 out of 7744 genes (8.12 pct)"
df.human <- df.human %>% rownames_to_column(var="ensembl") %>% as_tibble()


### Export
file.out <- here("out/wgcna/modules.mousebrain_bmicelltypes.specificity_human.csv.gz")
df.human %>% write_csv(file.out)



### Make module 'metadata' file
# file.out.meta <- here("out/wgcna/modules.mousebrain_bmicelltypes.metadata.csv.gz") 
# use ---> results/modules--metadata.txt instead


