############### SYNOPSIS ###################
### Cell-type BMI gene list enrichment 

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/geneset_enrichment"))

# ======================================================================= #
# ================================== INFO =============================== #
# ======================================================================= #


# Wilcoxon test analytical p-values.

### April 24th, run on FINAL/CORRECTED ESmu data in out/es/
# resultatfilerne er opdaterede (samme navn som foer)
# Campbell BMI: `/projects/jonatan/applied/17-holst-hsl/tables/campbell_GSA_SEMvsBMI_wilcoxon_analyt_1_mat_wilcoxonPval_genesetTests.txt`
# Mousebrain BMI: `/projects/jonatan/applied/18-mousebrain_7/tables/mb_ALL_GSA_SEMall_BMI_1.5_analyt_newgenelists_mat_wilcoxonPval_genesetTests.txt`
# Mousebrain WGCNA: `/projects/jonatan/applied/18-mousebrain_7/tables/mb_GSA_SEMall_WGCNAtop_2_analyt_mat_wilcoxonPval_genesetTests.txt`

# ======================================================================= #
# =============================== LOAD DATA ============================== #
# ======================================================================= #

### Data format: geneset x cell-types

### Mousebrain
file.data <- "/projects/jonatan/applied/18-mousebrain_7/tables/mb_ALL_GSA_SEMall_BMI_1.5_analyt_newgenelists_mat_wilcoxonPval_genesetTests.txt" # nearest_genes,associated_genes seperate.
# # df <- read.table(file.data, sep="\t", header=T) %>% rownames_to_column("geneset") %>% as_tibble()
df.mb <- read.table(file.data, sep="\t", header=T) %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("annotation") %>% 
  as_tibble()

# Campbell:
file.data <- "/projects/jonatan/applied/17-holst-hsl/tables/campbell_GSA_SEMvsBMI_wilcoxon_analyt_1_mat_wilcoxonPval_genesetTests.txt"
df.cb <- read.table(file.data, sep="\t", header=T) %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("annotation") %>% 
  as_tibble()

# ======================================================================= #
# =============================== MERGE DATA ============================== #
# ======================================================================= #
list.bind <- list("Mousebrain"=df.mb,
                  "Campbell"=df.cb)

df.join <- bind_rows(list.bind, .id="dataset")

# ======================================================================= #
# =============================== EXPORT DATA ============================== #
# ======================================================================= #

file.out <- here("results/es_enrichment--bmi_gene_lists.pvals.csv")
df.join %>% write_csv(file.out)

