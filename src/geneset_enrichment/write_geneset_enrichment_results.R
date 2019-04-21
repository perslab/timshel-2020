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
# =============================== LOAD DATA ============================== #
# ======================================================================= #

### Data format: geneset x cell-types

### Mousebrain
### JT April 8th 2019: all mousebrain celltypes ES scores vs BMI genesets (including two new ones). Wilcoxon test analytical p-values.
file.data <- "/projects/jonatan/applied/18-mousebrain_7/tables/mb_ALL_GSA_SEMall_BMI_1.5_analyt_newgenelists_mat_wilcoxonPval_genesetTests.txt" # nearest_genes,associated_genes seperate.
# # df <- read.table(file.data, sep="\t", header=T) %>% rownames_to_column("geneset") %>% as_tibble()
df.mb <- read.table(file.data, sep="\t", header=T) %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("annotation") %>% 
  as_tibble()

# Campbell:
file.data <- "/projects/jonatan/applied/17-holst-hsl/tables/cambpbell_GSA_SEMvsBMI_wilcoxon_analyt_1_mat_wilcoxonPval_genesetTests.txt"
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

