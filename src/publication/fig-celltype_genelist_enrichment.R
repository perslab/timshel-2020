

############### SYNOPSIS ###################
### Module cell-type enrichment plot

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

library(RColorBrewer)

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/publication"))


# ======================================================================= #
# =============================== LOAD DATA ============================== #
# ======================================================================= #

### JT April 8th 2019
### all mousebrain celltypes ES scores vs BMI genesets (including two new ones). Wilcoxon test analytical p-values.
file.data <- "/projects/jonatan/applied/18-mousebrain_7/tables/mb_ALL_GSA_SEMall_BMI_1.5_analyt_newgenelists_mat_wilcoxonPval_genesetTests.txt" # nearest_genes,associated_genes seperate.
df <- read.table(file.data, sep="\t", header=T) %>% rownames_to_column("module_id") %>% as_tibble
df



