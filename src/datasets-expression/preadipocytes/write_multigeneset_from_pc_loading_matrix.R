############### SYNOPSIS ###################
# Write LDSC pipeline multi_geneset_file for adipocyte principal component loading matrix
# incl. normalizating matrix

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)
library(here)


source(here("src/lib/load_functions.R")) # load sc-genetics library


here("src/GE-adipocyte")


# ======================================================================= #
# ================================== XXXXX ============================== #
# ======================================================================= #

# Here is the path to the PCA loading matrix for PC's 1-10 (tab seperated values): 
file.data <- "/projects/timshel/sc-scheele_lab_adipose_fluidigm_c1/data-preadipocytes_developing/10x-180831-geneloadings-PC1-PC10.tab"

df <- read_tsv(file.data, skip=0, col_types = NULL, trim_ws=F)
df <- read_tsv(file.data, col_types=cols(gene="c",.default = col_guess()))

# ======================================================================= #
# ================================== XXXXX ============================== #
# ======================================================================= #

# ======================================================================= #
# ================================== XXXXX ============================== #
# ======================================================================= #


# ======================================================================= #
# ================== Write multi_geneset_file for LDSC ================= #
# ======================================================================= #

# dataset_prefix <- "campbell_lvl1"
# load(file=sprintf("%s.sem_obj.RData", dataset_prefix))

# df_multi_geneset <- write_multi_geneset_file(sem_obj, dataset_prefix=dataset_prefix)
# df_multi_geneset
