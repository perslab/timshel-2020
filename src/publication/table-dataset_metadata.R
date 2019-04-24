############### SYNOPSIS ###################
### AIM: Make meta-data tables

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
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ============================== LOAD METADATA ============================= #
# ======================================================================= #

dataset_prefixes <- c("tabula_muris", "mousebrain_all", "campbell_lvl2")

for (dataset_prefix in dataset_prefixes) {
  df.metadata <- get_metadata(dataset_prefix)
  file.n_es <- here("src/publication/tables", sprintf("table-n_es_genes.%s.csv", dataset_prefix))
  df.n_es <- read_csv(file.n_es)
  df.metadata <- df.metadata %>% left_join(df.n_es, by="annotation")
  df.metadata <- df.metadata %>% arrange(annotation) # sort
  file.out <- here("src/publication/tables", sprintf("table-metadata_combined.%s.csv", dataset_prefix))
  df.metadata %>% write_csv(file.out)
}



