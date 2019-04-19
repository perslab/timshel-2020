############### SYNOPSIS ###################
### AIM: XXX top ESmu genes for selected annotations

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
# ================================ LOAD DATA ============================ #
# ======================================================================= #

dataset_prefix <- "campbell_lvl2"

### Get ES data
load(here(sprintf("out/es/%s.es_obj.RData", dataset_prefix)))

# []Check Campbell [n29.Nr5a1/Bdnf cell-type] for Cckar (Cck1), Cckbr (Cck2) --> Schwartz thinks this is an important cell-type.

# ======================================================================= #
# ===================== SOM table (csv): all genes ES =================== #
# ======================================================================= #

# filter.annotations <- sem_obj[["annotations"]]

### ES table for all genes (SOM)
df.all_table <- get_annotation_es.table(sem_obj, annotations=filter.annotations, es_metric="es_mu")
file.out <- sprintf("tables/table-es_mu.%s.csv.gz", dataset_prefix)
df.all_table %>% write_csv(file.out)

