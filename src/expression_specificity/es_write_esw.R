############### SYNOPSIS ###################
### AIM: Write ES 

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

setwd(here("src/expression_specificity"))


# ======================================================================= #
# ================================ PARAMS ============================ #
# ======================================================================= #

dataset_prefix <- "mousebrain_all"
path_out <- here("out/es")

# ======================================================================= #
# ================================ LOAD DATA ============================ #
# ======================================================================= #

load(here(sprintf("out/es/%s.es_obj.RData", dataset_prefix)))

# ======================================================================= #
# ================================ Write SEMs ================================= #
# ======================================================================= #


write_sems(sem_obj, slot="sem_transformed", dataset_prefix=dataset_prefix, dir_out=path_out) 
write_sems(sem_obj, slot="sem", dataset_prefix=dataset_prefix, dir_out=path_out) 
write_sems(sem_obj, slot="sem_pvalues", dataset_prefix=dataset_prefix, dir_out=path_out) 
write_sems(sem_obj, slot="null", dataset_prefix=dataset_prefix, dir_out=path_out) 
# sem_meta = mean, median, sd
# sem_transformed = ESw*
# sem = ESw
