############### SYNOPSIS ###################
### AIM: write_multi_geneset_file for LDSC

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
# ================================ SAVE/LOAD ================================= #
# ======================================================================= #

### load data
# load("/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain.sem_obj.RData")
load("/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain-190306.sem_obj.RData")

# ======================================================================= #
# ================================ Write SEMs ================================= #
# ======================================================================= #

# write_sems(sem_obj, slot="sem_meta", name.dataset=dataset_prefix)

# ======================================================================= #
# ================== Write multi_geneset_file for LDSC ================= #
# ======================================================================= #


df_multi_geneset <- write_multi_geneset_file(sem_obj, dataset_prefix="mousebrain_all_190306_es_fix",
                                             es_mean_only=T)
df_multi_geneset



