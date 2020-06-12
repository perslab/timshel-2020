############### SYNOPSIS ###################
# AIM: show that the prioritization results are not confounded by scRNA-seq associated factors
# N-cells
# N-genes
# ...

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)

library(corr)

# library(corrplot) # https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ========================== READ META DATA ============================= #
# ======================================================================= #

cols_metadata_keep <- c("ClusterName",
                        "Class",
                        "Description",
                        "NCells",
                        "Comment"
)

### Read
file.metadata <- here("data/expression/mousebrain/mousebrain_agg_L5.metadata_full.csv") # exported via loompy from "L5_All.agg.loom"
df.metadata <- read_csv(file.metadata)

### Select metadata
df.metadata <- df.metadata %>% 
  select(cols_metadata_keep) %>% 
  rename(annotation = ClusterName)


# ======================================================================= #
# ========================== XXXXXXXXXXXXXX ============================= #
# ======================================================================= #





# ======================================================================= #
# ========================== XXXXXXXXXXXXXX ============================= #
# ======================================================================= #



# ======================================================================= #
# ========================== XXXXXXXXXXXXXX ============================= #
# ======================================================================= #



# ======================================================================= #
# ========================== XXXXXXXXXXXXXX ============================= #
# ======================================================================= #



# ======================================================================= #
# ========================== XXXXXXXXXXXXXX ============================= #
# ======================================================================= #

