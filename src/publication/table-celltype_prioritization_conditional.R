############### SYNOPSIS ###################
### AIM: Write conditional MB cell-type prioritization result tables

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
# ================================= RUN ================================== #
# ======================================================================= #



file.data <- here("results/prioritization_celltypes_conditional--mousebrain_all.BMI_UKBB_Loh2018.csv.gz")
df <- read_csv(file.data)

df.spread <- df %>% 
  select(annotation, condition, p.value) %>%
  spread(key="condition", value="p.value") %>%
  arrange(baseline)

# file.out <- here("src/publication/tables/", sprintf("table-celltype_ldsc_conditional_results.%s.csv", dataset_prefix))
file.out <- here("src/publication/tables/table-celltype_ldsc_conditional_results.mousebrain.csv")
df.spread %>% write_csv(file.out)




