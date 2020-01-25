############### SYNOPSIS ###################
### AIM: Export MAGMA results

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

library(GGally)
library(corrr) # devtools::install_github("drsimonj/corrr")
library(gghighlight)


source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ================================ MAGMA ================================ #
# ======================================================================= #
file.magma <- here("results/cellect_magma/out.CELLECT_MAGMA.mousebrain.BMI_UKBB_Loh2018_no_mhc.es_mu.csv")
df.magma.raw <- read_csv(file.magma)
df.magma <- df.magma.raw %>% select(annotation, p.value, estimate, std.error, `FDR < 0.05`=fdr_significant)
df.magma <- df.magma %>% arrange(p.value)

file.out <- sprintf("tables/table_celltype_priori_magma.csv")
df.magma %>% write_csv(file.out)
