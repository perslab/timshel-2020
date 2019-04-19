############### SYNOPSIS ###################
### AIM: Make hypothalamus plot

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
# ================================ MOUSEBRAIN ================================ #
# ======================================================================= #

dataset_prefix <- "mousebrain_all"
filter.gwas <- "BMI_UKBB_Loh2018"

# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #

### Read LDSC results
file.results <- here("results/prioritization_celltypes--mousebrain_all.multi_gwas.csv.gz")
df.ldsc_cts <- read_csv(file.results)

# =========================== FILTER GWAS =========================== #
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == filter.gwas)

# =========================== ADD METADATA =========================== #

df.metadata <- get_metadata(dataset_prefix)
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation")

# ======================== FILTER REGION =========================== #

filter.annotations <- get_annotations.mousebrain.hypothalamus()
df.ldsc_cts <- df.ldsc_cts %>% filter(annotation %in% filter.annotations)

# ======================================================================= #
# ============================ CAMPBELL LVL2 ============================ #
# ======================================================================= #

dataset_prefix <- "campbell_lvl2"
filter.gwas <- "BMI_UKBB_Loh2018"

# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #

### Read LDSC results
file.results <- here("results/prioritization_celltypes--campbell_lvl2.multi_gwas.csv.gz")
df.ldsc_cts <- read_csv(file.results)

# =========================== FILTER GWAS =========================== #
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == filter.gwas)

# =========================== ADD METADATA =========================== #
df.metadata <- get_metadata(dataset_prefix)
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation")
df.ldsc_cts
# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #


# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #


# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #