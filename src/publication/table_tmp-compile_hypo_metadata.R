############### SYNOPSIS ###################
# Compile hypothalamus meta-data from all data sets
# We will use this list of compiled annotations as the basis for manually cleaning the annotation names
# SUMMMARY: THIS WILL ONLY OUTPUT A TEMPORARY FILE THAT IS USED FOR LATER MANUAL CLEANING.

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
# ========================= GET ALL HYPO ANNOTATIONS ==================== #
# ======================================================================= #


### Read LDSC results
file.results <- here("results/cellect_ldsc/prioritization.csv")
df.ldsc_cts <- read_csv(file.results)
df.ldsc_cts <- format_cellect_ldsc_results(df.ldsc_cts)

### Get unique hypo annotations
data_prefixes <- get_scrna_seq_dataset_prefixes("hypo")
df.anno <- df.ldsc_cts %>% 
  filter(specificity_id %in% data_prefixes) %>% 
  distinct(specificity_id, annotation)
df.anno
# df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == "BMI_UKBB_Loh2018")

# ======================================================================= #
# =======================  ROUGH CLEAN / REPLACE CHARS ================== #
# ======================================================================= #

### Init
df.anno <- df.anno %>% mutate(annotation_fmt=annotation)

### Clean Campbel annotation names (some meta-data for Campbell was not correctly provided by the authors)
# recode: x=new_name; name(x)=old_name
vector_rename <- c("n03"="n03.Th/Sst",
                   "n07"="n07.Arx/Nr5a2",
                   "n08"="n08.Th/Slc6a3",
                   "n27"="n27.Tbx19",
                   "n29.Nr5a1-Adcyap1"="n29.Nr5a1/Bdnf")
df.anno <- df.anno %>% mutate(annotation_fmt=recode(annotation_fmt, !!!vector_rename))


### Relace
df.anno <- df.anno %>% mutate(annotation_fmt = str_replace_all(annotation_fmt, "\\.", " "))
df.anno <- df.anno %>% mutate(annotation_fmt = str_replace_all(annotation_fmt, "_", " "))
df.anno <- df.anno %>% mutate(annotation_fmt = str_replace_all(annotation_fmt, "-", "/"))
df.anno


### Sort 
df.anno <- df.anno %>% arrange(specificity_id, annotation_fmt)

### Prefix
# -Campbell2016=ARC-ME
# -kim2019=VMH
# -chen2017=HYP-Chen
# -mikkelsen2019=LHA
# -moffitt2018=POA
# -romanov2017=HYP-Romanov
df.anno <- df.anno %>% mutate(annotation_prefix = case_when(
  specificity_id=="campbell2017_lvl2" ~ "ARCME",
  specificity_id=="romanov2017" ~ "HYPR",
  specificity_id=="chen2017" ~ "HYPC",
  specificity_id=="kimVMH2019_smartseq" ~ "VMH",
  specificity_id=="mikkelsen2019" ~ "LHA",
  specificity_id=="moffitt2018" ~ "POA",
  TRUE ~ "ERROR IN PREFIX MAPPING"
))


### ADD TAXONOMY
df.metadata.campbell <- get_metadata("campbell2017_lvl2")
df.anno <- df.anno %>% left_join(df.metadata.campbell, by="annotation")

### Final: add prefix
df.anno <- df.anno %>% mutate(annotation_fmt_prefix = paste0(annotation_prefix, "-", annotation_fmt))


### Export
df.anno %>% write_csv(here("data/expression/hypothalamus/hypothalamus_metadata.unformated.csv"))


