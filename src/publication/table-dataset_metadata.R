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


dataset_prefixes <- c("tabula_muris", "mousebrain", "hypothalamus")

for (dataset_prefix in dataset_prefixes) {
  df.metadata <- get_metadata(dataset_prefix)
  if (dataset_prefix == "tabula_muris") {
    df.metadata <- df.metadata %>% mutate(annotation_fmt=utils.rename_annotations.tabula_muris(annotation, style="tissue - celltype", check_all_matches=T))
    df.metadata <- df.metadata %>% mutate(tissue = case_when(
      tissue == "Brain_Myeloid" ~ "Brain",
      tissue == "Brain_Non-Myeloid" ~ "Brain",
      TRUE ~ as.character(tissue))
    )
    df.metadata <- df.metadata %>% mutate(tissue = stringr::str_replace_all(tissue, pattern="_", replacement=" "))
    df.metadata <- df.metadata %>% mutate(cell_type = stringr::str_replace_all(cell_type, pattern="_", replacement=" "))
  }
  if (dataset_prefix == "hypothalamus") {
    specificity_ids_tmp <- get_scrna_seq_dataset_prefixes("hypo")
    file.n_es <- here("src/publication/tables", paste0("table-annotation_summary.", specificity_ids_tmp, ".csv"))
    list.df.n_es <- file.n_es %>% map(read_csv)
    names(list.df.n_es) <- specificity_ids_tmp
    df.n_es <- bind_rows(list.df.n_es, .id="specificity_id")
    df.n_es <- df.n_es %>% mutate(annotation_fmt = utils.rename_annotations.hypothalamus(annotation, specificity_id, check_all_matches=F)) # check_all_matches False because MB anno included
    df.metadata <- df.metadata %>% left_join(df.n_es %>% select(-annotation, -specificity_id), by="annotation_fmt")
    df.metadata <- df.metadata %>% arrange(annotation_prefix, annotation_fmt) # sort | *WITH annotation_prefix*
  } else {
    file.n_es <- here("src/publication/tables", sprintf("table-annotation_summary.%s.csv", dataset_prefix))
    df.n_es <- read_csv(file.n_es)
    df.metadata <- df.metadata %>% left_join(df.n_es, by="annotation")
    df.metadata <- df.metadata %>% arrange(annotation) # sort
  }
  
  file.out <- here("src/publication/tables", sprintf("table-metadata_combined.%s.csv", dataset_prefix))
  df.metadata %>% write_csv(file.out)
}



