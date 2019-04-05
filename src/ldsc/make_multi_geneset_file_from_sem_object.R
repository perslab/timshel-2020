############### SYNOPSIS ###################
# Make "all genes" multi_geneset file from SEM objects

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

source(here("src/lib/load_functions.R")) # load sc-genetics library


setwd(here("src/ldsc/"))



# ======================================================================= #
# ========== MAKE multi_geneset_file mousebrain FDR sign. cell-type ======= #
# ======================================================================= #

load("/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain.sem_obj.RData")
# annotations <- c("MEGLU10","MEGLU1","DEINH3","TEGLU23","MEINH2","MEGLU11","DEGLU5","TEGLU17")
annotations <- c("MEGLU10","DEINH3","MEGLU1","MEINH2","TEGLU23","DEGLU5")
sem_obj.sub <- subset_annotations(sem_obj, annotations)
sem_obj.sub <- calc_sem_wrapper(sem_obj.sub)
sem_obj.sub <- calc_empirical_pvalues_wrapper(sem_obj.sub)
sem_obj.sub <- transform_sems(sem_obj.sub, method="rank_normalize")
sem_obj.sub <- set_group_by_annotation_slots(sem_obj.sub)
sem_obj.sub <- calc_sem_meta(sem_obj.sub)

### export
df_multi_geneset <- write_multi_geneset_file(sem_obj.sub, dataset_prefix="mousebrain_hierarchical_fdr_sign_only_190114", use_raw_sem_values=F)

# ======================================================================= #
# ================ MAKE multi_geneset_file for each dataset ============= #
# ======================================================================= #

list.datasets <- list("mousebrain_all"="/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain.sem_obj.RData",
                      "tabula_muris"="/projects/timshel/sc-genetics/sc-genetics/src/GE-maca/tabula_muris.sem_obj.RData",
                      "campbell_lvl1"="/projects/timshel/sc-genetics/sc-genetics/src/GE-hypothalamus/campbell_lvl1.sem_obj.RData",
                      "campbell_lvl2"="/projects/timshel/sc-genetics/sc-genetics/src/GE-hypothalamus/campbell_lvl2.sem_obj.RData")

for (dataset_prefix in names(list.datasets)) {
  print(dataset_prefix)
  file.sem_obj <- list.datasets[[dataset_prefix]]
  load(file.sem_obj)
  ### Write SEMs
  write_sems(sem_obj, slot="sem_meta", dataset_prefix)
  ### Write multi_geneset_file for LDSC
  # df_multi_geneset <- write_multi_geneset_file(sem_obj, dataset_prefix, use_raw_sem_values=F)
  df_multi_geneset <- write_multi_geneset_file(sem_obj, dataset_prefix=sprintf("%s_raw_sems", dataset_prefix), use_raw_sem_values=T) # *OBS*: only set use_raw_sem_values=T if you want raw (untransformed) SEM values exported.
}

### TMP FOR DEBUGGING 
# dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
# source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library
# load("/projects/timshel/sc-genetics/sc-genetics/src/GE-hypothalamus/campbell_lvl1.sem_obj.RData")
# df.multi_geneset <- write_multi_geneset_file(sem_obj, dataset_prefix, use_raw_sem_values=T)
# df.multi_geneset %>% filter(.sem_name =="tstat_top10pct_binary") %>% count(.annotation)
# df.multi_geneset %>% filter(.sem_name =="specificity_top10pct_binary") %>% count(.annotation)

# ======================================================================= #
# =============================== MAKE ALL GENES multi_geneset_file ================================= #
# ======================================================================= #

list.datasets <- list("mousebrain"="/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain.sem_obj.RData",
                      "tabula_muris"="/projects/timshel/sc-genetics/sc-genetics/src/GE-maca/tabula_muris.sem_obj.RData",
                      "campbell"="/projects/timshel/sc-genetics/sc-genetics/src/GE-hypothalamus/campbell_lvl1.sem_obj.RData") # OBS: we treat campbell_lvl1 and campbell_lvl2 identical.

### Get genes for each dataset
list.dfs <- list()
for (dataset_prefix in names(list.datasets)) {
  print(dataset_prefix)
  file.sem_obj <- list.datasets[[dataset_prefix]]
  load(file.sem_obj)
  df.tmp <- tibble(geneset_name=sprintf("all_genes_in_dataset.%s",dataset_prefix), gene=sem_obj[["genes"]], value=1) # value=1 is how we encode all genes
  print(nrow(df.tmp))
  list.dfs[[dataset_prefix]] <- df.tmp
}

### Combine
df_multi_geneset <- bind_rows(list.dfs)
df_multi_geneset

### Write file
file.out <- "multi_geneset.all_genes_in_dataset.txt"
df_multi_geneset %>% write_tsv(file.out, col_names=F) # no header!

# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #



# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #








