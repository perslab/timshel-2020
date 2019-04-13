############### SYNOPSIS ###################
# Calculate ES

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== USAGE ================================ #
# ======================================================================= #

# SEE run_es_calculation.sh

# ======================================================================= #
# =============================== SETUP ================================ #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library

# ======================================================================= #
# =============================== CONSTANTS ================================ #
# ======================================================================= #

# TODO: test that outdir writeable
path_out <- here("out/es")


# ======================================================================= #
# ================================== FUNCTION ============================== #
# ======================================================================= #

wrapper_calc_es <- function(path_prefix.es_precalc, path_out, dataset_prefix, version_stamp, type_mouse_gene_ids) {
  
  # ================================== Mouse ============================== #
  sem_obj_mouse <- create_sem_object(path_prefix.es_precalc)
  sem_obj_mouse <- exclude_sporadic_expressed_genes(sem_obj_mouse)
  # sem_obj_mouse <- calc_sem_wrapper(sem_obj_mouse)
  # sem_obj_mouse <- bin_sems(sem_obj_mouse, n_bins=101, threshold_bin_zero=0)
  
  # print("Saving object...")
  # file.out <- sprintf("%s/%s.es_obj_mouse-%s.RData", path_out, dataset_prefix, version_stamp)
  # save(sem_obj_mouse, file=file.out)
  
  # ============================== Human =============================== #
  
  sem_obj <- map_to_human(sem_obj_mouse, type_mouse_gene_ids=type_mouse_gene_ids)
  sem_obj <- calc_sem_wrapper(sem_obj)
  sem_obj <- calc_empirical_pvalues_wrapper(sem_obj)
  sem_obj <- transform_sems(sem_obj, method="rank_normalize")
  sem_obj <- set_group_by_annotation_slots(sem_obj)
  sem_obj <- calc_sem_meta(sem_obj)
  
  print("Saving object...")
  file.out <- sprintf("%s/%s.es_obj-%s.RData", path_out, dataset_prefix, version_stamp)
  save(sem_obj, file=file.out)
  print("Done saving")
  
  # ================================ Write SEMs ================================= #
  
  write_sems(sem_obj, slot="sem_meta", dataset_prefix=dataset_prefix, dir_out=path_out) # sem_meta writes out mean, median, sd
  print("Done with function")
}


# ======================================================================= #
# =============================== ARGS ================================ #
# ======================================================================= #

### CMD Args
dataset_prefix <- commandArgs(trailingOnly=TRUE)[1]
version_stamp <- commandArgs(trailingOnly=TRUE)[2]
print(sprintf("RUNNING dataset_prefix=%s, version_stamp=%s", dataset_prefix, version_stamp))

### DEV | MANUAL Args
# version_stamp <- "999999"
# dataset_prefix2run <- c("tabula_muris","mousebrain_all","campbell_lvl1","campbell_lvl2")
# dataset_prefix <- "campbell_lvl1"
# dataset_prefix <- "campbell_lvl2"
# dataset_prefix <- "tabula_muris"
# dataset_prefix <- "mousebrain_all"


if (dataset_prefix == "mousebrain_all") {
  path_prefix.es_precalc <- here("data/expression/mousebrain/mousebrain")
  type_mouse_gene_ids <- "ensembl"
} else if (dataset_prefix == "tabula_muris") {
  path_prefix.es_precalc <- here("data/expression/tabula_muris/tabula_muris") 
  type_mouse_gene_ids <- "mgi"
} else if (dataset_prefix %in% c("campbell_lvl1", "campbell_lvl2")) {
  path_prefix.es_precalc <- here("data/expression/hypothalamus_campbell/", dataset_prefix) 
  type_mouse_gene_ids <- "mgi"
} else {
  stop(sprintf("Got wrong dataset_prefix %s", dataset_prefix))
}

# ======================================================================= #
# =============================== TMP DEV ================================ #
# ======================================================================= #

# source(here("src/lib/load_functions.R")) # load sc-genetics library
# sem_obj_mouse <- create_sem_object(path_prefix.es_precalc)
# sem_obj_mouse[["annotations"]]
# sem_obj_mouse[["data"]][["mean"]]
# stringr::str_detect(sem_obj_mouse[["annotations"]], pattern="/")
# stringr::str_replace_all(sem_obj_mouse[["annotations"]], c("/"="-", "\\s+"="_"))

# ======================================================================= #
# =============================== RUN ================================ #
# ======================================================================= #

wrapper_calc_es(path_prefix.es_precalc, path_out, dataset_prefix, version_stamp, type_mouse_gene_ids)

print("Script done!")

# ======================================================================= #
# ============================== MAKE ES HEATMAPS ========================== #
# ======================================================================= #

### Run functions
# make_es_heatmap.per_anno.genesXes_metrics(sem_obj, dataset_prefix)
# make_es_heatmap.es_meta.genesXannotations(sem_obj, dataset_prefix)
# print("Script done")


# ======================================================================= #
# ======================== Hierarchical analysis ======================== #
# ======================================================================= #

# SEE MOUSEBRAIN SCRIPT

