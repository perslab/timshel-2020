############### SYNOPSIS ###################
# Calculate Tabula Muris SEM models
# USAGE: time Rscript tabula_muris-SEM_models.R |& tee tabula_muris-SEM_models.181211.out.R

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


# wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-hypothalamus/"
# setwd(wd)

# ======================================================================= #
# ================================== Mouse ============================== #
# ======================================================================= #

dataset_prefix <- "tabula_muris"

sem_obj_mouse <- create_sem_object(sprintf("/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/tabula_muris/%s", dataset_prefix))
sem_obj_mouse <- exclude_sporadic_expressed_genes(sem_obj_mouse)
# [1] "Number of genes with NA in anova pvalue: 390. (This number is not used for anything. We just report it for your information)"
# [1] "Number of sporadic_expressed_genes: 2365"
# [1] "Number of zero-mean genes: 390"
# [1] "Total number of genes marked for exclusion: 2755"
# [1] "Number of genes in object: 20586"

# sem_obj_mouse <- calc_sem_wrapper(sem_obj_mouse)
# sem_obj_mouse <- bin_sems(sem_obj_mouse, n_bins=101, threshold_bin_zero=0)

# print("Saving object...")
# save(sem_obj_mouse, file=sprintf("%s.sem_obj_mouse.RData", dataset_prefix))

# # ======================================================================= #
# # ============================== Human =============================== #
# # ======================================================================= #
# 
sem_obj <- map_to_human(sem_obj_mouse, type_mouse_gene_ids="mgi")
# [1] "Number of genes NOT mapped from MGI to Ensembl ID: 2063 out of 20586 genes (10.02 pct)"
# [1] "Number of genes NOT mapped to human ortholog: 6283 out of 18523 genes (33.92 pct)"
# [1] "Number of mapped genes in output: 14303 out of 20586 input genes (total mapping rate = 69.48 pct)"
# [1] "Total number of genes marked for exclusion: 6283"
# [1] "Number of genes in object: 14303"
sem_obj <- calc_sem_wrapper(sem_obj)
sem_obj <- calc_empirical_pvalues_wrapper(sem_obj)
sem_obj <- transform_sems(sem_obj, method="rank_normalize")
sem_obj <- set_group_by_annotation_slots(sem_obj)
sem_obj <- calc_sem_meta(sem_obj)

print("Saving object...")
save(sem_obj, file=sprintf("%s.sem_obj.RData", dataset_prefix))

# ======================================================================= #
# ============================== LOAD =============================== #
# ======================================================================= #

# ======================================================================= #
# ======================== Hierarchical analysis ======================== #
# ======================================================================= #
