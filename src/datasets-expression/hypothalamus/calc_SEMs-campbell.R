############### SYNOPSIS ###################
# Calculate Campbell SEM models

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library


wd <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-hypothalamus/"
setwd(wd)

# ======================================================================= #
# ================================== Mouse ============================== #
# ======================================================================= #

# dataset_prefix <- "campbell_lvl1"
dataset_prefix <- "campbell_lvl2"

sem_obj_mouse <- create_sem_object(sprintf("/projects/timshel/sc-genetics/sc-genetics/data/expression/hypothalamus_campbell/%s", dataset_prefix))
sem_obj_mouse <- exclude_sporadic_expressed_genes(sem_obj_mouse)
# sem_obj_mouse <- calc_sem_wrapper(sem_obj_mouse)
# sem_obj_mouse <- bin_sems(sem_obj_mouse, n_bins=101, threshold_bin_zero=0)

### LVL1
# [1] "Number of genes with NA in anova pvalue: 410. (This number is not used for anything. We just report it for your information)"
# [1] "Number of sporadic_expressed_genes: 10055"
# [1] "Number of zero-mean genes: 410"
# [1] "Total number of genes marked for exclusion: 10465"
# [1] "Number of genes in object: 16309"

### LVL2
# [1] "Number of genes with NA in anova pvalue: 410. (This number is not used for anything. We just report it for your information)"
# [1] "Number of sporadic_expressed_genes: 10129"
# [1] "Number of zero-mean genes: 410"
# [1] "Total number of genes marked for exclusion: 10539"
# [1] "Number of genes in object: 16235"

# print("Saving object...")
# save(sem_obj_mouse, file=sprintf("%s.sem_obj_mouse.RData", dataset_prefix))

# # ======================================================================= #
# # ============================== Human =============================== #
# # ======================================================================= #
# 
sem_obj <- map_to_human(sem_obj_mouse, type_mouse_gene_ids="mgi")
# [1] "Number of genes NOT mapped from MGI to Ensembl ID: 655 out of 16309 genes (4.02 pct)"
# [1] "Number of genes NOT mapped to human ortholog: 3775 out of 15654 genes (24.12 pct)"
# [1] "Number of mapped genes in output: 12534 out of 16309 input genes (total mapping rate = 76.85 pct)"
# [1] "Total number of genes marked for exclusion: 3775"
# [1] "Number of genes in object: 12534"
sem_obj <- calc_sem_wrapper(sem_obj) # NEW
sem_obj <- calc_empirical_pvalues_wrapper(sem_obj) # NEW
sem_obj <- transform_sems(sem_obj, method="rank_normalize")
sem_obj <- set_group_by_annotation_slots(sem_obj)
sem_obj <- calc_sem_meta(sem_obj)

print("Saving object...")
save(sem_obj, file=sprintf("%s.sem_obj.RData", dataset_prefix))

# ======================================================================= #
# ================== Write multi_geneset_file for LDSC ================= #
# ======================================================================= #

# dataset_prefix <- "campbell_lvl1"
# load(file=sprintf("%s.sem_obj.RData", dataset_prefix))

# df_multi_geneset <- write_multi_geneset_file(sem_obj, dataset_prefix=dataset_prefix)
# df_multi_geneset


# ======================================================================= #
# ============================== LOAD =============================== #
# ======================================================================= #


# ======================================================================= #
# ======================== Hierarchical analysis ======================== #
# ======================================================================= #
