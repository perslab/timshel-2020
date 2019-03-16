############### SYNOPSIS ###################
### AIM: calculate ES. 

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

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library

setwd(here("src/GE-mousebrain"))

# ======================================================================= #
# ================================== Mouse ============================== #
# ======================================================================= #

sem_obj_mouse <- create_sem_object("/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain")
sem_obj_mouse <- exclude_sporadic_expressed_genes(sem_obj_mouse)
# sem_obj_mouse <- calc_sem_wrapper(sem_obj_mouse)
# sem_obj_mouse <- bin_sems(sem_obj_mouse, n_bins=101, threshold_bin_zero=0)

print("Saving object...")
save(sem_obj_mouse,file="mousebrain-190306.sem_obj_mouse.RData")

# # ======================================================================= #
# # ============================== Human =============================== #
# # ======================================================================= #
# 
sem_obj <- map_to_human(sem_obj_mouse, type_mouse_gene_ids="ensembl")
sem_obj <- calc_sem_wrapper(sem_obj) # NEW
sem_obj <- calc_empirical_pvalues_wrapper(sem_obj) # NEW
sem_obj <- transform_sems(sem_obj, method="rank_normalize")
sem_obj <- set_group_by_annotation_slots(sem_obj)
sem_obj <- calc_sem_meta(sem_obj)

print("Saving object...")
save(sem_obj, file="mousebrain-190306.sem_obj.RData")
print("Done saving")

# ======================================================================= #
# ============================== TMP STOP POINT =============================== #
# ======================================================================= #

stop("QUITING SCRIPT - DONE")


# ======================================================================= #
# ============================== MAKE ES HEATMAPS ========================== #
# ======================================================================= #

### Load function
# source("lib-es_heatmap.R")

### Run functions
# make_es_heatmap.per_anno.genesXes_metrics(sem_obj)
# make_es_heatmap.es_meta.genesXannotations(sem_obj)
# print("Script done")


# ======================================================================= #
# ============================== LOAD =============================== #
# ======================================================================= #

# load(file="mousebrain.sem_obj.RData")

# ======================================================================= #
# ======================== Hierarchical analysis ======================== #
# ======================================================================= #

### Load metadata
cols_keep.metadata <- c("ClusterName","Class","NCells","Description","Probable_location","Region","Developmental_compartment","TaxonomyRank1","TaxonomyRank2","TaxonomyRank3","TaxonomyRank4")
file.metadata <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain-agg_L5.metadata.csv"
df.metadata <- read_csv(file.metadata) %>% select(cols_keep.metadata)
df.metadata

# GROUP BY: Class
annotations.Class <- with(df.metadata, split(ClusterName, Class))
names(annotations.Class)
sem_obj.hier.Class <- hierarchical_sem(object=sem_obj, list.annotations=annotations.Class)
save(sem_obj.hier.Class, file="mousebrain.sem_obj.hier.Class.RData")

# GROUP BY: TaxonomyRank2
annotations.TaxonomyRank2 <- with(df.metadata, split(ClusterName, TaxonomyRank2))
names(annotations.TaxonomyRank2)
sem_obj.hier.TaxonomyRank2 <- hierarchical_sem(object=sem_obj, list.annotations=annotations.TaxonomyRank2)
save(sem_obj.hier.TaxonomyRank2, file="mousebrain.sem_obj.hier.TaxonomyRank2.RData")


# ======================================================================= #
# ============================= Hierarchical: SUBSET ALL NEURONS ================================= #
# ======================================================================= #

file.metadata <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain-agg_L5.metadata.csv"
df.metadata <- read_csv(file.metadata) %>% select(ClusterName, Class)
annotations.neurons <- df.metadata %>% filter(Class=="Neurons") %>% pull(ClusterName)
length(annotations.neurons) # 214

### Subset
sem_obj.sub <- subset_annotations(sem_obj, annotations=annotations.neurons)
sem_obj.sub <- calc_sem_wrapper(sem_obj.sub) # NEW
sem_obj.sub <- calc_empirical_pvalues_wrapper(sem_obj.sub) # NEW
sem_obj.sub <- transform_sems(sem_obj.sub, method="rank_normalize")
sem_obj.sub <- set_group_by_annotation_slots(sem_obj.sub)
sem_obj.sub <- calc_sem_meta(sem_obj.sub)
# save(sem_obj.sub, file="mousebrain_neurons_190306.sem_obj.RData")

### Write file
df_multi_geneset <- write_multi_geneset_file(sem_obj.sub, 
                                             dataset_prefix="mousebrain_neurons",
                                             es_mean_only=T)
df_multi_geneset

# ======================================================================= #
# ============================= Hierarchical: SUBSET FDR CELL-TYPES  ================================= #
# ======================================================================= #

### Subset
annotations.fdr_sign <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12") # BMI_UKBB_Loh2018 FDR sign cell-types mousebrain
sem_obj.sub <- subset_annotations(sem_obj, annotations=annotations.fdr_sign)
sem_obj.sub <- calc_sem_wrapper(sem_obj.sub) # NEW
sem_obj.sub <- calc_empirical_pvalues_wrapper(sem_obj.sub) # NEW
sem_obj.sub <- transform_sems(sem_obj.sub, method="rank_normalize")
sem_obj.sub <- set_group_by_annotation_slots(sem_obj.sub)
sem_obj.sub <- calc_sem_meta(sem_obj.sub)

### Write file
df_multi_geneset <- write_multi_geneset_file(sem_obj.sub, 
                                             dataset_prefix="mousebrain_bmi_loh2018_11fdr_celltypes",
                                             es_mean_only=T)
df_multi_geneset


