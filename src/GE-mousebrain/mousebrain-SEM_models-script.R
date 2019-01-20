############### SYNOPSIS ###################
# SCRIPT to calculate SEMs on mousebrain

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library

# ======================================================================= #
# ================================== Mouse ============================== #
# ======================================================================= #

sem_obj_mouse <- create_sem_object("/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain")
sem_obj_mouse <- exclude_sporadic_expressed_genes(sem_obj_mouse)
# sem_obj_mouse <- calc_sem_wrapper(sem_obj_mouse)
# sem_obj_mouse <- bin_sems(sem_obj_mouse, n_bins=101, threshold_bin_zero=0)

print("Saving object...")
save(sem_obj_mouse,file="mousebrain.sem_obj_mouse.RData")

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
save(sem_obj, file="mousebrain.sem_obj.RData")

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
# ================================ XXXX ================================= #
# ======================================================================= #

library(tidyverse)
library(ComplexHeatmap)

# ======================= meta SEM - per anno ======================= #

for (anno_name in names(sem_obj[["group_by_annotation.sem_transformed"]])) {
  print(anno_name)
  df.plot <- sem_obj[["group_by_annotation.sem_transformed"]][[anno_name]] %>% 
    mutate(gene=sem_obj[["genes"]]) %>% 
    # slice(1:10) %>%
    column_to_rownames(var="gene") %>% as.data.frame() # this step must come last to keep the rownames.
  col_fun <- circlize::colorRamp2(c(0, 100), c("white", "red"))
  
  dummy <- tryCatch({ 
    message(sprintf("running ID=%s", anno_name))
    h1 <- Heatmap(df.plot, 
                  cluster_rows = T, 
                  cluster_columns = FALSE, 
                  show_row_names = F,
                  col = col_fun,
                  show_row_dend = F, # do not show row dendrogram | this might solve 'Error: node stack overflow'
                  use_raster = TRUE, 
                  raster_device = "png")
    file.out <- sprintf("heatmap_complex.sem_transformed.%s.pdf", anno_name)
    pdf(file.out, width = 8, height = 20) # {ncols=5}{w=8}
    time_elapsed <- system.time(draw(h1))
    print(time_elapsed)
    dev.off()
  }, warning = function(w) {
    message(sprintf("warning: %s", w))
  }, error = function(e) {
    message(sprintf("************* error *************: %s", e))
  }, finally = {
    message(sprintf("done with ID=%s", anno_name))
  }) # end tryCatch
}



# ======================= meta SEM ======================= #

for (sem_meta_name in names(sem_obj[["sem_meta"]])) {
  # sem_meta_name <- "mean"
  print(sem_meta_name)
  df.plot <- sem_obj[["sem_meta"]][[sem_meta_name]] %>% 
    mutate(gene=sem_obj[["genes"]]) %>% 
    # slice(1:10) %>%
    column_to_rownames(var="gene") %>% as.data.frame() # this step must come last to keep the rownames.
  col_fun <- circlize::colorRamp2(c(0, 100), c("white", "red"))
  h1 <- Heatmap(df.plot, 
                cluster_rows = T, 
                cluster_columns = FALSE, 
                show_row_names = F,
                col = col_fun,
                use_raster = TRUE, 
                raster_device = "png")
  file.out <- sprintf("heatmap_complex.sem_meta.%s.pdf", sem_meta_name)
  pdf(file.out, width = 40, height = 20)
  time_elapsed <- system.time(draw(h1))
  print(time_elapsed)
  dev.off()
}



