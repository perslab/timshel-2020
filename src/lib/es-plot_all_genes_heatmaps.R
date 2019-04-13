############### SYNOPSIS ###################
### AIM: Make ES heatmaps for ALL genes
### OBS: runtime and be many hours because of Clustering and saving of large .pdfs


### ComplexHeatmap
### Install
# source("https://bioconductor.org/biocLite.R")
# biocLite("ComplexHeatmap")
### ComplexHeatmap REFs:
# http://www.bioconductor.org/packages/3.6/bioc/html/ComplexHeatmap.html (release 3.6 for R 3.4)
# https://github.com/jokergoo/ComplexHeatmap
# https://jokergoo.github.io/ComplexHeatmap-reference/book/


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

library(ComplexHeatmap)

# setwd(here("src/expression_specificity"))

# ======================================================================= #
# ================================ FUNCTIONS ================================= #
# ======================================================================= #


# ======================= ES metrics (transformed) per anno ======================= #
# heatmap: genes x ES_metrics

make_es_heatmap.per_anno.genesXes_metrics <- function(sem_obj, dataset_prefix) {
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
      file.out <- sprintf("heatmap_complex.sem_transformed.%s__%s.pdf", dataset_prefix, anno_name)
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
}


# ======================= meta SEM (ALL ANNOTATIONS) ======================= #
### heatmap: genes x annotations
### THIS PLOT IS TOO BIG TO REALLY MAKE USE OF IT.

make_es_heatmap.es_meta.genesXannotations <- function(sem_obj, dataset_prefix) {
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
    file.out <- sprintf("heatmap_complex.sem_meta.%s__%s.pdf", dataset_prefix, sem_meta_name)
    pdf(file.out, width = 40, height = 20)
    time_elapsed <- system.time(draw(h1))
    print(time_elapsed)
    dev.off()
  }
}





##############################################################################################################
#############################################  RUN AS STAND-ALONE  ###########################################
##############################################################################################################

### Load data
# load(file="mousebrain.sem_obj.RData") # human

### Run functions
# make_es_heatmap.per_anno.genesXes_metrics(sem_obj)
# make_es_heatmap.es_meta.genesXannotations(sem_obj)
# print("Script done")


##############################################################################################################
######################################  Jon CODE ######################################
##############################################################################################################

# # ======================= Jon code ======================= #
# ht1 <- Heatmap(cellModEmbed_mat[,], 
#                cluster_rows = F,
#                #row_order = ,
#                cluster_columns = T, 
#                #show_row_dend = F, 
#                show_column_dend = F, 
#                show_heatmap_legend = T, 
#                show_row_names = F, 
#                show_column_names = T,
#                heatmap_legend_param = list(title = "Expression"), use_raster = T)
# pdf(sprintf("%s%s_%s_%s_cellModEmbed_cellOrder.pdf", dir_plots, data_prefix, run_prefix, fuzzyModMembership),h=max(10, nrow(cellModEmbed_mat) %/% 1000), w=max(8, ncol(cellModEmbed_mat) %/% 2))



##############################################################################################################
######################################  heatmaply : INTERACTIVE HEATMAP ######################################
##############################################################################################################

# ***PROBLEM*** ---> THIS DOES NOT SCALE TO ALL GENES

# library(tidyverse)
# # library(superheat) # https://rlbarter.github.io/superheat/basic-usage.html
# # library(pheatmap) # https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
# # library(d3heatmap) # interactive, https://cran.r-project.org/web/packages/d3heatmap/vignettes/Introduction.html
# library(heatmaply) # https://github.com/talgalili/heatmaply | Install this lib before you can save figs --> webshot::install_phantomjs()

# load(file="mousebrain.sem_obj.RData") # human

# for (sem_meta_name in names(sem_obj[["sem_meta"]])) {
#   print(sem_meta_name)
#   df.plot <- sem_obj[["sem_meta"]][[sem_meta_name]] %>% 
#     mutate(gene=sem_obj[["genes"]]) %>% 
#     slice(1:10) %>%
#     column_to_rownames(var="gene") %>% as.data.frame() # this step must come last to keep the rownames.
#   file.out <- sprintf("heatmap.sem_meta.%s.html", sem_meta_name)
#   heatmaply(df.plot,
#             main=sem_meta_name,
#             showticklabels = c(FALSE, FALSE),
#             file=file.out)
# }
