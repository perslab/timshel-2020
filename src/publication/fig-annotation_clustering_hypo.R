############### SYNOPSIS ###################
# AIM: Expression specificity dendrogram for hypothalamus

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


library(ggraph)
library(tidygraph)
library(corrr)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))



# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #


get_hypothalamus_integrated_dendrogram <- function(save_fig, clustering_method=c("correlation", "wards")) {
  clustering_method <- match.arg(clustering_method)
  # wards: bool param. If True, then Wards-Eucledian cluster will be used. Else Average_linkage-Correalation_dist will be used
  
  setwd(here("src/publication"))
  # ======================================================================= #
  # ================================== LOAD DATA ========================= #
  # ======================================================================= #
  
  df.metadata <- get_metadata("hypothalamus")
  
  df.es <- get_combined_hypo_es(merge="inner") # *inner*: only overlapping genes will be used for clustering
  
  # ======================================================================= #
  # ================================== MAKE META-DATA ========================= #
  # ======================================================================= #
  
  # df.metadata.join <- tibble(
  #   annotation=colnames(df.es)[-1] # exclude gene
  # )
  
  df.metadata.join <- df.metadata %>% mutate(dataset=annotation_prefix,
                                             annotation=annotation_fmt)
  
  # ======================================================================= #
  # ======================== CALCULATE CORR/DENDROGRAM ===================== #
  # ======================================================================= #
  
  ### Compute distances and hierarchical clustering
  if (clustering_method=="correlation") {
    cormat.es <- cor(df.es %>% select(-gene), method="pearson") 
    dd.corr <- as.dist(1-cormat.es)
    hc <- hclust(dd.corr, method = "average") # Hierarchical clustering 
  } else if (clustering_method=="wards") { 
    dd.eucledian <- dist(t(df.es %>% select(-gene)), method = "euclidean") # computes distances between the rows of a data matrix
    hc <- hclust(dd.eucledian, method = "ward.D2") # Hierarchical clustering 
  } else {
    stop("error")
  }
  
  # single linkage = nearest neighbour. 
  # complete linkage = farthest neighbour
  # ward distance does not have easy interpretable dedrogram 'height' for correlation distances.
  dend <- as.dendrogram(hc) # Turn the object into a dendrogram.
  # plot(dend)
  
  
  # ======================================================================= #
  # ========= plot dendrogram: plot_es_dendrogram.multi_dataset =========== #
  # ======================================================================= #
  # plot_es_dendrogram.multi_dataset() does not highlight cell-types
  
  ### 'Linear' - Ward dist
  # p.linear <- plot_es_dendrogram.multi_dataset(dend, df.metadata.join, circular=FALSE, var_aes="taxonomy_lvl1")
  # ggsave(filename=sprintf("figs/fig_clustering.hypothalamus.dendrogram.wards_tax1.linear.pdf"), width=16, height=8)
  # p.linear <- plot_es_dendrogram.multi_dataset(dend, df.metadata.join, circular=FALSE, var_aes="taxonomy_lvl2")
  # ggsave(filename=sprintf("figs/fig_clustering.hypothalamus.dendrogram.wards_tax2.linear.pdf"), width=16, height=8)
  # p.linear <- plot_es_dendrogram.multi_dataset(dend, df.metadata.join, circular=FALSE, var_aes="dataset")
  # ggsave(filename=sprintf("figs/fig_clustering.hypothalamus.dendrogram.wards_dataset.linear.pdf"), width=16, height=8)
  ### 'Linear' - Cor dist
  if (save_fig) {
    p.linear <- plot_es_dendrogram.multi_dataset(dend, df.metadata.join, circular=FALSE, var_aes="taxonomy_lvl1")
    ggsave(filename=sprintf("figs/fig_clustering.hypothalamus.dendrogram.%s.tax1.linear.pdf", clustering_method), width=16, height=8)
    p.linear <- plot_es_dendrogram.multi_dataset(dend, df.metadata.join, circular=FALSE, var_aes="taxonomy_lvl2")
    ggsave(filename=sprintf("figs/fig_clustering.hypothalamus.dendrogram.%s.tax2.linear.pdf", clustering_method), width=16, height=8)
    p.linear <- plot_es_dendrogram.multi_dataset(dend, df.metadata.join, circular=FALSE, var_aes="dataset")
    ggsave(filename=sprintf("figs/fig_clustering.hypothalamus.dendrogram.%s.dataset.linear.pdf", clustering_method), width=16, height=8)
  }
  # ======================================================================= #
  # =============== plot dendrogram: plot_es_dendrogram() ================= #
  # ======================================================================= #
  # Highlight cell-types
  ### 'Linear'
  # p.linear <- plot_es_dendrogram(dend, df.metadata.join, dataset_prefix="hypothalamus", label_only_prioritized=F, circular=FALSE, show_legend=T)
  p.linear <- plot_es_dendrogram(dend, df.metadata.join, dataset_prefix="hypothalamus", label_only_prioritized=F, circular=FALSE, show_legend=F)
  if (save_fig) {
    file.out <- sprintf("figs/fig_clustering.hypothalamus.dendrogram.%s.linear.pdf", clustering_method)
    ggsave(plot=p.linear, filename=file.out, width=12, height=5)
  }
  
  
  list.res <- list(plot=p.linear, dendrogram=dend)
  return(list.res)
}

# interactive(): returns TRUE when R is being used interactively and FALSE otherwise
# REF: https://stackoverflow.com/a/2968404/6639640
if (!interactive()) { # function will run if script is sourced when running non-interactively (e.g. called via Rscript)
  setwd(here("src/publication"))
  get_hypothalamus_integrated_dendrogram(save_fig=TRUE, clustering_method="correlation")
}


