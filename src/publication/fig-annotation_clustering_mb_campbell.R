############### SYNOPSIS ###################
# AIM: Expression specificity dendrogram for mousebrain hypothalamus

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


get_mousebrain_campbell_integrated_dendrogram <- function(save_fig) {
  
  # ======================================================================= #
  # ================================== MOUSEBRAIN ========================= #
  # ======================================================================= #
  dataset_prefix <- "mousebrain_all"
  df.metadata.mb <- get_metadata(dataset_prefix)
  ### Annotation ESmu
  file.es <- here(sprintf("out/es/%s.mean.csv.gz", dataset_prefix)) # New ESmu values
  df.es <- read_csv(file.es) # genes x cell-types
  
  ### Filter
  filter.annotations <- get_annotations.mousebrain.hypothalamus()
  df.es.mb <- df.es %>% select(gene, filter.annotations)
  
  # ======================================================================= #
  # ================================== CAMPBELL ========================= #
  # ======================================================================= #
  dataset_prefix <- "campbell_lvl2"
  df.metadata.cb <- get_metadata(dataset_prefix)
  
  ### Annotation ESmu
  file.es <- here(sprintf("out/es/%s.mean.csv.gz", dataset_prefix)) # New ESmu values
  df.es.cb <- read_csv(file.es) # genes x cell-types
  
  ### Filter | *CONSIDER ONLY INCLUDING CAMPBELL NEURONS*
  # df.es.cb %>% select(gene, filter.annotations)
  
  # ======================================================================= #
  # ================================== MERGE DATA FRAMES ========================= #
  # ======================================================================= #
  
  # dim(df.es.cb) # 12509    65
  # dim(df.es.mb) # 15071    18
  df.es <- inner_join(df.es.cb, df.es.mb, by="gene")
  
  # ======================================================================= #
  # ================================== MAKE META-DATA ========================= #
  # ======================================================================= #
  
  df.metadata.join <- tibble(
    annotation=colnames(df.es)[-1] # exclude gene
  )
  
  df.metadata.join <- df.metadata.join %>% mutate(dataset=if_else(annotation %in% df.metadata.mb$annotation, "Mousebrain", "Campbell"))
  df.metadata.join
  
  # ======================================================================= #
  # ======================== CALCULATE CORR/DENDROGRAM ===================== #
  # ======================================================================= #
  
  # ============================ Calculate correlation =============================== #
  cormat.es <- cor(df.es %>% select(-gene), method="pearson") 
  
  # ========================= [base] Calculate dendrogram ============================= #
  
  ### Compute distances and hierarchical clustering
  # dd.eucledian <- dist(t(df.es %>% select(-gene)), method = "euclidean") # computes distances between the rows of a data matrix
  dd.corr <- as.dist(1-cormat.es)
  hc <- hclust(dd.corr, method = "average") # Hierarchical clustering 
  # single linkage = nearest neighbour. 
  # complete linkage = farthest neighbour
  # ward distance does not have easy interpretable dedrogram 'height' for correlation distances.
  dend <- as.dendrogram(hc) # Turn the object into a dendrogram.
  # plot(dend)
  
  
  # ======================================================================= #
  # =========================== plot dendrogram =========================== #
  # ======================================================================= #
  
  ### 'Linear'
  p.linear <- plot_es_dendrogram.mb_campbell(dend, df.metadata.join, circular=FALSE)
  if (save_fig) {
    file.out <- sprintf("figs/fig_clustering.%s.dendrogram.linear.pdf", dataset_prefix="integrated_campbell_mousebrain")
    ggsave(plot=p.linear, filename=file.out, width=9, height=4)
  }
  
  # ### Circular
  # p.circular <- plot_es_dendrogram.mb_campbell(dend, df.metadata.join, circular=TRUE)
  # file.out <- sprintf("figs/fig_%s.dendrogram.circular.pdf", dataset_prefix="integrated_campbell_mousebrain")
  # ggsave(plot=p.circular, filename=file.out, width=6, height=6)
  
  
  list.res <- list(plot=p.linear, dendrogram=dend)
  return(list.res)
}

# interactive(): returns TRUE when R is being used interactively and FALSE otherwise
# REF: https://stackoverflow.com/a/2968404/6639640
if (!interactive()) { # function will run if script is sourced when running non-interactively (e.g. called via Rscript)
  setwd(here("src/publication"))
  get_mousebrain_campbell_integrated_dendrogram(save_fig=TRUE)
}


