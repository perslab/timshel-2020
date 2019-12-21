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


get_hypothalamus_integrated_dendrogram <- function(save_fig) {
  setwd(here("src/publication"))
  # ======================================================================= #
  # ================================== LOAD DATA ========================= #
  # ======================================================================= #
  
  ### Meta-data
  df.metadata <- get_metadata("hypothalamus")
  df.metadata <- df.metadata %>% mutate(annotation_uniq = paste0(specificity_id, "__", annotation)) # needed because some hyp annotations are duplicated across datasets
  
  ### ESmu matrix
  file.es <- here("out/es", paste0(get_scrna_seq_dataset_prefixes("hypo"), ".mu.csv.gz"))
  list.df.es <- map(file.es, read_csv) # genes x cell-types
  ### Rename ESmu annotations to avoid problem with duplicates during merge 
  names(list.df.es) <- get_scrna_seq_dataset_prefixes("hypo") # names is garantueed to be in same order as above
  list.df.es.prefixed <- list.df.es %>% imap(.f=function(df,specificity_id){colnames(df)[-1] <- paste0(specificity_id, "__", colnames(df)[-1]); df}) # genes x cell-types.
  # ^ x[-1]: first column is "gene"

  ### Merge (*inner*: only overlapping genes will be used for clustering)
  df.es <- plyr::join_all(list.df.es.prefixed, by="gene", type="inner", match="first") %>% as.tibble() # REF: https://stackoverflow.com/a/32066419/6639640. match argument should not matter (only speed)
  
  ### Map hypothalamus annotations to annotation_fmt
  annotation_uniq <- colnames(df.es)[-1] # exclude 'gene' column
  annotation_fmt <- df.metadata$annotation_fmt[match(annotation_uniq, df.metadata$annotation_uniq)]
  # df.tmp <- tibble(x=annotation_fmt, y=annotation_uniq) # test to show that the matching works
  colnames(df.es)[-1] <- annotation_fmt # ALTERNATIVE: use deframe() and rename(!!!x)
  
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
  
  # ============================ Calculate correlation =============================== #
  cormat.es <- cor(df.es %>% select(-gene), method="pearson") 
  
  # ========================= [base] Calculate dendrogram ============================= #
  
  ### Compute distances and hierarchical clustering
  # dd.eucledian <- dist(t(df.es %>% select(-gene)), method = "euclidean") # computes distances between the rows of a data matrix
  dd.corr <- as.dist(1-cormat.es)
  hc <- hclust(dd.corr, method = "average") # Hierarchical clustering 
  # hc <- hclust(dd.eucledian, method = "ward.D2") # Hierarchical clustering 
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
    ggsave(filename=sprintf("figs/fig_clustering.hypothalamus.dendrogram.cor_tax1.linear.pdf"), width=16, height=8)
    p.linear <- plot_es_dendrogram.multi_dataset(dend, df.metadata.join, circular=FALSE, var_aes="taxonomy_lvl2")
    ggsave(filename=sprintf("figs/fig_clustering.hypothalamus.dendrogram.cor_tax2.linear.pdf"), width=16, height=8)
    p.linear <- plot_es_dendrogram.multi_dataset(dend, df.metadata.join, circular=FALSE, var_aes="dataset")
    ggsave(filename=sprintf("figs/fig_clustering.hypothalamus.dendrogram.cor_dataset.linear.pdf"), width=16, height=8)
  }
  # ======================================================================= #
  # =============== plot dendrogram: plot_es_dendrogram() ================= #
  # ======================================================================= #
  # Highlight cell-types
  ### 'Linear'
  # p.linear <- plot_es_dendrogram(dend, df.metadata.join, dataset_prefix="hypothalamus", label_only_prioritized=F, circular=FALSE, show_legend=T)
  p.linear <- plot_es_dendrogram(dend, df.metadata.join, dataset_prefix="hypothalamus", label_only_prioritized=F, circular=FALSE, show_legend=F)
  if (save_fig) {
    file.out <- sprintf("figs/fig_clustering.%s.dendrogram.linear.pdf", dataset_prefix="hypothalamus")
    ggsave(plot=p.linear, filename=file.out, width=16, height=8)
  }
  
  
  list.res <- list(plot=p.linear, dendrogram=dend)
  return(list.res)
}

# interactive(): returns TRUE when R is being used interactively and FALSE otherwise
# REF: https://stackoverflow.com/a/2968404/6639640
if (!interactive()) { # function will run if script is sourced when running non-interactively (e.g. called via Rscript)
  setwd(here("src/publication"))
  get_hypothalamus_integrated_dendrogram(save_fig=TRUE)
}


