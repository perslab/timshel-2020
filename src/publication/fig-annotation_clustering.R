############### SYNOPSIS ###################
# AIM: Explore and plot espression specificity correlations
# Dendrogram
# Correlation heatmap

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

setwd(here("src/publication"))


# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #

### Tabula muris
# dataset_prefix <- "tabula_muris"
# filter.annotations <- c("Brain_Non-Myeloid.neuron","Brain_Non-Myeloid.oligodendrocyte_precursor_cell")

### Mousebrain
dataset_prefix <- "mousebrain_all"
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")

# ======================================================================= #
# ======================== LOAD DATA AND PRE-PROCESS ===================== #
# ======================================================================= #

file.es <- here(sprintf("data/genes_cell_type_specific/%s.mean.csv.gz", dataset_prefix))
df.es <- read_csv(file.es) # genes x cell-types

# ============================ Calculate correlation =============================== #
cormat.es <- cor(df.es %>% select(-gene), method="pearson") 

# ========================= [base] Calculate dendrogram ============================= #

### Compute distances and hierarchical clustering
# dd.eucledian <- dist(t(df.es %>% select(-gene)), method = "euclidean") # computes distances between the rows of a data matrix
dd.corr <- as.dist(1-cormat.es)
hc <- hclust(dd.corr, method = "ward.D2") # Hierarchical clustering 
dend <- as.dendrogram(hc) # Turn the object into a dendrogram.
plot(dend)

# ============================ Export clustering order =============================== #
# file_out.clustering_order <- here(sprintf("data/genes_cell_type_specific/%s.hclust_order.csv", dataset_prefix))
# df.clustering_order <- tibble(cluster_order=seq_along(labels(dend)), annotation=labels(dend))
# df.clustering_order %>% write_csv(file_out.clustering_order)

# ======================================================================= #
# ====================== [gggraph] plot dendrogram ==================== #
# ======================================================================= #

library(ggraph) # 2019-02-19: CRAN only has ggraph_1.0.2 | use devtools::install_github('thomasp85/ggraph') to install ggraph v1.1 since it is build on tidygraph
# devtools::install_github('thomasp85/ggraph')
library(tidygraph)

# TODO: write code here. 
# SEE EVERNOTE.
# SEE fig-module_network

file.out <- sprintf("figs/fig_%s.dendrogram_gggraph.pdf", dataset_prefix)

# ======================================================================= #
# ====================== [dendextend] plot dendrogram ==================== #
# ======================================================================= #

library(dendextend)

dend %>% plot()

n_labels <- length(labels(dend))
labels_col <- rep("gray", n_labels)
labels_col[labels(dend) %in% filter.annotations] <- "red"
leaves_cex <- rep(0, n_labels)
leaves_cex[labels(dend) %in% filter.annotations] <- 1
leaves_col <- rep("gray", n_labels)
leaves_col[labels(dend) %in% filter.annotations] <- "red"

pdf(sprintf("figs/fig_%s.dendrogram_dendextend.pdf",dataset_prefix), width=45, height=12)
dend %>% 
  set("labels_col", labels_col) %>%
  set("branches_k_color", k=20) %>% # branch color (k-means clustering)
  set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", leaves_cex) %>%  # node point size
  set("leaves_col", leaves_col) %>% #node point color
  plot()
dev.off()


# ======================================================================= #
# ====================== [corrplot] Correlation plots ==================== #
# ======================================================================= #

### https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
library(corrplot)

### All cell-types
pdf(sprintf("figs/fig_%s.corrplot.hclust.pdf",dataset_prefix), width=20, height=20)
corrplot(cormat.es, order = "hclust", addrect = 10)
dev.off()

### Selected cell-types
tmp.cor <- cor(df.es %>% select(one_of(filter.annotations), -gene), method="pearson")  # one_of(): variables in character vector. You need this when mixing unquoted names and character vector
pdf(sprintf("figs/fig_%s.corrplot.mixed_hclust_fdr_celltypes_only.pdf",dataset_prefix), width=15, height=15)
corrplot.mixed(tmp.cor, order = "hclust")
dev.off()


### Does not look good for many cell-types - need to adjust text size and more, and I don't want to.
# pdf("plot.corrplot.mixed.pdf", width=20, height=20)
# corrplot.mixed(cormat.es)
# dev.off()


### THIS version of 'selected FDR' GIVES WRONG RESULT BECAUSE focus() does not return diagonal matrix
# cormat.es.seleced <- as_cordf(cormat.es) %>% focus(filter.annotations, mirror=T) %>% as_matrix()
# corrplot.mixed(cormat.es.seleced, is.corr=FALSE) # is.corr=FALSE convert to correlation matrix

