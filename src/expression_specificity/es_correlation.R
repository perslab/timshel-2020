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


dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library

setwd(here("src/expression_specificity"))


# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #

### Tabula muris
dataset_prefix <- "tabula_muris"
filter.celltypes <- c("Brain_Non-Myeloid.neuron","Brain_Non-Myeloid.oligodendrocyte_precursor_cell")

### Mousebrain
# dataset_prefix <- "mousebrain_all"
# filter.celltypes <- c("TEGLU23","DEINH3","MEGLU1","MEINH2","DEGLU5","MEGLU10","TEGLU17","MEGLU11","TEGLU4","DEGLU4","TEINH12") # BMI_UKBB_Loh2018 FDR sign.

# ======================================================================= #
# ============================ LOAD DATA =============================== #
# ======================================================================= #

file.es <- here(sprintf("data/genes_cell_type_specific/%s.mean.csv.gz", dataset_prefix))
df.es <- read_csv(file.es) # genes x cell-types

# ======================================================================= #
# ============================ Calculate correlation =============================== #
# ======================================================================= #

cormat.es <- cor(df.es %>% select(-gene), method="pearson") 

# ======================================================================= #
# ============================ [base] Calculate dendrogram =============================== #
# ======================================================================= #

### Compute distances and hierarchical clustering
# dd.eucledian <- dist(t(df.es %>% select(-gene)), method = "euclidean") # computes distances between the rows of a data matrix
dd.corr <- as.dist(1-cormat.es)
hc <- hclust(dd.corr, method = "ward.D2") # Hierarchical clustering 
dend <- as.dendrogram(hc) # Turn the object into a dendrogram.
plot(dend)

# ======================================================================= #
# ============================ Export clustering order =============================== #
# ======================================================================= #

file_out.clustering_order <- here(sprintf("data/genes_cell_type_specific/%s.hclust_order.csv", dataset_prefix))
df.clustering_order <- tibble(cluster_order=seq_along(labels(dend)), annotation=labels(dend))
df.clustering_order %>% write_csv(file_out.clustering_order)

# ======================================================================= #
# ============================ [dendextend] plot dendrogram =============================== #
# ======================================================================= #

library(dendextend)

dend %>% plot()

n_labels <- length(labels(dend))
labels_col <- rep("gray", n_labels)
labels_col[labels(dend) %in% filter.celltypes] <- "red"
leaves_cex <- rep(0, n_labels)
leaves_cex[labels(dend) %in% filter.celltypes] <- 1
leaves_col <- rep("gray", n_labels)
leaves_col[labels(dend) %in% filter.celltypes] <- "red"

pdf(sprintf("plot.%s.dendrogram.pdf",dataset_prefix), width=45, height=12)
dend %>% 
  set("labels_col", labels_col) %>%
  set("branches_k_color", k=20) %>% # branch color (k-means clustering)
  set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", leaves_cex) %>%  # node point size
  set("leaves_col", leaves_col) %>% #node point color
  plot()
dev.off()


# ======================================================================= #
# ============================ [corrplot] Correlation plots =============================== #
# ======================================================================= #

### https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
library(corrplot)

### Does not look good for many cell-types - need to adjust text size and more, and I don't want to.
# pdf("plot.corrplot.mixed.pdf", width=20, height=20)
# corrplot.mixed(cormat.es)
# dev.off()

pdf(sprintf("plot.%s.corrplot.hclust.pdf",dataset_prefix), width=20, height=20)
corrplot(cormat.es, order = "hclust", addrect = 10)
dev.off()

### selected FDR using corrplot
tmp.cor <- cor(df.es %>% select(one_of(filter.celltypes), -gene), method="pearson")  # one_of(): variables in character vector. You need this when mixing unquoted names and character vector
pdf(sprintf("plot.%s.corrplot.mixed_hclust_fdr_only.pdf",dataset_prefix), width=15, height=15)
corrplot.mixed(tmp.cor, order = "hclust")
dev.off()


### selected FDR using corrplot --> GIVES WRONG RESULT BECAUSE focus() does not return diagonal matrix
# cormat.es.seleced <- as_cordf(cormat.es) %>% focus(filter.celltypes, mirror=T) %>% as_matrix()
# corrplot.mixed(cormat.es.seleced, is.corr=FALSE) # is.corr=FALSE convert to correlation matrix


# ======================================================================= #
# ============================ [corrr] Correlation plots =============================== #
# ======================================================================= #

### corrr: exploring correlations
### Code: https://github.com/drsimonj/corrr/
### Blog: https://drsimonj.svbtle.com/exploring-correlations-in-r-with-corrr
library(corrr)
# *OBS*: focus() does NOT return a diagonal matrix. That is, it is not 'ordered'.

df.es_cordf <- as_cordf(cormat.es) # instead of calling correlate(df.es %>% select(-gene))

### selecting cell-types
# df.es_cordf %>% focus(matches("^MSN"), mirror=T)
# df.es_cordf %>% focus(filter.celltypes, mirror=T) %>% rearrange(method="HC") %>% shave() %>% rplot(print_cor=T) # shave() looks weird when the matrix is not diagonal
df.es_cordf %>% focus(filter.celltypes, mirror=T) %>% rearrange(method="HC") %>% rplot(print_cor=T)
ggsave(sprintf("plot.%s.corrr_hc.pdf",dataset_prefix), w=10, h=10)

# df.es_cordf %>% focus(filter.celltypes, mirror=T) %>% select(sort(current_vars())) # ... trying to order the colnames but then rowname is included, so it becomes to much of a mess to solve

### INDEX
# as_cordf	Coerce lists and matrices to correlation data frames
# as_matrix	Convert a correlation data frame to matrix format
# correlate	Correlation Data Frame
# fashion	Fashion a correlation data frame for printing.
# first_col	Add a first column to a data.frame
# focus	Focus on section of a correlation data frame.
# focus_	Focus on section of a correlation data frame.
# focus_if	Conditionally focus correlation data frame
# network_plot	Network plot of a correlation data frame
# pair_n	Number of pairwise complete cases.
# rearrange	Re-arrange a correlation data frame
# rplot	Plot a correlation data frame.
# shave	Shave off upper/lower triangle.
# stretch	Stretch correlation data frame into long format.

# ======================================================================= #
# ============================ [ggcor] Correlation plots =============================== #
# ======================================================================= #
library(GGally)
# ggcorr (used for lira lab stats):
# part of GGally
# intro: https://briatte.github.io/ggcorr/ 
# Repo: https://github.com/briatte/ggcorr

# ggcorr(data = NULL, cor_matrix = cormat.es)
# ggcorr(df.lab, geom = "circle", hjust = 1, size = 2, layout.exp = 5) # label = TRUE
# ggcorr(df.lab, label = TRUE, label_size = 1.5, label_round = 2, label_alpha = TRUE, hjust = 1, size = 2, layout.exp = 5)



# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


