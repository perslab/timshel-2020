############### SYNOPSIS ###################
# AIM: Plot expression specificity:
# - dendrogram 
# - corrologram (correlation heatmap)

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
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #

### Tabula muris
# dataset_prefix <- "tabula_muris"
# filter.annotations <- get_prioritized_annotations_bmi(dataset="tabula_muris")

### Mousebrain
dataset_prefix <- "mousebrain_all"
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")

# ======================================================================= #
# ======================== LOAD DATA AND PRE-PROCESS ===================== #
# ======================================================================= #

### Annotation metadata
df.metadata <- get_metadata(dataset_prefix)

### Annotation ESmu
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

# source(here("src/publication/lib-load_pub_lib_functions.R"))

### Circular
p <- plot_es_dendrogram(dend, df.metadata, dataset_prefix, circular=TRUE)
p <- p + theme(plot.margin = unit(c(3,3,3,3), "cm")) # (t, r, b, l) all margin
file.out <- sprintf("figs/fig_%s.dendrogram.circular.pdf", dataset_prefix)
ggsave(plot=p, filename=file.out, width=6, height=6)

### 'Linear'
p <- plot_es_dendrogram(dend, df.metadata, dataset_prefix, circular=FALSE)
p <- p + theme(plot.margin = unit(c(1,1,4,1), "cm")) # (t, r, b, l) widen bottom margin
file.out <- sprintf("figs/fig_%s.dendrogram.linear.pdf", dataset_prefix)
ggsave(plot=p, filename=file.out, width=9, height=4)

# if (dataset_prefix == "mousebrain_all") {
#   p <- p + theme(plot.margin = unit(c(1,1,2,1), "cm")) # (t, r, b, l) widen bottom margin
# } else if (dataset_prefix == "tabula_muris") {
#   p <- p + theme(plot.margin = unit(c(1,1,4,1), "cm")) # (t, r, b, l) widen bottom margin
# }

p
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

