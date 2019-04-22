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

library(corr)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #

### Tabula muris
dataset_prefix <- "tabula_muris"
filter.annotations <- get_prioritized_annotations_bmi(dataset="tabula_muris")

### Mousebrain
# dataset_prefix <- "mousebrain_all"
# filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")

# ======================================================================= #
# ======================== LOAD DATA AND PRE-PROCESS ===================== #
# ======================================================================= #

### Annotation metadata
df.metadata <- get_metadata(dataset_prefix)

### Annotation ESmu
# file.es <- here(sprintf("data/genes_cell_type_specific/%s.mean.csv.gz", dataset_prefix))
file.es <- here(sprintf("out/es/%s.mean.csv.gz", dataset_prefix)) # New ESmu values
df.es <- read_csv(file.es) # genes x cell-types

# ============================ Format annotation names =============================== #
### DOWNSTREAM ANNOTATION SELECTIONS FAIL WHEN formatting names because filter.annotations and colormaps are tied to non-formatted names
# NB: this also works for mousebrain

# ### Metadata
# df.metadata <- df.metadata %>% mutate(annotation=utils.rename_annotations.tabula_muris(annotation, style="tissue - celltype", check_all_matches=F))
# 
# ### df.es 
# tmp.oldnames <- colnames(df.es)
# tmp.newnames <- utils.rename_annotations.tabula_muris(tmp.oldnames, style="tissue - celltype", check_all_matches=F)
# df.es <- df.es %>% rename_at(vars(tmp.oldnames), ~ tmp.newnames) # REF https://stackoverflow.com/questions/20987295/rename-multiple-columns-by-names

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
p.circular <- plot_es_dendrogram(dend, df.metadata, dataset_prefix, circular=TRUE)
file.out <- sprintf("figs/fig_clustering.%s.dendrogram.circular.pdf", dataset_prefix)
ggsave(plot=p.circular, filename=file.out, width=6, height=6)

### 'Linear'
p.linear <- plot_es_dendrogram(dend, df.metadata, dataset_prefix, circular=FALSE)
file.out <- sprintf("figs/fig_clustering.%s.dendrogram.linear.pdf", dataset_prefix)
ggsave(plot=p.linear, filename=file.out, width=9, height=4)

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

plot_corrplot_all <- function(){
  corrplot(cormat.es, 
           order = "hclust", 
           hclust.method = "ward.D2",
           addrect = 10,
           addgrid.col=NA, # don't add grid
           tl.cex=0.2, # Numeric, for the size of text label (variable names).
           tl.col = "black")
}

plot_corrplot_selected <- function(){
  corrplot.mixed(tmp.cor, 
                 order = "hclust",
                 hclust.method = "ward.D2",
                 tl.col = "black",
                 tl.cex=0.6 # Numeric, for the size of text label (variable names)
                 )
}

### All cell-types
pdf(sprintf("figs/fig_clustering.%s.corrplot.all.pdf",dataset_prefix), width=15, height=15)
plot_corrplot_all()
dev.off()

### Selected cell-types
tmp.cor <- cor(df.es %>% select(one_of(filter.annotations), -gene), method="pearson")  # one_of(): variables in character vector. You need this when mixing unquoted names and character vector
pdf(sprintf("figs/fig_clustering.%s.corrplot.bmi_celltypes.pdf",dataset_prefix), width=15, height=15)
plot_corrplot_selected()
dev.off()

### COMBINING PLOTS
pdf(sprintf("figs/fig_clustering.%s.corrplot.combined.pdf",dataset_prefix), width=15, height=15)
par(mfrow=c(1,2), # mfrow=c(nrows, ncols)
    oma=c(2,7,7,2) # bottom, left, top, and right.
    ) 
plot_corrplot_all()
plot_corrplot_selected()
dev.off() # Every time a new device is opened par() will reset

### PAR DOCS
# mar: A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. The default is c(5, 4, 4, 2) + 0.1.
# oma: A vector of the form c(bottom, left, top, right) giving the size of the outer margins in lines of text.

# ======================================================================= #
# ========================= COMBINING PLOTS - LEFTOVER ========================= #
# ======================================================================= #

### layout | not ggplots
# layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
#        widths=c(1,1), heights=c(1,1))
# plot_corrplot_all()
# plot_corrplot_selected()

### cowplot ref REF: https://bioinfo.iric.ca/introduction-to-cowplot/
### ---> displays baseplots weird
# library(cowplot)
# row1 <- plot_grid(p.linear) 
# row2 <- plot_grid(plot_corrplot_all, plot_corrplot_selected, nrow=1, rel_widths=c(1, 1)) # notice we pass the function, but do not call it
# row12 <- plot_grid(row1, row2, ncol=1)
# row12

# ======================================================================= #
# ================================ LEFTOVERS ============================ #
# ======================================================================= #


### Does not look good for many cell-types - need to adjust text size and more, and I don't want to.
# pdf("plot.corrplot.mixed.pdf", width=20, height=20)
# corrplot.mixed(cormat.es)
# dev.off()

### THIS version of 'selected FDR' GIVES WRONG RESULT BECAUSE focus() does not return diagonal matrix
# cormat.es.seleced <- as_cordf(cormat.es) %>% focus(filter.annotations, mirror=T) %>% as_matrix()
# corrplot.mixed(cormat.es.seleced, is.corr=FALSE) # is.corr=FALSE convert to correlation matrix

