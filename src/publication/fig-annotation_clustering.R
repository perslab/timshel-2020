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

source(here("src/publication/lib-load_pub_lib_functions.R"))
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

### Annotation metadata
df.metadata <- get_metadata(dataset_prefix)


if (dataset_prefix == "mousebrain_all") {
  df.metadata <- df.metadata %>% mutate(Region = case_when(
    Region == "Midbrain dorsal" ~ "Midbrain",
    Region == "Midbrain dorsal,Midbrain ventral" ~ "Midbrain",
    Region == "Hippocampus,Cortex" ~ "Hippocampus/Cortex",
    TRUE ~ as.character(Region))
  )
}

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

library(ggraph) # 2019-02-19: CRAN only has ggraph_1.0.2 | use devtools::install_github('thomasp85/ggraph') to install ggraph v1.1 since it is build on tidygraph
# devtools::install_github('thomasp85/ggraph')
library(tidygraph)


### make tidygraph
gr <- as_tbl_graph(dend)
gr

### add meta-data
gr <- gr %>% 
  activate(nodes) %>% 
  left_join(df.metadata, by=c("label"="annotation"))

### add flag for prioritized cell-types
gr <- gr %>% 
  activate(nodes) %>% 
  mutate(flag_prioritized = if_else(label %in% filter.annotations, TRUE, FALSE))

### Get annotation colors
colormap.annotations <- get_color_mapping.prioritized_annotations_bmi(dataset="mousebrain")
colormap.annotations
### Create TaxonomyRank color mapping 
# [only needed because we manually want to set colors of annotations, and geom_node_text() and geom_node_point() use the same 'color scale']
colormap.taxonomy <- scales::hue_pal()(n_distinct(df.metadata$TaxonomyRank4_reduced1)) # these are default ggplot colors
names(colormap.taxonomy) <- sort(unique(df.metadata$TaxonomyRank4_reduced1)) # we sort to ensure consistency of results
### Combine color mappings
colormap.nodes <- c(colormap.annotations, colormap.taxonomy)
colormap.nodes

# ====================== "Circular"/radial dendrogram ==================== #
### REF rotate/angle geom_node_text for ggraph circular plot : https://stackoverflow.com/questions/43153004/how-to-read-a-text-label-in-ggraph-radial-graph
layout <- create_layout(gr, layout = "dendrogram", circular=TRUE) # circular
# head(as.tibble(layout)) # ---> contains x,y coordinates.
p <- ggraph(layout) + 
  # geom_edge_elbow() + # draw all lines [don't use this with geom_edge_elbow2() as it will draw the same lines twice]
  geom_edge_elbow2(aes(color=node.TaxonomyRank4_reduced1)) + # draw all lines | color 'leaf edges' | REF: see under "2-variant" | https://www.data-imaginist.com/2017/ggraph-introduction-edges/
  geom_node_point(aes(filter=(flag_prioritized==TRUE), color=label), size=1.5) + # color prioritized annotation leaf nodes points
  # geom_node_point(aes(filter=(flag_prioritized==TRUE)), color="red") + # ALT: use this if you want to color prioritized annotations *red*
  geom_node_text(aes(filter=((leaf==TRUE) & (flag_prioritized!=TRUE)), # leaf node labels
                     color=TaxonomyRank4_reduced1, label=label,
                     x=x*1.05, y=y*1.05,
                     angle = -((-node_angle(x, y)+90)%%180)+90
                     ),
                 hjust='outward',
                 size=rel(1),
                 show.legend=F) +
  geom_node_text(aes(filter=(flag_prioritized==TRUE), # leaf node labels, prioritized cell-types
                     color=TaxonomyRank4_reduced1, label=label,
                     x=x*1.05, y=y*1.05,
                     angle = -((-node_angle(x, y)+90)%%180)+90
                     ),
                 hjust='outward',
                 size=rel(1.4),
                 fontface = "bold",
                 show.legend=F) +
  ### Color nodes points by main figure colors
  ### Apparently, geom_node_text() and geom_node_point() use the same scale_color_*, so the argument to values=X need to contain color mapping for both coloring
  scale_color_manual(values=colormap.nodes, guide=FALSE) + 
  ### Guides: don't display
  guides(node_color="none",
         edge_color="none") + 
  theme_graph() # base_family="Helvetica"
p <- p + coord_fixed() # equal x and y axis scale is appropriate for circular plot
p

file.out <- sprintf("figs/fig_%s.dendrogram_ggraph.circular.pdf", dataset_prefix)
ggsave(plot=p, filename=file.out, width=6, height=6)

# ====================== "Linear" dendrogram [WORKS, but should be updated] ==================== #
### TODO: update geom_*(filter=...) to reflect the lastest version in the circular version
### The linear and circular plot only differ by their geom_node_text()
### And coord_fixed()
# layout <- create_layout(gr, layout = "dendrogram", circular=FALSE)
# p <- ggraph(layout) + 
#   geom_edge_elbow() +
#   geom_edge_elbow2(aes(color = node.TaxonomyRank4_reduced1)) +
#   geom_node_point(aes(filter=(flag_prioritized==TRUE)), color="red") + 
#   geom_node_text(aes(filter=(leaf==TRUE), color=TaxonomyRank4_reduced1, label=label), 
#                  angle=90, hjust=1, size=rel(0.8), nudge_y=-0.2, show.legend=F) +
#   geom_node_text(aes(filter=(flag_prioritized==TRUE), color=TaxonomyRank4_reduced1, label=label), 
#                  angle=90, hjust=1, size=rel(1), nudge_y=-0.2, show.legend=F) +
#   guides(node_color="none",
#          edge_color="none") + 
#   theme_graph()
# p
# file.out <- sprintf("figs/fig_%s.dendrogram_ggraph.linear.pdf", dataset_prefix)


# ======================================================================= #
# ====================== [dendextend] plot dendrogram ==================== #
# ======================================================================= #

library(dendextend)

### Simple base plot
# dend %>% plot()

if (FALSE) { # WORKS and looks ok
  ### Prep
  n_labels <- length(labels(dend))
  labels_col <- rep("gray", n_labels)
  labels_col[labels(dend) %in% filter.annotations] <- "red"
  leaves_cex <- rep(0, n_labels)
  leaves_cex[labels(dend) %in% filter.annotations] <- 1
  leaves_col <- rep("gray", n_labels)
  leaves_col[labels(dend) %in% filter.annotations] <- "red"
  
  ### Plot
  pdf(sprintf("figs/fig_%s.dendrogram_dendextend.pdf",dataset_prefix), width=45, height=12)
  dend %>% 
    set("labels_col", labels_col) %>%
    set("branches_k_color", k=20) %>% # branch color (k-means clustering)
    set("leaves_pch", 19) %>%  # node point type
    set("leaves_cex", leaves_cex) %>%  # node point size
    set("leaves_col", leaves_col) %>% #node point color
    plot()
  dev.off()
}

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

