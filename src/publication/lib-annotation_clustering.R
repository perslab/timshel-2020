############### SYNOPSIS ###################
# Helper functions plotting expression specificity dendrogram 


# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)

library(ggraph) # 2019-02-19: CRAN only has ggraph_1.0.2 | use devtools::install_github('thomasp85/ggraph') to install ggraph v1.1 since it is build on tidygraph
# devtools::install_github('thomasp85/ggraph')
library(tidygraph)


# ======================================================================= #
# ====================== [gggraph] plot dendrogram ==================== #
# ======================================================================= #
plot_es_dendrogram <- function(dend, df.metadata, dataset_prefix, circular) {
  
  # ====================== Get dataset-specific params (color mappings) ==================== #
  if (dataset_prefix == "mousebrain_all") {
    filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
    colormap.annotations_highlight <- get_color_mapping.prioritized_annotations_bmi(dataset="mousebrain")
    sym_var.node_color <- sym("TaxonomyRank4_reduced1")
    sym_var.node_color_edge_elbow2 <- sym("node.TaxonomyRank4_reduced1")
  } else if (dataset_prefix == "tabula_muris") {
    filter.annotations <- get_prioritized_annotations_bmi(dataset="tabula_muris")
    colormap.annotations_highlight <- get_color_mapping.prioritized_annotations_bmi(dataset="tabula_muris")
    sym_var.node_color <- sym("tissue")
    sym_var.node_color_edge_elbow2 <- sym("node.tissue")
  } else {
    stop(sprintf("Wrong dataset_prefix: %s", dataset_prefix))
  }
  
  
  ### Create color mapping for 'all nodes'
  # [only needed because we manually want to set colors of annotations, and geom_node_text() and geom_node_point() use the same 'color scale']
  colormap.all_nodes <- scales::hue_pal()(n_distinct(df.metadata %>% pull(!!sym_var.node_color))) # these are default ggplot colors
  names(colormap.all_nodes) <- sort(unique(df.metadata %>% pull(!!sym_var.node_color))) # we sort to ensure consistency of results
  ### Combine color mappings
  colormap.nodes <- c(colormap.annotations_highlight, colormap.all_nodes)
  colormap.nodes
  
  
  # ====================== make tidygraph ==================== #
  gr <- as_tbl_graph(dend)
  ### add meta-data
  gr <- gr %>% 
    activate(nodes) %>% 
    left_join(df.metadata, by=c("label"="annotation")) # OK --> Warning message: Column `label`/`annotation` joining factor and character vector, coercing into character vector 
  ### add flag for prioritized cell-types
  gr <- gr %>% 
    activate(nodes) %>% 
    mutate(flag_prioritized = if_else(label %in% filter.annotations, TRUE, FALSE))
  
  
  
  # ====================== Dendrogram ==================== #
  ### REF rotate/angle geom_node_text for ggraph circular plot : https://stackoverflow.com/questions/43153004/how-to-read-a-text-label-in-ggraph-radial-graph
  layout <- create_layout(gr, layout = "dendrogram", circular=circular) # circular
  # head(as.tibble(layout)) # ---> contains x,y coordinates.
  p <- ggraph(layout) + 
    # geom_edge_elbow() + # draw all lines [don't use this with geom_edge_elbow2() as it will draw the same lines twice]
    geom_edge_elbow2(aes(color=!!sym_var.node_color_edge_elbow2)) + # draw all lines | color 'leaf edges' | REF: see under "2-variant" | https://www.data-imaginist.com/2017/ggraph-introduction-edges/
    geom_node_point(aes(filter=(flag_prioritized==TRUE), color=label), size=1.5) # color prioritized annotation leaf nodes points
    # geom_node_point(aes(filter=(flag_prioritized==TRUE)), color="red") + # ALT: use this if you want to color prioritized annotations *red*
    
  ### Circular
  ### The linear and circular plot only differ by their geom_node_text()
  ### And coord_fixed()
  if (circular) {
    p <- p + geom_node_text(aes(filter=((leaf==TRUE) & (flag_prioritized!=TRUE)), # leaf node labels
                       color=!!sym_var.node_color, label=label,
                       x=x*1.05, y=y*1.05,
                       angle = -((-node_angle(x, y)+90)%%180)+90),
                  hjust='outward',
                  size=rel(1),
                  show.legend=F) +
    geom_node_text(aes(filter=(flag_prioritized==TRUE), # leaf node labels, prioritized cell-types
                       color=!!sym_var.node_color, label=label,
                       x=x*1.05, y=y*1.05,
                       angle = -((-node_angle(x, y)+90)%%180)+90),
                    hjust='outward',
                    size=rel(1.4),
                    fontface = "bold",
                    show.legend=F)
    p <- p + coord_fixed( # equal x and y axis scale is appropriate for circular plot
      clip = 'off' # This keeps the labels from disappearing. It allows drawing of data points anywhere on the plot, including in the plot margins.
    )
  } else {
    p <- p + 
      geom_node_text(aes(filter=((leaf==TRUE) & (flag_prioritized!=TRUE)), # leaf node labels
                              color=!!sym_var.node_color, label=label),
                      size=rel(0.8), 
                      angle=90, hjust=1, nudge_y=-0.2, 
                      show.legend=F) +
      geom_node_text(aes(filter=(flag_prioritized==TRUE), # leaf node labels, prioritized cell-types
                             color=!!sym_var.node_color, label=label),
                     size=rel(1), 
                     angle=90, hjust=1, nudge_y=-0.2, 
                     show.legend=F)
    p <- p + coord_cartesian(
      clip = 'off' # This keeps the labels from disappearing. It allows drawing of data points anywhere on the plot, including in the plot margins.
    )
  }
  p <- p + 
    ### Color nodes points by main figure colors
    ### Apparently, geom_node_text() and geom_node_point() use the same scale_color_*, so the argument to values=X need to contain color mapping for both coloring
    scale_color_manual(values=colormap.nodes, guide=FALSE) + 
  ### Guides: don't display
    guides(node_color="none",
           edge_color="none") + 
    theme_graph() # base_family="Helvetica"
  return(p)
}

           
  

# ======================================================================= #
# ====================== [dendextend] plot dendrogram ==================== #
# ======================================================================= #

### WORKS and looks ok, but ggraph is better

if (FALSE) {
  library(dendextend)
  
  ### Simple base plot
  # dend %>% plot()
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


