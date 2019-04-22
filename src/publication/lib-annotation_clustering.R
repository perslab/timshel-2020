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

library(RColorBrewer)


# ======================================================================= #
# ====================== [gggraph] plot dendrogram ==================== #
# ======================================================================= #

# ============= MB+CAMPBELL PLOT ==================== #
# SIMILAR TO plot_es_dendrogram() but with the approximate following changes
# - no "flag_priorized"  and highligthig
# - increased text size
# - added color label
# base_family = 'Helvetica' to avoid ggsave problem
plot_es_dendrogram.mb_campbell <- function(dend, df.metadata, circular) {
  
  # ====================== Get dataset-specific params (color mappings) ==================== #
  sym_var.node_color <- sym("dataset")
  sym_var.node_color_edge_elbow2 <- sym("node.dataset")
  n_color_categories <- n_distinct(df.metadata %>% pull(!!sym_var.node_color))
  
  ### Create color mapping for 'all nodes'
  # [only needed because we manually want to set colors of annotations, and geom_node_text() and geom_node_point() use the same 'color scale']
  # colormap.nodes <- scales::hue_pal()(n_distinct(df.metadata %>% pull(!!sym_var.node_color))) # these are default ggplot colors
  colormap.nodes <- brewer.pal(name="Set1", n=9) # n=9 max for Set1 | "Set2"/"Dark2" ok for colorblind
  if (n_color_categories < length(n_color_categories)) {
    stop(sprintf("Too few colors in current colormap for n=%s categories", n_color_categories))
  }
  colormap.nodes <- colormap.nodes[1:n_color_categories] # select the colors we need
  names(colormap.nodes) <- sort(unique(df.metadata %>% pull(!!sym_var.node_color))) # we sort to ensure consistency of results

  # ====================== make tidygraph ==================== #
  gr <- as_tbl_graph(dend)
  ### add meta-data
  gr <- gr %>% 
    activate(nodes) %>% 
    left_join(df.metadata, by=c("label"="annotation")) # OK --> Warning message: Column `label`/`annotation` joining factor and character vector, coercing into character vector 

  # ====================== Dendrogram ==================== #
  layout <- create_layout(gr, layout = "dendrogram", height=height, circular=circular) # circular
  p <- ggraph(layout)
  p <- p + geom_edge_elbow2(aes(color=!!sym_var.node_color_edge_elbow2)) 
  ### Circular
  if (circular) {
    p <- p + geom_node_text(aes(filter=(leaf==TRUE), # leaf node labels
                                color=!!sym_var.node_color, label=label,
                                x=x*1.05, y=y*1.05,
                                angle = -((-node_angle(x, y)+90)%%180)+90),
                            hjust='outward',
                            size=rel(1.2),
                            show.legend=F)
    p <- p + coord_fixed( 
      clip = 'off'
    )
    p <- p + theme_graph(base_family="Helvetica")
    p <- p + theme(plot.margin = unit(c(3,3,3,3), "cm")) # (t, r, b, l) extra margins
  } else { ### Linear
    p <- p + geom_node_text(aes(filter=(leaf==TRUE), # leaf node labels
                         color=!!sym_var.node_color, label=label),
                     size=rel(1.2), 
                     angle=90, hjust=1, nudge_y=-0.05, 
                     show.legend=F)
    p <- p + labs(y=expression(Distance~(rho)))
    p <- p + coord_cartesian(clip="off", 
                               ylim=c(0,max(layout$height)+0.1), 
                               expand=F 
                               ) 
    ### ensure that annotation are placed some distance away from the y-axis.
    n_leaf_nodes <- length(order.dendrogram(dend)) # stats::order.dendrogram: A vector with length equal to the number of leaves in the dendrogram is returned
    p <- p + xlim(c(-1,n_leaf_nodes+1)) 
    p <- p + theme_classic(base_family="Helvetica")
    p <- p + theme(axis.line.x=element_blank(), # remove all x-axis 
               axis.title.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.text.x=element_blank())
    p  <- p  + theme(plot.margin = unit(c(1,1,4,1), "cm")) # (t, r, b, l) widen bottom margin
  }
  ### Add points (after drawing edges)
  p <- p + geom_node_point(aes(filter=(leaf==TRUE), color=!!sym_var.node_color), size=1) # color prioritized annotation leaf nodes points
  ### Guides
  p <- p + 
    ### Color nodes points by main figure colors
    ### Apparently, geom_node_text() and geom_node_point() use the same scale_color_*, so the argument to values=X need to contain color mapping for both coloring
    scale_color_manual(values=colormap.nodes) + 
    guides(edge_color="none", node.color="none") # node_color="none", node_color=guide_legend()
  ### Legend
  p <- p + 
    theme(legend.position="top")
  return(p)
}




# ====================== MOUSEBRAIN PLOT ==================== #

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
  ### FIX HEIGHT issues: ggraph(dend, layout = 'dendrogram', height = height) | REF: https://github.com/thomasp85/ggraph/issues/175 
  layout <- create_layout(gr, layout = "dendrogram", height=height, circular=circular) # circular
  # head(as.tibble(layout)) # ---> contains x,y coordinates.
  p <- ggraph(layout)
  p <- p + geom_edge_elbow2(aes(color=!!sym_var.node_color_edge_elbow2)) # draw all lines | color 'leaf edges' | REF: see under "2-variant" | https://www.data-imaginist.com/2017/ggraph-introduction-edges/
  ### Circular
  ### The linear and circular plot differ by: aes(y=node.height), geom_node_text() and coord_fixed()
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
    p <- p + theme_graph(base_family="Helvetica")
    p <- p + theme(plot.margin = unit(c(3,3,3,3), "cm")) # (t, r, b, l) extra margins
  } else {
    # p <- p + geom_edge_elbow2(aes(y=node.height, color=!!sym_var.node_color_edge_elbow2)) # gives correct height if 'height' is not set in ggraph
    p <- p + geom_node_text(aes(filter=((leaf==TRUE) & (flag_prioritized!=TRUE)), # leaf node labels
                                color=!!sym_var.node_color, label=label),
                        size=rel(1),
                        angle=90, hjust=1, nudge_y=-0.05,
                        show.legend=F)
    p <- p + geom_node_text(aes(filter=(flag_prioritized==TRUE), # leaf node labels, prioritized cell-types
                               color=!!sym_var.node_color, label=label),
                       size=rel(1.3),
                       angle=90, hjust=1, nudge_y=-0.05,
                       show.legend=F)
    p <- p + labs(y=expression(Distance~(rho)))
    p <- p + coord_cartesian(clip="off", # This keeps the labels from disappearing. It allows drawing of data points anywhere on the plot, including in the plot margins.
                             ylim=c(0,max(layout$height)+0.1), # set y-limits. For some reason ggraph get's this completely wrong on its own...
                             expand=F # don't add extra 'expansion'.
    )
    ### ensure that annotation are placed some distance away from the y-axis.
    n_leaf_nodes <- length(order.dendrogram(dend)) # stats::order.dendrogram: A vector with length equal to the number of leaves in the dendrogram is returned
    p <- p + xlim(c(-1,n_leaf_nodes+1)) 
    # p <- p + scale_x_continuous(expand=expand_scale(add = 2, mult = c(4, .1)))
    # p <- p + scale_x_discrete(expand=expand_scale(add = 2, mult = c(4, .1)))
    # ^ ggraph does not respond do expand argument for scale_x_*() ---> MAYBE THIS IS BECAUSE coord_cartesian(expand=F) is set. 
    # TODO: try setting: coord_cartesian(expand=T) + scale_y_continuous(expand=0).
    # ^	Vector of range expansion constants used to add some padding around the data, to ensure that they are placed some distance away from the axes.  
    # Default: c(0, 0.6) for discrete variables.
    # Default: c(0.05, 0) for continuous variables
    # ^ REF: https://stackoverflow.com/a/52922323/6639640
    # ^ REF: https://ggplot2.tidyverse.org/reference/scale_discrete.html
    p <- p + theme_classic(base_family="Helvetica")
    p <- p + theme(axis.line.x=element_blank(), # remove all x-axis 
                   axis.title.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.text.x=element_blank())
    p  <- p  + theme(plot.margin = unit(c(1,1,4,1), "cm")) # (t, r, b, l) widen bottom margin
  }
  ### Add points (after drawing edges)
  p <- p + geom_node_point(aes(filter=((flag_prioritized==TRUE) & (leaf==TRUE)), color=label), size=1) # color prioritized annotation leaf nodes points
  ### Guides
  p <- p + 
    ### Color nodes points by main figure colors
    ### Apparently, geom_node_text() and geom_node_point() use the same scale_color_*, so the argument to values=X need to contain color mapping for both coloring
    scale_color_manual(values=colormap.nodes, guide=FALSE) + 
    guides(node_color="none",
           edge_color="none")
  ### Potential errors:
  ### ^ without base_family = 'Helvetica' REF: FIX ggsave() problem: https://github.com/thomasp85/ggraph/issues/152
  # Error in grid.Call.graphics(C_text, as.graphicsAnnot(x$label), x$x, x$y,  : 
  # invalid font type
  # In addition: There were 50 or more warnings (use warnings() to see the first 50)
  return(p)
}


### NON-CIRCULAR - BEFORE geom_edge_elbow2(aes(y=node.height,...)
#   p <- p + 
#     geom_node_text(aes(filter=((leaf==TRUE) & (flag_prioritized!=TRUE)), # leaf node labels
#                             color=!!sym_var.node_color, label=label),
#                     size=rel(0.8), 
#                     angle=90, hjust=1, nudge_y=-0.2, 
#                     show.legend=F) +
#     geom_node_text(aes(filter=(flag_prioritized==TRUE), # leaf node labels, prioritized cell-types
#                            color=!!sym_var.node_color, label=label),
#                    size=rel(1), 
#                    angle=90, hjust=1, nudge_y=-0.2, 
#                    show.legend=F)
#   p <- p + coord_cartesian(
#     clip = 'off' # This keeps the labels from disappearing. It allows drawing of data points anywhere on the plot, including in the plot margins.
#   )
# }
  

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


