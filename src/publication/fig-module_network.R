############### SYNOPSIS ###################
### AIM: create a network module plot

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(tidyverse)
library(here)

library(corrr)
library(ggraph) # 2019-02-19: CRAN only has ggraph_1.0.2 | use devtools::install_github('thomasp85/ggraph') to install ggraph v1.1 since it is build on tidygraph
# devtools::install_github('thomasp85/ggraph')
library(tidygraph)
library(igraph)



library(RColorBrewer)

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/wgcna_modules"))

# ======================================================================= #
# ================================ PARAMETERS ================================ #
# ======================================================================= #

gwas <- "BMI_UKBB_Loh2018"
genomic_annotation_prefix <- "wgcna.mousebrain-190213.fdr_sign_celltypes.continuous"


### Set LDSC specific files
# file.ldsc_cts <- sprintf("/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/%s__%s.cell_type_results.txt", genomic_annotation_prefix, gwas)
file.ldsc_cts <- here("results/prioritization_modules--mousebrain.BMI_UKBB_Loh2018.csv.gz")
# file.module_geneset <- sprintf("/scratch/sc-ldsc/%s/log.%s.multi_geneset.txt", genomic_annotation_prefix, genomic_annotation_prefix)
file.module_geneset <- here("results/modules--metadata.txt")

# ======================================================================= #
# ================================ READ DATA ================================ #
# ======================================================================= #

### LOAD kME
### this file correspons to RUN_ID="wgcna.mousebrain-190213.fdr_sign_celltypes.continuous"
# file.kme <- "/projects/jonatan/mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_kMs_full_join.csv.gz" 
# file.kme <- "/projects/jonatan/applied/18-mousebrain_7/tables/Neurons_sub_ClusterName_7.2_run1_kMs_full_join.csv.gz" # NEW APRIL 2019
file.kme <- here("results/modules-kme_table.csv.gz")
df.kme <- read_csv(file.kme)

### LOAD LDSC CTS RESULTS
df.ldsc_cts <- read_tsv(file.ldsc_cts)
df.ldsc_cts <- df.ldsc_cts %>% 
  mutate(module_id = str_split_fixed(Name, "__", n=Inf)[,2]) %>%
  select(module_id, p.value=Coefficient_P_value)
df.ldsc_cts


### Module origin metadata
file.module_origin_metadata <- here("data/expression/mousebrain/mousebrain-agg_L5.metadata.csv")
df.module_origin_metadata <- read_csv(file.module_origin_metadata) %>% select(module_origin=annotation, Region)
### Module origin mapping file (log.%.multi_geneset.txt)
df.module_geneset <- read_tsv(file.module_geneset)
df.module_geneset <- df.module_geneset %>% rename(ensembl_gene_id=gene, module_id=annotation, module_origin=cell_cluster)
df.module_geneset <- df.module_geneset %>% mutate(module_origin = stringr::str_replace_all(module_origin, pattern="\\s+", replacement="_"))  # *HACK FOR TABULA MURIS*
df.module_geneset <- df.module_geneset %>% rename(pkME=annotation_value) # *NEW*

# ======================================================================= #
# ========================== CREATE MODULE METADATA ===================== #
# ======================================================================= #

### Create the metadata we need/want
df.module_metadata <- df.module_geneset %>% group_by(module_origin, module_id) %>% summarise(n_genes=n()) # get distinct modules and count number of genes
df.module_metadata <- df.module_metadata %>% left_join(df.ldsc_cts, by="module_id") %>% mutate(mlog10p=-log10(p.value)) # add LDSC pvalue
df.module_metadata <- df.module_metadata %>% left_join(df.module_origin_metadata, by="module_origin") # add module origin metadata
head(df.module_metadata)

### Rename Region
df.module_metadata <- df.module_metadata %>% mutate(Region = case_when(
  Region == "Midbrain dorsal" ~ "Midbrain",
  Region == "Midbrain dorsal,Midbrain ventral" ~ "Midbrain",
  Region == "Hippocampus,Cortex" ~ "Hippocampus/Cortex",
  TRUE ~ as.character(Region))
)

### Make sure module_id is first column in matrix (need when used as 'node' data)
df.module_metadata <- df.module_metadata %>% select(module_id, everything())

# ======================================================================= #
# ================================ Correlation ================================ #
# ======================================================================= #

### Compute correlation
df.corrr <- df.kme %>% select(-genes) %>% correlate(use = "pairwise.complete.obs", method="pearson")

### Convert to long format
df.corrr.long <- df.corrr %>%
  shave() %>% # # Convert the upper or lower triangle of a correlation data frame (cor_df) to missing values.
  stretch(na.rm = TRUE) # convert to long format. na.rm = TRUE: matrix diagonal (self-correlation) have NA values and will be dropped.
# x             y                      r
# <chr>         <chr>              <dbl>
#   1 antiquewhite3 antiquewhite3    NA     
# 2 antiquewhite3 aquamarine2      -0.375 
# 3 antiquewhite3 blue1            -0.134 


# ======================================================================= #
# ================================ ggraph + tidygraph ================================ #
# ======================================================================= #

### make tidygraph
# gr <- tbl_graph(nodes=df.module_metadata, edges=df.corrr.long, directed=FALSE)
# gr <- tbl_graph(nodes=df.module_metadata, edges=df.corrr.long %>% filter(abs(r)>=0.5), directed=FALSE)
gr <- tbl_graph(nodes=df.module_metadata, edges=df.corrr.long %>% filter(r>=0), directed=FALSE)
gr


fdr_threshold <- 0.05/nrow(df.module_metadata)
fdr_threshold

### Get color map for Region
colormap.region <- get_color_mapping.mb.region()


### plot (simple)
set.seed(1)
layout <- create_layout(gr, layout = 'fr', weights=gr %>% activate(edges) %>% pull(r)) # weigthed Fruchterman-Reingold layout
# layout <- create_layout(gr, layout = 'fr') # un-weigthed FR
# layout <- create_layout(gr, layout = 'kk')
# layout <- create_layout(gr, layout = 'mds') 
# layout <- create_layout(gr, layout = 'auto') 
p <- ggraph(layout) + 
  geom_edge_link(aes(filter=r>=0.3, edge_alpha = abs(r), edge_width = abs(r), color = r)) +
  scale_edge_width_continuous(range = c(0, 1)) +
  scale_edge_colour_gradientn(limits = c(0.3, 1), colors = c("white", "dodgerblue2")) +
  # scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("firebrick2", "dodgerblue2")) + # positive and negative correlation
  geom_node_point(aes(size=mlog10p, fill=Region), shape=21) +  # size=pval. Try also mlog10p^3
  # & shapes 21-24 have both stroke colour and a fill. The size of the filled part is controlled by size, the size of the stroke is controlled by stroke. Each
  geom_node_point(aes(filter=mlog10p>-log10(fdr_threshold), size=mlog10p, fill=Region), stroke=2, color="red", shape=21, show.legend=F) + 
  # geom_node_point(aes(size=n_genes, color=Region)) +  # size=n_genes
  geom_node_text(aes(filter=mlog10p>-log10(fdr_threshold), label = module_id), repel = TRUE) + # REF: https://stackoverflow.com/a/50365886/6639640: tidygraph has filter property which can be used in various geoms for filtering nodes, edges, etc.
  # labs(color="Cell-type origin", size=) +
  scale_fill_manual(values=colormap.region) + 
  guides(
    edge_width="none", 
    edge_alpha="none",
    edge_color=guide_legend(title=expression(rho)),
    color=guide_legend(title="Module cell-type origin"),
    size = guide_legend(title=expression(-log[10](P[S-LDSC])))) +
  theme_graph(base_family="Helvetica") # theme_graph() --> font family 'Arial Narrow' not found in PostScript font database
# sans,serif,Helvetica,
# missing: Arial
p 
file.out <- "figs/fig_module_network_plot.pdf"
ggsave(plot=p, filename=file.out, width=10, height=8)


### Correlation graph REFs:
# https://drsimonj.svbtle.com/how-to-create-correlation-network-plots-with-corrr-and-ggraph
# http://rstudio-pubs-static.s3.amazonaws.com/278083_11ffb2dbd47a4ee1b68f29e41c53efce.html





# ======================================================================= #
# ==================== igraph [works, but not used] ===================== #
# ======================================================================= #

if (FALSE) { # Don't run this
  gr <- graph_from_data_frame(d = df.corrr.long, vertices = df.module_metadata, directed = FALSE)
  # ^ ### If vertices is NULL, then the first two columns of d are used as a symbolic edge list and additional columns as edge attributes. 
  # ^ The names of the attributes are taken from the names of the columns.
  # ^ ### If vertices is *not NULL*, then it must be a data frame giving vertex metadata. 
  # ^ The first column of vertices is assumed to contain symbolic vertex names, this will be added to the graphs as the ‘name’ vertex attribute
  # ^ Other columns will be added as additional vertex attributes. 
  # ^ If vertices is not NULL then the symbolic edge list given in d is checked to contain only vertex names listed in vertices.
  gr
  
  
  ### Map color | REF: https://www.r-graph-gallery.com/249-igraph-network-map-a-color/
  vertex_attr(gr)
  node_color_map <- brewer.pal(length(unique(V(gr)$Region)), "Set1") # Make a palette of 3 colors
  node_color <- node_color_map[as.numeric(as.factor(V(gr)$Region))] # Create a vector of color
  
  ### Plot
  l <- layout_with_fr(gr, weights = E(gr)$r) #  The weight edge attribute is used by default, if present.
  plot(gr, layout = l, asp = 1, edge.width = E(gr)$r*4, vertex.label =  "",
       vertex.color = adjustcolor(node_color, alpha.f = 0.9), 
       vertex.size = V(gr)$mlog10p * 12, vertex.frame.width = 2)
  
}

# ======================================================================= #
# ================================ LAYOUT DOCS ================================ #
# ======================================================================= #

# igraph_layouts <- c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl', 'lgl')

# layout_with_fr(graph, coords = NULL, dim = 2, niter = 500,
#                start.temp = sqrt(vcount(graph)), grid = c("auto", "grid", "nogrid"),
#                weights = NULL, minx = NULL, maxx = NULL, miny = NULL, maxy = NULL,
#                minz = NULL, maxz = NULL, coolexp, maxdelta, area, repulserad, maxiter)
# ---> weights is a vector.

### fr: Places nodes according to the force-directed algorithm of Fruchterman and Reingold. See with_fr
# Fruchterman-Reingold is one of the most used force-directed layout algorithms out there.
# Force-directed layouts try to get a nice-looking graph where edges are similar
# in length and cross each other as little as possible. They simulate the graph
# as a physical system. Nodes are electrically charged particles that repulse
# each other when they get too close. The edges act as springs that attract
# connected nodes closer together. As a result, nodes are evenly distributed
# through the chart area, and the layout is intuitive in that nodes which share
# more connections are closer to each other. The disadvantage of these
# algorithms is that they are rather slow and therefore less often used in
# graphs larger than ~1000 vertices. You can set the “weight” parameter which
# increases the attraction forces among nodes connected by heavier edges.
### kk: Uses the spring-based algorithm by Kamada and Kawai to place nodes. See with_kk
### nicely: Tries to pick an appropriate layout. See nicely for a description of the simpe decision tree it uses

