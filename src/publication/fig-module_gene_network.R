############### SYNOPSIS ###################
### AIM: create a network gene module plot

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
library(ggraph)
library(tidygraph)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ============================== PARAMETERS ============================= #
# ======================================================================= #

n_max_genes <- 50

dirFiles = "/projects/jonatan/applied/18-mousebrain_7/tables/networkplot_args/"
prefixOut ="mb_7.2_190408"

celltype <- "MEINH2" 
#celltype <- "TEGLU23"

module <-  "lavenderblush" 
#module <- "lightpink3"

# ======================================================================= #
# ================================ LOAD DATA ================================ #
# ======================================================================= #

### Module information
file.module_geneset <- here("results/modules--metadata.txt")
df.module_geneset <- read_tsv(file.module_geneset)

### Gene expression
datExpr = read.csv(paste0(dirFiles, prefixOut,"_", module, "_", celltype, "_datExprScaled.csv.gz"), quote="", stringsAsFactors = F, row.names = 1)
# ^ cells x genes

### MAGMA
file.magma_gwas <- here("out/magma/gene_based_scores/BMI_UKBB_Loh2018_no_mhc.resid_correct_all.gsa.genes.mapped.out") # 100 KB window
df.magma <- read_tsv(file.magma_gwas)

# ======================================================================= #
# ========================== PROCESS MODULE DATA ======================== #
# ======================================================================= #

### Filter and rename
df.module <- df.module_geneset %>% 
  filter(annotation == module) %>%
  select(annotation, kme=annotation_value, gene_symbol=hgnc) # hgnc contain Jon's mouse gene symbols

# ======================================================================= #
# ======================= EXTRACT DATA - OLD using JT kMEs =============== #
# ======================================================================= #
# ---> problem is that Jon has not mapped module to human

# moduleColors = readRDS(paste0(dirFiles, prefixOut,"_", module, "_", celltype, "_moduleColors.RDS.gz"))
# ^ Named chr [1:5000]. Maps genes to modules

# kMEs =read.csv(paste0(dirFiles, prefixOut,"_", module, "_", celltype, "_kMEs.csv.gz"), quote="", stringsAsFactors = F, row.names=1)
# # ^ genes x modules
# 
# # Get the idx and kME values of top module genes ordered by decreasing kME
# gene_idx <- order(kMEs[,module], decreasing = T)[1:n_max_genes]
# names(gene_idx) <- rownames(kMEs)[gene_idx]
# # mod_kMEs = kMEs[gene_idx, module] # not used - legacy
# 
# ### Create df.kme
# df.kme <- kMEs %>% select(!!module) %>% 
#   rownames_to_column(var="gene_symbol") %>% 
#   mutate(gene_symbol=toupper(gene_symbol)) %>%
#   as_tibble()

# ======================================================================= #
# ================================= EXTRACT DATA ======================== #
# ======================================================================= #

genes_to_plot <- df.module %>% pull(gene_symbol)
genes_to_plot <- genes_to_plot[1:min(n_max_genes,length(genes_to_plot))]

# ======================================================================= #
# ================================ Correlation ================================ #
# ======================================================================= #

### Compute correlation
# df.corrr <- datExpr %>% select(genes_to_plot) %>% correlate(use = "pairwise.complete.obs", method="pearson")
df.corrr <- datExpr %>% as_tibble() %>% select(one_of(genes_to_plot)) %>% correlate(use = "pairwise.complete.obs", method="pearson")
# Unknown columns: `Pvrl3`
# ^ REF: dplyr::select() with some variables that may not exist in the data frame?: https://stackoverflow.com/a/51529352/6639640
df.corrr

### Convert to long format
df.corrr.long <- df.corrr %>%
  shave() %>% # # Convert the upper or lower triangle of a correlation data frame (cor_df) to missing values.
  stretch(na.rm = TRUE) # convert to long format. na.rm = TRUE: matrix diagonal (self-correlation) have NA values and will be dropped.
df.corrr.long
# # A tibble: 1,225 x 3
# x     y          r
# 1 Foxp2 Ndnf   0.363
# 2 Foxp2 Ntng1  0.481
# 3 Foxp2 Pde1c  0.322


# ======================================================================= #
# ========================= PREPARE NODE METADATA ======================== #
# ======================================================================= #

### Join MAGMA data with 
df.node_metadata <- df.magma %>% 
  mutate(gene_symbol = toupper(gene_symbol)) %>%
  filter(gene_symbol %in% toupper(genes_to_plot)) %>% # strictly not needed
  full_join(tibble(gene_symbol=toupper(genes_to_plot)), by="gene_symbol") %>% # *IMPORTANT*: add any genes that are not found in MAGMA. These genes will have NA vaules in ZSTAT
  select(gene_symbol, ZSTAT)

### Add kME data
df.node_metadata <- df.node_metadata %>% left_join(df.module %>% mutate(gene_symbol=toupper(gene_symbol)), by="gene_symbol")

### Convert gene symbols to uppercase
df.corrr.long <- df.corrr.long %>% mutate(x=toupper(x), y=toupper(y))
df.corrr.long
# ======================================================================= #
# ========================= ggraph + tidygraph ========================== #
# ======================================================================= #

### make tidygraph [positive correlations]
gr <- tbl_graph(nodes=df.node_metadata, edges=df.corrr.long %>% filter(r>=0), directed=FALSE)
gr


### plot (simple)
set.seed(1)
layout <- create_layout(gr, layout = 'fr', weights=gr %>% activate(edges) %>% pull(r)) # weigthed Fruchterman-Reingold layout
p <- ggraph(layout) + 
  geom_edge_link(aes(filter=r>=0, edge_width = abs(r), color = r)) +
  scale_edge_width_continuous(range = c(0, 1)) +
  scale_edge_colour_gradientn(limits = c(0, 1), colors = c("white", "dodgerblue2")) +
  # geom_node_circle(aes(r=kme/1.5, fill=ZSTAT), show.legend=T) +  
  # geom_node_text(aes(label=gene_symbol, size=kme^3), color="white", show.legend=F) + 
  geom_node_point(aes(size=kme, fill=ZSTAT), shape=21) +  
  geom_node_text(aes(label=gene_symbol, size=kme), repel=T, show.legend=F) + # IF USING geom_node_point()
  # & shapes 21-24 have both stroke colour and a fill. The size of the filled part is controlled by size, the size of the stroke is controlled by stroke. Each
  # labs(color="Cell-type origin", size=) +
  guides(
    edge_width="none", 
    edge_alpha="none",
    edge_color=guide_legend(title=expression(rho)),
    fill=guide_legend(title="MAGMA Z-Stat"),
    size = guide_legend(title="kME")
    ) +
  theme_graph(base_family="Helvetica") + # theme_graph() --> font family 'Arial Narrow' not found in PostScript font database
  theme(legend.position="right")
  # sans,serif,Helvetica,
# missing: Arial
p 
file.out <- "figs/fig_module_gene_networkplot.pdf"
ggsave(plot=p, filename=file.out, width=7, height=8)

