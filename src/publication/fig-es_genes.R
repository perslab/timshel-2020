############### SYNOPSIS ###################
### AIM: plot expression specificity plots for selected annotations and genes

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

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

library(patchwork)

setwd(here("src/publication"))

# ======================================================================= #
# ============================== LOAD DATA ============================== #
# ======================================================================= #

load(here("out/es/mousebrain_all.es_obj.RData"))
sem_obj.mb <- sem_obj
load(here("out/es/tabula_muris.es_obj.RData"))
sem_obj.tm <- sem_obj

# ======================================================================= #
# ============================= MAIN FIG ================================ #
# ======================================================================= #


### SETTINGS
annotations_highlight <- get_prioritized_annotations_bmi(dataset="mousebrain")
colormap.annotation <- get_color_mapping.prioritized_annotations_bmi(dataset="mousebrain")
sem_obj <- sem_obj.mb

### MAIN PLOT: ESmu
genes_select <- c("POMC", "AGRP", "LEPR", "MC4R")
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, 
                                           annotations_highlight=annotations_highlight,
                                           annotations_colormap=colormap.annotation)
p <- p + labs(y=expression(ES[mu]))
p
file.out <- "figs/fig_es.main.mb.pomc_agrp_lepr_mc4r.pdf"
ggsave(plot=p, filename=file.out, width=8.5, height=5)


### SOM: Mean expression
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="expr_mean")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, 
                                           annotations_highlight=annotations_highlight, 
                                           annotations_colormap=colormap.annotation,
                                           scale.zero_to_one=F,
                                           scale.log=F)
p <- p + labs(y="log avg UMI (counts per 10,000)")
file.out <- "figs/fig_es.mb.pomc_agrp_lepr_mc4r.avg_expr.pdf"
ggsave(plot=p, filename=file.out, width=8.5, height=5)
p

# ======================================================================= #
# =================== SOM FIG: INDIVUDAL GENES *COMBINED* ================ #
# ======================================================================= #

### SETTINGS
annotations_highlight <- get_prioritized_annotations_bmi(dataset="mousebrain")
colormap.annotation <- get_color_mapping.prioritized_annotations_bmi(dataset="mousebrain")
sem_obj <- sem_obj.mb


### MOUSEBRAIN
genes_select <- c("MC1R", "MC2R", "MC3R", "MC4R", "MC5R", "ADCY3")
genes_select <- c(genes_select, "FTO", "IRX3", "IRX5")
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, 
                                             annotations_highlight=annotations_highlight,
                                             annotations_colormap=colormap.annotation)
p <- p + labs(y=expression(ES[mu]))
p
ggsave(filename="figs/fig_es.mb.mcxr_fto_genes.pdf", plot=p, width=8, height=8) # A4 8.3 x 11.7


genes_select <- c("DRD2", "DRD1", "DRD3")
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, 
                                           annotations_highlight=annotations_highlight,
                                           annotations_colormap=colormap.annotation)
p <- p + labs(y=expression(ES[mu]))
p
ggsave(filename="figs/fig_es.mb.drdx_genes.pdf", plot=p, width=8, height=8) # A4 8.3 x 11.7



# ======================= Standard model genes [appetite regulation] ======================= #

genes_select <- c("Npy", "Agrp", "Pomc", "Cartpt", "Lepr", "Insr", "Trh", "Cck", "Glp1r") 
genes_select <- c(genes_select, "Bdnf", "Cadm2", "Negr1") # known GWAS genes
genes_select <- c(genes_select, "Gdf15", "Gfral", "Hdac5", "Fam46a") # Hindbrain
genes_select <- c(genes_select, "Cckar", "Cckbr") # Cckar (Cck1), Cckbr (Cck2) # Schwartz thinks this is an important cell-type. ---> Cckbr is specifically expressed in "n29.Nr5a1/Bdnf" cell-type 
genes_select <- c(genes_select, "Sim1", "Pdyn") # Li et al., 2019, Neuron (https://doi.org/10.1016/j.neuron.2019.02.028)
genes_select <- toupper(genes_select)
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, 
                                           annotations_highlight=annotations_highlight,
                                           annotations_colormap=colormap.annotation)
p <- p + labs(y=expression(ES[mu]))
p
ggsave(filename="figs/fig_es.mb.standard_model_genes_and_more.pdf", plot=p, width=10, height=12) # A4 8.3 x 11.7


# ======================================================================= #
# =============== SOM: ES metrics comparisons, top_n cell-types ========== #
# ======================================================================= #



list.es_metrics <- list(
  "DET"="tstat", 
  "GES"="ges", 
  "NSI"="si", 
  "EP"="specificity",
  "ESmu"="es_mu", 
  "log avg UMI\n(counts per 10,000)"="expr_mean")

### MOUSEBRAIN
sem_obj <- sem_obj.mb
genes_select <- c("POMC", "AGRP")
list.res <- list()
for (i in seq(list.es_metrics)) {
  es_metric <- list.es_metrics[[i]]
  df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric=es_metric)
  p <- plot_es.gene_centric.single_es_metric(df.es.gene, 
                                             n_top_annotations_highlight=5, 
                                             scale.zero_to_one=F,
                                             size.highlight.text=rel(1.5),
                                             size.highlight.points=rel(1.5))
  p <- p + labs(y=names(list.es_metrics)[i])
  list.res[[es_metric]] <- p
}
p.patch <- patchwork::wrap_plots(list.res[-1], ncol=1) & theme(strip.text=element_blank())
p.patch <- list.res[[1]] + p.patch + plot_layout(ncol=1, heights=c(heights=c(1/7, 1))) & labs(x="")
ggsave(filename = "figs/fig_es.mb.es_metrics.pomc_agrp.pdf", plot=p.patch, width=8, height=8) # A4 8.3 x 11.7


### TABULA MURIS
sem_obj <- sem_obj.tm
genes_select <- c("GCG", "ADIPOQ")
list.res <- list()
for (i in seq(list.es_metrics)) {
  es_metric <- list.es_metrics[[i]]
  df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric=es_metric)
  p <- plot_es.gene_centric.single_es_metric(df.es.gene, 
                                             n_top_annotations_highlight=5, 
                                             scale.zero_to_one=F,
                                             size.highlight.text=rel(1.5),
                                             size.highlight.points=rel(1.5))
  p <- p + labs(y=names(list.es_metrics)[i])
  list.res[[es_metric]] <- p
}
p.patch <- patchwork::wrap_plots(list.res[-1], ncol=1) & theme(strip.text=element_blank())
p.patch <- list.res[[1]] + p.patch + plot_layout(ncol=1, heights=c(1/10, 1)) & labs(x="")
ggsave(filename = "figs/fig_es.tm.es_metrics.gcg_adipoq.pdf", plot=p.patch, width=8, height=8) # A4 8.3 x 11.7



### LEFTOVER: SIMPLE pathwork
# p.patch <- patchwork::wrap_plots(list.res, ncol=1) & 
#   theme(strip.text=element_blank()) &
#   labs(x="")



# ======================================================================= #
# =================== SOM FIG: INDIVUDAL GENES *SINGLE FIGS* ================ #
# ======================================================================= #

# ### SETTINGS
# annotations_highlight <- get_prioritized_annotations_bmi(dataset="mousebrain")
# colormap.annotation <- get_color_mapping.prioritized_annotations_bmi(dataset="mousebrain")
# sem_obj <- sem_obj.mb
# 
# # ======================= MCXR ======================= #
# ### MC1R-MC5R (MC3R)
# genes_select <- c("MC1R", "MC2R", "MC3R", "MC4R", "MC5R", "ADCY3")
# df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
# p <- plot_es.gene_centric.single_es_metric(df.es.gene, 
#                                            annotations_highlight=annotations_highlight,
#                                            annotations_colormap=colormap.annotation)
# p <- p + labs(y=expression(ES[mu]))
# file.out <- "figs/fig_es.som.mb_mcxr.pdf"
# ggsave(plot=p, filename=file.out, width=8.5, height=5)
# 
# 
# # ======================= DRDX ======================= #
# genes_select <- c("DRD2", "DRD1")
# df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
# p <- plot_es.gene_centric.single_es_metric(df.es.gene, 
#                                            annotations_highlight=annotations_highlight,
#                                            annotations_colormap=colormap.annotation)
# p <- p + labs(y=expression(ES[mu]))
# file.out <- "figs/fig_es.som.mb_drd.pdf"
# ggsave(plot=p, filename=file.out, width=7, height=4)
# 
# 
# # ======================= FTO ======================= #
# genes_select <- c("FTO", "IRX3", "IRX5")
# df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
# p <- plot_es.gene_centric.single_es_metric(df.es.gene, 
#                                            annotations_highlight=annotations_highlight,
#                                            annotations_colormap=colormap.annotation)
# p <- p + labs(y=expression(ES[mu]))
# file.out <- "figs/fig_es.som.mb_fto.pdf"
# ggsave(plot=p, filename=file.out, width=7, height=4)



# # ======================================================================= #
# # =================== SOM FIG: INDIVUDAL GENES *COMBINED* ================ #
# # ======================================================================= #
# ----> MAKES NO SENSE TO DO THIS INSTEAD OF FACET_WRAP

# ### SETTINGS
# annotations_highlight <- get_prioritized_annotations_bmi(dataset="mousebrain")
# colormap.annotation <- get_color_mapping.prioritized_annotations_bmi(dataset="mousebrain")
# sem_obj <- sem_obj.mb
# 
# source(here("src/lib/load_functions.R")) # load sc-genetics library
# 
# ### MOUSEBRAIN
# genes_select <- c("MC1R", "MC2R", "MC3R", "MC4R", "MC5R", "ADCY3")
# genes_select <- c(genes_select, "DRD2", "DRD1")
# genes_select <- c(genes_select, "FTO", "IRX3", "IRX5")
# list.genes <- as.list(genes_select)
# list.res <- list()
# for (i in seq(list.genes)) {
#   gene_now <- list.genes[[i]]
#   df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=gene_now, es_metric="es_mu")
#   if (nrow(df.es.gene) == 0) {
#     print("Gene not found..")
#     next
#   }
#   p <- plot_es.gene_centric.single_es_metric(df.es.gene, 
#                                              annotations_highlight=annotations_highlight,
#                                              annotations_colormap=colormap.annotation,
#                                              scale.zero_to_one=F,
#                                              size.highlight.text=rel(1.5),
#                                              size.highlight.points=rel(1.5))
#   list.res[[gene_now]] <- p
# }
# p.patch <- patchwork::wrap_plots(list.res) & labs(x="", y="")
# ggsave(filename="figs/fig_es.som.es_metrics.mb_ALL_genes.pdf", plot=p.patch, width=8, height=8) # A4 8.3 x 11.7

