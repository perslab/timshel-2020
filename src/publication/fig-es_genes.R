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

setwd(here("src/publication"))

# ======================================================================= #
# ================================ MOUSEBRAIN PLOT ================================ #
# ======================================================================= #

load(here("src/datasets-expression/mousebrain/mousebrain.sem_obj.RData"))

annotations_highlight <- get_prioritized_annotations_bmi(dataset="mousebrain")


# ======================================================================= #
# ============================= MAIN FIG ================================ #
# ======================================================================= #

### ESmu
genes_select <- c("POMC", "AGRP", "LEPR", "MC4R")
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, show_only_nonzero_es=show_only_nonzero_es)
p

### Mean expression
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="mean")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, 
                                           show_only_nonzero_es=show_only_nonzero_es, 
                                           scale.es_mu=F,
                                           scale.log=F)
p

# ======================================================================= #
# ================================= SOM FIG ============================= #
# ======================================================================= #

# ======================= BMI gene lists ======================= #
### TODO
# 1) PLOT MENDELIAN + RARE VARIANT GENES
# 2) PLOT LOCKE2015 Genes

### COOKBOOK
# load gene lists
# map from Ensembl to gene symbol
# call plot functions

# ======================= Standard model genes ======================= #
genes_select <- c("Npy", "Agrp", "Pomc", "Cartpt", "Lepr", "Insr", "Trh", "Cck", "Glp1r") 
# genes_select <- c("Gdf15", "Gfral", "Glp1r", "Hdac5", "Fam46a") 
genes_select <- toupper(genes_select)
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, show_only_nonzero_es=show_only_nonzero_es)
p


# ======================= MCXR ======================= #
### MC1R-MC5R (MC3R)
genes_select <- c("MC1R", "MC2R", "MC3R", "MC4R", "MC5R")
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, show_only_nonzero_es=show_only_nonzero_es)
p


# ======================= DRDX ======================= #
genes_select <- c("DRD2", "DRD1")
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, show_only_nonzero_es=show_only_nonzero_es)
p


