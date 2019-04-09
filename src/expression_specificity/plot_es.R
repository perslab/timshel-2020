############### SYNOPSIS ###################
### AIM:
# create gene expression specificity plots

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

# library(plotly)

library(gghighlight)
library(ggrepel)
library(viridis)

source(here("src/lib/load_functions.R")) # load sc-genetics library


setwd(here("src/expression_specificity"))

# ======================================================================= #
# ================================ FUNCTIONS ================================ #
# ======================================================================= #

source(here("src/lib/plot_expression_specificity.R"))


# ======================================================================= #
# ================================ MOUSEBRAIN PLOT ================================ #
# ======================================================================= #

load(here("src/datasets-expression/mousebrain/mousebrain.sem_obj.RData"))

### Cell-types to plot
# AGRP ENSG00000159723 | DEINH6	Arcuate hypothalamic nucleus (probably, Npy, Agrp cells)
# DRD2 ENSG00000149295 | MSN2, MSN3, MSN4
# DRD1 ENSG00000184845 | MSN1, MSN5, MSN6


### MSN1
annotation <- "MSN2"
gene_highlight <- "DRD2"
df.es <- get_es.annotation_centric(sem_obj, annotation)
p <- plot_es.annotation_centric(df.es, genes_highlight=gene_highlight)
p + labs(title=annotation)
file.out <- sprintf("out.plot_expression_specificity_annotation_centric.%s.pdf", annotation)
# ggsave(file.out, width=10, height=6)


df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=gene_highlight, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=c("MSN2", "MSN3", "MSN4", "MSN1", "MSN5", "MSN6"))
p
# file.out <- sprintf("out.plot_expression_specificity_gene_centric.%s.pdf", gene_highlight)
# ggsave(file.out, width=12, height=6)


### DEINH6
annotation <- "DEINH6"
gene_highlight <- "AGRP"
df.es <- get_es.annotation_centric(sem_obj, annotation)
p <- plot_es.annotation_centric(df.es, genes_highlight=gene_highlight) # c("AGRP", "NPY")
p + labs(title=annotation)
file.out <- sprintf("out.plot_expression_specificity_annotation_centric.%s.pdf", annotation)
# ggsave(file.out, width=10, height=6)


df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=gene_highlight, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotation)
p
# file.out <- sprintf("out.plot_expression_specificity_gene_centric.%s.pdf", gene_highlight)
# ggsave(file.out, width=12, height=6)



# ======================================================================= #
# ================================ gene centric plots ================================ #
# ======================================================================= #

annotations_highlight <- get_prioritized_annotations_bmi(dataset="mousebrain")

### Main plot
genes_select <- c("POMC", "AGRP", "LEPR", "FTO", "BDNF", "MC4R", "BBS4", "DRD2", "DRD1")
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, show_only_nonzero_es=show_only_nonzero_es)
p


### PDYN / MC4R
# PVH-PDYN and PVH-MC4R neurons use separate efferent circuits to regulate appetite
# Li et al., 2019, Neuron (https://doi.org/10.1016/j.neuron.2019.02.028)
# "This expands the CNS satiety circuitry to include two non-overlapping PVH to hindbrain circuits."
genes_select <- c("PDYN", "MC4R", "SIM1")
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, show_only_nonzero_es=show_only_nonzero_es)
p

### TESTING other es_metric
genes_select <- c("PDYN", "MC4R", "SIM1")
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="tstat")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, 
                                           show_only_nonzero_es=show_only_nonzero_es, 
                                           scale.es_mu=F,
                                           scale.log=F)
p


### ADCY3 / MC4R
genes_select <- c("ADCY3", "MC4R")
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, show_only_nonzero_es=show_only_nonzero_es)
p


### Ciliopathy-related proteins
genes_select <- c("MRAP2","ALMS1","MKS1","CEP290","FRITZ","C2ORF86","SDCCAG8","CEP19","ANKRD26")
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, show_only_nonzero_es=show_only_nonzero_es)
p

genes_select <- c("BBS1","BBS2","BBS3","BBS4","BBS5","BBS6","BBS7","BBS8","BBS9","BBS10","BBS11","BBS12")
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, show_only_nonzero_es=show_only_nonzero_es)
p




### MC1R-MC5R (MC3R)
genes_select <- c("MC1R", "MC2R", "MC3R", "MC4R", "MC5R")
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, show_only_nonzero_es=show_only_nonzero_es)
p



### Other
annotations_highlight <- c("HYPEP5", "HYPEP2")
genes_select <- c("MC4R")
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, show_only_nonzero_es=show_only_nonzero_es)
p


genes_select <- c("SLC17A6")#VGLUT2 | excitatory
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, show_only_nonzero_es=show_only_nonzero_es)
p

### CCK
# "One mechanism proposed to explain how leptin re- duces food intake is by enhancing the CNS response to gut-derived satiety peptides, 
# such as cholecystokinin (CCK) and glucagon-like peptide-1 (GLP-1) (196) that are released upon food ingestion and activate 
# vagal afferent fibers that terminate in the NTS" (http://physrev.physiology.org/content/91/2/389.short p. 401ff)
### GDF15 seems to bind to its receptor (GFRAL) in the area postrema
### Hdac5 (which encodes histone deacetylase 5) is expressed in specific hypothalamic neurons that are activated by leptin.
### FAM46A is a trans regulator for leptin

genes_select <- c("Npy", "Agrp", "Pomc", "Cartpt", "Lepr", "Insr", "Trh", "Cck") 
genes_select <- c("Gdf15", "Gfral", "Glp1r", "Hdac5", "Fam46a") 
genes_select <- c("Pomc","Cartpt","Agrp","Anxa2","Gpr50","Sst","Th","Trh","Ttr","Glp1r","Oxt","Apoe") # long list from Arc-lira
genes_select <- c("Fos", "Jun", "Junb", "Egr1")

genes_select <- toupper(genes_select)
show_only_nonzero_es <- FALSE
df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=genes_select, es_metric="es_mu")
p <- plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotations_highlight, show_only_nonzero_es=show_only_nonzero_es)
p



# ======================================================================= #
# ================================ TABULA MURIS PLOT ================================ #
# ======================================================================= #
 
load(here("src/datasets-expression/tabula_muris/tabula_muris.sem_obj.RData"))

### Cell-types to plot
# - [INS NOT IN DATA - multiple orthologs] Beta-cells (insulin+). See http://www.ensembl.org/Homo_sapiens/Gene/Compara_Ortholog?db=core;g=ENSG00000254647;r=11:2181009-2182571
# - [OK] Alpha-cells (glucagon+)
# - [OK] Delta-cells (somatostatin+)
# - [OK] Gamma-cells (pancreatic polypeptide+)
# - [CELL-TYPE NOT IN DATA] Epsilon cells (ghrelin+)

### A
annotation <- "Pancreas.pancreatic A cell" # GCG
gene_highlight <- "GCG"
df.es <- get_es.annotation_centric(sem_obj, annotation)
p <- plot_es.annotation_centric(df.es, genes_highlight=gene_highlight)
p + labs(title="Pancreatic Alpha-cells (glucagon+)")
# file.out <- sprintf("out.plot_expression_specificity_annotation_centric.%s.pdf", str_replace_all(annotation, " ", "_"))
# ggsave(file.out, width=10, height=6)

df.es.gene <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=gene_highlight, es_metric="es_mu")
plot_es.gene_centric.single_es_metric(df.es.gene, annotations_highlight=annotation)
# file.out <- sprintf("out.plot_expression_specificity_gene_centric.%s.pdf", gene_highlight)
# ggsave(file.out, width=12, height=6)


### .....

### D
annotation <- "Pancreas.pancreatic D cell" # SST
df.es <- get_es.annotation_centric(sem_obj, annotation)
p <- plot_es.annotation_centric(df.es, genes_highlight="SST")
p + labs(title="Pancreatic Delta-cells (somatostatin+)")
# file.out <- sprintf("out.plot_expression_specificity.%s.pdf", str_replace_all(annotation, " ", "_"))
# ggsave(file.out, width=10, height=6)


### PP 
annotation <- "Pancreas.pancreatic PP cell" # PPY
df.es <- get_es.annotation_centric(sem_obj, annotation)
p <- plot_es.annotation_centric(df.es, genes_highlight="PPY")
p + labs(title="Pancreatic Gamma-cells (pancreatic polypeptide+)")
# file.out <- sprintf("out.plot_expression_specificity.%s.pdf", str_replace_all(annotation, " ", "_"))
# ggsave(file.out, width=10, height=6)




# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #
