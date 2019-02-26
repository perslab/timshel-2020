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

library(gghighlight)
library(ggrepel)
library(viridis)

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


setwd(here("src/expression_specificity"))

# ======================================================================= #
# ================================ FUNCTIONS ================================ #
# ======================================================================= #


get_es_annotation_centric <- function(sem_obj, annotation) {
  if (is.null(sem_obj[["group_by_annotation.sem"]][[annotation]])) {
    stop(sprintf("Annotation = '%s' not found in data", annotation))
  }
  ### Get data
  df.es <- sem_obj[["group_by_annotation.sem"]][[annotation]] %>% 
    mutate(gene=sem_obj$genes) %>% 
    hs_add_gene_symbol_from_ensembl_ids(colname_geneids_from="gene", colname_geneids_to="gene_name") # map gene
  
  ### Add SEM mean
  df.es <- df.es %>% 
    mutate(
    `mean nES` = sem_obj$sem_meta$mean[[annotation]] # returns a (unnamed) numeric vector. This works even though sem_obj$sem_meta$mean is a tibble
    )
  ### Wrangle data
  df.es.mod <- df.es %>% 
    mutate(gene_idx_fixed = seq_len(n())) %>%
    gather(key="es_metric", value="es_weight", -gene, -gene_name, -gene_idx_fixed) %>%
    group_by(es_metric) %>%
    mutate(gene_idx_sorted_within_metric = rank(es_weight, ties.method = "first")) %>%
    ungroup()
  
  ### rename
  df.es.mod <- df.es.mod %>% mutate(es_metric = case_when(
    es_metric == "ges" ~ "GES",
    es_metric == "specificity" ~ "Specificity",
    es_metric == "si" ~ "SI",
    es_metric == "tstat" ~ "DE t-test",
    TRUE ~ as.character(es_metric))
  )
  return(df.es.mod)
}

es_plot_annotation_centric <- function(df.es, genes_highlight) {
  # df.es <- df.es %>% mutate(gene_label = paste0(gene_name, "(", round(es_weight,2), ")") ) # looks too complication
  df.es <- df.es %>% mutate(gene_label = gene_name)
  p <- ggplot(df.es, aes(x=gene_idx_fixed, y=es_weight)) + 
    geom_point(color="red") + 
    gghighlight(gene_name %in% genes_highlight, label_key=gene_label) + 
    facet_wrap(~es_metric, scales = "free") + 
    labs(x="Gene", y="ES weight") + 
    # theme_light() # this looks ok as well.
    theme_classic()
  return(p)
}


es_plot_gene_centric <- function(sem_obj, gene_select, annotation_highlight) {
  ### OBS: this function uses 'sem_meta$mean' but could use any data.
  ### Get data
  df.es_mean <- sem_obj$sem_meta$mean %>% 
    mutate(gene=sem_obj$genes) %>% 
    hs_add_gene_symbol_from_ensembl_ids(colname_geneids_from="gene", colname_geneids_to="gene_name") %>% # map gene
    select(gene_name, everything(), -gene)
  # df.es_mean %>% select(gene_name, everything(), -gene)
  
  df.es_mean.gene_filter <- df.es_mean %>% filter(gene_name == gene_select)
  df.es_mean.gene_filter <- df.es_mean.gene_filter %>% select(-gene_name) %>% gather(key="annotation", value=es_weight)
  ### Plot
  p <- ggplot(df.es_mean.gene_filter, aes(x=fct_reorder(annotation, es_weight), y=es_weight, label=annotation)) + 
    geom_point(size=3) + 
    geom_segment(aes(x=annotation, xend=annotation, y=0, yend=es_weight), color="grey") +
    # geom_label_repel(data=df.es_mean.gene_filter %>% filter(es_weight > 0)) + # OK, but easier to use gghighlight
    # gghighlight(es_weight > 0, label_key=annotation) +  # OK, but a but too much labels when the x-axis is ordered by es_weight (the labels clutter)
    # gghighlight(annotation %in% annotation_highlight, label_key=NULL) + # WORKS - no label
    gghighlight(annotation %in% annotation_highlight, label_key=annotation) + # WORKS - labels
    # geom_label(data=df.es_mean.gene_filter %>% filter(annotation %in% annotation_highlight)) + # give same results as gg_repel |  vjust="inward",hjust="inward"
    # gghighlight(annotation %in% annotation_highlight, label_key=annotation, label_params = list(hjust=0, vjust=0, direction="y")) + # cannot make 'label_params = list(hjust=0, vjust=0, direction="y")' to really do anything
    # gghighlight(annotation %in% annotation_highlight, label_key=annotation, label_params = list(nudge_x=1, nudge_y=1, position=NULL)) + # Error: Specify either `position` or `nudge_x`/`nudge_y`
    labs(x="Cell-type", y="mean nES weight", title=gene_select) + 
    # scale_color_viridis("plasma", direction=-1) +
    # coord_flip() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 75, hjust = 1))
  # theme(plot.margin = unit(c(1,3,1,1), "cm")) # t, r, b, l 
  # theme(axis.text.x=element_blank()) # hide x-axis labels
  return(p)
}


# p <- es_plot_gene_centric(sem_obj, gene_select=gene_highlight, annotation_highlight=annotation)
# p

# ======================================================================= #
# ================================ TABULA MURIS PLOT ================================ #
# ======================================================================= #
 
load("/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-maca/tabula_muris.sem_obj.RData")

### Cell-types to plot
# - [INS NOT IN DATA - multiple orthologs] Beta-cells (insulin+). See http://www.ensembl.org/Homo_sapiens/Gene/Compara_Ortholog?db=core;g=ENSG00000254647;r=11:2181009-2182571
# - [OK] Alpha-cells (glucagon+)
# - [OK] Delta-cells (somatostatin+)
# - [OK] Gamma-cells (pancreatic polypeptide+)
# - [CELL-TYPE NOT IN DATA] Epsilon cells (ghrelin+)

### A
annotation <- "Pancreas.pancreatic A cell" # GCG
gene_highlight <- "GCG"
df.es <- get_es_annotation_centric(sem_obj, annotation)
p <- es_plot_annotation_centric(df.es, genes_highlight=gene_highlight)
p + labs(title="Pancreatic Alpha-cells (glucagon+)")
file.out <- sprintf("out.plot_expression_specificity_annotation_centric.%s.pdf", str_replace_all(annotation, " ", "_"))
ggsave(file.out, width=10, height=6)

es_plot_gene_centric(sem_obj, gene_select=gene_highlight, annotation_highlight=annotation)
file.out <- sprintf("out.plot_expression_specificity_gene_centric.%s.pdf", gene_highlight)
ggsave(file.out, width=12, height=6)


### .....

### D
annotation <- "Pancreas.pancreatic D cell" # SST
df.es <- get_es_annotation_centric(sem_obj, annotation)
p <- es_plot_annotation_centric(df.es, genes_highlight="SST")
p + labs(title="Pancreatic Delta-cells (somatostatin+)")
file.out <- sprintf("out.plot_expression_specificity.%s.pdf", str_replace_all(annotation, " ", "_"))
ggsave(file.out, width=10, height=6)


### PP 
annotation <- "Pancreas.pancreatic PP cell" # PPY
df.es <- get_es_annotation_centric(sem_obj, annotation)
p <- es_plot_annotation_centric(df.es, genes_highlight="PPY")
p + labs(title="Pancreatic Gamma-cells (pancreatic polypeptide+)")
file.out <- sprintf("out.plot_expression_specificity.%s.pdf", str_replace_all(annotation, " ", "_"))
ggsave(file.out, width=10, height=6)



# ======================================================================= #
# ================================ MOUSEBRAIN PLOT ================================ #
# ======================================================================= #

load("/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain.sem_obj.RData")

### Cell-types to plot
# AGRP ENSG00000159723 | DEINH6	Arcuate hypothalamic nucleus (probably, Npy, Agrp cells)
# DRD2 ENSG00000149295 | MSN2, MSN3, MSN4
# DRD1 ENSG00000184845 | MSN1, MSN5, MSN6


### MSN1
annotation <- "MSN2"
gene_highlight <- "DRD2"
df.es <- get_es_annotation_centric(sem_obj, annotation)
p <- es_plot_annotation_centric(df.es, genes_highlight=gene_highlight)
p + labs(title=annotation)
file.out <- sprintf("out.plot_expression_specificity_annotation_centric.%s.pdf", annotation)
ggsave(file.out, width=10, height=6)


p <- es_plot_gene_centric(sem_obj, gene_select=gene_highlight, annotation_highlight=c("MSN2", "MSN3", "MSN4", "MSN1", "MSN5", "MSN6"))
p
file.out <- sprintf("out.plot_expression_specificity_gene_centric.%s.pdf", gene_highlight)
ggsave(file.out, width=12, height=6)


### DEINH6
annotation <- "DEINH6"
gene_highlight <- "AGRP"
df.es <- get_es_annotation_centric(sem_obj, annotation)
p <- es_plot_annotation_centric(df.es, genes_highlight=gene_highlight) # c("AGRP", "NPY")
p + labs(title=annotation)
file.out <- sprintf("out.plot_expression_specificity_annotation_centric.%s.pdf", annotation)
ggsave(file.out, width=10, height=6)


p <- es_plot_gene_centric(sem_obj, gene_select=gene_highlight, annotation_highlight=annotation)
p
file.out <- sprintf("out.plot_expression_specificity_gene_centric.%s.pdf", gene_highlight)
ggsave(file.out, width=12, height=6)




# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #

# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #



