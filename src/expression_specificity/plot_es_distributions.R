############### SYNOPSIS ###################
### AIM: Plot distributions of ES metrics.

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

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


setwd(here("src/expression_specificity"))

# ======================================================================= #
# ================================ FUNCTIONS ================================ #
# ======================================================================= #

source(here("src/lib/plot_expression_specificity.R"))

# ======================================================================= #
# ============================ PUBLICATION PLOT ========================= #
# ======================================================================= #

load("/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain/mousebrain.sem_obj.RData")
sem_obj.mb <- sem_obj
load("/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-maca/tabula_muris.sem_obj.RData")
sem_obj.tm <- sem_obj

# ============================ ANNO DENSITY PLOT ========================= #

df.es <- get_es_annotation_centric(sem_obj.mb, annotation="DEINH6")
df.es <- df.es %>% mutate(
  es_weight_pseudo_log=sign(es_weight)*log10(abs(es_weight+1)),
  es_weight_log=log10(es_weight)
)
# df.es <- df.es %>% filter(es_metric=="tstat")
# df.es <- df.es %>% filter(es_metric == "mean nES")
df.es <- df.es %>% filter(es_weight > 0)
p <- ggplot(df.es, aes(x=es_weight)) + 
  geom_density() + 
  scale_x_log10() +
  facet_wrap(~es_metric, scales = "free")
p


# ======================= TMP COLOR density ======================= #
### REF: https://stackoverflow.com/questions/3494593/shading-a-kernel-density-plot-between-two-points
set.seed(1)
draws <- rnorm(100)^2
q75 <- quantile(draws, .75)
q95 <- quantile(draws, .95)
dens <- density(draws)
dd <- with(dens,data.frame(x,y)) # 'dens' contains elements 'x' and 'y'
qplot(x,y,data=dd,geom="line")+
  geom_ribbon(data=subset(dd,x>q75 & x<q95),aes(ymax=y),ymin=0,
              fill="red",colour=NA,alpha=0.5)



# ======================================================================= #
# =================== GENE TILE PLOT | **DID NOT WORK WELL** ================ #
# ======================================================================= #


### gene - DEINH6, AGRP
df.es.gene <- get_es_gene_centric(sem_obj.mb, genes_select="AGRP")
df.es.gene
p <- ggplot(df.es.gene, aes(x=gene_idx_sorted_within_metric, y=es_metric)) + 
  geom_tile(color="white", fill="gray", width=.9, height=.9) +
  geom_tile(data=df.es.gene %>% filter(annotation=="DEINH6"), width=.9, height=.9, fill="black") # highlight
p

### gene - MSN2, DRD2
df.es.gene <- get_es_gene_centric(sem_obj.mb, genes_select="DRD2")
df.es.gene
p <- ggplot(df.es.gene, aes(x=gene_idx_sorted_within_metric, y=es_metric)) + 
  geom_tile(color="white", fill="gray", width=.9, height=.9) +
  geom_tile(data=df.es.gene %>% filter(annotation=="MSN2"), width=.9, height=.9, fill="black") # highlight
p


### gene - DEINH2, NPY
df.es.gene <- get_es_gene_centric(sem_obj.mb, genes_select="NPY")
df.es.gene
p <- ggplot(df.es.gene, aes(x=gene_idx_sorted_within_metric, y=es_metric)) + 
  geom_tile(color="white", fill="gray", width=.9, height=.9) #+
geom_tile(data=df.es.gene %>% filter(annotation=="DEINH6"), width=.9, height=.9, fill="black") # highlight
p



