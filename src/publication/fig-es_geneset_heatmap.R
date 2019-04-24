############### SYNOPSIS ###################
### AIM: ES geneset heatmap

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

setwd(here("src/publication"))


# ======================================================================= #
# ==================================== TODO ============================= #
# ======================================================================= #
# TODO: implement dotplot. 
# dim = genes x selected cell-types
# size = ESmu
# color = blue

### TODO: BMI gene lists
#  PLOT LOCKE2015 NEAREST Genes



# ======================================================================= #
# ================================ LOAD ES DATA ============================ #
# ======================================================================= #

load(here("out/es/mousebrain_all.es_obj.RData"))
# sem_obj.mb <- sem_obj
# load(here("out/es/tabula_muris.es_obj.RData"))
# sem_obj.tm <- sem_obj

# ======================================================================= #
# ============================= LOAD GENE LISTS ========================= #
# ======================================================================= #

file.geneset <- here("data/genes_obesity/combined_gene_list.rare_and_mendelian_obesity_genes.txt")
df.geneset <- read_tsv(file.geneset)
df.geneset

# ======================================================================= #
# ================================ GET ES DF ============================ #
# ======================================================================= #

filter.genes <- df.geneset %>% pull(gene_symbol) %>% toupper()
print(length(filter.genes)) # 50

### Get ESmu
df.es <- get_es.gene_centric.single_es_metric(sem_obj, genes_select=filter.genes, es_metric="es_mu")
### "The following genes could not be found: BBS1;BBS5;ZBTB7B"
# gene_name annotation es_weight
# 1 ALMS1     ABC            0    
# 2 ARL6      ABC            0 

### Filter
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
df.es <- df.es %>% filter(annotation %in% filter.annotations)


# ======================================================================= #
# ================================ CLUSTER GENES ============================ #
# ======================================================================= #

### Convert to wide/'matrix' format
df.es.spread <- df.es %>% spread(key="annotation", value="es_weight") # genes x annotations
### Compute distances and hierarchical clustering
dd <- dist(df.es.spread %>% as.data.frame() %>% column_to_rownames(var="gene_name"), method = "euclidean") # dist() calcultes distances between rows
hc <- hclust(dd, method = "ward.D2")
# summary(hc)


# ======================================================================= #
# ================================ PLOT ============================ #
# ======================================================================= #
df.plot <- df.es

### Reorder gene_symbol by clustering
df.plot <- df.plot %>% mutate(gene_name = factor(gene_name, levels=hc$labels[hc$order]))

p <- ggplot(df.plot, aes(x=annotation, y=gene_name, fill=es_weight))
p <- p + geom_tile()
p <- p + labs(x="", y="", fill=expression(ES[mu]))
p <- p + colorspace::scale_fill_continuous_sequential(palette="Greens 2", rev=TRUE, na.value = "white", limits=c(0,1)) # "Blues 2","Blues 3","Purples 3"
p <- p + theme_classic()
p <- p + theme(axis.text.x=element_text(angle=45, hjust=1))
p

file.out <- sprintf("figs/fig_es.heatmap.mendelian_rare_genes.pdf")
ggsave(p, filename=file.out, width=5, height=8)


# ======================================================================= #
# ================================ XXXXXXXXX ============================ #
# ======================================================================= #


