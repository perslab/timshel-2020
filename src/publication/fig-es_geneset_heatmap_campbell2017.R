############### SYNOPSIS ###################
### AIM: ES geneset heatmap | campbell2017_lvl2

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:


# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

.libPaths(c(.libPaths(),"~/R/x86_64-pc-linux-gnu-library/"))
library("vctrs", lib.loc = "~/R/x86_64-pc-linux-gnu-library/")
library("rlang", lib.loc = "~/R/x86_64-pc-linux-gnu-library/")
library(tidyverse)

library(here)

library(patchwork)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))

# ======================================================================= #
# ================================ LOAD ES DATA ============================ #
# ======================================================================= #

load(here("out/es/campbell2017_lvl2.es_obj.RData")) # es_obj
#df.es = read_csv(here("out","es","campbell2017_lvl2.mu.csv.gz"))
# ======================================================================= #
# ================================ LOAD PVAL DATA ============================ #
# ======================================================================= #

# ================== LOAD LDSC RESULTS (multi GWAS) ================ #
file.results <- here("results/cellect_ldsc/prioritization.csv")
df.ldsc_cts <- read_csv(file.results)
df.ldsc_cts <- format_cellect_ldsc_results(df.ldsc_cts)
df.ldsc_cts <- df.ldsc_cts %>% filter(specificity_id == "campbell2017_lvl2", gwas =="BMI_UKBB_Loh2018")

# ========================== CELL-TYPE BMI ENRICHMENT =================== #
# ### READ: enrichment results for MB+hypothalamus cell-types
#file.enrich <- here("src/publication/tables/table-es_enrichment.combined.csv") # con
file.enrich_geneset23 = here("out","es_enrichment_test","per_geneset","campbell2017_lvl2_BMI_23_genes_wilcoxon_genetestOuts_200613.csv.gz")
file.enrich_geneset12 = here("out","es_enrichment_test","per_geneset","campbell2017_lvl2_mono_extr_12_genes_wilcoxon_genetestOuts_200613.csv.gz")

df.enrich_geneset23 <- read_csv(file.enrich_geneset23)
df.enrich_geneset12 <- read_csv(file.enrich_geneset12)

#df.enrich <- df.enrich %>% filter(str_detect(annotation_fmt, "^ARCME-")) # *IMPORTANT* to filter to avoid join problems downstream

# ========================== JOIN =================== #

df.tmp1 <- df.ldsc_cts %>% mutate(p.value.mlog10.ldsc = -log10(p.value)) %>% select(annotation, p.value.mlog10.ldsc)
df.tmp2 <- df.enrich_geneset23 %>% mutate(p.value.mlog10.enrich = -log10(p.value)) %>% select(annotation, p.value.mlog10.enrich)
df.tmp3 <- df.enrich_geneset12 %>% mutate(p.value.mlog10.enrich = -log10(p.value)) %>% select(annotation, p.value.mlog10.enrich)

all(df.tmp1$annotation %in% df.tmp2$annotation)
# [1] TRUE
all(df.tmp2$annotation %in% df.tmp1$annotation)
# [1] TRUE

df.pvals23 <- full_join(df.tmp1, df.tmp2, by="annotation")
df.pvals23 <- df.pvals23 %>% mutate(annotation = utils.rename_annotations.hypothalamus(annotation, "campbell2017_lvl2", check_all_matches=T))

df.pvals12 <- full_join(df.tmp1, df.tmp3, by="annotation")
df.pvals12 <- df.pvals12 %>% mutate(annotation = utils.rename_annotations.hypothalamus(annotation, "campbell2017_lvl2", check_all_matches=T))

anyNA(df.pvals23) # --> FALSE
anyNA(df.pvals12) # --> FALSE

# only keep neurons
df.pvals23 <- df.pvals23[grep("NEURO",df.pvals23$annotation),]
df.pvals12 <- df.pvals12[grep("NEURO",df.pvals12$annotation),]


#df.pvals23 <- df.pvals23 %>% mutate(annotation=annotation_fmt) %>% select(-annotation_fmt)
#df.pvals12 <- df.pvals12 %>% mutate(annotation=annotation_fmt) %>% select(-annotation_fmt)

# ======================================================================= #
# ================================ GENE LISTS =========================== #
# ======================================================================= #

file.geneset23 <- here("data/genes_obesity/prot_mono_early_extr_obesity_23genes.csv")
file.geneset12 <- here("data/genes_obesity/mono_extr_obesity_12genes.csv")

df.geneset23 <- read_csv(file.geneset23)
df.geneset12 <- read_csv(file.geneset12)

filter.genes23 <- df.geneset23 %>% pull(gene_symbol) %>% toupper()
filter.genes12 <- df.geneset12 %>% pull(gene_symbol) %>% toupper()
# filter.genes <- filter.genes[1:10]
#print(length(filter.genes)) # 50

# ======================================================================= #
# ================================ GET ES DF ============================ #
# ======================================================================= #



### Get ESmu
#df.es23 <- df.es %>% filter(gene %in% df.geneset23$ensembl_gene_ID) 
df.es23 <- get_es.gene_centric.single_es_metric(es_obj, genes_select=filter.genes23, es_metric="es_mu")
#df.es12 <- df.es %>% filter(gene %in% df.geneset12$ensembl_gene_ID)
df.es12 <- get_es.gene_centric.single_es_metric(es_obj, genes_select=filter.genes12, es_metric="es_mu")
# gene_name annotation es_weight
# 1 ALMS1     ABC            0    
# 2 ARL6      ABC            0 

### Filter
# filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
# df.es <- df.es %>% filter(annotation %in% filter.annotations)

### Change to clean names
df.es23 <- df.es23 %>% mutate(annotation = utils.rename_annotations.hypothalamus(annotation, "campbell2017_lvl2", check_all_matches=T))
df.es12 <- df.es12 %>% mutate(annotation = utils.rename_annotations.hypothalamus(annotation, "campbell2017_lvl2", check_all_matches=T))

# only keep neurons
df.es23 <- df.es23[grep("NEURO", df.es23$annotation),]
df.es12 <- df.es12[grep("NEURO", df.es12$annotation),]
# ======================================================================= #
# ========================== RETRIEVE PROP EXPR TAB ===================#
# ======================================================================= #

# generated by here("src", "publication", "table-expression-stat.R")
# annotation x gene *(ensembl)

df.prop_expr = read_csv(here("out", "qc_checks", "campbell_obesitygenes_prop_expr_200613.csv"))

# update cell type names
df.prop_expr <- df.prop_expr %>% mutate(annotation = utils.rename_annotations.hypothalamus(annotation, "campbell2017_lvl2", check_all_matches=T))

# transpose
mat_prop_expr = df.prop_expr[,2:ncol(df.prop_expr)] %>% as.matrix
rownames(mat_prop_expr) = df.prop_expr$annotation
mat_prop_expr = t(mat_prop_expr)
rownames(mat_prop_expr) = toupper(rownames(mat_prop_expr))

df.prop_expr = tibble("gene_name"=rownames(mat_prop_expr), as.tibble(mat_prop_expr)) 
# ======================================================================= #
# ============================= PLOT HEATMAP ==============================#
# ======================================================================= #
#### 23 GENES ### 

df.plot.heat <- df.es23

df.prop_expr_long = df.prop_expr %>% pivot_longer(cols=colnames(df.prop_expr)[-1], names_to = "annotation", values_to = "prop_expr")
#df.prop_expr_long <- df.prop_expr_long %>% arrange(gene_name, annotation)

# keep only neurons
df.prop_expr_long <- df.prop_expr_long[grep("NEURO",df.prop_expr_long$annotation),]

# add genes which were filtered out during CELLEX preprocessing to df.plot.heat, with 0 esmu values
# note that three of the 23 obesity genes are missing already from Campbell
# genes missing from df.prop_expr_long: "Lep"     "Zbtb7b"  "Rapgef3"
# genes missing from df.plot.head: "Lep"     "Wnt10b"  "Zbtb7b"  "Rapgef3" "Znf169"
vec_genes_filtered = setdiff(unique(df.prop_expr_long$gene_name), unique(df.plot.heat$gene_name))
vec_genes_filtered

# ignore the missing genes
# tbbl_add = tibble("gene_name"=unlist(lapply(vec_genes_filtered, rep, length(unique(df.prop_expr_long$annotation)))),
#                 "annotation"=rep(unique(df.prop_expr_long$annotation), 2),
#                 "es_weight"=0)
# df.plot.heat = bind_rows(df.plot.heat, tbbl_add)
df.prop_expr_long <- df.prop_expr_long %>% filter(!gene_name %in% vec_genes_filtered)

nrow(df.prop_expr_long)==nrow(df.plot.heat)
# [1] TRUE

# replace es values with NA for genes with low proportion of expression in cluster (>= 0.1)
df.plot.heat$gene_name <- as.character(df.plot.heat$gene_name)
df.plot.heat$annotation <- as.character(df.plot.heat$annotation)

df.plot.heat <- df.plot.heat %>% arrange(gene_name, annotation)
df.prop_expr_long <- df.prop_expr_long %>% arrange(gene_name, annotation)

all.equal(as.character(df.plot.heat$gene_name), df.prop_expr_long$gene_name)
all.equal(as.character(df.plot.heat$annotation), df.prop_expr_long$annotation)
# both TRUE

df.plot.heat$es_weight <- ifelse(df.prop_expr_long$prop_expr>0.1,df.plot.heat$es_weight, NA_real_)

### ORDER annotations by pvals (enrich) [*REQUIRED* to match barplot order]
tmp.order <- df.pvals23 %>% arrange(p.value.mlog10.enrich) %>% pull(annotation)
df.plot.heat <- df.plot.heat %>% mutate(annotation = factor(annotation, levels=tmp.order))

### CLUSTER GENES for better order [OPTIONAL] ### 
# Convert to wide/'matrix' format
df.es.spread <- df.es23 %>% spread(key="annotation", value="es_weight") # genes x annotations

df.es.spread <- df.es23 %>% spread(key="annotation", value="es_weight") # genes x annotations
# Compute distances and hierarchical clustering
dd <- dist(df.es.spread %>% as.data.frame() %>% column_to_rownames(var="gene_name"), method = "euclidean") # dist() calcultes distances between rows
hc <- hclust(dd, method = "ward.D2")
# summary(hc)
### Reorder gene_symbol by clustering
df.plot.heat <- df.plot.heat %>% mutate(gene_name = factor(gene_name, levels=hc$labels[hc$order]))
#df.plot.heat <- df.plot.heat %>% arrange(gene_name, annotation)
  
### 
p.heat <- ggplot(df.plot.heat, aes(x=annotation, y=gene_name, fill=es_weight))
p.heat <- p.heat + geom_tile()
p.heat <- p.heat + labs(x="", y="", fill=expression(ES[mu]))
#p.heat <- p.heat + colorspace::scale_fill_continuous_sequential(palette="Greens 2", rev=TRUE, na.value = "white", limits=c(0,1)) # "Blues 2","Blues 3","Purples 3"
# p.heat <- p.heat + scale_fill_binned() # requires ggplot2-3-3 REF: https://www.tidyverse.org/blog/2020/03/ggplot2-3-3-0/ + https://ggplot2.tidyverse.org/reference/scale_colour_continuous.html
p.heat <- p.heat + scale_fill_stepsn(colors=c("#edf8fb", "#b2e2e2", "#66c2a4", "#238b45"), 
                                     na.value = "white"
                                     # breaks=c(0, 0.25, 0.5, 0.75, 1),
                                     # limits=c(0,1)
                                     ) 
p.heat <- p.heat + theme_classic()
p.heat <- p.heat + theme(axis.text.x=element_text(angle=45, hjust=1))
#p.heat <- p.heat + theme(axis.text.y=element_text(size=10))
p.heat <- p.heat + coord_flip()
p.heat


# file.out <- sprintf("figs/fig_es.heatmap.campbell2017_lvl2.pdf")
# ggsave(p.heat, filename=file.out, width=10, height=4)


# ======================================================================= #
# ================================ BARPLOT ============================ #
# ======================================================================= #

df.plot.pvals <- df.pvals23 %>% 
  rename(Enrichment=p.value.mlog10.enrich, `S-LDSC`=p.value.mlog10.ldsc) %>%  # rename
  gather(key="method", value="value", -annotation)
df.plot.pvals


### ORDER annotations by pvals (enrich) [*REQUIRED* to match heatmap order]
tmp.order <- df.pvals23 %>% arrange(p.value.mlog10.enrich) %>% pull(annotation)
df.plot.pvals <- df.plot.pvals %>% mutate(annotation = factor(annotation, levels=tmp.order))


p.bar <- ggplot(df.plot.pvals, aes(x=annotation, y=value))
p.bar <- p.bar + geom_col(aes(fill=method), position=position_dodge2())
p.bar <- p.bar + coord_flip()
p.bar <- p.bar + labs(y=expression(-log[10](P)), x="", fill="")
p.bar <- p.bar + theme_classic()
p.bar <- p.bar + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank())
p.bar <- p.bar + scale_fill_manual(values=c("Enrichment"="#6baed6", "S-LDSC"="#02818a"))
p.bar

# ======================================================================= #
# ================================ PATCHWORK ============================ #
# ======================================================================= #

### Patch

p.patch <-(p.heat+theme(legend.position="top")) + 
  (p.bar+theme(legend.position="top")) + 
  plot_layout(nrow=1, widths=c(1, 0.2))
p.patch <- p.patch + plot_annotation(theme = theme(plot.margin = margin(r=30)))

### Save
file.out <- sprintf("figs/fig_es.heatmap.campbell2017_lvl2.patch.pdf")
ggsave(plot=p.patch, filename=file.out, width=10, height=8)

# ======================================================================= #
# ====================== version with just 12 genes  ==================== #
# ======================================================================= #

vec_12genes = c("BDNF","CEP19","KSR2","LEP","LEPR","MC4R","NTRK2","PCSK1","POMC","SH2B1","SIM1","WNT10B")

df.plot.heat <- df.plot.heat %>% filter(gene_name %in% vec_12genes)
df.prop_expr_long <- df.prop_expr_long %>% filter(gene_name %in% vec_12genes)


### 
pheat12 <- ggplot(df.plot.heat, aes(x=annotation, y=gene_name, fill=es_weight))
pheat12 <- pheat12 + geom_tile()
pheat12 <- pheat12 + labs(x="", y="", fill=expression(ES[mu]))
#pheat12 <- pheat12 + colorspace::scale_fill_continuous_sequential(palette="Greens 2", rev=TRUE, na.value = "white", limits=c(0,1)) # "Blues 2","Blues 3","Purples 3"
# pheat12 <- pheat12 + scale_fill_binned() # requires ggplot2-3-3 REF: https://www.tidyverse.org/blog/2020/03/ggplot2-3-3-0/ + https://ggplot2.tidyverse.org/reference/scale_colour_continuous.html
pheat12 <- pheat12 + scale_fill_stepsn(colors=c("#edf8fb", "#b2e2e2", "#66c2a4", "#238b45"), 
                                     na.value = "white"
                                     # breaks=c(0, 0.25, 0.5, 0.75, 1),
                                     # limits=c(0,1)
) 
pheat12 <- pheat12 + theme_classic()
pheat12 <- pheat12 + theme(axis.text.x=element_text(angle=45, hjust=1))
#pheat12 <- pheat12 + theme(axis.text.y=element_text(size=10))
pheat12 <- pheat12 + coord_flip()
pheat12


# file.out <- sprintf("figs/fig_es.heatmap.campbell2017_lvl2.pdf")
# ggsave(pheat12, filename=file.out, width=10, height=4)


# ======================================================================= #
# ================================ BARPLOT ============================ #
# ======================================================================= #

# df.plot.pvals <- df.pvals23 %>% 
#   rename(Enrichment=p.value.mlog10.enrich, `S-LDSC`=p.value.mlog10.ldsc) %>%  # rename
#   gather(key="method", value="value", -annotation)
# df.plot.pvals


### ORDER annotations by pvals (enrich) [*REQUIRED* to match heatmap order]
# tmp.order <- df.pvals23 %>% arrange(p.value.mlog10.enrich) %>% pull(annotation)
# df.plot.pvals <- df.plot.pvals %>% mutate(annotation = factor(annotation, levels=tmp.order))


# p.bar <- ggplot(df.plot.pvals, aes(x=annotation, y=value))
# p.bar <- p.bar + geom_col(aes(fill=method), position=position_dodge2())
# p.bar <- p.bar + coord_flip()
# p.bar <- p.bar + labs(y=expression(-log[10](P)), x="", fill="")
# p.bar <- p.bar + theme_classic()
# p.bar <- p.bar + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank())
# p.bar <- p.bar + scale_fill_manual(values=c("Enrichment"="#6baed6", "S-LDSC"="#02818a"))
# p.bar

# ======================================================================= #
# ================================ PATCHWORK ============================ #
# ======================================================================= #

### Patch

p.patch12 <-(pheat12+theme(legend.position="top")) + 
  (p.bar+theme(legend.position="top")) + 
  plot_layout(nrow=1, widths=c(1, 0.2))
p.patch12 <- p.patch12 + plot_annotation(theme = theme(plot.margin = margin(r=50)))

### Save
file.out <- sprintf("figs/fig_es.heatmap.campbell2017_lvl2.patch_12genes.pdf")
ggsave(plot=p.patch12, filename=file.out, width=6, height=8)
