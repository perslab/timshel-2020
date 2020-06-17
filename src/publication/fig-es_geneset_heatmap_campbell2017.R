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
file.enrich <- here("src/publication/tables/table-es_enrichment.combined.csv") # con
#file.enrich_geneset23 = here("out","es_enrichment_test","per_geneset","campbell2017_lvl2_BMI_23_genes_wilcoxon_genetestOuts_200613.csv.gz")
#file.enrich_geneset12 = here("out","es_enrichment_test","per_geneset","campbell2017_lvl2_mono_extr_12_genes_wilcoxon_genetestOuts_200613.csv.gz")

df.enrich = read_csv(file.enrich)
# df.enrich_geneset23 <- read_csv(file.enrich_geneset23)
# df.enrich_geneset12 <- read_csv(file.enrich_geneset12)

df.enrich <- df.enrich %>% filter(str_detect(annotation_fmt, "^ARCME-")) # *IMPORTANT* to filter to avoid join problems downstream

# ========================== CELL-TYPE MARKER GENES =================== #

df.arc_markergenes = get_metadata("hypothalamus") %>% select(annotation_fmt, annotation_marker)

# ========================== JOIN =================== #

df.tmp1 <- df.ldsc_cts %>% mutate(p.value.mlog10.ldsc = -log10(p.value)) %>% select(annotation, p.value.mlog10.ldsc)
df.tmp2 <- df.enrich %>% mutate(p.value.mlog10.enrich = -log10(p.value)) %>% select(annotation, p.value.mlog10.enrich, annotation_fmt)
#df.tmp2 <- df.enrich_geneset23 %>% mutate(p.value.mlog10.enrich = -log10(p.value)) %>% select(annotation, p.value.mlog10.enrich)
#df.tmp3 <- df.enrich_geneset12 %>% mutate(p.value.mlog10.enrich = -log10(p.value)) %>% select(annotation, p.value.mlog10.enrich)

# all(df.tmp1$annotation %in% df.tmp2$annotation)
# # [1] TRUE
# all(df.tmp2$annotation %in% df.tmp1$annotation)
# # [1] TRUE

df.pvals <- full_join(df.tmp1, df.tmp2, by="annotation")
#df.pvals <- df.pvals %>% mutate(annotation = utils.rename_annotations.hypothalamus(annotation, "campbell2017_lvl2", check_all_matches=T))

# df.pvals12 <- full_join(df.tmp1, df.tmp3, by="annotation")
# df.pvals12 <- df.pvals12 %>% mutate(annotation = utils.rename_annotations.hypothalamus(annotation, "campbell2017_lvl2", check_all_matches=T))
anyNA(df.pvals) # --> false
df.pvals <- df.pvals %>% mutate(annotation=annotation_fmt) %>% select(-annotation_fmt)
df.pvals
# anyNA(df.pvals) # --> FALSE
# anyNA(df.pvals12) # --> FALSE

#df.pvals12 <- df.pvals12[grep("NEURO",df.pvals12$annotation),]


#df.pvals <- df.pvals %>% mutate(annotation=annotation_fmt) %>% select(-annotation_fmt)
#df.pvals12 <- df.pvals12 %>% mutate(annotation=annotation_fmt) %>% select(-annotation_fmt)

# ======================================================================= #
# ================================ GENE LISTS =========================== #
# ======================================================================= #

file.geneset <- here("data/genes_obesity/prot_mono_early_extr_obesity_23genes.csv")
#file.geneset12 <- here("data/genes_obesity/mono_extr_obesity_12genes.csv")

df.geneset <- read_csv(file.geneset)
#df.geneset12 <- read_csv(file.geneset12)

filter.genes <- df.geneset %>% pull(gene_symbol) %>% toupper()
#filter.genes12 <- df.geneset12 %>% pull(gene_symbol) %>% toupper()
# filter.genes <- filter.genes[1:10]
#print(length(filter.genes)) # 50

# ======================================================================= #
# ================================ GET ES DF ============================ #
# ======================================================================= #
### Get ESmu
#df.es23 <- df.es %>% filter(gene %in% df.geneset23$ensembl_gene_ID) 
df.es <- get_es.gene_centric.single_es_metric(es_obj, genes_select=filter.genes, es_metric="es_mu")
#df.es12 <- df.es %>% filter(gene %in% df.geneset12$ensembl_gene_ID)
#df.es12 <- get_es.gene_centric.single_es_metric(es_obj, genes_select=filter.genes12, es_metric="es_mu")
# gene_name annotation es_weight
# 1 ALMS1     ABC            0    
# 2 ARL6      ABC            0 

### Filter
# filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
# df.es <- df.es %>% filter(annotation %in% filter.annotations)

### Change to clean names
df.es <- df.es %>% mutate(annotation = utils.rename_annotations.hypothalamus(annotation, "campbell2017_lvl2", check_all_matches=T))
#df.es12 <- df.es12 %>% mutate(annotation = utils.rename_annotations.hypothalamus(annotation, "campbell2017_lvl2", check_all_matches=T))

# ======================================================================= #
# ========================== RETRIEVE PROP EXPR TAB ===================#
# ======================================================================= #

# generated by here("src", "publication", "table-expression-stat.R")
# annotation x gene *(ensembl)

df.prop_expr = read_csv(here("out", "post_enrichment_analysis", "campbell_obesitygenes_prop_expr_200613.csv"))

# update cell type names
df.prop_expr <- df.prop_expr %>% mutate(annotation = utils.rename_annotations.hypothalamus(annotation, "campbell2017_lvl2", check_all_matches=T))

# transpose
mat_prop_expr = df.prop_expr[,2:ncol(df.prop_expr)] %>% as.matrix
rownames(mat_prop_expr) = df.prop_expr$annotation
mat_prop_expr = t(mat_prop_expr)
rownames(mat_prop_expr) = toupper(rownames(mat_prop_expr))

df.prop_expr = tibble("gene_name"=rownames(mat_prop_expr), as.tibble(mat_prop_expr)) 

df.prop_expr_long = df.prop_expr %>% pivot_longer(cols=colnames(df.prop_expr)[-1], names_to = "annotation", values_to = "prop_expr")
#df.prop_expr_long <- df.prop_expr_long %>% arrange(gene_name, annotation)

# ======================================================================= #
# =================================== PLOTS ============================== #
# ======================================================================= #

# ======================================================================= #
# ================ PLOT HEATMAP, ALL CELL TYPES ========================= #
# ======================================================================= #
#### 23 GENES ### 

# note that three of the 23 obesity genes are missing already from Campbell
# genes missing from df.prop_expr_long: "Lep"     "Zbtb7b"  "Rapgef3"

# NB two more genes were filtered out during CELLEX preprocessing 
# genes missing from df.plot.head: "Lep"     "Wnt10b"  "Zbtb7b"  "Rapgef3" "Znf169"

df.plot.heat <- df.es

vec_genes_filtered = setdiff(unique(df.prop_expr_long$gene_name), unique(df.plot.heat$gene_name))
vec_genes_filtered
#
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

# add marker genes
vec_markers1 = pull(df.arc_markergenes, annotation_marker)[match(df.plot.heat$annotation, df.arc_markergenes$annotation_fmt)]
df.plot.heat <- df.plot.heat %>% mutate(annotation_marker = paste0(annotation," (",vec_markers1, ")"))

vec_markers2 = pull(df.arc_markergenes, annotation_marker)[match(df.pvals$annotation, df.arc_markergenes$annotation_fmt)]
df.pvals <- df.pvals %>% mutate(annotation_marker = paste0(annotation," (",vec_markers2, ")"))

### ORDER annotations by pvals (enrich) [*REQUIRED* to match barplot order]
tmp.order <- df.pvals %>% arrange(p.value.mlog10.enrich) %>% pull(annotation_marker)
df.plot.heat <- df.plot.heat %>% mutate(annotation_marker = factor(annotation_marker, levels=tmp.order))

### CLUSTER GENES for better order [OPTIONAL] ### 
# Convert to wide/'matrix' format
df.es.spread <- df.es %>% spread(key="annotation", value="es_weight") # genes x annotations
# Compute distances and hierarchical clustering
dd <- dist(df.es.spread %>% as.data.frame() %>% column_to_rownames(var="gene_name"), method = "euclidean") # dist() calcultes distances between rows
hc <- hclust(dd, method = "ward.D2")
# summary(hc)
### Reorder gene_symbol by clustering
df.plot.heat <- df.plot.heat %>% mutate(gene_name = factor(gene_name, levels=hc$labels[hc$order]))
#df.plot.heat <- df.plot.heat %>% arrange(gene_name, annotation)


### 
p.heat <- ggplot(df.plot.heat, aes(x=annotation_marker, y=gene_name, fill=es_weight))
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
p.heat <- p.heat + theme(
  axis.text.x=element_text(angle=35, hjust=1), 
  axis.text.y=element_text(face="italic")
)
#p.heat <- p.heat + theme(axis.text.y=element_text(size=10))
#p.heat <- p.heat + coord_flip()
p.heat


# file.out <- sprintf("figs/fig_es.heatmap.campbell2017_lvl2.pdf")
# ggsave(p.heat, filename=file.out, width=10, height=4)

# ======================================================================= #
# ================================ BARPLOT ============================ #
# ======================================================================= #

df.plot.pvals <- df.pvals %>% 
  rename(Enrichment=p.value.mlog10.enrich, `S-LDSC`=p.value.mlog10.ldsc) %>%  # rename
  gather(key="method", value="value", -annotation, -annotation_marker)
df.plot.pvals


### ORDER annotations by pvals (enrich) [*REQUIRED* to match heatmap order]
tmp.order <- df.pvals %>% arrange(p.value.mlog10.enrich) %>% pull(annotation_marker)
df.plot.pvals <- df.plot.pvals %>% mutate(annotation_marker = factor(annotation_marker, levels=tmp.order))


p.bar <- ggplot(df.plot.pvals, aes(x=annotation_marker, y=value))
p.bar <- p.bar + geom_col(aes(fill=method), position=position_dodge2())
p.bar <- p.bar + labs(y=expression(-log[10](P)), x="", fill="")
#p.bar <- p.bar + coord_flip()

p.bar <- p.bar + theme_classic()
p.bar <- p.bar + theme(axis.text.x=element_blank())  #element_text(angle=45, hjust=1))
#p.bar <- p.bar + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank())
p.bar <- p.bar + scale_fill_manual(values=c("Enrichment"="#6baed6", "S-LDSC"="#02818a"))
p.bar <- p.bar + geom_hline(yintercept = -log10(0.05/nrow(df.pvals)), linetype="dashed", color="gray") 
p.bar

# ======================================================================= #
# ================================ PATCHWORK ============================ #
# ======================================================================= #

### Patch

#vertically
p.patch <- (p.bar+theme(legend.position="right")) / (p.heat+theme(legend.position="right")) + plot_layout(nrow=2)#, widths=c(1, 0.2))
p.patch <- p.patch + plot_annotation(theme = theme(plot.margin = margin(l=75,r=30)))

# p.patch <-(p.heat+theme(legend.position="top")) + 
#   (p.bar+theme(legend.position="top")) + 
#   plot_layout(nrow=1, widths=c(1, 0.2))
# p.patch <- p.patch + plot_annotation(theme = theme(plot.margin = margin(r=30)))

### Save
file.out1 <- sprintf("figs/fig_es.heatmap.campbell2017_lvl2_allcells.patch.pdf")
ggsave(plot=p.patch, filename=file.out1, width=18, height=8)

# ======================================================================= #
# ======================= PLOT HEATMAP, NEURONS ONLY =====================#
# ======================================================================= #
#### 23 GENES ### 

df.plot.heat_neur <- df.es

# keep only neurons, including n33,n34
df.prop_expr_long_neur <- df.prop_expr_long[grep("NEURO|OTHER1|OTHER2",df.prop_expr_long$annotation),]
df.pvals_neur <- df.pvals[grep("NEURO|OTHER1|OTHER2",df.pvals$annotation),]
df.plot.heat_neur <- df.plot.heat_neur[grep("NEURO|OTHER1|OTHER2", df.plot.heat_neur$annotation),]


df.prop_expr_long_neur <- df.prop_expr_long_neur %>% filter(!gene_name %in% vec_genes_filtered)

nrow(df.prop_expr_long_neur)==nrow(df.plot.heat_neur)
# [1] TRUE

# replace es values with NA for genes with low proportion of expression in cluster (>= 0.1)
df.plot.heat_neur$gene_name <- as.character(df.plot.heat_neur$gene_name)
df.plot.heat_neur$annotation <- as.character(df.plot.heat_neur$annotation)

df.plot.heat_neur <- df.plot.heat_neur %>% arrange(gene_name, annotation)
df.prop_expr_long_neur <- df.prop_expr_long_neur %>% arrange(gene_name, annotation)

all.equal(as.character(df.plot.heat_neur$gene_name), df.prop_expr_long_neur$gene_name)
all.equal(as.character(df.plot.heat_neur$annotation), df.prop_expr_long_neur$annotation)
# both TRUE

df.plot.heat_neur$es_weight <- ifelse(df.prop_expr_long_neur$prop_expr>0.1,df.plot.heat_neur$es_weight, NA_real_)

# add marker genes
vec_markers1_neur = pull(df.arc_markergenes, annotation_marker)[match(df.plot.heat_neur$annotation, df.arc_markergenes$annotation_fmt)]
df.plot.heat_neur <- df.plot.heat_neur %>% mutate(annotation_marker = paste0(annotation," (",vec_markers1_neur, ")"))

vec_markers2_neur = pull(df.arc_markergenes, annotation_marker)[match(df.pvals_neur$annotation, df.arc_markergenes$annotation_fmt)]
df.pvals_neur <- df.pvals_neur %>% mutate(annotation_marker = paste0(annotation," (",vec_markers2_neur, ")"))

### ORDER annotations by pvals (enrich) [*REQUIRED* to match barplot order]
tmp.order_neur <- df.pvals_neur %>% arrange(p.value.mlog10.enrich) %>% pull(annotation_marker)
df.plot.heat_neur <- df.plot.heat_neur %>% mutate(annotation_marker = factor(annotation_marker, levels=tmp.order_neur))

### CLUSTER GENES for better order [OPTIONAL] ### 
# Convert to wide/'matrix' format
df.es.spread_neur <- df.es %>% filter(grepl("NEURO|OTHER1|OTHER2",annotation)) %>% spread(key="annotation", value="es_weight") # genes x annotations
                                      
# Compute distances and hierarchical clustering
dd_neur <- dist(df.es.spread_neur %>% as.data.frame() %>% column_to_rownames(var="gene_name"), method = "euclidean") # dist() calcultes distances between rows
hc_neur <- hclust(dd_neur, method = "ward.D2")
# summary(hc)
### Reorder gene_symbol by clustering
df.plot.heat_neur <- df.plot.heat_neur %>% mutate(gene_name = factor(gene_name, levels=hc_neur$labels[hc_neur$order]))
#df.plot.heat_neur <- df.plot.heat_neur %>% arrange(gene_name, annotation)


### 
p.heat_neur <- ggplot(df.plot.heat_neur, aes(x=annotation_marker, y=gene_name, fill=es_weight))
p.heat_neur <- p.heat_neur + geom_tile()
p.heat_neur <- p.heat_neur + labs(x="", y="", fill=expression(ES[mu]))
#p.heat_neur <- p.heat_neur + colorspace::scale_fill_continuous_sequential(palette="Greens 2", rev=TRUE, na.value = "white", limits=c(0,1)) # "Blues 2","Blues 3","Purples 3"
# p.heat_neur <- p.heat_neur + scale_fill_binned() # requires ggplot2-3-3 REF: https://www.tidyverse.org/blog/2020/03/ggplot2-3-3-0/ + https://ggplot2.tidyverse.org/reference/scale_colour_continuous.html
p.heat_neur <- p.heat_neur + scale_fill_stepsn(colors=c("#edf8fb", "#b2e2e2", "#66c2a4", "#238b45"), 
                                     na.value = "white"
                                     # breaks=c(0, 0.25, 0.5, 0.75, 1),
                                     # limits=c(0,1)
                                     ) 
p.heat_neur <- p.heat_neur + theme_classic()
p.heat_neur <- p.heat_neur + theme(
  axis.text.x=element_text(angle=35, hjust=1), 
  axis.text.y=element_text(face="italic")
  )
#p.heat_neur <- p.heat_neur + theme(axis.text.y=element_text(size=10))
#p.heat_neur <- p.heat_neur + coord_flip()
p.heat_neur


# file.out <- sprintf("figs/fig_es.heatmap.campbell2017_lvl2.pdf")
# ggsave(p.heat_neur, filename=file.out, width=10, height=4)

# ======================================================================= #
# ========================= BARPLOT, NEURONS ONLY ======================= #
# ======================================================================= #

df.plot.pvals_neur <- df.pvals_neur %>% 
  rename(Enrichment=p.value.mlog10.enrich, `S-LDSC`=p.value.mlog10.ldsc) %>%  # rename
  gather(key="method", value="value", -annotation, -annotation_marker)
df.plot.pvals_neur


### ORDER annotations by pvals (enrich) [*REQUIRED* to match heatmap order]
tmp.order_neur <- df.pvals_neur %>% arrange(p.value.mlog10.enrich) %>% pull(annotation_marker)
df.plot.pvals_neur <- df.plot.pvals_neur %>% mutate(annotation_marker = factor(annotation_marker, levels=tmp.order_neur))


p.bar_neur <- ggplot(df.plot.pvals_neur, aes(x=annotation_marker, y=value))
p.bar_neur <- p.bar_neur + geom_col(aes(fill=method), position=position_dodge2())
p.bar_neur <- p.bar_neur + labs(y=expression(-log[10](P)), x="", fill="")
#p.bar_neur <- p.bar_neur + coord_flip()

p.bar_neur <- p.bar_neur + theme_classic()
p.bar_neur <- p.bar_neur + theme(axis.text.x=element_blank())  #element_text(angle=45, hjust=1))
#p.bar_neur <- p.bar_neur + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank())
p.bar_neur <- p.bar_neur + scale_fill_manual(values=c("Enrichment"="#6baed6", "S-LDSC"="#02818a"))
p.bar_neur <- p.bar_neur + geom_hline(yintercept = -log10(0.05/nrow(df.pvals_neur)), linetype="dashed", color="gray") 
p.bar_neur

# ======================================================================= #
# ========================PATCHWORK, NEURONS ============================ #
# ======================================================================= #

### Patch

#vertically
p.patch_neur <- (p.bar_neur+theme(legend.position="right")) / (p.heat_neur+theme(legend.position="right")) + plot_layout(nrow=2)#, widths=c(1, 0.2))
p.patch_neur <- p.patch_neur + plot_annotation(theme = theme(plot.margin = margin(l=75,r=30)))

# p.patch_neur <-(p.heat_neur+theme(legend.position="top")) + 
#   (p.bar_neur+theme(legend.position="top")) + 
#   plot_layout(nrow=1, widths=c(1, 0.2))
# p.patch_neur <- p.patch_neur + plot_annotation(theme = theme(plot.margin = margin(r=30)))

### Save
file.out2 <- sprintf("figs/fig_es.heatmap.campbell2017_lvl2_neur.patch.pdf")
ggsave(plot=p.patch_neur, filename=file.out2, width=11, height=8)

# ======================================================================= #
# #######==================== NEURONS, JUST 3 GENES  ==================== #
# ======================================================================= #

vec_3genes = c("LEPR", "POMC", "MC4R")
#vec_3genes = c("BDNF","CEP19","KSR2","LEP","LEPR","MC4R","NTRK2","PCSK1","POMC","SH2B1","SIM1","WNT10B")

df.plot.heat_neur <- df.plot.heat_neur %>% filter(gene_name %in% vec_3genes)
df.prop_expr_long_neur <- df.prop_expr_long_neur %>% filter(gene_name %in% vec_3genes)


### 
pheat3 <- ggplot(df.plot.heat_neur, aes(x=annotation_marker, y=gene_name, fill=es_weight))
pheat3 <- pheat3 + geom_tile()
pheat3 <- pheat3 + labs(x="", y="", fill=expression(ES[mu]))
#pheat3 <- pheat3 + scale_y_discrete(expand = c(0,0))
#pheat3 <- pheat3 + colorspace::scale_fill_continuous_sequential(palette="Greens 2", rev=TRUE, na.value = "white", limits=c(0,1)) # "Blues 2","Blues 3","Purples 3"
# pheat3 <- pheat3 + scale_fill_binned() # requires ggplot2-3-3 REF: https://www.tidyverse.org/blog/2020/03/ggplot2-3-3-0/ + https://ggplot2.tidyverse.org/reference/scale_colour_continuous.html
pheat3 <- pheat3 + scale_fill_stepsn(colors=c("#edf8fb", "#b2e2e2", "#66c2a4", "#238b45"), 
                                     na.value = "white"
                                     # breaks=c(0, 0.25, 0.5, 0.75, 1),
                                     # limits=c(0,1)
) 
pheat3 <- pheat3 + theme_classic()
pheat3 <- pheat3 + theme(
  axis.text.x=element_text(angle=35, hjust=1),
  axis.text.y=element_text(face="italic")
  )
#pheat3 <- pheat3 + theme(axis.text.y=element_text(size=10))
#pheat3 <- pheat3 + coord_flip()
pheat3

# ======================================================================= #
# ================================ PATCHWORK ============================ #
# ======================================================================= #

### Patch

p.patch3 <- (p.bar_neur+theme(legend.position="right")) / (pheat3+theme(legend.position="right")) + plot_layout(nrow=2)#, widths=c(1, 0.2))
p.patch3 <- p.patch3 + plot_annotation(theme = theme(plot.margin = margin(l=50,r=50)))
p.patch3 <- p.patch3 + plot_layout(heights = unit(c(4, 1.5), c('null')))
### Save
file.out3 <- sprintf("figs/fig_es.heatmap.campbell2017_lvl2_neur_3genes.patch.pdf")
ggsave(plot=p.patch3, filename=file.out3, width=11, height=4)
