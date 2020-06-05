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
df.enrich <- read_csv(file.enrich)
df.enrich <- df.enrich %>% filter(str_detect(annotation_fmt, "^ARCME-")) # *IMPORTANT* to filter to avoid join problems downstream


# ========================== JOIN =================== #

df.tmp1 <- df.ldsc_cts %>% mutate(p.value.mlog10.ldsc = -log10(p.value)) %>% select(annotation, p.value.mlog10.ldsc)
df.tmp2 <- df.enrich %>% mutate(p.value.mlog10.enrich = -log10(p.value)) %>% select(annotation, p.value.mlog10.enrich, annotation_fmt)
df.pvals <- full_join(df.tmp1, df.tmp2, by="annotation")
anyNA(df.pvals) # --> false
df.pvals <- df.pvals %>% mutate(annotation=annotation_fmt) %>% select(-annotation_fmt)
df.pvals

# ======================================================================= #
# ================================ GENE LISTS =========================== #
# ======================================================================= #

### TMP PLACEHOLDER
file.geneset <- here("data/genes_obesity/combined_gene_list.rare_and_mendelian_obesity_genes.txt")
df.geneset <- read_tsv(file.geneset)
filter.genes <- df.geneset %>% pull(gene_symbol) %>% toupper()
# filter.genes <- filter.genes[1:10]
print(length(filter.genes)) # 50

# ======================================================================= #
# ================================ GET ES DF ============================ #
# ======================================================================= #

### Get ESmu
df.es <- get_es.gene_centric.single_es_metric(es_obj, genes_select=filter.genes, es_metric="es_mu")
# gene_name annotation es_weight
# 1 ALMS1     ABC            0    
# 2 ARL6      ABC            0 

### Filter
# filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
# df.es <- df.es %>% filter(annotation %in% filter.annotations)

### Change to clean names
df.es <- df.es %>% mutate(annotation = utils.rename_annotations.hypothalamus(annotation, "campbell2017_lvl2", check_all_matches=T))

# ======================================================================= #
# ================================ PLOT HEATMAP ============================ #
# ======================================================================= #

df.plot.heat <- df.es

### ORDER annotations by pvals (enrich) [*REQUIRED* to match barplot order]
tmp.order <- df.pvals %>% arrange(p.value.mlog10.enrich) %>% pull(annotation)
df.plot.heat <- df.plot.heat %>% mutate(annotation = factor(annotation, levels=tmp.order))

### CLUSTER GENES for better order [OPTIONAL] ### 
# Convert to wide/'matrix' format
df.es.spread <- df.es %>% spread(key="annotation", value="es_weight") # genes x annotations
# Compute distances and hierarchical clustering
dd <- dist(df.es.spread %>% as.data.frame() %>% column_to_rownames(var="gene_name"), method = "euclidean") # dist() calcultes distances between rows
hc <- hclust(dd, method = "ward.D2")
# summary(hc)
### Reorder gene_symbol by clustering
df.plot.heat <- df.plot.heat %>% mutate(gene_name = factor(gene_name, levels=hc$labels[hc$order]))

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
p.heat <- p.heat + coord_flip()
p.heat


file.out <- sprintf("figs/fig_es.heatmap.campbell2017_lvl2.pdf")
ggsave(p.heat, filename=file.out, width=10, height=4)


# ======================================================================= #
# ================================ BARPLOT ============================ #
# ======================================================================= #

df.plot.pvals <- df.pvals %>% 
  rename(Enrichment=p.value.mlog10.enrich, `S-LDSC`=p.value.mlog10.ldsc) %>%  # rename
  gather(key="method", value="value", -annotation)
df.plot.pvals


### ORDER annotations by pvals (enrich) [*REQUIRED* to match heatmap order]
tmp.order <- df.pvals %>% arrange(p.value.mlog10.enrich) %>% pull(annotation)
df.plot.pvals <- df.plot.pvals %>% mutate(annotation = factor(annotation, levels=tmp.order))


p.bar <- ggplot(df.plot.pvals, aes(x=annotation, y=value))
p.bar <- p.bar + geom_col(aes(fill=method), position=position_dodge2())
p.bar <- p.bar + coord_flip()
p.bar <- p.bar + labs(x="", y="", fill="")
p.bar <- p.bar + theme_classic()
p.bar <- p.bar + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank())
p.bar <- p.bar + scale_fill_manual(values=c("Enrichment"="#6baed6", "S-LDSC"="#02818a"))
p.bar

# ======================================================================= #
# ================================ PATCHWORK ============================ #
# ======================================================================= #

### Patch
p.patch <-(p.heat+theme(legend.position="top")) + (p.bar+theme(legend.position="top")) + plot_layout(nrow=1, widths=c(1, 0.2))
p.patch
# & theme(plot.margin = unit(c(3,1,1,1), "cm")) # (t, r, b, l) widen margin --> does not work for pathwork

### Save
# file.out <- sprintf("figs/fig_es.heatmap.campbell2017_lvl2.patch.pdf")
# ggsave(plot=p.patch, filename=file.out, width=10, height=5)

# ======================================================================= #
# ================================ XXXXXXXXX ============================ #
# ======================================================================= #


# ======================================================================= #
# ================================ XXXXXXXXX ============================ #
# ======================================================================= #

