############### SYNOPSIS ###################
### AIM: Make hypothalamus plot + BMI geneset cell-type enrichment
# Point jitter version

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

library(ggbeeswarm)

source(here("src/lib/load_functions.R")) # load sc-genetics library
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))


# ======================================================================= #
# ========================== BMI GENESET ENRICHMENT =================== #
# ======================================================================= #

### ---> NOT USED (but should work)
# ### READ: enrichment results for MB+hypothalamus cell-types
# file.enrich <- here("src/publication/tables/table-es_enrichment.combined.csv") # con
# df.enrich <- read_csv(file.enrich)
# df.enrich <- df.enrich %>% mutate(annotation_fmt=annotation) # copy column | only needed for legacy reasons because of _join() operation

# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
file.results <- here("results/cellect_ldsc/prioritization.csv")
df.ldsc_cts <- read_csv(file.results)
df.ldsc_cts <- format_cellect_ldsc_results(df.ldsc_cts)
df.ldsc_cts <- df.ldsc_cts %>% filter(specificity_id %in% get_scrna_seq_dataset_prefixes("hypo"))
df.ldsc_cts %>% count(specificity_id)

# ======================================================================= #
# ============================ BMI point plot =============================== #
# ======================================================================= #

### Create 'plotting' data frame
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == "BMI_UKBB_Loh2018")  # filter BMI results

### Add fdr_significant flag (within GWAS)
df.ldsc_cts <- df.ldsc_cts %>% mutate(fdr_significant = if_else(p.value <= 0.05/n(), true=T, false=F))

### Add columns
df.ldsc_cts <- df.ldsc_cts %>% mutate(p.value.mlog10 = -log10(p.value))

# ======================================================================= #
# ============================ PROCESS: add meta-data ====================== #
# ======================================================================= #

### Init
df.join <- df.ldsc_cts

### Add clean names
df.join <- df.join %>% mutate(annotation_fmt = utils.rename_annotations.hypothalamus(annotation, specificity_id, check_all_matches=T))

### Add meta-data
df.metadata <- get_metadata("hypothalamus")
df.join <- df.join %>% left_join(df.metadata %>% select(-annotation, -specificity_id), by="annotation_fmt")


### Add extra columns
df.join <- df.join %>% mutate(
  #category=taxonomy_lvl1, # *DETERMINES WHAT IS COLORED BY*
  category=annotation_prefix, # *DETERMINES WHAT IS COLORED BY*
  text=annotation_fmt
  # text=annotation_marker
  
)


# ======================================================================= #
# ============================ GET DENDROGRAM =========================== #
# ======================================================================= #

# source("fig-annotation_clustering_hypo.R")
# list.res <- get_hypothalamus_integrated_dendrogram(save_fig=FALSE)
# # p.dendro.orig <- list.res[["plot"]]
# dend <- list.res[["dendrogram"]]

# ======================================================================= #
# ============================ PLOT =============================== #
# ======================================================================= #

### Init
df.plot <- df.join

### Order annotations by tax1
# df.plot <- df.plot %>% arrange(taxonomy_lvl1, annotation_fmt) %>% mutate(annotation_fmt = factor(annotation_fmt, levels=annotation_fmt))
### Order annotations by annotation_prefix
df.plot <- df.plot %>% arrange(annotation_prefix, annotation_fmt) %>% mutate(annotation_fmt = factor(annotation_fmt, levels=annotation_fmt))


### Order annotations by dendrogram
# order.annotations <- labels(dend) # labels(dend) returns ordered labels | order.dendrogram(dend) returns ordered index.
# df.plot <- df.plot %>% mutate(annotation_fmt = factor(annotation_fmt, levels=order.annotations))

### Colormap
names.category <- sort(unique(df.plot$category))
# colormap.category <- brewer.pal(name="Set1", n=9) # n=9 max for Set1 | "Set2"/"Dark2" ok for colorblind
# colormap.category <- colormap.category[1:length(names.category)] # select the colors we need
colormap.category <- colorspace::qualitative_hcl(n=length(names.category), palette = "Dark 3")
names(colormap.category) <- names.category
colormap.category


fdr_threshold <- 0.05/nrow(df.plot)
p.priori <- ggplot() +
  ### fdr line
  geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color="darkgray") + 
  # p + geom_violin(show.legend=F)
  ggbeeswarm::geom_quasirandom(data=df.plot, aes(x=category, y=p.value.mlog10, fill=category, color=category), shape=21, stroke=0.1, show.legend=F) + # varwidth = TRUE
  # geom_text(aes(label=sample_id), position=ggbeeswarm::position_quasirandom(), size=1, show.legend=F, color="darkgray") # https://github.com/eclarke/ggbeeswarm/issues/35
  ### points - no enrichment
  #**SEE ME* geom_boxplot(data=df.plot, aes(x=category, y=p.value.mlog10, fill=category, color=category), shape=21, stroke=0.1) +
  ### points - w. enrichment
  # geom_point(data=df.plot, aes(x=annotation_fmt, y=p.value.mlog10, fill=category, size=-log10(combined_rare_mendelian_obesity), color=category), shape=21, stroke=0.1) +
  # ^ shapes 21-24 have both stroke colour and a fill. The size of the filled part is controlled by size, the size of the stroke is controlled by stroke.
  # geom_point(data=df.plot %>% filter(combined_rare_mendelian_obesity<fdr_threshold), aes(annotation_fmt, y=p.value.mlog10, size=-log10(combined_rare_mendelian_obesity), fill=category), stroke=0.8, color="black", shape=21, show.legend=F) + 
  ### text/description
  geom_text(data=df.plot %>% filter(fdr_significant), aes(x=category, y=p.value.mlog10, label=text), position=ggbeeswarm::position_quasirandom(), size=rel(2), hjust=0, show.legend=F) +  # https://github.com/eclarke/ggbeeswarm/issues/35
  # geom_text(data=df.plot, aes(x=annotation_fmt, y=p.value.mlog10, label=text), size=rel(2), angle=90, vjust=0.22, hjust=-0.08, show.legend=F) +
  # ^ this works when no using coord_flip()
  ### axes
  labs(x="", y=expression(-log[10](P[S-LDSC])), size=expression(-log[10](P[enrichment])), fill="") +
  # coord
  coord_flip(clip = 'off') + # This keeps the labels from disappearing 
  ### color
  scale_color_manual(values=colormap.category) +
  scale_fill_manual(values=colormap.category) +
  ### guides
  guides(color=F,
         fill=guide_legend(order = 1),
         size=guide_legend(order = 2)
  ) + 
  ### theme
  theme_classic() + 
  theme(axis.text.y=element_text(size=rel(1)),
        axis.ticks.y=element_blank()) +
  # theme(axis.text.x=element_text(angle=65, hjust=1, size=rel(0.6))) + # no coord_flip
  theme(legend.position="top")
### Add margin to plot if displaying it directly (without pathwork)
p.priori <- p.priori + theme(plot.margin = unit(c(1,3,1,1), "cm")) # (t, r, b, l) widen margin
p.priori

# The default behavior of beeswarm has changed in version 0.6.0. 
# In versions <0.6.0, this plot would have been dodged on the y-axis.  
# In versions >=0.6.0, grouponX=FALSE must be explicitly set to group on y-axis. 
# Please set grouponX=TRUE/FALSE to avoid this warning and ensure proper axis choice.

file.out <- "figs/fig_celltypepriori.hypothalamus.priori.by_study.jitter.pdf"
ggsave(plot=p.priori, filename=file.out, width=7, height=4)


# ======================================================================= #
# ============================ XXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ XXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ LEFTOVERS =============================== #
# ======================================================================= #
