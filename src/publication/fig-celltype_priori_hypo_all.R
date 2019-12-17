############### SYNOPSIS ###################
### AIM: Make hypothalamus plot + BMI geneset cell-type enrichment

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
# ========================== BMI GENESET ENRICHMENT =================== #
# ======================================================================= #


### READ: enrichment results for MB+hypothalamus cell-types
file.enrich <- here("src/publication/tables/table-es_enrichment.combined.csv") # con
df.enrich <- read_csv(file.enrich)
df.enrich <- df.enrich %>% mutate(annotation_clean=annotation) # copy column | only needed for legacy reasons because of _join() operation

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
df.ldsc_cts.tmp.summary <- df.ldsc_cts %>% group_by(specificity_id) %>% summarise(n_obs=n())
df.ldsc_cts <- left_join(df.ldsc_cts, df.ldsc_cts.tmp.summary, by="specificity_id")
df.ldsc_cts <- df.ldsc_cts %>% mutate(fdr_significant = if_else(p.value <= 0.05/n_obs, true=T, false=F))
df.ldsc_cts



# ======================================================================= #
# ============================ MERGE AND PREP LDSC DATA ====================== #
# ======================================================================= #
list.bind <- list("Mousebrain"=df.mb %>% select(annotation, annotation_clean, p.value, category=Class, text, n_cells=NCells), 
                  "Campbell"=df.cb %>% select(annotation, annotation_clean, p.value, category=taxonomy_lvl2, text, n_cells=n))
# TODO: make clean text column for each dataset before joining.

df.join <- bind_rows(list.bind, .id="dataset")

df.join %>% count(dataset)
# dataset        n
# 1 Campbell      64
# 2 Mousebrain    17
# --> Total = 81

### Add pvalue
df.join <- df.join %>% mutate(p.value.mlog10 = -log10(p.value))

### Add enrichment data
df.join <- df.join %>% left_join(df.enrich %>% select(annotation_clean, combined_rare_mendelian_obesity), by="annotation_clean")

# ======================================================================= #
# ============================ GET DENDROGRAM =========================== #
# ======================================================================= #

source("fig-annotation_clustering_mb_campbell.R")
list.res <- get_mousebrain_campbell_integrated_dendrogram(save_fig=FALSE)
p.dendro.orig <- list.res[["plot"]]
dend <- list.res[["dendrogram"]]


# ======================================================================= #
# ============================ XXXXXXXXX =============================== #
# ======================================================================= #

# ======================================================================= #
# ============================ XXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ XXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ LEFTOVERS =============================== #
# ======================================================================= #


# ### Add pvalue
# df.plot <- df.ldsc_cts
# df.plot <- df.plot %>% mutate(p.value.mlog10 = -log10(p.value))
# df.plot <- df.plot %>% arrange(specificity_id) %>% mutate(annotation=factor(annotation, levels=unique(annotation)))
# p <- ggplot(df.plot, aes(x=annotation, y=p.value.mlog10, color=specificity_id)) + 
#   geom_point(aes(size=fdr_significant))
# p

# pp <- ggplotly(p)
# pp
# htmlwidgets::saveWidget(pp, here("src/publication/figs/tmp-prioritization.html"))
# library(plotly)

