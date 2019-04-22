############### SYNOPSIS ###################
### AIM: Mousebrain BMI geneset cell-type enrichment

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
# ========================== CELL-TYPE BMI ENRICHMENT =================== #
# ======================================================================= #

### READ: all MB + Campbell cell-types
file.enrich <- here("results/es_enrichment--bmi_gene_lists.pvals.csv")
df.enrich <- read_csv(file.enrich)

# ======================================================================= #
# ============================ MOUSEBRAIN LDSC =========================== #
# ======================================================================= #

dataset_prefix <- "mousebrain_all"
filter.gwas <- "BMI_UKBB_Loh2018"

# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #

### Read LDSC results
file.results <- here("results/prioritization_celltypes--mousebrain_all.multi_gwas.csv.gz")
df.ldsc_cts <- read_csv(file.results)

# =========================== FILTER GWAS =========================== #
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == filter.gwas)

# =========================== FILTER GWAS =========================== #
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == filter.gwas)


# =========================== ADD METADATA =========================== #

df.metadata <- get_metadata(dataset_prefix)
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation")

# ======================================================================= #
# ================================= PLOT ================================ #
# ======================================================================= #

df.plot <- df.ldsc_cts

### Filter
filter.annotations <- get_prioritized_annotations_bmi(dataset="mousebrain")
df.plot <- df.plot %>% filter(annotation %in% filter.annotations)
# df.plot <- df.plot %>% mutate(flag_highlight = if_else(annotation %in% filter.annotations, TRUE, FALSE))

### Add enrichment data
df.plot <- df.plot %>% left_join(df.enrich %>% select(annotation, combined_rare_mendelian_obesity), by="annotation")
df.plot <- df.plot %>% mutate(enrichment = -log10(combined_rare_mendelian_obesity))

### Get annotation colors
colormap.annotations <- get_color_mapping.prioritized_annotations_bmi(dataset="mousebrain")


### Plot
p <- ggplot(df.plot, aes(x=annotation, y=enrichment, label=annotation))
p <- p + geom_segment(aes(x=annotation, xend=annotation, y=0, yend=enrichment), color="grey", alpha=0.3)
p <- p + geom_point(aes(color=annotation), size=3)
p <- p + labs(x=" ", y=expression(-log[10](P[enrichment])))
p <- p + scale_color_manual(values=colormap.annotations)
p <- p + guides(color=F) # hide legend
p <- p + coord_flip()
### Theme
p <- p + theme_classic()
p
file.out <- sprintf("figs/fig_celltype_geneset_enrichment.mb.bmi_celltypes.pdf")
ggsave(p, filename=file.out, width=4, height=4)

