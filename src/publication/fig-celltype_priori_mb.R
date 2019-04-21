############### SYNOPSIS ###################
# Make mousebrain main fig


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ================== IMPORTANT LEGEND INFORMATION ======================= #
# ======================================================================= #

# 1. Pval for heatmap is 'clipped' at 5 (all values about 5 get the same max color)
# 2. h2 barplot is on liability scale for case-control traits
# 3. FDR correction is done per-trait.


# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library


library(RColorBrewer)
library(patchwork)

library(colorspace)

source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))

# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #



# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
file.results <- here("results/prioritization_celltypes--mousebrain_all.multi_gwas.csv.gz")
df.ldsc_cts <- read_csv(file.results)

# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (single GWAS) ================ #
# ======================================================================= #

# gwas <- "BMI_UKBB_Loh2018"
# dataset_prefix <- "mousebrain_all"
# genomic_annotation_prefix <- get_genomic_annotation_prefix(dataset_prefix)
# 
# ### Loading BMI data
# file.ldsc_cts <- sprintf("/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/%s__%s.cell_type_results.txt", genomic_annotation_prefix, gwas)
# df.ldsc_cts <- load_ldsc_cts_results(file.ldsc_cts, dataset_prefix)
# df.ldsc_cts <- df.ldsc_cts %>% filter(sem=="sem_mean")
# 


# ======================================================================= #
# ========================= LOAD CELL-TYPE METADATA ===================== #
# ======================================================================= #


### Read annotation metadata
df.metadata <- get_metadata(dataset_prefix="mousebrain_all")

### Add meta data
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation") # add meta data
# df.ldsc_cts <- df.ldsc_cts %>% left_join(df.taxonomy_metadata, by="TaxonomyRank4") # add taxonomy meta data
df.ldsc_cts

# ======================================================================= #
# ============================ PROCESS DATA =============================== #
# ======================================================================= #

### Order annotations
df.ldsc_cts <- df.ldsc_cts %>% mutate(annotation=factor(annotation, levels=unique(annotation[order(tax_order_idx_mb_fig1c, as.character(annotation))]))) # order annotations by their taxonomy (from MB Fig1c), and secondary order by their annotation name
# ^*IMPORTANT*:  we here order annotation factor in the SAME ORDER 'TaxonomyRank4_reduced1' is ordered later in df.tax_text_position
# ^ if df.ldsc_cts only contain 1 gwas, then you don't need the levels=unique().
df.ldsc_cts$annotation %>% levels()

### Add fdr_significant flag (within GWAS)
df.ldsc_cts.tmp.summary <- df.ldsc_cts %>% group_by(gwas) %>% summarise(n_obs_gwas=n())
df.ldsc_cts <- left_join(df.ldsc_cts, df.ldsc_cts.tmp.summary, by="gwas")
df.ldsc_cts <- df.ldsc_cts %>% mutate(fdr_significant = if_else(p.value <= 0.05/n_obs_gwas, true=T, false=F))
df.ldsc_cts

# ======================================================================= #
# ============================ BMI point plot =============================== #
# ======================================================================= #

### Create 'plotting' data frame
df.plot <- df.ldsc_cts %>% filter(gwas == "BMI_UKBB_Loh2018")  # filter BMI results

### Add pvalue
df.plot <- df.plot %>% mutate(p.value.mlog10 = -log10(p.value))

### Rename
df.plot <- df.plot %>% mutate(Region = case_when(
  Region == "Midbrain dorsal" ~ "Midbrain",
  Region == "Midbrain dorsal,Midbrain ventral" ~ "Midbrain",
  Region == "Hippocampus,Cortex" ~ "Hippocampus/Cortex",
  TRUE ~ as.character(Region))
)


### Get tax text
df.tax_text_position <- get_celltype_taxonomy_text_position.mb(df.metadata, df.plot)

### PLOT
p.main <- get_celltype_priori_base_tax_plot.mb(df.plot, df.tax_text_position)
### add more data
p.main <- p.main + geom_point(data=df.plot, aes(x=annotation, y=-log10(p.value)), color="gray")
p.main <- p.main + geom_point(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), color=Region))
p.main <- p.main + ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=Region), hjust = 0, nudge_x = 1.5, show.legend=F)
p.main <- p.main + scale_color_manual(values=colormap.region)
p.main <- p.main + theme(legend.position="bottom")
p.main

file.out <- sprintf("figs/fig_celltypepriori.mb.pdf")
ggsave(p.main, filename=file.out, width=9, height=8)


# ======================================================================= #
# =================== HEATMAP + H2 BARPLOT =================== #
# ======================================================================= #

set_multi_gwas_heatmap_plot(df.ldsc_cts)
p.multi_gwas
set_h2_barplot()


# ======================================================================= #
# =================== COMBINE PLOTS: MAIN + HEATMAP =================== #
# ======================================================================= #

### Make blank plot to funciton as a 'dummy margin plot'
# REASON: SEE my notes on "Patchwork 'outer' margin issue"
p.blank <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) + 
  geom_blank() + 
  theme(panel.grid = element_blank(),axis.line = element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.background = element_blank())


### Patch - with barplot
p.patch <- 
  (
    # LEFT
    (p.h2.blank / 
              (p.blank + p.main + plot_layout(widths = c(6,10)))
    + plot_layout(heights = c(1,10))
    ) |
     # ) & theme(plot.margin =  unit(c(1,1,1,1), "pt")) | # dont know why "& theme()" does not recurse into nested levels
    # RIGHT
    (p.h2_no_x /
              p.multi_gwas 
    + plot_layout(heights = c(1,10))
    )
  # RIGHT + LEFT
  ) + plot_layout(widths = c(10,3)) & theme(plot.margin =  unit(c(1,1,1,1), "pt"))
p.patch

### Save
file.out <- sprintf("figs/fig_celltypepriori.mb.with_heatmap_and_barplot.pdf")
ggsave(p.patch, filename=file.out, width=12, height=10)


### Patch (no h2 barplot)
# p.patch <- p.blank + p.main + p.multi_gwas + plot_layout(ncol = 3, widths = c(6,10,5))
# p.patch


# # ======================================================================= #
# # ============================ CLASS annotation [WORKS] =============================== #
# # ======================================================================= #
# 
# df.plot.anno_class <- df.plot # copy (important to retain annotation factor order)
# ### Rename
# df.plot.anno_class <- df.plot.anno_class %>% mutate(Class = case_when(
#   Class == "Oligoes" ~ "Oligodendrocytes",
#   Class == "PeripheralGlia" ~ "Peripheral Glia",
#   TRUE ~ as.character(Class))
# )
# 
# ### Class
# p_anno_class <- ggplot(df.plot.anno_class, aes(x=annotation, y=0, fill=Class)) + 
#   geom_tile() + 
#   coord_flip() + 
#   theme_minimal() + 
#   theme(legend.position="bottom") +
#   theme(panel.grid = element_blank(), # remove background, axis ticks and texts
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         panel.background = element_blank()) 
# # p_anno_class
# 
# 
# # ======================================================================= #
# # =================== COMBINE PLOTS: MAIN + ANNO_CLASS [WORKS] ================== #
# # ======================================================================= #
# 
# p_patch <- p_anno_class + p.main + plot_layout(ncol = 2, widths = c(1,10)) # guides = "collect" # "auto", "collect", "keep"
# p_patch
# # plot_annotation(tag_levels = "A")
# 
# file.out <- sprintf("fig_celltypepriori_mb_with_class_anno.pdf")
# # ggsave(p_patch, filename=file.out, width=10, height=8)
# 
# 
# # ======================================================================= #
# # ============================ COLORS =============================== #
# # ======================================================================= #


### Colors
# n_colors <- 12
# display.brewer.pal(n = n_colors, name = 'Paired') # # View a single RColorBrewer palette by specifying its name
# brewer.pal(n = n_colors, name = "Paired") # # Hexadecimal color specification

# colormap.region <- c("Cortex"="#1B9E77",
#                      "Thalamus"="#D95F02",
#                      "Midbrain"="#7570B3",
#                      "Hypothalamus"="#E7298A",
#                      "Hippocampus"="#66A61E",
#                      "Hippocampus/Cortex"="#E6AB02")
# colormap.region <- c("Cortex"="#1F78B4",
#                      "Hippocampus/Cortex"="#A6CEE3",
#                      "Hippocampus"="#33A02C",
#                      "Thalamus"="#FDBF6F",
#                      "Midbrain"="#FF7F00",
#                      "Hypothalamus"="#E31A1C"
#                      )
# colormap.region <- c("Cortex"="#1F78B4",
#                      "Hippocampus/Cortex"="#CAB2D6",
#                      "Hippocampus"="#6A3D9A",
#                      "Thalamus"="#B15928",
#                      "Midbrain"="#FF7F00",
#                      "Hypothalamus"="#E31A1C"
# )


# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #
