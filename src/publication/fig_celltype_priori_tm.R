############### SYNOPSIS ###################
# Make Tabula Muris main fig


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ================== IMPORTANT LEGEND INFORMATION ======================= #
# ======================================================================= #

# SEE MOUSEBRAIN

# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library


library(RColorBrewer)
library(patchwork)

library(colorspace)

source(here("src/publication/lib-fig_celltype.R"))

setwd(here("src/publication"))

# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #



# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
file.results <- here("results/primary-tabula_muris.multi_gwas.csv.gz")
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
df.metadata <- get_metadata(dataset_prefix="tabula_muris")
df.metadata

### Add meta data
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation") # add meta data
df.ldsc_cts

# ======================================================================= #
# ============================ PROCESS DATA =============================== #
# ======================================================================= #

### Order annotations
df.ldsc_cts <- df.ldsc_cts %>% mutate(annotation=factor(annotation, levels=unique(annotation[order(tissue, cell_type)]))) # order annotations by their taxonomy (from MB Fig1c), and secondary order by their annotation name
# ^ if df.ldsc_cts only contain 1 gwas, then you don't need the levels=unique().
df.ldsc_cts$annotation %>% levels()

### Add fdr_significant flag (within GWAS)
df.ldsc_cts.tmp.summary <- df.ldsc_cts %>% group_by(gwas) %>% summarise(n_obs_gwas=n())
df.ldsc_cts <- left_join(df.ldsc_cts, df.ldsc_cts.tmp.summary, by="gwas")
df.ldsc_cts <- df.ldsc_cts %>% mutate(fdr_significant = if_else(p.value <= 0.05/n_obs_gwas, true=T, false=F))
df.ldsc_cts


# ======================================================================= #
# ======================== PRE-PROCESS ANNOTATIONS ====================== #
# ======================================================================= #

# Self-defined label formatting function
tm_annotation_formatter <- function(annotation_tissue_celltype) {
  # INPUT: annotation_tissue_celltype, string e.g. 'Marrow.Slamf1-negative_multipotent_progenitor_cell'
  # USAGE 1: df %>% mutate
  # USAGE 2: ggplot() + scale_y_discrete(label=tm_annotation_formatter)
  label <- stringr::str_split_fixed(annotation_tissue_celltype, pattern="\\.", n=Inf)[,2]
  label <- stringr::str_replace_all(label, pattern="-", " - ") # add space to any hyphens
  label <- stringr::str_replace_all(label, pattern="_", " ") # convert _ to space
  # label <- stringr::str_to_sentence(label) # title case
  substr(label, 1, 1) <- toupper(substr(label, 1, 1)) # capitalize first character
  return(label)
}


### Create 'plotting' data frame
df.plot.tax_order <- df.ldsc_cts %>%filter(gwas == "BMI_UKBB_Loh2018")  # filter BMI results

### annotation: create clean names
df.plot.tax_order <- df.plot.tax_order %>% mutate(annotation_label_fmt = tm_annotation_formatter(annotation))
df.plot.tax_order

### tissue: rename + remove some text because they contain too few cell-types
df.plot.tax_order %>% count(tissue, sort=T) # ---> potentially filter : n_annotations_in_tax >= 5
df.plot.tax_order <- df.plot.tax_order %>% mutate(tissue = case_when(
  tissue == "Brain_Myeloid" ~ "Brain",
  tissue == "Brain_Non-Myeloid" ~ "Brain",
  TRUE ~ stringr::str_replace(tissue, pattern="_", replacement=" ") # TRUE ~ as.character(tissue)
  )
)

# ======================================================================= #
# ============================ BMI point plot =============================== #
# ======================================================================= #


### Create data frame with taxonomy text data
df.tax_text_position <- df.plot.tax_order %>% 
  group_by(tissue) %>% 
  summarize(n_annotations_in_tax=n()) %>% # count how many annotations in each tax
  left_join(df.metadata %>% select(tissue) %>% distinct(), by="tissue") %>%
  # ^ OBS: we use 'distinct()' because df.taxonomy_metadata contains (intentional) duplicated rows for some combinations of {tissue, order_idx_mb_taxonomy_fig1c}.
  mutate(tissue = factor(tissue, levels=tissue[order(tissue)])) # IMPORTANT: order factor by the SAME ORDER as 'annotation' column is ordered by in df.plot.tax_order
### Add information of the text's position in the plot
df.tax_text_position <- df.tax_text_position %>% 
  arrange(tissue) %>%
  mutate(pos_mid=cumsum(n_annotations_in_tax)-n_annotations_in_tax/2,
         pos_start=cumsum(n_annotations_in_tax)-n_annotations_in_tax,
         pos_end=cumsum(n_annotations_in_tax),
         idx=1:n()) # idx is used for identifying every n'th tax
df.tax_text_position <- df.tax_text_position %>% mutate(flag_draw_rect=if_else(idx %% 2 == 0, TRUE, FALSE)) 


### Set 'tissue_display_name' name
# Some tissues have too few annotations for us to display their name
df.plot.tax_order %>% count(tissue, sort=T) 
df.tax_text_position <- df.tax_text_position %>% mutate(tissue_display_name = case_when(
  tissue == "Tongue" ~ "", 
  TRUE ~ as.character(tissue)
  )
)

### Reverse tax order [only relevant for gganimate direction. x axis is plottet by annotation].
df.plot.tax_order <- df.plot.tax_order %>% mutate(tissue = fct_rev(tissue))
df.tax_text_position <- df.tax_text_position %>% mutate(tissue = fct_rev(tissue))


fdr_threshold <- 0.05/nrow(df.plot.tax_order)
p.main <- ggplot() +
  ### add ggplot 'baselayer'. This makes our 'canvas' and needs to go first (for some reason I don't fully understand...)
  geom_point(data=df.plot.tax_order, aes(x=annotation, y=-log10(p.value)), color="gray") +
  ### add tax info and gray rects (add gray rect before the data points, so avoid them 'covering' the points)
  geom_text(data=df.tax_text_position, aes(x=pos_mid, y=-3, label=tissue_display_name, group=tissue), hjust="right", size=rel(3)) + # group=tissue solves gganimate problem with 'jumps' in the position. It is not needed for a static ggplot
  geom_rect(data=df.tax_text_position %>% filter(flag_draw_rect), aes(xmin=pos_start, xmax=pos_end, ymin=-5, ymax=Inf, fill=tissue), color="gray", alpha=0.1) +
  ### cell-types
  geom_point(data=df.plot.tax_order %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), color=tissue)) +
  ggrepel::geom_text_repel(data=df.plot.tax_order %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation_label_fmt, color=tissue, group=tissue), hjust = 0, nudge_x = 1.5, show.legend=F) + # group is to try to improve gganimate
  ### extra
  geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color="darkgray") + 
  ### axes
  labs(x="", y=expression(-log[10](P))) +
  # coord
  coord_flip(ylim = c( 0, max(-log10(df.plot.tax_order$p.value)) ), # This focuses the y-axis on the range of interest
             clip = 'off') +   # This keeps the labels from disappearing 
  # ^ clip = 'off': it allows drawing of data points anywhere on the plot, including in the plot margins. If limits are set via xlim and ylim and some data points fall outside those limits, then those data points may show up in places such as the axes, the legend, the plot title, or the plot margins.
  # ^ clip = 'off': disable cliping so df.tax_text_position text annotation can be outside of plot | REF https://stackoverflow.com/a/51312611/6639640
  ### guides
  #...
  ### scale
  scale_x_discrete(label=tm_annotation_formatter) +
  ### Color
  #scale_fill_brewer(palette="Paired") +
  # scale_fill_brewer(palette="Dark2") +
  # scale_fill_brewer(palette="Set2") +
  ### theme
  theme_classic() + 
  theme(axis.text.y=element_text(size=rel(0.5))) + # REF: https://stackoverflow.com/a/47144823/6639640
  theme(legend.position="bottom")

### Remove legend
p.main <- p.main + guides(fill=F, color=F)
### Add margin to plot if displaying it directly (without pathwork)
p.main.margin <- p.main + theme(plot.margin = unit(c(1,1,1,5), "cm")) # (t, r, b, l) widen left margin
p.main.margin

# file.out <- sprintf("fig_celltypepriori_tm.pdf")
# ggsave(p.main.margin, filename=file.out, width=9, height=8)


# ======================================================================= #
# =================== HEATMAP + H2 BARPLOT =================== #
# ======================================================================= #

set_multi_gwas_heatmap_plot(df.ldsc_cts)
set_h2_barplot()


# ======================================================================= #
# ============================ TOP3 cell-types ========================= #
# ======================================================================= #

# df.top <- df.ldsc_cts %>% 
#   filter(gwas %in% filter.gwas) %>% 
#   group_by(gwas) %>% 
#   top_n(n=3, wt=-log10(p.value))
# # df.top

# ======================================================================= #
# =================== COMBINE PLOTS: MAIN + HEATMAP =================== #
# ======================================================================= #

### Make blank plot to funciton as a 'dummy margin plot'
p.blank <- ggplot(df.plot.tax_order, aes(x=annotation, y=-log10(p.value))) + 
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
file.out <- sprintf("fig_celltypepriori_tm_with_heatmap_and_barplot.pdf")
# ggsave(p.patch, filename=file.out, width=12, height=10)


### Patch (no h2 barplot)
# p.patch <- p.blank + p.main + p.multi_gwas + plot_layout(ncol = 3, widths = c(6,10,5))
# p.patch
