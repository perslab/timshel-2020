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

setwd(here("src/publication"))

# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #



# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
file.results <- here("results/primary-mousebrain_all.multi_gwas.csv.gz")
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
# ============ Calculate tau normalized (LOAD and ADD GWAS h2) ========== #
# ======================================================================= #

M_NUMBER_SNPS <- 5961159 # EUR 1000G, MAF>=5%
### reference SNPs for European LD score estimation were the set of 9,997,231 SNPs with a minor allele count >= 5 from 489 unrelated European individuals in Phase 3 of 1000 Genomes.  ===> This is the number of SNPs in the annot files and in the .bim/.bed files. These are the SNPs used to calculate the ldscores (only some of the SNPs are written because of the --print-snps argument).
### Heritability was partitioned for the set of 5,961,159 SNPs with MAF >= 0.05 ===> This is the number of SNPs in the LDSC baseline model (number of SNPs that have written ldscores for)
### Regression coefficient estimation was performed with 1,217,312 HapMap3 SNPs (SNPs in HapMap3 are used because they are generally well-imputed). ===> regression weights
### SEE ALSO EVERNOTE "SOFTWARE | LDSC - understanding | SNPs in LDSC (weights, ldscores, summarystats...)"

### Load h2 (observed scaled) estimates
file.h2 <- here("results/h2_trait.multi_gwas.csv")
df.h2 <- read_csv(file.h2)

### Join
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.h2, by="gwas")

### Calculate tau norm
### normalized tau = tau/(h2g/M), where h2g/M is the"mean per-SNP heritability"
### [LDSD-SEG]:    The normalized tau can be interpreted as the proportion by which the per-SNP heritability of an average SNP would increase if tau_k were added to it.
### [Benchmarker]: The normalized tau can be interpreted as the change in heritability associated with the value of the annotation increasing from 0 to 1.
df.ldsc_cts <- df.ldsc_cts %>% mutate(tau_norm = estimate/(h2/M_NUMBER_SNPS))


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
df.plot.tax_order <- df.ldsc_cts %>%filter(gwas == "BMI_UKBB_Loh2018")  # filter BMI results

### Rename
df.plot.tax_order <- df.plot.tax_order %>% mutate(Region = case_when(
  Region == "Midbrain dorsal" ~ "Midbrain",
  Region == "Midbrain dorsal,Midbrain ventral" ~ "Midbrain",
  Region == "Hippocampus,Cortex" ~ "Hippocampus/Cortex",
  TRUE ~ as.character(Region))
)

### Create data frame with taxonomy text data
df.tax_text_position <- df.plot.tax_order %>% 
  group_by(TaxonomyRank4_reduced1) %>% 
  summarize(n_annotations_in_tax=n()) %>% # count how many annotations in each tax
  left_join(df.taxonomy_metadata %>% select(TaxonomyRank4_reduced1, tax_order_idx_mb_fig1c) %>% distinct(), by="TaxonomyRank4_reduced1") %>% # add 'tax_order_idx_mb_fig1c' to be able to sort factor correctly.
    # ^ OBS: we use 'distinct()' because df.taxonomy_metadata contains (intentional) duplicated rows for some combinations of {TaxonomyRank4_reduced1, tax_order_idx_mb_fig1c}.
  mutate(TaxonomyRank4_reduced1 = factor(TaxonomyRank4_reduced1, levels=TaxonomyRank4_reduced1[order(tax_order_idx_mb_fig1c)])) # IMPORTANT: order factor by the SAME ORDER as 'annotation' column is ordered by in df.plot.tax_order
### Add information of the text's position in the plot
df.tax_text_position <- df.tax_text_position %>% 
  arrange(TaxonomyRank4_reduced1) %>%
  mutate(pos_mid=cumsum(n_annotations_in_tax)-n_annotations_in_tax/2,
         pos_start=cumsum(n_annotations_in_tax)-n_annotations_in_tax,
         pos_end=cumsum(n_annotations_in_tax),
         idx=1:n()) # idx is used for identifying every n'th tax
df.tax_text_position <- df.tax_text_position %>% mutate(flag_draw_rect=if_else(idx %% 2 == 0, TRUE, FALSE)) 

### Remove some text because they contain too few cell-types
df.tax_text_position %>% arrange(n_annotations_in_tax) # ---> potentially filter : n_annotations_in_tax >= 5
df.tax_text_position <- df.tax_text_position %>% mutate(TaxonomyRank4_reduced1 = case_when(
  TaxonomyRank4_reduced1 == "Other CNS neurons" ~ "",
  # TaxonomyRank4_reduced1 == "Other glia" ~ "",
  TRUE ~ as.character(TaxonomyRank4_reduced1))
)


### Colors
n_colors <- 12
display.brewer.pal(n = n_colors, name = 'Paired') # # View a single RColorBrewer palette by specifying its name
brewer.pal(n = n_colors, name = "Paired") # # Hexadecimal color specification

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
colormap.region <- c("Cortex"="#FF7F00",
                     "Hippocampus/Cortex"="#B15928",
                     "Hippocampus"="#E31A1C",
                     "Thalamus"="#6A3D9A",
                     "Midbrain"="#1F78B4",
                     "Hypothalamus"="#33A02C"
)

fdr_threshold <- 0.05/nrow(df.plot.tax_order)
p.main <- ggplot() +
  ### add ggplot 'baselayer'. This makes our 'canvas' and needs to go first (for some reason I don't fully understand...)
  geom_point(data=df.plot.tax_order, aes(x=annotation, y=-log10(p.value)), color="gray") +
  ### add tax info and gray rects (add gray rect before the data points, so avoid them 'covering' the points)
  geom_text(data=df.tax_text_position, aes(x=pos_mid, y=-0.5, label=TaxonomyRank4_reduced1), hjust="right", size=rel(3)) +
  geom_rect(data=df.tax_text_position %>% filter(flag_draw_rect), aes(xmin=pos_start, xmax=pos_end, ymin=-3, ymax=Inf), color="gray", alpha=0.1) +
  ### cell-types
  geom_point(data=df.plot.tax_order %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), color=Region)) +
  ggrepel::geom_text_repel(data=df.plot.tax_order %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=Region), hjust = 0, nudge_x = 1.5, show.legend=F) +
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
  ### color
  scale_color_manual(values=colormap.region)  +
  # scale_color_brewer(palette="Dark2") +
  # scale_color_viridis_d() + # viridis ---> does not work well
  ### theme
  theme_classic() + 
  theme(axis.text.y=element_text(size=rel(0.3))) + # REF: https://stackoverflow.com/a/47144823/6639640
  theme(legend.position="bottom")
### Add margin to plot if displaying it directly (without pathwork)
p.main.margin <- p.main + theme(plot.margin = unit(c(1,1,1,10), "cm")) # (t, r, b, l) widen left margin
p.main.margin

# file.out <- sprintf("fig_celltypepriori_mb.pdf")
# ggsave(p, filename=file.out, width=9, height=8)



# ======================================================================= #
# ============================ MULTI-GWAS PLOT =============================== #
# ======================================================================= #

filter.gwas <- c("BMI_UKBB_Loh2018",
#"AN_PGC_Duncan2017",
"EA3_Lee2018",
"INTELLIGENCE_Savage2018",
"SCZ_Pardinas2018",
"INSOMNIA_Jansen2018",
"MS_Patsopoulos2011",
"WHRadjBMI_UKBB_Loh2018",
"HEIGHT_UKBB_Loh2018",
"LIPIDS_LDL_Willer2013",
"RA_Okada2014")
# ALTERNATIVE Immune: LUPUS_2015

### Create 'plotting' data frame
df.plot.multi_gwas <- df.ldsc_cts %>% filter(gwas %in% filter.gwas)
df.plot.multi_gwas <- df.plot.multi_gwas %>% mutate(gwas=factor(gwas, levels=filter.gwas)) # Order GWAS
df.plot.multi_gwas


### Plot
p.multi_gwas <- ggplot(df.plot.multi_gwas, aes(x=gwas, y=annotation, fill=-log10(p.value))) + 
# p.multi_gwas <- ggplot(df.plot.multi_gwas, aes(x=annotation, y=gwas, fill=tau_norm)) + 
  geom_tile() + 
  geom_text(data=df.plot.multi_gwas %>% filter(fdr_significant), label="*", hjust=0.5, vjust=0.75) + # add asterisk if fdr significant
  # ^ hjust/vjust: https://stackoverflow.com/questions/7263849/what-do-hjust-and-vjust-do-when-making-a-plot-using-ggplot
  # ^ hjust="center", vjust="middle"
  # scale_fill_viridis_c(option="magma", direction=-1) + 
  # scale_fill_distiller(palette="Greys", direction=1) + # "Greys" "Blues" "Greens" "BluGn" "Reds"
  scale_fill_distiller(palette="Blues", direction=1, limits=c(0,5), oob=scales::squish) + # "Greys" "Blues" "Greens" "BluGn" "Reds"
  # ^ oob: Function that handles limits outside of the scale limits (out of bounds). scales::squish "squishes" values into the range specified by limits
  # scale_fill_distiller(palette="Greys", direction=1, limits=c(0,1.8), na.value = "white") + # tau norm plot
  # coord_flip() + 
  labs(fill=expression(-log[10](P))) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(legend.position="bottom") +
  theme(axis.text.y=element_blank(), 
        axis.ticks.y = element_blank(),
        axis.title = element_blank())
p.multi_gwas





# ======================================================================= #
# ================================ h2 barplot =========================== #
# ======================================================================= #

filter.gwas <- c("BMI_UKBB_Loh2018",
                 "EA3_Lee2018",
                 "INTELLIGENCE_Savage2018",
                 "SCZ_Pardinas2018_liability_scale", # OBS
                 "INSOMNIA_Jansen2018_liability_scale", # OBS
                 "MS_Patsopoulos2011_liability_scale", # OBS
                 "WHRadjBMI_UKBB_Loh2018",
                 "HEIGHT_UKBB_Loh2018",
                 "LIPIDS_LDL_Willer2013",
                 "RA_Okada2014_liability_scale") # OBS

### Create 'plotting' data frame
df.plot.h2 <- df.h2 %>% filter(gwas %in% filter.gwas)
df.plot.h2 <- df.plot.h2 %>% mutate(gwas=factor(gwas, levels=filter.gwas)) # Order GWAS
df.plot.h2

### Plot
p.h2 <- ggplot(df.plot.h2, aes(x=gwas, y=h2)) + 
  geom_col(fill="gray") +
  labs(y=expression(h[SLDSC]^{2})) +
  theme(axis.text.x=element_text(angle=45, hjust=1), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line = element_blank())
p.h2
p.h2_no_x <- p.h2 + theme(axis.text.x=element_blank())

### Make blank
p.h2.blank <- ggplot(df.plot.h2, aes(x=gwas, y=h2)) + 
  geom_blank() +
  theme(panel.grid = element_blank(),axis.line = element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.background = element_blank())
p.h2.blank
# ======================================================================= #
# =================== COMBINE PLOTS: MAIN + HEATMAP =================== #
# ======================================================================= #

### Make blank plot to funciton as a 'dummy margin plot'
p.blank <- ggplot(df.plot.tax_order, aes(x=annotation, y=-log10(p.value))) + 
  geom_blank() + 
  theme(panel.grid = element_blank(),axis.line = element_blank(),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),panel.background = element_blank())

### Patch (no h2 barplot)
p.patch <- p.blank + p.main + p.multi_gwas + plot_layout(ncol = 3, widths = c(6,10,5))
p.patch


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
file.out <- sprintf("fig_celltypepriori_mb_with_heatmap_and_barplot.pdf")
# ggsave(p.patch, filename=file.out, width=12, height=10)


### Simple syntax (but needs the power of plot_layout()
# (p.h2.blank / (p.blank | p.main) + plot_layout(heights = c(1,10))) | (p.h2 / p.multi_gwas)

### Patch (with h2 barplot) - works but not aligned
# p.patch <- p.blank + p.main + 
#   (p.h2.blank + p.multi_gwas + plot_layout(ncol = 1, heights = c(1,10)) ) +
#   plot_layout(ncol = 3, widths = c(6,10,5)) & theme(plot.margin =  unit(c(1,1,1,1), "pt"))
# p.patch


### WORKS BUT WRONG
# p.patch <- ( ( (p.blank + p.main + plot_layout(ncol = 2, widths = c(6,10)))
#                + p.h2.blank + plot_layout(ncol = 1, heights = c(1,10)) ) +
#                (p.h2_no_x + p.multi_gwas + plot_layout(ncol = 1, heights = c(1,10)) ) ) +
#   plot_layout(ncol = 2, widths = c(10,5)) & theme(plot.margin =  unit(c(1,1,1,1), "pt"))
# p.patch



### Patchwork 'outer' margin issue
### For some reason, patchwork ignores the 'outside' plot margin for the plots. 
### The only plot.margin patchwork seems to respond to, is the 'inner' margin between plots (e.g. right margin for left plot or left margin for right plot).
### (Top and bottom margins are also ignored.)
### This means that adjusting plot.margin does not help us.
# pm <- p.main + theme(plot.margin = unit(c(1,1,1,200), "pt"))
# pm_multi_gwas <- p.multi_gwas + theme(plot.margin = unit(c(1,1,1,1), "pt"))
# pm + pm_multi_gwas + plot_layout(ncol = 2, widths = c(3,1))


# ======================================================================= #
# ============================ CLASS annotation =============================== #
# ======================================================================= #

df.plot.anno_class <- df.plot.tax_order # copy (important to retain annotation factor order)
### Rename
df.plot.anno_class <- df.plot.anno_class %>% mutate(Class = case_when(
  Class == "Oligoes" ~ "Oligodendrocytes",
  Class == "PeripheralGlia" ~ "Peripheral Glia",
  TRUE ~ as.character(Class))
)

### Class
p_anno_class <- ggplot(df.plot.anno_class, aes(x=annotation, y=0, fill=Class)) + 
  geom_tile() + 
  coord_flip() + 
  theme_minimal() + 
  theme(legend.position="bottom") +
  theme(panel.grid = element_blank(), # remove background, axis ticks and texts
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank()) 
# p_anno_class


# ======================================================================= #
# =================== COMBINE PLOTS: MAIN + ANNO_CLASS ================== #
# ======================================================================= #

p_patch <- p_anno_class + p.main + plot_layout(ncol = 2, widths = c(1,10)) # guides = "collect" # "auto", "collect", "keep"
p_patch
# plot_annotation(tag_levels = "A")

file.out <- sprintf("fig_celltypepriori_mb_with_class_anno.pdf")
# ggsave(p_patch, filename=file.out, width=10, height=8)




# ======================================================================= #
# ============================ COLORS =============================== #
# ======================================================================= #
display.brewer.all(colorblindFriendly = TRUE)
n_colors <- 6
display.brewer.pal(n = n_colors, name = 'Dark2') # # View a single RColorBrewer palette by specifying its name
brewer.pal(n = n_colors, name = "Dark2") # # Hexadecimal color specification

# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #


# ======================================================================= #
# ============================ XXXXXXXXXX =============================== #
# ======================================================================= #
