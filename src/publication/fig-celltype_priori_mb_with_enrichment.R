############### SYNOPSIS ###################
### AIM: Make hypothalamus plot: Campbell+Mousebrain + BMI geneset cell-type enrichment

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
# ================================ MOUSEBRAIN LDSC ================================ #
# ======================================================================= #

dataset_prefix <- "mousebrain_all"
filter.gwas <- "BMI_UKBB_Loh2018"
var_color_by <- sym("TaxonomyRank4_reduced1")

# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #

### Read LDSC results
file.results <- here("results/prioritization_celltypes--mousebrain_all.multi_gwas.csv.gz")
df.ldsc_cts <- read_csv(file.results)

# =========================== FILTER GWAS =========================== #
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == filter.gwas)

# =========================== ADD METADATA =========================== #

df.metadata <- get_metadata(dataset_prefix)
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation")

### Add 'text' column
text.description <- str_split(df.ldsc_cts$Description, pattern=",") %>% purrr::map_chr(1)
# text.annotation <- df.ldsc_cts$annotation
df.ldsc_cts <- df.ldsc_cts %>% mutate(text=paste(annotation, text.description, sep="\n"))

### Add 'clean name' annotation (dummy for mb)
# df.ldsc_cts <- df.ldsc_cts %>% mutate(annotation=annotation)

# ======================================================================= #
# ================================ PLOT LDSC ================================ #
# ======================================================================= #

df.plot <- df.ldsc_cts
### Add pvalue
df.plot <- df.plot %>% mutate(p.value.mlog10 = -log10(p.value))
### Add enrichment data
df.plot <- df.plot %>% left_join(df.enrich %>% select(annotation, combined_rare_mendelian_obesity), by="annotation")

### Order annotation
df.plot <- df.plot %>% mutate(annotation=factor(annotation, levels=unique(annotation[order(tax_order_idx_mb_fig1c, as.character(annotation))]))) # order annotations by their taxonomy (from MB Fig1c), and secondary order by their annotation name
# df.plot <- df.plot %>% mutate(TaxonomyRank4_reduced1 = factor(TaxonomyRank4_reduced1, levels=TaxonomyRank4_reduced1[order(tax_order_idx_mb_fig1c)]))

### Colormap
factor_levels <- sort(unique(df.plot %>% pull(!!var_color_by)))
colormap <- colorspace::qualitative_hcl(n=length(factor_levels), palette = "Dark 3")
names(colormap) <- factor_levels


### Get tax text
df.tax_text_position <- get_celltype_taxonomy_text_position.mb(df.metadata, df.plot)

### PLOT
fdr_threshold <- 0.05/nrow(df.plot)
p.main <- get_celltype_priori_base_tax_plot.mb(df.plot, df.tax_text_position)
### points
p.main <- p.main + geom_point(data=df.plot, aes(x=annotation, y=p.value.mlog10, fill=!!var_color_by, size=-log10(combined_rare_mendelian_obesity), color=!!var_color_by), shape=21, stroke=0.1)
# ^ shapes 21-24 have both stroke colour and a fill. The size of the filled part is controlled by size, the size of the stroke is controlled by stroke.
p.main <- p.main + geom_point(data=df.plot %>% filter(combined_rare_mendelian_obesity<fdr_threshold), aes(annotation, y=p.value.mlog10, size=-log10(combined_rare_mendelian_obesity), fill=!!var_color_by), stroke=0.8, color="black", shape=21, show.legend=F)
### text/description
p.main <- p.main + ggrepel::geom_text_repel(data=df.plot %>% filter(combined_rare_mendelian_obesity<fdr_threshold), 
                                            aes(x=annotation, y=p.value.mlog10, label=text), 
                                            size=rel(2.5), show.legend=F, 
                                            min.segment.length=0.1)# vjust=-2, hjust=1 #  box.padding=0
# p.main <- p.main + geom_text(data=df.plot %>% filter(combined_rare_mendelian_obesity<fdr_threshold), aes(x=annotation, y=p.value.mlog10, label=text), size=rel(2), hjust=0, nudge_y=0.08, show.legend=F)
# ^ no repel: works semi-ok
### axes
p.main <- p.main + labs(x="", y=expression(-log[10](P[S-LDSC])), size=expression(-log[10](P[enrichment])))
### guides
p.main <- p.main + guides(color=F, fill=F)
### color
p.main <- p.main + scale_color_manual(values=colormap) + scale_fill_manual(values=colormap)
### legend
p.main <- p.main + theme(legend.position="top")
p.main
file.out <- sprintf("figs/fig_celltypepriori.mb.with_enrichment.pdf")
ggsave(p.main, filename=file.out, width=11, height=10)



# =================================================================== #
# =========================== PLOT ENRICHMENT ======================= #
# =================================================================== #

### FILTER to match cell-types in df.plot
# df.enrich.filter <- df.enrich %>% filter(annotation %in% df.plot$annotation)
# df.enrich.filter

fdr_threshold <- 0.05/nrow(df.plot)
p.enrich <- ggplot() +
  geom_col(data=df.plot, aes(x=annotation, y=-log10(combined_rare_mendelian_obesity), fill=!!var_color_by)) +
  ### extra
  geom_hline(yintercept=-log10(fdr_threshold), linetype="dashed", color="darkgray") + 
  ### axes
  labs(x="", y=expression(-log[10](P[enrichment]))) +
  # coord
  coord_flip(clip = 'off') + # This keeps the labels from disappearing 
  ### guides
  guides(fill=F) +
  ### color
  scale_fill_manual(values=colormap) +
  ### theme
  theme_classic() + 
  theme(axis.text.y=element_text(size=rel(0.6)))

p.enrich
file.out <- "figs/fig_celltype_geneset_enrichment.mb.pdf"
ggsave(plot=p.enrich, filename=file.out, width=8, height=8)


### Add margin to plot if displaying it directly (without pathwork)
p.enrich.mod <- p.enrich +
  scale_y_reverse() + 
  theme(axis.ticks.y = element_blank(), # wipe away all x-axis elements
        axis.text.y = element_blank(),
        axis.line.y = element_blank())
p.enrich.mod


# ======================================================================= #
# ================================ PATCHWORK ================================ #
# ======================================================================= #

p.patch <- p.enrich.mod + p.main
p.patch

file.out <- "figs/fig_celltypepriori.mb.priori_with_enrichment.pdf"
ggsave(plot=p.patch, filename=file.out, width=12, height=10)


# ======================================================================= #
# ================================ XXXXX ================================ #
# ======================================================================= #

