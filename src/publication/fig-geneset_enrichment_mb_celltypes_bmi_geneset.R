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

library(ggpubr)


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

dataset_prefix <- "mousebrain"
filter.gwas <- "BMI_UKBB_Loh2018"

# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #

### Read LDSC results
file.results <- here("results/prioritization_celltypes--mousebrain.multi_gwas.csv.gz")
df.ldsc_cts <- read_csv(file.results)

# =========================== FILTER GWAS =========================== #
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == filter.gwas)

# =========================== FILTER GWAS =========================== #
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == filter.gwas)


# =========================== ADD METADATA =========================== #

df.metadata <- get_metadata(dataset_prefix)
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation")


# ======================================================================= #
# ========= EXPORT ENRICHMENT RESULT TABLE WITH MB META-DATA ============ #
# ======================================================================= #

df.export <- df.ldsc_cts %>% left_join(df.enrich, by="annotation")
df.export <- df.export %>% select(annotation, 
                                  p.value.enrichment_combined_rare_mendelian_obesity=combined_rare_mendelian_obesity,
                                  p.value.ldsc=p.value, 
                                  Region, Probable_location, Description)
df.export <- df.export %>% mutate(fdr_significant_enrichment=if_else(p.value.enrichment_combined_rare_mendelian_obesity < 0.05/n(), TRUE, FALSE))
df.export <- df.export %>% arrange(p.value.enrichment_combined_rare_mendelian_obesity)
# file.out <- "tables/table-celltype_bmi_geneset_enrichment.mb.csv"
# df.export %>% write_csv(file.out)

# ======================================================================= #
# ============ SCATTERPLOT - S-LDSC VS OBESITY ENRICHMENT ============ #
# ======================================================================= #

df.plot <- df.export

fdr_threshold <- 0.05/nrow(df.plot)
fdr_threshold.mlog10 <- -log10(fdr_threshold)

df.plot <- df.plot %>% mutate(
  significance_group = case_when(
    ((p.value.ldsc < fdr_threshold) & (p.value.enrichment_combined_rare_mendelian_obesity < fdr_threshold)) ~ "Geneset and S-LDSC",
    (p.value.ldsc < fdr_threshold) ~ "S-LDSC (common obesity variation)",
    (p.value.enrichment_combined_rare_mendelian_obesity < fdr_threshold) ~ "Geneset (rare obesity variation)",
    TRUE ~ "ns"
  )
)

df.plot %>% count(significance_group)
# significance_group     n
# 1 both                   2
# 2 enrich                13
# 3 ldsc                   9
# 4 ns                   241
df.plot

p <- ggplot(df.plot, aes(x=-log10(p.value.ldsc), y=-log10(p.value.enrichment_combined_rare_mendelian_obesity)))
p <- p + geom_point(color="gray")
p <- p + geom_point(data=df.plot %>% filter(significance_group != "ns"), aes(color=significance_group))
p <- p + geom_text_repel(data=df.plot %>% filter(significance_group != "ns"), aes(label=annotation), show.legend=F)
p <- p + geom_abline()
p <- p + geom_hline(yintercept=fdr_threshold.mlog10, linetype="dashed", color="gray")
p <- p + geom_vline(xintercept=fdr_threshold.mlog10, linetype="dashed", color="gray")
p <- p + ggpubr::stat_cor(method = "pearson")
p <- p + labs(y=expression(-log[10](P[enrichment])), x=expression(-log[10](P[S-LDSC])), color="")
### Theme
p <- p + theme_classic()
p <- p + theme(legend.position = "top")
p
file.out <- sprintf("figs/fig_celltype_geneset_enrichment_vs_ldsc.mb.pdf")
ggsave(p, filename=file.out, width=7, height=5)

# p <- p + geom_text_repel(data=df.export %>% filter((-log10(p.value.ldsc) > fdr_threshold.mlog10) | (-log10(p.value.enrichment_combined_rare_mendelian_obesity) > fdr_threshold.mlog10)), 
#                          aes(label=annotation), color="black")

# ======================================================================= #
# ============ ENRICHMENT BARPLOT - LDSC FDR CELL-TYPES ONLY ============ #
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




