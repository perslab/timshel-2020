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

### READ: MB + hypothalamus enrichment results
# file.enrich <- here("results/es_enrichment--bmi_gene_lists.pvals.csv")
# df.enrich <- read_csv(file.enrich)
file.enrich <- here("src/publication/tables/table-es_enrichment.combined.csv")
df.enrich <- read_csv(file.enrich)
df.enrich <- df.enrich %>% filter(dataset=="MSN") # *IMPORTANT* to filter to avoid join problems downstream

# ======================================================================= #
# ============================ MOUSEBRAIN LDSC =========================== #
# ======================================================================= #

### Read LDSC results
file.results <- here("results/cellect_ldsc/prioritization.csv")
df.ldsc_cts <- read_csv(file.results)
df.ldsc_cts <- format_cellect_ldsc_results(df.ldsc_cts)
df.ldsc_cts <- df.ldsc_cts %>% filter(specificity_id == "mousebrain")
df.ldsc_cts <- df.ldsc_cts %>% filter(gwas == "BMI_UKBB_Loh2018")

# ======================================================================= #
# ======================== JOIN LDSC + ENRICH DATA ======================= #
# ======================================================================= #

df.join <- df.ldsc_cts %>% rename(p.value.ldsc=p.value) %>%
  left_join(df.enrich %>% rename(p.value.enrich=p.value), by="annotation")


# ======================================================================= #
# ========= EXPORT ENRICHMENT RESULT TABLE WITH MB META-DATA ============ #
# ======================================================================= #
# *NOT DONE IN 2020*. The code is now slightly incompatible with the 2020 updated code.

# ADD METADATA
# df.metadata <- get_metadata("mousebrain")
# df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation")


# df.export <- df.ldsc_cts %>% left_join(df.enrich, by="annotation")
# df.export <- df.export %>% select(annotation, 
#                                   p.value.enrich=combined_rare_mendelian_obesity,
#                                   p.value.ldsc=p.value, 
#                                   Region, Probable_location, Description)
# df.export <- df.export %>% mutate(fdr_significant_enrichment=if_else(p.value.enrich < 0.05/n(), TRUE, FALSE))
# df.export <- df.export %>% arrange(p.value.enrich)
# file.out <- "tables/table-celltype_bmi_geneset_enrichment.mb.csv"
# df.export %>% write_csv(file.out)

# ======================================================================= #
# ============ SCATTERPLOT - S-LDSC VS OBESITY ENRICHMENT ============ #
# ======================================================================= #

df.plot <- df.join

fdr_threshold <- 0.05/nrow(df.plot)
fdr_threshold.mlog10 <- -log10(fdr_threshold)

df.plot <- df.plot %>% mutate(
  significance_group = case_when(
    ((p.value.ldsc < fdr_threshold) & (p.value.enrich < fdr_threshold)) ~ "Geneset and S-LDSC",
    (p.value.ldsc < fdr_threshold) ~ "S-LDSC (common obesity variation)",
    (p.value.enrich < fdr_threshold) ~ "Geneset (rare obesity variation)",
    TRUE ~ "ns"
  )
)

df.plot %>% count(significance_group) %>% arrange(n)
# significance_group                    n
# 1 Geneset and S-LDSC                    4
# 2 Geneset (rare obesity variation)     11
# 3 S-LDSC (common obesity variation)    18
# 4 ns                                  232
df.plot

p <- ggplot(df.plot, aes(x=-log10(p.value.ldsc), y=-log10(p.value.enrich)))
p <- p + geom_point(color="gray")
p <- p + geom_point(data=df.plot %>% filter(significance_group != "ns"), aes(color=significance_group))
p <- p + geom_text_repel(data=df.plot %>% filter(significance_group != "ns"), aes(label=annotation), show.legend=F, size=2)
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

# p <- p + geom_text_repel(data=df.export %>% filter((-log10(p.value.ldsc) > fdr_threshold.mlog10) | (-log10(p.value.enrich) > fdr_threshold.mlog10)), 
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
df.plot <- df.plot %>% left_join(df.enrich %>% select(annotation, p.value.enrich=p.value), by="annotation")
df.plot <- df.plot %>% mutate(enrichment = -log10(p.value.enrich))

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




