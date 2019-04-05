############### SYNOPSIS ###################
# TMP SCRIPT FOR PLOTTING LDSC AND MAGMA RESULTS

### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)

source(here("src/lib/load_functions.R")) # load sc-genetics library

wd <- "/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain"
setwd(wd)

dataset_prefix <- "mousebrain"

# ======================================================================= #
# ============================== READ METADATA =============================== #
# ======================================================================= #

cols_metadata_keep <- c("ClusterName",
                        "Class",
                        "Description",
                        "NCells",
                        "Neurotransmitter",
                        "Probable_location",
                        "Region",
                        "TaxonomyRank1",
                        "TaxonomyRank2",
                        "TaxonomyRank3",
                        "TaxonomyRank4")

file.metadata <- "/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain-agg_L5.metadata.csv"
df.metadata <- read_csv(file.metadata) %>% select(cols_metadata_keep) %>% rename(annotation = ClusterName)

# df.metadata <- df.metadata %>% mutate(color_by_variable = Class)

# ======================================================================= #
# ================================ LOAD LDSC ================================= #
# ======================================================================= #

# The last column gives a P-value from a one-sided test that the coefficient is greater than zero. 
file.ldsc_cts <- "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/celltypes.mousebrain.all.mean.BMI_Yengo2018.cell_type_results.txt"
df.ldsc_cts <- read_tsv(file.ldsc_cts)
tmp_split <- stringr::str_split_fixed(df.ldsc_cts$Name,pattern="\\.", n=Inf)
df.ldsc_cts <- df.ldsc_cts %>% mutate(sem=tmp_split[,length(tmp_split[1,])],
                                      annotation=tmp_split[,length(tmp_split[1,])-1],
                                      dataset=tmp_split[,length(tmp_split[1,])-2]) %>% 
  select(-Name) # drop Name col now
df.ldsc_cts <- df.ldsc_cts %>% rename(estimate=Coefficient, std.error=Coefficient_std_error, p.value=Coefficient_P_value)
df.ldsc_cts <- df.ldsc_cts %>% mutate(fdr_significant = if_else(p.value <= 0.05/n(), true=T, false=F))
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation") # add meta data
df.ldsc_cts

### export
# file.out <- sprintf("out.cell_prioritization.%s.%s.sem_meta_mean.csv", dataset_prefix, "ldsc_cts")
# df.ldsc_cts %>% write_csv(file.out)

# ======================================================================= #
# ================================ SAVE/LOAD ================================= #
# ======================================================================= #

### load data
load(file="mousebrain.sem_obj.RData") # human

# ======================================================================= #
# ================================ LOAD MAGMA ================================= #
# ======================================================================= #

file.magma <- "/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/BMI_Yengo2018.resid.correct_all.gsa.genes.out"
df.magma <- read_table(file.magma, comment = "#")

### add ensembl ids
df.magma <- add_ensembl_ids_from_entrez(df.magma, colname_geneids_from="GENE", colname_geneids_to="gene") %>% 
  filter(!is.na(gene)) %>% # exclude un-mapped genes
  select(-GENE)

# ======================================================================= #
# ============================== Fit SEM ================================ #
# ======================================================================= #

df.model_sumstats.all_anno <- fit_sems(object=sem_obj, slot="mean", df.magma, df.metadata, exclude_bin_zero=F)
df.model_sumstats.all_anno

df.model_sumstats.all_anno.bin_zero_excl <- fit_sems(object=sem_obj, slot="mean", df.magma, df.metadata, exclude_bin_zero=T)
df.model_sumstats.all_anno.bin_zero_excl

df.model_sumstats.all_anno.tstat <- fit_sems_tstat(object=sem_obj, slot="mean", df.magma, df.metadata)
df.model_sumstats.all_anno.tstat

### boxplot
boxplot_bin_zero(object=sem_obj, slot="mean", df.magma, df.metadata, annotation=NULL)
boxplot_bin_zero(object=sem_obj, slot="mean", df.magma, df.metadata, annotation="TEGLU10")

# ======================================================================= #
# ================================ EXPORT results ================================= #
# ======================================================================= #

list.model_sumstats <- list("all_anno.lm_with_bin_zero" = df.model_sumstats.all_anno,
                            "all_anno.lm_without_bin_zero" = df.model_sumstats.all_anno.bin_zero_excl,
                            "all_anno.ttest" = df.model_sumstats.all_anno.tstat,
                            "all_annot.ldsc" = df.ldsc_cts)


for (name_elem in names(list.model_sumstats)) {
  print(name_elem)
  file.out <- sprintf("out.cell_prioritization.%s.%s.sem_meta_mean.csv", dataset_prefix, name_elem)
  list.model_sumstats[[name_elem]] %>% write_csv(file.out)
}



# ======================================================================= #
# ================================ PLOT: cell priori ================================= #
# ======================================================================= #

color_by_variable <- rlang::sym("Class")
for (name_elem in names(list.model_sumstats)) {
  print(name_elem)
  df.plot <- list.model_sumstats[[name_elem]]
  # df.plot <- list.model_sumstats[["all_annot.ldsc"]]
  
  fdr_threshold <- 0.05/nrow(df.plot)
  p <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) +
    geom_point(aes(color=!!color_by_variable, size=estimate)) + 
    ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=!!color_by_variable), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
    geom_hline(yintercept=-log10(fdr_threshold), color="red") + 
    labs(x="Cell Type", y="-log10(P-value)") + # expression(-log[10]("P-value"))
    ggtitle(sprintf("%s - %s", "BMI_Yengo2018", name_elem)) +
    theme_classic() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  file.out <- sprintf("out.plot.cell_prioritization.%s.%s.sem_meta_mean.color_by_class.pdf", dataset_prefix, name_elem)
  ggsave(p, filename=file.out, width=20, height=8)
}



### LDSC plot | tau estimates
df.plot <- df.ldsc_cts
color_by_variable <- rlang::sym("Class")
fdr_threshold <- 0.05/nrow(df.plot)
p <- ggplot(df.plot, aes(x=annotation, y=estimate)) +
  geom_point(aes(color=!!color_by_variable, size=-log10(p.value))) + 
  ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=estimate, label=annotation, color=!!color_by_variable), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
  #geom_hline(yintercept=-log10(fdr_threshold), color="red") + 
  labs(x="Cell Type", y="Estimate") +
  ggtitle(sprintf("%s - %s", "BMI_Yengo2018", name_elem)) +
  theme_classic() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p

# ======================================================================= #
# ================================ Compare LDSC and MAGMA ================================= #
# ======================================================================= #

df.ldsc_vs_magma <- full_join(df.ldsc_cts %>% select(annotation, p.value, estimate, std.error, fdr_significant), 
                              # df.model_sumstats.all_anno %>% select(-c(r.squared, adj.r.squared, sigma, statistic, df, logLik, AIC, BIC, deviance, df.residual, p.value.orig)), 
                              df.model_sumstats.all_anno %>% select(annotation, p.value, estimate, std.error, fdr_significant),
                              by="annotation", suffix=c(".ldsc", ".magma"))
df.ldsc_vs_magma <- df.ldsc_vs_magma %>% mutate(fdr_significant.both = if_else(fdr_significant.ldsc & fdr_significant.magma, T, F))
df.ldsc_vs_magma <- df.ldsc_vs_magma %>% left_join(df.metadata, by="annotation") # add meta data

### export
file.out <- sprintf("out.cell_prioritization.%s.%s.sem_meta_mean.csv", dataset_prefix, "LDSC_VS_MAGMA")
df.ldsc_vs_magma %>% write_csv(file.out)

### cor
with(df.ldsc_vs_magma, cor.test(x=p.value.ldsc, y=p.value.magma, method="spearman"))

### plot
ggplot(df.ldsc_vs_magma, aes(x=p.value.ldsc, y=p.value.magma)) + geom_point() + geom_abline()
ggplot() + 
  geom_point(data=df.ldsc_vs_magma, aes(x=p.value.ldsc, y=p.value.magma)) + 
  geom_abline() + 
  scale_y_log10() + scale_x_log10() +
  ggrepel::geom_text_repel(data=df.ldsc_vs_magma %>% filter(fdr_significant.ldsc), aes(x=p.value.ldsc, y=p.value.magma, label=annotation, color=!!color_by_variable), hjust = 0, nudge_x = 0) # OBS geom_text_repel


# ======================================================================= #
# ======================== ENRICHMENT TESTS ============================= #
# ======================================================================= #


enrichment_category <- function(category.now, df, colname_test_for_enrichment, colname_statistic) {
  # print(category.now)
  df.now <- df %>% mutate(category.tmp = if_else(!!rlang::sym(colname_test_for_enrichment)==category.now, category.now, "Other"))
  x <- df.now %>% filter(category.tmp==category.now) %>% pull(!!rlang::sym(colname_statistic))
  y <- df.now %>% filter(category.tmp=="Other") %>% pull(!!rlang::sym(colname_statistic))
  wtest <- wilcox.test(x,y , alternative="greater") # alternative hypothesis: x > y.
  # wtest <- wilcox.test(bt_value ~ category.tmp, data = df.now, alternative="greater") # *OBS*: using formula is DANGEROUS (when using a one-sided test) because the levels in 'category.tmp' are ordered alphabetically
  return(wtest$p.value)
}


df.dataset <- df.ldsc_cts
df.dataset <- df.dataset %>% mutate(p.value.mlog10 = -log10(p.value))
colname_statistic <- "p.value.mlog10"
# SWITCH: 'Class' or 'category' - or any other column
# colname_test_for_enrichment <- "Class" 
colname_test_for_enrichment <- "Region"
colname_test_for_enrichment <- "TaxonomyRank2"
colname_test_for_enrichment <- "TaxonomyRank3"
colname_test_for_enrichment <- "TaxonomyRank4"
# colname_test_for_enrichment <- "Neurotransmitter"

df.wtest_res <- data.frame(enrichment=sapply(df.dataset %>% distinct(!!rlang::sym(colname_test_for_enrichment)) %>% pull(), 
                                         enrichment_category, 
                                         df=df.dataset, 
                                         colname_test_for_enrichment,
                                         colname_statistic))
df.wtest_res <- df.wtest_res %>% 
  rownames_to_column(var="category") %>% 
  as.tibble() %>%
  arrange(enrichment) %>% 
  left_join(df.metadata %>% count(!!sym(colname_test_for_enrichment)), by=c(category=colname_test_for_enrichment))
  
df.wtest_res

### PLOT enrichment
# p <- ggplot(df.wtest_res, aes(x=category, y=-log10(enrichment))) + 
#   geom_col() + 
#   geom_col(data=df.wtest_res %>% filter(category %in% c("Medulla", "Pons", "Hypothalamus")), aes(fill=category)) +  # HIGHLIGHT only for category
#   theme_classic() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# p
# file.out = "plot.BMI_yengo.mousebrain-enrichment_region_highlight.pdf"
# # ggsave(p, filename=file.out, width=8, height=5)