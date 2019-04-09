############### SYNOPSIS ###################
### AIM: Perform enrichment tests to answer the question:
# Does certain mousebrain metadata (cell-type Region, Neurotransmitter, ...) have non-random genetic prioritization score/rank/values (LDSC results)?


### OUTPUT: 
# ....

### REMARKS:
# Most content is copied and modified from "src/datasets-expression/mousebrain/analysis-genetic_prioritization_enrichment_tests_metadata.R"

### REFERENCE:

# ======================================================================= #
# ================================ SETUP ================================ #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/publication"))

# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
file.results <- here("results/prioritization_celltypes--mousebrain_all.multi_gwas.csv.gz")
df.ldsc_cts <- read_csv(file.results)

filter.gwas <- c("BMI_UKBB_Loh2018", "SCZ_Pardinas2018")

df.ldsc_cts <- df.ldsc_cts %>% filter(gwas %in% filter.gwas)

# ======================================================================= #
# =========================== READ METADATA =========================== #
# ======================================================================= #

### Read
file.metadata <- here("data/expression/mousebrain/mousebrain-agg_L5.metadata.csv")
df.metadata <- read_csv(file.metadata)

### Join with LDSC
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation")

# ======================================================================= #
# ======================== TEST SETTINGS ============================= #
# ======================================================================= #

### *IMPORTANT*: only do test for Neurons
df.dataset <- df.dataset %>% filter(Class=="Neurons") # subset Neurons only

# ======================================================================= #
# ======================== ENRICHMENT TEST ============================= #
# ======================================================================= #


### ===== Run enrichment test =====
cols2test <- c("Neurotransmitter_class")
df.dataset <- df.ldsc_cts
list.res.gwas <- list()
for (name_gwas in filter.gwas) {
  list.res.metadata <- list()
  print(name_gwas)
  for (colname_test_for_enrichment in cols2test) {
    print(colname_test_for_enrichment)
    ### Modify data frame
    df.dataset <- df.dataset %>% filter(gwas==name_gwas) 
    df.dataset <- df.dataset %>% mutate(p.value.mlog10 = -log10(p.value))
    df.dataset <- df.dataset %>% filter(!is.na(!!rlang::sym(colname_test_for_enrichment))) # make sure we don't get any NA values
    colname_statistic <- "p.value.mlog10"
    df.wtest_res <- run_enrichment(df.dataset, colname_test_for_enrichment, colname_statistic, do_empirical_pvalue=F, run_parallel=F) # no empirical pval
    list.res.metadata[[colname_test_for_enrichment]] <- df.wtest_res
  }
  df.wtest_multi_metadata <- bind_rows(list.res.metadata, .id="metadata_col") # bind col2test
  list.res.gwas[[name_gwas]] <- df.wtest_multi_metadata
}
df.wtest_multi_gwas <- bind_rows(list.res.gwas, .id="gwas") # bind gwas
df.wtest_multi_gwas


# ======================================================================= #
# =========================== PLOT ENRICHMENT ============================ #
# ======================================================================= #


### PLOT enrichment_val | FACET GRID
p <- ggplot(df.wtest_multi_gwas, aes(x=category, y=-log10(enrichment_val))) +
  geom_col() +
  facet_grid(gwas~metadata_col)
p

for (name_gwas in filter.gwas) {
  for (name_metadata_col in cols2test) {
    df.plot <- df.wtest_multi_gwas %>% filter(gwas==name_gwas & metadata_col==name_metadata_col)
    # df.plot <- df.plot %>% filter(category!="Ex/In") # ***************** REMOVE THIS LATER *******************
    p <- ggplot(df.plot, aes(x=category, y=-log10(enrichment_val))) +
      geom_col() +
      geom_hline(yintercept=-log10(0.05)) + 
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    file.out <- sprintf("plot_enrichment__%s__%s.pdf", name_metadata_col, name_gwas)
    ggsave(filename=file.out, plot=p)
  }
}


### PLOT enrichment_val | SINGLE
p <- ggplot(df.wtest_res, aes(x=category, y=-log10(enrichment_val))) +
  geom_col() +
  # geom_col(data=df.wtest_res %>% filter(category %in% c("Medulla", "Pons", "Hypothalamus")), aes(fill=category)) +  # HIGHLIGHT only for category
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p



# ======================================================================= #
# ======================== PLOT MB RESULTS ============================= #
# ======================================================================= #

name_gwas <- "BMI_UKBB_Loh2018"

name_gwas <- "SCZ_Pardinas2018"
df.plot <- df.ldsc_cts %>% filter(gwas == name_gwas)
for (sorted_plot in c(T,F)) {
  SELECTION <- c("Ex", "In")
  fdr_threshold <- 0.05/nrow(df.plot)
  if (sorted_plot) {
    p <- ggplot(df.plot, aes(x=fct_reorder(annotation, log10(p.value)), y=-log10(p.value))) # x-axis ranked by p-value
  } else {
    p <- ggplot(df.plot, aes(x=fct_reorder(annotation, tax_order_idx_mb_fig1c), y=-log10(p.value))) # order by tax_order_idx_mb_fig1c
  }
  p <- p + geom_point(aes(size=fdr_significant), color="gray50") +
    geom_point(data=df.plot %>% filter(Neurotransmitter_class %in% SELECTION), aes(size=fdr_significant, color=Neurotransmitter_class)) +
    ggrepel::geom_text_repel(data=df.plot %>% filter(Neurotransmitter_class %in% SELECTION, fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=Neurotransmitter_class), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
    geom_hline(yintercept=-log10(fdr_threshold), color="gray") + 
    labs(x="Cell Type", color="Neurotransmitter class", y=expression(-log[10]("P-value"))) + 
    scale_size_manual(values=c("TRUE"=3, "FALSE"=1)) +
    guides(size=F) +  # NO FDR/size legeng
    theme_classic() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) # remove x-labels
  p
  print(file.out)
  file.out <- sprintf("plot_transmitter.%s%s.pdf", name_gwas, if_else(sorted_plot, ".sorted_plot", ""))
  ggsave(p, filename=file.out, width=15, height=8)
  
}


