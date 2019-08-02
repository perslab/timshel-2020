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
source(here("src/publication/lib-load_pub_lib_functions.R"))

setwd(here("src/publication"))

# ======================================================================= #
# ================== LOAD LDSC CTS RESULTS (multi GWAS) ================ #
# ======================================================================= #

### Read LDSC results
file.results <- here("results/prioritization_celltypes--mousebrain_all.multi_gwas.csv.gz")
df.ldsc_cts <- read_csv(file.results)

# filter.gwas <- c("BMI_UKBB_Loh2018", "SCZ_Pardinas2018")
filter.gwas <- c("BMI_UKBB_Loh2018")

df.ldsc_cts <- df.ldsc_cts %>% filter(gwas %in% filter.gwas)

# ======================================================================= #
# =========================== READ METADATA =========================== #
# ======================================================================= #

### Read
file.metadata <- here("data/expression/mousebrain/mousebrain-agg_L5.metadata.csv")
df.metadata <- read_csv(file.metadata)


### df.metadata %>% count(Neurotransmitter_class)
# 1 Acetylcholine                7
# 2 Ex/In                        4
# 3 Excitatory                 102
# 4 Inhibitory                  78
# 5 Missing neurotransmitter     5
# 6 Monoamines                  12
# 7 Nitric oxide                 6
# 8 Non-neuron                  51

### Process Neurotransmitter_class 
df.metadata <- df.metadata %>% mutate(Neurotransmitter_class = case_when(
  is.na(Neurotransmitter_class) ~ "Undefined neurotransmitter",
  Neurotransmitter_class == "Ex/In" ~ "Undefined neurotransmitter",
  Neurotransmitter_class == "Missing neurotransmitter" ~ "Undefined neurotransmitter",
  TRUE ~ as.character(Neurotransmitter_class)
))

### Join with LDSC
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation")


# ======================================================================= #
# ======================== ENRICHMENT TEST ============================= #
# ======================================================================= #

### SETTINGS
df.dataset <- df.ldsc_cts
### *IMPORTANT*: only do test for Neurons
df.dataset <- df.dataset %>% filter(Class=="Neurons") # subset Neurons only


### ===== Run enrichment test =====
cols2test <- c("Neurotransmitter_class")

list.res.gwas <- list()
for (name_gwas in filter.gwas) {
  list.res.metadata <- list()
  print(name_gwas)
  for (colname_test_for_enrichment in cols2test) {
    print(colname_test_for_enrichment)
    ### Modify data frame
    df.dataset.tmp <- df.dataset
    df.dataset.tmp <- df.dataset.tmp %>% filter(gwas==name_gwas) 
    df.dataset.tmp <- df.dataset.tmp %>% mutate(p.value.mlog10 = -log10(p.value))
    df.dataset.tmp <- df.dataset.tmp %>% filter(!is.na(!!rlang::sym(colname_test_for_enrichment))) # make sure we don't get any NA values
    colname_statistic <- "p.value.mlog10"
    df.wtest_res <- run_enrichment(df.dataset.tmp, colname_test_for_enrichment, colname_statistic, do_empirical_pvalue=F, run_parallel=F) # no empirical pval
    list.res.metadata[[colname_test_for_enrichment]] <- df.wtest_res
  }
  df.wtest_multi_metadata <- bind_rows(list.res.metadata, .id="metadata_col") # bind col2test
  list.res.gwas[[name_gwas]] <- df.wtest_multi_metadata
}
df.wtest_multi_gwas <- bind_rows(list.res.gwas, .id="gwas") # bind gwas
df.wtest_multi_gwas


# ======================================================================= #
# =========================== EXPORT ENRICHMENT ============================ #
# ======================================================================= #

file.out <- "tables/table-transmitter_enrichment.csv"
df.wtest_multi_gwas %>% write_csv(file.out)

# ======================================================================= #
# =========================== PLOT ENRICHMENT ============================ #
# ======================================================================= #

name_gwas <- "BMI_UKBB_Loh2018"
colname_test_for_enrichment <- "Neurotransmitter_class"

df.plot <- df.wtest_multi_gwas %>% filter(gwas==name_gwas & metadata_col==colname_test_for_enrichment)
p <- ggplot(df.plot, aes(x=category, y=-log10(enrichment_val))) +
  geom_col() +
  geom_hline(yintercept=-log10(0.05/nrow(df.plot)), linetype="dashed", color="gray") + 
  labs(x="", y=expression(-log[10](P[enrichment]))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
file.out <- "figs/fig_transmitter_enrichment.barplot.pdf"
ggsave(filename=file.out, plot=p, width=7, height=4)


# =============================================================================== #
# == PLOT RESULTS [ALL TRANSMITTER CLASSES, SORTED BY CLASS] + NON-NEURONS ====== #
# =============================================================================== #

name_gwas <- "BMI_UKBB_Loh2018"
# name_gwas <- "SCZ_Pardinas2018"

### Init
df.plot <- df.ldsc_cts %>% filter(gwas == name_gwas)

### Order Neurotransmitter_class
# unique(df.plot$Neurotransmitter_class)
order.x <- c("Excitatory",
  "Inhibitory",
  "Acetylcholine",
  "Monoamines",
  "Nitric oxide",
  "Undefined neurotransmitter",
  "Non-neuron")
df.plot <- df.plot %>% mutate(Neurotransmitter_class = factor(Neurotransmitter_class, levels=order.x))
df.plot.labels <- df.plot %>% 
  count(Neurotransmitter_class) %>% # this df will be ordered by factor levels of Neurotransmitter_class
  mutate(n_label = paste0(Neurotransmitter_class, " (", n, ")"))
df.plot.labels

### Order annotations
annotations_highlight <- get_prioritized_annotations_bmi(dataset="mousebrain")
df.plot <- df.plot %>% arrange(Neurotransmitter_class, p.value) %>% mutate(annotation = factor(annotation, levels=annotation))
df.plot <- df.plot %>% mutate(flag_highlight = if_else(annotation %in% annotations_highlight, TRUE, FALSE))
fdr_threshold <- 0.05/nrow(df.plot)

p <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) +
  geom_point(aes(color=Neurotransmitter_class)) +
  geom_text_repel(data=df.plot %>% filter(flag_highlight), aes(label=annotation), size=rel(2)) +
  geom_hline(yintercept=-log10(fdr_threshold), color="gray", linetype="dashed") + 
  labs(x="Cell-type", y=expression(-log[10](P[S-LDSC])), color="") + # color="Neurotransmitter class"
  scale_color_brewer(palette="Set2", labels=df.plot.labels$n_label) + 
  theme_classic() + 
  theme(legend.position="bottom") + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) # remove x-labels
p
file.out <- "figs/fig_transmitter_enrichment.celltypes.pdf"
ggsave(p, filename=file.out, width=10, height=4)



# ======================================================================= #
# ============ PLOT MB RESULTS [ONLY EX/IN, SORTED BY P-VAL] ============= #
# ======================================================================= #

# name_gwas <- "BMI_UKBB_Loh2018"
# # name_gwas <- "SCZ_Pardinas2018"
# df.plot <- df.ldsc_cts %>% filter(gwas == name_gwas)
# annotations_highlight <- get_prioritized_annotations_bmi(dataset="mousebrain")
# df.plot <- df.plot %>% mutate(flag_highlight = if_else(annotation %in% annotations_highlight, TRUE, FALSE))
# SELECTION <- c("Ex", "In")
# fdr_threshold <- 0.05/nrow(df.plot)
# 
# p <- ggplot(df.plot, aes(x=fct_reorder(annotation, tax_order_idx_mb_fig1c), y=-log10(p.value))) # order by tax_order_idx_mb_fig1c
# p <- p + geom_point(aes(size=flag_highlight), color="gray50") +
#   geom_point(data=df.plot %>% filter(Neurotransmitter_class %in% SELECTION), aes(size=flag_highlight, color=Neurotransmitter_class)) +
#   ggrepel::geom_text_repel(data=df.plot %>% filter(Neurotransmitter_class %in% SELECTION, flag_highlight), aes(x=annotation, y=-log10(p.value), label=annotation, color=Neurotransmitter_class), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
#   geom_hline(yintercept=-log10(fdr_threshold), color="gray", linetype="dashed") + 
#   labs(x="Cell Type", color="Neurotransmitter class", y=expression(-log[10]("P-value"))) + 
#   theme_classic() + 
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) # remove x-labels
# p

# file.out <- sprintf("plot_transmitter.%s.pdf", name_gwas)
# ggsave(p, filename=file.out, width=15, height=8)


