############### SYNOPSIS ###################
# FIT SEM models
# Export Mousebrain enrichment scores



### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# ==============================  SETUP  =============================== #
# ======================================================================= #

library(tidyverse)

dir.sc_genetics_lib <- "/projects/timshel/sc-genetics/sc-genetics/src/lib/"
source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library


# ======================================================================= #
# ============================== CONSTANTS =============================== #
# ======================================================================= #


wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-hypothalamus"
setwd(wd)

# dataset_prefix <- "campbell_lvl1"
dataset_prefix <- "campbell_lvl2"

# ======================================================================= #
# ============================== READ METADATA =============================== #
# ======================================================================= #

# file.metadata <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/tabula_muris/tabula_muris_facs.tissue_celltype.celltype_metadata.csv"

if (dataset_prefix == "campbell_lvl1") {
  file.metadata <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/hypothalamus_campbell/campbell_lvl1.cell_type_metadata.csv"
  df.metadata <- read_csv(file.metadata) %>% rename(annotation = cell_type_all_lvl1)
} else if (dataset_prefix == "campbell_lvl2") {
  file.metadata <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/hypothalamus_campbell/campbell_lvl2.cell_type_metadata.csv"
  df.metadata <- read_csv(file.metadata) %>% rename(annotation = cell_type_all_lvl2)
} else {
  stop("wrong dataset_prefix")
}

df.metadata <- df.metadata %>% mutate(color_by_variable = taxonomy_lvl1)


# ======================================================================= #
# ================================ SAVE/LOAD ================================= #
# ======================================================================= #

### load data
load(file=sprintf("%s.sem_obj.RData", dataset_prefix)) # human

# ======================================================================= #
# ================================ Adjust binning ================================= #
# ======================================================================= #

# sem_obj.bins11 <- bin_sems(sem_obj, n_bins=11, threshold_pval=0.05)
# sem_obj.bins11 <- set_group_by_annotation_slots(sem_obj.bins11)
# sem_obj.bins11 <- calc_sem_meta_from_bins(sem_obj.bins11, median_only=T)

# ======================================================================= #
# ================================ LOAD MAGMA ================================= #
# ======================================================================= #

file.magma <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/BMI_Yengo2018.resid.correct_all.gsa.genes.out"
df.magma <- read_table(file.magma, comment = "#")

### add ensembl ids
df.magma <- add_ensembl_ids_from_entrez(df.magma, colname_geneids_from="GENE", colname_geneids_to="gene") %>% 
  filter(!is.na(gene)) %>% # exclude un-mapped genes
  select(-GENE)
# [1] "Number of genes mapped: 17467"
# [1] "Number of genes not mapped: 158"
# [1] "Number of genes with a NON-unique mapping (genes with duplicated ensembl gene IDs after mapping): 33"
# [1] "Total mapping stats: 191 genes have no mapping (not mapped + duplicates) out of 17625 input genes."
# [1] "Total genes mapped (non NA genes): 17434"
# [1] "Returning tibble with the column 'gene' added where all gene identifiers unique. Unmapped genes have NA values"

# ======================================================================= #
# ============================== Fit SEM ================================ #
# ======================================================================= #

source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library

df.model_sumstats.all_anno <- fit_sems(object=sem_obj, df.magma, df.metadata, exclude_bin_zero=F)
df.model_sumstats.all_anno

df.model_sumstats.all_anno.bin_zero_excl <- fit_sems(object=sem_obj, df.magma, df.metadata, exclude_bin_zero=T)
df.model_sumstats.all_anno.bin_zero_excl

df.model_sumstats.all_anno.tstat <- fit_sems_tstat(object=sem_obj, df.magma, df.metadata)
df.model_sumstats.all_anno.tstat

### boxplot
# boxplot_bin_zero(object=sem_obj, df.magma, df.metadata, annotation=NULL)
# boxplot_bin_zero(object=sem_obj, df.magma, df.metadata, annotation="TEGLU10")

# ======================================================================= #
# ================================ EXPORT results ================================= #
# ======================================================================= #

list.model_sumstats <- list("all_anno.lm_with_bin_zero" = df.model_sumstats.all_anno,
                            "all_anno.lm_without_bin_zero" = df.model_sumstats.all_anno.bin_zero_excl,
                            "all_anno.ttest" = df.model_sumstats.all_anno.tstat)


for (name_elem in names(list.model_sumstats)) {
  print(name_elem)
  file.out <- sprintf("out.cell_prioritization.%s.%s.sem_meta_median.csv", dataset_prefix, name_elem)
  list.model_sumstats[[name_elem]] %>% write_csv(file.out)
}

# ======================================================================= #
# ================================ PLOT: cell priori ================================= #
# ======================================================================= #

for (name_elem in names(list.model_sumstats)) {
  print(name_elem)
  df.plot <- list.model_sumstats[[name_elem]]
  
  fdr_threshold <- 0.05/nrow(df.plot)
  p <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) +
    geom_point(aes(color=color_by_variable, size=estimate)) + 
    ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=color_by_variable), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
    geom_hline(yintercept=-log10(fdr_threshold), color="red") + 
    labs(x="Cell Type", y="-log10(P-value)") + # expression(-log[10]("P-value"))
    ggtitle(sprintf("%s - %s", "BMI_Yengo2018", name_elem)) +
    theme_classic() + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  file.out <- sprintf("out.plot.cell_prioritization.%s.%s.sem_meta_median.color_by_class.pdf", dataset_prefix, name_elem)
  ggsave(p, filename=file.out, width=20, height=8)
  
}


# ======================================================================= #
# ================================ PLOT: SEM vs ZSTAT ================================= #
# ======================================================================= #

df.sem_meta <- get_sem_meta(sem_obj, df.magma)

x <- "DEINH3"

# p <- ggplot(df.sem_meta, aes(x=!!sym(x), y=ZSTAT)) + labs(title=x)
p <- ggplot(df.sem_meta %>% filter(!!sym(x)!=0), aes(x=!!sym(x), y=ZSTAT)) + labs(title=x) # exclude bin zero from plot
p <- p + geom_point(alpha=0.3) + geom_density_2d()
p
p <- p + geom_smooth(method = "lm", se=T) # auto", "lm", "glm", "gam", "loess" or a function, e.g. MASS::rlm or mgcv::gam, base::lm, or base::loess.
p <- p + geom_smooth(method = "loess", se=T) # auto", "lm", "glm", "gam", "loess" or a function, e.g. MASS::rlm or mgcv::gam, base::lm, or base::loess.
p


# ======================================================================= #
# ================================ SEM Post-processing ================================= #
# ======================================================================= #


### Empirical pval distribution
df.plot <- get_empirical_distribution(sem_obj, sem_name="ges", annotation=NULL) # ACBG (median=1)
df.plot.gather <- df.plot %>% gather(key="distribution", value="sem_value")
df.plot.gather %>% ggplot(aes(x=sem_value, fill=distribution)) + geom_density(alpha=0.2) + xlim(c(-10,10))
df.plot.gather %>% ggplot(aes(x=sem_value, fill=distribution)) + geom_density(alpha=0.2) + xlim(c(10,1000))

### qvals summary
df.summary <- get_empirical_pvalues_summary(sem_obj, threshold=0.05, slot="sem_qvalues")

### pvals summary
df.summary <- get_empirical_pvalues_summary(sem_obj, threshold=0.05, slot="sem_pvalues")
df.summary.stats <- sapply(df.summary %>% select(-annotation), summary) %>% # matrix with rownames from summary (Min., 1st Qu., Median, ..)
  as.data.frame() %>%
  rownames_to_column(var="statistic") %>%
  as.tibble()
df.summary.stats
df.summary.stats %>% select(-statistic) %>% summarise_all(.funs=funs(ratio = max(.)/min(.))) # min/max ratio


### meta-SEM median: counts per bin
# df.sem_meta.bin_counts <- sem_obj.hier.Class[["Neurons"]][["sem_meta"]][["median"]] %>% gather(key="annotation", value="sem_meta_bin") %>% group_by(annotation) %>% count(sem_meta_bin)
df.sem_meta.bin_counts <- sem_obj[["sem_meta"]][["median"]] %>% gather(key="annotation", value="sem_meta_bin") %>% group_by(annotation) %>% count(sem_meta_bin)
df.sem_meta.bin_counts <- df.sem_meta.bin_counts %>% filter(!sem_meta_bin == 0) # remove zero-bin
# plot: single-cell type
ggplot(df.sem_meta.bin_counts %>% filter(annotation == "DEINH3"), aes(x=sem_meta_bin, y=n)) + geom_col() # single cell-type
ggplot(df.sem_meta.bin_counts %>% filter(annotation == "MEINH2"), aes(x=sem_meta_bin, y=n)) + geom_col() # single cell-type
# plot: number of non-zero genes per cell-type (with meta-data)
df.sem_meta.bin_counts.summary <- df.sem_meta.bin_counts %>% 
  group_by(annotation) %>% 
  summarise(sum = sum(n)) %>% 
  left_join(df.metadata, by=c("annotation"))
df.sem_meta.bin_counts.summary %>% group_by(Class) %>% summarise(mean = mean(sum))
df.sem_meta.bin_counts.summary %>% 
  ggplot(aes(x=annotation, y=sum, fill=Class)) + geom_col() + geom_hline(yintercept = mean(df.sem_meta.bin_counts.summary$sum), color="red") + labs(y="number of non-zero bin genes per cell-type")

# plot: mean number of genes per bin
df.sem_meta.bin_counts %>% group_by(sem_meta_bin) %>% summarise(mean = mean(n)) %>% ggplot(aes(x=sem_meta_bin, y=mean)) + geom_col() + labs(y="mean number of genes in bin across all cell-types")



# ======================================================================= #
# ================================ XXXX ================================= #
# ======================================================================= #


# ======================================================================= #
# ================================ XXXX ================================= #
# ======================================================================= #


