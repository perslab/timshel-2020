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

wd <- "/raid5/projects/timshel/sc-genetics/sc-genetics/src/GE-mousebrain"
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

file.metadata <- "/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain-agg_L5.metadata.csv"
df.metadata <- read_csv(file.metadata) %>% select(cols_metadata_keep) %>% rename(annotation = ClusterName)


# ======================================================================= #
# ============================== TESTING =============================== #
# ======================================================================= #

# sem_obj.sub <- subset_annotations(sem_obj, c("EPMB","HYPEN","EPEN","EPSC","CHOR"))
# sem_obj.sub <- calc_sem(sem_obj.sub, "ges")
# sem_obj.sub <- calc_sem_wrapper(sem_obj.sub)
# sem_obj.sub <- bin_sems(sem_obj.sub, n_bins=101, threshold_bin_zero=0)
# sem_obj.sub.human <- map_to_human(sem_obj.sub, type_mouse_gene_ids="ensembl")
# sem_obj.sub.human <- set_group_by_annotation_slots(sem_obj.sub.human)
# sem_obj.sub.human <- calc_sem_meta_from_bins(sem_obj.sub.human)
# 
# write_sems(sem_obj, slot="sem", name.dataset="mousebrain.mouse.all")
# write_sems(sem_obj, slot="sem_bin", name.dataset="mousebrain.mouse.all")
# write_sems(sem_obj.human, slot="sem", name.dataset="mousebrain.human.all")
# write_sems(sem_obj.human, slot="sem_bin", name.dataset="mousebrain.human.all")

# ======================================================================= #
# ======================== NEW WORKFLOW (empirical) ===================== #
# ======================================================================= #

# source(sprintf("%s/load_functions.R", dir.sc_genetics_lib)) # load sc-genetics library
# 
# sem_obj_mouse <- create_sem_object("/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain")
# sem_obj_mouse <- exclude_sporadic_expressed_genes(sem_obj_mouse)
# sem_obj.mouse.sub <- subset_annotations(sem_obj_mouse, c("EPMB","HYPEN","EPEN","EPSC","CHOR"))
# sem_obj.sub <- map_to_human(sem_obj.mouse.sub, type_mouse_gene_ids="ensembl")
# # sem_obj.sub <- calc_sem(sem_obj.sub, "ges")
# sem_obj.sub <- calc_sem_wrapper(sem_obj.sub)
# sem_obj.sub <- calc_empirical_pvalues_wrapper(sem_obj.sub, threshold_pval=0.05)

# ======================================================================= #
# ================================ SAVE/LOAD ================================= #
# ======================================================================= #

### load data
load(file="mousebrain.sem_obj.RData") # human

# ======================================================================= #
# ================================ Write SEMs ================================= #
# ======================================================================= #

# write_sems(sem_obj, slot="sem_meta", name.dataset=dataset_prefix)

# ======================================================================= #
# ================== Write multi_geneset_file for LDSC ================= #
# ======================================================================= #


# df_multi_geneset <- write_multi_geneset_file(sem_obj, dataset_prefix="mousebrain_all")
# df_multi_geneset

# ======================================================================= #
# ================================ Adjust binning ================================= #
# ======================================================================= #

# sem_obj.bins11 <- bin_sems(sem_obj, n_bins=11, threshold_pval=0.05)
# sem_obj.bins11 <- set_group_by_annotation_slots(sem_obj.bins11)
# sem_obj.bins11 <- calc_sem_meta_from_bins(sem_obj.bins11, median_only=T)

# ======================================================================= #
# ================================ Post-processing ================================= #
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


### meta-SEM mean: counts per bin
# df.sem_meta.bin_counts <- sem_obj.hier.Class[["Neurons"]][["sem_meta"]][["mean"]] %>% gather(key="annotation", value="sem_meta_bin") %>% group_by(annotation) %>% count(sem_meta_bin)
df.sem_meta.bin_counts <- sem_obj[["sem_meta"]][["mean"]] %>% gather(key="annotation", value="sem_meta_bin") %>% group_by(annotation) %>% count(sem_meta_bin)
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

df.model_sumstats.all_anno <- fit_sems(object=sem_obj, slot="mean", df.magma, df.metadata, exclude_bin_zero=F)
df.model_sumstats.all_anno

### you get the same p-values when multiplying all values by 100
# sem_obj.scaled100 <- sem_obj
# # sem_obj.scaled100[["sem_meta"]][["mean_orig"]] <- sem_obj.scaled100[["sem_meta"]][["mean"]]
# sem_obj.scaled100[["sem_meta"]][["mean"]] <- sem_obj.scaled100[["sem_meta"]][["mean"]]*100
# df.model_sumstats.all_anno.scaled100 <- fit_sems(object=sem_obj, df.magma, df.metadata, exclude_bin_zero=F)
# df.model_sumstats.all_anno.scaled100

df.model_sumstats.all_anno.bin_zero_excl <- fit_sems(object=sem_obj, slot="mean", df.magma, df.metadata, exclude_bin_zero=T)
df.model_sumstats.all_anno.bin_zero_excl

df.model_sumstats.all_anno.tstat <- fit_sems_tstat(object=sem_obj, slot="mean", df.magma, df.metadata)
df.model_sumstats.all_anno.tstat

### boxplot
boxplot_bin_zero(object=sem_obj, slot="mean", df.magma, df.metadata, annotation=NULL)
boxplot_bin_zero(object=sem_obj, slot="mean", df.magma, df.metadata, annotation="TEGLU10")


# ======================================================================= #
# ============================== TMP CALCULATE EMPIRICAL BETA PVAL ================================ #
# ======================================================================= #

### Add genetic data to sem_meta
# df.regression <- get_genetic_and_sem_data(sem_obj, slot="mean", df.magma)

### Run regressions
# print("Running regressions...")
# list.res <- list()
# annotation <- "TEGLU10"
# df.tmp <- df.sem_meta %>% select(ZSTAT, !!rlang::sym(annotation)) # two-column
# for (idx.perm in 1:100) {
#   # <DO SHUFFLING HERE>
#   list.res[[annotation]][[i]] <- lm(as.formula(paste0("ZSTAT ~", annotation)), data = df.tmp)
# }

# ======================================================================= #
# ============================== Fit SEM - HIERARCHICAL ================================ #
# ======================================================================= #

load(file="mousebrain.sem_obj.hier.Class.RData") # sem_obj.hier.Class

list.df.model_sumstats <- list()
list.df.model_sumstats.tstat <- list()
for (name_elem in names(sem_obj.hier.Class)) {
  print(name_elem)
  list.df.model_sumstats[[name_elem]] <- fit_sems(object=sem_obj.hier.Class[[name_elem]], df.magma, df.metadata)
  list.df.model_sumstats.tstat[[name_elem]] <- fit_sems_tstat(object=sem_obj.hier.Class[[name_elem]], df.magma, df.metadata)
}

list.df.model_sumstats[["Neurons"]]
df <- list.df.model_sumstats.tstat[["Neurons"]]

# ======================================================================= #
# ============================== INVESTIGATE HIERARCHICAL ================================ #
# ======================================================================= #

selected_annotation <- "SCINH3"
selected_annotation <- "DEINH3"
df.hier <- sem_obj.hier.Class[["Neurons"]][["sem_meta"]][["mean"]] %>% mutate(gene=sem_obj.hier.Class[["Neurons"]][["genes"]]) %>% select(gene, hier=!!sym(selected_annotation))
df.global <- sem_obj[["sem_meta"]][["mean"]] %>% mutate(gene=sem_obj[["genes"]]) %>% select(gene, global=!!sym(selected_annotation))
print(nrow(df.hier)); print(nrow(df.global)) # check that they have similar number of genes.
df.cmp <- inner_join(df.hier, df.global, by="gene") # *OBS* inner join.
df.cmp
df.cmp <- df.cmp %>% mutate(
  cases = case_when( # Like an if statement, the arguments are evaluated in order, so you must proceed from the most specific to the most general.
    (hier == 0) & (global == 0) ~ "both.zero",
    (hier > 0) & (global) ~ "both.non_zero",
    (hier > 0) ~ "hier.non_zero",
    (global > 0) ~ "global.non_zero"
  ))
df.cmp %>% count(cases)
# DEINH3
# 1 both.non_zero    1884
# 2 both.zero       10753
# 3 global.non_zero  1886
# 4 hier.non_zero     522

library("ggpubr")
ggplot(df.cmp, aes(x=hier, y=global)) + geom_point() + 
  geom_rug(col=rgb(.5,0,0,alpha=.2)) + 
  geom_smooth(method="lm") +
  geom_abline() +
  stat_cor(method = "pearson") +
  # geom_density_2d() + # does not work: Computation failed in `stat_density2d()`: bandwidths must be strictly positive
  labs(title=selected_annotation)

# ======================================================================= #
# ================================ GSEA ================================= #
# ======================================================================= #


do_gprofiler <- function(input, background){
  df.res = tryCatch({ # to avoid gprofiler error "transfer closed with outstanding read data remaining"
    df.res <- gProfileR::gprofiler(input, organism="hsapiens", ordered_query=F, significant=T, custom_bg=background)
    df.res <- df.res %>% mutate(call_status="Successful")
    return(df.res)
  }, warning = function(w) {
    message(sprintf("warning: %s", w))
  }, error = function(e) {
    message(sprintf("error: %s", e))
    return(data.frame(call_status="Failed")) # return single-row data frame
  }, finally = {
    message("done with function do_gprofiler")
  })
  return(df.res)
}

### GO test
# TODO: add function
df.cmp %>% filter(cases == "hier.non_zero") %>% pull(gene) %>% cat() # non neuronal specific
df.cmp %>% filter(cases == "global.non_zero") %>% pull(gene) %>% cat() # synaptic signaling, synapse
x <- df.cmp %>% filter(cases == "both.non_zero") %>% pull(gene) # synaptic signaling, synapse

df.go <- do_gprofiler(input=x, background=sem_obj[["genes"]])

# ======================================================================= #
# ================================ EXPORT results ================================= #
# ======================================================================= #

# df.model_sumstats.all_anno %>% write_csv("out.cell_prioritization.mousebrain.all_annotations.sem_meta_mean.csv")

for (name_elem in names(list.df.model_sumstats)) {
  print(name_elem)
  list.df.model_sumstats[[name_elem]] %>% write_csv(sprintf("out.cell_prioritization.mousebrain.%s.sem_meta_mean.csv", name_elem))
}


# ======================================================================= #
# ================================ PLOT: SEM vs ZSTAT ================================= #
# ======================================================================= #

df.sem_meta <- get_sem_meta(sem_obj, df.magma)

x <- "SCINH2"
x <- "CBGRC"
x <- "VLMC1" # lowest sign.
# x <- "MFOL1" 
x <- "DEINH3"

# p <- ggplot(df.sem_meta, aes(x=!!sym(x), y=ZSTAT)) + labs(title=x)
p <- ggplot(df.sem_meta %>% filter(!!sym(x)!=0), aes(x=!!sym(x), y=ZSTAT)) + labs(title=x) # exclude bin zero from plot
p <- p + geom_point(alpha=0.3) + geom_density_2d()
p
p <- p + geom_smooth(method = "lm", se=T) # auto", "lm", "glm", "gam", "loess" or a function, e.g. MASS::rlm or mgcv::gam, base::lm, or base::loess.
p <- p + geom_smooth(method = "loess", se=T) # auto", "lm", "glm", "gam", "loess" or a function, e.g. MASS::rlm or mgcv::gam, base::lm, or base::loess.
p




# ======================================================================= #
# ================================ PLOT: cell priori ================================= #
# ======================================================================= #


df.plot <- df.model_sumstats.all_anno.tstat
df.plot <- df.model_sumstats.all_anno
df.plot <- df.model_sumstats.all_anno.bin_zero_excl
df.plot <- list.df.model_sumstats.tstat[["Neurons"]]
fdr_threshold <- 0.05/nrow(df.plot)
p <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) +
  geom_point(aes(color=Class, size=estimate)) + 
  # geom_text(data=df.plot %>% filter(fdr_significant), aes(x=name, y=-log10(bp_value), label=ClusterName, color=Class), hjust = 0, nudge_x = 1.5) + 
  ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=Class), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
  geom_hline(yintercept=-log10(fdr_threshold), color="red") + 
  # labs(title=sprintf("%s - %s", name.gwas.clean, name.expression.clean), x="Cell Type", y="-log10(P-value)") + # expression(-log[10]("P-value"))
  labs(x="Cell Type", y="-log10(P-value)") + # expression(-log[10]("P-value"))
  ggtitle(sprintf("%s - %s", "BMI", "Mousebrain")) +
  # scale_size_manual(values=c("TRUE"=3, "FALSE"=1)) +
  
  # guides(size=F, color=guide_legend(title="Tissue category")) + 
  # guides(size=F, color=F) +  # NO LEGEND at all
  # guides(size=F) +  # NO FDR/size legeng
  theme_classic() + 
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #      panel.background = element_blank(), axis.line = element_line(colour = "black")) + # alternative to theme_classic()
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) # remove x-
p
file.out <- sprintf("out.plot.cell_prioritization.mousebrain.%s.sem_meta_mean.color_by_class.pdf", "all_annotation")
ggsave(p, filename=file.out, width=20, height=8)




for (name_elem in names(list.df.model_sumstats)) {
  print(name_elem)
  
  df.plot <- list.df.model_sumstats[[name_elem]]
  fdr_threshold <- 0.05/nrow(df.plot)
  
  
  p <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) +
    geom_point(aes(color=Class, size=estimate)) + 
    # geom_text(data=df.plot %>% filter(fdr_significant), aes(x=name, y=-log10(bp_value), label=ClusterName, color=Class), hjust = 0, nudge_x = 1.5) + 
    ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=Class), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
    geom_hline(yintercept=-log10(fdr_threshold), color="red") + 
    # labs(title=sprintf("%s - %s", name.gwas.clean, name.expression.clean), x="Cell Type", y="-log10(P-value)") + # expression(-log[10]("P-value"))
    labs(x="Cell Type", y="-log10(P-value)") + # expression(-log[10]("P-value"))
    ggtitle(sprintf("%s - %s", "BMI", "Mousebrain")) +
    # scale_size_manual(values=c("TRUE"=3, "FALSE"=1)) +
    
    # guides(size=F, color=guide_legend(title="Tissue category")) + 
    # guides(size=F, color=F) +  # NO LEGEND at all
    # guides(size=F) +  # NO FDR/size legeng
    theme_classic() + 
    # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    #      panel.background = element_blank(), axis.line = element_line(colour = "black")) + # alternative to theme_classic()
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) # remove x-
  file.out <- sprintf("out.plot.cell_prioritization.mousebrain.%s.sem_meta_mean.color_by_class.pdf", name_elem)
  ggsave(p, filename=file.out, width=20, height=8)
  
}

# ======================================================================= #
# ================================ XXXX ================================= #
# ======================================================================= #

# ======================================================================= #
# ================================ XXXX ================================= #
# ======================================================================= #


# ======================================================================= #
# ======================== Hierarchical analysis ======================== #
# ======================================================================= #

# TaxonomyRank2              n
# <chr>                  <int>
#   1 CNS neurons              181
# 2 PNS neurons               33
# 3 CNS glia                  23
# 4 Neural crest-like glia    13
# 5 Vascular cells            10
# 6 Immune cells               5

# TaxonomyRank1      n
# <chr>          <int>
#   1 Neurons          214
# 2 Glia              36
# 3 Vascular cells    10
# 4 Immune cells       5

# Class              n
# <chr>          <int>
#   1 Neurons          214
# 2 Vascular          11
# 3 Astrocytes        10
# 4 Oligos            10
# 5 PeripheralGlia    10
# 6 Ependymal          5
# 7 Immune             5

# ======================================================================= #
# ================= Compare to Linnarson original GES values ============ #
# ======================================================================= #
### Conclusion
# 1: My calculation works. There is a high correlation between my GES and Linnarson GES. The difference is due to normalization
# 2: My GES values are distributed 'close' to normal but with a long tail. The distribution is 'nicer' (less 'spiky', less variance) compared to Linnarson GES.
# 3: The median of my GES for each cell-type seems to be VERY close to 1. This is a nice attribute.

df.ges.linnarson <- read_file_fast("/raid5/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain-agg_L5.enrichment_values.csv.gz")
df.ges.my <- df.ges %>% mutate(gene=sem.ex[["genes"]])
df.ges.my.join <- df.ges.my %>% left_join(df.ges.linnarson, by="gene")
ggplot(df.ges.my.join, aes(x=ENT8.x, y=ENT8.y)) + geom_point()
with(df.ges.my.join, cor.test(x=ENT9.x, y=ENT9.y, method="spearman"))

ggplot(df.ges.my.join, aes(x=ENT8.x)) + geom_density() # my GES are more 'gently' (less spike) distributed
ggplot(df.ges.my.join, aes(x=ENT8.y)) + geom_density() + xlim(0,100)
summary(df.ges.my.join %>% pull(ENT8.x)) # ---> all cell-types have median GES VERY close to 1. (my GES)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.09231  0.95519  0.99979  1.16738  1.23329 17.68845 
summary(df.ges.my.join %>% pull(ENT8.y))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000    0.951    2.185    7.261    6.102 6017.554

df <- df.ges.my.join %>% summarise_if(is.numeric, .funs=funs(mean, median, min, max))
df.x <- as.data.frame(df %>% t()) %>% rownames_to_column(var="stat")



