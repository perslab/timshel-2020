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
# sem_obj_mouse <- create_sem_object("/projects/timshel/sc-genetics/sc-genetics/data/expression/mousebrain/mousebrain")
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
# ================================ Adjust binning ================================= #
# ======================================================================= #

# sem_obj.bins11 <- bin_sems(sem_obj, n_bins=11, threshold_pval=0.05)
# sem_obj.bins11 <- set_group_by_annotation_slots(sem_obj.bins11)
# sem_obj.bins11 <- calc_sem_meta_from_bins(sem_obj.bins11, median_only=T)


# ======================================================================= #
# ================================ LOAD MAGMA ================================= #
# ======================================================================= #

file.magma <- "/projects/timshel/sc-genetics/sc-genetics/src/magma-fit_models/BMI_Yengo2018.resid.correct_all.gsa.genes.out"
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


