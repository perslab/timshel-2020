############### SYNOPSIS ###################
# FIT ES models
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
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library


setwd(here("src/magma/"))

# ======================================================================= #
# ================================ PARAMETERS ================================ #
# ======================================================================= #


# gwas_name <- "BMI_UKBB_Loh2018_no_mhc"
# gwas_name <- "BMI_UPDATE_Yengo2018_no_mhc"
# dataset_prefix <- "mousebrain"
# dataset_prefix <- "tabula_muris"

for (dataset_prefix in c("mousebrain", "tabula_muris")) {
  for (gwas_name in c("BMI_UKBB_Loh2018_no_mhc", "BMI_UPDATE_Yengo2018_no_mhc") ) {

  # ======================================================================= #
  # ================================ LOAD ES ================================= #
  # ======================================================================= #
  
  ### load data
  if (dataset_prefix == "mousebrain") {
    load(file=here("src/GE-mousebrain/mousebrain.es_obj.RData"))
  } else if (dataset_prefix == "tabula_muris") {
    load(file=here("src/GE-maca/tabula_muris.es_obj.RData"))
  } else {
    stop("Wrong dataset_prefix")
  }
  
  # ======================================================================= #
  # ============== LOAD MAGMA GWAS gene-based p-values (corrected) ======== #
  # ======================================================================= #
  
  file.magma <- file.path(here("src/magma/out_magma_BMI"), sprintf("%s.resid_correct_all.gsa.genes.out", gwas_name))
  df.magma <- read_table(file.magma, comment = "#")
  
  ### add ensembl ids
  df.magma <- add_ensembl_ids_from_entrez(df.magma, colname_geneids_from="GENE", colname_geneids_to="gene") %>% 
    filter(!is.na(gene)) %>% # exclude un-mapped genes
    select(-GENE)
  
  
  # ======================================================================= #
  # ============================== Fit ES ================================ #
  # ======================================================================= #
  
  df.model_sumstats.all_anno <- fit_ess(object=es_obj, slot="mean", df.magma, df.metadata=NULL, exclude_bin_zero=F)
  df.model_sumstats.all_anno
  
  # df.model_sumstats.all_anno.bin_zero_excl <- fit_ess(object=es_obj, slot="mean", df.magma, df.metadata, exclude_bin_zero=T)
  # df.model_sumstats.all_anno.bin_zero_excl
  # 
  # df.model_sumstats.all_anno.tstat <- fit_ess_tstat(object=es_obj, slot="mean", df.magma, df.metadata)
  # df.model_sumstats.all_anno.tstat
  
  # ======================================================================= #
  # ================================ EXPORT results ================================= #
  # ======================================================================= #
  
  file.out <- sprintf("out.cell_prioritization.%s.%s.es_meta_mean.csv", dataset_prefix, gwas_name)
  df.model_sumstats.all_anno %>% write_csv(file.out)

  }
}
# ======================================================================= #
# ================================ PLOT: cell priori ================================= #
# ======================================================================= #
# REQUIRES metadata (color=color_by_variable)

# df.plot <- df.model_sumstats.all_anno
# fdr_threshold <- 0.05/nrow(df.plot)
# 
# p <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) +
#     geom_point(aes(color=color_by_variable, size=estimate)) + 
#     ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=color_by_variable), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
#     geom_hline(yintercept=-log10(fdr_threshold), color="red") + 
#     labs(x="Cell Type", y="-log10(P-value)") + # expression(-log[10]("P-value"))
#     ggtitle(sprintf("%s - %s", "BMI_Yengo2018", gwas_name)) +
#     theme_classic() + 
#     theme(axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())
# 
# file.out <- sprintf("out.plot.cell_prioritization.%s.%s.es_meta_mean.color_by_class.pdf", dataset_prefix, gwas_name)
# ggsave(p, filename=file.out, width=20, height=8)


