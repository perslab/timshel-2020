############### SYNOPSIS ###################
# Analyze and export LDSC results for multiple datasets (tabula muris, mousebrain, ...)
# CONTITIONAL!!!


### OUTPUT: 
# ....

### REMARKS:
# ....

### REFERENCE:

# ======================================================================= #
# =============================== SETUP ================================= #
# ======================================================================= #

library(tidyverse)
library(here)

source(here("src/lib/load_functions.R")) # load sc-genetics library

setwd(here("src/ldsc"))


# ======================================================================= #
# ============================ PARAMETERS =============================== #
# ======================================================================= #


### PARAMS
dataset_prefix <- "mousebrain"
# dataset_prefix <- "tabula_muris"

genomic_annotation_prefix <- get_genomic_annotation_prefix(dataset_prefix)

# ======================================================================= #
# =============================== LOAD LDSC CTS RESULTS ================================= #
# ======================================================================= #


### Load - CONDITIONAL
dir.data <- "/projects/timshel/sc-genetics/sc-genetics/out/out.ldsc/"
filenames <- list.files(path=dir.data,  pattern=sprintf("%s__BMI_UKBB_Loh2018(.*).cell_type_results.txt", genomic_annotation_prefix))
# [12] "celltypes.mousebrain.all__BMI_UKBB_Loh2018__CONDITIONAL__wgcna.mousebrain-190111-dodgerblue.cell_type_results.txt"
# [13] "celltypes.mousebrain.all__BMI_UKBB_Loh2018_no_mhc_max_chisq_720.cell_type_results.txt"                            
# [14] "celltypes.mousebrain.all__BMI_UKBB_Loh2018_no_mhc_max_chisq_80.cell_type_results.txt"                             
# [15] "celltypes.mousebrain.all__BMI_UKBB_Loh2018.cell_type_results.txt"
filenames <- filenames[!grepl(pattern="no_mhc", filenames, perl=T)] # *TMP HACKY* exclude any *no_mhc_max_chisq_* results
list.dfs <- lapply(file.path(dir.data, filenames), load_ldsc_cts_results, dataset_prefix)
# stringr::str_match(filenames, pattern=sprintf("%s__(.*?)(__CONDITIONAL__)?(.*?).cell_type_results.txt", genomic_annotation_prefix))
tmp <- stringr::str_match(filenames, pattern=sprintf("%s__(.*?)(__CONDITIONAL__)(.*?)\\.cell_type_results.txt", genomic_annotation_prefix))[,4]
stopifnot(length(tmp[is.na(tmp)])==1) # we only allow for one match
tmp[is.na(tmp)] <- "baseline"
names(list.dfs) <- tmp
names(list.dfs)
df.ldsc_cts <- list.dfs %>% bind_rows(.id="condition")
df.ldsc_cts <- df.ldsc_cts %>% filter(es=="es_mean")

df.ldsc_cts %>% count(condition)

### format file [*not* compatible with multi-CONDITION plotting]
df.ldsc_cts.export <- df.ldsc_cts %>% select(condition, p.value, annotation) %>% spread(key=condition, value=p.value) # format file

# ======================================================================= #
# =========================== EXPORT to results ========================= #
# ======================================================================= #

### Export (selected columns)
file.out <- here("results", sprintf("prioritization_celltypes_conditional--%s.BMI_UKBB_Loh2018.csv.gz", dataset_prefix))
file.out
# df.ldsc_cts %>% select(-es, -dataset, -n_obs_es, -fdr_significant, -p.value.adj) %>% write_csv(file.out)


# ======================================================================= #
# =============================== METADATA ================================= #
# ======================================================================= #

df.metadata <- get_metadata(dataset_prefix)
df.ldsc_cts <- df.ldsc_cts %>% left_join(df.metadata, by="annotation") # add meta data
df.ldsc_cts.export <- df.ldsc_cts.export %>% left_join(df.metadata, by="annotation") # add meta data

# ======================================================================= #
# =============================== EXPORT ================================= #
# ======================================================================= #

### Multi CONDITION
df.ldsc_cts.export %>% write_csv(sprintf("out.tmp.190212.%s.es_mean.CONDITIONAL_BMI_UKBB_Loh2018.csv", dataset_prefix))

# ======================================================================= #
# ================================ PLOT: cell priori [FACET WRAP MULTI-CONDITIONAL] ================================= #
# ======================================================================= #


df.plot <- df.ldsc_cts
fdr_threshold <- 0.05/nrow(df.plot %>% filter(condition==(df.plot %>% pull(condition))[1])) # *OBS*: TMP DIRTY CODE FOR CONDITIONAL
p <- ggplot(df.plot, aes(x=annotation, y=-log10(p.value))) +
  geom_point(aes(color=color_by_variable, size=estimate)) + 
  ggrepel::geom_text_repel(data=df.plot %>% filter(fdr_significant), aes(x=annotation, y=-log10(p.value), label=annotation, color=color_by_variable), hjust = 0, nudge_x = 1.5) + # OBS geom_text_repel
  geom_hline(yintercept=-log10(fdr_threshold), color="red") + 
  facet_wrap(~condition, ncol = 2) + 
  labs(x="Cell Type", y="-log10(P-value)") + # expression(-log[10]("P-value"))
  # ggtitle(sprintf("%s - %s", "BMI_UKBB_Loh2018", name_elem)) +
  # ggtitle(sprintf("%s - %s", dataset_prefix, name_elem)) +
  theme_classic() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p
file.out <- sprintf("out.tmp.190212.plot.cell_prioritization.%s.%s.es_mean.color_by_class.pdf", dataset_prefix, "CONDITIONAL_BMI_UKBB_Loh2018")
#ggsave(p, filename=file.out, width=15, height=20)
ggsave(p, filename=file.out, width=15, height=20)


### https://github.com/tidyverse/ggplot2/issues/2402
# Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : 
#                      Viewport has zero dimension(s)

# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #



# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #




# ======================================================================= #
# =============================== XXXXX ================================= #
# ======================================================================= #

